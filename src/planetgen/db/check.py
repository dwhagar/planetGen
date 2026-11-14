# planetgen/db/check.py

"""
The database check (DB.8)
=========================

`planetgen check-db` (and the Admin dashboard's "Check the database"
button, which runs it as a job) looks for damage in a galaxy database and
changes nothing. See `docs/design/operations/db-check-and-parity-repair.md` section 6.

Each check is a command of its own on the work queue (`check-revision`,
`check-tables`, ... `queue/commands.py`, PERF.126) whose unit the engine runs
as a list of pieces (OPS.204, `queue/engine.py`): a table, a foreign key, a
sector, a star system (`Check.pieces`, `Check.piece`). The unit's result is the
`CheckResult` that `Check.combine` makes of what the pieces found (`Problem`s,
none when the check passes): a check that itself failed to run is an "error"
result (could not check, which is not the same as damaged), never an exception.
Every piece opens its own read-only connection. `planetgen check-db`
only builds the list of the checks asked for and hands it to the wrapper; `Report`
turns the results of the list into the lines and the exit code.

The checks (each one command):

1. `revision`: the Alembic revision, `schema_migrations` and the code's
   `store.SCHEMA_VERSION` agree, and the tables match `db/models.py`.
2. `tables`: `CHECK TABLE ... QUICK` on every table; `tables-extended` is the
   slow `EXTENDED` one (a check of its own, not a flag).
3. `orphans`: every foreign key's child rows have a parent (the keys come
   from `information_schema`), plus the pointers no key covers
   (`nearest_systems` and `sector_paths` objects).
4. `ids`: `id_blocks` is ahead of every id in the tables it hands out.
5. `values`: columns hold plausible numbers (a radius above zero, a phase
   in 0 to 360 degrees, an eccentricity below 1).
6. `counts`: a warning, not damage: `sector_stats.actual_systems` differs
   from the star systems filed in the sector (an orbit refile moves a
   system without updating the stats).
7. `version keys`: every sector has a well formed version key.
8. `systems` (a check of its own, DB.21): every star system is loaded and run
   through `validation.check_star_system` (orbit spacing, the stable-orbit
   ceiling, each planet's moons). Slow, so `deep_estimate` tells how long
   first, from the speed this server recorded for it, and the unit reports its
   progress.

The scope (`--shell` and `--sector`) belongs to the value, count, version-key and
systems checks, which look only at those sectors; the others always look at the
whole database. The systems check also takes `wide_binaries`. A galaxy too big
to check in one go is checked a region at a time.
"""

import re
import time
import warnings
from dataclasses import dataclass, field

import sqlalchemy
from alembic.autogenerate import compare_metadata
from alembic.migration import MigrationContext

from planetgen.db import alembic_runner, models, storage_checksum, store
from planetgen.generation import steps, validation
from planetgen.util import log

MAX_ROWS_SHOWN = 10
"""int: How many offending rows a problem lists before "and N more"."""

EXIT_OK = 0
EXIT_DAMAGED = 1
EXIT_UNCHECKED = 2
"""Exit codes: nothing wrong (warnings allowed), damage found, or a check
could not run (and found no damage)."""

PASS, FAIL, WARN, ERROR = "pass", "fail", "warn", "error"

VALUE_RULES = (
    # (table, column, low, high, low_inclusive, high_inclusive)
    ("stars", "mass_kg", 0, None, False, None),
    ("stars", "radius_km", 0, None, False, None),
    ("stars", "temperature_k", 0, None, False, None),
    ("stars", "age_gy", 0, None, True, None),
    ("stars", "galactic_orbital_phase_deg", 0, 360, True, True),
    ("planets", "radius_km", 0, None, False, None),
    ("planets", "mass_kg", 0, None, False, None),
    ("planets", "period_years", 0, None, False, None),
    ("planets", "surface_temperature_k", 0, None, False, None),
    ("planets", "orbital_inclination_deg", 0, 180, True, True),
    ("moons", "radius_km", 0, None, False, None),
    ("moons", "mass_kg", 0, None, False, None),
    ("moons", "period_years", 0, None, False, None),
    ("moons", "surface_temperature_k", 0, None, False, None),
    ("moons", "orbital_inclination_deg", 0, 180, True, True),
    ("planets", "dipole_tilt_deg", 0, 90, True, True),
    ("moons", "dipole_tilt_deg", 0, 90, True, True),
    ("planets", "core_ni_fraction", 0, 1, True, True),
    ("planets", "core_light_fraction", 0, 1, True, True),
    ("planets", "core_s_fraction", 0, 1, True, True),
    ("planets", "core_si_fraction", 0, 1, True, True),
    ("planets", "core_o_fraction", 0, 1, True, True),
    ("planets", "core_c_fraction", 0, 1, True, True),
    ("planets", "core_h_fraction", 0, 1, True, True),
    ("moons", "core_ni_fraction", 0, 1, True, True),
    ("moons", "core_light_fraction", 0, 1, True, True),
    ("moons", "core_s_fraction", 0, 1, True, True),
    ("moons", "core_si_fraction", 0, 1, True, True),
    ("moons", "core_o_fraction", 0, 1, True, True),
    ("moons", "core_c_fraction", 0, 1, True, True),
    ("moons", "core_h_fraction", 0, 1, True, True),
    ("star_systems", "binary_eccentricity", 0, 1, True, False),
    ("star_systems", "binary_mutual_orbital_phase_deg", 0, 360, True, True),
    ("sectors", "galactic_radius_pc", 0, None, True, None),
)
"""tuple: `(table, column, low, high, low_inclusive, high_inclusive)`; a
NULL value is not checked, a `None` bound is open. Deliberately loose:
each is a value no generator can write, not one that is merely unusual."""

_SECTOR_OF = {
    "sectors": "t.id",
    "star_systems": "t.sector_id",
    "stars": "(SELECT s.sector_id FROM star_systems s WHERE s.id = t.star_system_id)",
    "planets": "(SELECT s.sector_id FROM star_systems s WHERE s.id = t.star_system_id)",
    "moons": "(SELECT s.sector_id FROM star_systems s WHERE s.id = t.star_system_id)",
}
"""dict: For each table with a value rule, the SQL giving its row's sector
id (the row is aliased `t`)."""

_VERSION_KEY = re.compile(r"^[0-9A-Fa-f]{22}$")


@dataclass
class Problem:
    """One thing wrong: `what` is a sentence, `rows` the offenders (each a
    short text), `warning` marks something that is not damage."""
    what: str
    rows: list = field(default_factory=list)
    warning: bool = False

    def as_dict(self):
        return {"what": self.what, "rows": list(self.rows), "warning": self.warning}

    @classmethod
    def from_dict(cls, data):
        return cls(data["what"], list(data["rows"]), bool(data["warning"]))

    def lines(self):
        shown = self.rows[:MAX_ROWS_SHOWN]
        out = [self.what] + [f"    {row}" for row in shown]
        if len(self.rows) > len(shown):
            out.append(f"    ... and {len(self.rows) - len(shown)} more")
        return out


@dataclass
class CheckResult:
    """The outcome of one check: `status` is pass, fail, warn or error."""
    name: str
    status: str
    problems: list = field(default_factory=list)
    error: str = ""

    def as_dict(self):
        """The result as plain JSON, as a unit of the work queue returns it."""
        return {"name": self.name, "status": self.status, "problems": [p.as_dict() for p in self.problems],
                "error": self.error}

    @classmethod
    def from_dict(cls, data):
        return cls(data["name"], data["status"], [Problem.from_dict(p) for p in data["problems"]], data["error"])

    def line(self):
        label = {PASS: "pass", FAIL: "FAIL", WARN: "warn", ERROR: "COULD NOT CHECK"}[self.status]
        return f"{label}: {self.name}" + (f" ({self.error})" if self.error else "")


@dataclass
class Report:
    """Every check's result."""
    results: list

    @property
    def damaged(self):
        return any(r.status == FAIL for r in self.results)

    @property
    def unchecked(self):
        return any(r.status == ERROR for r in self.results)

    def exit_code(self):
        if self.damaged:
            return EXIT_DAMAGED
        return EXIT_UNCHECKED if self.unchecked else EXIT_OK

    def lines(self):
        out = []
        for result in self.results:
            for problem in result.problems:
                out.append(f"[{result.name}] " + ("warning: " if problem.warning else "") + problem.what)
                out.extend(problem.lines()[1:])
        out.extend(r.line() for r in self.results)
        if self.damaged:
            out.append("The database is DAMAGED.")
        elif self.unchecked:
            out.append("No damage found, but some checks could not run.")
        else:
            out.append("The database passed every check.")
        return out


@dataclass
class Scope:
    """Which sectors the sector-scoped checks look at: whole `shells` and
    single `addresses` (`(shell, band, wedge)`); neither means all.
    `wide_binaries` narrows the `--deep` star system pass to the wide
    binaries (GEN.228); it does not change which sectors are looked at."""
    shells: tuple = ()
    addresses: tuple = ()
    wide_binaries: bool = False

    @property
    def everything(self):
        return not self.shells and not self.addresses

    def where(self, sector_sql):
        """Returns `(sql, params)`: a condition on `sector_sql` (empty for
        everything)."""
        if self.everything:
            return "", ()
        clauses, params = [], []
        if self.shells:
            clauses.append("shell_index IN (%s)" % ",".join(["%s"] * len(self.shells)))
            params.extend(self.shells)
        for address in self.addresses:
            clauses.append("(shell_index = %s AND band_index = %s AND wedge_index = %s)")
            params.extend(address)
        return f" AND {sector_sql} IN (SELECT id FROM sectors WHERE {' OR '.join(clauses)})", tuple(params)


@dataclass
class _Context:
    conn: object
    config: object
    scope: Scope

    def rows(self, sql, params=()):
        return self.conn.execute(sql, params).fetchall()

    def tables(self):
        return [row["t"] for row in self.rows(
            "SELECT table_name AS t FROM information_schema.tables"
            " WHERE table_schema = DATABASE() AND table_type = 'BASE TABLE' ORDER BY 1")]


def _table_exists(ctx, name):
    return name in {t.lower() for t in ctx.tables()}


def _dicts(problems):
    """What a piece hands the engine: its problems as plain JSON."""
    return [problem.as_dict() for problem in problems]


def _problems_of(outputs):
    """The usual `combine` of a check: the problems of every piece, in the order of the pieces."""
    return [Problem.from_dict(data) for output in outputs for data in output]


# --- revision: the records of the version, then the tables against the models ------------------------------------------------------

def _revision_pieces(ctx):
    return ["versions", "models"]


def _revision_piece(ctx, piece):
    """The three records of the schema version agree with the code (`versions`), and the tables match the models (`models`)."""
    if piece == "models":
        drift = _model_drift(ctx.config)
        return _dicts([Problem("The tables differ from db/models.py:", drift)] if drift else [])
    problems = []
    head = store.SCHEMA_VERSION
    expected = alembic_runner.revision_id(head)
    tables = {t.lower() for t in ctx.tables()}
    recorded = None
    if "schema_migrations" in tables:
        row = ctx.rows("SELECT MAX(version) AS v FROM schema_migrations")[0]
        recorded = row["v"]
    revision = None
    if "alembic_version" in tables:
        found = [row["version_num"] for row in ctx.rows("SELECT version_num FROM alembic_version")]
        revision = found[0] if len(found) == 1 else found
    if recorded is None:
        problems.append(Problem("schema_migrations has no version recorded."))
    elif recorded != head:
        problems.append(Problem(f"schema_migrations says version {recorded}; this code is version {head}."))
    # A database built whole from schema.sql was never stamped: an empty
    # table is normal, a different or doubled revision is not.
    if revision not in (None, [], expected):
        problems.append(Problem(f"alembic_version says {revision!r}; version {head} is revision {expected!r}."))
    return _dicts(problems)


def _model_drift(config):
    engine = alembic_runner._engine(config)
    try:
        with engine.connect() as connection:
            context = MigrationContext.configure(connection, opts={"compare_type": True})
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", sqlalchemy.exc.SAWarning)  # the star/remnant table cycle
                diffs = compare_metadata(context, models.metadata)
    finally:
        engine.dispose()
    return [_diff_text(diff) for diff in _flatten(diffs)]


def _flatten(diffs):
    for diff in diffs:
        if isinstance(diff, list):
            yield from diff
        else:
            yield diff


def _diff_text(diff):
    kind, *rest = diff
    parts = []
    for item in rest:
        name = getattr(item, "name", None)
        table = getattr(item, "table", None)
        parts.append(f"{getattr(table, 'name', '')}.{name}".strip(".") if name else str(item))
    return f"{kind} {' '.join(parts)}"


# --- tables: one table a piece -----------------------------------------------------------------------------------------------------

def _tables_pieces(ctx):
    return ctx.tables()


def _tables_piece(ctx, table):
    """`CHECK TABLE ... QUICK` on one table."""
    return _check_table(ctx, table, "QUICK")


def _tables_extended_piece(ctx, table):
    """`CHECK TABLE ... EXTENDED` on one table: slow, so a check of its own."""
    return _check_table(ctx, table, "EXTENDED")


def _check_table(ctx, table, mode):
    problems = []
    for row in ctx.rows(f"CHECK TABLE `{table}` {mode}"):
        kind = str(row.get("Msg_type", "")).lower()
        if kind in ("error", "corrupt") or (kind == "status" and str(row.get("Msg_text")).upper() != "OK"):
            problems.append(Problem(f"CHECK TABLE {table}: {row.get('Msg_text')}"))
    return _dicts(problems)


# --- orphans: one foreign key a piece, then the pointers no key covers --------------------------------------------------------------

POINTER_TABLES = ("nearest_systems", "sector_paths")
"""tuple: The tables whose rows name an object by `object_table` and `object_id`, which no foreign key covers (the target table varies)."""


def _orphans_pieces(ctx):
    keys = ctx.rows(
        "SELECT constraint_name AS c, table_name AS t, column_name AS col,"
        " referenced_table_name AS rt, referenced_column_name AS rcol"
        " FROM information_schema.key_column_usage"
        " WHERE table_schema = DATABASE() AND referenced_table_name IS NOT NULL"
        " ORDER BY constraint_name, ordinal_position")
    grouped = {}
    for key in keys:
        grouped.setdefault((key["c"], key["t"], key["rt"]), []).append([key["col"], key["rcol"]])
    return [{"key": [name, table, parent, pairs]} for (name, table, parent), pairs in grouped.items()] \
        + [{"pointer": table} for table in POINTER_TABLES]


def _orphans_piece(ctx, piece):
    """Child rows of one foreign key whose parent row is gone, or the rows of one pointer table that name an object that is gone."""
    if "pointer" in piece:
        return _dicts(_pointer_orphans(ctx, piece["pointer"], "object_table", "object_id"))
    name, table, parent, pairs = piece["key"]
    join = " AND ".join(f"c.`{a}` = p.`{b}`" for a, b in pairs)
    present = " AND ".join(f"c.`{a}` IS NOT NULL" for a, _ in pairs)
    orphans = ctx.rows(
        f"SELECT c.`{pairs[0][0]}` AS v FROM `{table}` c LEFT JOIN `{parent}` p ON {join}"
        f" WHERE {present} AND p.`{pairs[0][1]}` IS NULL LIMIT {MAX_ROWS_SHOWN + 1}")
    if not orphans:
        return []
    total = ctx.rows(
        f"SELECT COUNT(*) AS n FROM `{table}` c LEFT JOIN `{parent}` p ON {join}"
        f" WHERE {present} AND p.`{pairs[0][1]}` IS NULL")[0]["n"]
    shown = [f"{table}.{pairs[0][0]} = {row['v']} (no row in {parent})" for row in orphans]
    return _dicts([Problem(f"{total} row(s) of {table} have no parent in {parent} ({name}).", shown + ["..."] * (total > len(shown)))])


def _pointer_orphans(ctx, table, table_column, id_column):
    """Rows of `table` naming an object `(table_column, id_column)` that
    does not exist. No foreign key covers it, as the target table varies."""
    if not _table_exists(ctx, table):
        return []
    problems = []
    tables = {t.lower() for t in ctx.tables()}
    for kind in [row["k"] for row in ctx.rows(f"SELECT DISTINCT `{table_column}` AS k FROM `{table}`")]:
        if kind not in tables:
            problems.append(Problem(f"{table} points into {kind!r}, which is not a table."))
            continue
        rows = ctx.rows(
            f"SELECT c.`{id_column}` AS v FROM `{table}` c LEFT JOIN `{kind}` p ON p.id = c.`{id_column}`"
            f" WHERE c.`{table_column}` = %s AND p.id IS NULL LIMIT {MAX_ROWS_SHOWN + 1}", (kind,))
        if rows:
            problems.append(Problem(f"{table} has rows for {kind} objects that are gone.",
                                    [f"{kind} id {row['v']}" for row in rows]))
    return problems


# --- ids: one table a piece --------------------------------------------------------------------------------------------------------

def _ids_pieces(ctx):
    return sorted(store.ID_BLOCK_TABLES)


def _ids_piece(ctx, table):
    """`id_blocks.next_id` is above the largest id of one table it serves."""
    top = ctx.rows(f"SELECT MAX(id) AS m FROM `{table}`")[0]["m"]
    if top is None:
        return []
    rows = ctx.rows("SELECT next_id AS n FROM id_blocks WHERE table_name = %s", (table,))
    if not rows:
        return _dicts([Problem(f"{table} has rows (largest id {top}) but no id_blocks entry.")])
    if rows[0]["n"] <= top:
        return _dicts([Problem(f"id_blocks.next_id for {table} is {rows[0]['n']}, not above its largest id {top}; the next row would clash.")])
    return []


# --- values: one table a piece -----------------------------------------------------------------------------------------------------

def _bad_value_sql(column, low, high, low_in, high_in):
    """The SQL condition under which a value of `column` (the row aliased `t`) breaks its `VALUE_RULES` rule."""
    bad = []
    if low is not None:
        bad.append(f"t.`{column}` {'<' if low_in else '<='} {low}")
    if high is not None:
        bad.append(f"t.`{column}` {'>' if high_in else '>='} {high}")
    return f"t.`{column}` IS NOT NULL AND ({' OR '.join(bad)})"


def _values_pieces(ctx):
    return list(dict.fromkeys(rule[0] for rule in VALUE_RULES))


def _values_piece(ctx, table):
    """Columns of one table holding numbers no generator can write. One scan of the table counts every rule on it (a conditional sum per rule,
    PERF.258); only a rule that finds rows is read again, for the rows it lists."""
    rules = [rule for rule in VALUE_RULES if rule[0] == table]
    scope_sql, params = ctx.scope.where(_SECTOR_OF[table])
    sums = ", ".join(f"COALESCE(SUM(CASE WHEN {_bad_value_sql(*rule[1:])} THEN 1 ELSE 0 END), 0) AS n{index}" for index, rule in enumerate(rules))
    counts = ctx.rows(f"SELECT {sums} FROM `{table}` t WHERE 1 = 1{scope_sql}", params)[0]
    problems = []
    for index, rule in enumerate(rules):
        column, total = rule[1], int(counts[f"n{index}"])
        if total:
            where = f"{_bad_value_sql(*rule[1:])}{scope_sql}"
            rows = ctx.rows(f"SELECT t.id AS id, t.`{column}` AS v FROM `{table}` t WHERE {where} LIMIT {MAX_ROWS_SHOWN}", params)
            problems.append(Problem(
                f"{total} row(s) of {table} have {column} out of range.",
                [f"{table} id {row['id']}: {column} = {row['v']}" for row in rows] + ["..."] * (total > len(rows))))
    return _dicts(problems)


# --- counts: the sectors' stats, then the galaxy's stored counts -------------------------------------------------------------------

def _counts_pieces(ctx):
    return ["sectors", "galaxy"] if ctx.scope.everything else ["sectors"]


def _counts_piece(ctx, piece):
    """A sector's stats count the systems it holds (`sectors`), and the galaxy's stored counts (DB.121) the rows the tables hold (`galaxy`); a warning only."""
    problems = []
    if piece == "sectors":
        scope_sql, params = ctx.scope.where("s.id")
        rows = ctx.rows(
            "SELECT s.shell_index AS r, s.band_index AS l, s.wedge_index AS p, st.actual_systems AS stated,"
            " (SELECT COUNT(*) FROM star_systems y WHERE y.sector_id = s.id) AS actual"
            " FROM sectors s JOIN sector_stats st ON st.shell_index = s.shell_index"
            " AND st.band_index = s.band_index AND st.wedge_index = s.wedge_index"
            f" WHERE st.actual_systems IS NOT NULL{scope_sql}"
            " HAVING stated <> actual ORDER BY r, l, p", params)
        if rows:
            problems.append(Problem(f"{len(rows)} sector(s) hold a different number of systems than their stats say.",
                                    [f"sector {row['r']}/{row['l']}/{row['p']}: stats {row['stated']}, found {row['actual']}" for row in rows], warning=True))
    else:
        # DB.121: the running totals and the sectors' facet rows against a full count of the tables (`planetgen fill-block-totals` repairs).
        from planetgen.db import galaxy_counts
        found = galaxy_counts.problems(ctx.conn)
        if found:
            problems.append(Problem("The stored galaxy counts differ from a full count.", found, warning=True))
    return _dicts(problems)


# --- object index: one kind of object a piece --------------------------------------------------------------------------------------

def _object_index_pieces(ctx):
    from planetgen.db import objectindex
    return list(objectindex.KINDS) if _table_exists(ctx, "object_index") else ["no table"]


def _object_index_piece(ctx, kind):
    """The object index (DB.79, PERF.140) holds a row for one kind of object for every object, no row for a gone object and what each object
    holds. An index nothing has built yet is a warning, not damage."""
    from planetgen.db import objectindex
    if kind == "no table":
        return _dicts([Problem("The object index table is not there: run the update.", warning=True)])
    found = objectindex.check(ctx.conn, [kind])[kind]
    if not found["indexed"] and found["live"]:
        return _dicts([Problem(f"The object index holds no {kind} yet ({found['live']} exist): run `planetgen warm object-index`.", warning=True)])
    if found["missing"] or found["extra"] or found["different"]:
        return _dicts([Problem(
            f"The object index differs from the {kind} table: {found['missing']} missing, {found['extra']} "
            f"extra, {found['different']} different (of {found['live']}). `planetgen warm object-index` rebuilds it.", warning=False)])
    return []


# --- version keys: a run of sector ids a piece ---------------------------------------------------------------------------------------

VERSION_KEY_SECTORS = 2000
"""int: Sectors a piece of `check-version-keys` reads. Defined by the program, not a physical value: a piece is one read, and the bar counts them."""


def _version_keys_pieces(ctx):
    scope_sql, params = ctx.scope.where("s.id")
    ids = [row["id"] for row in ctx.rows(f"SELECT s.id AS id FROM sectors s WHERE 1 = 1{scope_sql} ORDER BY s.id", params)]
    return [[ids[start], ids[min(start + VERSION_KEY_SECTORS, len(ids)) - 1]] for start in range(0, len(ids), VERSION_KEY_SECTORS)]


def _version_keys_piece(ctx, span):
    """Every sector of ids `span` (in scope) points at a generation run that carries a 22 hex digit version key (DB.16): the sectors that do not, as text."""
    scope_sql, params = ctx.scope.where("s.id")
    rows = ctx.rows("SELECT s.id AS id, r.version_key AS k FROM sectors s LEFT JOIN generation_runs r ON r.id = s.run_id"
                    f" WHERE s.id BETWEEN %s AND %s{scope_sql}", (span[0], span[1], *params))
    return [f"sector id {row['id']}: {row['k']!r}" for row in rows if not row["k"] or not _VERSION_KEY.match(str(row["k"]))]


def _version_keys_combine(outputs):
    bad = [text for output in outputs for text in output]
    return [Problem(f"{len(bad)} sector(s) have a missing or malformed version key.", bad)] if bad else []


# --- storage: one sector a piece -----------------------------------------------------------------------------------------------------

def _storage_pieces(ctx):
    scope_sql, params = ctx.scope.where("s.id")
    return [row["id"] for row in ctx.rows(f"SELECT s.id AS id FROM sectors s WHERE 1 = 1{scope_sql} ORDER BY s.id", params)]


def _storage_piece(ctx, sector_id):
    """
    A sector whose `parity_dirty` flag is clear still hashes to its stored checksum (DB.45, DB.8): the rows are as they were when
    the parity record was made. Gives `["damaged", id]` for a mismatch, `["unprotected", id]` for a sector whose flag is set or that was never sealed (a writer changed it since,
    or nothing has sealed it yet: not checked), else `None`.
    """
    row = ctx.rows("SELECT s.parity_dirty AS dirty, s.storage_checksum AS sum FROM sectors s WHERE s.id = %s", (sector_id,))
    if not row:
        return None
    if row[0]["dirty"] or row[0]["sum"] is None:
        return ["unprotected", sector_id]
    return ["damaged", sector_id] if storage_checksum.check_sector(ctx.conn, sector_id) == storage_checksum.DAMAGED else None


def _storage_combine(outputs):
    damaged = [f"sector id {found[1]}" for found in outputs if found and found[0] == "damaged"]
    unprotected = [found[1] for found in outputs if found and found[0] == "unprotected"]
    problems = []
    if damaged:
        problems.append(Problem(f"{len(damaged)} sector(s) no longer hash to their stored checksum.", damaged))
    if unprotected:
        problems.append(Problem(f"{len(unprotected)} sector(s) have no current storage checksum (a write since, or never "
                                "sealed; not damage).", [], warning=True))
    return problems


DEEP_STATS_KIND = "db-check"
"""str: The `generation_stats` kind the per-system pass records its speed under (a unit is a system)."""

DEEP_FALLBACK_SECONDS_PER_SYSTEM = 0.04
"""float: Seconds a system takes to load and validate when the server has no speed recorded: a measured
figure on a mid-size server with a margin, so a first estimate errs long."""


@dataclass
class DeepEstimate:
    """What the `--deep` per-system pass will take: `systems` to validate, `seconds`, and whether `measured`
    (from this server's recorded speed) or a fallback guess."""
    systems: int
    seconds: float
    measured: bool

    def summary(self):
        from planetgen.generation import stats as generation_stats

        if self.measured:
            return (f"Deep check: {self.systems:,} star systems to validate, about "
                    f"{generation_stats.format_duration(self.seconds)} (from the speed this server recorded).")
        return (f"Deep check: {self.systems:,} star systems to validate, about "
                f"{generation_stats.time_text(self.seconds, False)}.")

    def as_dict(self):
        return {"deep_check": True, "systems": self.systems, "seconds": self.seconds, "measured": self.measured,
                "what": "the deep database check", "summary": self.summary()}


def _system_filter(scope):
    """`(sql, params)`: the condition that picks the star systems the `--deep` pass validates."""
    where, params = scope.where("sector_id")
    if scope.wide_binaries:
        where += " AND binary_configuration = 'wide'"
    return where, params


def _system_ids(conn, scope):
    where, params = _system_filter(scope)
    return [row["id"] for row in conn.execute(
        f"SELECT id FROM star_systems WHERE 1 = 1{where} ORDER BY id", params).fetchall()]


def deep_estimate(conn, config, scope=None):
    """The `DeepEstimate` for the systems check over every star system in `scope` (all by default), from the recorded
    speed of the pass (`DEEP_STATS_KIND`), else `DEEP_FALLBACK_SECONDS_PER_SYSTEM` a system."""
    from planetgen.generation import run_common

    where, params = _system_filter(scope or Scope())
    count = conn.execute(f"SELECT COUNT(*) AS n FROM star_systems WHERE 1 = 1{where}", params).fetchone()["n"]
    stats = run_common._stats_for_config(config)
    seconds = steps.predict(stats, DEEP_STATS_KIND, count)
    if seconds is None:
        return DeepEstimate(count, count * DEEP_FALLBACK_SECONDS_PER_SYSTEM, False)
    return DeepEstimate(count, seconds, True)


def _systems_pieces(ctx):
    return _system_ids(ctx.conn, ctx.scope)


def _systems_piece(ctx, system_id):
    """One star system in scope loads and passes `validation.check_star_system`. Gives `{"seconds", "failing"}`: the seconds it took (the speed of the pass is recorded
    from them, `DEEP_STATS_KIND`) and, for a system that fails, the text that says what is wrong (a system that cannot be loaded says so), else `None`."""
    started = time.monotonic()
    try:
        system = store.load_star_system(ctx.conn, system_id)
        found = validation.check_star_system(system)
        lines = [f"{problem.body}: {problem.message}" for problem in found]
        name = system.name
    except Exception as exc:  # noqa: BLE001 -- a system that will not load is the finding
        log.swallowed("db/check", "_systems_piece", exc)
        lines, name = [f"could not be loaded or checked: {type(exc).__name__}: {exc}"], "?"
    failing = f"star system {name!r} (id {system_id}): " + "; ".join(lines) if lines else None
    return {"seconds": time.monotonic() - started, "failing": failing}


def _systems_combine(outputs):
    failing = [output["failing"] for output in outputs if output["failing"]]
    return [Problem(f"{len(failing)} star system(s) fail validation.", failing)] if failing else []


def _record_systems_speed(config, outputs):
    """Records how long the pass took a system (`DEEP_STATS_KIND`), so the next estimate (`deep_estimate`) knows this server's speed; never fails a check."""
    from planetgen.generation import run_common
    seconds = sum(output["seconds"] for output in outputs)
    if not outputs or seconds <= 0:
        return
    try:
        stats = run_common._stats_for_config(config)
        stats.record(DEEP_STATS_KIND, 0.0, seconds, systems=len(outputs), stars=len(outputs), workers=1)
        stats.flush()
    except Exception as exc:  # noqa: BLE001 -- statistics never fail a check
        log.debug(f"Database check: could not record the speed of the systems pass: {exc}")


@dataclass(frozen=True)
class Check:
    """One check: its `name` (the word `planetgen check-db NAME` takes, and the tail of its command `check-NAME`), the `label` its result line
    carries, and how the engine runs it (OPS.204): `pieces(ctx)` lists the pieces of the one unit (JSON-able: a table, a key, a sector, a system),
    `piece(ctx, piece)` does one and gives what it found (JSON-able), and `combine(outputs)` makes the problems of the pieces' outputs
    (`_problems_of`, the pieces' problems one after another, unless the check says another). `scoped` marks the checks that take the
    sector list, `wide_binaries` the one that takes that."""
    name: str
    label: str
    pieces: object
    piece: object
    combine: object = _problems_of
    scoped: bool = False
    wide_binaries: bool = False


CHECKS = (
    Check("revision", "revision", _revision_pieces, _revision_piece),
    Check("tables", "tables", _tables_pieces, _tables_piece),
    Check("tables-extended", "tables (extended)", _tables_pieces, _tables_extended_piece),
    Check("orphans", "orphans", _orphans_pieces, _orphans_piece),
    Check("ids", "ids", _ids_pieces, _ids_piece),
    Check("values", "values", _values_pieces, _values_piece, scoped=True),
    Check("counts", "counts", _counts_pieces, _counts_piece, scoped=True),
    Check("version-keys", "version keys", _version_keys_pieces, _version_keys_piece, _version_keys_combine, scoped=True),
    Check("storage", "storage checksums", _storage_pieces, _storage_piece, _storage_combine, scoped=True),
    Check("systems", "systems", _systems_pieces, _systems_piece, _systems_combine, scoped=True, wide_binaries=True),
    Check("object-index", "object index", _object_index_pieces, _object_index_piece),
)
"""tuple: Every check, in the order a list of them runs and reports."""

BY_NAME = {check.name: check for check in CHECKS}

_CHECK_ERRORS = (sqlalchemy.exc.SQLAlchemyError, store.pymysql.err.MySQLError, KeyError)
"""tuple: What a check that fails to run raises: "could not check", which is not the same as damage."""


def scope_of(args):
    """The `Scope` a check's unit arguments name (`shells`, `addresses`, `wide_binaries`)."""
    return Scope(tuple(args.get("shells") or ()), tuple(tuple(a) for a in args.get("addresses") or ()),
                 bool(args.get("wide_binaries")))


def _config_of(args):
    return store.MySQLConfig.from_args(args["database"])


def _read_only(args, work):
    """`work(_Context)` on a connection of the database the run named (`wrapper.database_args`; the password is this worker's settings'), opened as the checks always did: without
    `ensure_schema`, so a check never mends the schema it is looking at, and closed after. A check writes nothing, so it has nothing to commit."""
    config = _config_of(args)
    conn = store.get_connection(config, ensure_schema=False)
    try:
        return work(_Context(conn, config, scope_of(args)))
    finally:
        conn.close()


def _pieces_of(name, args, cursor):
    """The pieces of the unit of the command `check-<name>` after the first `cursor` of them. A check that cannot list its pieces gives one piece that says so."""
    def list_them(ctx):
        try:
            return list(BY_NAME[name].pieces(ctx))
        except _CHECK_ERRORS as exc:
            return [{"could_not_check": f"{type(exc).__name__}: {exc}"}]
    return _read_only(args, list_them)[cursor:]


def check_piece(args, piece, ctx):
    """
    One piece of a check (`do` of every `check-NAME` command; the command is `ctx.kind`): runs it on a connection of its own and gives `{"found": what it found}`, or
    `{"could_not_check": the error}` when the database refused a query (a check that cannot run is "could not check", never an exception and never damage).
    """
    if isinstance(piece, dict) and "could_not_check" in piece:
        return piece
    check = BY_NAME[ctx.kind[len("check-"):]]
    try:
        return {"found": _read_only(args, lambda context: check.piece(context, piece))}
    except _CHECK_ERRORS as exc:
        return {"could_not_check": f"{type(exc).__name__}: {exc}"}


def _result_of(name, args, results):
    """The unit's result: the `CheckResult` of the check `name` as JSON, from what its pieces found."""
    check = BY_NAME[name]
    for result in results:
        if "could_not_check" in result:
            return CheckResult(check.label, ERROR, error=result["could_not_check"]).as_dict()
    outputs = [result["found"] for result in results]
    if name == "systems":
        _record_systems_speed(_config_of(args), outputs)
    problems = check.combine(outputs)
    if any(not p.warning for p in problems):
        status = FAIL
    else:
        status = WARN if problems else PASS
    return CheckResult(check.label, status, problems).as_dict()


def check_revision_pieces(args, cursor):
    """The pieces of `check-revision`: the records of the version, then the tables against the models (PERF.126)."""
    return _pieces_of("revision", args, cursor)


def check_revision_gather(args, results):
    """The result of `check-revision`."""
    return _result_of("revision", args, results)


def check_tables_pieces(args, cursor):
    """The pieces of `check-tables`: one table each, `CHECK TABLE ... QUICK`."""
    return _pieces_of("tables", args, cursor)


def check_tables_gather(args, results):
    """The result of `check-tables`."""
    return _result_of("tables", args, results)


def check_tables_extended_pieces(args, cursor):
    """The pieces of `check-tables-extended`: one table each, `CHECK TABLE ... EXTENDED`."""
    return _pieces_of("tables-extended", args, cursor)


def check_tables_extended_gather(args, results):
    """The result of `check-tables-extended`."""
    return _result_of("tables-extended", args, results)


def check_orphans_pieces(args, cursor):
    """The pieces of `check-orphans`: one foreign key each, then the pointer tables no key covers."""
    return _pieces_of("orphans", args, cursor)


def check_orphans_gather(args, results):
    """The result of `check-orphans`: child rows whose parent is gone."""
    return _result_of("orphans", args, results)


def check_ids_pieces(args, cursor):
    """The pieces of `check-ids`: one table each of those `id_blocks` serves."""
    return _pieces_of("ids", args, cursor)


def check_ids_gather(args, results):
    """The result of `check-ids`: `id_blocks` is ahead of every id."""
    return _result_of("ids", args, results)


def check_values_pieces(args, cursor):
    """The pieces of `check-values`: one table each of those with a value rule, in the scope the arguments name."""
    return _pieces_of("values", args, cursor)


def check_values_gather(args, results):
    """The result of `check-values`: columns holding numbers no generator can write."""
    return _result_of("values", args, results)


def check_counts_pieces(args, cursor):
    """The pieces of `check-counts`: the sectors' stats, then (with no scope) the galaxy's stored counts."""
    return _pieces_of("counts", args, cursor)


def check_counts_gather(args, results):
    """The result of `check-counts`: a sector's stats count its systems."""
    return _result_of("counts", args, results)


def check_version_keys_pieces(args, cursor):
    """The pieces of `check-version-keys`: a run of `VERSION_KEY_SECTORS` sectors each, in the scope the arguments name."""
    return _pieces_of("version-keys", args, cursor)


def check_version_keys_gather(args, results):
    """The result of `check-version-keys`: every sector in scope has a well formed version key."""
    return _result_of("version-keys", args, results)


def check_storage_pieces(args, cursor):
    """The pieces of `check-storage`: one sector each, in the scope the arguments name."""
    return _pieces_of("storage", args, cursor)


def check_storage_gather(args, results):
    """The result of `check-storage`: every sealed sector in scope still hashes to its storage checksum (DB.45)."""
    return _result_of("storage", args, results)


def check_systems_pieces(args, cursor):
    """The pieces of `check-systems`: one star system each, in the scope the arguments name (the wide binaries only with `wide_binaries`)."""
    return _pieces_of("systems", args, cursor)


def check_systems_gather(args, results):
    """The result of `check-systems`: every star system in scope validates; its speed is recorded for the next estimate."""
    return _result_of("systems", args, results)


def check_object_index_pieces(args, cursor):
    """The pieces of `check-object-index`: one kind of object each (PERF.140)."""
    return _pieces_of("object-index", args, cursor)


def check_object_index_gather(args, results):
    """The result of `check-object-index`: the object index matches the tables."""
    return _result_of("object-index", args, results)


def expand_check(params):
    """The one unit of a check's stage: the arguments the step gave, unchanged."""
    from planetgen.queue import commands
    return commands.Expansion([(dict(params), 1.0)], {})


def expand_check_systems(params):
    """The one unit of the systems check, its cost the number of systems it will validate (the bar and the ETA)."""
    from planetgen.queue import commands
    config = _config_of(params)
    conn = store.get_connection(config, ensure_schema=False)
    try:
        where, values = _system_filter(scope_of(params))
        count = conn.execute(f"SELECT COUNT(*) AS n FROM star_systems WHERE 1 = 1{where}", values).fetchone()["n"]
    finally:
        conn.close()
    return commands.Expansion([(dict(params), float(max(count, 1)))], {"systems": count})
