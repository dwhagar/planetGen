# planetgen/db/check.py

"""
The database check (DB.8)
=========================

`planetgen check-db` (and the Admin dashboard's "Check the database"
button, which runs it as a job) looks for damage in a galaxy database and
changes nothing. See `docs/design/db-check-and-parity-repair.md` section 6.

Each check is a function of a `_Context` that returns a list of `Problem`s
(empty when the check passes); `run_checks` runs them all, turns a check
that itself failed to run into an "error" result (could not check, which is
not the same as damaged), and returns a `Report`.

The checks:

1. `revision`: the Alembic revision, `schema_migrations` and the code's
   `store.SCHEMA_VERSION` agree, and the tables match `db/models.py`.
2. `tables`: `CHECK TABLE` on every table (`--no-check-table` skips it,
   `--deep` makes it EXTENDED).
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
8. `systems` (`--deep` only, DB.21): every star system is loaded and run
   through `validation.check_star_system` (orbit spacing, the stable-orbit
   ceiling, each planet's moons). Slow, so `deep_estimate` tells how long
   first, from the speed this server recorded for it, and the pass draws
   its own progress bar.

`--ring` and `--sector` scope the value, count and version-key checks to
those sectors; the others always look at the whole database. A galaxy too
big to check in one go is checked a region at a time.
"""

import re
import warnings
from dataclasses import dataclass, field

import sqlalchemy
from alembic.autogenerate import compare_metadata
from alembic.migration import MigrationContext

from planetgen.db import alembic_runner, models, store
from planetgen.generation import steps, validation

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
    ("planets", "orbital_phase_deg", 0, 360, True, True),
    ("planets", "orbital_inclination_deg", 0, 180, True, True),
    ("moons", "radius_km", 0, None, False, None),
    ("moons", "mass_kg", 0, None, False, None),
    ("moons", "period_years", 0, None, False, None),
    ("moons", "surface_temperature_k", 0, None, False, None),
    ("moons", "orbital_phase_deg", 0, 360, True, True),
    ("moons", "orbital_inclination_deg", 0, 180, True, True),
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
    """Which sectors the sector-scoped checks look at: whole `rings` and
    single `addresses` (`(ring, layer, slot)`); neither means all."""
    rings: tuple = ()
    addresses: tuple = ()

    @property
    def everything(self):
        return not self.rings and not self.addresses

    def where(self, sector_sql):
        """Returns `(sql, params)`: a condition on `sector_sql` (empty for
        everything)."""
        if self.everything:
            return "", ()
        clauses, params = [], []
        if self.rings:
            clauses.append("ring_index IN (%s)" % ",".join(["%s"] * len(self.rings)))
            params.extend(self.rings)
        for address in self.addresses:
            clauses.append("(ring_index = %s AND layer_index = %s AND ring_slot_index = %s)")
            params.extend(address)
        return f" AND {sector_sql} IN (SELECT id FROM sectors WHERE {' OR '.join(clauses)})", tuple(params)


@dataclass
class _Context:
    conn: object
    config: object
    scope: Scope
    deep: bool
    check_table: bool

    def rows(self, sql, params=()):
        return self.conn.execute(sql, params).fetchall()

    def tables(self):
        return [row["t"] for row in self.rows(
            "SELECT table_name AS t FROM information_schema.tables"
            " WHERE table_schema = DATABASE() AND table_type = 'BASE TABLE' ORDER BY 1")]


def _table_exists(ctx, name):
    return name in {t.lower() for t in ctx.tables()}


def check_revision(ctx):
    """The three records of the schema version agree with the code, and the
    tables match the models."""
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
    detected = store.detect_schema_version(ctx.conn)
    if detected != head:
        problems.append(Problem(f"The tables look like schema version {detected}, not {head}."))
    drift = _model_drift(ctx.config)
    if drift:
        problems.append(Problem("The tables differ from db/models.py:", drift))
    return problems


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


def check_tables(ctx):
    """`CHECK TABLE` on every table (skipped with `--no-check-table`)."""
    if not ctx.check_table:
        return []
    mode = "EXTENDED" if ctx.deep else "QUICK"
    problems = []
    for table in ctx.tables():
        for row in ctx.rows(f"CHECK TABLE `{table}` {mode}"):
            kind = str(row.get("Msg_type", "")).lower()
            if kind in ("error", "corrupt") or (kind == "status" and str(row.get("Msg_text")).upper() != "OK"):
                problems.append(Problem(f"CHECK TABLE {table}: {row.get('Msg_text')}"))
    return problems


def check_orphans(ctx):
    """Child rows whose parent row is gone."""
    problems = []
    keys = ctx.rows(
        "SELECT constraint_name AS c, table_name AS t, column_name AS col,"
        " referenced_table_name AS rt, referenced_column_name AS rcol"
        " FROM information_schema.key_column_usage"
        " WHERE table_schema = DATABASE() AND referenced_table_name IS NOT NULL"
        " ORDER BY constraint_name, ordinal_position")
    grouped = {}
    for key in keys:
        grouped.setdefault((key["c"], key["t"], key["rt"]), []).append((key["col"], key["rcol"]))
    for (name, table, parent), pairs in grouped.items():
        join = " AND ".join(f"c.`{a}` = p.`{b}`" for a, b in pairs)
        present = " AND ".join(f"c.`{a}` IS NOT NULL" for a, _ in pairs)
        orphans = ctx.rows(
            f"SELECT c.`{pairs[0][0]}` AS v FROM `{table}` c LEFT JOIN `{parent}` p ON {join}"
            f" WHERE {present} AND p.`{pairs[0][1]}` IS NULL LIMIT {MAX_ROWS_SHOWN + 1}")
        if orphans:
            total = ctx.rows(
                f"SELECT COUNT(*) AS n FROM `{table}` c LEFT JOIN `{parent}` p ON {join}"
                f" WHERE {present} AND p.`{pairs[0][1]}` IS NULL")[0]["n"]
            shown = [f"{table}.{pairs[0][0]} = {row['v']} (no row in {parent})" for row in orphans]
            problems.append(Problem(f"{total} row(s) of {table} have no parent in {parent} ({name}).",
                                    shown + ["..."] * (total > len(shown))))
    problems.extend(_pointer_orphans(ctx, "nearest_systems", "object_table", "object_id"))
    problems.extend(_pointer_orphans(ctx, "sector_paths", "object_table", "object_id"))
    return problems


def _pointer_orphans(ctx, table, table_column, id_column):
    """Rows of `table` naming an object `(table_column, id_column)` that
    does not exist. No foreign key covers it, as the target table varies."""
    if not _table_exists(ctx, table):
        return []
    problems = []
    for kind in [row["k"] for row in ctx.rows(f"SELECT DISTINCT `{table_column}` AS k FROM `{table}`")]:
        if kind not in {t.lower() for t in ctx.tables()}:
            problems.append(Problem(f"{table} points into {kind!r}, which is not a table."))
            continue
        rows = ctx.rows(
            f"SELECT c.`{id_column}` AS v FROM `{table}` c LEFT JOIN `{kind}` p ON p.id = c.`{id_column}`"
            f" WHERE c.`{table_column}` = %s AND p.id IS NULL LIMIT {MAX_ROWS_SHOWN + 1}", (kind,))
        if rows:
            problems.append(Problem(f"{table} has rows for {kind} objects that are gone.",
                                    [f"{kind} id {row['v']}" for row in rows]))
    return problems


def check_ids(ctx):
    """`id_blocks.next_id` is above the largest id of each table it serves."""
    problems = []
    blocks = {row["t"]: row["n"] for row in ctx.rows("SELECT table_name AS t, next_id AS n FROM id_blocks")}
    for table in sorted(store.ID_BLOCK_TABLES):
        top = ctx.rows(f"SELECT MAX(id) AS m FROM `{table}`")[0]["m"]
        if top is None:
            continue
        if table not in blocks:
            problems.append(Problem(f"{table} has rows (largest id {top}) but no id_blocks entry."))
        elif blocks[table] <= top:
            problems.append(Problem(f"id_blocks.next_id for {table} is {blocks[table]}, "
                                    f"not above its largest id {top}; the next row would clash."))
    return problems


def check_values(ctx):
    """Columns holding numbers no generator can write."""
    problems = []
    for table, column, low, high, low_in, high_in in VALUE_RULES:
        bad = []
        if low is not None:
            bad.append(f"t.`{column}` {'<' if low_in else '<='} {low}")
        if high is not None:
            bad.append(f"t.`{column}` {'>' if high_in else '>='} {high}")
        scope_sql, params = ctx.scope.where(_SECTOR_OF[table])
        where = f"t.`{column}` IS NOT NULL AND ({' OR '.join(bad)}){scope_sql}"
        total = ctx.rows(f"SELECT COUNT(*) AS n FROM `{table}` t WHERE {where}", params)[0]["n"]
        if total:
            rows = ctx.rows(f"SELECT t.id AS id, t.`{column}` AS v FROM `{table}` t WHERE {where}"
                            f" LIMIT {MAX_ROWS_SHOWN}", params)
            problems.append(Problem(
                f"{total} row(s) of {table} have {column} out of range.",
                [f"{table} id {row['id']}: {column} = {row['v']}" for row in rows] + ["..."] * (total > len(rows))))
    return problems


def check_counts(ctx):
    """A sector's stats count the systems it holds (a warning only)."""
    scope_sql, params = ctx.scope.where("s.id")
    rows = ctx.rows(
        "SELECT s.ring_index AS r, s.layer_index AS l, s.ring_slot_index AS p, st.actual_systems AS stated,"
        " (SELECT COUNT(*) FROM star_systems y WHERE y.sector_id = s.id) AS actual"
        " FROM sectors s JOIN sector_stats st ON st.ring_index = s.ring_index"
        " AND st.layer_index = s.layer_index AND st.ring_slot_index = s.ring_slot_index"
        f" WHERE st.actual_systems IS NOT NULL{scope_sql}"
        " HAVING stated <> actual ORDER BY r, l, p", params)
    if not rows:
        return []
    return [Problem(f"{len(rows)} sector(s) hold a different number of systems than their stats say.",
                    [f"sector {row['r']}/{row['l']}/{row['p']}: stats {row['stated']}, found {row['actual']}"
                     for row in rows], warning=True)]


def check_version_keys(ctx):
    """Every sector carries a 22 hex digit version key."""
    scope_sql, params = ctx.scope.where("s.id")
    rows = ctx.rows(f"SELECT s.id AS id, s.version_key AS k FROM sectors s WHERE 1 = 1{scope_sql}", params)
    bad = [row for row in rows if not row["k"] or not _VERSION_KEY.match(str(row["k"]))]
    if not bad:
        return []
    return [Problem(f"{len(bad)} sector(s) have a missing or malformed version key.",
                    [f"sector id {row['id']}: {row['k']!r}" for row in bad])]


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

        basis = ("from the speed this server recorded" if self.measured
                 else "a rough guess: this server has not recorded the speed yet")
        return (f"Deep check: {self.systems:,} star systems to validate, about "
                f"{generation_stats.format_duration(self.seconds)} ({basis}), plus the table checks.")

    def as_dict(self):
        return {"deep_check": True, "systems": self.systems, "seconds": self.seconds, "measured": self.measured,
                "what": "the deep database check", "summary": self.summary()}


def _system_ids(conn, scope):
    where, params = scope.where("sector_id")
    return [row["id"] for row in conn.execute(
        f"SELECT id FROM star_systems WHERE 1 = 1{where} ORDER BY id", params).fetchall()]


def deep_estimate(conn, config, scope=None):
    """The `DeepEstimate` for validating every star system in `scope` (all by default), from the recorded
    speed of the pass (`DEEP_STATS_KIND`), else `DEEP_FALLBACK_SECONDS_PER_SYSTEM` a system."""
    from planetgen.generation import run_common

    where, params = (scope or Scope()).where("sector_id")
    count = conn.execute(f"SELECT COUNT(*) AS n FROM star_systems WHERE 1 = 1{where}", params).fetchone()["n"]
    stats = run_common._stats_for_config(config)
    seconds = steps.predict(stats, DEEP_STATS_KIND, count)
    if seconds is None:
        return DeepEstimate(count, count * DEEP_FALLBACK_SECONDS_PER_SYSTEM, False)
    return DeepEstimate(count, seconds, True)


def check_systems(ctx):
    """Every star system in scope loads and passes `validation.check_star_system`. One problem lists each
    failing system with what is wrong (a system that cannot be loaded says so). Only run with `--deep`."""
    from planetgen.generation import run_common

    ids = _system_ids(ctx.conn, ctx.scope)
    failing = []
    stats = run_common._stats_for_config(ctx.config)
    with steps.step("Validating star systems", DEEP_STATS_KIND, len(ids), stats=stats, own=True) as bar:
        for system_id in ids:
            try:
                system = store.load_star_system(ctx.conn, system_id)
                found = validation.check_star_system(system)
                lines = [f"{problem.body}: {problem.message}" for problem in found]
                name = system.name
            except Exception as exc:  # noqa: BLE001 -- a system that will not load is the finding
                lines, name = [f"could not be loaded or checked: {type(exc).__name__}: {exc}"], "?"
            if lines:
                failing.append(f"star system {name!r} (id {system_id}): " + "; ".join(lines))
            bar.advance()
    if not failing:
        return []
    return [Problem(f"{len(failing)} star system(s) fail validation.", failing)]


CHECKS = (
    ("revision", check_revision),
    ("tables", check_tables),
    ("orphans", check_orphans),
    ("ids", check_ids),
    ("values", check_values),
    ("counts", check_counts),
    ("version keys", check_version_keys),
)
"""tuple: `(name, function)` in the order they run and report."""


def run_checks(conn, config, scope=None, deep=False, check_table=True, on_progress=None):
    """
    Runs every check and returns the `Report`. Nothing is written. A check
    that raises is reported as "could not check" with its error.

    Args:
        conn (Connection): An open connection (no schema changes needed).
        config (MySQLConfig): The same database, for the model comparison.
        scope (Scope, optional): The sectors the sector-scoped checks use.
        deep (bool): EXTENDED table checks.
        check_table (bool): Run `CHECK TABLE` at all.
        on_progress (callable, optional): `on_progress(name)` before each.
    """
    ctx = _Context(conn, config, scope or Scope(), deep, check_table)
    results = []
    for name, function in CHECKS + ((("systems", check_systems),) if deep else ()):
        if on_progress is not None:
            on_progress(name)
        try:
            problems = function(ctx)
        except (sqlalchemy.exc.SQLAlchemyError, store.pymysql.err.MySQLError, KeyError) as exc:
            results.append(CheckResult(name, ERROR, error=f"{type(exc).__name__}: {exc}"))
            continue
        if any(not p.warning for p in problems):
            status = FAIL
        else:
            status = WARN if problems else PASS
        results.append(CheckResult(name, status, problems))
    return Report(results)
