# planetgen/db/fingerprint.py

"""
A galaxy's fingerprint (GEN.58)
===============================

A canonical SHA-256 digest of each sector's generated content, and one for
a region, so two builds can be compared in the sense
`docs/design/reproducible-galaxies.md` section 2 defines: every generated
object has the same address, position, properties and name.

A sector's content is its `sectors` row, every object filed in it (star
systems with their stars, planets, moons, belts, comets and facilities;
phenomena; stand-alone facilities), those objects' child rows
(compositions, paragraphs, spectra, nebula shape balls), and the
address-keyed rows of its cell (`bright_stars`, `phenomenon_scatter`).

Left out, as the design says: row ids and timestamps; each row's clock
(`epoch_unix`, `next_update_due`, GEN.106), which records when it was
last moved, not what was made; and what is rebuilt from positions
(`nearest_systems`, the `location` text written from it, sector paths,
sector stats). Population data (species, polities, owners) is rebuilt
later, so it is left out too, as are bookkeeping tables (run history, id
blocks, name registries).

A foreign key would carry a row id, so it is replaced by what it points
at: a sector by its address, an object by its unique ID (`uid`, GEN.170), a
system configuration by the digest of its content. Each row becomes one
line of canonical JSON (keys sorted, floats in Python's shortest
round-trip form, -0.0 as 0.0, bytes in hex), and a sector's digest is the
SHA-256 of its lines sorted, so the order rows were saved in doesn't
count. A region's digest is the SHA-256 of its sectors' lines, in address
order, with the plan's (the `galaxy_*` tables) first for the whole galaxy.

This is the galaxy as stored now. Replaying the settings file's edits and
regenerations onto a rebuild before comparing (GEN.59) comes later; the
same functions back GEN.57's order test, TEST.77 and OPS.12.
"""

import decimal
import hashlib
import json
import math
from dataclasses import dataclass

SECTOR_TABLES = ("star_systems", "black_holes", "neutron_stars", "nebulae", "supernova_remnants", "quasars",
                 "rogue_planets", "interstellar_comets", "asteroid_fields", "facilities")
"""tuple: Tables whose rows are filed in a sector by `sector_id`."""

SYSTEM_TABLES = ("stars", "planets", "moons", "asteroid_belts", "comets", "facilities")
"""tuple: Tables whose rows belong to a star system by `star_system_id`."""

STAR_TABLES = ("black_holes", "neutron_stars")
"""tuple: Tables whose rows can belong to a star by `star_id` (a compact
remnant in a system) with no `sector_id` of their own."""

CHILD_TABLES = (
    ("planet_evolutionary_paragraphs", "planet_id", "planets"),
    ("planet_reflection_spectrum", "planet_id", "planets"),
    ("moon_evolutionary_paragraphs", "moon_id", "moons"),
    ("moon_reflection_spectrum", "moon_id", "moons"),
    ("asteroid_belt_composition", "belt_id", "asteroid_belts"),
    ("comet_composition", "comet_id", "comets"),
    ("interstellar_comet_composition", "comet_id", "interstellar_comets"),
    ("asteroid_field_composition", "field_id", "asteroid_fields"),
    ("nebula_shape_balls", "nebula_id", "nebulae"),
)
"""tuple: `(table, foreign key, parent table)` of the child rows of an
object a sector holds."""

ADDRESS_TABLES = ("bright_stars", "phenomenon_scatter")
"""tuple: Tables keyed by a sector address (`ring_index`, `layer_index`,
`ring_slot_index`) rather than a sector row."""

PLAN_TABLES = ("galaxy_shape", "galaxy_layer", "galaxy_column")
"""tuple: The galaxy's plan, compared for the whole galaxy."""

CONTENT_TABLES = frozenset(SECTOR_TABLES + SYSTEM_TABLES + STAR_TABLES + ADDRESS_TABLES + PLAN_TABLES
                           + tuple(table for table, _fk, _parent in CHILD_TABLES) + ("sectors",))
"""frozenset: Every table a fingerprint reads."""

LEFT_OUT_TABLES = frozenset({
    "schema_migrations", "alembic_version", "orbit_simulation_state", "nearest_systems", "sector_paths",
    "sector_path_knots", "sector_stats", "sector_name_registry", "system_name_registry", "generation_runs",
    "generation_run_arguments", "id_blocks", "id_counters", "phenomenon_scatter_classes", "sector_system_counts", "species", "polities", "system_owners", "population_state",
    "system_configs", "system_config_slots",
})
"""frozenset: Tables a fingerprint doesn't read as content (see the module
docstring; a system configuration counts through the systems that use
it)."""

LEFT_OUT_COLUMNS = frozenset({"id", "epoch_unix", "next_update_due", "binary_epoch_unix", "binary_next_update_due"})
"""frozenset: Columns left out of every table: the row id and the clocks."""

LEFT_OUT_TABLE_COLUMNS = {"star_systems": frozenset({"location"})}
"""dict: Columns left out of one table: `location` is written from the
nearest systems."""

_TIME_TYPES = frozenset({"datetime", "timestamp", "date", "time"})

_IN_BATCH = 1000


def canonical_value(value):
    """`value` as a JSON-safe canonical form: floats in shortest round-trip
    form (with -0.0 as 0.0, and NaN and the infinities as strings),
    decimals as their normalized text, bytes in hex."""
    if isinstance(value, float):
        if math.isnan(value):
            return "nan"
        if math.isinf(value):
            return "inf" if value > 0 else "-inf"
        return 0.0 if value == 0 else value
    if isinstance(value, decimal.Decimal):
        return format(value.normalize(), "f")
    if isinstance(value, (bytes, bytearray, memoryview)):
        return "0x" + bytes(value).hex()
    return value


def _line(table, row):
    return json.dumps([table, row], sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def _digest(lines):
    return hashlib.sha256("\n".join(lines).encode("utf-8")).hexdigest()


@dataclass
class Fingerprint:
    """A region's fingerprint: `sectors` is `[(label, digest)]` in address
    order (a placed sector's label is `"ring layer slot"`, an unplaced
    one's `"unplaced <name>"`), `plan` the plan's digest (the whole galaxy
    only, else `None`), `region` the digest over both."""
    sectors: list
    plan: object
    region: str

    def lines(self):
        """The printed form, one line a sector, then the plan and region."""
        out = [f"{label} {digest}" for label, digest in self.sectors]
        if self.plan is not None:
            out.append(f"plan {self.plan}")
        out.append(f"region {self.region} ({len(self.sectors)} sector(s))")
        return out


class _Reader:
    """Reads and canonicalizes rows of one database, caching the schema's
    columns and foreign keys and the keys of the rows they point at."""

    def __init__(self, conn):
        self.conn = conn
        self._columns = {}
        self._foreign = None
        self._keys = {}

    def columns(self, table):
        if table not in self._columns:
            rows = self.conn.execute(
                "SELECT column_name AS name, data_type AS type FROM information_schema.columns"
                " WHERE table_schema = DATABASE() AND table_name = ? ORDER BY ordinal_position", (table,)).fetchall()
            skipped = LEFT_OUT_COLUMNS | LEFT_OUT_TABLE_COLUMNS.get(table, frozenset())
            self._columns[table] = [row["name"] for row in rows
                                    if row["name"] not in skipped and row["type"].lower() not in _TIME_TYPES]
        return self._columns[table]

    def foreign_keys(self, table):
        """`{column: referenced table}` of `table`."""
        if self._foreign is None:
            self._foreign = {}
            for row in self.conn.execute(
                    "SELECT table_name AS t, column_name AS c, referenced_table_name AS r"
                    " FROM information_schema.key_column_usage"
                    " WHERE table_schema = DATABASE() AND referenced_table_name IS NOT NULL").fetchall():
                self._foreign.setdefault(row["t"], {})[row["c"]] = row["r"]
        return self._foreign.get(table, {})

    def select(self, table, where, params=()):
        """The rows of `table` matching `where`: the compared columns, and
        `id` when the table has one (to tell rows apart, never compared)."""
        names = [f"`{c}`" for c in self.columns(table)]
        if "id" in self._all_columns(table):
            names.insert(0, "id")
        return self.conn.execute(f"SELECT {', '.join(names)} FROM `{table}` WHERE {where}", tuple(params)).fetchall()

    def _all_columns(self, table):
        key = ("all", table)
        if key not in self._columns:
            self._columns[key] = {row["name"] for row in self.conn.execute(
                "SELECT column_name AS name FROM information_schema.columns"
                " WHERE table_schema = DATABASE() AND table_name = ?", (table,)).fetchall()}
        return self._columns[key]

    def select_in(self, table, column, ids):
        ids = sorted(set(ids))
        rows = []
        for first in range(0, len(ids), _IN_BATCH):
            chunk = ids[first:first + _IN_BATCH]
            rows.extend(self.select(table, f"`{column}` IN ({', '.join('?' * len(chunk))})", chunk))
        return rows

    def line(self, table, row):
        """`row` of `table` as its canonical line, foreign keys replaced by
        the keys of what they point at."""
        foreign = self.foreign_keys(table)
        canon = {}
        for column in self.columns(table):
            value = row[column]
            if column in foreign:
                value = None if value is None else self.key(foreign[column], value)
            canon[column] = canonical_value(value)
        return _line(table, canon)

    def key(self, table, row_id):
        """A stable key for row `row_id` of `table`: never its id."""
        cache_key = (table, row_id)
        if cache_key in self._keys:
            return self._keys[cache_key]
        if table == "sectors":
            row = self.conn.execute("SELECT ring_index, layer_index, ring_slot_index, name FROM sectors WHERE id = ?",
                                    (row_id,)).fetchone()
            key = None if row is None else "sector:" + _sector_label(row)
        elif table == "system_configs":
            key = "config:" + self._config_digest(row_id)
        elif "uid" in self._all_columns(table):
            row = self.conn.execute(f"SELECT uid FROM `{table}` WHERE id = ?", (row_id,)).fetchone()
            key = None if row is None else f"{table}:{_uid_text(row['uid'])}"
        else:
            key = f"{table}:unkeyed"
        self._keys[cache_key] = key
        return key

    def prime(self, table, rows):
        """Caches the keys of `rows` of `table` (read with their `id`), so a
        sector's own references cost no extra query."""
        if "uid" not in self._all_columns(table):
            return
        for row in rows:
            if "id" not in row:
                continue
            self._keys[(table, row["id"])] = f"{table}:{_uid_text(row['uid'])}"

    def _config_digest(self, config_id):
        lines = [self.line("system_configs", row) for row in self.select("system_configs", "id = ?", (config_id,))]
        lines += [self.line("system_config_slots", row)
                  for row in self.select("system_config_slots", "config_id = ?", (config_id,))]
        return _digest(sorted(lines))


def _uid_text(uid):
    if uid is None:
        return "none"
    if isinstance(uid, (bytes, bytearray, memoryview)):
        return bytes(uid).hex()
    return format(int(uid), "x")


def _sector_label(row):
    if row["ring_index"] is None:
        return f"unplaced {row['name']}"
    return f"{row['ring_index']} {row['layer_index']} {row['ring_slot_index']}"


def _sector_lines(reader, sector_id, address):
    """The canonical lines of one sector's content: `sector_id` its row's
    id (or `None` for a cell with only address-keyed rows), `address`
    `(ring, layer, slot)` (or `None` when unplaced)."""
    owned = {}

    def keep(table, rows):
        # A facility is found by its sector and its system, a compact
        # remnant by its sector and its star: keep each row once.
        kept = owned.setdefault(table, {})
        for row in rows:
            kept[row.get("id", ("row", len(kept)))] = row
        if table != "sectors":
            reader.prime(table, rows)

    if sector_id is not None:
        keep("sectors", reader.select("sectors", "id = ?", (sector_id,)))
        for table in SECTOR_TABLES:
            keep(table, reader.select(table, "sector_id = ?", (sector_id,)))
        systems = list(owned.get("star_systems", {}))
        for table in SYSTEM_TABLES:
            keep(table, reader.select_in(table, "star_system_id", systems))
        stars = list(owned.get("stars", {}))
        for table in STAR_TABLES:
            keep(table, reader.select_in(table, "star_id", stars))
        for table, column, parent in CHILD_TABLES:
            keep(table, reader.select_in(table, column, owned.get(parent, {})))
    if address is not None:
        for table in ADDRESS_TABLES:
            keep(table, reader.select(table, "ring_index = ? AND layer_index = ? AND ring_slot_index = ?", address))
    return sorted(reader.line(table, row) for table, rows in owned.items() for row in rows.values())


def sector_digest(conn, sector_id=None, address=None):
    """
    The fingerprint of one sector: its row `sector_id`, or the cell at
    `address` `(ring, layer, slot)` (with its sector row, if there is one).

    Returns:
        str: The hex SHA-256 of the sector's sorted canonical lines.
    """
    if sector_id is None and address is None:
        raise ValueError("give a sector id or an address")
    if sector_id is None:
        row = conn.execute("SELECT id FROM sectors WHERE ring_index = ? AND layer_index = ? AND ring_slot_index = ?",
                           tuple(address)).fetchone()
        sector_id = None if row is None else row["id"]
    elif address is None:
        row = conn.execute("SELECT ring_index, layer_index, ring_slot_index FROM sectors WHERE id = ?",
                           (sector_id,)).fetchone()
        if row is not None and row["ring_index"] is not None:
            address = (row["ring_index"], row["layer_index"], row["ring_slot_index"])
    return _digest(_sector_lines(_Reader(conn), sector_id, address))


def plan_digest(conn):
    """The fingerprint of the galaxy's plan (`PLAN_TABLES`)."""
    reader = _Reader(conn)
    lines = [reader.line(table, row) for table in PLAN_TABLES for row in reader.select(table, "1 = 1")]
    return _digest(sorted(lines))


def _cells(conn, rings=None, addresses=None):
    """`[(label, sector id or None, address or None)]` of every sector and
    address-keyed cell in the region, in address order (unplaced sectors
    last, by name)."""
    cells = {}
    for row in conn.execute("SELECT id, ring_index, layer_index, ring_slot_index, name FROM sectors").fetchall():
        if row["ring_index"] is None:
            cells[(1, row["name"] or "", row["id"])] = (_sector_label(row), row["id"], None)
        else:
            address = (row["ring_index"], row["layer_index"], row["ring_slot_index"])
            cells[(0,) + address] = (_sector_label(row), row["id"], address)
    for table in ADDRESS_TABLES:
        for row in conn.execute(f"SELECT DISTINCT ring_index, layer_index, ring_slot_index FROM `{table}`").fetchall():
            address = (row["ring_index"], row["layer_index"], row["ring_slot_index"])
            cells.setdefault((0,) + address, (" ".join(str(n) for n in address), None, address))
    chosen = []
    for order in sorted(cells):
        label, sector_id, address = cells[order]
        if rings is not None or addresses is not None:
            if address is None:
                continue
            if not ((rings and address[0] in rings) or (addresses and address in addresses)):
                continue
        chosen.append((label, sector_id, address))
    return chosen


def region_fingerprint(conn, rings=None, addresses=None, on_progress=None):
    """
    The fingerprint of a region: the sectors on `rings` and at `addresses`
    (`(ring, layer, slot)` tuples), or with neither the whole galaxy, plan
    included.

    Args:
        conn (Connection): An open connection.
        rings (iterable of int, optional): Whole rings to include.
        addresses (iterable of tuple, optional): Single sectors to include.
        on_progress (callable, optional): `on_progress(label, done, total)`
            after each sector.

    Returns:
        Fingerprint: The per-sector digests, the plan's and the region's.
    """
    rings = None if rings is None else set(rings)
    addresses = None if addresses is None else {tuple(a) for a in addresses}
    reader = _Reader(conn)
    cells = _cells(conn, rings, addresses)
    sectors = []
    for done, (label, sector_id, address) in enumerate(cells, start=1):
        sectors.append((label, _digest(_sector_lines(reader, sector_id, address))))
        if on_progress is not None:
            on_progress("sectors", done, len(cells))
    plan = plan_digest(conn) if rings is None and addresses is None else None
    lines = ([f"plan {plan}"] if plan is not None else []) + [f"{label} {digest}" for label, digest in sectors]
    return Fingerprint(sectors=sectors, plan=plan, region=_digest(lines))
