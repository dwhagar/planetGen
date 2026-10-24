# planetgen/db/objectindex.py

"""
The object index builder (PERF.140, docs/design/operations/cache-warm-up.md; the table is DB.79).

`object_index` has one row per object of every kind the site lists. This module fills it from the live tables and keeps
it in step:

- `build(conn, kinds)` writes every object of the kinds, row by row, in chunks of `CHUNK` source ids, and removes the rows
  of objects that are gone. A row is written with `INSERT ... ON DUPLICATE KEY UPDATE` that changes `indexed_at` only when
  a column of the row changed, so a build over unchanged data rewrites no row (the server leaves an identical row alone).
- `refresh(conn, kinds, since)` writes only the objects changed since `since` (a database time): a system, phenomenon
  or facility is changed when its `modified_at` is, and the bodies of a system when the system's is (a change to a child
  row bumps its system, `schema.sql` v27). It then counts each kind against its table and, when the counts differ, sweeps
  the rows of objects that are gone. What the orbital update changes (which sector holds an object) leaves `modified_at`
  alone, so a refresh after it is `build`.
- `check(conn, kinds)` compares the index with the live tables, both ways, and says what is missing, extra or different.
- `status(conn)` gives the size and age the admin page shows.

The state of the last run of each kind (when it ran, what it covered and the time `refresh` starts from next) is kept in
Redis, `planetgen:objindex:<database>`; without it (a Redis that lost its data) the next refresh is a build.

Kinds are `KINDS`: `system`, `star`, `planet`, `moon`, `belt`, `comet`, `facility` and the phenomenon kinds. A phenomenon
is indexed when it stands alone (a black hole or neutron star that is a normal system's own star is on that system's
page already, `query._phenomena_union`).
"""

import json
import time

import redis

from planetgen.db import query
from planetgen.galaxy import objectref
from planetgen.physics import constants
from planetgen.queue import redisclient

CHUNK = 50_000
"""int: Source ids one statement covers; a chunk commits before the next starts."""

SAFETY_SECONDS = 2.0
"""float: How far before a run's start `refresh` goes back, for a change that was committing while the run started."""

KM_PER_LY = constants.LY_TO_AU * (constants.AU_M / 1000.0)
"""float: Kilometres in a light year."""

SOLAR_MASS_KG = constants.SOLAR_MASS_TO_KG

COLUMNS = ("kind", "uid", "source_id", "name", "sector_id", "system_id", "parent_uid", "object_class", "mass_kg",
           "radius_km", "placed")
"""tuple: The columns a source writes, in this order; `indexed_at` is the table's own."""

BODY_KINDS = objectref.BODY_KINDS
PHENOMENON_KINDS = objectref.PHENOMENON_KINDS

_SYSTEM_STAR_TYPE = ("(SELECT st.star_type FROM stars st WHERE st.star_system_id = t.id"
                 " ORDER BY (st.role = 'primary') DESC, st.id LIMIT 1)")
_SYSTEM_MASS = "(SELECT SUM(st.mass_kg) FROM stars st WHERE st.star_system_id = t.id)"


def _select(kind, table, name, sector, system, parent, object_class, mass, radius, placed, joins="", where="",
            scope="own"):
    """A kind's source: the index columns as `SELECT` expressions over `table` (alias `t`) and its `joins`, an extra
    `WHERE`, and how a refresh finds what changed (`scope`: `own`, the row's `modified_at`; `system`, its system's)."""
    return {"kind": kind, "table": table, "scope": scope, "where": where, "from": f"{table} t {joins}".strip(),
            "columns": (f"'{kind}' AS kind, t.uid AS uid, t.id AS source_id, {name} AS name, {sector} AS sector_id, "
                        f"{system} AS system_id, {parent} AS parent_uid, {object_class} AS object_class, "
                        f"{mass} AS mass_kg, {radius} AS radius_km, {placed} AS placed")}


def _phenomenon_sources():
    descriptors = {type_label: (table, descriptor, radius) for table, type_label, descriptor, radius in
                   query._PHENOMENON_TABLES}
    sources = {}
    for kind in PHENOMENON_KINDS:
        table, descriptor, radius = descriptors[kind]
        mass_column = query._PHENOMENON_MASS_COLUMNS.get(table)
        sources[kind] = _select(
            kind, table, "t.name", "t.sector_id", "NULL", "NULL", descriptor,
            f"t.{mass_column} * {SOLAR_MASS_KG!r}" if mass_column else "NULL",
            f"t.radius_ly * {KM_PER_LY!r}" if radius == "radius_ly" else "NULL", "(t.center_x_pc IS NOT NULL)",
            where="t.star_id IS NULL" if table in ("black_holes", "neutron_stars") else "")
    return sources


def _sources():
    ss = "JOIN star_systems ss ON ss.id = t.star_system_id"
    placed = "(ss.galaxy_x_pc IS NOT NULL)"
    sources = {
        "system": _select("system", "star_systems", "t.name", "t.sector_id", "t.id", "NULL", _SYSTEM_STAR_TYPE,
                          _SYSTEM_MASS, "NULL", "(t.galaxy_x_pc IS NOT NULL)"),
        "star": _select("star", "stars", "t.name", "ss.sector_id", "ss.id", "ss.uid", "t.star_type", "t.mass_kg",
                        "t.radius_km", placed, joins=ss, scope="system"),
        "planet": _select("planet", "planets", "t.name", "ss.sector_id", "ss.id", "ss.uid", "t.planet_class",
                          "t.mass_kg", "t.radius_km", placed, joins=ss, scope="system"),
        "moon": _select("moon", "moons", "t.name", "ss.sector_id", "ss.id", "pl.uid", "t.planet_class", "t.mass_kg",
                        "t.radius_km", placed, joins=f"{ss} JOIN planets pl ON pl.id = t.planet_id", scope="system"),
        "belt": _select("belt", "asteroid_belts", "CONCAT('Asteroid belt ', t.orbital_index + 1)", "ss.sector_id",
                        "ss.id", "ss.uid", "t.density", "t.mass_kg", "NULL", placed, joins=ss, scope="system"),
        "comet": _select("comet", "comets", "t.name", "ss.sector_id", "ss.id", "ss.uid", "t.period_class", "NULL",
                         "t.nucleus_diameter_km / 2", placed, joins=ss, scope="system"),
        "facility": _select("facility", "facilities", "t.name", "t.sector_id", "t.star_system_id", "NULL", "t.kind",
                            "t.mass_kg", "NULL", "(t.center_x_pc IS NOT NULL)"),
    }
    sources.update(_phenomenon_sources())
    return sources


SOURCES = _sources()
"""dict: Kind -> its source (`kind`, `table`, `columns`, `from`, `where`, `scope`)."""

KINDS = tuple(SOURCES)
"""tuple: Every kind the index holds."""


class UnknownKind(ValueError):
    """A kind the index does not hold."""


def _kinds(kinds):
    kinds = list(kinds or KINDS)
    unknown = [kind for kind in kinds if kind not in SOURCES]
    if unknown:
        raise UnknownKind(f"no such kind of object {unknown[0]!r}; one of {', '.join(KINDS)}")
    return kinds


def _columns_changed():
    """The `ON DUPLICATE KEY UPDATE` list: `indexed_at` moves only when a column of the row differs (it comes first, so
    it compares the old values), then each column takes the new one from `fresh`, the derived table the rows come from."""
    assigned = COLUMNS[2:]
    same = " AND ".join(f"object_index.{column} <=> fresh.{column}" for column in assigned)
    return (f"object_index.indexed_at = IF({same}, object_index.indexed_at, CURRENT_TIMESTAMP(3)), "
            + ", ".join(f"object_index.{column} = fresh.{column}" for column in assigned))


def _upsert_sql(source, scope_sql):
    where = [part for part in (source["where"], scope_sql, "t.id BETWEEN %s AND %s", "t.uid IS NOT NULL") if part]
    return (f"INSERT INTO object_index ({', '.join(COLUMNS)}) SELECT * FROM (SELECT {source['columns']}"
            f" FROM {source['from']} WHERE {' AND '.join(where)}) AS fresh ON DUPLICATE KEY UPDATE {_columns_changed()}")


def _live(source, scope_sql=""):
    """The objects of `source` as a derived table's `SELECT`, `scope_sql` an extra condition on the alias `t`."""
    where = " AND ".join(part for part in (source["where"], scope_sql, "t.uid IS NOT NULL") if part)
    return f"SELECT {source['columns']} FROM {source['from']} WHERE {where}"


def _write(conn, kind, scope_sql="", params=()):
    """Writes the objects of `kind` that `scope_sql` (an extra condition on the alias `t`) selects, in chunks of source
    ids. Returns how many rows the server changed or added."""
    source = SOURCES[kind]
    row = conn.execute(f"SELECT MIN(id) AS low, MAX(id) AS high FROM {source['table']}").fetchone()
    if row["low"] is None:
        return 0
    sql = _upsert_sql(source, scope_sql)
    changed = 0
    for start in range(row["low"], row["high"] + 1, CHUNK):
        cursor = conn.execute(sql, (*params, start, min(start + CHUNK - 1, row["high"])))
        changed += max(cursor.rowcount, 0)
        conn.commit()
    return changed


def _live_count(conn, kind):
    source = SOURCES[kind]
    return conn.execute(f"SELECT COUNT(*) AS n FROM ({_live(source)}) live").fetchone()["n"]


def _indexed_count(conn, kind):
    return conn.execute("SELECT COUNT(*) AS n FROM object_index WHERE kind = %s", (kind,)).fetchone()["n"]


def _sweep(conn, kind):
    """Deletes the rows of `kind` whose object is gone (or has no `uid` any more). Returns how many."""
    cursor = conn.execute(
        f"DELETE oi FROM object_index oi LEFT JOIN ({_live(SOURCES[kind])}) live ON live.uid = oi.uid"
        " WHERE oi.kind = %s AND live.uid IS NULL", (kind,))
    conn.commit()
    return max(cursor.rowcount, 0)


def _now(conn):
    """The database's clock, as the text a `modified_at` comparison takes, `SAFETY_SECONDS` back."""
    return conn.execute("SELECT DATE_FORMAT(NOW(3) - INTERVAL %s SECOND, '%%Y-%%m-%%d %%H:%%i:%%s.%%f') AS now",
                        (SAFETY_SECONDS,)).fetchone()["now"]


def build(conn, kinds=None, database=None):
    """
    Writes every object of `kinds` (default all) and removes the rows of objects that are gone.

    Returns:
        dict: Kind -> `{"changed", "removed", "rows", "seconds"}`.
    """
    result = {}
    for kind in _kinds(kinds):
        started = time.monotonic()
        since = _now(conn)
        changed = _write(conn, kind)
        removed = _sweep(conn, kind)
        rows = _indexed_count(conn, kind)
        result[kind] = {"changed": changed, "removed": removed, "rows": rows, "seconds": time.monotonic() - started}
        _remember(conn, database, kind, "build", since, rows, result[kind])
    return result


def refresh(conn, kinds=None, database=None):
    """
    Writes the objects of `kinds` (default all) that changed since the kind's last run, and sweeps a kind whose rows
    are not as many as its objects. A kind with no record of a run is built.

    Returns:
        dict: As `build`, with `mode` (`build` or `refresh`) for each kind.
    """
    result = {}
    for kind in _kinds(kinds):
        state = recorded(database).get(kind)
        if not state or not state.get("since"):
            result[kind] = {**build(conn, [kind], database)[kind], "mode": "build"}
            continue
        started = time.monotonic()
        since = _now(conn)
        source = SOURCES[kind]
        if source["scope"] == "own":
            scope_sql = "t.modified_at > %s"
        else:
            scope_sql = "t.star_system_id IN (SELECT m.id FROM star_systems m WHERE m.modified_at > %s)"
        changed = _write(conn, kind, scope_sql, (state["since"],))
        removed = 0
        if _live_count(conn, kind) != _indexed_count(conn, kind):
            removed = _sweep(conn, kind)
            if _live_count(conn, kind) != _indexed_count(conn, kind):
                # An object that changed under the watermark was not seen: write the kind whole.
                changed += _write(conn, kind)
        rows = _indexed_count(conn, kind)
        result[kind] = {"changed": changed, "removed": removed, "rows": rows, "mode": "refresh",
                        "seconds": time.monotonic() - started}
        _remember(conn, database, kind, "refresh", since, rows, result[kind])
    return result


def check(conn, kinds=None):
    """
    Compares the index with the live tables.

    Returns:
        dict: Kind -> `{"live", "indexed", "missing", "extra", "different"}`: the objects there are, the rows there
        are, objects with no row, rows with no object, and rows that do not hold what their object does.
    """
    report = {}
    for kind in _kinds(kinds):
        live = f"({_live(SOURCES[kind])}) live"
        missing = conn.execute(
            f"SELECT COUNT(*) AS n FROM {live} LEFT JOIN object_index oi ON oi.kind = live.kind AND oi.uid = live.uid"
            " WHERE oi.uid IS NULL").fetchone()["n"]
        compared = " OR ".join(f"NOT (oi.{column} <=> live.{column})" for column in COLUMNS[2:])
        different = conn.execute(
            f"SELECT COUNT(*) AS n FROM {live} JOIN object_index oi ON oi.kind = live.kind AND oi.uid = live.uid"
            f" WHERE {compared}").fetchone()["n"]
        extra = conn.execute(
            f"SELECT COUNT(*) AS n FROM object_index oi LEFT JOIN {live} ON oi.kind = live.kind AND oi.uid = live.uid"
            " WHERE oi.kind = %s AND live.uid IS NULL", (kind,)).fetchone()["n"]
        report[kind] = {"live": _live_count(conn, kind), "indexed": _indexed_count(conn, kind), "missing": missing,
                        "extra": extra, "different": different}
    return report


def size_bytes(conn):
    """The bytes `object_index` takes on disk (data and indexes), 0 when the table is not there."""
    row = conn.execute("SELECT COALESCE(DATA_LENGTH + INDEX_LENGTH, 0) AS size FROM information_schema.TABLES"
                       " WHERE TABLE_SCHEMA = DATABASE() AND TABLE_NAME = 'object_index'").fetchone()
    return int(row["size"]) if row else 0


def status(conn, database=None):
    """
    What the admin page shows: the table's size and, for each kind, its rows and when it was last run.

    Returns:
        dict: `{"size_bytes", "kinds": [{"kind", "rows", "mode", "built_at", "seconds"}]}`.
    """
    counts = {row["kind"]: row["n"] for row in
              conn.execute("SELECT kind, COUNT(*) AS n FROM object_index GROUP BY kind").fetchall()}
    state = recorded(database)
    kinds = []
    for kind in KINDS:
        last = state.get(kind) or {}
        kinds.append({"kind": kind, "rows": counts.get(kind, 0), "mode": last.get("mode"),
                      "built_at": last.get("built_at"), "seconds": last.get("seconds")})
    return {"size_bytes": size_bytes(conn), "kinds": kinds}


# --- what the last runs did ---------------------------------------------------------------------------------------

def state_key(database):
    """The Redis hash of the last run of each kind of `database`."""
    return f"planetgen:objindex:{database}"


def _client():
    return redisclient.client(timeout=0.5, retries=1, decode_responses=True)


SIZE_FIELD = "_size_bytes"
"""str: The hash field that holds the table's size when the last run ended."""


def recorded(database):
    """Kind -> the record of its last run (`{"mode", "since", "built_at", "rows", "seconds"}`); `{}` without Redis."""
    if database is None:
        return {}
    try:
        raw = _client().hgetall(state_key(database))
    except (redis.RedisError, OSError):
        return {}
    state = {}
    for kind, text in raw.items():
        if kind == SIZE_FIELD:
            continue
        try:
            state[kind] = json.loads(text)
        except ValueError:
            continue
    return state


def recorded_size(database):
    """The bytes the table took when the last run ended (`0` without a record)."""
    try:
        return int(_client().hget(state_key(database), SIZE_FIELD) or 0)
    except (redis.RedisError, OSError, ValueError):
        return 0


def reserved_bytes():
    """
    The disk the object index takes of the cache budget (`cache.max_mb`), all databases together: the index is kept
    first, so the tile cache's budget is what is left (`tilecache.max_cache_bytes`, PERF.140).
    """
    try:
        client = _client()
        return sum(int(client.hget(key, SIZE_FIELD) or 0) for key in client.scan_iter(match=state_key("*"), count=100))
    except (redis.RedisError, OSError, ValueError):
        return 0


def _remember(conn, database, kind, mode, since, rows, outcome):
    """Keeps the record of a kind's run, and the table's size (which the tile cache's budget leaves out)."""
    if database is None:
        return
    record = {"mode": mode, "since": since, "built_at": time.time(), "rows": rows,
              "seconds": round(outcome["seconds"], 2)}
    try:
        pipe = _client().pipeline(transaction=True)
        pipe.hset(state_key(database), kind, json.dumps(record))
        pipe.hset(state_key(database), SIZE_FIELD, size_bytes(conn))
        pipe.execute()
    except (redis.RedisError, OSError):
        pass


def forget(database):
    """Drops the record of the runs, so the next `refresh` builds."""
    try:
        _client().delete(state_key(database))
    except (redis.RedisError, OSError):
        pass
