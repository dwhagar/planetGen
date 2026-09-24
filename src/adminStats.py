# src/adminStats.py

"""
Read-only statistics about one content database, for the admin stats page
(`html/adminstats.py`, via `GET /api/admin/stats` and `GET
/api/admin/duplicate-names` in `html/api/admin.py`).

Kept out of `queryDb.py` because none of it is a query a visitor's page
runs: it's bookkeeping about the database itself (table sizes, schema
version, row timestamps, how many names the uniqueness rules had to
decorate).

Everything here is written to stay cheap on a galaxy-scale database,
because the API serves every page from one small thread pool and a slow
admin request would stall the whole site (see `docs/apache-deployment.md`):

- Per-table row counts come from `information_schema.tables`, which is an
  estimate for InnoDB (and can lag by up to a day on MySQL 8, see
  `information_schema_stats_expiry`). Only `sectors` and `star_systems`
  get an exact `COUNT(*)`, the same two `/api/databases` already counts.
- Newest/last-modified times use the primary key and the v27
  `idx_*_modified_at` indexes, never a scan of `created_at`.
- Duplicate-name lookups look up exact candidate names through each
  table's `idx_*_name` index instead of pattern-matching the whole table.
"""

from stellarObjects import _db
from stellarObjects.names import COMPANION_SUFFIXES, DIMINUTIVE_PREFIXES, GREEK_LETTERS, ROMAN_NUMERALS_BY_VALUE
from stellarObjects.nameUniqueness import strip_decoration

NAME_REGISTRIES = (
    ("sector", "sector_name_registry", "occurrence_count > 1"),
    ("system", "system_name_registry", "(occurrence_count > 1 OR diminutive_index IS NOT NULL)"),
    ("body", "body_name_registry", "(occurrence_count > 1 OR suffix_index IS NOT NULL)"),
)
"""tuple: `(level, registry table, collided condition)` for each level of
the naming hierarchy -- see `stellarObjects/nameUniqueness.py`. Every
name gets a registry row when it's first used, so only the rows matching
the condition had to be decorated: a sector or system name used more than
once (Greek/Roman), a system named after a sector (a diminutive), or a
planet/moon name that collided with anything (a companion suffix).

None of these columns is indexed, so counting them scans each registry
once. That's one row per name, which is fine for a page only admins load,
but it's why the stats page doesn't auto-refresh."""

_NAMED_TABLES = (
    ("sector", "sectors"),
    ("system", "star_systems"),
    ("planet", "planets"),
    ("moon", "moons"),
)

_CANDIDATE_CHUNK = 1000
"""int: Most names bound into one `IN (...)` list."""


def _format_time(value):
    return value.strftime("%Y-%m-%d %H:%M:%S") if value is not None else None


def server_info(conn):
    """
    The MySQL server's version, uptime and connection count. Uptime and
    connections come from `SHOW GLOBAL STATUS`, which some hosted setups
    restrict; those come back `None` rather than failing the page.

    Returns:
        dict: `version`, `uptime_seconds`, `threads_connected`.
    """
    info = {"version": None, "uptime_seconds": None, "threads_connected": None}
    try:
        info["version"] = conn.execute("SELECT VERSION() AS v").fetchone()["v"]
    except Exception:
        pass
    try:
        rows = conn.execute(
            "SHOW GLOBAL STATUS WHERE Variable_name IN ('Uptime', 'Threads_connected')"
        ).fetchall()
        status = {row["Variable_name"]: row["Value"] for row in rows}
        if "Uptime" in status:
            info["uptime_seconds"] = int(status["Uptime"])
        if "Threads_connected" in status:
            info["threads_connected"] = int(status["Threads_connected"])
    except Exception:
        pass
    return info


def schema_version(conn):
    """The database's applied schema version, or `None` before `migrateDb.py`
    has ever run against it (no `schema_migrations` table yet)."""
    try:
        row = conn.execute("SELECT MAX(version) AS version FROM schema_migrations").fetchone()
    except Exception:
        return None
    return row["version"] if row else None


def table_stats(conn):
    """
    One entry per base table in the connection's database, sorted by name.

    Returns:
        list[dict]: `name`, `approx_rows` (`information_schema`'s estimate),
            `data_bytes`, `index_bytes`.
    """
    rows = conn.execute(
        "SELECT table_name AS name, table_rows AS approx_rows, "
        "data_length AS data_bytes, index_length AS index_bytes "
        "FROM information_schema.tables "
        "WHERE table_schema = DATABASE() AND table_type = 'BASE TABLE' "
        "ORDER BY table_name"
    ).fetchall()
    return [
        {
            "name": row["name"],
            "approx_rows": int(row["approx_rows"] or 0),
            "data_bytes": int(row["data_bytes"] or 0),
            "index_bytes": int(row["index_bytes"] or 0),
        }
        for row in rows
    ]


def exact_counts(conn):
    """Exact row counts for the two tables small enough to always count
    (`sectors`, `star_systems`), each `None` if the table is missing."""
    counts = {}
    for table in ("sectors", "star_systems"):
        try:
            counts[table] = conn.execute(f"SELECT COUNT(*) AS n FROM {table}").fetchone()["n"]
        except Exception:
            counts[table] = None
    return counts


def timestamp_stats(conn):
    """
    For every table with v27 row timestamps (`_db.TIMESTAMPED_TABLES`):
    when its newest row was created and when any row last changed. A
    database still on a schema older than v27 reports `None` for both.

    Returns:
        list[dict]: `table`, `newest_created_at`, `last_modified_at`
            (`"YYYY-MM-DD HH:MM:SS"` strings, or `None`).
    """
    result = []
    for table in _db.TIMESTAMPED_TABLES:
        newest = last_modified = None
        try:
            row = conn.execute(f"SELECT created_at FROM {table} ORDER BY id DESC LIMIT 1").fetchone()
            newest = row["created_at"] if row else None
            last_modified = conn.execute(f"SELECT MAX(modified_at) AS t FROM {table}").fetchone()["t"]
        except Exception:
            pass
        result.append({
            "table": table,
            "newest_created_at": _format_time(newest),
            "last_modified_at": _format_time(last_modified),
        })
    return result


def name_collision_summary(conn):
    """
    How many base names had to be made unique at each level of the
    naming hierarchy, plus how many distinct base names that is across
    all three (one name can collide at more than one level).

    Returns:
        dict: `distinct_base_names`, and per level (`sector`/`system`/
            `body`) the number of base names that collided at that level.
    """
    summary = {}
    selects = []
    for level, table, collided in NAME_REGISTRIES:
        try:
            summary[level] = conn.execute(f"SELECT COUNT(*) AS n FROM {table} WHERE {collided}").fetchone()["n"]
            selects.append(f"SELECT base_name FROM {table} WHERE {collided}")
        except Exception:
            summary[level] = None
    if selects:
        summary["distinct_base_names"] = conn.execute(
            f"SELECT COUNT(*) AS n FROM ({' UNION '.join(selects)}) AS names"
        ).fetchone()["n"]
    else:
        summary["distinct_base_names"] = None
    return summary


def decorated_name_candidates(base_name):
    """
    Every name the uniqueness rules can give a row whose base name is
    `base_name`: the bare name, each Greek-letter prefix, the "Alpha
    <name> <roman>" overflow forms, each diminutive prefix (a system that
    collided with a sector) and each companion suffix (a planet or moon).
    `strip_decoration` maps each of these back to `base_name`.

    Returns:
        list[str]
    """
    names = [base_name]
    names.extend(f"{letter} {base_name}" for letter in GREEK_LETTERS)
    names.extend(f"{GREEK_LETTERS[0]} {base_name} {numeral}" for numeral in ROMAN_NUMERALS_BY_VALUE.values())
    names.extend(f"{prefix} {base_name}" for prefix in DIMINUTIVE_PREFIXES)
    names.extend(f"{base_name} {suffix}" for suffix in COMPANION_SUFFIXES)
    return names


def _registry_union_sql():
    parts = [
        f"SELECT base_name, '{level}' AS level FROM {table} WHERE {collided}"
        for level, table, collided in NAME_REGISTRIES
    ]
    return " UNION ALL ".join(parts)


def _rows_named(conn, kind, table, names):
    """Rows of one named table whose name is in `names`, with enough
    context to link to them (planets and moons have no page of their own,
    so they carry their system's id and name)."""
    if kind == "sector":
        select = "SELECT t.id, t.name FROM sectors t"
    elif kind == "system":
        select = "SELECT t.id, t.name, t.sector_id FROM star_systems t"
    elif kind == "planet":
        select = (
            "SELECT t.id, t.name, s.id AS star_system_id, s.name AS system_name "
            "FROM planets t JOIN star_systems s ON s.id = t.star_system_id"
        )
    else:
        select = (
            "SELECT t.id, t.name, s.id AS star_system_id, s.name AS system_name "
            "FROM moons t JOIN planets p ON p.id = t.planet_id "
            "JOIN star_systems s ON s.id = p.star_system_id"
        )
    names = list(names)
    rows = []
    for start in range(0, len(names), _CANDIDATE_CHUNK):
        chunk = names[start:start + _CANDIDATE_CHUNK]
        placeholders = ", ".join("?" for _ in chunk)
        rows.extend(conn.execute(f"{select} WHERE t.name IN ({placeholders})", tuple(chunk)).fetchall())
    return rows


def duplicate_names(conn, limit=100, offset=0):
    """
    One page of the base names the uniqueness rules had to decorate
    (see `NAME_REGISTRIES`, alphabetical), each
    with every sector, system, planet and moon that now carries a form of
    it.

    Args:
        conn: An open, read-only connection.
        limit (int): Base names per page.
        offset (int): Base names to skip.

    Returns:
        dict: `total` (distinct base names overall) and `items`, each
            `{"base_name", "levels": [...], "rows": [{"kind", "id",
            "name", ...}]}`. `levels` says which registries hold the base
            name (`sector`/`system`/`body`). A system row carries
            `sector_id`; a planet or moon row carries `star_system_id`
            and `system_name`.
    """
    union_sql = _registry_union_sql()
    total = conn.execute(
        f"SELECT COUNT(DISTINCT base_name) AS n FROM ({union_sql}) AS names"
    ).fetchone()["n"]
    page = conn.execute(
        f"SELECT base_name, GROUP_CONCAT(DISTINCT level ORDER BY level) AS levels "
        f"FROM ({union_sql}) AS names GROUP BY base_name ORDER BY base_name LIMIT ? OFFSET ?",
        (limit, offset),
    ).fetchall()

    items = []
    by_base = {}
    candidates = set()
    for row in page:
        item = {"base_name": row["base_name"], "levels": row["levels"].split(","), "rows": []}
        items.append(item)
        by_base[row["base_name"].casefold()] = item
        candidates.update(decorated_name_candidates(row["base_name"]))

    if candidates:
        for kind, table in _NAMED_TABLES:
            for row in _rows_named(conn, kind, table, sorted(candidates)):
                item = by_base.get(strip_decoration(row["name"]).casefold())
                if item is None:
                    continue
                entry = {"kind": kind, "id": row["id"], "name": row["name"]}
                if kind == "system":
                    entry["sector_id"] = row["sector_id"]
                elif kind in ("planet", "moon"):
                    entry["star_system_id"] = row["star_system_id"]
                    entry["system_name"] = row["system_name"]
                item["rows"].append(entry)

    kind_order = {kind: index for index, (kind, _table) in enumerate(_NAMED_TABLES)}
    for item in items:
        item["rows"].sort(key=lambda r: (kind_order[r["kind"]], r["name"].casefold(), r["id"]))
    return {"total": total, "limit": limit, "offset": offset, "items": items}
