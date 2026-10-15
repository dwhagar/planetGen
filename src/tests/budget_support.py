# tests/budget_support.py

"""
The big-galaxy query budget (TEST.133): a database of millions of systems, and a way to count the rows a page's
queries read.

`ensure_big_galaxy` builds (once, then reuses) a database of `PLANETGEN_BUDGET_SYSTEMS` systems (default 2,000,000),
one star each, spread over `PLANETGEN_BUDGET_SECTORS` sectors (default 150,000), by cloning one real saved system.
That takes about eight minutes and 3 GB, which is why the tier is opt-in (`scripts/test.py budget`) and why the
database is kept between runs, named for its size and the schema version so a new schema rebuilds it.

`rows_read(conn, call)` runs one page's query function on a connection and returns how many rows the server read to
answer it (the handler read counters of that session), so the budget test can fail a page that reads far more rows
than it shows, however fast the box happens to be.
"""

import os
import time

import pymysql

from planetgen.db import store
from planetgen.db.store import MySQLConfig

DEFAULT_SYSTEMS = 2_000_000
DEFAULT_SECTORS = 150_000
PREFIX = "planetgen_budget_"
HANDLER_COUNTERS = ("Handler_read_first", "Handler_read_key", "Handler_read_last", "Handler_read_next",
                    "Handler_read_prev", "Handler_read_rnd", "Handler_read_rnd_next")
_STAR_TYPE_SQL = ("CASE WHEN (s.seq * 2654435761) % 100 < 76 THEN 'M5V Red Dwarf' WHEN (s.seq * 2654435761) % 100 < 88 THEN 'K5V Orange Dwarf' "
                  "WHEN (s.seq * 2654435761) % 100 < 96 THEN 'G2V Yellow Main Sequence Star' WHEN (s.seq * 2654435761) % 100 < 99 THEN 'F5V Yellow-White Dwarf' "
                  "ELSE 'A1V White Main Sequence Star' END")
"""str: A star type for a cloned star, in about the mix of a real galaxy (three in four red dwarfs), so a filter by type is as selective as it is there."""
FULLTEXT_INDEXES = (("sectors", "ft_sectors_name"), ("star_systems", "ft_star_systems_name"), ("stars", "ft_stars_name"))


def sizes():
    return (int(os.environ.get("PLANETGEN_BUDGET_SYSTEMS", DEFAULT_SYSTEMS)),
            int(os.environ.get("PLANETGEN_BUDGET_SECTORS", DEFAULT_SECTORS)))


def database_name(systems, sectors):
    return f"{PREFIX}{systems}_{sectors}_v{store.SCHEMA_VERSION}"


def _server(kwargs, database=None):
    return pymysql.connect(**kwargs, database=database, autocommit=True, read_timeout=3600, write_timeout=3600)


def _exists(kwargs, name):
    conn = _server(kwargs)
    try:
        with conn.cursor() as cur:
            cur.execute("SELECT 1 FROM information_schema.SCHEMATA WHERE SCHEMA_NAME = %s", (name,))
            return cur.fetchone() is not None
    finally:
        conn.close()


def _is_complete(kwargs, name, systems):
    conn = _server(kwargs, name)
    try:
        with conn.cursor() as cur:
            cur.execute("SELECT 1 FROM information_schema.TABLES WHERE TABLE_SCHEMA = %s AND TABLE_NAME = 'budget_built'", (name,))
            if cur.fetchone() is None:
                return False
            cur.execute("SELECT systems FROM budget_built")
            row = cur.fetchone()
            return row is not None and row[0] == systems
    except pymysql.MySQLError:
        return False
    finally:
        conn.close()


def ensure_big_galaxy(server_kwargs, systems=None, sectors=None, log=print):
    """The `MySQLConfig` of the big database, built if it is not there (or was left half-built). Other budget
    databases of another size or schema version are dropped, to keep the disk."""
    default_systems, default_sectors = sizes()
    systems = systems or default_systems
    sectors = sectors or default_sectors
    name = database_name(systems, sectors)
    config = MySQLConfig(database=name, **server_kwargs)
    if os.environ.get("PLANETGEN_BUDGET_REBUILD") != "1" and _exists(server_kwargs, name) and _is_complete(server_kwargs, name, systems):
        return config
    build(server_kwargs, name, systems, sectors, log=log)
    return config


def drop_other_budget_databases(server_kwargs, keep):
    conn = _server(server_kwargs)
    try:
        with conn.cursor() as cur:
            cur.execute("SELECT SCHEMA_NAME FROM information_schema.SCHEMATA WHERE SCHEMA_NAME LIKE %s", (PREFIX + "%",))
            for (found,) in cur.fetchall():
                if found != keep:
                    cur.execute(f"DROP DATABASE `{found}`")
    finally:
        conn.close()


def _seed_one_system(config):
    """One real saved sector holding one real single-star system, to clone."""
    import math

    from planetgen.galaxy.sector import SpaceSector
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem

    sector = SpaceSector("Seed Sector", edge_ly=10.0)
    system_config = SystemConfig()
    system_config.STAR_TYPE = "G2V"
    system_config.PLANETS = False
    system_config.BINARY_SYSTEM = False
    sector.add_system(StarSystem(system_config=system_config), position=(0.0, 0.0, 0.0), system_config=system_config)
    center = (5.0, 5.0, 5.0)
    position = {"center_x_pc": center[0], "center_y_pc": center[1], "center_z_pc": center[2],
                "galactic_radius_pc": math.dist(center, (0.0, 0.0, 0.0))}
    store.save_sector(sector, config=config, galaxy_position=position)


def build(server_kwargs, name, systems, sectors, log=print, drop_others=True):
    """Creates `name` with the project's schema and `systems` cloned systems over `sectors` cloned sectors."""
    started = time.monotonic()
    if drop_others:
        drop_other_budget_databases(server_kwargs, keep=None)
    conn = _server(server_kwargs)
    with conn.cursor() as cur:
        cur.execute(f"CREATE DATABASE `{name}`")
    conn.close()
    config = MySQLConfig(database=name, **server_kwargs)
    store.get_connection(config).close()
    _seed_one_system(config)
    conn = _server(server_kwargs, name)
    cur = conn.cursor()

    def run(sql, what):
        t = time.monotonic()
        cur.execute(sql)
        log(f"budget: {time.monotonic() - t:6.1f} s  {what}")

    def columns(table):
        cur.execute("SELECT COLUMN_NAME FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = %s AND TABLE_NAME = %s "
                    "AND EXTRA NOT LIKE %s ORDER BY ORDINAL_POSITION", (name, table, "%GENERATED%"))
        return [row[0] for row in cur.fetchall()]

    try:
        run("SET SESSION foreign_key_checks = 0, unique_checks = 0", "relaxed checks for the bulk fill")
        for table, index in FULLTEXT_INDEXES:
            run(f"ALTER TABLE {table} DROP INDEX {index}", f"drop {index} while filling")
        run("CREATE TABLE budget_digits (d TINYINT PRIMARY KEY)", "digits")
        cur.execute("INSERT INTO budget_digits VALUES (0),(1),(2),(3),(4),(5),(6),(7),(8),(9)")
        top = max(systems, sectors)
        joins = " JOIN ".join(f"budget_digits d{i}" for i in range(1, 8))
        run("CREATE TABLE budget_numbers (seq INT PRIMARY KEY)", "numbers")
        run("INSERT INTO budget_numbers SELECT seq FROM (SELECT d1.d + 10 * d2.d + 100 * d3.d + 1000 * d4.d + 10000 * d5.d "
            f"+ 100000 * d6.d + 1000000 * d7.d + 1 AS seq FROM {joins}) n "
            f"WHERE seq <= {top}", f"numbers 1 to {top}")
        cur.execute("SELECT id FROM system_configs LIMIT 1")
        config_id = cur.fetchone()[0]

        def clone(table, where, overrides, count):
            names = columns(table)
            select = ", ".join(overrides.get(c, f"t.`{c}`") for c in names)
            run(f"INSERT INTO {table} ({', '.join('`' + c + '`' for c in names)}) SELECT {select} FROM "
                f"(SELECT * FROM {table} WHERE {where} LIMIT 1) t JOIN budget_numbers s ON s.seq <= {count} AND s.seq >= 2",
                f"clone {count - 1} rows into {table}")

        clone("sectors", "center_x_pc IS NOT NULL", {
            "id": "s.seq + 1000000000", "name": "CONCAT('Sector ', CHAR(65 + (s.seq * 7919) % 26), ' ', (s.seq * 2654435761) % 1000003)",
            "ring_index": "s.seq DIV 1000", "layer_index": "0", "ring_slot_index": "s.seq MOD 1000", "uid": "s.seq + 1000000000"}, sectors)
        clone("star_systems", "is_binary = 0 AND sector_id IS NOT NULL", {
            "id": "s.seq + 1000000000", "sector_id": f"1000000000 + 2 + (s.seq * 7) % {sectors - 1}", "system_config_id": str(config_id),
            "name": "CONCAT(ELT(1 + s.seq % 6, 'Alder', 'Brae', 'Cinder', 'Dune', 'Ember', 'Fell'), ' ', (s.seq * 2654435761) % 9000017)",
            "quadrant": "ELT(1 + s.seq % 8, 'I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII')",
            "uid": "UNHEX(LPAD(HEX(s.seq + 1000000000), 20, '0'))", "inside_nebula_id": "NULL", "inside_remnant_id": "NULL"}, systems)
        clone("stars", "role = 'single'", {
            "id": "s.seq + 1000000000", "star_system_id": "s.seq + 1000000000", "name": "CONCAT('Star ', s.seq)",
            "star_type": _STAR_TYPE_SQL,
            "uid": "UNHEX(LPAD(HEX(s.seq + 1000000000), 20, '0'))"}, systems)
        for table, index in FULLTEXT_INDEXES:
            run(f"ALTER TABLE {table} ADD FULLTEXT KEY {index} (name)", f"rebuild {index}")
        run("DROP TABLE budget_numbers", "drop numbers")
        run("DROP TABLE budget_digits", "drop digits")
        run("ANALYZE TABLE star_systems, stars, sectors", "statistics")
        cur.execute("CREATE TABLE budget_built (systems INT)")
        cur.execute("INSERT INTO budget_built VALUES (%s)", (systems,))
        log(f"budget: built {name} in {time.monotonic() - started:.0f} s")
    finally:
        cur.close()
        conn.close()
    store.close_pool(config)


def _handler_counts(conn):
    rows = conn.execute("SHOW SESSION STATUS WHERE Variable_name LIKE ?", ("Handler_read%",)).fetchall()
    values = {row["Variable_name"]: int(row["Value"]) for row in rows}
    return sum(values.get(name, 0) for name in HANDLER_COUNTERS)


def rows_read(conn, call):
    """(rows the server read to run `call(conn)`, seconds). Counts every read of the connection's session, so
    `call` must use only `conn`."""
    before = _handler_counts(conn)
    started = time.monotonic()
    call(conn)
    seconds = time.monotonic() - started
    return _handler_counts(conn) - before, seconds
