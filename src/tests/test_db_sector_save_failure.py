# tests/test_db_sector_save_failure.py

"""
TEST.15 "Sector save fails halfway": `store.save_sector`/`insert_sector`
failing after a sector's systems are written but before its phenomena or
its neighbour links. Nothing may persist (no rows, no orphans, no stale
name reservations), the neighbour named lock must be free again, and the
same sector must save afterwards under its generated names. Also covers
running out of deadlock retries, and a lock wait timeout on the neighbour
lock that clears once the other holder lets go.

Every test here takes the `mysql_config` fixture (see `conftest.py`) --
skipped, not failed, when no MySQL test server is configured/reachable.
"""

import pymysql
import pytest

from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.rogue import RoguePlanet
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.system import StarSystem
from tests.conftest import _test_server_kwargs

# schema_migrations holds the bootstrap version row; id_blocks reservations
# commit on their own connection by design (PERF.13), so a rolled-back save
# leaves a gap there, never rows anywhere else.
_BOOKKEEPING_TABLES = {"schema_migrations", "id_blocks"}

_POSITION = {"center_x_pc": 1.0, "center_y_pc": 2.0, "center_z_pc": 0.5, "galactic_radius_pc": 2.29}
_NEIGHBOR_POSITION = {"center_x_pc": 4.526, "center_y_pc": 2.0, "center_z_pc": 0.5, "galactic_radius_pc": 4.97}

# Repeated names, so the save really renames things (Greek letters) and a
# name left over from a failed attempt would show.
_SYSTEM_NAMES = ["Vestara", "Vestara", "Corvain"]
_ROGUE_NAMES = ["Vestara", "Drifter"]
_SAVED_SYSTEMS = ["Alpha Vestara", "Beta Vestara", "Corvain"]
_SAVED_ROGUES = ["Drifter", "Gamma Vestara"]


def _sector(name, system_names=_SYSTEM_NAMES, rogue_names=_ROGUE_NAMES):
    sector = SpaceSector(name, edge_ly=11.5)
    for index, system_name in enumerate(system_names):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.BINARY_SYSTEM = False
        system = StarSystem(system_config=cfg)
        system.name = system_name
        sector.add_system(system, position=(float(index) + 1.0, 1.0, 1.0), system_config=cfg)
    for index, rogue_name in enumerate(rogue_names):
        rogue = RoguePlanet(SystemConfig(), name=rogue_name)
        sector.add_phenomenon(rogue, "rogue-planet", position=(1.0, float(index) + 2.5, 1.0))
    return sector


def _names(sector):
    return [sector.name] + [e.star_system.name for e in sector.entries] + [e.phenomenon.name for e in sector.phenomena]


def _query(config, sql, params=()):
    conn = store.get_connection(config)
    try:
        return conn.execute(sql, params).fetchall()
    finally:
        conn.close()


def _content_row_counts(config):
    tables = _query(config, "SELECT table_name AS name FROM information_schema.tables "
                            "WHERE table_schema = ? AND table_type = 'BASE TABLE'", (config.database,))
    names = sorted(row["name"] for row in tables if row["name"] not in _BOOKKEEPING_TABLES)
    assert "star_systems" in names and "system_name_registry" in names
    return {name: _query(config, f"SELECT COUNT(*) AS n FROM `{name}`")[0]["n"] for name in names}


def _orphans(config):
    """`[(constraint, child rows)]` for every foreign key with a child row
    pointing at a missing parent, plus `nearest_systems` rows (whose
    `object_id` has no foreign key: it names a row in `object_table`)."""
    columns = _query(config, """
        SELECT constraint_name AS name, table_name AS child, column_name AS col,
               referenced_table_name AS parent, referenced_column_name AS ref
        FROM information_schema.key_column_usage
        WHERE table_schema = ? AND referenced_table_name IS NOT NULL
        ORDER BY table_name, constraint_name, ordinal_position
    """, (config.database,))
    keys = {}
    for row in columns:
        keys.setdefault((row["child"], row["name"], row["parent"]), []).append((row["col"], row["ref"]))
    found = []
    for (child, name, parent), pairs in keys.items():
        joined = " AND ".join(f"p.`{ref}` = c.`{col}`" for col, ref in pairs)
        present = " AND ".join(f"c.`{col}` IS NOT NULL" for col, _ref in pairs)
        n = _query(config, f"SELECT COUNT(*) AS n FROM `{child}` c LEFT JOIN `{parent}` p ON {joined} "
                           f"WHERE {present} AND p.`{pairs[0][1]}` IS NULL")[0]["n"]
        if n:
            found.append((name, n))
    for row in _query(config, "SELECT DISTINCT object_table AS t FROM nearest_systems"):
        n = _query(config, f"SELECT COUNT(*) AS n FROM nearest_systems ns LEFT JOIN `{row['t']}` o "
                           f"ON o.id = ns.object_id WHERE ns.object_table = ? AND o.id IS NULL", (row["t"],))[0]["n"]
        if n:
            found.append((f"nearest_systems.{row['t']}", n))
    return found


def _lock_name(config):
    conn = store.get_connection(config)
    try:
        return store._neighbor_lock_name(conn)
    finally:
        conn.close()


def _lock_is_free(config, name):
    """Asked from a connection outside the pool, as another writer would."""
    other = pymysql.connect(**_test_server_kwargs())
    try:
        with other.cursor() as cur:
            cur.execute("SELECT IS_FREE_LOCK(%s)", (name,))
            return cur.fetchone()[0] == 1
    finally:
        other.close()


def _assert_nothing_persisted(config):
    assert {t: n for t, n in _content_row_counts(config).items() if n} == {}
    assert _query(config, "SELECT base_name FROM system_name_registry") == []
    assert _query(config, "SELECT base_name FROM sector_name_registry") == []
    assert _orphans(config) == []


def _assert_saved_once(config, sector_id):
    systems = _query(config, "SELECT sector_id, name FROM star_systems ORDER BY name")
    assert [(r["sector_id"], r["name"]) for r in systems] == [(sector_id, n) for n in _SAVED_SYSTEMS]
    rogues = _query(config, "SELECT sector_id, name FROM rogue_planets ORDER BY name")
    assert [(r["sector_id"], r["name"]) for r in rogues] == [(sector_id, n) for n in _SAVED_ROGUES]
    registry = {r["base_name"]: r["occurrence_count"] for r in _query(
        config, "SELECT base_name, occurrence_count FROM system_name_registry")}
    assert registry == {"Vestara": 3, "Corvain": 1, "Drifter": 1}
    assert _orphans(config) == []


def _fail_on_second_rogue(monkeypatch, config, lock_name):
    real = store._PHENOMENON_INSERTERS["rogue-planet"]
    calls = []

    def inserter(conn, phenomenon, **kwargs):
        calls.append(phenomenon.name)
        if len(calls) == 1:
            return real(conn, phenomenon, **kwargs)
        conn.flush()
        # The systems (and the first phenomenon) really are written.
        assert conn.execute("SELECT COUNT(*) AS n FROM star_systems").fetchone()["n"] == 3
        assert conn.execute("SELECT COUNT(*) AS n FROM rogue_planets").fetchone()["n"] == 1
        raise pymysql.err.IntegrityError(1452, "Cannot add or update a child row (injected)")

    monkeypatch.setitem(store._PHENOMENON_INSERTERS, "rogue-planet", inserter)


def _fail_linking_neighbours(monkeypatch, config, lock_name):
    def link(conn, sector_id):
        assert conn.execute("SELECT COUNT(*) AS n FROM star_systems").fetchone()["n"] == 3
        assert conn.execute("SELECT COUNT(*) AS n FROM rogue_planets").fetchone()["n"] == 2
        assert not _lock_is_free(config, lock_name)
        raise pymysql.err.IntegrityError(1452, "Cannot add or update a child row (injected)")

    monkeypatch.setattr(store, "_add_sector_to_nearest", link)


@pytest.mark.parametrize("inject", [_fail_on_second_rogue, _fail_linking_neighbours],
                         ids=["phenomena", "neighbours"])
def test_a_save_failing_after_the_systems_leaves_nothing_behind(mysql_config, monkeypatch, inject):
    lock_name = _lock_name(mysql_config)
    sector = _sector("Halfway Sector")
    generated = _names(sector)
    with monkeypatch.context() as patch:
        inject(patch, mysql_config, lock_name)
        with pytest.raises(pymysql.err.IntegrityError):
            store.save_sector(sector, config=mysql_config, galaxy_position=_POSITION)

    assert _names(sector) == generated
    _assert_nothing_persisted(mysql_config)
    assert _lock_is_free(mysql_config, lock_name)

    sector_id = store.save_sector(sector, config=mysql_config, galaxy_position=_POSITION)
    assert _query(mysql_config, "SELECT name FROM sectors")[0]["name"] == "Halfway Sector"
    _assert_saved_once(mysql_config, sector_id)
    assert _query(mysql_config, "SELECT COUNT(*) AS n FROM nearest_systems")[0]["n"] > 0
    assert _lock_is_free(mysql_config, lock_name)


def test_a_deadlock_on_every_attempt_gives_up_after_the_last_one(mysql_config, monkeypatch):
    lock_name = _lock_name(mysql_config)
    sector = _sector("Deadlocked Sector")
    generated = _names(sector)
    real_insert = store.insert_sector
    seen, sleeps = [], []

    def insert(conn, sector_arg, galaxy_position=None, **kwargs):
        seen.append(_names(sector_arg))
        return real_insert(conn, sector_arg, galaxy_position=galaxy_position, **kwargs)

    def deadlock(conn, sector_id):
        raise pymysql.err.OperationalError(1213, "Deadlock found when trying to get lock")

    monkeypatch.setattr(store, "insert_sector", insert)
    monkeypatch.setattr(store, "_add_sector_to_nearest", deadlock)
    monkeypatch.setattr(store.time, "sleep", sleeps.append)
    with pytest.raises(pymysql.err.OperationalError) as raised:
        store.save_sector(sector, config=mysql_config, galaxy_position=_POSITION)

    assert raised.value.args[0] == 1213
    assert seen == [generated] * store.SECTOR_SAVE_ATTEMPTS
    assert len(sleeps) == store.SECTOR_SAVE_ATTEMPTS - 1
    assert _names(sector) == generated
    _assert_nothing_persisted(mysql_config)
    assert _lock_is_free(mysql_config, lock_name)


def test_a_lock_wait_timeout_on_the_neighbour_lock_retries_until_it_is_free(mysql_config, monkeypatch):
    neighbor_id = store.save_sector(_sector("Older Sector", ["Ilmaren", "Quorra"], []),
                                  config=mysql_config, galaxy_position=_NEIGHBOR_POSITION)
    lock_name = _lock_name(mysql_config)
    holder = pymysql.connect(**_test_server_kwargs())
    try:
        with holder.cursor() as cur:
            cur.execute("SELECT GET_LOCK(%s, 0)", (lock_name,))
            assert cur.fetchone()[0] == 1

        real_insert = store.insert_sector
        attempts, sleeps = [], []

        def insert(conn, sector_arg, galaxy_position=None, **kwargs):
            attempts.append(_names(sector_arg))
            return real_insert(conn, sector_arg, galaxy_position=galaxy_position, **kwargs)

        def sleep(seconds):
            sleeps.append(seconds)
            with holder.cursor() as cur:
                cur.execute("SELECT RELEASE_LOCK(%s)", (lock_name,))

        # The default timeout is bound when `lock_until_commit` is defined.
        monkeypatch.setattr(store.Connection.lock_until_commit, "__defaults__", (0.2,))
        monkeypatch.setattr(store, "insert_sector", insert)
        monkeypatch.setattr(store.time, "sleep", sleep)
        sector = _sector("Waiting Sector")
        generated = _names(sector)
        sector_id = store.save_sector(sector, config=mysql_config, galaxy_position=_POSITION)
    finally:
        holder.close()

    assert len(sleeps) == 1
    assert attempts == [generated, generated]
    assert _lock_is_free(mysql_config, lock_name)
    systems = _query(mysql_config, "SELECT name FROM star_systems WHERE sector_id = ? ORDER BY name", (sector_id,))
    assert [r["name"] for r in systems] == _SAVED_SYSTEMS
    assert _query(mysql_config, "SELECT COUNT(*) AS n FROM star_systems")[0]["n"] == 5
    assert _query(mysql_config, "SELECT COUNT(*) AS n FROM rogue_planets")[0]["n"] == 2
    assert _orphans(mysql_config) == []
    # The stored neighbour lists are exactly what a full recompute finds,
    # and the older sector's lists picked up the new sector's systems.
    conn = store.get_connection(mysql_config)
    try:
        assert store.refresh_nearest_systems(conn, [neighbor_id, sector_id]) == set()
        conn.rollback()
        linked = conn.execute(
            "SELECT COUNT(*) AS n FROM nearest_systems ns JOIN star_systems s ON s.id = ns.neighbor_system_id "
            "WHERE ns.sector_id = ? AND s.sector_id = ?", (neighbor_id, sector_id)).fetchone()["n"]
        assert linked > 0
    finally:
        conn.close()


def test_the_neighbour_lock_holder_gives_up_a_row_its_waiter_holds(mysql_config):
    """PERF.21: a save waiting for the neighbour lock keeps the rows it
    already wrote (a name collision renames an existing system), and the
    holder's nearest-neighbour rows can need one of them. InnoDB can't see
    the wait on the named lock, so both used to sit until its 50 s
    timeout, every time two workers met like that. The holder now gives
    up within `LOCK_HOLDER_ROW_WAIT_S` (a 1205 that `save_sector`
    retries), and the waiter goes on."""
    import threading
    import time

    setup = store.get_connection(mysql_config)
    try:
        setup.execute("CREATE TABLE lock_probe (id INT PRIMARY KEY, v INT) ENGINE=InnoDB")
        setup.execute("INSERT INTO lock_probe VALUES (1, 0)")
        setup.commit()
    finally:
        setup.close()
    lock_name = "planetgen.test.lock-probe"
    holder = store.get_connection(mysql_config)
    waiter = store.get_connection(mysql_config)
    waited = {}

    def wait_for_lock():
        waiter.execute("UPDATE lock_probe SET v = 1 WHERE id = 1")
        waiter.lock_until_commit(lock_name)
        waiter.commit()
        waited["done"] = True

    try:
        holder.lock_until_commit(lock_name)
        thread = threading.Thread(target=wait_for_lock, daemon=True)
        thread.start()
        time.sleep(0.5)  # the waiter has its row and waits for the lock
        started = time.monotonic()
        with pytest.raises(pymysql.err.OperationalError) as raised:
            holder.execute("UPDATE lock_probe SET v = 2 WHERE id = 1")
        assert raised.value.args[0] == 1205
        assert time.monotonic() - started < store.LOCK_HOLDER_ROW_WAIT_S + 5
        holder.rollback()
        thread.join(timeout=20)
        assert waited.get("done")
        # The pooled connection's own row wait is back to normal.
        row = holder.execute("SELECT @@SESSION.innodb_lock_wait_timeout AS t").fetchone()
        assert int(row["t"]) > store.LOCK_HOLDER_ROW_WAIT_S
    finally:
        holder.close()
        waiter.close()
