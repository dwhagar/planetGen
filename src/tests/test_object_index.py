# tests/test_object_index.py

"""PERF.140: the object index builder (`db/objectindex.py`): it matches the live tables after a build, a refresh writes only what
changed, a gone object loses its row and `check` says when the index and the tables differ."""

import pytest

from planetgen.db import objectindex, store
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.nebula import Nebula
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.system import StarSystem
from tests.publicids import pid

pytestmark = pytest.mark.fresh_database


@pytest.fixture
def galaxy(mysql_config):
    """A sector of three systems with their planets, moons, belts and comets, a nebula and a facility-less galaxy."""
    sector = SpaceSector("Index Sector", edge_ly=10.0)
    for index, star_type in enumerate(("G2V", "K5V", "M5V")):
        config = SystemConfig()
        config.STAR_TYPE = star_type
        config.PLANETS = True
        config.BINARY_SYSTEM = False
        sector.add_system(StarSystem(system_config=config), position=(index + 1.0, 1.0, 1.0), system_config=config)
    sector_id = store.save_sector(sector, config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        nebula_id = store.insert_nebula(conn, Nebula(SystemConfig(), name="Indexed Cloud"), placement={
            "center_x_pc": 5.0, "center_y_pc": 1.0, "center_z_pc": 2.0, "galactic_radius_pc": 5.5})
        conn.commit()
        pid("nebula", nebula_id, conn)   # issues the object ID a nebula is indexed by
    finally:
        conn.close()
    return mysql_config, sector_id


def _conn(config):
    return store.get_connection(config)


def _rows(conn, kind=None):
    where = " WHERE kind = %s" if kind else ""
    return conn.execute(f"SELECT * FROM object_index{where} ORDER BY kind, uid", (kind,) if kind else ()).fetchall()


def _assert_matches(conn):
    report = objectindex.check(conn)
    for kind, found in report.items():
        assert (found["missing"], found["extra"], found["different"]) == (0, 0, 0), (kind, found)
        assert found["live"] == found["indexed"], (kind, found)


def test_a_build_indexes_every_object_and_nothing_differs(galaxy):
    config, _sector = galaxy
    conn = _conn(config)
    try:
        assert _rows(conn) == []
        built = objectindex.build(conn)
        _assert_matches(conn)
        assert built["system"]["rows"] == 3 and built["nebula"]["rows"] == 1 and built["star"]["rows"] >= 3
        counts = {kind: conn.execute(f"SELECT COUNT(*) AS n FROM {objectindex.SOURCES[kind]['table']}").fetchone()["n"]
                  for kind in ("star", "planet", "moon", "belt", "comet")}
        assert {kind: built[kind]["rows"] for kind in counts} == counts
        system = _rows(conn, "system")[0]
        assert system["sector_id"] is not None and system["placed"] in (0, 1) and system["object_class"]
        assert _rows(conn, "nebula")[0]["name"] == "Indexed Cloud"
    finally:
        conn.close()


def test_a_second_build_changes_no_row(galaxy):
    config, _sector = galaxy
    conn = _conn(config)
    try:
        objectindex.build(conn)
        before = _rows(conn)
        again = objectindex.build(conn)
        assert _rows(conn) == before
        assert all(found["changed"] == 0 and found["removed"] == 0 for found in again.values())
    finally:
        conn.close()


def test_a_refresh_after_one_edit_rewrites_only_that_objects_row(galaxy):
    config, _sector = galaxy
    conn = _conn(config)
    try:
        objectindex.build(conn, database=config.database)
        before = {(row["kind"], bytes(row["uid"])): row for row in _rows(conn)}
        target = conn.execute("SELECT id, uid FROM star_systems ORDER BY id LIMIT 1").fetchone()
        conn.execute("UPDATE star_systems SET name = 'Renamed Primary' WHERE id = %s", (target["id"],))
        conn.commit()
        refreshed = objectindex.refresh(conn, database=config.database)
        after = {(row["kind"], bytes(row["uid"])): row for row in _rows(conn)}
        changed = [key for key in after if after[key] != before.get(key)]
        assert changed == [("system", bytes(target["uid"]))]
        assert after[changed[0]]["name"] == "Renamed Primary"
        assert all(found["mode"] == "refresh" for found in refreshed.values())
        _assert_matches(conn)
    finally:
        conn.close()
        objectindex.forget(config.database)


def test_an_edit_to_a_planet_reaches_its_row_through_its_system(galaxy):
    config, _sector = galaxy
    conn = _conn(config)
    try:
        planet = conn.execute("SELECT id, star_system_id FROM planets ORDER BY id LIMIT 1").fetchone()
        if planet is None:
            pytest.skip("this draw made no planet")
        objectindex.build(conn, database=config.database)
        conn.execute("UPDATE planets SET name = 'Edited World' WHERE id = %s", (planet["id"],))
        store.touch_star_system(conn, planet["star_system_id"])
        conn.commit()
        objectindex.refresh(conn, database=config.database)
        assert [row["name"] for row in _rows(conn, "planet") if row["source_id"] == planet["id"]] == ["Edited World"]
        _assert_matches(conn)
    finally:
        conn.close()
        objectindex.forget(config.database)


def test_a_deleted_object_loses_its_row_and_a_new_one_gains_one(galaxy):
    config, _sector = galaxy
    conn = _conn(config)
    try:
        objectindex.build(conn, database=config.database)
        gone = conn.execute("SELECT id FROM nebulae ORDER BY id LIMIT 1").fetchone()["id"]
        conn.execute("DELETE FROM nebulae WHERE id = %s", (gone,))
        conn.commit()
        objectindex.refresh(conn, ["nebula"], database=config.database)
        assert _rows(conn, "nebula") == []
        second = store.insert_nebula(conn, Nebula(SystemConfig(), name="Second Cloud"), placement={
            "center_x_pc": 6.0, "center_y_pc": 1.0, "center_z_pc": 2.0, "galactic_radius_pc": 6.5})
        conn.commit()
        pid("nebula", second, conn)
        objectindex.refresh(conn, ["nebula"], database=config.database)
        assert [row["name"] for row in _rows(conn, "nebula")] == ["Second Cloud"]
        _assert_matches(conn)
    finally:
        conn.close()
        objectindex.forget(config.database)


def test_a_kind_with_no_record_of_a_run_is_built(galaxy):
    config, _sector = galaxy
    conn = _conn(config)
    try:
        objectindex.forget(config.database)
        result = objectindex.refresh(conn, ["system", "nebula"], database=config.database)
        assert {kind: found["mode"] for kind, found in result.items()} == {"system": "build", "nebula": "build"}
        again = objectindex.refresh(conn, ["system", "nebula"], database=config.database)
        assert {found["mode"] for found in again.values()} == {"refresh"}
    finally:
        conn.close()
        objectindex.forget(config.database)


def test_check_says_what_is_missing_extra_and_different(galaxy):
    config, _sector = galaxy
    conn = _conn(config)
    try:
        objectindex.build(conn, ["system"])
        conn.execute("DELETE FROM object_index WHERE kind = 'system' LIMIT 1")
        conn.execute("UPDATE object_index SET name = 'Wrong' WHERE kind = 'system' LIMIT 1")
        conn.execute("INSERT INTO object_index (kind, uid, source_id, name) VALUES ('system', %s, 99, 'Ghost')",
                     (b"\x80" * 12,))
        conn.commit()
        found = objectindex.check(conn, ["system"])["system"]
        assert (found["missing"], found["extra"], found["different"]) == (1, 1, 1)
        objectindex.build(conn, ["system"])
        assert objectindex.check(conn, ["system"])["system"]["extra"] == 0
        _assert_matches_kind = objectindex.check(conn, ["system"])["system"]
        assert (_assert_matches_kind["missing"], _assert_matches_kind["different"]) == (0, 0)
    finally:
        conn.close()


def test_an_unknown_kind_is_refused(galaxy):
    config, _sector = galaxy
    conn = _conn(config)
    try:
        with pytest.raises(objectindex.UnknownKind, match="no such kind of object 'planetoid'"):
            objectindex.build(conn, ["planetoid"])
    finally:
        conn.close()


def test_the_status_gives_the_size_and_each_kinds_rows(galaxy):
    config, _sector = galaxy
    conn = _conn(config)
    try:
        objectindex.build(conn, ["system"], database=config.database)
        status = objectindex.status(conn, config.database)
        assert status["size_bytes"] > 0
        by_kind = {found["kind"]: found for found in status["kinds"]}
        assert by_kind["system"]["rows"] == 3 and by_kind["system"]["mode"] == "build" and by_kind["system"]["built_at"]
        assert by_kind["nebula"]["rows"] == 0 and by_kind["nebula"]["mode"] is None
    finally:
        conn.close()
        objectindex.forget(config.database)


# --- the warm-up command, the check and the shared budget -----------------------------------------------------

@pytest.fixture
def indexed_database(galaxy, redis_server):
    config, _sector = galaxy
    objectindex.forget(config.database)
    yield config
    objectindex.forget(config.database)


def _unit(config, mode=None):
    from planetgen.web import warmup
    args = {"database": config.database, "cap": None}
    if mode:
        args["index_mode"] = mode
    return warmup.warm_object_index_unit(args, None)


def test_the_warm_command_builds_the_first_time_and_then_writes_only_what_changed(indexed_database):
    first = _unit(indexed_database)
    assert first["family"] == "object-index" and first["built"] == first["warmed"] > 0 and first["cached"] == 0
    second = _unit(indexed_database)
    assert second["built"] == 0 and second["cached"] == second["warmed"] == first["warmed"]
    assert "0 rows written" in second["line"] and second["line"].startswith("object-index: ")
    conn = _conn(indexed_database)
    try:
        _assert_matches(conn)
    finally:
        conn.close()


def test_the_build_mode_writes_every_row_again_and_a_refresh_after_it_finds_all_current(indexed_database):
    _unit(indexed_database)
    conn = _conn(indexed_database)
    try:
        conn.execute("DELETE FROM object_index WHERE kind = 'planet'")
        conn.commit()
    finally:
        conn.close()
    rebuilt = _unit(indexed_database, "build")
    assert rebuilt["built"] > 0
    assert _unit(indexed_database)["built"] == 0


def test_the_index_is_a_phase_of_the_opening_and_the_changes_that_move_objects_build_it_whole():
    from planetgen.queue import warmstart
    from planetgen.web import warmup
    assert warmup.phase_of("object-index") == "opening" and "object-index" in warmup.FAMILIES
    assert set(warmstart.FULL_INDEX_TAGS) == {"positions", "epoch"}


def test_the_index_takes_its_size_off_the_tile_cache_budget(indexed_database, monkeypatch):
    from planetgen.web.lib import tilecache
    budget = tilecache.budget_bytes()
    assert tilecache.max_cache_bytes() == budget
    _unit(indexed_database)
    reserved = objectindex.reserved_bytes()
    assert reserved > 0 and objectindex.recorded_size(indexed_database.database) == reserved
    assert tilecache.max_cache_bytes() == max(1, budget - reserved)
    monkeypatch.setattr(objectindex, "reserved_bytes", lambda: budget * 2)
    assert tilecache.max_cache_bytes() == 1, "a budget the index fills leaves the tile cache full, not off"


def test_the_check_command_warns_before_the_first_build_and_fails_on_a_differing_row(indexed_database):
    from planetgen.db import check
    from planetgen.queue import commands
    from tests.test_db_check import _problem_text, _run, _status
    assert "object-index" not in commands.DEFAULT_CHECKS
    report = _run(indexed_database, ("object-index",))
    assert _status(report) == {"object index": check.WARN}
    assert "holds no system yet" in _problem_text(report, "object index")
    _unit(indexed_database)
    assert _status(_run(indexed_database, ("object-index",))) == {"object index": check.PASS}
    conn = _conn(indexed_database)
    try:
        conn.execute("UPDATE object_index SET name = 'wrong' WHERE kind = 'system' LIMIT 1")
        conn.commit()
    finally:
        conn.close()
    report = _run(indexed_database, ("object-index",))
    assert _status(report) == {"object index": check.FAIL}
    assert "1 different" in _problem_text(report, "object index")
