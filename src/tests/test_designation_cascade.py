# tests/test_designation_cascade.py

"""
GEN.1141: a system that changes its ID takes its bodies with it (each body keeps its number under the new ID and gets its own history row), a belt
member's or captured body's designation follows its parent without a row, a delete writes tombstones so every old link answers 410, and the history
of a body follows it when its system's content is written again.
"""

import pytest

from planetgen.api import ids
from planetgen.db import designation_history, store
from planetgen.galaxy import designations, object_uid
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.util import draw
from tests.publicids import pid
from tests.test_admin_edits import _saved_nebula, _saved_system, admin, client  # noqa: F401 -- fixtures
from tests.test_designation_history import PLACE, _planetary_system, stored  # noqa: F401 -- a fixture

pytestmark = pytest.mark.db

BODY_TABLES = designation_history.BODY_TABLES


def _bodies(conn, system_id):
    return {(table, row["id"]): object_uid.from_bytes(row["uid"]) for table in BODY_TABLES for row in conn.execute(
        f"SELECT id, uid FROM {table} WHERE star_system_id = ?", (system_id,)).fetchall()}


def _free_system_uid(conn, near):
    for serial_step in range(1, 100):
        candidate = object_uid.with_body(near + (serial_step << object_uid.LAYOUT.body_bits), 0)
        if not designation_history.is_live_uid(conn, candidate) and not designation_history.is_retired_uid(conn, candidate):
            return candidate
    pytest.fail("no free system ID")


def test_a_systems_bodies_follow_its_new_id(stored):
    conn, system_id = stored
    before = _bodies(conn, system_id)
    old_system = object_uid.from_bytes(conn.execute("SELECT uid FROM star_systems WHERE id = ?", (system_id,)).fetchone()["uid"])
    new_system = _free_system_uid(conn, old_system)
    with conn:
        designation_history.renumber(conn, "star_systems", system_id, new_system)
    after = _bodies(conn, system_id)
    assert before.keys() == after.keys() and len(before) > 2
    for key, old in before.items():
        assert after[key] == object_uid.with_body(new_system, object_uid.body_of(old))
        assert designation_history.is_retired_uid(conn, old) and not designation_history.is_live_uid(conn, old)
        rows = designation_history.history(conn, key[0], key[1])
        assert [(r.designation, r.reason) for r in rows] == [(old, "parent change")]
    assert [(r.designation, r.reason) for r in designation_history.history(conn, "star_systems", system_id)] == [(old_system, "admin renumber")]
    with pytest.raises(designations.DesignationError):
        designation_history.renumber(conn, "star_systems", system_id, old_system)


def test_a_cascade_that_cannot_finish_changes_nothing(stored):
    conn, system_id = stored
    old_system = object_uid.from_bytes(conn.execute("SELECT uid FROM star_systems WHERE id = ?", (system_id,)).fetchone()["uid"])
    bodies = _bodies(conn, system_id)
    taken = next(iter(bodies.values()))
    new_system = _free_system_uid(conn, old_system)
    clash = object_uid.with_body(new_system, object_uid.body_of(taken))
    with conn:
        conn.execute("UPDATE designations SET id = id WHERE 1 = 0")
        conn.execute("INSERT INTO designations (designation_kind, uid, object_kind, object_row_id, retired_epoch, reason)"
                     " VALUES ('uid', ?, 'moons', 999999, 1, 'rename')", (object_uid.to_bytes(clash),))
    with pytest.raises(designations.DesignationError, match="reserved for ever"):
        designation_history.renumber(conn, "star_systems", system_id, new_system)
    assert _bodies(conn, system_id) == bodies


def test_a_body_takes_only_an_id_of_its_own_system(stored):
    conn, system_id = stored
    moon = conn.execute("SELECT id, uid FROM moons WHERE star_system_id = ? LIMIT 1", (system_id,)).fetchone()
    other = _free_system_uid(conn, object_uid.from_bytes(moon["uid"]))
    with pytest.raises(designations.DesignationError, match="not an ID of the system"):
        designation_history.renumber(conn, "moons", moon["id"], object_uid.with_body(other, 5))


def _belted_system(mysql_config):
    cfg = SystemConfig()
    cfg.PLANETS = True
    cfg.ASTEROID_BELT = True
    cfg.BINARY_SYSTEM = False
    for seed in range(200):
        with draw.bound(seed):
            system = StarSystem(system_config=cfg)
        if any(p.body_type == "a" for p in system.planets):
            sector = SpaceSector("Belt Sector", edge_ly=11.5)
            sector.add_system(system, position=(0.0, 0.0, 0.0), system_config=cfg)
            store.save_sector(sector, config=mysql_config, galaxy_position=PLACE)
            return
    pytest.fail("no belted system in 200 seeds")


def test_a_belt_members_old_designation_names_the_member_under_the_new_belt_id(mysql_config):
    _belted_system(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        belt = conn.execute("SELECT id, uid FROM asteroid_belts LIMIT 1").fetchone()
        old = object_uid.from_bytes(belt["uid"])
        free = next(object_uid.with_body(old, body) for body in range(2000, 4000)
                    if not designation_history.is_live_uid(conn, object_uid.with_body(old, body))
                    and not designation_history.is_retired_uid(conn, object_uid.with_body(old, body)))
        with conn:
            designation_history.renumber(conn, "asteroid_belts", belt["id"], free)
        found, current = designation_history.resolve_member(conn, object_uid.designation(old, 3))
        assert found.status == designations.REDIRECT and current == object_uid.designation(free, 3)
        live, same = designation_history.resolve_member(conn, object_uid.designation(free, 3))
        assert live.status == designations.CURRENT and same == object_uid.designation(free, 3)
        gone_found, gone_text = designation_history.resolve_member(conn, object_uid.designation(old + (1 << 40), 3))
        assert gone_found.status == designations.UNKNOWN and gone_text is None
    finally:
        conn.close()


def test_a_deleted_systems_ids_answer_410(admin, mysql_config):
    sector_id, system_id = _saved_system(mysql_config)
    printed = pid("system", system_id, mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        before = _bodies(conn, system_id)
        system_uid = object_uid.from_bytes(conn.execute("SELECT uid FROM star_systems WHERE id = ?", (system_id,)).fetchone()["uid"])
    finally:
        conn.close()
    assert admin.delete(f"/api/systems/{printed}").status_code == 200
    response = admin.get(f"/api/systems/{printed}")
    assert response.status_code == 410 and response.get_json()["gone"] is True
    assert admin.get(f"/system/{printed}").status_code == 410
    conn = store.get_connection(mysql_config)
    try:
        for old in [system_uid, *before.values()]:
            assert designation_history.is_retired_uid(conn, old)
        reasons = {row["reason"] for row in conn.execute("SELECT reason FROM designations").fetchall()}
        assert reasons == {"deleted"}
        assert conn.execute("SELECT COUNT(*) AS n FROM designations WHERE gone_epoch IS NULL").fetchone()["n"] == 0
    finally:
        conn.close()
    assert admin.get(f"/api/systems/{pid('system', 999999, mysql_config) if False else 'FFFFFFFFFFFF-FFFFFFF-FFF'}").status_code == 404


def test_a_deleted_phenomenon_and_a_deleted_sector_answer_410(admin, mysql_config):
    _sector_id, nebula_id = _saved_nebula(mysql_config)
    printed = pid("nebula", nebula_id, mysql_config)
    assert admin.delete(f"/api/phenomena/nebula/{printed}").status_code == 200
    assert admin.get(f"/api/phenomena/nebula/{printed}").status_code == 410
    sector_id, _system_id = _saved_system(mysql_config)
    system_printed = pid("system", _system_id, mysql_config)
    assert admin.delete(f"/api/sectors/{pid('sector', sector_id, mysql_config)}/contents").status_code == 200
    assert admin.get(f"/api/systems/{system_printed}").status_code == 410


def test_a_names_history_follows_its_body_when_the_system_is_written_again(mysql_config):
    system, cfg = _planetary_system()
    sector = SpaceSector("Carry Sector", edge_ly=11.5)
    sector.add_system(system, position=(0.0, 0.0, 0.0), system_config=cfg)
    store.save_sector(sector, config=mysql_config, galaxy_position=PLACE)
    conn = store.get_connection(mysql_config)
    try:
        system_id = conn.execute("SELECT id FROM star_systems").fetchone()["id"]
        planet = conn.execute("SELECT id, name, uid FROM planets WHERE star_system_id = ? ORDER BY id LIMIT 1", (system_id,)).fetchone()
        with conn:
            store.rename_body(conn, "planets", planet["id"], "Hand-named")
        fresh, fresh_cfg = _planetary_system()
        with conn:
            assert store.replace_star_system_content(conn, system_id, fresh, fresh_cfg)
        now = conn.execute("SELECT id, name FROM planets WHERE uid = ?", (planet["uid"],)).fetchone()
        assert now is not None and now["id"] != planet["id"]
        assert [(r.designation, r.reason) for r in designation_history.history(conn, "planets", now["id"])] == [(planet["name"], "rename")]
        assert designation_history.history(conn, "planets", planet["id"]) == []
        # no stand-in name of the rebuilt system was kept as an object's
        stand_ins = conn.execute("SELECT COUNT(*) AS n FROM designations WHERE designation_kind = 'name' AND reason = 'parent change'").fetchone()["n"]
        assert stand_ins == 0
    finally:
        conn.close()
