# tests/test_designation_tombstones.py

"""GEN.1141 (part e, tombstones): deleting an object keeps its ID (reserved for ever, `designations.reason` `deleted`) and its name, and stamps `gone_epoch`
on every row that names it, so a link to it answers 410 and the ID is never issued again. A system's deletion takes its bodies' IDs with it; an edit that drops a planet drops its moons' IDs too."""

import pytest

from planetgen.db import designation_rows as rows
from planetgen.db import edits, store
from planetgen.galaxy import designations as d
from planetgen.galaxy import object_uid as ou
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.rogue import RoguePlanet
from planetgen.generation.system import StarSystem
from planetgen.util import draw
from tests.test_api import admin_client, client, default_admin_client, first_admin_password  # noqa: F401

pytestmark = pytest.mark.fresh_database


@pytest.fixture
def galaxy(mysql_config, admin_client):  # noqa: F811
    sector = SpaceSector("Tombstone Sector", edge_ly=20.0)
    config = SystemConfig()
    config.PLANETS = True
    config.MAX_PLANETS = True
    config.BINARY_SYSTEM = False
    made = 0
    for seed in range(80):      # two systems, the first with a moon
        with draw.bound(seed):
            system = StarSystem(system_config=config)
        if made == 0 and not any(getattr(p, "moons", None) for p in system.planets):
            continue
        sector.add_system(system, position=(made + 1.0, 1.0, 1.0), system_config=config)
        made += 1
        if made == 2:
            break
    assert made == 2
    sector.add_phenomenon(RoguePlanet(SystemConfig()), "rogue-planet")
    store.save_sector(sector, config=mysql_config)
    return mysql_config, admin_client


def _bodies(conn, system_id):
    found = []
    for table in ("stars", "planets", "moons", "asteroid_belts", "comets"):
        for row in conn.execute(f"SELECT id, uid FROM {table} WHERE star_system_id = ? AND uid IS NOT NULL", (system_id,)).fetchall():
            found.append((table, row["id"], ou.from_bytes(row["uid"])))
    return found


def test_a_deleted_system_keeps_its_id_and_each_bodys_and_answers_410(galaxy):
    config, client = galaxy
    conn = store.get_connection(config)
    try:
        system = conn.execute("SELECT id, uid, name FROM star_systems ORDER BY id LIMIT 1").fetchone()
        other = conn.execute("SELECT id FROM star_systems ORDER BY id DESC LIMIT 1").fetchone()
        old = ou.from_bytes(system["uid"])
        bodies = _bodies(conn, system["id"])
        assert len(bodies) > 2
    finally:
        conn.close()
    response = client.delete(f"/api/systems/{ou.format_id(old)}")
    assert response.status_code == 200, response.get_json()
    conn = store.get_connection(config)
    try:
        assert conn.execute("SELECT 1 FROM star_systems WHERE id = ?", (system["id"],)).fetchone() is None
        assert conn.execute("SELECT 1 FROM star_systems WHERE id = ?", (other["id"],)).fetchone() is not None
        for table, row_id, uid in [("star_systems", system["id"], old), *bodies]:
            kept = rows.retired_uid(conn, uid)
            assert kept is not None and kept.reason == "deleted" and kept.gone_epoch is not None
            assert (kept.object_kind, kept.object_row_id) == (table, row_id)
            assert rows.live_holder(conn, uid) is None
        kept_name = rows.retired_names(conn, system["name"])
        assert [(r.object_kind, r.reason) for r in kept_name if r.object_kind == "star_systems"] == [("star_systems", "deleted")]
        assert all(r.gone_epoch is not None for r in kept_name if r.object_kind == "star_systems")
        assert conn.execute("SELECT COUNT(*) AS n FROM designations WHERE object_kind = 'star_systems' AND object_row_id = ? AND gone_epoch IS NULL", (other["id"],)).fetchone()["n"] == 0
    finally:
        conn.close()
    answer = client.get(f"/api/systems/{ou.format_id(old)}")
    assert answer.status_code == 410 and answer.get_json()["gone_epoch"] is not None
    assert client.get(f"/system/{ou.format_id(old)}").status_code == 410


def test_a_deleted_id_cannot_be_given_to_another_object(galaxy):
    config, client = galaxy
    conn = store.get_connection(config)
    try:
        system = conn.execute("SELECT id, uid FROM star_systems ORDER BY id LIMIT 1").fetchone()
        old = ou.from_bytes(system["uid"])
        other = conn.execute("SELECT id FROM star_systems ORDER BY id DESC LIMIT 1").fetchone()
    finally:
        conn.close()
    assert client.delete(f"/api/systems/{ou.format_id(old)}").status_code == 200
    conn = store.get_connection(config)
    try:
        with pytest.raises(d.DesignationError):
            with conn:
                rows.renumber(conn, "star_systems", other["id"], old, "admin renumber")
    finally:
        conn.close()


def test_a_deleted_phenomenon_keeps_its_id_and_answers_410(galaxy):
    config, client = galaxy
    conn = store.get_connection(config)
    try:
        rogue = conn.execute("SELECT id, uid FROM rogue_planets ORDER BY id LIMIT 1").fetchone()
        old = ou.from_bytes(rogue["uid"])
        with conn:
            assert edits.delete_phenomenon(conn, "rogue_planet", rogue["id"])
        kept = rows.retired_uid(conn, old)
        assert kept is not None and kept.reason == "deleted" and kept.gone_epoch is not None and kept.object_kind == "rogue_planets"
    finally:
        conn.close()


def test_an_edit_that_drops_a_planet_keeps_the_ids_of_the_planet_and_its_moons(galaxy):
    config, _client = galaxy
    conn = store.get_connection(config)
    try:
        system_id = conn.execute("SELECT s.id FROM star_systems s JOIN moons m ON m.star_system_id = s.id ORDER BY s.id LIMIT 1").fetchone()["id"]
        system = store.load_star_system(conn, system_id)
        victim = next(p for p in system.planets if getattr(p, "moons", None))
        planet_uid = ou.from_bytes(conn.execute("SELECT uid FROM planets WHERE id = ?", (victim.db_id,)).fetchone()["uid"])
        moon_uids = [ou.from_bytes(conn.execute("SELECT uid FROM moons WHERE id = ?", (moon.db_id,)).fetchone()["uid"]) for moon in victim.moons]
        system.planets = [p for p in system.planets if p is not victim]
        with conn:
            edits.save_system_edits(conn, system_id, system)
        assert conn.execute("SELECT 1 FROM planets WHERE id = ?", (victim.db_id,)).fetchone() is None
        for uid in [planet_uid, *moon_uids]:
            kept = rows.retired_uid(conn, uid)
            assert kept is not None and kept.reason == "deleted" and kept.gone_epoch is not None
        # the planets that stayed keep what they had
        for planet in system.planets:
            assert not rows.history(conn, "planets", planet.db_id)
    finally:
        conn.close()
