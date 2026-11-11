# tests/test_system_scene.py

"""
MAP.69: `GET /api/systems/<id>/scene` and `web/maps/systemscene.build_scene`
-- every star, planet, moon, belt and comet of a system with its ref,
radius, colour, orbit elements and position at the epoch. Needs MySQL
(`mysql_config`).
"""

import json
import pytest

from planetgen.api.config import Config
from planetgen.db import store
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.web.app import create_app
from planetgen.web.maps.systemscene import build_scene
from tests.publicids import pid


def _system(want_binary):
    for _ in range(400):
        cfg = SystemConfig()
        cfg.PLANETS = True
        cfg.MAX_PLANETS = True
        cfg.BINARY_SYSTEM = want_binary
        system = StarSystem(system_config=cfg)
        if want_binary and len(system.stars) < 2:
            continue
        if any(getattr(p, "moons", None) for p in system.planets) and system.comets:
            return system, cfg
    pytest.skip("no suitable system generated")


@pytest.fixture
def saved(mysql_config):
    sector = SpaceSector("Scene Sector", edge_ly=11.5)
    single, cfg_a = _system(False)
    sector.add_system(single, position=(1.0, 1.0, 1.0), system_config=cfg_a)
    try:
        binary, cfg_b = _system(True)
        sector.add_system(binary, position=(-2.0, 0.5, 3.0), system_config=cfg_b)
    except pytest.skip.Exception:
        pass
    store.save_sector(sector, config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        ids = [r["id"] for r in conn.execute("SELECT id FROM star_systems ORDER BY id").fetchall()]
    finally:
        conn.close()
    return mysql_config, ids


def test_the_scene_has_every_body_with_orbit_elements(saved):
    config, ids = saved
    conn = store.get_connection(config)
    try:
        scene = build_scene(conn, ids[0])
        detail_counts = {
            "planets": conn.execute("SELECT COUNT(*) AS n FROM planets WHERE star_system_id = ?", (ids[0],)).fetchone()["n"],
            "moons": conn.execute("SELECT COUNT(*) AS n FROM moons WHERE star_system_id = ?", (ids[0],)).fetchone()["n"],
            "comets": conn.execute("SELECT COUNT(*) AS n FROM comets WHERE star_system_id = ?", (ids[0],)).fetchone()["n"],
        }
    finally:
        conn.close()
    assert scene["system"]["ref"] == f"system:{ids[0]}"
    assert len(scene["planets"]) == detail_counts["planets"]
    assert sum(len(p["moons"]) for p in scene["planets"]) == detail_counts["moons"]
    assert len(scene["comets"]) == detail_counts["comets"]
    star = scene["stars"][0]
    assert star["ref"].startswith("star:") and star["position_km"] == [0.0, 0.0, 0.0] and star["look"]["color"].startswith("#")
    planet = scene["planets"][0]
    assert planet["ref"].startswith("planet:") and planet["orbit"]["around"] == star["ref"]
    assert set(planet["orbit"]) >= {"distance_km", "period_years", "inclination_deg", "ascending_node_deg",
                                    "position_km", "velocity_kms", "epoch_unix"}
    assert "phase_deg" not in planet["orbit"]  # GEN.289: a planet is its state vector, not a phase
    assert planet["orbit"]["position_km"] == planet["position_km"] and len(planet["orbit"]["velocity_kms"]) == 3
    assert len(planet["position_km"]) == 3 and planet["radius_km"] > 0 and planet["look"]["color"].startswith("#")
    moon = next(m for p in scene["planets"] for m in p["moons"])
    assert moon["ref"].startswith("moon:") and moon["orbit"]["around"] == moon["parent"]
    assert star["look"]["body"] == "star" and star["look"]["seed"] and "color" not in star
    assert planet["look"]["body"] in ("solid", "gas") and planet["look"]["seed"] and "color" not in planet
    assert moon["look"]["kind"] == "moon" and moon["look"]["body"] == "solid" and "color" not in moon
    # MAP.343: the fields the map's reasons to look read, as the rows hold them
    for body in (planet, moon):
        assert {"life_microbial", "life_flora", "life_fauna", "life_sentient", "life_stage", "planet_class"} <= set(body), sorted(body)
    json.dumps(scene, allow_nan=False)    # a look record never carries a NaN or an infinity
    comet = scene["comets"][0]
    assert comet["ref"].startswith("comet:") and comet["orbit"]["kepler"]["eccentricity"] >= 0
    assert len(comet["position_km"]) == 3


def test_a_binary_pair_is_two_positioned_stars_with_a_mutual_orbit(saved):
    config, ids = saved
    if len(ids) < 2:
        pytest.skip("no binary generated")
    conn = store.get_connection(config)
    try:
        scene = build_scene(conn, ids[1])
    finally:
        conn.close()
    assert len(scene["stars"]) == 2 and scene["system"]["is_binary"]
    secondary = next(s for s in scene["stars"] if s["role"] == "secondary")
    assert secondary["orbit"]["distance_km"] > 0 and secondary["orbit"]["around"]
    assert secondary["position_km"] != [0.0, 0.0, 0.0]


def test_the_endpoint_returns_the_scene_and_404s_for_a_missing_system(saved):
    config, ids = saved

    class TestConfig(Config):
        MYSQL_CONFIG = config
        WRITE_MYSQL_CONFIG = config
        CONTROL_MYSQL_CONFIG = config

    app = create_app(TestConfig)
    app.testing = True
    client = app.test_client()
    body = client.get(f"/api/systems/{pid('system', ids[0])}/scene").get_json()
    assert body["system"]["id"] == pid("system", ids[0]) and body["planets"] and "epoch" in body
    assert client.get("/api/systems/999999999/scene").status_code == 404


def test_a_belt_carries_its_heaviest_bodies_with_their_looks_and_places(mysql_config):
    """MAP.414: each belt of the scene lists the heaviest bodies it stores (at most `BELT_MEMBERS_DRAWN`) by rank, each with the designation the body's page prints, its radius, a
    `member_look` whose seed is the designation's, and where its stored orbit puts it, a distance from the belt's star inside the orbit's reach."""
    import math

    from planetgen.galaxy import object_uid
    from planetgen.web.maps import bodylook, systemscene
    from tests.fuzz_support import first_seed_with

    def build(seed):
        cfg = SystemConfig()
        cfg.ASTEROID_BELT = True
        return StarSystem(system_config=cfg), cfg

    _seed, (system, cfg) = first_seed_with(build, lambda built: any(getattr(b, "members", None) for b in built[0].planets if b.body_type == "a"), what="a belt that stores bodies")
    sector = SpaceSector("Belt Scene Sector", edge_ly=11.5)
    sector.add_system(system, position=(1.0, 1.0, 1.0), system_config=cfg)
    store.save_sector(sector, config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        system_id = conn.execute("SELECT id FROM star_systems LIMIT 1").fetchone()["id"]
        scene = build_scene(conn, system_id)
        stored = {row["id"]: row["uid"] for row in conn.execute("SELECT id, uid FROM asteroid_belts WHERE star_system_id = ?", (system_id,)).fetchall()}
        kept = {row["id"]: conn.execute("SELECT COUNT(*) AS n FROM belt_members WHERE belt_id = ? AND status <> 'moon'", (row["id"],)).fetchone()["n"] for row in scene["belts"]}
    finally:
        conn.close()
    drawn = [belt for belt in scene["belts"] if belt["members"]]
    assert drawn
    for belt in scene["belts"]:
        assert len(belt["members"]) == min(kept[belt["id"]], systemscene.BELT_MEMBERS_DRAWN)
        assert [m["rank"] for m in belt["members"]] == sorted(m["rank"] for m in belt["members"])
        for member in belt["members"]:
            assert member["designation"] == object_uid.designation(object_uid.from_bytes(stored[belt["id"]]), member["rank"])
            assert member["radius_km"] > 0.0 and member["look"]["kind"] == "belt_member"
            assert member["look"]["seed"] == bodylook.seed_for("belt_member", member["designation"])
            assert 0.0 < math.dist(member["position_km"], (0.0, 0.0, 0.0)) < 10.0 * belt["outer_km"]


class _Rows:
    def __init__(self, rows):
        self.rows = rows

    def fetchall(self):
        return self.rows


class _MembersOnly:
    """A connection that answers the one query `_belt_members` asks, with the rows given and the SQL kept."""

    def __init__(self, rows):
        self.rows, self.asked = rows, []

    def execute(self, sql, params=()):
        self.asked.append((sql, params))
        return _Rows(self.rows)


def _stored_row(rank, **changes):
    row = {"member_rank": rank, "mass_kg": 1.0e20 / rank, "radius_km": 200.0 / rank, "distance_km": 3.0e8 + 1.0e6 * rank, "eccentricity": 0.1, "inclination_deg": 5.0 * rank,
           "ascending_node_deg": 40.0, "arg_periapsis_deg": 70.0, "mean_anomaly_deg": 30.0 * rank, "taxon": "C", "density_g_cm3": 1.5, "rock_fraction": 0.7, "metal_fraction": 0.1,
           "water_fraction": 0.2, "ice_fraction": 0.0, "status": "minor body"}
    row.update(changes)
    return row


def test_the_members_of_a_belt_are_the_heaviest_stored_ones_placed_by_their_orbits():
    """MAP.414 without a database: the query asks for the belt's bodies by rank, never the ones that became moons, up to the cap; the designation is the belt ID's and the rank's, the
    look is the body page's, and the position is on the stored orbit (distance from the star between periapsis and apoapsis, the plane tilted by the inclination)."""
    import math

    from planetgen.galaxy import object_uid
    from planetgen.physics import constants as pc
    from planetgen.web.maps import bodylook, systemscene

    belt_id = object_uid.parse_id("0000FE81000A2B-0000005-003")
    uid = object_uid.to_bytes(belt_id)
    rows = [_stored_row(rank) for rank in (1, 2, 3)]
    conn = _MembersOnly(rows)
    members = systemscene._belt_members(conn, {"id": 7, "uid": uid}, limit=2)
    sql, params = conn.asked[0]
    assert "ORDER BY member_rank" in sql and "LIMIT 2" in sql and "status <> 'moon'" in sql and params == (7,)
    assert len(members) == 3      # the stub returns what it holds; the LIMIT is the database's
    for member, row in zip(members, rows):
        rank = row["member_rank"]
        assert member["rank"] == rank and member["designation"] == object_uid.designation(belt_id, rank)
        assert member["radius_km"] == row["radius_km"]
        assert member["look"]["kind"] == "belt_member" and member["look"]["seed"] == bodylook.seed_for("belt_member", member["designation"])
        distance = math.dist(member["position_km"], (0.0, 0.0, 0.0))
        a_km = row["distance_km"]
        assert a_km * (1.0 - row["eccentricity"]) * (1.0 - 1e-9) <= distance <= a_km * (1.0 + row["eccentricity"]) * (1.0 + 1e-9)
        sine = member["position_km"][2] / distance
        assert abs(sine) <= math.sin(math.radians(row["inclination_deg"])) * (1.0 + 1e-9)
    assert systemscene._belt_members(_MembersOnly(rows), {"id": 7, "uid": uid}, limit=0) == []
    assert pc.AU_TO_KM > 1.0e8


def test_a_belt_with_no_stored_id_still_names_its_bodies():
    from planetgen.web.maps import systemscene

    members = systemscene._belt_members(_MembersOnly([_stored_row(1)]), {"id": 9, "uid": None})
    assert members[0]["designation"] == "9.0001"
