# tests/test_system_scene.py

"""
MAP.69: `GET /api/systems/<id>/scene` and `web/maps/systemscene.build_scene`
-- every star, planet, moon, belt and comet of a system with its ref,
radius, colour, orbit elements and position at the epoch. Needs MySQL
(`mysql_config`).
"""

import pytest

from planetgen.api.config import Config
from planetgen.db import store
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.web.app import create_app
from planetgen.web.maps.systemscene import build_scene


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
    assert star["ref"].startswith("star:") and star["position_km"] == [0.0, 0.0, 0.0] and star["color"].startswith("#")
    planet = scene["planets"][0]
    assert planet["ref"].startswith("planet:") and planet["orbit"]["around"] == star["ref"]
    assert set(planet["orbit"]) >= {"distance_km", "period_years", "inclination_deg", "ascending_node_deg", "phase_deg"}
    assert len(planet["position_km"]) == 3 and planet["radius_km"] > 0 and planet["color"].startswith("#")
    moon = next(m for p in scene["planets"] for m in p["moons"])
    assert moon["ref"].startswith("moon:") and moon["orbit"]["around"] == moon["parent"]
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
    body = client.get(f"/api/systems/{ids[0]}/scene").get_json()
    assert body["system"]["id"] == ids[0] and body["planets"] and "epoch" in body
    assert client.get("/api/systems/999999999/scene").status_code == 404
