# tests/test_nav_system.py

"""
NAV.27: `GET /api/nav?adjusted=1` inside a star system -- two ends in one system are a `within` leg, a course between systems has an `out` and
an `into` leg; each is planned again around the keep-out spheres of the system's bodies (`api/nav_system.adjust_legs`). The scene is the saved
system's, with its first two planets put on either side of the star so that the star is in the way.
"""

import pytest

from planetgen import tuning
from planetgen.api import nav_system
from planetgen.db import obstacles, store
from planetgen.galaxy.sector import SpaceSector
from planetgen.physics.body_positions import positions_at
from planetgen.physics.constants import AU_TO_KM as AU_KM
from tests.publicids import pid
from tests.test_api import client  # noqa: F401  (the fixture)
from tests.test_body_positions import state_orbit
from tests.test_system_scene import _system
from planetgen.web.maps import systemscene


@pytest.fixture
def saved(mysql_config, monkeypatch):
    sector = SpaceSector("Nav System Sector", edge_ly=11.5)
    system, cfg = _system(False)
    sector.add_system(system, position=(1.0, 1.0, 1.0), system_config=cfg)
    store.save_sector(sector, config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        system_id = conn.execute("SELECT id FROM star_systems ORDER BY id").fetchone()["id"]
        planets = [r["id"] for r in conn.execute("SELECT id FROM planets WHERE star_system_id = ? ORDER BY id", (system_id,)).fetchall()]
        star_id = conn.execute("SELECT id FROM stars WHERE star_system_id = ? ORDER BY id", (system_id,)).fetchone()["id"]
    finally:
        conn.close()
    if len(planets) < 2:
        pytest.skip("the generated system has fewer than two planets")
    real_build = systemscene.build_scene

    def scene_with_the_star_between(conn, sid):
        scene = real_build(conn, sid)
        around = f"star:{star_id}"
        scene["epoch_unix"] = 1_760_000_000.0
        for planet, phase in zip(scene["planets"][:2], (0.0, 180.0)):
            planet["orbit"] = state_orbit(around, 1.0, 1.0, 0.0, 0.0, phase, scene["epoch_unix"])
        return scene

    monkeypatch.setattr(nav_system, "build_scene", scene_with_the_star_between)
    monkeypatch.setattr(obstacles, "system_keep_outs", lambda conn, scene: {
        f"star:{star_id}": 0.2 * AU_KM, **{p["ref"]: 1.0e5 for p in scene["planets"]}})
    return system_id, planets, star_id


def test_two_planets_either_side_of_the_star_are_one_leg_bent_round_it(saved, client):  # noqa: F811
    system_id, planets, star_id = saved
    response = client.get(f"/api/nav?from=planet:{pid('planet', planets[0])}&to=planet:{pid('planet', planets[1])}&adjusted=1")
    assert response.status_code == 200, response.get_json()
    shown = response.get_json()
    assert shown["scope"] == "system"
    adjusted = shown["adjusted"]
    assert (adjusted["frame"], adjusted["unit"]) == ("system", "km")
    assert adjusted["avoided"] == [f"star:{star_id}"] and adjusted["extra"] > 0.0
    assert adjusted["names"][f"star:{star_id}"]
    assert adjusted["extra_ly"] == pytest.approx(adjusted["length_ly"] - adjusted["straight_length_ly"])
    assert adjusted["warp_times"] and adjusted["fold_times"]
    assert adjusted["moving_end"] is True and adjusted["speed_kms"] > 0.0 and adjusted["arrival_years"] > 0.0
    [leg] = shown["legs"]
    assert leg["kind"] == "within" and leg["adjusted"]["avoided"] == adjusted["avoided"]


def test_nothing_is_planned_inside_a_system_unless_asked_for(saved, client):  # noqa: F811
    _system_id, planets, _star_id = saved
    shown = client.get(f"/api/nav?from=planet:{pid('planet', planets[0])}&to=planet:{pid('planet', planets[1])}").get_json()
    assert "adjusted" not in shown and all("adjusted" not in leg for leg in shown["legs"])


def test_the_legs_out_of_and_into_a_system_are_planned_too(saved, mysql_config):
    system_id, planets, star_id = saved
    conn = store.get_connection(mysql_config)
    try:
        where = positions_at(nav_system.build_scene(conn, system_id), 0.0)
        first, second = where[f"planet:{planets[0]}"], where[f"planet:{planets[1]}"]
        edge = [-3.0e10 if first[0] > 0 else 3.0e10, 0.0, 0.0]        # the heliopause side beyond the star from the first planet
        result = {"scope": "galaxy", "legs": [
            {"kind": "out", "from": f"planet:{planets[0]}", "to": f"system:{system_id}", "system_id": system_id,
             "from_km": first, "to_km": edge},
            {"kind": "between", "from": f"system:{system_id}", "to": "system:99"},
            {"kind": "into", "from": f"system:{system_id}", "to": f"planet:{planets[1]}", "system_id": system_id,
             "from_km": [-edge[0], 0.0, 0.0], "to_km": second}]}
        nav_system.adjust_legs(conn, result)
    finally:
        conn.close()
    out, between, into = result["legs"]
    assert "adjusted" not in between and "adjusted" not in result          # the leg between systems is NAV.6's
    # The star's sphere is 1e-3 of the way out to the heliopause: below the planner's resolution, stepped round and listed apart.
    assert out["adjusted"]["ignored_below_resolution"] == [f"star:{star_id}"] and out["adjusted"]["moving_end"] is False
    assert any(point["kind"] == "tangent_in" for point in out["adjusted"]["waypoints"])
    assert out["adjusted"]["waypoints"][-1]["p"] == pytest.approx(edge)
    assert into["adjusted"]["moving_end"] is True and into["adjusted"]["ignored_below_resolution"] == [f"star:{star_id}"]


def test_the_planning_speed_is_the_slowest_warp(saved):
    assert tuning.COURSE_SYSTEM_WARP_FACTOR == 1
