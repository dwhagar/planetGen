# tests/test_body_positions.py

"""
MAP.70: `physics/body_positions.positions_at` -- where every body of a
system is at any time, from the scene `GET /api/systems/<id>/scene`
returns. `sample_scene` is a hand-built scene (a close pair with a planet,
moon and both kinds of comet) that `test_js_unit.py` also hands to
`static/orbitpositions.js` to check the two copies against each other.
"""

import math

import pytest

from planetgen.physics.body_positions import AU_KM, positions_at


def circular(around, distance_au, period, incl, node, phase, **extra):
    return dict({"around": around, "distance_km": distance_au * AU_KM, "period_years": period,
                 "inclination_deg": incl, "ascending_node_deg": node, "phase_deg": phase}, **extra)


def sample_scene():
    return {
        "stars": [
            {"ref": "star:1", "orbit": None},
            {"ref": "star:2", "orbit": circular("barycenter", 0.3, 0.2, 20.0, 40.0, 100.0, secondary_mass_fraction=0.25)},
        ],
        "planets": [
            {"ref": "planet:1", "orbit": circular("barycenter", 1.5, 1.8, 3.0, 80.0, 10.0),
             "moons": [{"ref": "moon:1", "orbit": circular("planet:1", 0.002, 0.03, 5.0, 10.0, 200.0)}]},
            {"ref": "planet:2", "orbit": circular("barycenter", 5.0, 11.0, 0.0, 0.0, 300.0), "moons": []},
        ],
        "comets": [
            {"ref": "comet:1", "orbit": {"around": "barycenter", "type": "elliptical", "kepler": {
                "perihelion_distance_km": 0.6 * AU_KM, "eccentricity": 0.97, "inclination_deg": 162.0,
                "arg_periapsis_deg": 111.0, "ascending_node_deg": 58.0, "period_years": 76.0,
                "mean_anomaly_deg": 35.0, "parabolic_mean_anomaly": None, "primary_mass_solar": 1.2}}},
            {"ref": "comet:2", "orbit": {"around": "barycenter", "type": "parabolic", "kepler": {
                "perihelion_distance_km": 2.0 * AU_KM, "eccentricity": 1.0, "inclination_deg": 70.0,
                "arg_periapsis_deg": 200.0, "ascending_node_deg": 15.0, "period_years": None,
                "mean_anomaly_deg": None, "parabolic_mean_anomaly": -3.0, "primary_mass_solar": 1.2}}},
        ],
    }


def norm(v):
    return math.sqrt(sum(c * c for c in v))


def sub(a, b):
    return [x - y for x, y in zip(a, b)]


def test_a_circular_orbit_keeps_its_radius_and_returns_after_a_period():
    scene = sample_scene()
    start = positions_at(scene, 0.0)
    later = positions_at(scene, 1.8)
    assert norm(start["planet:1"]) == pytest.approx(1.5 * AU_KM, rel=1e-9)
    assert norm(positions_at(scene, 0.7)["planet:1"]) == pytest.approx(1.5 * AU_KM, rel=1e-9)
    assert later["planet:1"] == pytest.approx(start["planet:1"], rel=1e-6, abs=1.0)
    assert positions_at(scene, 1.8 / 2)["planet:1"] != pytest.approx(start["planet:1"], rel=1e-3)


def test_a_moon_is_relative_to_its_planet():
    scene = sample_scene()
    here = positions_at(scene, 0.4)
    assert norm(sub(here["moon:1"], here["planet:1"])) == pytest.approx(0.002 * AU_KM, rel=1e-9)


def test_a_close_pair_balances_on_the_barycenter():
    scene = sample_scene()
    for years in (0.0, 0.05, 3.0):
        p = positions_at(scene, years)
        gap = sub(p["star:2"], p["star:1"])
        assert norm(gap) == pytest.approx(0.3 * AU_KM, rel=1e-9)
        # mass fraction 0.25 for the second star: weighted centre sits at the origin
        balance = [0.75 * a + 0.25 * b for a, b in zip(p["star:1"], p["star:2"])]
        assert norm(balance) < 1e-3 * AU_KM * 1e-6


def test_a_wide_pair_orbits_the_first_star():
    scene = sample_scene()
    scene["stars"][1]["orbit"] = circular("star:1", 30.0, 100.0, 0.0, 0.0, 0.0)
    p = positions_at(scene, 0.0)
    assert p["star:1"] == [0.0, 0.0, 0.0]
    assert norm(p["star:2"]) == pytest.approx(30.0 * AU_KM, rel=1e-9)


def test_an_elliptical_comet_is_fast_near_perihelion_and_periodic():
    scene = sample_scene()
    near = norm(positions_at(scene, 0.0)["comet:1"])
    assert near >= 0.6 * AU_KM * (1 - 1e-9)
    assert positions_at(scene, 76.0)["comet:1"] == pytest.approx(positions_at(scene, 0.0)["comet:1"], rel=1e-6, abs=1e3)
    far = norm(positions_at(scene, 38.0 - 35.0 / 360.0 * 76.0)["comet:1"])  # aphelion: mean anomaly 180 degrees
    assert far == pytest.approx(0.6 * AU_KM * (1 + 0.97) / (1 - 0.97), rel=1e-6)


def test_a_parabolic_comet_passes_perihelion_at_the_right_time():
    scene = sample_scene()
    mu = 4 * math.pi ** 2 * 1.2
    rate = math.sqrt(mu / (2 * 2.0 ** 3))
    years_to_perihelion = 3.0 / rate
    assert norm(positions_at(scene, years_to_perihelion)["comet:2"]) == pytest.approx(2.0 * AU_KM, rel=1e-9)
    assert norm(positions_at(scene, 0.0)["comet:2"]) > 2.0 * AU_KM


def test_the_scene_positions_at_the_epoch_match_the_stored_ones(mysql_config):
    from planetgen.db import store
    from planetgen.galaxy.sector import SpaceSector
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.web.maps.systemscene import build_scene

    sector = SpaceSector("Position Sector", edge_ly=11.5)
    for _ in range(400):
        cfg = SystemConfig()
        cfg.PLANETS = True
        cfg.MAX_PLANETS = True
        system = StarSystem(system_config=cfg)
        if any(getattr(p, "moons", None) for p in system.planets) and system.comets and len(system.stars) > 1:
            break
    sector.add_system(system, position=(1.0, 1.0, 1.0), system_config=cfg)
    store.save_sector(sector, config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        scene = build_scene(conn, conn.execute("SELECT id FROM star_systems LIMIT 1").fetchone()["id"])
    finally:
        conn.close()
    here = positions_at(scene, 0.0)
    for star in scene["stars"]:
        if star["orbit"] is not None and star["orbit"]["around"] == "barycenter":
            assert here[star["ref"]] == pytest.approx(star["position_km"], rel=1e-6, abs=1e4)
    for planet in scene["planets"]:
        origin = here[planet["orbit"]["around"]] if planet["orbit"]["around"] != "barycenter" else [0.0, 0.0, 0.0]
        assert sub(here[planet["ref"]], origin) == pytest.approx(planet["position_km"], rel=1e-6, abs=1e4)
        for moon in planet["moons"]:
            assert sub(here[moon["ref"]], here[planet["ref"]]) == pytest.approx(moon["position_km"], rel=1e-6, abs=1e4)
    for comet in scene["comets"]:
        assert here[comet["ref"]] == pytest.approx(comet["position_km"], rel=1e-5, abs=1e6)
