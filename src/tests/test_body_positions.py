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

from planetgen.physics import constants
from planetgen.physics.body_positions import positions_at
from planetgen.physics.constants import AU_TO_KM as AU_KM
from planetgen.physics.orbits import circular_orbital_velocity_au_per_year, orbital_position_au

KMS_PER_AU_YR = constants.AU_TO_KM / constants.SECONDS_PER_YEAR


def circular(around, distance_au, period, incl, node, phase, **extra):
    """A circular orbit by its phase: how a pair's second star is described."""
    return dict({"around": around, "distance_km": distance_au * AU_KM, "period_years": period,
                 "inclination_deg": incl, "ascending_node_deg": node, "phase_deg": phase}, **extra)


def state_orbit(around, distance_au, period, incl, node, phase, epoch_unix=None):
    """A planet's or moon's orbit as the scene gives it (GEN.289): the elements it keeps and the state vector it is,
    here the one a circular orbit has at `phase`."""
    return {"around": around, "distance_km": distance_au * AU_KM, "period_years": period,
            "inclination_deg": incl, "ascending_node_deg": node,
            "position_km": [c * AU_KM for c in orbital_position_au(distance_au, incl, node, phase)],
            "velocity_kms": [c * KMS_PER_AU_YR for c in
                             circular_orbital_velocity_au_per_year(distance_au, incl, node, phase, period)],
            "epoch_unix": epoch_unix}


SCENE_EPOCH = 1_760_000_000.0
"""float: The sample scene's epoch (unix seconds): the time its bodies' stored states are taken to hold at."""

LATE_BY_S = 40 * 86400.0
"""float: How long before the scene epoch the sample scene's second planet and moon were last stepped."""


def sample_scene():
    return {
        "epoch_unix": SCENE_EPOCH,
        "tree": {
            "nodes": [{"node_no": 0, "parent_node_no": None, "slot": 0, "mass_solar": 1.0, "separation_km": 0.3 * AU_KM,
                       "eccentricity": 0.0, "inclination_deg": 20.0, "ascending_node_deg": 40.0,
                       "arg_periapsis_deg": 0.0, "mean_anomaly_deg": 100.0, "epoch_unix": SCENE_EPOCH,
                       "period_years": 0.2, "kind": "close", "interaction": "detached", "flags": 0}],
            "members": [{"ref": "star:1", "member_no": 0, "parent_node_no": 0, "slot": 0, "mass_solar": 0.75},
                        {"ref": "star:2", "member_no": 1, "parent_node_no": 0, "slot": 1, "mass_solar": 0.25}],
            "origin": None},
        "stars": [
            {"ref": "star:1", "orbit": None},
            {"ref": "star:2", "orbit": circular("barycenter", 0.3, 0.2, 20.0, 40.0, 100.0, secondary_mass_fraction=0.25)},
        ],
        "planets": [
            {"ref": "planet:1", "orbit": state_orbit("barycenter", 1.5, 1.8, 3.0, 80.0, 10.0),
             "moons": [{"ref": "moon:1", "orbit": state_orbit("planet:1", 0.002, 0.03, 5.0, 10.0, 200.0, SCENE_EPOCH - 3 * 86400.0)}]},
            {"ref": "planet:2", "orbit": state_orbit("barycenter", 5.0, 11.0, 0.0, 0.0, 300.0, SCENE_EPOCH - LATE_BY_S), "moons": []},
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


def _wide_scene(separation_au, eccentricity):
    scene = sample_scene()
    scene["stars"][1]["orbit"] = circular("star:1", separation_au, 100.0, 0.0, 0.0, 0.0)
    node = scene["tree"]["nodes"][0]
    node.update(separation_km=separation_au * AU_KM, eccentricity=eccentricity, inclination_deg=0.0,
                ascending_node_deg=0.0, mean_anomaly_deg=0.0, period_years=100.0, kind="wide")
    scene["tree"]["origin"] = "star:1"
    return scene


def eccentric_wide_scene():
    """MAP.255: a wide pair on an eccentric, tilted orbit with the second star's planet, which `test_js_unit.py` also hands
    to `static/orbitpositions.js`: the browser draws the pair on the stored ellipse and not on a circle."""
    scene = _wide_scene(30.0, 0.6)
    scene["tree"]["nodes"][0].update(inclination_deg=25.0, ascending_node_deg=70.0, arg_periapsis_deg=130.0, mean_anomaly_deg=40.0)
    scene["planets"] = [{"ref": "planet:3", "orbit": state_orbit("star:2", 1.2, 1.4, 5.0, 20.0, 60.0, SCENE_EPOCH), "moons": []}]
    scene["comets"] = []
    return scene


def triple_scene(epoch_unix=SCENE_EPOCH):
    """MAP.255: an eccentric close pair about the centre of a third star on a wide orbit (`scene_epoch` None: the scene has
    none, the tree's first node says where time zero is)."""
    scene = sample_scene()
    scene["epoch_unix"] = epoch_unix
    scene["stars"].append({"ref": "star:3", "orbit": None})
    inner, outer = scene["tree"]["nodes"][0], dict(scene["tree"]["nodes"][0])
    inner.update(parent_node_no=1, slot=0, eccentricity=0.35, arg_periapsis_deg=75.0, mean_anomaly_deg=200.0)
    outer.update(node_no=1, parent_node_no=None, slot=0, mass_solar=2.0, separation_km=40.0 * AU_KM, eccentricity=0.45,
                 period_years=150.0, kind="wide", mean_anomaly_deg=310.0, inclination_deg=12.0, ascending_node_deg=95.0,
                 arg_periapsis_deg=20.0, epoch_unix=SCENE_EPOCH - 7 * 86400.0)
    scene["tree"]["nodes"].append(outer)
    scene["tree"]["members"].append({"ref": "star:3", "member_no": 2, "parent_node_no": 1, "slot": 1, "mass_solar": 1.0})
    scene["planets"], scene["comets"] = [], []
    return scene


def test_a_wide_pair_orbits_the_first_star():
    scene = _wide_scene(30.0, 0.0)
    p = positions_at(scene, 0.0)
    assert p["star:1"] == [0.0, 0.0, 0.0]
    assert norm(p["star:2"]) == pytest.approx(30.0 * AU_KM, rel=1e-9)


def test_an_eccentric_wide_pair_is_closest_at_periapsis_and_farthest_half_a_period_later():
    """The tree's orbit is a Kepler orbit, not a circle: at mean anomaly 0 the pair is at periapsis a(1 - e)."""
    scene = _wide_scene(30.0, 0.6)
    assert norm(positions_at(scene, 0.0)["star:2"]) == pytest.approx(30.0 * 0.4 * AU_KM, rel=1e-9)
    assert norm(positions_at(scene, 50.0)["star:2"]) == pytest.approx(30.0 * 1.6 * AU_KM, rel=1e-6)


@pytest.mark.parametrize("phase", [0.0, 37.0, 100.0, 215.5, 359.0])
@pytest.mark.parametrize("years", [0.0, 0.07, 1.3])
def test_a_circular_node_gives_the_positions_the_stored_phase_did(phase, years):
    """GEN.311: reading a pair's node moves the stars where the old circular orbit moved them."""
    scene = sample_scene()
    scene["tree"]["nodes"][0]["mean_anomaly_deg"] = phase
    scene["stars"][1]["orbit"]["phase_deg"] = phase
    p = positions_at(scene, years)
    here = orbital_position_au(0.3, 20.0, 40.0, (phase + 360.0 * years / 0.2) % 360.0)
    fraction = 0.25
    assert p["star:2"] == pytest.approx([c * AU_KM * (1 - fraction) for c in here], rel=1e-9, abs=1.0)
    assert p["star:1"] == pytest.approx([-c * AU_KM * fraction for c in here], rel=1e-9, abs=1.0)


def test_a_triple_places_the_close_pair_about_the_outer_centre():
    scene = sample_scene()
    scene["stars"].append({"ref": "star:3", "orbit": None})
    scene["tree"]["nodes"][0].update(parent_node_no=1, slot=0, mass_solar=1.0)
    scene["tree"]["nodes"].append(dict(scene["tree"]["nodes"][0], node_no=1, parent_node_no=None, slot=0,
                                       mass_solar=2.0, separation_km=40.0 * AU_KM, period_years=150.0, kind="wide",
                                       mean_anomaly_deg=0.0, inclination_deg=0.0, ascending_node_deg=0.0))
    scene["tree"]["members"].append({"ref": "star:3", "member_no": 2, "parent_node_no": 1, "slot": 1, "mass_solar": 1.0})
    p = positions_at(scene, 0.0)
    pair_centre = [0.75 * a + 0.25 * b for a, b in zip(p["star:1"], p["star:2"])]
    # The pair's centre of mass is 1/3 of the way from the pair to the third star, which sits at 40 AU
    assert norm(sub(p["star:3"], pair_centre)) == pytest.approx(40.0 * AU_KM, rel=1e-9)
    total = [(1.0 * c + 1.0 * t) for c, t in zip(pair_centre, p["star:3"])]
    assert norm(total) < 1e-6 * AU_KM


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


def test_a_body_last_stepped_before_the_scene_epoch_is_carried_on_to_it():
    scene = sample_scene()
    stamped = scene["planets"][1]["orbit"]
    assert stamped["epoch_unix"] == SCENE_EPOCH - LATE_BY_S
    here = positions_at(scene, 0.0)["planet:2"]
    unstamped = sample_scene()
    unstamped["planets"][1]["orbit"]["epoch_unix"] = None
    carried = positions_at(unstamped, LATE_BY_S / constants.SECONDS_PER_YEAR)["planet:2"]
    assert here == pytest.approx(carried, rel=1e-9, abs=1.0)
    assert here != pytest.approx(positions_at(unstamped, 0.0)["planet:2"], rel=1e-3)
    assert norm(here) == pytest.approx(5.0 * AU_KM, rel=1e-9)


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
        cfg.MOONS = True
        cfg.COMETS = True  # a comet by construction, not by the draw (TEST.130)
        system = StarSystem(system_config=cfg)
        if any(getattr(p, "moons", None) for p in system.planets) and system.comets and len(system.stars) > 1:
            break
    else:
        pytest.fail("400 systems built and none had a planet with moons, a comet and a second star")
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
        # Like a planet's, a comet's stored position is relative to what it goes round (GEN.165: a comet
        # round the second star of a wide pair was compared against the barycenter and missed by 64 AU).
        origin = here[comet["orbit"]["around"]] if comet["orbit"]["around"] != "barycenter" else [0.0, 0.0, 0.0]
        assert sub(here[comet["ref"]], origin) == pytest.approx(comet["position_km"], rel=1e-5, abs=1e6)


@pytest.mark.parametrize("orbit_type", ["elliptical", "parabolic"])
def test_a_comets_position_at_the_epoch_is_the_one_it_was_generated_with(orbit_type):
    """TEST.130: the scene's comet path (`positions_at`) and the position a comet is generated and stored with
    are the same function of the same orbit, for 2,000 comets of each type (60,000 drawn once, none differed)."""
    from planetgen.generation.comet import Comet
    from planetgen.generation.config import SystemConfig
    from planetgen.physics.body_positions import _comet_relative_km
    from planetgen.util import draw

    draw.set_run_seed(130)
    for _ in range(2000):
        mass = draw.uniform(0.08, 5.0)
        comet = Comet(SystemConfig(), mass, orbit_type=orbit_type)
        orbit = {"type": comet.orbit_type, "kepler": {
            "perihelion_distance_km": comet.perihelion_distance_au * AU_KM, "eccentricity": comet.eccentricity,
            "inclination_deg": comet.inclination_deg, "arg_periapsis_deg": comet.arg_periapsis_deg,
            "ascending_node_deg": comet.ascending_node_deg, "period_years": comet.orbital_period_years,
            "mean_anomaly_deg": comet.mean_anomaly_deg, "parabolic_mean_anomaly": comet.parabolic_mean_anomaly,
            "primary_mass_solar": mass}}
        stored = [comet.position_x_au * AU_KM, comet.position_y_au * AU_KM, comet.position_z_au * AU_KM]
        assert _comet_relative_km(orbit, 0.0) == pytest.approx(stored, rel=1e-5, abs=1e6), vars(comet)


def test_a_planet_of_a_pair_below_a_further_star_goes_round_the_pairs_centre_not_the_barycenter():
    """A close pair that hangs below another orbit has its planets round `node:N`, the pair's centre of mass, which is
    away from the barycenter the third star pulls it from."""
    scene = triple_scene()
    scene["planets"] = [{"ref": "planet:1", "orbit": state_orbit("node:0", 1.5, 1.8, 3.0, 80.0, 10.0, SCENE_EPOCH), "moons": []}]
    here = positions_at(scene, 0.0)
    members = [m for m in scene["tree"]["members"] if m["parent_node_no"] == 0]
    total = sum(m["mass_solar"] for m in members)
    centre = [sum(m["mass_solar"] * here[m["ref"]][axis] for m in members) / total for axis in range(3)]
    away = [here["planet:1"][axis] - centre[axis] for axis in range(3)]
    assert norm(away) == pytest.approx(1.5 * AU_KM, rel=1e-9)
    assert norm(centre) > 1.0 * AU_KM, "the pair's centre is far from the barycenter"
