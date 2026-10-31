# tests/test_gen_hosting_wiring.py

"""
GEN.761: the moon draw asks the hosting gate (`physics.hosting_gate`, GEN.579). `generate_moons` stops at a host the gate
refuses, a small body (under `HOSTING_SMALL_BODY_MAX_RADIUS_KM`) is named for what it is to the gate
(`planets.hosting_population`) and has moons at the share of the gate's `satellite-share` table for its population and size,
and every other planet keeps the even chance.
"""

import pathlib
import types

import pytest

from planetgen import tuning
from planetgen.generation.config import SystemConfig
from planetgen.generation.planet import Planet
from planetgen.generation.star import Star
from planetgen.physics import constants, hosting_gate, planets as planet_physics
from planetgen.util import draw


@pytest.fixture(scope="module")
def sun_like():
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    return Star(cfg)


# What the caller gives to get each sort of host: the planet's kind and size follow from its mass.
GIVEN_MASS_KG = {"giant": 300 * constants.EARTH_MASS_TO_KG, "small": 2.0e19, "earth": constants.EARTH_MASS_TO_KG}
SMALL_RADIUS_KM = tuning.HOSTING_SMALL_BODY_MAX_RADIUS_KM / 2.0


def _planet(star, kind, zone, **kwargs):
    inner, outer = star.habitable_zone
    distance = {"h": inner * 0.5, "e": (inner + outer) / 2.0, "c": outer * 2.0}[zone]
    return Planet(SystemConfig(), star, star.habitable_zone, distance, mass=GIVEN_MASS_KG[kind], **kwargs)


def test_the_moon_draw_stops_at_a_host_the_gate_refuses(sun_like, monkeypatch):
    asked = []

    def refuse(host):
        asked.append(host)
        return hosting_gate.Verdict(False, 2, "no stable orbit")

    with draw.bound(5):
        assert _planet(sun_like, "giant", "c", moon_count=3).moons, "the real gate lets a gas giant have its moons"
    monkeypatch.setattr(hosting_gate, "can_host", refuse)
    with draw.bound(5):
        refused = _planet(sun_like, "giant", "c", moon_count=3)
    assert refused.moons == [] and len(asked) == 1
    host = asked[0]
    assert host.population == "planet" and host.fluid and host.parent_mass_kg == pytest.approx(sun_like.mass)
    assert host.diameter_m == pytest.approx(2000.0 * refused.radius) and host.mass_kg == refused.mass


def test_a_directed_moon_count_still_goes_through_the_gate(sun_like):
    counts = []
    for seed in range(24):
        with draw.bound(seed):
            planet = _planet(sun_like, "earth", "e", moon_count=2)
        assert hosting_gate.can_host(planet_physics.hosting_host(planet)).allowed
        counts.append(len(planet.moons))
    # the directive holds where there is room: never more than asked, and a planet whose orbits cannot take the second moon keeps one (about 1 in 12)
    assert max(counts) == 2 and min(counts) >= 1 and counts.count(2) >= 0.8 * len(counts)


@pytest.mark.parametrize("zone,distance_au,inclination,expected", [
    ("h", 0.3, 0.0, "inner-small-body"),
    ("e", 1.0, 0.0, "inner-small-body"),
    ("c", 2.8, 0.0, "belt-small-body"),
    ("c", 9.0 * 2.7, 30.0, "belt-small-body"),
    ("c", 11.0 * 2.7, 2.0, "trans-neptunian-cold"),
    ("c", 11.0 * 2.7, 20.0, "trans-neptunian-hot"),
])
def test_a_small_body_is_named_for_its_place_and_orbit(zone, distance_au, inclination, expected):
    body = types.SimpleNamespace(radius=SMALL_RADIUS_KM, zone=zone, distance=distance_au, orbital_inclination_deg=inclination,
                                 star=types.SimpleNamespace(snow_line_au=2.7))
    assert planet_physics.hosting_population(body) == expected


def test_every_other_planet_is_just_a_planet_to_the_gate():
    body = types.SimpleNamespace(radius=60000.0, zone="c", distance=5.0, orbital_inclination_deg=0.0,
                                 star=types.SimpleNamespace(snow_line_au=2.7))
    assert planet_physics.hosting_population(body) == "planet"


def test_a_small_body_has_the_share_of_the_table_and_a_planet_has_none(sun_like):
    with draw.bound(3):
        small, big = _planet(sun_like, "small", "h"), _planet(sun_like, "earth", "e")
    assert planet_physics.satellite_share(small) == pytest.approx(0.05)       # the small body row
    assert planet_physics.satellite_share(big) is None


@pytest.mark.parametrize("share,expected", [(0.2, 0.2), (None, 0.5)])
def test_the_chance_of_moons_is_the_tables_share_else_an_even_chance(monkeypatch, share, expected):
    monkeypatch.setattr(planet_physics, "satellite_share", lambda planet: share)
    body = types.SimpleNamespace(system_config=SystemConfig())
    runs = 6000
    with draw.bound(11):
        hits = sum(Planet._draw_moons(body) for _ in range(runs))
    sigma = (runs * expected * (1.0 - expected)) ** 0.5
    assert abs(hits - runs * expected) < 4.0 * sigma


def test_the_hosting_gate_has_a_caller_in_the_generator():
    root = pathlib.Path(__file__).resolve().parents[1] / "planetgen"
    users = [path.relative_to(root).as_posix() for path in root.rglob("*.py")
             if path.name != "hosting_gate.py" and "hosting_gate" in path.read_text(encoding="utf-8")]
    assert "physics/planets.py" in users
