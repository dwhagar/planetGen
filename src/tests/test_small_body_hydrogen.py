# tests/test_small_body_hydrogen.py

"""
GEN.1199: a small body far from its star does not keep a hydrogen air. The escape filter of the drawn mixture (`mixture_escape`) reads the star's mass to find the
XUV; the subject built from a body once left it out, so the filter saw no XUV, its energy budget was nothing and no gas left any body. A 0.02 Earth-mass moon of a giant at
5 AU then had a hydrogen air and the class T of a gas dwarf.
"""

import math

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.star import Star
from planetgen.physics import activity, constants, thermosphere
from planetgen.util import draw
from tests import class_support

AIR = {"H2": 0.95, "N2": 0.04, "Ar": 0.01}
SUN_W = constants.SOLAR_LUMINOSITY


def kept_hydrogen(mass_earth, teq_k, distance_au, age_gy=4.5, density=2500.0, pressure_pa=1.0e4):
    """The share of a hydrogen air left at `age_gy` on a body of rock and ice at `distance_au` from the Sun."""
    mass = mass_earth * constants.EARTH_MASS_TO_KG
    radius = (3.0 * mass / (4.0 * math.pi * density)) ** (1.0 / 3.0)
    flux = activity.coronal_xuv_at(1.0, age_gy, SUN_W) / (4.0 * math.pi * (distance_au * constants.AU_TO_M) ** 2) / constants.EARTH_XUV_FLUX_W_M2
    kept, _exobase = thermosphere.kept_shares(AIR, mass, radius, teq_k, pressure_pa, age_gy, flux, 1.0, SUN_W)
    return kept["H2"]


@pytest.mark.parametrize("distance_au,teq_k,mass_earth", [(1.0, 255.0, 5.0), (5.2, 110.0, 5.0), (5.2, 110.0, 0.0226)])
def test_a_body_this_small_this_near_its_star_loses_its_hydrogen(distance_au, teq_k, mass_earth):
    assert kept_hydrogen(mass_earth, teq_k, distance_au) < 1e-3


@pytest.mark.parametrize("distance_au,teq_k,mass_earth", [(30.0, 45.0, 2.0), (50.0, 35.0, 0.5), (100.0, 25.0, 0.1)])
def test_a_body_far_from_its_star_keeps_its_hydrogen_above_a_mass_that_falls_with_distance(distance_au, teq_k, mass_earth):
    """The mass that keeps half of a hydrogen air over 4.5 Gy: 23 Earth masses at 1 AU, 14 at 5 AU, 8 at 10, 2.2 at 20, 0.8 at 30, 0.2 at 50, 0.04 at 100 (icy rock, 10 kPa)."""
    assert kept_hydrogen(mass_earth, teq_k, distance_au) > 0.5


def test_no_small_moon_of_a_giant_carries_a_hydrogen_air_or_class_t():
    config = SystemConfig()
    config.STAR_TYPE = "G2V"
    with draw.bound("small-body-hydrogen:star"):
        star = Star(config)
    giant = 300.0 * constants.EARTH_MASS_TO_KG
    seen = 0
    for seed in range(12):
        planet = class_support.draw_planet(star, 5.2, seed, config, mass=giant, zone_override="c", moon_count=6)
        for moon in planet.moons:
            if moon.mass / constants.EARTH_MASS_TO_KG >= 0.2 or moon.body_type != "t":
                continue
            seen += 1
            partials = getattr(moon, "gas_partials_pa", None) or {}
            total = sum(partials.values())
            assert not total or partials.get("H2", 0.0) / total < 0.5, (seed, moon.planet_class, partials)
            assert moon.planet_class != "T", (seed, moon.mass / constants.EARTH_MASS_TO_KG)
    assert seen > 10
