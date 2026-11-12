# tests/test_bare_body.py

"""
GEN.1548 (built ahead by Lane 07 for Lane 05): the stand-ins of the backward edges of the forward build: what a body has before its air is solved.
"""

import math
import types

import pytest

from planetgen.generation.config import SystemConfig
from planetgen.generation.star import Star
from planetgen.physics import bare_body, constants, interior_heat
from planetgen.util import draw
from tests import class_support

SOLAR_L_W = 3.828e26
EARTH_MASS_KG = 5.972e24


@pytest.fixture(scope="module")
def star():
    config = SystemConfig()
    config.STAR_TYPE = "G2V"
    with draw.bound("bare-body-star"):
        return Star(config)


def planet_at(star, au, mass_earth, seed=2):
    config = SystemConfig()
    config.STAR_TYPE = "G2V"
    distance = au * math.sqrt(star.luminosity / SOLAR_L_W)
    with draw.bound(f"bare-body-planet:{au}:{mass_earth}"):
        return class_support.draw_planet(star, distance, seed, config, mass=mass_earth * EARTH_MASS_KG, moon_count=0)


@pytest.fixture(scope="module")
def earth_like(star):
    return planet_at(star, 1.0, 1.0)


def test_the_bare_values_draw_nothing_and_repeat(earth_like):
    with draw.bound(draw.Stream("bare-body-run")):
        before = draw.getstate()
        first = (bare_body.scale_height_km(earth_like), bare_body.equilibrium_k(earth_like), bare_body.reference_albedo(earth_like))
        after = draw.getstate()
    assert before == after
    assert first == (bare_body.scale_height_km(earth_like), bare_body.equilibrium_k(earth_like), bare_body.reference_albedo(earth_like))


def test_an_earth_twins_bare_scale_height_and_temperature_are_the_airs_to_within_a_fifth(earth_like):
    """The bare height is for the margin that keeps a moon clear of the air: within 25 percent of the solved air's at 1 au."""
    assert earth_like.scale_height
    assert bare_body.scale_height_km(earth_like) == pytest.approx(earth_like.scale_height, rel=0.25)
    assert bare_body.equilibrium_k(earth_like) == pytest.approx(255.0, abs=15.0)


def test_a_bare_rocky_body_uses_the_molar_mass_of_the_stand_in(earth_like):
    expected = constants.R * bare_body.equilibrium_k(earth_like) / (bare_body.BARE_ROCKY_MOLAR_MASS_KG_MOL * bare_body.gravity_m_s2(earth_like)) / 1000.0
    assert bare_body.scale_height_km(earth_like) == pytest.approx(expected)


def test_a_giant_uses_the_molar_mass_of_a_hydrogen_helium_envelope(star):
    giant = planet_at(star, 5.2, 318.0)
    assert giant.body_type == "g"
    expected = constants.R * bare_body.equilibrium_k(giant) / (2.35e-3 * bare_body.gravity_m_s2(giant)) / 1000.0
    assert bare_body.scale_height_km(giant) == pytest.approx(expected)


def test_the_surface_temperature_is_the_bodys_own_once_it_has_one_and_else_its_equilibrium(star):
    far = planet_at(star, 5.0, 1.0)
    bare = types.SimpleNamespace(**{**vars(far), "surface_temperature": None})
    assert bare_body.surface_k(bare) == bare_body.equilibrium_k(bare)
    assert abs(bare_body.surface_k(bare) - interior_heat.SURFACE_TEMPERATURE_K) > 50.0, "at 5 au the old 273 K fallback is far off"
    own = types.SimpleNamespace(**{**vars(far), "surface_temperature": 300.0})
    assert bare_body.surface_k(own) == 300.0


def test_a_body_with_no_size_has_no_scale_height(earth_like):
    nothing = types.SimpleNamespace(**{**vars(earth_like), "radius": None})
    assert bare_body.scale_height_km(nothing) == 0.0


def test_the_thermal_tide_is_known_only_once_the_planet_has_a_pressure():
    assert not bare_body.thermal_tide_is_known(types.SimpleNamespace(atmospheric_pressure=None))
    assert not bare_body.thermal_tide_is_known(types.SimpleNamespace(atmospheric_pressure=0.0))
    assert bare_body.thermal_tide_is_known(types.SimpleNamespace(atmospheric_pressure=1.0e5))
