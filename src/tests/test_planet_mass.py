"""The planet mass law without a class (GEN.650, physics/planet_mass.py)."""

import statistics

import pytest

from planetgen import tuning
from planetgen.physics import constants, formation, planet_mass
from planetgen.util import draw

SUN_KG = constants.SOLAR_MASS_TO_KG
EARTH_KG = constants.EARTH_MASS_TO_KG
SNOW_LINE_AU = 2.7


def masses(distance_au, count=400, star_kg=SUN_KG, fe_h=0.0, seed=1):
    stream = draw.Stream(seed)
    return [planet_mass.planet_mass(distance_au, star_kg, SNOW_LINE_AU, fe_h, stream) for _ in range(count)]


def median_earth(draws):
    return statistics.median(d.core_kg for d in draws) / EARTH_KG


def test_same_stream_same_planet():
    first = planet_mass.planet_mass(5.0, SUN_KG, SNOW_LINE_AU, 0.0, draw.Stream(7))
    second = planet_mass.planet_mass(5.0, SUN_KG, SNOW_LINE_AU, 0.0, draw.Stream(7))
    assert first == second


def test_the_core_is_the_isolation_mass_times_the_assembly_factor_at_zero_scatter(monkeypatch):
    monkeypatch.setattr(tuning, "PLANET_ASSEMBLY_SCATTER_DEX", 0.0)
    density = formation.mmsn_surface_density_gcm2(1.0, SNOW_LINE_AU, tuning.PLANET_DISC_SOLIDS_ENHANCEMENT)
    expected = formation.isolation_mass_kg(1.0, density, SUN_KG) / tuning.GIANT_IMPACT_SURVIVAL_FRACTION
    got = planet_mass.assembled_core_kg(1.0, SUN_KG, SNOW_LINE_AU, 0.0, draw.Stream(3))
    assert got == pytest.approx(expected)


def test_solids_grow_outward_and_jump_at_the_snow_line():
    inner, just_inside, just_outside, far = (median_earth(masses(d)) for d in (0.3, 2.6, 2.8, 20.0))
    assert inner < just_inside < just_outside < far
    assert just_outside > 3.0 * just_inside   # the ice boost of the snow line


def test_metal_rich_discs_build_heavier_planets():
    poor = median_earth(masses(1.0, fe_h=-0.5))
    rich = median_earth(masses(1.0, fe_h=0.3))
    assert rich > 3.0 * poor


def test_no_runaway_inside_the_snow_line():
    assert not any(d.runaway for d in masses(2.0, count=600))


def test_runaway_giants_are_heavy_and_beyond_the_line_only():
    draws = masses(15.0, count=800)
    giants = [d for d in draws if d.runaway]
    assert giants
    assert all(d.mass_kg >= tuning.GIANT_RUNAWAY_MIN_MASS_EARTH * EARTH_KG for d in giants)
    assert all(d.core_kg >= tuning.GIANT_CORE_CRITICAL_MASS_EARTH * EARTH_KG for d in giants)
    assert all(d.mass_kg == d.core_kg for d in draws if not d.runaway)


def test_metal_poor_stars_make_fewer_giants():
    rich = sum(d.runaway for d in masses(15.0, count=800, fe_h=0.2))
    poor = sum(d.runaway for d in masses(15.0, count=800, fe_h=-1.0))
    assert poor < rich


def test_a_planet_is_never_heavier_than_the_disc_allows():
    light_star = 0.1 * SUN_KG
    for d in masses(20.0, count=300, star_kg=light_star):
        assert d.mass_kg <= tuning.PLANET_MAX_STAR_MASS_RATIO * light_star * (1 + 1e-9) or \
            d.mass_kg <= tuning.GIANT_RUNAWAY_MIN_MASS_EARTH * EARTH_KG * (1 + 1e-9)


def test_the_runaway_mass_function_favours_lighter_giants():
    stream = draw.Stream(11)
    values = [planet_mass.runaway_mass_kg(SUN_KG, stream) for _ in range(2000)]
    low = tuning.GIANT_RUNAWAY_MIN_MASS_EARTH * EARTH_KG
    assert min(values) >= low
    assert max(values) <= 13.0 * constants.JUPITER_MASS_TO_KG * (1 + 1e-9)
    assert sum(v < 3 * low for v in values) > sum(v > 30 * low for v in values)
