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
    return [planet_mass.planet_mass(distance_au, star_kg, SNOW_LINE_AU, fe_h, 10.0 ** fe_h, stream) for _ in range(count)]


def median_earth(draws):
    return statistics.median(d.core_kg for d in draws) / EARTH_KG


def test_same_stream_same_planet():
    first = planet_mass.planet_mass(5.0, SUN_KG, SNOW_LINE_AU, 0.0, 1.0, draw.Stream(7))
    second = planet_mass.planet_mass(5.0, SUN_KG, SNOW_LINE_AU, 0.0, 1.0, draw.Stream(7))
    assert first == second


def _no_scatter(monkeypatch):
    """The disc and the spacing at their medians: nothing is drawn but the median."""
    monkeypatch.setattr(tuning, "DISC_SOLIDS_SCATTER_DEX", 0.0)
    monkeypatch.setattr(tuning, "OLIGARCH_SEPARATION_SIGMA_LN", 0.0)


def test_the_core_is_the_isolation_mass_of_the_median_spacing_at_zero_scatter(monkeypatch):
    _no_scatter(monkeypatch)
    density = formation.mmsn_surface_density_gcm2(1.0, SNOW_LINE_AU, tuning.PLANET_DISC_SOLIDS_ENHANCEMENT)
    expected = formation.isolation_mass_kg(1.0, density, SUN_KG, tuning.OLIGARCH_SEPARATION_MEDIAN)
    got = planet_mass.assembled_core_kg(1.0, SUN_KG, SNOW_LINE_AU, 1.0, draw.Stream(3))
    assert got == pytest.approx(expected)


def test_solids_grow_outward_and_jump_at_the_snow_line():
    def assembled(distance):
        stream = draw.Stream(1)
        return statistics.median(planet_mass.assembled_core_kg(distance, SUN_KG, SNOW_LINE_AU, 1.0, stream)
                                 for _ in range(400)) / EARTH_KG
    inner, just_inside, just_outside, far = (assembled(d) for d in (0.3, 2.6, 2.8, 5.0))
    assert inner < just_inside < just_outside < far
    assert just_outside > 3.0 * just_inside   # the ice boost of the snow line


def test_the_slow_outer_disc_builds_less_than_the_isolation_mass_alone(monkeypatch):
    _no_scatter(monkeypatch)
    far = 4 * tuning.PLANET_OUTER_GROWTH_AU
    tapered = planet_mass.assembled_core_kg(far, SUN_KG, SNOW_LINE_AU, 0.0, draw.Stream(3))
    monkeypatch.setattr(tuning, "PLANET_OUTER_GROWTH_AU", 1.0e6)
    untapered = planet_mass.assembled_core_kg(far, SUN_KG, SNOW_LINE_AU, 0.0, draw.Stream(3))
    assert tapered == pytest.approx(untapered * 0.25 ** tuning.PLANET_OUTER_GROWTH_EXPONENT)


def test_a_solid_core_is_pressed_softly_under_its_cap():
    cap = 10.0 * EARTH_KG
    assert planet_mass.saturate_kg(0.1 * cap, cap) == pytest.approx(0.1 * cap, rel=1e-3)
    pressed = [planet_mass.saturate_kg(f * cap, cap) for f in (0.5, 1.0, 2.0, 10.0, 1000.0)]
    assert pressed == sorted(pressed) and len(set(pressed)) == len(pressed)
    assert all(value < cap for value in pressed)


def test_a_core_inside_the_snow_line_stays_under_the_inner_cap():
    draws = masses(1.0, count=800, fe_h=0.5)
    assert max(d.core_kg for d in draws) < tuning.PLANET_INNER_CORE_MAX_EARTH * EARTH_KG


def test_cooler_stars_keep_more_of_their_solids_than_the_dust_law_gives(monkeypatch):
    _no_scatter(monkeypatch)
    light = 0.3 * SUN_KG
    got = planet_mass.assembled_core_kg(1.0, light, SNOW_LINE_AU, 1.0, draw.Stream(3))
    dust_scale = formation.disk_surface_density_scale(light, 1.0) * tuning.PLANET_DISC_SOLIDS_ENHANCEMENT
    dust = formation.isolation_mass_kg(1.0, formation.mmsn_surface_density_gcm2(1.0, SNOW_LINE_AU, dust_scale), light)
    assert got > formation.isolation_mass_kg(1.0, formation.mmsn_surface_density_gcm2(1.0, SNOW_LINE_AU, dust_scale), light,
                                             tuning.OLIGARCH_SEPARATION_MEDIAN)


def test_metal_rich_discs_build_heavier_planets():
    poor = median_earth(masses(1.0, fe_h=-0.5))
    rich = median_earth(masses(1.0, fe_h=0.3))
    assert rich > 3.0 * poor


def test_no_runaway_inside_the_snow_line():
    assert not any(d.runaway for d in masses(2.0, count=600))


def test_runaway_giants_are_heavy_and_beyond_the_line_only():
    draws = masses(4.0, count=800)
    giants = [d for d in draws if d.runaway]
    assert giants
    assert all(d.mass_kg >= tuning.GIANT_RUNAWAY_MIN_MASS_EARTH * EARTH_KG for d in giants)
    assert all(d.core_kg >= tuning.GIANT_CORE_CRITICAL_MASS_EARTH * EARTH_KG for d in giants)
    assert all(d.mass_kg == d.core_kg for d in draws if not d.runaway)


def test_a_merged_core_past_the_deuterium_limit_makes_a_giant_no_heavier_than_the_limit():
    """Two cores that merged (`planet_architecture.sculpt`) can sum to more than 13 Jupiter masses; the giant they make is held to the limit (a heavier
    body is a brown dwarf and the planet code refuses it)."""
    limit_kg = planet_mass._BROWN_DWARF_LIMIT_KG
    stream = draw.Stream(3)
    outcomes = [planet_mass.gas_outcome(1.4 * limit_kg, 5.0, SUN_KG, SNOW_LINE_AU, 0.5, stream) for _ in range(200)]
    giants = [outcome for outcome in outcomes if outcome.runaway]
    assert len(giants) > 10
    assert all(outcome.mass_kg <= limit_kg and outcome.core_kg <= limit_kg for outcome in giants)


def test_metal_poor_stars_make_fewer_giants():
    rich = sum(d.runaway for d in masses(4.0, count=800, fe_h=0.2))
    poor = sum(d.runaway for d in masses(4.0, count=800, fe_h=-1.0))
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


def test_a_core_that_finds_no_gas_stops_at_an_ice_giants_mass():
    draws = masses(4.0, count=800)
    capped = [d for d in draws if not d.runaway and d.core_kg > 10 * EARTH_KG]
    assert capped
    assert max(d.core_kg for d in capped) < tuning.GIANT_NO_GAS_CORE_MAX_EARTH * EARTH_KG


def test_giants_become_rarer_beyond_the_cut_distance():
    near = sum(d.runaway for d in masses(4.0, count=1500, seed=5))
    far = sum(d.runaway for d in masses(40.0, count=1500, seed=5))
    assert far < near / 5


def test_a_sun_like_star_has_about_a_tenth_of_a_giant_between_a_third_and_ten_au():
    """The tuning target (tier L3): [R1479] give 0.10 to 0.15 giants per Sun-like star in 0.3 to 10 AU."""
    stream = draw.Stream(21)
    slots = [0.3 * 1.5 ** index for index in range(12) if 0.3 * 1.5 ** index <= 10.0]
    giants = 0
    for _ in range(400):
        for a in slots:
            found = planet_mass.planet_mass(a, SUN_KG, SNOW_LINE_AU, 0.0, 1.0, stream)
            giants += found.runaway and found.mass_kg >= 100 * EARTH_KG
    assert 0.03 < giants / 400 < 0.4


def test_the_inner_disc_is_steeper_than_the_nebula_inside_one_au(monkeypatch):
    _no_scatter(monkeypatch)
    got = planet_mass.assembled_core_kg(0.1, SUN_KG, SNOW_LINE_AU, 1.0, draw.Stream(3))
    scale = tuning.PLANET_DISC_SOLIDS_ENHANCEMENT
    plain = formation.isolation_mass_kg(0.1, formation.mmsn_surface_density_gcm2(0.1, SNOW_LINE_AU, scale), SUN_KG,
                                        tuning.OLIGARCH_SEPARATION_MEDIAN)
    assert got == pytest.approx(plain * (0.1 ** -tuning.PLANET_INNER_DISC_STEEPENING) ** 1.5)
