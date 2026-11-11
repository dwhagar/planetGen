# tests/test_close_companion_gas.py

"""GEN.1475: a star with a companion inside 47 AU keeps its disk's gas in only 0.34 of the cases (Kraus et al. 2012, 2016): the chance a core finds gas is multiplied by that, and every other star draws as before."""

from types import SimpleNamespace

import pytest

from planetgen import tuning
from planetgen.physics import constants, disk_truncation, planet_architecture, planet_mass, planet_kind
from planetgen.util import draw

SUN_KG = constants.SOLAR_MASS_TO_KG
EARTH_KG = constants.EARTH_MASS_TO_KG
SNOW_LINE_AU = 2.7


def star_with(*axes_au, separation=None):
    star = SimpleNamespace(_companions=[SimpleNamespace(axis_au=axis) for axis in axes_au])
    if separation is not None:
        star.binary_separation_au = separation
    return star


def test_a_companion_inside_47_au_cuts_the_gas_share_to_0_34_and_no_other_star_is_touched():
    assert disk_truncation.gas_survival_factor(star_with(20.0)) == 0.34
    assert disk_truncation.gas_survival_factor(star_with(300.0, 46.9)) == 0.34       # the tightest companion decides
    assert disk_truncation.gas_survival_factor(star_with(47.0)) == 1.0
    assert disk_truncation.gas_survival_factor(star_with(120.0)) == 1.0
    assert disk_truncation.gas_survival_factor(star_with()) == 1.0                      # a lone star
    assert disk_truncation.gas_survival_factor(SimpleNamespace()) == 1.0                # a star the system attached nothing to
    assert disk_truncation.gas_survival_factor(star_with(separation=0.12)) == 0.34      # the proxy of a close pair


def runaways(factor, count=3000, seed=11):
    stream = draw.Stream(seed)
    critical = tuning.GIANT_CORE_CRITICAL_MASS_EARTH * EARTH_KG * 1.5
    return sum(planet_mass.gas_outcome(critical, 4.0, SUN_KG, SNOW_LINE_AU, 0.0, stream, factor).runaway for _ in range(count))


def test_the_factor_scales_the_chance_a_core_finds_gas():
    full, cut = runaways(1.0), runaways(0.34)
    assert full > 300 and cut == pytest.approx(0.34 * full, rel=0.2)
    assert runaways(0.0) == 0


def test_a_factor_of_one_draws_exactly_what_was_drawn_before():
    stream_a, stream_b = draw.Stream(5), draw.Stream(5)
    for distance in (0.3, 1.0, 4.0, 9.0, 20.0):
        assert planet_mass.planet_mass(distance, SUN_KG, SNOW_LINE_AU, 0.0, 1.0, stream_a) == planet_mass.planet_mass(distance, SUN_KG, SNOW_LINE_AU, 0.0, 1.0, stream_b, 1.0)
    assert planet_kind.draw_formation(5.0, SUN_KG, SNOW_LINE_AU, 0.0, 1.0, draw.Stream(3)) == planet_kind.draw_formation(5.0, SUN_KG, SNOW_LINE_AU, 0.0, 1.0, draw.Stream(3), 1.0)


def test_an_architecture_built_with_the_factor_makes_fewer_giants_and_the_rocky_planets_stay():
    def built(factor):
        giants = rocky = 0
        for seed in range(200):
            architecture = planet_architecture.build(SUN_KG, SNOW_LINE_AU, 1.0, 0.0, 60.0, stream=draw.Stream(seed), gas_factor=factor)
            giants += sum(1 for planet in architecture.planets if planet.runaway)
            rocky += sum(1 for planet in architecture.planets if not planet.runaway)
        return giants, rocky
    full_giants, full_rocky = built(1.0)
    cut_giants, cut_rocky = built(0.34)
    assert full_giants > 20 and cut_giants < 0.6 * full_giants
    assert cut_rocky >= full_rocky              # a core that finds no gas stays a solid or ice-giant core: no planet is lost
