"""GEN.1611: the snow line a star's disk was built round (`physics/formation_snow_line.py`), and that a star stores it."""

import math
import types

import pytest

from planetgen import tuning
from planetgen.generation.config import SystemConfig
from planetgen.generation.star import Star
from planetgen.physics import birth_luminosity, constants, formation, formation_snow_line as fsl
from planetgen.util import draw


def test_the_luminosity_is_the_birth_luminosity_at_the_disk_age():
    for mass in (0.3, 1.0, 2.5, 8.0):
        assert fsl.formation_luminosity_sol(mass, 0.0) == birth_luminosity.birth_luminosity_sol(mass, 0.0, tuning.FORMATION_DISK_AGE_GY)
    assert tuning.FORMATION_DISK_AGE_GY == 0.01


def test_a_sun_like_disk_saw_less_light_than_the_sun_gives_and_an_m_dwarfs_more_than_it_gives():
    assert fsl.formation_luminosity_sol(1.0, 0.0) < 1.0
    assert fsl.formation_snow_line_au(1.0, 0.0) < formation.snow_line_au(constants.SOLAR_LUMINOSITY)
    assert fsl.formation_luminosity_sol(0.3, 0.0) > 0.03        # 0.02 of the Sun on the main sequence


@pytest.mark.parametrize("mass, low, high", [
    (4.36, 20.0, 80.0),        # a white dwarf born at 4.36 solar masses: its disk saw a B star (was 0.02 AU)
    (1.43, 2.0, 8.0),
    (1.04, 0.8, 4.0),          # a K giant born at 1.04 solar masses (was 21.8 AU)
    (1.73, 3.0, 12.0),
])
def test_an_evolved_star_gets_the_snow_line_of_its_progenitors_birth_not_of_today(mass, low, high):
    assert low < fsl.formation_snow_line_au(mass, 0.0) < high


def _star_like(initial_mass_sol, mass_kg, fe_h=0.0, yerkes="VII"):
    return types.SimpleNamespace(initial_mass_sol=initial_mass_sol, mass=mass_kg, metallicity_fe_h=fe_h, yerkes_class=yerkes)


def test_a_star_with_no_stored_initial_mass_uses_the_mass_its_planets_formed_round():
    star = types.SimpleNamespace(initial_mass_sol=None, mass=1.0 * constants.SOLAR_MASS_TO_KG, metallicity_fe_h=0.0, yerkes_class="V")
    assert fsl.initial_mass_sol_of(star) == pytest.approx(1.0)
    assert fsl.initial_mass_sol_of(_star_like(4.36, 0.6 * constants.SOLAR_MASS_TO_KG)) == 4.36


def test_the_drawn_line_is_inside_the_scatter_cut_of_the_birth_position_and_draws_from_its_own_stream():
    star = _star_like(1.0, constants.SOLAR_MASS_TO_KG, yerkes="V")
    center = fsl.formation_snow_line_au(1.0, 0.0)
    lines = [fsl.draw_for_star(star, draw.Stream(f"t:{index}")) for index in range(200)]
    limit = tuning.SNOW_LINE_SCATTER_DEX * tuning.SNOW_LINE_SCATTER_MAX_SIGMA
    assert all(abs(math.log10(line / center)) <= limit + 1e-9 for line in lines)
    assert lines[0] == fsl.draw_for_star(star, draw.Stream("t:0"))
    draw.set_run_seed(7)
    first = draw.random()
    draw.set_run_seed(7)
    fsl.draw_for_star(star, draw.Stream("t:0"))
    assert draw.random() == first


def test_a_generated_white_dwarf_stores_the_line_of_its_progenitors_disk():
    seen = 0
    for seed in range(200):
        with draw.bound(seed):
            star = Star(SystemConfig())
        if star.yerkes_class not in ("VII", "D"):
            continue
        seen += 1
        center = fsl.formation_snow_line_au(fsl.initial_mass_sol_of(star), star.metallicity_fe_h)
        assert abs(math.log10(star.snow_line_au / center)) <= tuning.SNOW_LINE_SCATTER_DEX * tuning.SNOW_LINE_SCATTER_MAX_SIGMA + 1e-9
        assert star.snow_line_au > 0.2, "not the 0.02 AU of a white dwarf's own light"
        if seen >= 4:
            break
    assert seen


# --- GEN.1678: the water line ---------------------------------------------------------------------------------------------------------------

PAPER_TABLE_1 = {1.0: (3.3, 1.6), 0.8: (2.6, 1.2), 0.6: (1.9, 0.9), 0.4: (1.2, 0.6), 0.2: (0.5, 0.2)}
"""dict: Mulders et al. 2015 ([R3193]) Table 1, the one-sigma middle snow line (au) with small (570 cm2/g) and large (20 cm2/g) grains, by stellar mass."""


def test_the_viscous_line_is_the_papers_at_both_grain_sizes():
    for mass, (small, large) in PAPER_TABLE_1.items():
        assert fsl.viscous_snow_line_au(mass, 570.0) == pytest.approx(small, rel=0.12, abs=0.08)
        assert fsl.viscous_snow_line_au(mass, 20.0) == pytest.approx(large, rel=0.12, abs=0.08)
        assert fsl.viscous_snow_line_au(mass, 20.0) < fsl.viscous_snow_line_au(mass) < fsl.viscous_snow_line_au(mass, 570.0)


def test_the_viscous_line_scales_with_the_papers_exponent_of_stellar_mass():
    assert fsl.viscous_snow_line_au(1.0) / fsl.viscous_snow_line_au(0.5) == pytest.approx(2.0 ** 1.14, rel=0.01)       # [R3193] sec. 2.4: R ~ M^1.14 (Eq. 2 and 4)
    assert fsl.viscous_snow_line_au(0.001) == fsl.viscous_snow_line_au(fsl.LUPUS_LOWEST_MASS_SOL)


def _lone(mass_sol, line_au, fe_h=0.0):
    return types.SimpleNamespace(initial_mass_sol=mass_sol, mass=mass_sol * constants.SOLAR_MASS_TO_KG, metallicity_fe_h=fe_h, snow_line_au=line_au, yerkes_class="V")


def test_the_water_line_moves_the_stars_own_line_by_the_viscous_over_the_luminosity_line_and_keeps_its_scatter():
    mean = fsl.formation_snow_line_au(0.4, 0.0)
    factor = fsl.viscous_snow_line_au(0.4) / mean
    assert fsl.water_snow_line_au(_lone(0.4, mean)) == pytest.approx(mean * factor)
    assert fsl.water_snow_line_au(_lone(0.4, 2.0 * mean)) == pytest.approx(2.0 * mean * factor)      # the star's scatter stays
    assert fsl.water_snow_line_au(_lone(0.4, None)) is None


def test_the_water_line_of_a_pair_adds_its_members_lines_as_the_stored_line_does():
    first, second = _lone(1.0, 1.6), _lone(0.3, 0.5)
    pair = types.SimpleNamespace(stars=[first, second])
    assert fsl.water_snow_line_au(pair) == pytest.approx(math.hypot(fsl.water_snow_line_au(first), fsl.water_snow_line_au(second)))


def test_the_delivered_water_reads_the_water_line():
    from planetgen.physics import bulk_draw
    star = _lone(0.4, fsl.formation_snow_line_au(0.4, 0.0))
    body = types.SimpleNamespace(star=star, formation_distance_au=0.5, distance=0.5)
    assert bulk_draw._snow_ratio(body) == pytest.approx(0.5 / fsl.water_snow_line_au(star))
    assert bulk_draw._snow_ratio(body) != pytest.approx(0.5 / star.snow_line_au)


def test_a_body_is_icy_by_the_water_line_not_the_luminosity_line():
    """GEN.1679: the icy/rocky test and the stored 'formed beyond the snow line' flag read the line the water is read against."""
    from planetgen.physics import bulk_draw, planets
    star = _lone(1.0, fsl.formation_snow_line_au(1.0, 0.0))
    water_line = fsl.water_snow_line_au(star)
    assert water_line > star.snow_line_au          # a Sun-like star: the viscous line lies outside the luminosity line
    between = types.SimpleNamespace(star=star, formation_distance_au=0.5 * (water_line + star.snow_line_au), distance=1.0)
    beyond = types.SimpleNamespace(star=star, formation_distance_au=1.01 * water_line, distance=1.0)
    assert bulk_draw.origin_of(between) == bulk_draw.ROCKY and not planets.formed_beyond_snow_line(between)
    assert bulk_draw.origin_of(beyond) == bulk_draw.ICY and planets.formed_beyond_snow_line(beyond)
    dim = _lone(0.2, fsl.formation_snow_line_au(0.2, 0.0))          # an M dwarf: the viscous line lies inside the luminosity line
    assert fsl.water_snow_line_au(dim) < dim.snow_line_au
    inside = types.SimpleNamespace(star=dim, formation_distance_au=0.5 * (fsl.water_snow_line_au(dim) + dim.snow_line_au), distance=1.0)
    assert bulk_draw.origin_of(inside) == bulk_draw.ICY
