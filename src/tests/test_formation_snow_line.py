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
