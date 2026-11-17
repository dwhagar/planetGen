# src/tests/test_coronal_wind.py

"""GEN.1625: one wind for a cool dwarf, from its X-ray luminosity (A to K) or a measured mass loss per area (M), read by the heliosphere, the black holes it feeds and
the pressure on a planet's magnetosphere; a white dwarf has none."""

import math
from types import SimpleNamespace

import pytest

from planetgen import tuning
from planetgen.generation.star import Star
from planetgen.physics import constants, magnetism, stellar_wind

SOLAR_MASS = 1.989e30
SUN_LOG_LX_LBOL = math.log10(constants.SOLAR_X_RAY_LUMINOSITY_W / constants.SOLAR_LUMINOSITY)
SUN_LOSS_KG_S = 2e-14 * SOLAR_MASS / constants.SECONDS_PER_YEAR


def _wind(mass_sol, luminosity_sol, radius_sol, star_type, yerkes_class="V", log_lx_lbol=SUN_LOG_LX_LBOL):
    return stellar_wind.wind(mass_sol * SOLAR_MASS, luminosity_sol * constants.SOLAR_LUMINOSITY, radius_sol * 695700.0, star_type, yerkes_class, log_lx_lbol)


def _star(*args, **kwargs):
    mass_sol, luminosity_sol, radius_sol, star_type = args
    return SimpleNamespace(mass=mass_sol * SOLAR_MASS, luminosity=luminosity_sol * constants.SOLAR_LUMINOSITY, radius=radius_sol * 695700.0, type=star_type,
                           yerkes_class=kwargs.get("yerkes_class", "V"), log_lx_lbol=kwargs.get("log_lx_lbol", SUN_LOG_LX_LBOL))


def test_the_sun_has_the_suns_wind():
    loss, speed = _wind(1.0, 1.0, 1.0, "G2V")
    assert loss == pytest.approx(SUN_LOSS_KG_S, rel=0.01)
    assert speed == pytest.approx(constants.SOLAR_WIND_VELOCITY, rel=0.01)


def test_an_a_to_k_dwarfs_wind_rises_with_its_surface_x_ray_flux_to_the_measured_power():
    quiet, active = _wind(0.8, 0.3, 0.8, "K1V", log_lx_lbol=-5.5)[0], _wind(0.8, 0.3, 0.8, "K1V", log_lx_lbol=-5.0)[0]
    assert active / quiet == pytest.approx(10 ** (0.5 * tuning.CORONAL_WIND_X_RAY_EXPONENT), rel=1e-9)


def test_the_flux_is_held_at_the_dividing_line_where_the_measured_winds_turn_weak():
    saturated = _wind(0.8, 0.3, 0.8, "K1V", log_lx_lbol=-3.0)[0]
    assert _wind(0.8, 0.3, 0.8, "K1V", log_lx_lbol=-2.0)[0] == pytest.approx(saturated)
    solar_flux_w_m2 = constants.SOLAR_X_RAY_LUMINOSITY_W / (4 * math.pi * constants.SOLAR_RADIUS_M ** 2)
    cap = tuning.CORONAL_WIND_X_RAY_FLUX_CAP_ERG_CM2_S * constants.ERG_CM2_S_TO_W_M2 / solar_flux_w_m2
    assert saturated == pytest.approx(SUN_LOSS_KG_S * 0.8 ** 2 * cap ** tuning.CORONAL_WIND_X_RAY_EXPONENT, rel=0.01)


def test_an_m_dwarfs_wind_is_a_measured_mass_loss_per_area_and_not_stronger_than_proxima_centauris_limit():
    """It was 2.8 times the Sun's for an M5 dwarf; [R5302] puts Proxima Centauri (M5.5) under 0.2 of it."""
    loss, speed = _wind(0.15, 0.003, 0.2, "M5V", log_lx_lbol=-3.0)
    assert loss / SUN_LOSS_KG_S == pytest.approx(tuning.M_DWARF_WIND_PER_AREA_SOLAR * 0.2 ** 2, rel=0.01)
    assert loss / SUN_LOSS_KG_S < 0.2
    assert 2e5 < speed < 6e5
    assert _wind(0.15, 0.003, 0.2, "M5V", log_lx_lbol=-5.0)[0] == pytest.approx(loss)        # it does not follow the X-ray flux


def test_a_white_dwarf_has_no_wind_and_a_star_with_a_corona_needs_its_x_ray_luminosity():
    assert _wind(0.6, 1e-3, 0.012, "DA", yerkes_class="VII", log_lx_lbol=None) == (0.0, 0.0)
    assert _wind(0.6, 1e-3, 0.012, "DA", yerkes_class="D", log_lx_lbol=None)[0] == 0.0
    with pytest.raises(ValueError):
        _wind(1.0, 1.0, 1.0, "G2V", log_lx_lbol=None)


def test_the_planets_wind_pressure_is_the_same_wind_as_the_heliospheres():
    sun = _star(1.0, 1.0, 1.0, "G2V")
    assert magnetism.wind_pressure_pa(sun, 1.0) == pytest.approx(magnetism.EARTH_WIND_PA, rel=0.01)
    m_dwarf = _star(0.15, 0.003, 0.2, "M5V", log_lx_lbol=-3.0)
    loss, speed = stellar_wind.wind(m_dwarf.mass, m_dwarf.luminosity, m_dwarf.radius, m_dwarf.type, m_dwarf.yerkes_class, m_dwarf.log_lx_lbol)
    sun_loss, sun_speed = stellar_wind.wind(sun.mass, sun.luminosity, sun.radius, sun.type, sun.yerkes_class, sun.log_lx_lbol)
    assert magnetism.wind_pressure_pa(m_dwarf, 1.0) / magnetism.wind_pressure_pa(sun, 1.0) == pytest.approx((loss * speed) / (sun_loss * sun_speed), rel=1e-9)
    assert magnetism.wind_pressure_pa(m_dwarf, 0.05) / magnetism.wind_pressure_pa(m_dwarf, 0.1) == pytest.approx(4.0)


def test_a_white_dwarfs_planet_sits_in_the_interstellar_pressure_and_its_bubble_is_the_floor_size():
    white_dwarf = _star(0.6, 1e-3, 0.012, "DA", yerkes_class="VII", log_lx_lbol=None)
    assert magnetism.wind_pressure_pa(white_dwarf, 1.0) == constants.ISM_PRESSURE      # it used to get the Earth's wind
    radius_au = Star._calculate_heliosphere_radius_static(white_dwarf.mass, white_dwarf.luminosity, white_dwarf.radius, white_dwarf.type, "VII", None)
    assert 0.0 <= radius_au < 1e-3


def test_a_star_with_no_light_has_no_wind_and_so_no_standoff():
    """A starless host (a rogue planet's) and a black hole nothing feeds give no light and no wind; the interstellar floor is for a quiet star, not for no star at all (GEN.1111)."""
    nothing = SimpleNamespace(mass=0.0, luminosity=0.0, radius=0.0, type="", yerkes_class="", log_lx_lbol=None)
    assert magnetism.wind_pressure_pa(nothing, 1.0) == 0.0
    planet = SimpleNamespace(magnetic_moment_a_m2=8.0e22, star=nothing, radius=6371.0)
    assert magnetism.magnetopause_rp(planet, 1.0) is None
    quiet = _star(0.3, 0.001, 0.3, "M5V", log_lx_lbol=-9.0)
    assert magnetism.wind_pressure_pa(quiet, 1.0) >= constants.ISM_PRESSURE > 0.0
