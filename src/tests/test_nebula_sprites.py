# tests/test_nebula_sprites.py

"""VIEW.28: the sprites of nebulae and supernova remnants: angular size, surface brightness and extinction, against hand cases."""

import math

import numpy as np
import pytest

from planetgen.physics import constants, light_travel, nebula_sprites as ns

LY = constants.LY_TO_M
YEAR = constants.SECONDS_PER_YEAR


def test_a_ball_seen_from_outside_has_the_angular_radius_of_its_silhouette():
    theta = ns.angular_radius_rad([10.0, 10.0], [100.0, 10000.0])
    assert theta[0] == pytest.approx(math.asin(0.1)) and theta[1] == pytest.approx(math.atan(1e-3), rel=1e-6)


def test_an_observer_inside_the_ball_sees_it_cover_the_sky():
    theta = ns.angular_radius_rad([50.0, 50.0], [20.0, 50.0])
    assert list(theta) == [pytest.approx(math.pi), pytest.approx(math.pi)]


def test_the_chord_extinction_is_the_diameter_at_the_centre_and_none_at_the_edge():
    dimming = ns.chord_extinction_mag(2.0, 10.0, [0.0, 5.0, 10.0, 12.0])
    assert dimming[0] == pytest.approx(2.0) and dimming[1] == pytest.approx(2.0 * math.sqrt(0.75)) and dimming[2] == 0.0 and dimming[3] == 0.0


def test_a_sprite_inside_its_ball_dims_by_half_the_diameter():
    assert list(ns.sprite_extinction_mag([2.0, 2.0], [10.0, 10.0], [5.0, 50.0])) == [1.0, 2.0]


def test_the_emission_measure_is_density_squared_over_the_mean_chord():
    pc_in_ly = constants.PARSEC_M / LY
    measure = ns.emission_measure_cm6_pc(100.0, pc_in_ly)         # a ball of radius 1 pc
    assert measure == pytest.approx(100.0 ** 2 * 4.0 / 3.0)


def test_a_rayleigh_of_hydrogen_alpha_is_about_thirty_magnitudes_per_square_arcsecond_and_ten_times_that_is_2_5_brighter():
    assert ns.rayleigh_to_v_mag(1.0) == pytest.approx(30.4, abs=0.1)
    assert ns.rayleigh_to_v_mag(1.0) - ns.rayleigh_to_v_mag(10.0) == pytest.approx(2.5)


def test_a_bright_h_ii_region_is_far_brighter_than_the_naked_eye_limit_and_a_diffuse_gas_is_not_visible():
    bright = ns.emission_surface_brightness(1.0e3, 3.0)           # a compact H II region
    faint = ns.emission_surface_brightness(0.3, 100.0)            # diffuse ionised gas
    assert bright < 18.0 and faint > 25.0


def test_a_larger_or_denser_gas_is_brighter_and_a_shock_shell_is_fainter_than_a_photoionised_gas():
    assert ns.emission_surface_brightness(1.0e3, 3.0) < ns.emission_surface_brightness(1.0e2, 3.0)
    assert ns.emission_surface_brightness(100.0, 30.0) < ns.emission_surface_brightness(100.0, 3.0)
    assert ns.remnant_surface_brightness(1.0, 50.0) > ns.emission_surface_brightness(1.0, 50.0)      # the same gas, a thin shell at half the light
    assert ns.remnant_surface_brightness(1.0, 50.0) == pytest.approx(ns.rayleigh_to_v_mag(
        ns.SHELL_HALPHA_SHARE * (ns.SHELL_COMPRESSION * 1.0) ** 2 * ns.SHELL_THICKNESS_OVER_RADIUS * 50.0 * ns.PC_PER_LY / ns.EM_PER_RAYLEIGH))


def test_a_bright_reflection_nebula_is_about_twenty_two_or_twenty_three():
    assert 21.5 < ns.reflection_surface_brightness(10.0, 1.5) < 24.0


def test_dark_and_diffuse_clouds_and_plerions_have_no_light_of_their_own():
    family = np.array(["dark", "diffuse", "plerion", "emission"])
    brightness = ns.surface_brightness(family, [1e4] * 4, [10.0] * 4, [10.0] * 4)
    assert list(np.isnan(brightness)) == [True, True, True, False]


def test_the_total_magnitude_is_the_surface_brightness_less_the_area_term():
    one_arcmin = 60.0 * math.pi / (180.0 * 3600.0)
    assert ns.total_magnitude(20.0, one_arcmin) == pytest.approx(20.0 - 2.5 * math.log10(math.pi * 3600.0))


def test_the_foreground_dust_dims_the_surface_brightness_and_the_limit_decides_visibility():
    near = ns.sprite_data(np.array(["emission"]), [3.0], [100.0], [1e3], [0.5])
    dimmed = ns.sprite_data(np.array(["emission"]), [3.0], [100.0], [1e3], [0.5], foreground_av=20.0)
    assert dimmed["surface_brightness"][0] == pytest.approx(near["surface_brightness"][0] + 20.0)
    assert near["visible"][0] and not dimmed["visible"][0]


def test_a_remnant_younger_than_its_light_time_is_not_seen_from_the_retarded_geometry():
    """A 500-year-old remnant 1,000 light years away: its light left 1,000 years before the sky's time, before the explosion. At 100 light years it has arrived."""
    offsets = np.array([[1000.0 * LY, 0.0, 0.0], [100.0 * LY, 0.0, 0.0]])
    apparent, tau = light_travel.retarded_offsets(offsets, np.zeros((2, 3)), (0.0, 0.0, 0.0), 0.0)
    distance_ly = np.sqrt((apparent ** 2).sum(axis=1)) / LY
    sprites = ns.sprite_data(np.array(["supernova-remnant"] * 2), [20.0, 20.0], distance_ly, [1.0, 1.0], [0.05, 0.05],
                             age_years=[500.0, 500.0], light_time_s=tau)
    assert list(sprites["seen"]) == [False, True] and list(sprites["visible"]) == [False, False]


def test_a_nebula_is_always_seen_and_a_plerion_remnant_follows_the_light_time_too():
    sprites = ns.sprite_data(np.array(["emission", "plerion"]), [3.0, 3.0], [100.0, 100.0], [1e3, 1.0], [0.5, 0.05],
                             age_years=[float("nan"), 50.0], light_time_s=[1e3 * YEAR, 1e3 * YEAR])
    assert list(sprites["seen"]) == [True, False]
