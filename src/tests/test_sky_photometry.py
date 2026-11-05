"""VIEW.20: photometry for the system sky on arrays (`planetgen/physics/sky_photometry.py`): apparent magnitude, the
bolometric correction for any temperature, and the galactic dust along a line of sight."""

import math

import numpy as np
import pytest

from planetgen.physics import photometry, sky_photometry as sky

AU_PC = 1.0 / 206264.806


# --- magnitude against published columns ------------------------------------------------------------------------------

def test_the_sun_is_minus_26_75_at_1_au_and_4_82_at_10_pc():
    assert sky.apparent_magnitude(1.0, 5772.0, AU_PC) == pytest.approx(-26.75, abs=0.02)         # published: -26.74
    assert sky.apparent_magnitude(1.0, 5772.0, 10.0) == pytest.approx(4.82, abs=0.01)            # published M_V 4.81 to 4.83


@pytest.mark.parametrize("name, lum, teff, distance, published", [
    ("Sirius A", 25.4, 9940.0, 2.637, -1.46),
    ("Alpha Centauri A", 1.519, 5790.0, 1.34, -0.01),
    ("Procyon A", 6.93, 6530.0, 3.51, 0.34),
])
def test_a_nearby_star_comes_out_within_a_tenth_of_a_magnitude_of_its_catalogue_value(name, lum, teff, distance, published):
    assert sky.apparent_magnitude(lum, teff, distance) == pytest.approx(published, abs=0.1), name


def test_a_hundred_times_the_distance_is_ten_magnitudes_fainter_and_dust_adds_straight_on():
    near = sky.apparent_magnitude(3.0, 6000.0, 5.0)
    assert sky.apparent_magnitude(3.0, 6000.0, 500.0) == pytest.approx(near + 10.0)
    assert sky.apparent_magnitude(3.0, 6000.0, 5.0, a_v=0.7) == pytest.approx(near + 0.7)


def test_arrays_broadcast_and_agree_with_one_star_at_a_time():
    lum = np.array([0.001, 1.0, 1e4])
    teff = np.array([3100.0, 5772.0, 25000.0])
    dist = np.array([2.0, 30.0, 900.0])
    many = sky.apparent_magnitude(lum, teff, dist)
    assert many.shape == (3,)
    for i in range(3):
        assert many[i] == pytest.approx(float(sky.apparent_magnitude(lum[i], teff[i], dist[i])))


@pytest.mark.parametrize("bad", [0.0, -1.0, math.nan, math.inf])
def test_a_distance_luminosity_or_temperature_that_is_not_positive_and_finite_raises(bad):
    with pytest.raises(ValueError):
        sky.apparent_magnitude(1.0, 5772.0, bad)
    with pytest.raises(ValueError):
        sky.apparent_magnitude(bad, 5772.0, 10.0)
    with pytest.raises(ValueError):
        sky.apparent_magnitude(1.0, bad, 10.0)


# --- the bolometric correction for any temperature -----------------------------------------------------------------

def test_inside_the_fit_it_is_the_fit_of_view_9():
    for temperature in (4100.0, 4500.0, 5772.0, 9000.0, 20000.0, 45000.0, 50000.0):
        assert sky.bolometric_correction(temperature) == pytest.approx(photometry.bolometric_correction_v(temperature))


def test_the_cool_rule_is_the_main_sequence_table_at_its_rows_and_between_them():
    for temperature, expected in sky.COOL_BC_TABLE_K_MAG:
        assert sky.bolometric_correction(temperature - 1e-9 if temperature == sky.COOL_BELOW_K else temperature) == \
            pytest.approx(expected, abs=0.01 if temperature == sky.COOL_BELOW_K else 1e-9)
    assert sky.bolometric_correction(3060.0) == pytest.approx(-3.11)                           # M5V
    middle = sky.bolometric_correction(3160.0)                                                # between M4.5V and M4V
    assert -2.84 < middle < -2.51


def test_the_cool_table_is_what_the_fit_gets_wrong_below_3500_k():
    assert photometry.bolometric_correction_v(3500.0) == pytest.approx(-2.31, abs=0.01)         # the fit
    assert sky.bolometric_correction(3500.0) == pytest.approx(-1.73, abs=0.05)                   # the dwarfs: 0.6 mag brighter in V
    assert sky.bolometric_correction(2000.0) == pytest.approx(-6.62)                              # held at the last row, not -24


def test_it_steps_by_under_a_tenth_of_a_magnitude_where_the_table_hands_over_to_the_fit():
    below = sky.bolometric_correction(sky.COOL_BELOW_K - 0.01)
    above = sky.bolometric_correction(sky.COOL_BELOW_K + 0.01)
    assert abs(above - below) < 0.1


def test_a_cooler_star_always_has_the_larger_correction_in_magnitude_from_4000_k_down():
    temperatures = np.linspace(2060.0, 4099.0, 200)
    assert np.all(np.diff(sky.bolometric_correction(temperatures)) >= 0)                        # warmer: less negative


def test_over_50000_k_it_follows_the_blackbody_tail_continuously():
    top = float(sky.bolometric_correction(50000.0))
    assert float(sky.bolometric_correction(50000.01)) == pytest.approx(top, abs=1e-3)
    assert np.all(np.diff(sky.bolometric_correction(np.geomspace(50000.0, 1e6, 100))) < 0)       # hotter: more negative


def test_the_tail_falls_toward_7_5_magnitudes_a_decade_as_the_v_band_goes_down_the_rayleigh_jeans_side():
    """The V flux of a blackbody goes as T and its bolometric flux as T^4: 2.5 x 3 = 7.5 mag per decade in the limit."""
    def slope(low):
        return float(sky.bolometric_correction(low) - sky.bolometric_correction(low * 10.0))
    assert 7.1 < slope(50000.0) < 7.3 and slope(1e6) == pytest.approx(7.5, abs=0.02)
    assert slope(1e8) == pytest.approx(7.5, abs=1e-3)


def test_an_unusable_temperature_raises():
    for bad in (0.0, -5.0, math.nan):
        with pytest.raises(ValueError):
            sky.bolometric_correction(bad)


# --- galactic dust ---------------------------------------------------------------------------------------------------

def test_a_line_straight_up_from_the_solar_circle_is_the_sech_squared_integral():
    h = sky.DUST_HEIGHT_PC + sky.DUST_HEIGHT_FLARE_PC_PER_KPC * (8000.0 - sky.DUST_FLARE_START_PC) / 1000.0
    expected = 1.0 / 1000.0 * h * math.tanh(1000.0 / h)
    assert sky.galactic_extinction_av([8000.0, 0.0, 0.0], [0.0, 0.0, 1.0], 1000.0) == pytest.approx(expected, rel=1e-4)


def test_a_line_outwards_in_the_plane_is_the_exponential_integral():
    expected = 1.0 / 1000.0 * sky.DUST_SCALE_LENGTH_PC * (1.0 - math.exp(-1000.0 / sky.DUST_SCALE_LENGTH_PC))
    assert sky.galactic_extinction_av([8000.0, 0.0, 0.0], [1.0, 0.0, 0.0], 1000.0) == pytest.approx(expected, rel=1e-4)


def test_the_local_plane_dims_one_magnitude_per_kiloparsec_for_a_short_chord():
    assert sky.galactic_extinction_av([8000.0, 0.0, 0.0], [0.0, 1.0, 0.0], 100.0) == pytest.approx(0.1, rel=0.01)


def test_the_dust_is_the_same_whichever_way_the_direction_is_scaled_and_zero_at_zero_distance():
    a = sky.galactic_extinction_av([8000.0, 0.0, 0.0], [1.0, 2.0, 3.0], 700.0)
    assert sky.galactic_extinction_av([8000.0, 0.0, 0.0], [5.0, 10.0, 15.0], 700.0) == pytest.approx(a)
    assert sky.galactic_extinction_av([8000.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.0) == 0.0


def test_many_lines_at_once_agree_with_one_at_a_time():
    rng = np.random.default_rng(7)
    directions = rng.normal(size=(5, 3))
    distances = rng.uniform(50.0, 6000.0, 5)
    many = sky.galactic_extinction_av([7000.0, 500.0, 20.0], directions, distances)
    assert many.shape == (5,)
    for i in range(5):
        assert many[i] == pytest.approx(sky.galactic_extinction_av([7000.0, 500.0, 20.0], directions[i], distances[i]))


def test_the_midpoint_steps_are_fine_enough_for_a_long_grazing_line():
    coarse = sky.galactic_extinction_av([8000.0, 0.0, 60.0], [-1.0, 0.0, 0.0], 25000.0)
    fine = sky.galactic_extinction_av([8000.0, 0.0, 60.0], [-1.0, 0.0, 0.0], 25000.0, steps=8192)
    assert coarse == pytest.approx(fine, rel=0.005)


def test_a_denser_law_or_a_smaller_galaxy_scales_the_extinction():
    base = sky.galactic_extinction_av([8000.0, 0.0, 0.0], [0.0, 1.0, 0.0], 500.0)
    assert sky.galactic_extinction_av([8000.0, 0.0, 0.0], [0.0, 1.0, 0.0], 500.0, av_per_kpc_at_sun=2.0) == pytest.approx(2 * base)


def test_a_zero_direction_or_a_negative_distance_raises():
    with pytest.raises(ValueError):
        sky.galactic_extinction_av([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 10.0)
    with pytest.raises(ValueError):
        sky.galactic_extinction_av([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], -10.0)


def test_a_star_seen_through_the_disk_is_fainter_than_the_same_star_seen_out_of_it():
    inside = sky.galactic_extinction_av([8000.0, 0.0, 0.0], [1.0, 0.0, 0.0], 3000.0)
    above = sky.galactic_extinction_av([8000.0, 0.0, 0.0], [0.0, 0.0, 1.0], 3000.0)
    assert inside > 2.5 * above
