# tests/test_band_winds.py

"""
GEN.1110: the mean surface winds by day and by night in each latitude band (`planetgen.physics.band_winds`): the Earth, Mars, Venus and Titan to their order of magnitude, the
heat-engine speed, the thermal wind and its smooth minimum, the Hadley edge, the regimes, the fade in thin air, the locked world's day-to-night flow and the inputs refused.
"""

import math

import pytest

from planetgen.physics import band_winds as bw

B = bw.Band
EARTH_OMEGA = 7.2921e-5
EARTH = dict(bands=[B(0, 300, 293), B(25, 295, 288), B(50, 284, 279), B(80, 260, 255)], pressure_pa=101325.0, mean_molar_mass_kg_mol=0.02897,
             radius_m=6.371e6, rotation_rate_rad_s=EARTH_OMEGA)
MARS = dict(bands=[B(0, 260, 190), B(30, 245, 180), B(60, 215, 165), B(85, 175, 150)], pressure_pa=600.0, mean_molar_mass_kg_mol=0.04334,
            radius_m=3.39e6, rotation_rate_rad_s=7.088e-5)
VENUS = dict(bands=[B(0, 737.5, 736.5), B(45, 736.5, 735.5), B(85, 735.0, 734.0)], pressure_pa=9.2e6, mean_molar_mass_kg_mol=0.04345, radius_m=6.05e6,
             rotation_rate_rad_s=-2.99e-7)
TITAN = dict(bands=[B(0, 94, 93.5), B(45, 92, 91.5), B(85, 90, 89.5)], pressure_pa=1.467e5, mean_molar_mass_kg_mol=0.0284, radius_m=2.575e6,
             rotation_rate_rad_s=4.56e-6)


def test_the_earth_has_surface_winds_of_a_few_to_ten_metres_a_second_calmer_at_night_and_weakest_at_the_pole():
    winds = bw.band_winds(**EARTH)
    assert all(2.0 < w.day_m_s < 15.0 for w in winds)
    assert 4.0 < winds[1].day_m_s < 9.0 and 3.0 < winds[2].day_m_s < 8.0             # the mid-latitudes: 5 to 8 m/s
    assert all(w.night_m_s < w.day_m_s for w in winds)
    assert 0.4 < sum(w.night_m_s for w in winds) / sum(w.day_m_s for w in winds) < 0.7
    assert winds[3].day_m_s < winds[1].day_m_s


def test_the_earths_bands_are_trades_westerlies_and_polar_winds():
    assert [w.regime for w in bw.band_winds(**EARTH)] == [bw.TRADES, bw.WESTERLIES, bw.WESTERLIES, bw.POLAR]


def test_mars_is_a_few_metres_a_second_and_venus_and_titan_are_calm():
    mars = bw.band_winds(**MARS)
    assert all(1.0 < w.day_m_s < 12.0 for w in mars)
    for case in (VENUS, TITAN):
        assert all(w.day_m_s < 3.0 for w in bw.band_winds(**case))
    assert all(w.day_m_s < 1.0 for w in bw.band_winds(**VENUS))
    assert all(w.regime == bw.OVERTURNING for w in bw.band_winds(**VENUS)) and all(w.regime == bw.OVERTURNING for w in bw.band_winds(**TITAN))


def test_the_heat_engine_speed_is_the_carnot_share_of_the_pressure_work():
    assert bw.heat_engine_speed_m_s(10.0, 288.0, 287.0) == pytest.approx(math.sqrt(10.0 / 288.0 * 287.0 * 10.0))
    assert bw.heat_engine_speed_m_s(-10.0, 288.0, 287.0) == bw.heat_engine_speed_m_s(10.0, 288.0, 287.0)
    assert bw.heat_engine_speed_m_s(0.0, 288.0, 287.0) == 0.0
    # doubling the contrast at a fixed mean doubles the speed (dT^2 / T under the root)
    assert bw.heat_engine_speed_m_s(20.0, 288.0, 287.0) == pytest.approx(2.0 * bw.heat_engine_speed_m_s(10.0, 288.0, 287.0))


def test_the_thermal_wind_falls_with_spin_and_with_latitude_and_is_infinite_at_the_equator():
    args = (40.0, 287.0, 6.371e6)
    assert math.isinf(bw.thermal_wind_m_s(0.0, *args, EARTH_OMEGA)) and math.isinf(bw.thermal_wind_m_s(45.0, *args, 0.0))
    assert bw.thermal_wind_m_s(45.0, *args, EARTH_OMEGA) == pytest.approx(2.0 * bw.thermal_wind_m_s(45.0, *args, 2.0 * EARTH_OMEGA))
    assert bw.thermal_wind_m_s(20.0, *args, EARTH_OMEGA) > bw.thermal_wind_m_s(60.0, *args, EARTH_OMEGA)
    # the Earth's mid-latitude shear is a few tens of metres a second
    assert 15.0 < bw.thermal_wind_m_s(45.0, *args, EARTH_OMEGA) < 50.0


def test_the_hadley_cell_reaches_the_pole_for_a_thermal_rossby_number_of_one():
    assert bw.hadley_edge_deg(1.0) == 90.0 and bw.hadley_edge_deg(5.0) == 90.0
    assert bw.hadley_edge_deg(0.0) == 0.0
    assert bw.hadley_edge_deg(0.12) == pytest.approx(math.degrees(math.sqrt(5.0 / 3.0 * 0.12)))
    ro = bw.thermal_rossby(39.0, 287.0, 6.371e6, EARTH_OMEGA)
    assert 0.03 < ro < 0.2 and 15.0 < bw.hadley_edge_deg(ro) < 40.0               # the Earth's cell ends near 30 degrees
    assert math.isinf(bw.thermal_rossby(39.0, 287.0, 6.371e6, 0.0))


def test_the_surface_wind_fades_in_thin_air_to_a_half_at_a_kilopascal_and_to_nothing_with_no_air():
    assert bw.pressure_fade(bw.PRESSURE_FADE_PA) == pytest.approx(0.5)
    assert bw.pressure_fade(1.0e5) > 0.99 and bw.pressure_fade(10.0) < 1e-3 and bw.pressure_fade(0.0) == 0.0
    fades = [bw.pressure_fade(p) for p in (1.0, 10.0, 100.0, 1000.0, 1e4, 1e5)]
    assert fades == sorted(fades)
    thin = bw.band_winds(**{**EARTH, "pressure_pa": 100.0})
    assert all(t.day_m_s < 0.05 * w.day_m_s for t, w in zip(thin, bw.band_winds(**EARTH)))
    assert all(w.day_m_s == 0.0 for w in bw.band_winds(**{**EARTH, "pressure_pa": 0.0}))


def test_a_larger_day_night_difference_is_a_stronger_wind_and_none_leaves_the_zonal_flow():
    flat = [B(0, 290, 290), B(60, 270, 270)]
    swing = [B(0, 310, 270), B(60, 280, 260)]
    calm = bw.band_winds(flat, 101325.0, 0.02897, 6.371e6, EARTH_OMEGA)
    windy = bw.band_winds(swing, 101325.0, 0.02897, 6.371e6, EARTH_OMEGA)
    assert all(w.day_night_m_s == 0.0 for w in calm) and all(w.day_night_m_s > 0.0 for w in windy)
    assert windy[0].day_m_s > calm[0].day_m_s


def test_the_diurnal_flow_is_cut_at_high_latitude_by_the_inertial_limit():
    swing = [B(0, 310, 270), B(89.9, 310, 270)]
    equator, pole = bw.band_winds(swing, 101325.0, 0.02897, 6.371e6, EARTH_OMEGA)
    assert pole.day_night_m_s == pytest.approx(equator.day_night_m_s / math.sqrt(5.0), rel=0.01)


def test_a_locked_world_has_one_flow_from_the_day_side_to_the_night_side_and_a_calmer_night():
    locked = bw.band_winds([B(0, 330, 200), B(45, 300, 215), B(85, 260, 230)], 1.0e5, 0.029, 6.4e6, 1.0e-6, locked=True)
    assert all(w.regime == bw.LOCKED for w in locked)
    assert locked[0].day_m_s > locked[1].day_m_s > locked[2].day_m_s > 0.0
    assert all(w.night_m_s == pytest.approx(bw.NIGHT_COUPLING * w.day_m_s) for w in locked)
    assert 5.0 < locked[0].day_m_s < 80.0


def test_a_world_that_does_not_turn_is_all_heat_engine():
    still = bw.band_winds(**{**EARTH, "rotation_rate_rad_s": 0.0})
    assert all(w.regime == bw.OVERTURNING for w in still)
    assert still[1].zonal_aloft_m_s == pytest.approx(still[2].zonal_aloft_m_s)


def test_speeds_are_unchanged_by_the_sense_of_the_spin_and_the_order_of_the_bands():
    forward = bw.band_winds(**EARTH)
    backward = bw.band_winds(**{**EARTH, "rotation_rate_rad_s": -EARTH_OMEGA})
    assert [w.day_m_s for w in forward] == pytest.approx([w.day_m_s for w in backward])
    shuffled = bw.band_winds(**{**EARTH, "bands": list(reversed(EARTH["bands"]))})
    assert [w.day_m_s for w in shuffled] == pytest.approx([w.day_m_s for w in reversed(forward)])


def test_inputs_the_model_cannot_use_are_refused():
    with pytest.raises(bw.BandWindError):
        bw.band_winds([B(0, 290, 280)], 1e5, 0.029, 6e6, 7e-5)
    with pytest.raises(bw.BandWindError):
        bw.band_winds([B(10, 290, 280), B(-10, 290, 280)], 1e5, 0.029, 6e6, 7e-5)
    with pytest.raises(bw.BandWindError):
        bw.band_winds([B(0, 0.0, 280), B(60, 270, 260)], 1e5, 0.029, 6e6, 7e-5)
    with pytest.raises(bw.BandWindError):
        bw.band_winds([B(0, 290, 280), B(95, 270, 260)], 1e5, 0.029, 6e6, 7e-5)
    with pytest.raises(bw.BandWindError):
        bw.band_winds(**{**EARTH, "pressure_pa": -1.0})
    with pytest.raises(bw.BandWindError):
        bw.band_winds(**{**EARTH, "mean_molar_mass_kg_mol": 0.0})
    with pytest.raises(bw.BandWindError):
        bw.band_winds(**{**EARTH, "radius_m": 0.0})
