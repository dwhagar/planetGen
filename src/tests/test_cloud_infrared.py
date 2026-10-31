# tests/test_cloud_infrared.py

"""GEN.1020: the infrared optics of the sulfuric acid cloud (`planetgen.physics.cloud_infrared`) and the minor-gas far wing of `co2_far_wing`."""

import numpy as np
import pytest

from planetgen.physics import cloud_infrared as ci
from planetgen.physics import co2_far_wing as fw
from planetgen.physics import line_opacity, sulfuric_cloud

CENTRES = line_opacity.CENTRES_CM1


def band(wavenumber):
    return int(np.argmin(np.abs(CENTRES - wavenumber)))


def test_the_table_has_a_row_for_every_radius_and_a_value_for_every_band():
    assert ci.DOCUMENT["bands"] == line_opacity.BANDS
    for percent in ci.PERCENTS:
        assert len(ci.DOCUMENT["optics"][str(percent)]) == len(ci.RADII_UM)
        for row in ci.DOCUMENT["optics"][str(percent)]:
            assert all(len(row[key]) == line_opacity.BANDS for key in ("ext_over_visible", "omega", "g"))


def test_the_extinction_rises_with_the_radius_in_the_thermal_infrared():
    ten_micron = band(1000.0)
    ratio = [ci.optics(r)[0][ten_micron] for r in (0.3, 0.7, 1.4, 3.85)]
    assert np.all(np.diff(ratio) > 0.0)


def test_the_cloud_is_a_clearer_window_near_four_micron_and_a_scatterer_there():
    """Sulfuric acid barely absorbs between 3.5 and 5 micron: the albedo of the droplets is near 1 there and low near 10 micron."""
    ratio, omega, g = ci.optics(1.0)
    assert omega[band(2400.0)] < omega[band(4400.0)] > 0.95 and omega[band(900.0)] < 0.5 and 0.0 < g.min() and g.max() < 1.0


def test_the_visible_efficiency_is_a_mode_two_value():
    row = ci.DOCUMENT["optics"]["81"][3]
    assert 2.0 < row["q_visible"] < 3.0 and ci.DOCUMENT["visible_real_index"] == pytest.approx(1.4345, abs=0.01)


def test_the_effective_depth_is_the_visible_depth_scaled_and_linear():
    ratio, omega, g = ci.optics(1.0)
    assert ci.effective_depth(10.0, 1.0) == pytest.approx(10.0 * ratio * (1.0 - omega * g))
    assert ci.effective_depth(20.0, 1.0) == pytest.approx(2.0 * ci.effective_depth(10.0, 1.0))


def test_the_acid_strength_is_read_between_the_two_tables():
    low, high, middle = (ci.optics(1.0, percent) for percent in (76.0, 81.0, 78.5))
    for a, b, c in zip(low, high, middle):
        assert c == pytest.approx(0.5 * (a + b))
    assert ci.optics(1.0, 70.0)[0] == pytest.approx(low[0]) and ci.optics(1.0, 90.0)[0] == pytest.approx(high[0])


def test_the_layers_hold_the_whole_depth_between_the_base_and_the_top():
    pressure = np.exp(np.linspace(np.log(1.0e3), np.log(9.2e6), 40))
    mass = np.gradient(pressure) / 8.87
    temperature = np.where(pressure < 3.0e5, 250.0, 700.0)
    depth = ci.layer_depth(pressure, mass, temperature, 25.0, 1.0)
    assert depth.shape == (40, line_opacity.BANDS)
    assert depth.sum(axis=0) == pytest.approx(ci.effective_depth(25.0, 1.0), rel=1e-9)
    inside = depth.sum(axis=1) > 0.0
    assert pressure[inside].max() < 3.0e5 and pressure[inside].min() >= 0.9 * 3.0e5 * ci.TOP_OVER_BASE_PRESSURE
    assert not depth[~inside].any()


def test_no_cloud_and_no_cool_air_give_no_depth():
    pressure = np.exp(np.linspace(np.log(1.0e3), np.log(9.2e6), 40))
    mass = np.gradient(pressure) / 8.87
    assert not ci.layer_depth(pressure, mass, np.full(40, 250.0 + 400.0 * (pressure > 3.0e5)), 0.0, 1.0).any()
    assert not ci.layer_depth(pressure, mass, np.full(40, 700.0), 25.0, 1.0).any()


def test_an_all_cool_air_puts_the_base_at_the_ground():
    pressure = np.exp(np.linspace(np.log(1.0e3), np.log(1.0e5), 40))
    mass = np.gradient(pressure) / 9.81
    depth = ci.layer_depth(pressure, mass, np.full(40, 300.0), 20.0, 1.0)
    assert depth[-1].any() and depth.sum(axis=0) == pytest.approx(ci.effective_depth(20.0, 1.0), rel=1e-9)


def test_the_top_ratio_is_the_published_pressure_ratio():
    assert ci.TOP_OVER_BASE_PRESSURE == pytest.approx(0.0267 / 1.391)
    assert sulfuric_cloud.number("base_temperature_k") == 367.0


EDGES = np.arange(0.0, 1100.0, 100.0)


def one_band(index, strength=1.0e-18):
    strengths = np.zeros((1, len(EDGES) - 1))
    strengths[0, index] = strength
    return strengths


def test_a_minor_gas_wing_is_zero_without_carbon_dioxide_and_grows_with_it():
    none = fw.minor_gas_cross_section_cm2("H2O", one_band(2), EDGES, [400.0], [5.0e6], [0.0])
    some = fw.minor_gas_cross_section_cm2("H2O", one_band(2), EDGES, [400.0], [5.0e6], [0.5])
    full = fw.minor_gas_cross_section_cm2("H2O", one_band(2), EDGES, [400.0], [5.0e6], [1.0])
    assert not none.any() and full.sum() == pytest.approx(2.0 * some.sum(), rel=1e-9)


def test_the_pure_lorentz_tail_is_the_upper_bracket_of_the_minor_wing(monkeypatch):
    with_chi = fw.minor_gas_cross_section_cm2("SO2", one_band(2), EDGES, [400.0], [5.0e6], [1.0])
    monkeypatch.setitem(fw.NUMBERS["minor_uses_chi"], "value", 0)
    without = fw.minor_gas_cross_section_cm2("SO2", one_band(2), EDGES, [400.0], [5.0e6], [1.0])
    assert np.all(without >= with_chi - 1e-30) and without.sum() > 1.5 * with_chi.sum()


def test_a_wider_line_gives_a_stronger_wing():
    h2o = fw.minor_gas_cross_section_cm2("H2O", one_band(2), EDGES, [300.0], [5.0e6], [1.0]).sum()
    co = fw.minor_gas_cross_section_cm2("CO", one_band(2), EDGES, [300.0], [5.0e6], [1.0]).sum()
    assert h2o > co
