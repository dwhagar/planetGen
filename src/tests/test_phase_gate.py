"""GEN.1221: the phase gate of a class family, the window of its liquid at the body's own pressure."""

import pytest

from planetgen.physics import hydrosphere, phase_gate


def test_water_is_liquid_between_freezing_and_the_boiling_point_at_the_pressure():
    low, high = phase_gate.window_k("water", 101325.0)
    assert low == pytest.approx(273.16, abs=0.1) and high == pytest.approx(hydrosphere.water_boiling_k(101325.0))
    assert phase_gate.holds("water", 288.0, 101325.0) and not phase_gate.holds("water", 380.0, 101325.0)
    assert phase_gate.window_k("water", 100.0) is None     # under its triple pressure (611 Pa) water has no liquid at any temperature


def test_carbon_dioxide_needs_518_kilopascals_and_becomes_a_fluid_above_the_critical_pressure():
    assert phase_gate.window_k("co2", 100000.0) is None and not phase_gate.holds("co2", 250.0, 100000.0)
    assert phase_gate.holds("co2", 225.0, 1.0e6)
    assert phase_gate.window_k("co2", 1.0e7)[1] == phase_gate.CO2_FLUID_MAX_K and phase_gate.holds("co2", 350.0, 1.0e7)


def test_the_methane_ethane_sea_lies_between_90_and_the_boiling_point_of_ethane():
    low, high = phase_gate.window_k("methane_ethane", 150000.0)
    assert 90.0 < low < 91.0 and 180.0 < high < 200.0
    assert phase_gate.holds("methane_ethane", 94.0, 150000.0) and not phase_gate.holds("methane_ethane", 94.0, 1.0e-3)


def test_sulfuric_acid_and_ammonia_water_have_their_freezing_points():
    assert phase_gate.window_k("sulfuric_acid", 101325.0)[0] == pytest.approx(283.5)
    assert phase_gate.window_k("ammonia_water", 101325.0)[0] == pytest.approx(176.0)


def test_the_thirds_name_the_positions_of_the_window():
    low, high = phase_gate.window_k("water", 101325.0)
    width = (high - low) / 3.0
    assert phase_gate.third("water", low + 0.5 * width, 101325.0) == "cold"
    assert phase_gate.third("water", low + 1.5 * width, 101325.0) == "ideal"
    assert phase_gate.third("water", low + 2.5 * width, 101325.0) == "hot"
    assert phase_gate.third("water", 200.0, 101325.0) is None


def test_an_unknown_family_and_no_pressure_are_handled():
    with pytest.raises(ValueError):
        phase_gate.window_k("silicon", 1.0e5)
    assert phase_gate.window_k("water", 0.0) is None
