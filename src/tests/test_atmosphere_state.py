"""
The air's density, scale height and greenhouse warming from its mixture (GEN.423, `planetgen.physics.atmosphere_state`).
Earth, Venus and Mars are test targets, never inputs of the code.
"""

import math

import pytest

from planetgen.physics import atmosphere_state as state
from planetgen.physics import constants, greenhouse

EARTH = {"N2": 78084.0, "O2": 20946.0, "Ar": 934.0, "CO2": 41.0, "H2O": 1000.0}
VENUS = {"CO2": 0.965 * 92.0e5, "N2": 0.035 * 92.0e5}
MARS = {"CO2": 600.0 * 0.953, "N2": 600.0 * 0.027, "Ar": 600.0 * 0.016}
G_EARTH, G_VENUS, G_MARS = 9.81, 8.87, 3.71


def _swap(mix, remove, add, pressure):
    out = {k: v for k, v in mix.items() if k != remove}
    out[add] = pressure
    return out


def test_earth_air_state():
    air = state.air_state(EARTH, 288.0, G_EARTH)
    assert air.mean_molar_mass_kg_mol == pytest.approx(0.0289, rel=0.01)
    assert air.density_kg_m3 == pytest.approx(1.22, rel=0.02)
    assert air.scale_height_m == pytest.approx(8.4e3, rel=0.03)
    assert air.column_mass_kg_m2 == pytest.approx(1.03e4, rel=0.02)


def test_the_barometric_formula_closes():
    for mix, temperature, gravity in ((EARTH, 288.0, G_EARTH), (VENUS, 737.0, G_VENUS), (MARS, 215.0, G_MARS)):
        air = state.air_state(mix, temperature, gravity)
        assert air.density_kg_m3 * gravity * air.scale_height_m == pytest.approx(air.total_pressure_pa, rel=1e-12)


def test_a_heavier_mix_has_a_lower_scale_height_and_a_higher_density():
    light = state.air_state(_swap(EARTH, "N2", "He", 78084.0), 288.0, G_EARTH)
    heavy = state.air_state(_swap(EARTH, "N2", "CO2", 78084.0), 288.0, G_EARTH)
    base = state.air_state(EARTH, 288.0, G_EARTH)
    assert light.scale_height_m > base.scale_height_m > heavy.scale_height_m
    assert light.density_kg_m3 < base.density_kg_m3 < heavy.density_kg_m3
    # the scale height follows 1 / M exactly at one temperature and gravity
    assert (heavy.scale_height_m * heavy.mean_molar_mass_kg_mol
            == pytest.approx(base.scale_height_m * base.mean_molar_mass_kg_mol, rel=1e-12))
    assert base.scale_height_m == pytest.approx(constants.R * 288.0 / (base.mean_molar_mass_kg_mol * G_EARTH), rel=1e-12)


def test_gas_columns_add_up_to_the_column_mass():
    columns = state.gas_column_mass_kg_m2(EARTH, G_EARTH)
    assert sum(columns.values()) == pytest.approx(state.air_state(EARTH, 288.0, G_EARTH).column_mass_kg_m2, rel=1e-12)
    # nitrogen weighs more per molecule than the water above it; the heavy trace gas outweighs its share of molecules
    assert columns["CO2"] / columns["H2O"] == pytest.approx(41.0 * 44.0095 / (1000.0 * 18.015), rel=0.01)


def test_doubling_a_gas_doubles_its_optical_depth():
    base = state.optical_depth_by_gas({"N2": 1.0e5, "CO2": 1.0e4}, G_EARTH)
    # both the gas and the total change, so compare the dilute limit: a trace gas on a big inert bed
    trace = state.optical_depth_by_gas({"N2": 1.0e5, "CO2": 10.0}, G_EARTH)
    doubled = state.optical_depth_by_gas({"N2": 1.0e5, "CO2": 20.0}, G_EARTH)
    assert doubled["CO2"] == pytest.approx(2.0 * trace["CO2"], rel=2e-3)
    assert base["CO2"] > trace["CO2"]
    assert set(base) == {"CO2"}


def test_optical_depth_follows_the_formula():
    mix = {"N2": 9.0e4, "CO2": 1.0e4}
    columns = state.gas_column_mass_kg_m2(mix, G_EARTH)
    assert state.infrared_optical_depth(mix, G_EARTH) == pytest.approx(state.mass_absorption_m2_kg("CO2") * columns["CO2"], rel=1e-12)


def test_less_gravity_means_more_column_per_pascal():
    assert (state.infrared_optical_depth(EARTH, G_EARTH / 2.0)
            == pytest.approx(2.0 * state.infrared_optical_depth(EARTH, G_EARTH), rel=1e-12))


def test_an_air_with_no_line_absorber_has_no_warming():
    warm = state.greenhouse_warming({"N2": 7.0e4, "O2": 2.0e4, "Ar": 1.0e3}, G_EARTH, 240.0)
    assert warm.optical_depth == 0.0
    assert warm.factor == pytest.approx(0.0, abs=1e-12)
    assert warm.surface_temperature_k == pytest.approx(warm.effective_temperature_k, rel=1e-12)


def test_the_warming_rises_with_the_greenhouse_gas():
    previous = -1.0
    for co2 in (1.0, 1.0e2, 1.0e3, 1.0e4, 1.0e5, 1.0e6, 1.0e7):
        warm = state.greenhouse_warming({"N2": 5.0e4, "CO2": co2}, G_EARTH, 240.0)
        assert warm.factor > previous
        previous = warm.factor


def test_changing_the_mix_moves_the_warming_as_the_formulas_say():
    dry = {"N2": 8.0e4, "O2": 2.0e4}
    wet = dict(dry, H2O=2000.0)
    tau_dry = state.greenhouse_warming(dry, G_EARTH, 240.0)
    tau_wet = state.greenhouse_warming(wet, G_EARTH, 240.0)
    assert tau_dry.optical_depth == 0.0 and tau_wet.optical_depth > 0.0
    expected = greenhouse.surface_temperature_k(tau_wet.optical_depth, 240.0, tau_wet.gamma, state.LAPSE_SCALING,
                                                state.PRESSURE_EXPONENT)[0]
    assert tau_wet.surface_temperature_k == pytest.approx(expected, rel=1e-12)
    assert tau_wet.factor == pytest.approx((expected / tau_wet.effective_temperature_k) ** 4 - 1.0, rel=1e-12)


def test_earth_venus_and_mars_land_where_the_world_does():
    # not inputs: the stand-in opacities come from Koll 2019's statements, not from these bodies
    earth = state.greenhouse_warming(EARTH, G_EARTH, 240.0)
    venus = state.greenhouse_warming(VENUS, G_VENUS, 160.0)
    mars = state.greenhouse_warming(MARS, G_MARS, 110.0)
    assert 0.2 < earth.factor < 1.0           # the measured 0.625; water vapour held at a fixed 1 kPa
    assert 60.0 < venus.factor < 150.0        # the measured 99.6
    assert mars.factor < 0.1                  # a thin dry world: almost none
    assert mars.factor < earth.factor < venus.factor


def test_heat_capacity_ratio_of_simple_gases():
    assert state.heat_capacity_ratio({"N2": 1.0e5}, 300.0) == pytest.approx(1.40, abs=0.01)
    assert state.heat_capacity_ratio({"He": 1.0e5}, 300.0) == pytest.approx(5.0 / 3.0, abs=0.005)
    assert state.heat_capacity_ratio({"CO2": 1.0e5}, 300.0) == pytest.approx(1.29, abs=0.01)
    # below the fit's range the cold end stands in instead of raising
    assert math.isfinite(state.heat_capacity_ratio({"CO2": 1.0e5}, 50.0))


def test_flags_decide_a_gas_with_no_entry():
    assert state.mass_absorption_m2_kg("CO2") == state.INFRARED_MASS_ABSORPTION_M2_KG["CO2"]
    assert state.mass_absorption_m2_kg("CH4") == state.GENERIC_BAND_M2_KG
    assert state.mass_absorption_m2_kg("N2") == 0.0      # collision-induced only
    assert state.mass_absorption_m2_kg("He") == 0.0


@pytest.mark.parametrize("mix,flux,gravity", [
    ({"H2": 8.0e4, "He": 2.0e4, "CH4": 20.0}, 1.0, 25.0),
    ({"H2O": 2.0e7, "CO2": 1.0e6}, 400.0, 9.81),
    ({"CO2": 1.0e-3}, 50.0, 1.0),
    ({"N2": 1.5e5, "CH4": 6.0e3}, 2.6, 1.35),
    ({"Na": 1.0, "O2": 1.0e-2}, 3000.0, 2.0),
])
def test_no_nan_for_any_family(mix, flux, gravity):
    air = state.air_state(mix, 300.0, gravity)
    warm = state.greenhouse_warming(mix, gravity, flux)
    for value in (air.density_kg_m3, air.scale_height_m, warm.optical_depth, warm.surface_temperature_k, warm.factor, warm.gamma):
        assert math.isfinite(value)
    assert warm.factor >= -1e-12


def test_bad_input_is_refused():
    with pytest.raises(ValueError):
        state.air_state({}, 288.0, G_EARTH)
    with pytest.raises(ValueError):
        state.air_state(EARTH, 0.0, G_EARTH)
    with pytest.raises(ValueError):
        state.air_state(EARTH, 288.0, -1.0)
    with pytest.raises(ValueError):
        state.greenhouse_warming(EARTH, G_EARTH, 0.0)
