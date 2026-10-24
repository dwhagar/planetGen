"""
The air's density, scale height and greenhouse warming from its mixture (GEN.423, `planetgen.physics.atmosphere_state`).
Earth, Venus and Mars are test targets, never inputs of the code.
"""

import math

import pytest

from planetgen.physics import atmosphere_state as state
from planetgen.physics import constants, greenhouse, infrared_opacity

EARTH = {"N2": 78084.0, "O2": 20946.0, "Ar": 934.0, "CO2": 41.0, "H2O": 1000.0}
EARTH_DRY_ISH = {"N2": 78084.0, "O2": 20946.0, "Ar": 934.0, "CO2": 41.0, "H2O": 400.0}
"""dict: Earth's air with the water vapour of a well-mixed column that holds the measured 25 kg/m2 (400 Pa)."""
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


def test_doubling_carbon_dioxide_warms_earth_by_the_amount_the_forcing_law_gives():
    """The measured forcing of a doubling is 3.7 W/m2 (Byrne and Goldblatt 2014a, 5.32 ln 2); with the water vapour held fixed that
    warms an Earth-like state by about 1.1 K in a full radiative-convective model. The grey equivalent depth (Koll's definition,
    which counts a band as black where a real one still radiates from the cold air above) gives the right sign and order but a third to
    a half of that, and a well-mixed water column masks the 15 um band further (GEN.703 limits)."""
    base = state.greenhouse_warming(EARTH_DRY_ISH, G_EARTH, 240.0)
    doubled = state.greenhouse_warming(dict(EARTH_DRY_ISH, CO2=2.0 * EARTH_DRY_ISH["CO2"]), G_EARTH, 240.0)
    rise = doubled.surface_temperature_k - base.surface_temperature_k
    assert 0.2 < rise < 2.5, rise
    # a band that saturates takes less than twice the optical depth from twice the gas
    assert base.optical_depth < doubled.optical_depth < 1.5 * base.optical_depth


def test_less_gravity_means_more_column_per_pascal_and_more_opacity():
    low = state.greenhouse_warming(EARTH_DRY_ISH, G_EARTH / 2.0, 240.0)
    high = state.greenhouse_warming(EARTH_DRY_ISH, G_EARTH, 240.0)
    assert low.optical_depth > high.optical_depth


def test_an_air_with_no_absorber_has_no_warming():
    warm = state.greenhouse_warming({"Ar": 9.0e4, "He": 1.0e4}, G_EARTH, 240.0)
    assert warm.optical_depth == 0.0
    assert warm.factor == pytest.approx(0.0, abs=1e-12)
    assert warm.surface_temperature_k == pytest.approx(warm.effective_temperature_k, rel=1e-12)


def test_a_dry_nitrogen_air_has_only_the_little_its_collisions_give():
    warm = state.greenhouse_warming({"N2": 7.0e4, "O2": 2.0e4, "Ar": 1.0e3}, G_EARTH, 240.0)
    assert 0.0 < warm.optical_depth < 0.1
    assert warm.factor < 0.05
    assert set(warm.by_gas) == set()


def test_the_warming_rises_with_the_greenhouse_gas():
    previous = -1.0
    for co2 in (1.0, 1.0e2, 1.0e3, 1.0e4, 1.0e5, 1.0e6, 1.0e7):
        warm = state.greenhouse_warming({"N2": 5.0e4, "CO2": co2}, G_EARTH, 240.0)
        assert warm.factor > previous
        previous = warm.factor


def test_changing_the_mix_moves_the_warming_as_the_formulas_say():
    dry = {"N2": 8.0e4, "O2": 2.0e4}
    wet = dict(dry, H2O=2000.0)
    warm_dry = state.greenhouse_warming(dry, G_EARTH, 240.0)
    warm_wet = state.greenhouse_warming(wet, G_EARTH, 240.0)
    assert warm_wet.optical_depth > warm_dry.optical_depth
    expected = greenhouse.surface_temperature_k(warm_wet.optical_depth, 240.0, warm_wet.gamma, warm_wet.lapse_scaling,
                                                infrared_opacity.PRESSURE_EXPONENT)[0]
    assert warm_wet.surface_temperature_k == pytest.approx(expected, rel=1e-12)
    assert warm_wet.factor == pytest.approx((expected / warm_wet.effective_temperature_k) ** 4 - 1.0, rel=1e-12)
    assert set(warm_wet.by_gas) == {"H2O"} and warm_wet.by_part["water_continuum"] > 0.0


def test_earth_venus_and_mars_land_where_the_world_does():
    """Not inputs: the band, wing and continuum numbers come from the research (part 03), not from these bodies. The measured
    targets are G 0.625 for Earth, 99.6 for Venus and 0.10 to 0.21 for Mars."""
    earth = state.greenhouse_warming(EARTH_DRY_ISH, G_EARTH, 240.0)
    venus = state.greenhouse_warming(VENUS, G_VENUS, 160.0)
    mars = state.greenhouse_warming(MARS, G_MARS, 110.0)
    assert earth.factor == pytest.approx(0.625, rel=0.3)
    assert venus.factor == pytest.approx(99.6, rel=0.4)
    assert 0.03 < mars.factor < 0.25
    assert mars.factor < earth.factor < venus.factor


def test_heat_capacity_ratio_of_simple_gases():
    assert state.heat_capacity_ratio({"N2": 1.0e5}, 300.0) == pytest.approx(1.40, abs=0.01)
    assert state.heat_capacity_ratio({"He": 1.0e5}, 300.0) == pytest.approx(5.0 / 3.0, abs=0.005)
    assert state.heat_capacity_ratio({"CO2": 1.0e5}, 300.0) == pytest.approx(1.29, abs=0.01)
    # below the fit's range the cold end stands in instead of raising
    assert math.isfinite(state.heat_capacity_ratio({"CO2": 1.0e5}, 50.0))


@pytest.mark.parametrize("mix,flux,gravity", [
    ({"H2": 8.0e4, "He": 2.0e4, "CH4": 20.0}, 1.0, 25.0),
    ({"H2O": 2.0e7, "CO2": 1.0e6}, 400.0, 9.81),
    ({"CO2": 1.0e-3}, 50.0, 1.0),
    ({"N2": 1.5e5, "CH4": 6.0e3}, 2.6, 1.35),
    ({"Na": 1.0, "O2": 1.0e-2}, 3000.0, 2.0),
    ({"H2O": 2.2e7}, 3000.0, 9.81),
    ({"CO2": 9.0e7, "N2": 3.0e6}, 160.0, 8.87),
    ({"Ar": 1.0e5}, 240.0, 9.81),
    ({"NH3": 5.0e3, "H2": 9.0e4}, 5.0, 10.0),
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
