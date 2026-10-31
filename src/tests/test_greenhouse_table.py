"""
The greenhouse table (GEN.850): the fitted stand-in for `radiative_convective.warming` that the surface-air fixed point reads. The table is held to the
model it replaces, to the numbers the model is checked against (Earth, Mars, Venus and Titan are test targets, never inputs of the code), to the
physics a table can break (more absorber or more light is never colder, a thin air has no greenhouse) and to the accuracy its file states. The
comparison with the model itself needs the two infrared datasets in the cache folder and is skipped without them.
"""

import math

import pytest

from planetgen.datasets import fetch, hitran_cia, ir_ktables
from planetgen.physics import constants, elements, greenhouse_table as gt

EARTH = {"N2": 78084.0, "O2": 20946.0, "Ar": 934.0, "CO2": 41.0, "H2O": 1000.0}
EARTH_FLUX = 240.0
MARS = {"CO2": 600.0 * 0.953, "N2": 600.0 * 0.027, "Ar": 600.0 * 0.016}
VENUS = {"CO2": 0.965 * 92.0e5, "N2": 0.035 * 92.0e5}
TITAN = {"N2": 1.47e5 * 0.985, "CH4": 1.47e5 * 0.015}

DATA_PRESENT = not any(fetch.locate(dataset).missing for dataset in (ir_ktables.DATASET, hitran_cia.DATASET))


def _co2(ppm):
    return {**EARTH, "CO2": ppm * 1.0e-6 * 1.0e5}


# --- the file ---------------------------------------------------------------------------------------------------------------------

def test_the_file_is_fitted_to_the_modules_features_and_the_catalogue_data():
    fitted = gt.table()
    assert len(fitted.mean) == len(fitted.scale) == len(fitted.low) == len(fitted.high) == len(gt.FEATURE_NAMES)
    assert fitted.weights[0].shape[0] == len(gt.FEATURE_NAMES) and fitted.weights[-1].shape[1] == 1
    assert all(w.shape[1] == b.shape[0] for w, b in zip(fitted.weights, fitted.biases))


def test_a_table_fitted_to_other_data_is_refused(monkeypatch):
    gt.table.cache_clear()
    monkeypatch.setattr(gt, "data_digest", lambda: "0" * 64)
    try:
        with pytest.raises(elements.DataError, match="fit it again"):
            gt.table()
    finally:
        gt.table.cache_clear()


def test_the_file_states_an_accuracy_the_fit_reached():
    """The validation is on airs the generator made from seeds the fit never saw; the numbers are the file's own. The cool airs, which the habitability and
    the life read, are held to a few percent; the hot dense airs (over 500 K) are the model's weakest and the file says how weak."""
    check = gt.table().summary["validation"]
    cool = check["by_surface_temperature_k"]["0 to 300"]
    assert cool["samples"] >= 500
    assert cool["relative_percentiles"]["50"] < 0.01
    assert cool["relative_percentiles"]["90"] < 0.05
    assert cool["kelvin_percentiles"]["99"] < 25.0
    assert check["relative_percentiles"]["50"] < 0.02


# --- the features -----------------------------------------------------------------------------------------------------------------

def test_the_features_follow_the_names_and_an_absent_gas_sits_at_the_floor():
    vector = gt.features({"N2": 1.0e5}, 9.81, 240.0)
    named = dict(zip(gt.FEATURE_NAMES, vector))
    assert len(vector) == len(gt.FEATURE_NAMES)
    assert named["ln_pressure"] == pytest.approx(math.log(1.0e5))
    assert named["ln_fraction_CO2"] == pytest.approx(math.log(gt.FRACTION_FLOOR))
    assert named["ln_fraction_N2"] == pytest.approx(math.log(1.0 + gt.FRACTION_FLOOR))
    assert named["aloft_share"] == 0.0


def test_gases_with_stand_in_bands_are_counted_by_the_strength_the_gas_table_flags_them_with():
    from planetgen.physics import gases

    strong = next(f for f in ("NH3", "H2S", "SO2") if gases.gas(f).ir == "ir")
    named = dict(zip(gt.FEATURE_NAMES, gt.features({"N2": 0.9e5, strong: 0.1e5}, 9.81, 240.0)))
    assert named["ln_fraction_ir"] == pytest.approx(math.log(0.1 + gt.FRACTION_FLOOR))
    assert named["ln_fraction_weak"] == pytest.approx(math.log(gt.FRACTION_FLOOR))


# --- the physics a table can break ------------------------------------------------------------------------------------------------

def test_a_thin_air_has_no_greenhouse():
    for pressure in (1.0e-6, 1.0e-2, 1.0):
        warming = gt.warming({"N2": pressure}, 9.81, EARTH_FLUX)
        assert warming.factor < 0.02, pressure
        assert warming.surface_temperature_k == pytest.approx(warming.effective_temperature_k, rel=0.01)


def test_more_light_is_never_colder_and_more_carbon_dioxide_is_never_colder():
    fluxes = [60.0, 120.0, 240.0, 480.0, 960.0]
    temperatures = [gt.warming(EARTH, 9.81, flux).surface_temperature_k for flux in fluxes]
    assert temperatures == sorted(temperatures)
    amounts = [gt.warming(_co2(ppm), 9.81, EARTH_FLUX).surface_temperature_k for ppm in (30, 100, 280, 1000, 5000, 30000)]
    assert amounts == sorted(amounts)


def test_the_factor_is_the_repos_variable():
    warming = gt.warming(EARTH, 9.81, EARTH_FLUX)
    assert warming.factor == pytest.approx((warming.surface_temperature_k / warming.effective_temperature_k) ** 4 - 1.0)
    assert warming.effective_temperature_k == pytest.approx((EARTH_FLUX / constants.STEFAN_BOLTZMANN_CONSTANT) ** 0.25)


def test_a_sunlight_taken_aloft_cools_the_surface_of_a_titan_like_air():
    clear = gt.warming(TITAN, 1.35, 2.61).surface_temperature_k
    hazy = gt.warming(TITAN, 1.35, 2.61, 0.6, 3.0e3).surface_temperature_k
    assert hazy < clear - 3.0


def test_an_input_outside_the_fitted_range_is_held_at_its_edge_and_stays_finite():
    for air, gravity, flux in (({"H2": 1.0e12}, 5.0, 1.0e6), ({"CO2": 1.0e-12}, 0.01, 1.0e-6), (VENUS, 80.0, 40.0)):
        warming = gt.warming(air, gravity, flux)
        assert math.isfinite(warming.surface_temperature_k)
        assert warming.effective_temperature_k <= warming.surface_temperature_k <= max(gt.HIGHEST_SURFACE_K, warming.effective_temperature_k)


def test_bad_inputs_are_refused():
    with pytest.raises(ValueError):
        gt.warming(EARTH, 9.81, 0.0)
    with pytest.raises(ValueError):
        gt.warming(EARTH, -1.0, 240.0)
    with pytest.raises(ValueError):
        gt.warming(EARTH, 9.81, 240.0, aloft_share=1.0)


# --- the numbers the model is checked against -------------------------------------------------------------------------------------

@pytest.mark.parametrize("name, air, gravity, flux, haze, model_k, tolerance_k", (
    ("Earth", EARTH, 9.81, EARTH_FLUX, 0.0, 284.0, 4.0),
    ("Mars", MARS, 3.71, 110.0, 0.0, 214.6, 4.0),
    ("Venus", VENUS, 8.87, 160.0, 0.0, 539.5, 40.0),
    ("Titan", TITAN, 1.35, 2.61, 0.0, 106.8, 4.0),
    ("Titan under its haze", TITAN, 1.35, 2.61, 0.6, 93.3, 8.0),
))
def test_the_table_gives_the_temperatures_the_model_is_checked_to(name, air, gravity, flux, haze, model_k, tolerance_k):
    """The model's own results (`docs/design/composition/line-opacity.md`, Checks): Earth with 1 kPa of vapour and the cloud share 284.0 K (275.6 K clear sky; measured 288 K), Mars 210 to 215 K,
    Venus 539.5 K in a clear sky with the sourced far wings of carbon dioxide (GEN.1019; 757 K once the sulfuric acid cloud of GEN.1020 is in the column, measured 735 K), Titan 106.8 K clear and 93.3 K under a haze that takes 60 percent of the sunlight at 30 mbar (measured 94 K)."""
    warming = gt.warming(air, gravity, flux, haze, 3.0e3 if haze else None)
    assert abs(warming.surface_temperature_k - model_k) < tolerance_k


def test_doubling_the_carbon_dioxide_warms_by_about_what_the_model_gives():
    """The model gives +1.72 K with the vapour free (+1.04 K with it held; the table has the vapour free, as the air the fixed point makes does)."""
    base = gt.warming(_co2(280.0), 9.81, EARTH_FLUX).surface_temperature_k
    doubled = gt.warming(_co2(560.0), 9.81, EARTH_FLUX).surface_temperature_k
    assert 1.0 < doubled - base < 2.6


# --- against the model ------------------------------------------------------------------------------------------------------------

@pytest.mark.skipif(not DATA_PRESENT, reason="the infrared datasets are not in this machine's cache folder")
@pytest.mark.parametrize("name, air, gravity, flux, tolerance", (
    ("Earth", EARTH, 9.81, EARTH_FLUX, 0.03),
    ("Earth with ten times the carbon dioxide", _co2(2800.0), 9.81, EARTH_FLUX, 0.03),
    ("Mars", MARS, 3.71, 110.0, 0.03),
    ("a 5 bar carbon dioxide air", {"CO2": 4.8e5, "N2": 2.0e4}, 9.0, 300.0, 0.05),
    ("Venus", VENUS, 8.87, 160.0, 0.06),                       # the hot dense airs are the table's weakest (the file's validation says by how much): 5.5 percent warm at 92 bar
    ("Titan", TITAN, 1.35, 2.61, 0.03),
    ("a cold 2 bar nitrogen-methane air", {"N2": 1.9e5, "CH4": 1.0e4}, 3.0, 8.0, 0.05),
    ("a hydrogen-helium air", {"H2": 8.0e5, "He": 2.0e5}, 12.0, 60.0, 0.05),
))
def test_the_table_agrees_with_the_model_on_airs_the_checks_use(name, air, gravity, flux, tolerance):
    from planetgen.physics import radiative_convective as rc

    model = rc.warming(air, gravity, flux)
    assert model.converged
    table = gt.warming(air, gravity, flux)
    assert table.surface_temperature_k == pytest.approx(model.surface_temperature_k, rel=tolerance), name
