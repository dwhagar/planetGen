# tests/test_sulfur_cloud_set.py

"""GEN.1425: the sulfur-cloud chemistry set V1 to V10 and its six gates (`planetgen.physics.sulfur_cloud_set`)."""

import json
import math

import pytest

from planetgen.physics import atmosphere, sulfur_cloud_set as sc, sulfuric_vapor

GAS = atmosphere.GAS_FIELD
VENUS_PRESSURE_PA = 92.0e5


def hot_volcanic(**changes):
    """A hot carbon-dioxide world with sulfur, a trace of water, carbon monoxide and chlorine under a bright star (a test air, not stored anywhere)."""
    pressure = VENUS_PRESSURE_PA
    gas = {"CO2": 0.965 * pressure, "N2": 0.035 * pressure, "SO2": 150e-6 * pressure, "H2O": 30e-6 * pressure, "CO": 20e-6 * pressure,
           "OCS": 0.3e-6 * pressure, "HCl": 0.45e-6 * pressure, "HF": 0.005e-6 * pressure}
    values = {GAS: gas, "surface_temperature_k": 735.0, "gravity_g": 0.904, "energy_flux_w_m2": 2600.0, "uv_c_share": 0.002, "body_type": "t",
              "hydrosphere": None}
    values.update(changes)
    return values


def cold_clean():
    pressure = 1.0e5
    return {GAS: {"N2": 0.78 * pressure, "O2": 0.21 * pressure, "H2O": 0.01 * pressure}, "surface_temperature_k": 288.0, "gravity_g": 1.0,
            "energy_flux_w_m2": 1361.0, "uv_c_share": 0.002, "body_type": "t", "hydrosphere": "surface ocean"}


def test_the_rule_file_is_the_ten_species_with_tiers_and_numbers_with_bases():
    document = json.load(open(sc.DATA_FILE, encoding="utf-8"))
    assert [row["id"] for row in document["species"]] == [f"V{n}" for n in range(1, 11)]
    assert all(row["default"] in ("on", "off") and row["tier"] and row["tweak"] and row["evidence"] and row["status"] for row in document["species"])
    assert {row["id"] for row in document["species"] if row["default"] == "off"} == {"V8", "V9"}
    assert all(entry["basis"] and entry["tier"] and entry["tweak"] for entry in document["numbers"].values())


def test_a_rule_file_entry_without_a_tweak_is_refused():
    with pytest.raises(sc.SulfurCloudSetError):
        sc._check("x", {"tier": "T3"})


def test_a_hot_volcanic_carbon_dioxide_world_under_a_bright_star_gets_the_acid_cloud_the_gases_and_the_sulfur_species():
    found = sc.chemistry(hot_volcanic())
    for species_id in ("V1", "V2", "V3", "V4", "V5"):
        assert found.has(species_id), (species_id, found.absent.get(species_id))
    assert found.cloud is not None and 345.0 <= found.cloud.base_temperature_k <= 372.0, "the published base is 345 to 367 K"
    assert found.cloud.weight_percent >= 95.0
    assert found.gates[:4] == (1, 2, 3, 4)


def test_a_cold_sulfur_poor_world_gets_none_of_the_set():
    found = sc.chemistry(cold_clean())
    assert found.species == () and found.cloud is None and found.gates == ()
    assert set(found.absent) == set(sc.SPECIES)


def test_an_airless_body_gets_none():
    found = sc.chemistry({GAS: {}, "surface_temperature_k": 300.0, "gravity_g": 1.0})
    assert found.species == () and set(found.absent.values()) == {"no air"}


def test_no_water_no_cloud_and_no_acid_cloud_means_no_acid_species():
    found = sc.chemistry(hot_volcanic(**{GAS: {**hot_volcanic()[GAS], "H2O": 1.0e-8 * VENUS_PRESSURE_PA}}))
    assert not found.has("V1") and not found.has("V6") and not found.has("V7") and not found.has("V10")
    assert found.has("V2"), "the sulfur dioxide is still there when the water floor of the cloud is not met"


def test_no_ultraviolet_no_lower_sulfur_oxides_and_no_sulfur_aerosol():
    found = sc.chemistry(hot_volcanic(uv_c_share=1.0e-6))
    assert found.has("V1") and not found.has("V4") and not found.has("V5")
    assert "gate 3" in found.absent["V4"]
    assert not sc.chemistry(hot_volcanic(energy_flux_w_m2=None)).has("V4"), "no flux, no ultraviolet to read"


def test_no_reductant_no_lower_sulfur_oxides():
    gas = {k: v for k, v in hot_volcanic()[GAS].items() if k not in ("CO", "OCS")}
    assert not sc.chemistry(hot_volcanic(**{GAS: gas})).has("V4")


def test_a_thin_sulfur_dioxide_air_gets_the_lower_oxides_but_not_the_sulfur_aerosol():
    gas = dict(hot_volcanic()[GAS])
    gas["SO2"] = 2.0e-6 * VENUS_PRESSURE_PA
    found = sc.chemistry(hot_volcanic(**{GAS: gas}))
    assert found.has("V4") and not found.has("V5")


def test_the_chlorine_and_iron_species_need_their_supplies():
    assert sc.chemistry(hot_volcanic()).has("V6") and sc.chemistry(hot_volcanic()).has("V7")
    no_chlorine = {k: v for k, v in hot_volcanic()[GAS].items() if k != "HCl"}
    found = sc.chemistry(hot_volcanic(**{GAS: no_chlorine}))
    assert not found.has("V6") and found.has("V7")
    buried = sc.chemistry(hot_volcanic(hydrosphere="surface ocean"))
    assert not buried.has("V6") and not buried.has("V7"), "no exposed iron under a global ocean"
    assert not sc.chemistry(hot_volcanic(body_type="g")).has("V7")


def test_the_single_source_species_are_off_by_default_and_need_their_supplies_when_on():
    values = hot_volcanic()
    assert not sc.chemistry(values).has("V8") and not sc.chemistry(values).has("V9")
    assert "off by default" in sc.chemistry(values).absent["V9"]
    on = sc.chemistry(values, include_off_by_default=True)
    assert on.has("V9") and not on.has("V8"), "the polycyclic hydrocarbons need cosmic dust as well"
    assert sc.chemistry(values, carbon_dust=True, include_off_by_default=True).has("V8")
    assert not sc.chemistry(values, carbon_dust=True).has("V8"), "off by default even with the dust"


def test_the_measured_absorber_comes_with_the_acid_cloud_and_is_replaced_by_the_published_spectrum_of_v8_or_v9():
    assert sc.chemistry(hot_volcanic()).has("V10")
    switched = sc.chemistry(hot_volcanic(), include_off_by_default=True)
    assert not switched.has("V10") and "replaced" in switched.absent["V10"]


def test_the_acid_cloud_base_agrees_with_the_vapor_curve_at_its_own_level():
    found = sc.chemistry(hot_volcanic()).cloud
    gas = hot_volcanic()[GAS]
    share = {k: v / sum(gas.values()) for k, v in gas.items()}
    acid = share["SO2"] * sc.number("acid_vapor_share_of_so2") * found.base_pressure_pa
    water = share["H2O"] * found.base_pressure_pa
    assert sulfuric_vapor.condenses(found.base_temperature_k, acid, water)
    assert not sulfuric_vapor.condenses(found.base_temperature_k + 10.0, acid, min(water, sulfuric_vapor.water_pressure_pa(found.base_temperature_k + 10.0)))


def test_the_absorber_band_is_the_measured_gaussian():
    center, fwhm = sc.number("absorber_center_um"), sc.number("absorber_fwhm_um")
    assert (center, fwhm) == (0.34, 0.14)
    peak = sc.absorber_depth(20.0, 1.0e-3, center)
    assert peak == pytest.approx(sc.number("absorber_depth_per_cloud_depth") * 20.0)
    assert sc.absorber_depth(20.0, 1.0e-3, center + fwhm / 2.0) == pytest.approx(0.5 * peak)
    assert sc.absorber_depth(20.0, 1.0e-3, center - fwhm / 2.0) == pytest.approx(0.5 * peak)
    assert sc.absorber_depth(20.0, 1.0e-3, 0.7) < 0.01 * peak
    assert sc.absorber_depth(20.0, 0.5 * sc.number("absorber_so2_reference_share"), center) == pytest.approx(0.5 * peak)
    assert sc.absorber_depth(0.0, 1.0e-3, center) == 0.0
    assert sc.absorber_spread_um() == (0.03, 0.01)


def test_the_gate_numbers_are_those_the_report_proposes_and_the_floors_are_hand_set():
    assert sc.number("chlorine_min_share") == pytest.approx(1.0e-7) and sc.number("cloud_water_min_share") == pytest.approx(1.0e-6)
    for name in ("acid_vapor_share_of_so2", "cloud_water_min_share", "uv_c_min_w_m2", "reductant_min_share", "chlorine_min_share", "absorber_depth_per_cloud_depth"):
        assert sc.NUMBERS[name]["tweak"] == "L3", name
    assert math.isclose(sc.number("acid_vapor_share_of_so2") * 150.0e-6, 4.5e-6)


def test_the_band_moves_inside_its_stated_uncertainties_with_the_stored_draw_and_not_beyond_two_sigma():
    wavelengths = [0.30 + 0.01 * n for n in range(30)]
    center_sigma, fwhm_sigma = sc.absorber_spread_um()
    peaks = []
    for draw in (0.01, 0.2, 0.5, 0.8, 0.99):
        band = sc.absorber_band(20.0, 1.0e-3, wavelengths, draw)
        peaks.append(wavelengths[band.index(max(band))])
        assert max(band) > 0.0
    assert max(peaks) - min(peaks) <= 4.0 * center_sigma + 0.01 + 1e-9, "the center stays within two sigma either way (to the grid step)"
    assert sc.absorber_band(20.0, 1.0e-3, wavelengths, 0.3) == sc.absorber_band(20.0, 1.0e-3, wavelengths, 0.3), "repeats from the draw"


def test_a_band_given_its_own_center_and_width_is_that_gaussian():
    assert sc.absorber_depth(20.0, 1.0e-3, 0.4, center_um=0.4, fwhm_um=0.1) == pytest.approx(sc.number("absorber_depth_per_cloud_depth") * 20.0)
