# tests/test_sulfur_cloud_set.py

"""GEN.1425: the sulfur-cloud chemistry set V1 to V10 and its six gates (`planetgen.physics.sulfur_cloud_set`)."""

import json
import math

import pytest

from types import SimpleNamespace

from planetgen.physics import atmosphere, elements, stellar_spectrum, sulfur_cloud_set as sc, sulfuric_vapor

GAS = atmosphere.GAS_FIELD
VENUS_PRESSURE_PA = 92.0e5


def chem(values, **kwargs):
    """`sulfur_cloud_set.chemistry` of `values`, whose `uv_c_w_m2` is the star's ultraviolet-C flux at the body (GEN.1872: `chemistry` takes it from the caller, `fill` works it out of the star)."""
    values = dict(values)
    return sc.chemistry(values, uv_c_w_m2=values.pop("uv_c_w_m2", None), **kwargs)


def hot_volcanic(**changes):
    """A hot carbon-dioxide world with sulfur, a trace of water, carbon monoxide and chlorine under a bright star (a test air, not stored anywhere)."""
    pressure = VENUS_PRESSURE_PA
    gas = {"CO2": 0.965 * pressure, "N2": 0.035 * pressure, "SO2": 150e-6 * pressure, "H2O": 30e-6 * pressure, "CO": 20e-6 * pressure,
           "OCS": 0.3e-6 * pressure, "HCl": 0.45e-6 * pressure, "HF": 0.005e-6 * pressure}
    values = {GAS: gas, "surface_temperature_k": 735.0, "gravity_g": 0.904, "energy_flux_w_m2": 2600.0, "uv_c_w_m2": 5.0, "body_type": "t",
              "hydrosphere": None}
    values.update(changes)
    return values


def cold_clean():
    pressure = 1.0e5
    return {GAS: {"N2": 0.78 * pressure, "O2": 0.21 * pressure, "H2O": 0.01 * pressure}, "surface_temperature_k": 288.0, "gravity_g": 1.0,
            "energy_flux_w_m2": 1361.0, "uv_c_w_m2": 5.0, "body_type": "t", "hydrosphere": "surface ocean"}


def test_the_rule_file_is_the_ten_species_with_tiers_and_numbers_with_bases():
    document = json.load(open(sc.DATA_FILE, encoding="utf-8"))
    assert [row["id"] for row in document["species"]] == [f"V{n}" for n in range(1, 11)]
    assert all(row["default"] in ("on", "off") and row["tier"] and row["tweak"] and row["evidence"] and row["status"] for row in document["species"])
    assert {row["id"] for row in document["species"] if row["default"] == "off"} == {"V8", "V9"}
    assert all(entry["basis"] and entry["tier"] and entry["tweak"] for entry in document["numbers"].values())


def test_a_rule_file_entry_without_a_tweak_is_refused():
    with pytest.raises(sc.SulfurCloudSetError):
        elements.check_rule_entry("x", {"tier": "T3"}, sc.SulfurCloudSetError, basis=False)


def test_a_hot_volcanic_carbon_dioxide_world_under_a_bright_star_gets_the_acid_cloud_the_gases_and_the_sulfur_species():
    found = chem(hot_volcanic())
    for species_id in ("V1", "V2", "V3", "V4", "V5"):
        assert found.has(species_id), (species_id, found.absent.get(species_id))
    assert found.cloud is not None and 345.0 <= found.cloud.base_temperature_k <= 372.0, "the published base is 345 to 367 K"
    assert found.cloud.weight_percent >= 95.0
    assert found.gates[:4] == (1, 2, 3, 4)


def test_a_cold_sulfur_poor_world_gets_none_of_the_set():
    found = chem(cold_clean())
    assert found.species == () and found.cloud is None and found.gates == ()
    assert set(found.absent) == set(sc.SPECIES)


def test_an_airless_body_gets_none():
    found = chem({GAS: {}, "surface_temperature_k": 300.0, "gravity_g": 1.0})
    assert found.species == () and set(found.absent.values()) == {"no air"}


def test_no_water_no_cloud_and_no_acid_cloud_means_no_acid_species():
    found = chem(hot_volcanic(**{GAS: {**hot_volcanic()[GAS], "H2O": 1.0e-8 * VENUS_PRESSURE_PA}}))
    assert not found.has("V1") and not found.has("V6") and not found.has("V7") and not found.has("V10")
    assert found.has("V2"), "the sulfur dioxide is still there when the water floor of the cloud is not met"


def test_no_ultraviolet_no_lower_sulfur_oxides_and_no_sulfur_aerosol():
    found = chem(hot_volcanic(uv_c_w_m2=0.01))
    assert found.has("V1") and not found.has("V4") and not found.has("V5")
    assert "gate 3" in found.absent["V4"]
    assert not chem(hot_volcanic(uv_c_w_m2=None)).has("V4"), "no flux, no ultraviolet to read"


def test_no_reductant_no_lower_sulfur_oxides():
    gas = {k: v for k, v in hot_volcanic()[GAS].items() if k not in ("CO", "OCS")}
    assert not chem(hot_volcanic(**{GAS: gas})).has("V4")


def test_a_thin_sulfur_dioxide_air_gets_the_lower_oxides_but_not_the_sulfur_aerosol():
    gas = dict(hot_volcanic()[GAS])
    gas["SO2"] = 2.0e-6 * VENUS_PRESSURE_PA
    found = chem(hot_volcanic(**{GAS: gas}))
    assert found.has("V4") and not found.has("V5")


def test_the_chlorine_and_iron_species_need_their_supplies():
    assert chem(hot_volcanic()).has("V6") and chem(hot_volcanic()).has("V7")
    no_chlorine = {k: v for k, v in hot_volcanic()[GAS].items() if k != "HCl"}
    found = chem(hot_volcanic(**{GAS: no_chlorine}))
    assert not found.has("V6") and found.has("V7")
    buried = chem(hot_volcanic(hydrosphere="surface ocean"))
    assert not buried.has("V6") and not buried.has("V7"), "no exposed iron under a global ocean"
    assert not chem(hot_volcanic(body_type="g")).has("V7")


def test_the_single_source_species_are_off_by_default_and_need_their_supplies_when_on():
    values = hot_volcanic()
    assert not chem(values).has("V8") and not chem(values).has("V9")
    assert "off by default" in chem(values).absent["V9"]
    on = chem(values, include_off_by_default=True)
    assert on.has("V9") and not on.has("V8"), "the polycyclic hydrocarbons need cosmic dust as well"
    assert chem(values, carbon_dust=True, include_off_by_default=True).has("V8")
    assert not chem(values, carbon_dust=True).has("V8"), "off by default even with the dust"


def test_the_measured_absorber_comes_with_the_acid_cloud_and_is_replaced_by_the_published_spectrum_of_v8_or_v9():
    assert chem(hot_volcanic()).has("V10")
    switched = chem(hot_volcanic(), include_off_by_default=True)
    assert not switched.has("V10") and "replaced" in switched.absent["V10"]


def test_the_acid_cloud_base_agrees_with_the_vapor_curve_at_its_own_level():
    found = chem(hot_volcanic()).cloud
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


def test_the_stored_text_round_trips_and_reads_as_words():
    """GEN.1425: a body stores the ids of its species in order; a page reads them as the names and the status of each."""
    found = chem(hot_volcanic())
    text = sc.encode(found)
    assert text == ",".join(found.species) and sc.species_of(text) == found.species
    words = sc.sentence(text)
    assert "sulfuric acid water droplets: observed in a planet (Venus)" in words and words.count(";") == len(found.species) - 1
    assert sc.encode(chem(cold_clean())) is None and sc.species_of(None) == () and sc.sentence(None) == ""


def test_a_given_cloud_decides_gate_two_so_it_agrees_with_the_aerosols():
    """GEN.1425: the cloud the body's aerosols were drawn with is gate 2; a body whose aerosols drew none gets no V1, V6, V7 or V10 though its values would open it."""
    values = hot_volcanic()
    assert chem(values).has("V1")
    with_none = chem(values, cloud=None)
    assert not with_none.has("V1") and not with_none.has("V10") and with_none.has("V2")
    cloud = sc.acid_cloud(values)
    assert chem(values, cloud=cloud).species == chem(values).species
    assert chem(cold_clean(), cloud=cloud).species == ()      # no sulfur dioxide: gate 1 stays shut whatever cloud is given


def test_the_gas_line_names_only_the_gases_above_their_floors():
    """UX.309: V3 opens on any one of its four gases, so a page names those the air holds, not the four of the set."""
    only_hcl = sc.sentence("V3", {"N2": 24.8e3, "HCl": 8.6e-3})
    assert only_hcl.startswith("hydrogen chloride gas:") and "fluoride" not in only_hcl and "monoxide" not in only_hcl
    two = sc.sentence("V3", {"N2": 24.8e3, "HCl": 8.6e-3, "CO": 1.0})
    assert two.startswith("hydrogen chloride and carbon monoxide gases:")
    assert sc.sentence("V3") == sc.sentence("V3", None) and "hydrogen fluoride" in sc.sentence("V3")


def _star(temperature, uv_c_share):
    """A star with only what `uv_c_flux_w_m2` reads (no mass or radius, so `activity.uv_shares_of` has no library spectrum and the share is the star's `uv_c_share`)."""
    return SimpleNamespace(temperature=temperature, uv_c_share=uv_c_share)


def test_the_ultraviolet_c_flux_is_the_blackbody_band_of_the_star_times_its_share_of_that_band():
    """GEN.1872: `uv_c_share` is the star's UV-C power as a multiple of its blackbody's, so the flux is the total flux times the blackbody fraction of 100 to 280 nm at the star's temperature
    times the share, as `body_radiation` gives the spectrum, and not the total flux times the share."""
    flux, temperature, share = 1361.0, 5772.0, 0.3
    band = stellar_spectrum.blackbody_fraction(100.0, 280.0, temperature)
    assert sc.uv_c_flux_w_m2(_star(temperature, share), flux) == pytest.approx(flux * band * share)
    assert sc.uv_c_flux_w_m2(_star(temperature, share), flux) < 0.1 * flux * share, "the total flux times the share (the old reading) is some tens of times too much for a Sun-like star"


def _m_dwarfs(count=20):
    """`count` M2V stars (each drawn on a stream of its own) with the flux they give at the warm edge of their habitable zone: `(star, flux_w_m2)`."""
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.star import Star
    from planetgen.util import draw

    found = []
    for seed in range(count):
        config = SystemConfig()
        config.STAR_TYPE, config.BINARY_SYSTEM = "M2V", False
        with draw.bound(f"sulfur-uv-c:{seed}"):
            star = Star(config)
        found.append((star, stellar_spectrum.bolometric_flux_w_m2(star.luminosity, star.habitable_zone[0])))
    return found


def test_a_quiet_m_dwarf_does_not_pass_gate_3_and_an_active_one_does_though_the_old_reading_passed_both():
    """GEN.1872: the old reading (total flux times `uv_c_share`) put every one of these at hundreds of W/m2; the star's UV-C band at its warm edge is some W/m2 for an active one and a tenth of that or less for a quiet one."""
    stars = _m_dwarfs()
    floor = sc.number("uv_c_min_w_m2")
    assert all(flux * star.uv_c_share >= floor for star, flux in stars if star.uv_c_share >= 0.5), "the old reading passes every star of a share of one half or more"
    quiet = [(star, flux) for star, flux in stars if sc.uv_c_flux_w_m2(star, flux) < floor]
    active = [(star, flux) for star, flux in stars if sc.uv_c_flux_w_m2(star, flux) >= floor]
    assert quiet and active, "both kinds occur among twenty M2V stars"
    assert any(flux * star.uv_c_share >= floor for star, flux in quiet), "a quiet star the old reading let through"
    star, flux = quiet[0]
    assert not chem(hot_volcanic(uv_c_w_m2=sc.uv_c_flux_w_m2(star, flux))).has("V4")
    star, flux = active[0]
    assert chem(hot_volcanic(uv_c_w_m2=sc.uv_c_flux_w_m2(star, flux))).has("V4")


def test_the_ultraviolet_c_flux_is_unknown_when_the_flux_the_temperature_or_the_share_is():
    assert sc.uv_c_flux_w_m2(_star(5772.0, 0.3), None) is None and sc.uv_c_flux_w_m2(_star(5772.0, 0.3), 0.0) is None
    assert sc.uv_c_flux_w_m2(_star(None, 0.3), 1361.0) is None and sc.uv_c_flux_w_m2(_star(5772.0, None), 1361.0) is None
    assert sc.uv_c_flux_w_m2(None, 1361.0) is None
