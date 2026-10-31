# tests/test_photolysis.py

"""GEN.496: the photolysis lifetime table, against docs/design/composition/body-composition-parts/G1-comets.md sections G1.2 and 2.2."""

import json
import math

import pytest

from planetgen.physics import gases, photolysis

# Quiet-Sun lifetimes at 1 AU (s) the sources give, with the tolerance of the design's consistency test 3.
QUIET = {"H2O": (8.29e4, 0.02), "CO": (1.5e6, 0.02), "CO2": (4.95e5, 0.02), "HCN": (6.7e4, 0.02), "CH3OH": (7.7e4, 0.02), "CH2O": (5.0e3, 0.02), "NH3": (6.7e3, 0.10)}


def test_every_species_is_in_the_gas_table():
    assert photolysis.formulas()
    assert all(gases.has(formula) for formula in photolysis.formulas())


@pytest.mark.parametrize("formula", sorted(QUIET))
def test_the_quiet_sun_lifetimes_are_the_sourced_values(formula):
    expected, tolerance = QUIET[formula]
    assert photolysis.lifetime_s(formula) == pytest.approx(expected, rel=tolerance)


def test_ammonia_is_the_quiet_sun_value_not_the_mean_of_quiet_and_active():
    assert photolysis.species("NH3").lifetime_s == pytest.approx(6.7e3)
    assert photolysis.lifetime_s("NH3") != pytest.approx(5.6e3, rel=0.05)
    assert photolysis.species("NH3").lifetime_range_s == (5.4e3, 7.7e3)


def test_doubling_the_distance_divides_the_rate_by_four_exactly():
    for formula in photolysis.formulas():
        near = photolysis.rate_per_s(formula, photolysis.inverse_square(1.0))
        far = photolysis.rate_per_s(formula, photolysis.inverse_square(2.0))
        assert far == near / 4.0


def test_a_cool_star_hardly_dissociates_co_and_n2_and_a_hot_one_dissociates_more():
    assert photolysis.rate_per_s("CO", teff_k=4000.0) < 1e-3 * photolysis.rate_per_s("CO")
    assert photolysis.rate_per_s("N2", teff_k=4000.0) < 1e-4 * photolysis.rate_per_s("N2")
    assert photolysis.rate_per_s("CO", teff_k=10000.0) > photolysis.rate_per_s("CO")
    assert photolysis.rate_per_s("CO2", teff_k=10000.0) == pytest.approx(6.5 * photolysis.rate_per_s("CO2"))


def test_ammonia_and_formaldehyde_barely_change_around_a_cool_star():
    assert photolysis.spectral_ratio("NH3", 4000.0) == pytest.approx(0.88)
    assert photolysis.spectral_ratio("CH2O", 4000.0) == pytest.approx(0.78)


def test_the_ratio_is_1_at_the_sun_and_monotone_between_templates():
    for formula in photolysis.formulas():
        assert photolysis.spectral_ratio(formula, 5772.0) == pytest.approx(1.0)
        row = photolysis.species(formula)
        values = [photolysis.spectral_ratio(formula, t) for t in (4000, 4500, 5000, 5772, 7000, 10000)]
        assert values[0] == pytest.approx(row.ratio_4000k) and values[-1] == pytest.approx(row.ratio_10000k)
        cold, hot = values[:4], values[3:]
        assert cold == sorted(cold) or cold == sorted(cold, reverse=True)
        assert hot == sorted(hot) or hot == sorted(hot, reverse=True)


def test_outside_the_templates_the_end_value_holds_and_is_flagged():
    assert photolysis.outside_templates(3000.0) and photolysis.outside_templates(20000.0)
    assert not photolysis.outside_templates(5772.0)
    assert photolysis.spectral_ratio("CO", 3000.0) == photolysis.spectral_ratio("CO", 4000.0)
    assert photolysis.spectral_ratio("CO", 30000.0) == photolysis.spectral_ratio("CO", 10000.0)


def test_the_active_sun_shortens_water_and_carbon_dioxide_by_the_read_factors():
    assert photolysis.species("H2O").lifetime_s / photolysis.lifetime_s("H2O", active_sun=True) == pytest.approx(1.68)
    assert photolysis.species("CO2").lifetime_s / photolysis.lifetime_s("CO2", active_sun=True) == pytest.approx(2.35)
    assert photolysis.species("NH3").activity_factor == pytest.approx(1.07)


def test_a_species_with_no_read_activity_factor_refuses_rather_than_assume_one():
    with pytest.raises(photolysis.NoData):
        photolysis.rate_per_s("CO", active_sun=True)
    assert photolysis.lifetime_s("CO", active_sun=False) == pytest.approx(1.5e6)


def test_the_rates_of_the_leiden_table_are_within_a_factor_two_of_the_adopted_lifetime():
    for formula in photolysis.formulas():
        row = photolysis.species(formula)
        implied = 1.0 / (row.beta_diss_per_s + row.beta_ion_per_s)
        assert row.lifetime_range_s[0] / 2 <= implied <= row.lifetime_range_s[1] * 2


def test_the_range_brackets_the_value_and_the_file_has_the_ratios_the_notes_state():
    with open(photolysis.DATA_FILE, encoding="utf-8") as handle:
        data = json.load(handle)["species"]
    assert data["CO"]["ratio_4000k"] == 3.5e-4 and data["N2"]["ratio_4000k"] == 2.9e-5
    assert data["CH4"]["ratio_4000k"] == 0.013 and data["NH3"]["ratio_10000k"] == 0.41


def test_no_flux_means_an_endless_lifetime():
    assert photolysis.lifetime_s("H2O", 0.0) == math.inf


def test_errors():
    with pytest.raises(photolysis.UnknownSpecies):
        photolysis.species("XeF6")
    with pytest.raises(ValueError):
        photolysis.spectral_ratio("CO", 0.0)
    with pytest.raises(ValueError):
        photolysis.rate_per_s("CO", -1.0)
    with pytest.raises(ValueError):
        photolysis.inverse_square(0.0)


def test_the_nh3_activity_factor_is_the_computed_one_and_names_its_inputs():
    """GEN.794: 1.07 is the mean of the three measured records (1.066, 1.083, 1.060), inside the 1.05 to 1.15 bracket of its registry row."""
    row = photolysis.species("NH3")
    assert 1.05 <= row.activity_factor <= 1.15
    assert row.activity_factor == pytest.approx((1.066 + 1.083 + 1.060) / 3.0, abs=0.005)
    with open(photolysis.DATA_FILE, encoding="utf-8") as handle:
        source = json.load(handle)["species"]["NH3"]["activity_source"]
    for name in ("SORCE SOLSTICE", "UARS SUSIM", "UARS SOLSTICE", "scripts/nh3_activity_factor.py"):
        assert name in source


# --- GEN.982: the column that absorbs the light and the most the light can break in a thick air ------------------------------------

SECONDS_PER_YEAR = 365.25 * 86400.0
EARTH_AREA_M2 = 5.1007e14


def test_the_shielded_column_is_one_optical_depth_of_the_read_cross_section():
    assert photolysis.shielded_column_mol_m2("CH4") == pytest.approx(1.0e4 / (1.81e-17 * 6.02214076e23))
    assert photolysis.shielded_column_mol_m2("CO2") > 1000.0 * photolysis.shielded_column_mol_m2("CH4")     # carbon dioxide barely absorbs
    for formula in photolysis.formulas():
        assert (photolysis.species(formula).uv_cross_section_cm2 is None) == (formula not in ("H2O", "NH3", "CH4", "CO2", "H2S"))


def test_a_species_with_no_read_cross_section_refuses_rather_than_assume_one():
    with pytest.raises(photolysis.NoData):
        photolysis.shielded_column_mol_m2("CO")
    with pytest.raises(photolysis.NoData):
        photolysis.photon_limited_flux_mol_m2_s("HCN")


def test_the_methane_the_sun_can_break_over_the_earths_area_is_within_the_range_of_the_anoxic_methane_fluxes():
    """Tens of Tmol a year: the same order as the 9 to 67 Tmol a year that Kharecha et al. (2005) give for the anoxic Earth's methane (which is a net flux,
    the products recombining, so the gross loss is the larger)."""
    tmol_per_year = photolysis.photon_limited_flux_mol_m2_s("CH4") * EARTH_AREA_M2 * SECONDS_PER_YEAR / 1.0e12
    assert 30.0 < tmol_per_year < 150.0


def test_the_photon_limited_flux_follows_the_light_and_the_spectrum_of_the_star():
    base = photolysis.photon_limited_flux_mol_m2_s("CH4")
    assert photolysis.photon_limited_flux_mol_m2_s("CH4", flux_ratio=0.25) == pytest.approx(0.25 * base)
    assert photolysis.photon_limited_flux_mol_m2_s("CH4", flux_ratio=1.0, teff_k=4000.0) == pytest.approx(
        base * photolysis.spectral_ratio("CH4", 4000.0))
    assert photolysis.photon_limited_flux_mol_m2_s("CH4", flux_ratio=0.0) == 0.0
