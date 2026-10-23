# tests/test_photolysis.py

"""GEN.496: the photolysis lifetime table, against docs/design/body-composition-parts/G1-comets.md sections G1.2 and 2.2."""

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
    assert 1.05 <= photolysis.species("NH3").activity_factor <= 1.10


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
