# tests/test_comet_species.py

"""GEN.1005, k7: the species only one comet has measured, from the host's ratios and a trapping factor each (`physics/comet_species.py`). 67P is the calibration, never the value."""

import pytest

from planetgen.physics import comet_species as cs

MEASURED_67P = {"N2": 8.9e-4, "Ar": 5.8e-6, "OCS": 4.1e-4, "SO2": 1.3e-3}      # Rubin et al. 2019, per water


def test_the_carbon_in_co_is_half_of_the_nebulas():
    assert cs.carbon_in_co_share() == pytest.approx(0.52, abs=0.01)


@pytest.mark.parametrize("species", sorted(MEASURED_67P))
def test_a_host_of_solar_composition_gives_the_calibration_comet_within_two_percent(species):
    assert cs.median_over_h2o(species) == pytest.approx(MEASURED_67P[species], rel=0.02)


def test_a_host_with_twice_the_nitrogen_gives_twice_the_n2_and_nothing_else_moves():
    assert cs.median_over_h2o("N2", n_over_c=2.0) == pytest.approx(2.0 * cs.median_over_h2o("N2"))
    assert cs.median_over_h2o("Ar", n_over_c=2.0) == cs.median_over_h2o("Ar")
    assert cs.median_over_h2o("SO2", n_over_c=2.0) == cs.median_over_h2o("SO2")


def test_sulfur_species_follow_the_hosts_sulfur_and_argon_its_argon():
    assert cs.median_over_h2o("OCS", s_over_c=3.0) == pytest.approx(3.0 * cs.median_over_h2o("OCS"))
    assert cs.median_over_h2o("SO2", s_over_c=0.5) == pytest.approx(0.5 * cs.median_over_h2o("SO2"))
    assert cs.median_over_h2o("Ar", ar_over_c=4.0) == pytest.approx(4.0 * cs.median_over_h2o("Ar"))


def test_o2_krypton_and_xenon_have_no_host_ratio_and_keep_one_value():
    for species, value in cs.UNSCALED_OVER_H2O.items():
        assert cs.median_over_h2o(species, n_over_c=3.0, s_over_c=3.0, ar_over_c=3.0) == value
    with pytest.raises(ValueError):
        cs.nebular_per_co("O2")


def test_the_width_sits_inside_the_scatter_of_the_four_species_a_population_measures_and_the_old_design_width():
    """The four k7 species measured in 9 to 18 comets scatter by 0.33 to 0.39 dex; the old width of the other seven was a design 0.6."""
    assert 0.33 <= cs.sigma_dex() <= 0.6


def test_a_host_ratio_that_is_not_above_zero_is_a_refusal():
    with pytest.raises(ValueError):
        cs.median_over_h2o("N2", n_over_c=0.0)


def test_every_species_of_row_k7_is_here():
    assert set(cs.SPECIES) == {"N2", "O2", "OCS", "SO2", "Ar", "Kr", "Xe"}
