# tests/test_comet_species.py

"""GEN.1005, GEN.1021, k7: the species few comets have measured, each between the comets that did (`physics/comet_species.py`). The comets measured are the tests, never the inputs."""

import pytest

from planetgen.physics import comet_species as cs


def test_the_carbon_in_co_is_half_of_the_nebulas():
    assert cs.carbon_in_co_share() == pytest.approx(0.52, abs=0.01)


def test_o2_krypton_and_xenon_have_no_host_ratio():
    for species in ("O2", "Kr", "Xe"):
        with pytest.raises(ValueError):
            cs.nebular_per_co(species)


def test_every_species_of_row_k7_is_here():
    assert set(cs.SPECIES) == {"N2", "O2", "OCS", "SO2", "Ar", "Kr", "Xe"}


# --- GEN.1021: the brackets between the comets measured --------------------------------------------------------

BRACKET_SPECIES = ("N2", "Ar", "OCS", "SO2", "O2", "Kr", "Xe")


@pytest.mark.parametrize("species", BRACKET_SPECIES)
def test_every_bracket_is_a_proper_range_above_zero(species):
    low, high = cs.bracket_over_h2o(species)
    assert 0.0 < low < high
    assert (cs.low_over_h2o(species), cs.high_over_h2o(species)) == (low, high)


@pytest.mark.parametrize("species, measured", [("N2", 8.9e-4), ("Ar", 5.8e-6), ("SO2", 1.3e-3), ("OCS", 7e-4), ("O2", 3.8e-2), ("O2", 3.7e-2)])
def test_a_comet_measured_falls_inside_its_species_bracket(species, measured):
    """Tests against measurement, never inputs: 67P (N2, argon, SO2, OCS, O2), Halley (O2)."""
    low, high = cs.bracket_over_h2o(species)
    assert low <= measured <= high


def test_hale_bopps_argon_and_c2016_r2s_nitrogen_fall_inside_and_reach_the_ends():
    low, high = cs.bracket_over_h2o("Ar")
    assert low <= 8.7e-3 <= high * 1.0001 and high / low > 1e3                      # 3.2 decades per water, the widest range
    n2_low, n2_high = cs.bracket_over_h2o("N2")
    assert n2_low == pytest.approx(0.0287 * cs.CO_OVER_H2O_MEDIAN) and n2_high == pytest.approx(0.09 * cs.CO_OVER_H2O_MEDIAN)
    assert n2_high / n2_low == pytest.approx(0.09 / 0.0287)


def test_argon_follows_the_water_not_the_co():
    """Hale-Bopp's argon per CO (0.043) is 2.7 times what the nebula holds per CO, so a share of the CO chain cannot reach it."""
    per_co_at_hale_bopp = 8.7e-3 / 0.2
    assert per_co_at_hale_bopp > 2.0 * cs.nebular_per_co("Ar")
    assert cs.high_over_h2o("Ar") > cs.nebular_per_co("Ar") * cs.CO_OVER_H2O_MEDIAN * 10


def test_krypton_and_xenon_are_argon_times_the_ratios_of_67p():
    ar_low, ar_high = cs.bracket_over_h2o("Ar")
    kr_low, kr_high = cs.bracket_over_h2o("Kr")
    xe_low, xe_high = cs.bracket_over_h2o("Xe")
    assert (kr_low, kr_high) == pytest.approx((ar_low * 0.089, ar_high * 0.089))
    assert (xe_low, xe_high) == pytest.approx((kr_low * 0.269, kr_high * 0.269))


def test_the_brackets_follow_the_hosts_ratios_and_o2_follows_nothing():
    assert cs.bracket_over_h2o("N2", n_over_c=2.0) == pytest.approx(tuple(2.0 * v for v in cs.bracket_over_h2o("N2")))
    assert cs.bracket_over_h2o("Ar", ar_over_c=0.5) == pytest.approx(tuple(0.5 * v for v in cs.bracket_over_h2o("Ar")))
    assert cs.bracket_over_h2o("Xe", ar_over_c=3.0) == pytest.approx(tuple(3.0 * v for v in cs.bracket_over_h2o("Xe")))
    assert cs.bracket_over_h2o("N2", ar_over_c=9.0, s_over_c=9.0) == cs.bracket_over_h2o("N2")
    assert cs.bracket_over_h2o("O2", n_over_c=3.0, ar_over_c=3.0, s_over_c=3.0) == cs.bracket_over_h2o("O2")


@pytest.mark.parametrize("species", ("OCS", "SO2"))
def test_a_sulfur_species_never_holds_more_than_the_nebular_sulfur_of_its_host(species):
    for s_over_c in (0.1, 1.0, 5.0):
        low, high = cs.bracket_over_h2o(species, s_over_c=s_over_c)
        assert high <= cs.nebular_per_co(species, s_over_c=s_over_c) * cs.CO_OVER_H2O_MEDIAN * (1 + 1e-12)
        assert 0.0 < low < high


def test_a_host_ratio_that_is_not_above_zero_or_a_species_not_here_is_a_refusal_for_the_brackets_too():
    with pytest.raises(ValueError):
        cs.bracket_over_h2o("Ar", ar_over_c=0.0)
    with pytest.raises(ValueError):
        cs.bracket_over_h2o("CH4")
