# tests/test_dry_carbon.py

"""GEN.1005: the carbon of a dry body's mantle from its building blocks and the core's share (`physics/dry_carbon.py`; [R0419]). The bulk silicate Earth is a check, never an input."""

import pytest

from planetgen.physics import dry_carbon as dc

BSE_CARBON_PPM = (100.0, 180.0)      # 140 +/- 40 of the mantle ([R0419], Li et al. 2021); a test point only


def test_a_body_that_took_only_undifferentiated_material_holds_what_the_table_says():
    assert dc.bulk_carbon_ppm(0.0) == dc.NC_UNDIFFERENTIATED_C_PPM == 1600.0
    assert dc.bulk_carbon_ppm(1.0) == dc.NC_DIFFERENTIATED_C_PPM
    assert dc.bulk_carbon_ppm(0.5) == pytest.approx(820.0)


def test_the_bracket_is_the_product_of_the_extremes_in_weight_percent_of_co2():
    low, high = dc.co2_wt_percent_bracket()
    assert low == pytest.approx(dc.mantle_carbon_ppm(0.60, 0.12) * dc.CO2_PER_C / 1e4)
    assert high == pytest.approx(dc.mantle_carbon_ppm(0.30, 0.26) * dc.CO2_PER_C / 1e4)
    assert 0.0 < low < high


def test_the_bracket_is_close_to_the_one_row_b1_drew_and_holds_the_earths_carbon():
    """The row drew 0.018 to 0.1 weight percent CO2, tuned to the bulk silicate Earth; the derived bracket is within 0.6 dex of both ends and holds the carbon of the Earth's mantle (as a test point)."""
    low, high = dc.co2_wt_percent_bracket()
    assert 0.018 / 3.0 < low < 0.018 * 4.0 and 0.1 / 4.0 < high < 0.1 * 3.0
    for ppm_of_mantle in BSE_CARBON_PPM:
        ppm_of_body = ppm_of_mantle * 0.68
        assert low * 1e4 / dc.CO2_PER_C / 2.0 <= ppm_of_body <= high * 1e4 / dc.CO2_PER_C * 2.0


def test_the_mantle_share_and_the_differentiated_share_are_ordered_brackets_inside_zero_and_one():
    for bracket in (dc.DIFFERENTIATED_SHARE, dc.MANTLE_SHARE):
        low, high = bracket
        assert 0.0 < low < high < 1.0


def test_a_share_outside_zero_and_one_is_a_refusal():
    with pytest.raises(ValueError):
        dc.bulk_carbon_ppm(1.5)
    with pytest.raises(ValueError):
        dc.mantle_carbon_ppm(0.5, -0.1)


def test_the_carbon_partitioning_of_the_source_gives_a_siderophile_element_that_the_core_still_leaves_in_the_mantle():
    """Eq. 4 of Blanchard et al.: D_C = 10^(0.3 + 3822/T) x D_Fe with D_Fe = 10^(-dIW/2); at 3800 K and dIW -1.06 that is about 70 (the paper's table has 141 and 8 to 157 across its runs),
    so a fully equilibrated core would take nearly all the carbon, and the mantle share of 0.12 to 0.26 is the incomplete equilibration the paper describes."""
    d_c = 10.0 ** (0.3 + 3822.0 / 3800.0) * 10.0 ** (1.06 / 2.0)
    assert 23 < d_c < 157
    core_fraction = 0.32
    fully_equilibrated_mantle_share = 1.0 / (1.0 + d_c * core_fraction / (1.0 - core_fraction))
    assert fully_equilibrated_mantle_share < min(dc.MANTLE_SHARE)
