# tests/test_sheet_star_rules.py

"""
GEN.788: the sheet columns a white dwarf leaves empty by design (no corona, no flares, no end to its stage) say so in a rule,
and every other star still owes them.
"""

import pytest

from planetgen.describe import sheet

COLUMNS = ("flare_alpha", "uv_c_share", "log_lx_lbol", "lifespan_gy")


@pytest.mark.parametrize("key", COLUMNS)
def test_a_white_dwarf_is_excused_the_activity_columns_and_any_other_star_owes_them(key):
    assert not sheet.applies("star", "F3VII", key, {"yerkes_class": "VII"})
    for yerkes in ("V", "IV", "III", "PMS"):
        assert sheet.applies("star", "G2V", key, {"yerkes_class": yerkes})
    assert sheet.applies("star", "G2V", key)
    assert ("star", key) not in sheet.OPEN_GAPS


def test_the_rule_gives_its_reason():
    [rule] = [r for r in sheet.APPLIES["star"] if "flare_alpha" in r[0]]
    assert "white dwarf" in rule[2] and "corona" in rule[2]


def test_gen789_the_critical_orbit_applies_to_a_wide_pair_and_a_tree_member_but_not_to_a_single_star_or_a_close_pair():
    key = "wide_binary_a_crit_km"
    assert not sheet.applies("star", "G2V", key, {"role": "single"})
    for role in ("primary", "secondary"):
        assert sheet.applies("star", "G2V", key, {"role": role, "binary_configuration": "wide"})
        assert not sheet.applies("star", "G2V", key, {"role": role, "binary_configuration": "close"})
    assert sheet.applies("star", "G2V", key, {"role": "member", "binary_configuration": "wide"})
    assert sheet.applies("star", "G2V", key, {"role": "member", "binary_configuration": "close"})
    assert sheet.applies("star", "G2V", key)
    assert (("star", key)) not in sheet.OPEN_GAPS


def test_gen790_a_comet_off_an_elliptical_orbit_has_no_period_and_one_on_it_owes_them():
    for key in ("mean_anomaly_deg", "orbital_period_years", "period_class"):
        assert sheet.applies("comet", "comet", key, {"orbit_type": "elliptical"})
        for orbit_type in ("parabolic", "hyperbolic"):
            assert not sheet.applies("comet", "comet", key, {"orbit_type": orbit_type})
        assert ("comet", key) not in sheet.OPEN_GAPS
