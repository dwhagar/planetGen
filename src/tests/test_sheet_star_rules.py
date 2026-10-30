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
