"""
The stand-in registry of the white dwarf binaries (GEN.707, docs/design/research/gap-filling-plan.md section 6.329): every row names a constant that exists and a
test that exists, and every stand-in constant of the GEN.271 block of `tuning.py` that is a draw range or share has a row.
"""

import pytest

from tests import standin_registry_support as support

MODULES = {"tuning": "planetgen.tuning"}

STANDINS = (
    "tuning.CATACLYSMIC_DONOR_MASS_RANGE_SOL", "tuning.AM_CVN_DONOR_MASS_RANGE_SOL", "tuning.ACCRETING_FILL_FACTOR",
    "tuning.HOT_SUBDWARF_POST_CE_SHARE", "tuning.BARIUM_GIANT_MAX_INITIAL_MASS_SOL", "tuning.BARIUM_PERIOD_DAYS", "tuning.SYMBIOTIC_PERIOD_DAYS",
)


def _rows():
    return support.rows_of("white-dwarf-binary-standins")


def test_the_registry_has_a_row_for_each_stand_in():
    assert len(_rows()) >= len(STANDINS)


@pytest.mark.parametrize("row", _rows(), ids=lambda row: row[0])
def test_a_row_names_a_constant_and_a_test_that_exist(row):
    support.assert_row_names_a_constant_and_a_test(row, MODULES.__getitem__)


def test_every_stand_in_constant_has_a_row():
    support.assert_every_stand_in_has_a_row(_rows(), STANDINS)


def test_the_row_ids_are_the_free_wk_ids():
    assert [row[0] for row in _rows()] == [f"WK{i}" for i in range(1, len(_rows()) + 1)]
