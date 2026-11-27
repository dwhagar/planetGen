"""
The stand-in registry of the blue stragglers (GEN.326, docs/design/research/gap-filling-plan.md section 6.331): every row names a constant that
exists and a test that exists, and every constant of the blue straggler block of `tuning.py` has a row.
"""

import pytest

from tests import standin_registry_support as support

MODULES = {"tuning": "planetgen.tuning"}

STANDINS = ("tuning.CLUSTER_BLUE_STRAGGLER_CORE_MASS_EXPONENT", "tuning.CLUSTER_BLUE_STRAGGLERS_AT_1E4_CORE_MASS",
            "tuning.BLUE_STRAGGLER_LIFETIME_GY", "tuning.BLUE_STRAGGLER_MASS_RANGE_SOL")


def _rows():
    return support.rows_of("cluster-blue-straggler-standins")


def test_the_registry_has_a_row_for_each_stand_in():
    assert len(_rows()) >= len(STANDINS)


@pytest.mark.parametrize("row", _rows(), ids=lambda row: row[0])
def test_a_row_names_a_constant_and_a_test_that_exist(row):
    support.assert_row_names_a_constant_and_a_test(row, MODULES.__getitem__)


def test_every_stand_in_constant_has_a_row():
    support.assert_every_stand_in_has_a_row(_rows(), STANDINS)


def test_the_row_ids_are_the_free_bn_ids():
    assert [row[0] for row in _rows()] == [f"BN{i}" for i in range(1, len(_rows()) + 1)]
