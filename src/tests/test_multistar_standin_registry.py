"""
The stand-in registry of the young-group draw (GEN.797, docs/design/research/gap-filling-plan.md section 6.11): every row names a
constant that exists and a test that exists, and every stand-in constant of the draw has a row.
"""

import pytest

from tests import standin_registry_support as support

MODULES = {"tuning": "planetgen.tuning", "system_types": "planetgen.generation.phenomena.system_types"}

STANDINS = (
    "system_types.SYSTEM_TYPE_RATES", "tuning.YOUNG_GROUP_AGE_MAX_YR", "tuning.YOUNG_GROUP_DECAY_AGE_RANGE_YR",
    "tuning.YOUNG_GROUP_SIZE_WEIGHTS", "tuning.YOUNG_GROUP_RADIUS_AU_RANGE", "tuning.YOUNG_GROUP_VELOCITY_DISPERSION_KMS",
    "tuning.YOUNG_GROUP_REMNANT_MAX_AXIS_AU",
)


def _rows():
    return support.rows_of("multistar-standins")


def test_the_registry_has_rows():
    assert len(_rows()) >= len(STANDINS)


@pytest.mark.parametrize("row", _rows(), ids=lambda row: row[0])
def test_a_row_names_a_constant_and_a_test_that_exist(row):
    support.assert_row_names_a_constant_and_a_test(row, MODULES.__getitem__)


def test_every_stand_in_constant_has_a_row():
    support.assert_every_stand_in_has_a_row(_rows(), STANDINS)
