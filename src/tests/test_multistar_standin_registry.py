"""
The stand-in registry of the young-group draw (GEN.797, docs/design/research/gap-filling-plan.md section 6.11): every row names a
constant that exists and a test that exists, and every stand-in constant of the draw has a row.
"""

import importlib
import pathlib
import re

import pytest

PLAN = pathlib.Path(__file__).resolve().parents[2] / "docs" / "design" / "research" / "gap-filling-plan.md"
TESTS = pathlib.Path(__file__).resolve().parent

MODULES = {"tuning": "planetgen.tuning", "system_types": "planetgen.generation.phenomena.system_types"}

STANDINS = (
    "system_types.SYSTEM_TYPE_RATES", "tuning.YOUNG_GROUP_AGE_MAX_YR", "tuning.YOUNG_GROUP_DECAY_AGE_RANGE_YR",
    "tuning.YOUNG_GROUP_SIZE_WEIGHTS", "tuning.YOUNG_GROUP_RADIUS_AU_RANGE", "tuning.YOUNG_GROUP_VELOCITY_DISPERSION_KMS",
    "tuning.YOUNG_GROUP_REMNANT_MAX_AXIS_AU",
)


def _rows():
    text = PLAN.read_text()
    block = text.split("<!-- multistar-standins:begin -->")[1].split("<!-- multistar-standins:end -->")[0]
    rows = [[cell.strip() for cell in line.strip().strip("|").split("|")] for line in block.strip().splitlines()]
    return [row for row in rows[2:] if len(row) == 10]


def test_the_registry_has_rows():
    assert len(_rows()) >= len(STANDINS)


@pytest.mark.parametrize("row", _rows(), ids=lambda row: row[0])
def test_a_row_names_a_constant_and_a_test_that_exist(row):
    fill_id, kind, tier, tweak, evidence, constant, knobs, validity, falsifier, replaced_by = row
    assert kind in ("D", "B", "T", "P")
    assert evidence and knobs and validity and replaced_by
    module, name = constant.split(".")
    assert hasattr(importlib.import_module(MODULES[module]), name), constant
    if not falsifier.startswith("none"):
        file, test = falsifier.split("::")
        assert re.search(rf"^def {re.escape(test)}\b", (TESTS / f"{file}.py").read_text(), re.M), falsifier
    else:
        assert len(falsifier) > len("none: ")


def test_every_stand_in_constant_has_a_row():
    named = {row[5] for row in _rows()}
    assert set(STANDINS) <= named
    assert len({row[0] for row in _rows()}) == len(_rows())
