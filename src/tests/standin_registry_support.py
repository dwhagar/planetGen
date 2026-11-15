# tests/standin_registry_support.py

"""
The checks the stand-in registry tests of the atmosphere draw and the young-group draw share (TEST.735): the rows of one marked block of
`gap-filling-plan.md`, a row that names a constant and a test that exist, and a constant list every member of which has a row.
"""

import importlib
import pathlib
import re

PLAN = pathlib.Path(__file__).resolve().parents[2] / "docs" / "design" / "research" / "gap-filling-plan.md"
TESTS = pathlib.Path(__file__).resolve().parent


def rows_of(block_name):
    """The ten-cell rows of the table between `<!-- {block_name}:begin -->` and `<!-- {block_name}:end -->`."""
    text = PLAN.read_text()
    block = text.split(f"<!-- {block_name}:begin -->")[1].split(f"<!-- {block_name}:end -->")[0]
    rows = [[cell.strip() for cell in line.strip().strip("|").split("|")] for line in block.strip().splitlines()]
    return [row for row in rows[2:] if len(row) == 10]


def assert_row_names_a_constant_and_a_test(row, module_of):
    """`module_of(name)` is the module a constant's module name stands for."""
    fill_id, kind, tier, tweak, evidence, constant, knobs, validity, falsifier, replaced_by = row
    assert kind in ("D", "B", "T", "P")
    assert evidence and knobs and validity and replaced_by
    module, name = constant.split(".")
    assert hasattr(importlib.import_module(module_of(module)), name), constant
    if not falsifier.startswith("none"):
        file, test = falsifier.split("::")
        assert re.search(rf"^def {re.escape(test)}\b", (TESTS / f"{file}.py").read_text(), re.M), falsifier
    else:
        assert len(falsifier) > len("none: ")


def assert_every_stand_in_has_a_row(rows, stand_ins):
    """Every constant of `stand_ins` is named by a row, and the row ids are unique."""
    named = {row[5] for row in rows}
    assert set(stand_ins) <= named
    assert len({row[0] for row in rows}) == len(rows)
