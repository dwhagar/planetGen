"""
Every stand-in of the accreting compact binaries block of `tuning.py` has a row in the stand-in registry (GEN.724), and the
registry row `millisecond-pulsar-companion` keeps the neutron star companion share the code keeps.
"""

import ast
import pathlib
import re

import pytest

from planetgen.generation import multistar
from planetgen.generation.phenomena import system_types

ROOT = pathlib.Path(__file__).resolve().parents[2]
TUNING = ROOT / "src" / "planetgen" / "tuning.py"
PLAN = ROOT / "docs" / "design" / "research" / "gap-filling-plan.md"

BLOCK_OPENS = "# Accreting compact binaries (GEN.272"
BLOCK_ENDS = "WIND_MASS_LOSS_MSUN_YR_AT_1E5_LSUN"


def _block_constants():
    text = TUNING.read_text(encoding="utf-8")
    lines = text.splitlines()
    first = next(i for i, line in enumerate(lines, 1) if line.startswith(BLOCK_OPENS))
    last = next(i for i, line in enumerate(lines, 1) if line.startswith(BLOCK_ENDS))
    return [node.targets[0].id for node in ast.parse(text).body
            if isinstance(node, ast.Assign) and first < node.lineno < last and isinstance(node.targets[0], ast.Name)]


def _registry_section():
    text = PLAN.read_text(encoding="utf-8")
    start = text.index("### 6.23 ")
    return text[start:text.index("\n## ", start)]


def test_the_block_is_found_and_holds_the_constants_of_the_item():
    names = _block_constants()
    assert {"SUPERGIANT_XRAY_PERIOD_DAYS", "BE_XRAY_PERIOD_DAYS", "LOW_MASS_XRAY_PERIOD_DAYS",
            "LOW_MASS_XRAY_DONOR_MASS_RANGE_SOL", "XRAY_SPAN_ROOM", "SPIDER_PERIOD_DAYS",
            "MICROQUASAR_SHARE_OF_XRAY_BINARIES", "ACCRETION_EFFICIENCY_BLACK_HOLE"} <= set(names)


@pytest.mark.parametrize("name", _block_constants())
def test_a_tuning_constant_of_the_block_has_a_registry_row(name):
    assert re.search(rf"tuning\.{name}\b", _registry_section()), f"{name} has no row in section 6.23 of the gap-filling plan"


def test_the_pulsar_row_keeps_the_neutron_star_share_the_code_keeps():
    row = system_types.SYSTEM_TYPE_RATES["millisecond-pulsar-companion"]
    code = multistar.MEMBER_STATE_BY_PLACED["recycled-pulsar"]["neutron-star"][0]
    assert row.default * row.by_class["neutron-star"] == pytest.approx(code, rel=0.05)
    assert row.by_class["neutron-star"] < 0.01
    assert sum(row.by_class.values()) == pytest.approx(1.0, abs=0.001)
