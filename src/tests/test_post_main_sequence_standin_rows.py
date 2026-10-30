"""
The stand-in registry rows of the stars after the main sequence (GEN.935, docs/design/research/gap-filling-plan.md section 6.51):
every row names a constant that exists and a test that exists, and every constant the post-main-sequence code chooses has a row.
"""

import pathlib
import re

import pytest

from planetgen import tuning
from planetgen.datasets import evolved_tracks

PLAN = pathlib.Path(__file__).resolve().parents[2] / "docs" / "design" / "research" / "gap-filling-plan.md"
TESTS = pathlib.Path(__file__).resolve().parent
NAME = re.compile(r"\b(tuning|evolved_tracks)\.([A-Z][A-Z0-9_]*)\b")
TEST = re.compile(r"\b(test_[a-z0-9_]+)::(test_[a-z0-9_]+)")
MODULES = {"tuning": tuning, "evolved_tracks": evolved_tracks}
CHOSEN = {"tuning.HELIUM_FLASH_MASS_SOL", "tuning.HORIZONTAL_BRANCH_DURATION_GY", "tuning.AGB_MS_FRACTION",
          "tuning.SUPERGIANT_MIN_MASS_SOL", "tuning.BRIGHT_GIANT_MIN_MASS_SOL", "tuning.SUPERGIANT_YERKES_THRESHOLDS_SOL",
          "tuning.CEPHEID_TEMPERATURE_RANGE_K", "tuning.CEPHEID_MIN_LUMINOSITY_SOL", "tuning.CEPHEID_MAX_INITIAL_MASS_SOL",
          "evolved_tracks.MASS_NODES_SOL", "tuning.RED_GIANT_BRANCH_REIMERS_ETA", "tuning.MIST_REIMERS_ETA",
          "tuning.WR_LUMINOSITY_RANGE_SOL"}
"""set: The constants of the post-main-sequence code that are choices, not data or structure."""


def _rows():
    block = PLAN.read_text().split("<!-- post-main-sequence-standins:begin -->")[1].split("<!-- post-main-sequence-standins:end -->")[0]
    rows = [[cell.strip() for cell in line.strip().strip("|").split("|")] for line in block.strip().splitlines()]
    return [row for row in rows[2:] if len(row) == 10]


def test_the_registry_has_a_row_for_every_choice_of_the_code():
    named = {f"{module}.{name}" for row in _rows() for module, name in NAME.findall(row[4])}
    assert CHOSEN <= named, sorted(CHOSEN - named)


@pytest.mark.parametrize("row", _rows(), ids=lambda row: row[0])
def test_a_row_names_constants_and_tests_that_exist(row):
    constants = NAME.findall(row[4])
    assert constants, row[0]
    for module, name in constants:
        assert hasattr(MODULES[module], name), (row[0], module, name)
    for file, test in TEST.findall(row[6] + " " + row[7]):
        assert re.search(rf"^def {test}\(", (TESTS / f"{file}.py").read_text(), re.M), (row[0], file, test)
    assert re.fullmatch(r"T[1-9]", row[2]) and re.fullmatch(r"L[0-4]", row[3]), row[0]
