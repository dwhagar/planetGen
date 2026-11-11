"""
TEST.715 (the gate TEST.455 names, from the GEN.945 census): a numeric constant added to the generator's own modules (physics, galaxy, generation, population,
`tuning.py`) has a stand-in registry row, a tag in its own docstring or comment, or a place in the constants module; the committed list of those that have none
(`baselines/unregistered_constants.txt`) may only shrink.
"""

import pathlib
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import constants_census as census  # noqa: E402


@pytest.fixture(scope="module")
def found():
    return set(census.unregistered())


@pytest.fixture(scope="module")
def baseline():
    return census.read_baseline()


def test_the_baseline_is_sorted_and_has_no_repeats(baseline):
    assert baseline == sorted(set(baseline))


def test_no_new_constant_lacks_a_registry_row_a_tag_or_a_place_in_the_constants_module(found, baseline):
    new = sorted(found - set(baseline))
    assert not new, ("numeric constants with no registry row (docs/design/research/stand-in-values.md), no tag in their docstring (stand-in, tier, [R....], derived, measured) "
                     f"and not in a constants module: {new}")


def test_a_baseline_name_that_is_registered_or_gone_is_removed_from_the_baseline(found, baseline):
    stale = sorted(set(baseline) - found)
    assert not stale, f"delete these lines from src/tests/baselines/unregistered_constants.txt (registered, tagged, moved or removed): {stale}"
