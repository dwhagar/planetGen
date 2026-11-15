"""
TEST.715 (the gate TEST.455 names, from the GEN.945 census): a numeric constant added to the generator's own modules (physics, galaxy, generation, population,
`tuning.py`) has a stand-in registry row, a tag in its own docstring or comment, or a place in the constants module; the committed list of those that have none
(`baselines/unregistered_constants.txt`) may only shrink.
"""

import pathlib
import sys
import warnings

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


def test_a_baseline_name_that_is_registered_or_gone_is_only_a_warning(found, baseline):
    """The rule is that no NEW name appears; a name that has since been registered, tagged, moved or removed leaves the baseline the next time it is written (`scripts/constants_census.py --write`), so a PR that registers a constant need not touch the file."""
    stale = sorted(set(baseline) - found)
    if stale:
        warnings.warn(f"{len(stale)} baseline names no longer qualify; run `python scripts/constants_census.py --write` to drop them: {stale[:5]}", stacklevel=1)
