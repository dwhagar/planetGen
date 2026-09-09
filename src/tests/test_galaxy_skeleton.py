# tests/test_galaxy_skeleton.py

"""
Tests for `stellarObjects.galaxySkeleton` -- the per-shell candidate-band
finder the galaxy-wide density "skeleton" (`galaxyPlan.py`) is built from.
See that module's own docstring for the reasoning: a shell's qualifying
region is found via an exact upper bound over `theta`, coarse-scanned and
bisected in `phi`, then converted to a slot-index range via the existing
exact `galaxyGeometry.slot_index_bounds_for_phi_range` inverse -- a safe
(possibly wider-than-exact) candidate range, not a per-slot membership
list.

The headline property under test is soundness, not tightness: every slot
that genuinely qualifies (checked via the real, exact
`galaxyDensity.relative_density` -- not the upper bound) must fall inside
one of `find_shell_bands`'s returned bands. A band that were too *narrow*
would silently drop real content; a band that's wider than strictly
necessary just means a few extra live checks at generation time, which
this design already expects and accepts (see the module docstring).
"""

import math

import pytest

from stellarObjects.galaxyDensity import build_galaxy_shape, relative_density
from stellarObjects.galaxyGeometry import sector_position_pc, shell_radius_pc, shell_sector_count
from stellarObjects.galaxySkeleton import (
    bound_relative_density_at_phi,
    expected_system_count_at_density_1,
    find_shell_bands,
    whole_shell_qualifies,
)

SHAPE = build_galaxy_shape(
    disk_scale_length_pc=40.0,
    disk_scale_height_pc=12.0,
    bulge_scale_radius_pc=10.0,
    bulge_amplitude=2.0,
    arm_count=2,
    pitch_angle_rad=math.radians(15),
    arm_amplitude=0.4,
)
EDGE_PC = 3.526


def test_expected_system_count_at_density_1_matches_spacesector():
    from stellarObjects.spaceSector import SpaceSector
    expected = SpaceSector(name="x", edge_ly=11.5).expected_system_count()
    assert expected_system_count_at_density_1(11.5) == pytest.approx(expected)


def test_whole_shell_qualifies_true_when_bulge_alone_clears_threshold():
    # At the galactic center, bulge_amplitude=2.0 -> raw bulge term = 2.0,
    # comfortably above a permissive threshold.
    r_k = shell_radius_pc(0, EDGE_PC)
    assert whole_shell_qualifies(SHAPE, r_k, threshold_rho=0.01)


def test_whole_shell_qualifies_false_far_out_with_strict_threshold():
    r_k = shell_radius_pc(500, EDGE_PC)  # far beyond both scale lengths
    assert not whole_shell_qualifies(SHAPE, r_k, threshold_rho=0.5)


def test_bound_relative_density_at_phi_is_symmetric_about_the_equator():
    r_k = shell_radius_pc(20, EDGE_PC)
    for phi in (0.1, 0.5, 1.0, 1.4):
        a = bound_relative_density_at_phi(SHAPE, r_k, math.pi / 2 - phi)
        b = bound_relative_density_at_phi(SHAPE, r_k, math.pi / 2 + phi)
        assert a == pytest.approx(b)


def test_bound_relative_density_at_phi_peaks_at_the_equator():
    # disk_scale_height (12) << disk_scale_length (40) in this SHAPE, so
    # the equator (phi=pi/2) should dominate over the poles at any shell
    # radius large enough for the disk term to matter at all -- see
    # galaxySkeleton's own module docstring for why this holds generally.
    r_k = shell_radius_pc(20, EDGE_PC)
    equator = bound_relative_density_at_phi(SHAPE, r_k, math.pi / 2)
    for phi in (0.0, 0.3, 0.8, math.pi - 0.3, math.pi):
        assert bound_relative_density_at_phi(SHAPE, r_k, phi) <= equator


def _exact_qualifying_slots(shape, edge_pc, shell_index, threshold_rho):
    """Brute-force ground truth: every slot whose *exact* relative_density
    (not the upper bound) clears threshold_rho."""
    n_k = shell_sector_count(shell_index)
    qualifying = set()
    for i in range(n_k):
        position = sector_position_pc(shell_index, i, edge_pc)
        if relative_density(position, shape) >= threshold_rho:
            qualifying.add(i)
    return qualifying


@pytest.mark.parametrize("shell_index,threshold_rho", [
    (5, 0.05), (5, 0.2), (15, 0.05), (15, 0.15), (30, 0.02), (60, 0.01),
])
def test_find_shell_bands_never_excludes_a_true_qualifier(shell_index, threshold_rho):
    bands = find_shell_bands(SHAPE, EDGE_PC, shell_index, threshold_rho)
    exact_qualifiers = _exact_qualifying_slots(SHAPE, EDGE_PC, shell_index, threshold_rho)

    for slot in exact_qualifiers:
        assert any(lo <= slot <= hi for lo, hi in bands), (
            f"slot {slot} genuinely qualifies (exact check) but falls outside every "
            f"candidate band {bands} for shell {shell_index}"
        )


def test_find_shell_bands_whole_shell_fast_path_returns_full_range():
    n_k = shell_sector_count(0)
    bands = find_shell_bands(SHAPE, EDGE_PC, 0, threshold_rho=0.01)
    assert bands == [(0, n_k - 1)]


def test_find_shell_bands_empty_far_beyond_the_edge():
    # At shell 5000, both the bulge and disk terms are astronomically
    # small for this toy-scale SHAPE (scale lengths of 10-40 pc) -- no
    # reasonable threshold should ever qualify there.
    bands = find_shell_bands(SHAPE, EDGE_PC, 5000, threshold_rho=1e-6)
    assert bands == []


def test_find_shell_bands_returns_bands_in_ascending_order_and_non_overlapping():
    bands = find_shell_bands(SHAPE, EDGE_PC, 15, threshold_rho=0.15)
    for (lo1, hi1), (lo2, hi2) in zip(bands, bands[1:]):
        assert hi1 < lo2 - 1  # merged already if touching/overlapping
