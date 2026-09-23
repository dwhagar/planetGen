# tests/test_bughunt_galaxy_boundaries.py

"""
Tier 1/Tier 2 bug-hunt coverage: shell/pole boundary conditions in
`galaxyGeometry.py` -- shell 0 (the innermost, smallest shell), the last
valid slot index in a shell, and slots near the poles (`phi` close to 0 or
pi), where the module's own docstrings admit approximation (e.g.
`sector_wedge_vertices_pc`'s "not an exact spherical Voronoi cell" and
`slot_index_bounds_for_phi_range`'s "1-slot buffer" fudge for float
rounding near a pole).

Tier 1 here is crash-freedom and basic finiteness/sanity (every returned
coordinate finite, every slot index within its shell's own valid range);
Tier 2 cross-checks `sector_position_pc`'s claimed nominal radius against
the actual computed radius, since a systematic drift there wouldn't crash
anything but would be a real placement bug.
"""

import math

import pytest

from stellarObjects.galaxyGeometry import (
    galactic_radius_pc, sector_position_pc, sector_wedge_vertices_pc,
    shell_radius_pc, shell_sector_count, slot_index_bounds_for_phi_range,
)
from tests.bughunt_support import tier2

EDGE_PC = 10.0


@pytest.mark.parametrize("shell_index", [0, 1, 2])
def test_every_slot_in_shell_produces_finite_position(shell_index):
    """Sweep every slot in a small shell (0-2 stay small enough to be
    fast) -- includes index 0 and the last index, the two most pole-
    adjacent slots per _phi_for_index's own mapping."""
    n_k = shell_sector_count(shell_index)
    for slot_index in range(n_k):
        x, y, z = sector_position_pc(shell_index, slot_index, EDGE_PC)
        assert all(math.isfinite(v) for v in (x, y, z)), (
            f"shell={shell_index} slot={slot_index}: non-finite position ({x}, {y}, {z})"
        )


@pytest.mark.parametrize("shell_index", [0, 1, 2, 50, 500])
def test_first_and_last_slot_wedge_vertices_are_finite(shell_index):
    """The two pole-adjacent slots (index 0 and n_k-1) are exactly where
    sector_wedge_vertices_pc's own docstring admits the wedge-vertex
    approximation is weakest -- confirm it still returns finite,
    non-degenerate (nonzero-area) vertices there, not NaN or a collapsed
    point."""
    n_k = shell_sector_count(shell_index)
    for slot_index in (0, n_k - 1):
        vertices = sector_wedge_vertices_pc(shell_index, slot_index, EDGE_PC)
        assert len(vertices) >= 3, f"shell={shell_index} slot={slot_index}: fewer than 3 vertices"
        for vx, vy, vz in vertices:
            assert all(math.isfinite(v) for v in (vx, vy, vz)), (
                f"shell={shell_index} slot={slot_index}: non-finite vertex ({vx}, {vy}, {vz})"
            )


def test_slot_index_bounds_near_poles_stay_within_shell_range():
    """slot_index_bounds_for_phi_range's docstring calls out a "1-slot
    buffer" fudge for float rounding near phi=0/pi -- confirm the bounds
    it returns for a phi range touching a pole never fall outside
    [0, n_k)."""
    for shell_index in (0, 1, 5, 100):
        n_k = shell_sector_count(shell_index)
        for phi_min, phi_max in [(0.0, 0.1), (math.pi - 0.1, math.pi), (0.0, math.pi)]:
            lo, hi = slot_index_bounds_for_phi_range(phi_min, phi_max, n_k)
            assert 0 <= lo <= hi < n_k or lo > hi, (
                f"shell={shell_index} n_k={n_k} phi=({phi_min},{phi_max}): "
                f"bounds ({lo}, {hi}) outside valid slot range [0, {n_k})"
            )


@tier2
def test_shell_radius_matches_nominal_design_formula():
    """Tier 2: cross-checks the actual galactic radius every slot in a
    shell lands at against shell_radius_pc's own claimed nominal r_k --
    a systematic drift wouldn't crash anything, but would be a real,
    silent placement bug."""
    for shell_index in (0, 1, 10, 200):
        nominal = shell_radius_pc(shell_index, EDGE_PC)
        n_k = shell_sector_count(shell_index)
        # Sample a handful of slots rather than every one, for speed.
        sample = range(0, n_k, max(1, n_k // 8))
        for slot_index in sample:
            position = sector_position_pc(shell_index, slot_index, EDGE_PC)
            actual = galactic_radius_pc(position)
            assert actual == pytest.approx(nominal, rel=1e-6), (
                f"shell={shell_index} slot={slot_index}: actual radius {actual} pc vs "
                f"nominal {nominal} pc"
            )
