"""
html/lib/galaxymap.py regression tests -- this module used to also build
the Galaxy Map's own flat SVG rendering (`render_galaxy_map_panel`),
replaced by a real 3D scene (`html/lib/galaxymap3d.py`, tested in
`test_galaxymap3d.py`) that `html/galaxy.py` renders directly now. What's
left here is the plain Quadrant/Ring classification math
(`sector_quadrant`/`sector_ring`/`ring_bounds_ly`) `galaxy.py`'s own data
tables and `sector.py`/`browse.py`'s "Quadrant N" links still depend on --
previously untested in isolation (the old version of this file only ever
exercised them indirectly, through the now-removed panel renderer).

Run with: pytest src/tests/test_galaxymap.py
"""
import math
import os
import sys

_SRC_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_SRC_DIR, "html", "lib"))

import pytest  # noqa: E402

from galaxymap import (  # noqa: E402
    QUADRANT_LABELS,
    RING_SHELL_WIDTH,
    RING_TARGET_LY,
    ring_bounds_ly,
    sector_quadrant,
    sector_ring,
)


def test_quadrant_labels_are_the_four_roman_numerals():
    assert QUADRANT_LABELS == ("I", "II", "III", "IV")


@pytest.mark.parametrize("x_pc,y_pc,expected", [
    (10.0, 1.0, "I"),      # just above +x, theta just past 0
    (1.0, 10.0, "I"),      # just short of +y, theta just under 90deg
    (-1.0, 10.0, "II"),
    (-10.0, -1.0, "III"),
    (1.0, -10.0, "IV"),
    (1.0, 0.0, "I"),       # theta == 0 exactly, the I/IV boundary's own side
])
def test_sector_quadrant_classifies_by_azimuth_from_plus_x(x_pc, y_pc, expected):
    assert sector_quadrant(x_pc, y_pc) == expected


def test_sector_quadrant_is_blind_to_z_by_construction():
    # sector_quadrant only ever takes x/y -- there's no z parameter to
    # pass in the first place, confirming the "purely azimuthal" design
    # the module docstring describes (a regression here would be a
    # different function signature, not a wrong return value).
    import inspect
    assert list(inspect.signature(sector_quadrant).parameters) == ["x_pc", "y_pc"]


def test_sector_ring_groups_shell_index_into_fixed_width_bands():
    assert sector_ring(0) == 0
    assert sector_ring(RING_SHELL_WIDTH - 1) == 0
    assert sector_ring(RING_SHELL_WIDTH) == 1
    assert sector_ring(2 * RING_SHELL_WIDTH + 1) == 2


def test_ring_shell_width_is_at_least_one():
    # RING_TARGET_LY / DEFAULT_SECTOR_EDGE_LY could round to 0 for a
    # pathologically large edge length -- max(1, ...) in the module
    # guards against a ZeroDivisionError in sector_ring.
    assert RING_SHELL_WIDTH >= 1


def test_ring_bounds_ly_are_contiguous_and_span_ring_target_ly():
    inner_0, outer_0 = ring_bounds_ly(0)
    inner_1, outer_1 = ring_bounds_ly(1)
    assert inner_0 == 0.0
    assert outer_0 == pytest.approx(inner_1)  # rings tile with no gap/overlap
    assert outer_0 == pytest.approx(RING_SHELL_WIDTH * (outer_0 / RING_SHELL_WIDTH))
    # Each ring's own width should land close to RING_TARGET_LY by
    # construction (RING_SHELL_WIDTH is derived to make this so).
    assert (outer_0 - inner_0) == pytest.approx(RING_TARGET_LY, rel=0.15)


def test_ring_bounds_ly_scale_linearly_with_ring_index():
    inner_5, outer_5 = ring_bounds_ly(5)
    inner_1, outer_1 = ring_bounds_ly(1)
    width = outer_1 - inner_1
    assert inner_5 == pytest.approx(5 * width)
    assert outer_5 == pytest.approx(6 * width)


def test_sector_ring_and_ring_bounds_ly_agree_with_each_other():
    # A shell_index's own ring, fed back into ring_bounds_ly, must bound
    # that same shell's real radial position -- these two functions
    # describe the same tiling from two directions and must stay
    # consistent with each other.
    from stellarObjects.galaxyGeometry import shell_radius_pc
    from stellarObjects.program_constants import DEFAULT_SECTOR_EDGE_LY
    from stellarObjects.utils import ly_to_pc, pc_to_ly

    edge_pc = ly_to_pc(DEFAULT_SECTOR_EDGE_LY)
    for shell_index in (0, 1, RING_SHELL_WIDTH, RING_SHELL_WIDTH * 3 + 2, 500):
        ring = sector_ring(shell_index)
        inner_ly, outer_ly = ring_bounds_ly(ring)
        shell_radius_ly = pc_to_ly(shell_radius_pc(shell_index, edge_pc))
        assert inner_ly <= shell_radius_ly <= outer_ly
