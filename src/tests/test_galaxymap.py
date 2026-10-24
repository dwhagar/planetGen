"""
planetgen/web/maps/galaxymap.py regression tests -- this module used to also build
the Galaxy Map's own flat SVG rendering (`render_galaxy_map_panel`),
replaced by a real 3D scene (`planetgen/web/maps/galaxymap3d.py`, tested in
`test_galaxymap3d.py`) that `html/galaxy.py` renders directly now. What's
left here is the plain Quadrant/Zone classification math
(`sector_quadrant`/`sector_zone`/`zone_bounds_ly`) `galaxy.py`'s own data
tables and `sector.py`/`browse.py`'s "Quadrant N" links still depend on --
previously untested in isolation (the old version of this file only ever
exercised them indirectly, through the now-removed panel renderer).

Run with: pytest src/tests/test_galaxymap.py
"""

import pytest  # noqa: E402

from planetgen.web.maps.galaxymap import (  # noqa: E402
    QUADRANT_LABELS,
    ZONE_SHELL_DEPTH,
    ZONE_TARGET_LY,
    sector_quadrant,
    sector_zone,
    zone_bounds_ly,
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


def test_sector_zone_groups_shell_index_into_fixed_depth_bands():
    assert sector_zone(0) == 0
    assert sector_zone(ZONE_SHELL_DEPTH - 1) == 0
    assert sector_zone(ZONE_SHELL_DEPTH) == 1
    assert sector_zone(2 * ZONE_SHELL_DEPTH + 1) == 2


def test_zone_shell_depth_is_at_least_one():
    # ZONE_TARGET_LY / DEFAULT_SECTOR_EDGE_LY could round to 0 for a
    # pathologically large edge length -- max(1, ...) in the module
    # guards against a ZeroDivisionError in sector_zone.
    assert ZONE_SHELL_DEPTH >= 1


def test_zone_bounds_ly_are_contiguous_and_span_zone_target_ly():
    inner_0, outer_0 = zone_bounds_ly(0)
    inner_1, outer_1 = zone_bounds_ly(1)
    assert inner_0 == 0.0
    assert outer_0 == pytest.approx(inner_1)  # zones tile with no gap/overlap
    assert (outer_1 - inner_1) == pytest.approx(ZONE_TARGET_LY, rel=0.15)


def test_zone_bounds_ly_step_by_whole_shells():
    inner_5, outer_5 = zone_bounds_ly(5)
    inner_1, outer_1 = zone_bounds_ly(1)
    width = outer_1 - inner_1
    assert outer_5 - inner_5 == pytest.approx(width)
    assert inner_5 - inner_1 == pytest.approx(4 * width)


def test_sector_zone_and_zone_bounds_ly_agree_with_the_shell_grid():
    # A shell's own zone, fed back into zone_bounds_ly, must contain that
    # shell's whole radial extent -- zones are whole shells.
    from planetgen.galaxy.spherical import shell_bounds_pc
    from planetgen.tuning import DEFAULT_SECTOR_EDGE_LY
    from planetgen.physics.units import ly_to_pc, pc_to_ly

    edge_pc = ly_to_pc(DEFAULT_SECTOR_EDGE_LY)
    for shell_index in (0, 1, ZONE_SHELL_DEPTH, ZONE_SHELL_DEPTH * 3 + 2, 500):
        inner_ly, outer_ly = zone_bounds_ly(sector_zone(shell_index))
        shell_inner_pc, shell_outer_pc = shell_bounds_pc(shell_index, edge_pc)
        assert inner_ly <= pc_to_ly(shell_inner_pc) + 1e-6
        assert pc_to_ly(shell_outer_pc) <= outer_ly + 1e-6
