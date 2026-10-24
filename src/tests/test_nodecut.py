# tests/test_nodecut.py

"""
MAP.261: the node cut (`planetgen.galaxy.nodecut`), which picks the sector blocks a camera needs by their size on
screen, and the weights that hand a block's place over to its children. The page runs the twin
(`static/galaxynodecut.js`, tests/js/galaxynodecut.test.mjs reads the cases from here), so the numbers matter.
"""

import math

import pytest

from planetgen import tuning
from planetgen.galaxy import nodecut
from planetgen.galaxy.drill import DRILL_LEVELS, DrillBlock, drill_chain_of, format_drill_key, parse_node_key
from planetgen.galaxy.spherical import sector_address_at

EDGE = float(tuning.DEFAULT_SECTOR_EDGE_PC)
FOCAL = nodecut.DEFAULT_FOCAL_PX
BOUNDS = ((-30000.0, -30000.0, -3000.0), (30000.0, 30000.0, 3000.0))

CAMERAS = (
    ((34000.0, -20000.0, 28000.0), (0.0, 0.0, 0.0)),    # the opening view, out of the galaxy
    ((8000.0, 0.0, 0.0), (8000.0, 0.0, -100.0)),        # in the disc, looking down
    ((8012.0, 21.0, 9.0), (8000.0, 4.0, 4.0)),          # close to a point
    ((0.0, 0.0, 0.0), (1.0, 1.0, 0.0)),                 # at the centre
    ((-5000.0, 3000.0, 700.0), (-5100.0, 3000.0, 700.0)),
)


def _cut(camera, look, **kwargs):
    return nodecut.node_cut(camera, look, EDGE, FOCAL, kwargs.pop("bounds_pc", BOUNDS), **kwargs)


def test_the_cut_holds_the_roots_and_every_parent_of_what_it_holds():
    for camera, look in CAMERAS:
        cut = _cut(camera, look)
        keys = set(cut.keys)
        roots = {format_drill_key(block) for block in nodecut.roots()}
        assert roots & keys, "a root is always in view of a camera inside the bounds"
        for key in keys:
            parent = nodecut.drill_parent(parse_node_key(key))
            assert parent is None or format_drill_key(parent) in keys, key
        assert cut.expanded <= keys


def test_the_cut_is_bounded_and_finer_near_than_far():
    for camera, look in CAMERAS:
        cut = _cut(camera, look)
        assert len(cut.keys) <= nodecut.MAX_NODES
        levels = {parse_node_key(key).m for key in cut.keys}
        assert len(levels) > 1 and min(levels) <= 2187, "near space is cut finer than the roots"
    close = _cut(*CAMERAS[2])
    assert min(parse_node_key(key).m for key in close.keys) < 243


def test_a_smaller_budget_gives_a_smaller_cut_that_is_still_closed():
    cut = _cut(*CAMERAS[1], max_nodes=100)
    assert len(cut.keys) <= 100
    for key in cut.keys:
        parent = nodecut.drill_parent(parse_node_key(key))
        assert parent is None or format_drill_key(parent) in cut.keys


def test_nothing_outside_the_bounds_is_in_the_cut():
    bounds = ((-200.0, -200.0, -200.0), (200.0, 200.0, 200.0))
    cut = _cut((100.0, 0.0, 0.0), (0.0, 0.0, 0.0), bounds_pc=bounds, half_angle_deg=180.0)
    for key in cut.keys:
        lo, hi = nodecut.block_box(parse_node_key(key), EDGE)
        assert all(hi[a] > bounds[0][a] and lo[a] <= bounds[1][a] for a in range(3)), key


def test_the_handover_runs_from_70_to_140_pixels_and_a_sector_gives_nothing():
    block = DrillBlock(729, 0, 0, 0)
    assert nodecut.handover(60.0, block) == 0.0
    assert nodecut.handover(70.0, block) == 0.0
    assert nodecut.handover(105.0, block) == pytest.approx(0.5)
    assert nodecut.handover(140.0, block) == 1.0
    assert nodecut.handover(500.0, DrillBlock(1, 0, 0, 0)) == 0.0


def test_the_weights_of_the_blocks_over_a_point_add_up_to_one():
    """A star in a block and in its children shows once from each, the two shares adding up to one."""
    for camera, look in CAMERAS:
        cut = _cut(camera, look, half_angle_deg=180.0, bounds_pc=None)
        weigh = nodecut.weight_calculator(camera, EDGE, FOCAL, cut.expanded)
        held = set(cut.keys)
        for point in ((8000.0, 100.0, 5.0), (-3000.0, 4000.0, 20.0), (50.0, 20.0, 1.0), (15000.0, -2000.0, 100.0)):
            chain = [block for block in drill_chain_of(*sector_address_at(point, EDGE)) if format_drill_key(block) in held]
            assert chain, point
            assert sum(weigh(block) for block in chain) == pytest.approx(1.0, abs=1e-9), (camera, point)


def test_a_block_the_cut_left_shut_stays_whole_and_its_children_show_nothing():
    camera, look = CAMERAS[1]
    cut = _cut(camera, look)
    shut = next(parse_node_key(key) for key in cut.keys if key not in cut.expanded and parse_node_key(key).m > 1)
    weigh = nodecut.weight_calculator(camera, EDGE, FOCAL, cut.expanded)
    assert weigh(shut) > 0.0
    for _shell, group in nodecut.drill_children(shut)[:1]:
        assert all(weigh(child) == 0.0 for child in group), "a block that gives nothing leaves its children no share"


def test_a_block_over_140_pixels_has_handed_everything_to_its_children():
    camera = (8000.0, 0.0, 0.0)
    block = next(block for block in nodecut.roots())
    size = nodecut.projected_size(block, camera, EDGE, FOCAL)
    assert size > 140.0
    expanded = {format_drill_key(block)}
    assert nodecut.node_weight(block, camera, EDGE, FOCAL, expanded) == 0.0
    assert nodecut.node_weight(block, camera, EDGE, FOCAL, set()) == 1.0


def test_the_projected_size_falls_with_distance():
    block = DrillBlock(DRILL_LEVELS[3], 0, 0, 0)
    near = nodecut.projected_size(block, (100.0, 0.0, 0.0), EDGE, FOCAL)
    far = nodecut.projected_size(block, (100000.0, 0.0, 0.0), EDGE, FOCAL)
    assert near > far > 0
    assert math.isfinite(far)
