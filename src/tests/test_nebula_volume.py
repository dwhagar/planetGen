"""
MAP.365: a nebula group's shape sampled from the model's density (`galaxy/nebula_volume.py`). The tests use a seeded stand-in for
`nebula_density.density_at` (smooth blobs and lattice noise keyed to the absolute position), so they hold the sampler, not the
model: the points sit on one absolute lattice, a nebula is not a ball, two groups that touch draw every place once, the same
group gives the same points, and the object does not matter, only where it is.
"""

import math

import numpy as np

from planetgen.galaxy import nebula_volume as rules
from tests import nebula_density_stand_in as stand_in

FAMILIES = ("diffuse", "emission", "reflection", "planetary", "dark", "supernova-remnant")
FORM = stand_in.FORM_LEVEL
MAX_DETAIL = stand_in.MAX_DETAIL


density = stand_in.density


def _nebula(x, y, z, r, family="dark"):
    return {"x": x, "y": y, "z": z, "radius_pc": r, "family": family}


def _sample(nebulae, inside=None, **kw):
    return rules.sample_group(nebulae, inside, density, FORM, MAX_DETAIL, FAMILIES, **kw)


def _positions(volume):
    return np.stack([volume.ix, volume.iy, volume.iz], axis=1).astype(np.float64) * volume.spacing_pc


def test_levels_follow_the_radius_and_the_lattice_is_a_power_of_two_times_the_base():
    assert rules.spacing_pc(0) == rules.LATTICE_BASE_PC
    assert rules.spacing_pc(3) == 8 * rules.LATTICE_BASE_PC
    assert rules.level_for_radius(1.5) == 0
    assert rules.level_for_radius(48.0) == 5
    assert rules.level_for_radius(1e9) == rules.MAX_LEVEL
    assert rules.level_for_radius(0.0) == rules.MAX_LEVEL


def test_points_sit_on_the_absolute_lattice_inside_the_nebula():
    nebula = _nebula(1003.3, -250.1, 12.0, 30.0)
    volume = _sample([nebula])
    assert len(volume.ix) > 100
    spots = _positions(volume)
    distance = np.linalg.norm(spots - np.array([nebula["x"], nebula["y"], nebula["z"]]), axis=1)
    assert distance.max() <= nebula["radius_pc"] * (1 + 1e-9)
    assert volume.strength.min() >= rules.MIN_STRENGTH
    assert volume.ix.dtype == np.int32 and volume.strength.dtype == np.uint8


def test_a_nebula_is_not_a_ball():
    """The points it holds fill well under a ball's share, and their outline is not a circle: the circularity (4 pi A / P^2) of the
    projection is under a stated bound."""
    nebula = _nebula(0.0, 0.0, 0.0, 60.0)
    volume = _sample([nebula])
    spots = _positions(volume)
    ball_points = rules._lattice_in_spheres([(0.0, 0.0, 0.0, 60.0)], volume.spacing_pc)
    assert len(volume.ix) < 0.8 * len(ball_points)
    grid = np.zeros((41, 41), dtype=bool)
    cell = 2 * 60.0 / 40
    for x, y in spots[:, :2]:
        grid[int((y + 60.0) / cell), int((x + 60.0) / cell)] = True
    area = grid.sum()
    padded = np.pad(grid, 1)
    perimeter = sum((padded & ~np.roll(padded, shift, axis)).sum() for axis in (0, 1) for shift in (1, -1))
    circularity = 4 * math.pi * area / perimeter ** 2
    assert circularity < 0.6, circularity


def test_the_same_group_gives_the_same_points_and_the_object_does_not_matter():
    members = [_nebula(0.0, 0.0, 0.0, 40.0), _nebula(55.0, 10.0, -5.0, 30.0)]
    first, second = _sample(members), _sample(members)
    assert all(np.array_equal(a, b) for a, b in zip(first[3:], second[3:]))
    renamed = [dict(m, id="other-object", name="Renamed") for m in members]
    third = _sample(renamed)
    assert all(np.array_equal(a, b) for a, b in zip(first[3:], third[3:]))


def test_two_blocks_that_share_a_nebula_draw_each_place_once_and_leave_no_gap():
    """A nebula crossing the border x = 20 is drawn half by each block; together the halves are the whole, and a point has the same
    value whichever block asked for it."""
    nebulae = [_nebula(20.0, 0.0, 0.0, 40.0), _nebula(75.0, 5.0, 0.0, 30.0)]
    west = _sample(nebulae, inside=lambda p: p[:, 0] < 20.0, level=3)
    east = _sample(nebulae, inside=lambda p: p[:, 0] >= 20.0, level=3)
    whole = _sample(nebulae, level=3)
    cells = lambda v: {(a, b, c): s for a, b, c, s in zip(v.ix.tolist(), v.iy.tolist(), v.iz.tolist(), v.strength.tolist())}
    assert len(west.ix) and len(east.ix)
    assert not (cells(west).keys() & cells(east).keys())
    assert {**cells(west), **cells(east)} == cells(whole)


def test_a_finer_level_keeps_the_coarse_blobs_where_they_were():
    nebula = _nebula(0.0, 0.0, 0.0, 28.0)
    coarse, fine = _sample([nebula], level=4), _sample([nebula], level=3)
    assert len(fine.ix) > len(coarse.ix)
    centre_coarse = _positions(coarse).mean(axis=0)
    centre_fine = _positions(fine).mean(axis=0)
    assert np.linalg.norm(centre_coarse - centre_fine) < 4.0


def test_a_shell_family_is_a_ring_not_a_filled_ball():
    ball = _sample([_nebula(0.0, 0.0, 0.0, 40.0, "dark")], level=3)
    shell = _sample([_nebula(0.0, 0.0, 0.0, 40.0, "supernova-remnant")], level=3)
    inner = lambda v: (np.linalg.norm(_positions(v), axis=1) < 12.0).sum() / max(1, len(v.ix))
    assert inner(shell) < 0.3 * inner(ball)


def test_a_group_with_no_member_is_empty_and_pack_round_trips():
    empty = _sample([])
    assert len(empty.ix) == 0
    volume = _sample([_nebula(10.0, 20.0, 30.0, 25.0)])
    blob = rules.pack(volume)
    assert blob[:4] == rules.MAGIC and blob[4] == rules.PACK_VERSION and blob[5] == volume.level
    n = int(np.frombuffer(blob[6:10], dtype="<u4")[0])
    assert n == len(volume.ix)
    assert len(blob) == 14 + n * (3 * 4 + 2)
    ix = np.frombuffer(blob[14:14 + 4 * n], dtype="<i4")
    assert np.array_equal(ix, volume.ix)


def test_big_groups_are_raised_to_a_coarser_level_to_hold_the_point_budget():
    many = [_nebula(float(i * 3), 0.0, 0.0, 30.0) for i in range(200)]
    volume = _sample(many, level=0)
    assert volume.level > 0


def test_inside_block_follows_the_ladders_cut_to_within_a_sector():
    from planetgen import tuning
    from planetgen.db import nebula_volumes as db_volumes, nodes
    from planetgen.galaxy import drill, spherical
    from planetgen.galaxy.nodecut import block_box

    edge_pc = float(tuning.DEFAULT_SECTOR_EDGE_PC)
    block = drill.DrillBlock(27, 3, 0, 5)
    lo, hi = block_box(block, edge_pc)
    rng = np.random.default_rng(5)
    points = rng.uniform(lo, hi, (400, 3))
    mask = db_volumes.inside_block(block, edge_pc)(points)
    memo = {}
    expected = [nodes.address_in_block(spherical.sector_address_at(tuple(p), edge_pc), block, memo) for p in points]
    # The cut is the block's analytic bounds; an object is in the block of its sector, a few parsecs of difference at the border.
    assert sum(a != b for a, b in zip(mask.tolist(), expected)) <= 0.03 * len(points)
    assert 0 < mask.sum() < len(points)
