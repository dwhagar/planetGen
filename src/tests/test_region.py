"""MAP.151: the round region as layer, ring and wrapping slot-arc ranges."""
import math

import pytest

from planetgen.galaxy.geometry import (
    enumerate_sectors_within_radius, ring_sector_count, sector_address_at, sector_position_pc,
)
from planetgen.galaxy.region import CYLINDER_EDGE_ALLOWANCE, Region
from planetgen.galaxy.skeleton import GalaxyBounds

EDGE = 4.0


def _walked(region):
    """The region by brute force: every sector of the sphere round it, cut to the can."""
    reach = math.hypot(region.horizontal_pc, (region.half_layers + 0.5) * region.edge_pc)
    low, high = region.layer_range
    cx, cy = region.center_pc[0], region.center_pc[1]
    return {
        (ring, layer, slot)
        for ring, layer, slot, x, y, _z, _d in enumerate_sectors_within_radius(region.center_pc, reach, region.edge_pc)
        if low <= layer <= high and math.hypot(x - cx, y - cy) <= region.horizontal_pc + 1e-9
    }


CENTRES = [
    (0, 0, 0),      # the axis
    (1, 0, 1),      # beside it
    (2, 1, 0),
    (12, 0, 5),
    (40, -2, 0),    # slot 0 on a bigger ring
    (40, 3, 63),
    (150, 0, 400),
]


@pytest.mark.parametrize("centre", CENTRES)
@pytest.mark.parametrize("sectors,layers", [(0.5, 0), (1, 1), (2, 1), (5, 2), (13, 3), (30, 0)])
def test_a_region_is_the_sectors_a_walk_of_the_sphere_finds(centre, sectors, layers):
    ring, layer, slot = centre
    slot %= ring_sector_count(ring)
    region = Region.around_sector(ring, layer, slot, sectors, layers, EDGE)
    listed = list(region.addresses())
    assert len(listed) == len(set(listed)) == region.count()
    assert set(listed) == _walked(region)
    assert all(region.contains(*address) for address in listed)


def test_the_centre_sector_is_always_in_and_a_far_sector_is_not():
    region = Region.around_sector(40, 1, 3, 4, 2, EDGE)
    assert region.contains(40, 1, 3)
    assert not region.contains(40, 4, 3)
    assert not region.contains(40, 1, 30)


def test_an_arc_that_crosses_slot_zero_wraps_and_splits_for_sql():
    region = Region.around_sector(40, 0, 0, 3, 0, EDGE)
    arc = next(a for a in region.ring_arcs() if a.ring == 40)
    n = ring_sector_count(40)
    assert not arc.whole and arc.first > arc.count and arc.first + arc.count > n
    low, high = region.slot_ranges(arc)
    assert low == (arc.first, n - 1) and high[0] == 0 and high[1] == arc.first + arc.count - n - 1
    assert sorted(slot for a, b in (low, high) for slot in range(a, b + 1)) == sorted(
        slot for _r, _l, slot in region.addresses() if _r == 40 and _l == 0)


def test_a_region_over_the_axis_holds_whole_rings():
    region = Region.around_sector(0, 0, 0, 3, 0, EDGE)
    arcs = region.ring_arcs()
    assert arcs[0].ring == 0 and all(a.whole for a in arcs[:3])
    assert [a.ring for a in arcs] == list(range(len(arcs)))
    # the ring just past the reach holds an arc or nothing, never a whole ring
    assert sum(a.count for a in arcs) == region.count()


def test_a_region_beside_the_axis_has_whole_inner_rings_and_arcs_outside():
    region = Region.around_sector(2, 0, 4, 6, 0, EDGE)
    arcs = region.ring_arcs()
    assert arcs[0].whole
    assert not arcs[-1].whole
    assert [a.ring for a in arcs] == sorted(a.ring for a in arcs)


def test_the_count_and_cells_follow_the_outline():
    bounds = GalaxyBounds({0: 6, 1: 3, -1: 3, 2: 1, -2: 1}, EDGE)
    region = Region.around_sector(2, 0, 3, 3, 2, EDGE)
    listed = list(region.addresses(bounds))
    assert region.count(bounds) == len(listed) < region.count()
    assert all(bounds.contains(ring, layer) for ring, layer, _slot in listed)
    # rings past a layer's edge are the only ones left out
    assert set(listed) == {a for a in region.addresses() if bounds.contains(a[0], a[1])}


def test_cells_carry_their_centre_and_distance():
    region = Region.around_sector(12, 0, 5, 2, 1, EDGE)
    centre = region.center_pc
    for ring, layer, slot, x, y, z, distance in region.cells():
        assert (x, y, z) == pytest.approx(sector_position_pc(ring, layer, slot, EDGE))
        assert distance == pytest.approx(math.dist((x, y, z), centre))


def test_spans_describe_the_same_sectors():
    region = Region.around_sector(40, 0, 0, 3, 1, EDGE)
    bounds = GalaxyBounds({layer: 100 for layer in range(-3, 4)}, EDGE)
    from_spans = set()
    for span in region.spans():
        from_spans.update(span.addresses(bounds))
    assert from_spans == set(region.addresses(bounds))
    assert sum(span.count(bounds) for span in region.spans()) == region.count(bounds)


def test_a_centre_off_a_sector_centre_uses_its_own_layer():
    position = (60.0, 25.0, 9.0)
    region = Region(position, 2 * EDGE, 1, EDGE)
    assert region.center_layer == sector_address_at(position, EDGE)[1]
    assert set(region.addresses()) == _walked(region)


def test_the_allowance_matches_the_fill_that_uses_it():
    from planetgen.generation import run_galaxy
    assert run_galaxy.cylinder_extent(1, 1, EDGE)[0] == pytest.approx((1 + CYLINDER_EDGE_ALLOWANCE) * EDGE)


@pytest.mark.parametrize("kwargs", [
    {"horizontal_pc": -1, "half_layers": 0, "edge_pc": 4}, {"horizontal_pc": 1, "half_layers": -1, "edge_pc": 4},
    {"horizontal_pc": 1, "half_layers": 0, "edge_pc": 0},
])
def test_a_region_refuses_nonsense(kwargs):
    with pytest.raises(ValueError):
        Region((0, 0, 0), **kwargs)
