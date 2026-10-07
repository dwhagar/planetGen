# tests/test_gen_grid_seams.py

"""
`docs/TODO.md` TEST.30 (grid seams and the nucleus), on the real sector
grid (`tuning.DEFAULT_SECTOR_EDGE_PC`, 4 pc):

- the angular seam: points at theta just under 2*pi (including a `y` so
  small that `atan2(y, x) % 2*pi` rounds to exactly 2*pi, the case
  `sector_address_at`'s clamp exists for) and at `-0.0`, each landing in
  a slot whose real cell holds it;
- the outermost planned ring and layer of a real plan (the Milky-Way
  defaults `planetgen plan` uses, and a small toy shape) against
  `GalaxyBounds`: the last ring and layer are inside, one more is not,
  and the faces between them fall on the documented side;
- the nucleus sector (ring 0, slot 0) on layers 0 and -1: the axis, the
  shared face between the two layers, neighbors, designations and the
  drill block holding both.

`test_fuzz_galaxy_grid.py` fuzzes the same functions over arbitrary edges;
this file pins the exact seam values on the edge production uses.
"""

import math
import random

import pytest

from planetgen.generation import run_galaxy
from planetgen.db import store
from planetgen import tuning
from planetgen.galaxy.density import build_galaxy_shape
from planetgen.galaxy.drill import parse_drill_key
from planetgen.galaxy.geometry import (
    SectorCell, galaxy_to_local_pc, layer_bounds_pc, layer_center_z_pc, neighbor_addresses,
    parse_sector_designation, provisional_sector_designation, ring_bounds_pc, ring_radius_pc,
    ring_sector_count, sector_address_at, sector_cell_vertices_pc, sector_position_pc, slot_angle_bounds,
)
from planetgen.galaxy.skeleton import (
    GalaxyBounds, bound_relative_density_at, build_layer_extents, column_extents,
    expected_system_count_at_density_1,
)

from tests.bughunt_support import mysql_argv, run_cli

EDGE = float(tuning.DEFAULT_SECTOR_EDGE_PC)
THRESHOLD = 1.0 / expected_system_count_at_density_1(tuning.DEFAULT_SECTOR_EDGE_LY)
TWO_PI = 2 * math.pi
JUST_UNDER_TWO_PI = math.nextafter(TWO_PI, 0.0)

# planetgen plan's own defaults (a real Milky-Way-scale outline, ~3,900
# rings) and test_galaxy_gen.py's small toy shape.
_MILKY_WAY = dict(disk_scale_length_pc=2800.0, disk_scale_height_pc=350.0, bulge_scale_radius_pc=200.0,
                  bulge_amplitude=1.0, arm_count=2, pitch_angle_deg=15.0, arm_amplitude=0.4)
_TOY = dict(disk_scale_length_pc=40.0, disk_scale_height_pc=12.0, bulge_scale_radius_pc=10.0,
            bulge_amplitude=2.0, arm_count=2, pitch_angle_deg=15.0, arm_amplitude=0.4)


def _shape(params):
    params = dict(params)
    return build_galaxy_shape(pitch_angle_rad=math.radians(params.pop("pitch_angle_deg")), **params)


_PLANS = {}


def _plan(name):
    if name not in _PLANS:
        shape = _shape(_MILKY_WAY if name == "milky_way" else _TOY)
        extents, outer, confirmed = build_layer_extents(shape, EDGE, THRESHOLD)
        assert confirmed and extents
        _PLANS[name] = shape, extents, outer
    return _PLANS[name]


def _cell_holds(address, point):
    """Whether `address`'s real cylindrical cell (in its own local frame)
    holds galaxy-frame `point`, and its ring/layer/angle bounds agree."""
    ring, layer, slot = address
    local = galaxy_to_local_pc(sector_position_pc(ring, layer, slot, EDGE), point)
    return SectorCell.for_ring(ring, EDGE).contains(local)


# ---------------------------------------------------------------------------
# The angular seam
# ---------------------------------------------------------------------------

SEAM_RINGS = [0, 1, 2, 7, 100, 3855, 99_999]


@pytest.mark.parametrize("ring", SEAM_RINGS)
def test_theta_just_under_two_pi_is_the_last_slot(ring):
    n = ring_sector_count(ring)
    r = ring_radius_pc(ring, EDGE)
    point = (r * math.cos(JUST_UNDER_TWO_PI), r * math.sin(JUST_UNDER_TWO_PI), 0.0)
    assert point[1] < 0
    address = sector_address_at(point, EDGE)
    assert address == (ring, 0, n - 1)
    assert _cell_holds(address, point)


@pytest.mark.parametrize("ring", SEAM_RINGS)
@pytest.mark.parametrize("y", [-1e-300, -1e-200, -1e-17])
def test_theta_that_rounds_to_two_pi_is_clamped_to_the_last_slot(ring, y):
    """`atan2(-tiny, x) % 2*pi` is exactly 2*pi in floating point, which
    would be slot `n`; it must be clamped to slot `n - 1`, never wrap
    or overflow."""
    n = ring_sector_count(ring)
    point = (ring_radius_pc(ring, EDGE), y, 0.0)
    address = sector_address_at(point, EDGE)
    assert address == (ring, 0, n - 1)
    assert _cell_holds(address, point)
    # Just the other side of the seam is slot 0.
    assert sector_address_at((point[0], -y, 0.0), EDGE) == (ring, 0, 0)


@pytest.mark.parametrize("ring", SEAM_RINGS)
def test_a_denormal_y_underflows_to_the_zero_side_of_the_seam(ring):
    """At `y = -5e-324` the angle itself underflows to -0.0, so the point
    is on the seam's slot-0 side -- still a cell that holds it."""
    point = (ring_radius_pc(ring, EDGE), -5e-324, 0.0)
    address = sector_address_at(point, EDGE)
    assert address == (ring, 0, 0)
    assert _cell_holds(address, point)


@pytest.mark.parametrize("ring", SEAM_RINGS)
def test_negative_zero_on_the_plus_x_axis_is_slot_zero_layer_zero(ring):
    r = ring_radius_pc(ring, EDGE)
    for point in ((r, -0.0, 0.0), (r, 0.0, -0.0), (r, -0.0, -0.0), (r, 0.0, 0.0)):
        address = sector_address_at(point, EDGE)
        assert address == (ring, 0, 0), point
        assert _cell_holds(address, point)


@pytest.mark.parametrize("ring", SEAM_RINGS)
def test_negative_zero_on_the_minus_x_axis_matches_positive_zero(ring):
    """theta = pi from either side of the -X axis (`y = +0.0` gives pi,
    `y = -0.0` gives -pi, the same angle mod 2*pi): one slot for both."""
    r = ring_radius_pc(ring, EDGE)
    above = sector_address_at((-r, 0.0, 0.0), EDGE)
    below = sector_address_at((-r, -0.0, 0.0), EDGE)
    assert above == below
    assert above[:2] == (ring, 0)
    assert _cell_holds(above, (-r, -0.0, 0.0))
    start, end = slot_angle_bounds(ring, above[2])
    assert start <= math.pi <= end


def test_seam_points_on_every_layer_face_keep_their_slot():
    """The seam slot is the same on the bottom face, midplane and just
    under the top face of a layer (the top face itself is the next layer)."""
    ring = 3855
    n = ring_sector_count(ring)
    r = ring_radius_pc(ring, EDGE)
    for layer in (-1, 0, 1, -317, 317):
        bottom, top = layer_bounds_pc(layer, EDGE)
        for z in (bottom, layer_center_z_pc(layer, EDGE), top - 1e-9):
            assert sector_address_at((r, -1e-300, z), EDGE) == (ring, layer, n - 1)
        assert sector_address_at((r, -1e-300, top), EDGE) == (ring, layer + 1, n - 1)


def test_one_ulp_under_a_layers_top_face_is_still_that_layer():
    """Every layer keeps a point one ulp under its top face, the plane's
    layer 0 included (GEN.31: `z / edge + 0.5` used to round up there)."""
    r = ring_radius_pc(5, EDGE)
    for layer in range(-400, 400):
        bottom, top = layer_bounds_pc(layer, EDGE)
        z = math.nextafter(top, -math.inf)
        assert sector_address_at((r, 0.0, z), EDGE)[1] == layer, layer
        assert sector_address_at((r, 0.0, bottom), EDGE)[1] == layer, layer


def test_one_ulp_under_the_planes_top_face_is_still_layer_zero():
    _bottom, top = layer_bounds_pc(0, EDGE)
    z = math.nextafter(top, -math.inf)
    assert sector_address_at((ring_radius_pc(5, EDGE), 0.0, z), EDGE)[1] == 0


def test_points_one_ulp_inside_every_face_stay_in_their_cell():
    """Rings and slots too: one ulp inside either face of a cell, and the
    face itself, give the cell whose own bounds hold the point."""
    for edge in (EDGE, 3.0, 0.7):
        for ring in (0, 1, 2, 5, 37, 400, 3855):
            inner, outer = ring_bounds_pc(ring, edge)
            for r in (inner, math.nextafter(outer, 0.0)):
                if r == 0.0:
                    continue
                assert sector_address_at((r, 0.0, 0.0), edge)[0] == ring, (edge, ring, r)
            n = ring_sector_count(ring)
            r = ring_radius_pc(ring, edge)
            for slot in (0, 1, n // 2, n - 1):
                start, end = slot_angle_bounds(ring, slot)
                for theta in (start, math.nextafter(end, 0.0)):
                    x, y = r * math.cos(theta), r * math.sin(theta)
                    got_ring, _layer, got_slot = sector_address_at((x, y, 0.0), edge)
                    low, high = slot_angle_bounds(got_ring, got_slot)
                    back = math.atan2(y, x) % (2 * math.pi)
                    assert low <= back < high or got_slot == n - 1, (edge, ring, slot, back, got_slot)


# ---------------------------------------------------------------------------
# The outermost planned ring and layer against GalaxyBounds
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("plan", ["milky_way", "toy"])
def test_each_layers_last_ring_is_inside_and_the_next_is_not(plan):
    shape, extents, outer = _plan(plan)
    bounds = GalaxyBounds(extents, EDGE)
    assert bounds.outer_ring_index == outer == dict(extents)[0]
    for layer, last in extents:
        assert bounds.contains(last, layer) and bounds.contains(0, layer)
        assert not bounds.contains(last + 1, layer)
        z = layer_center_z_pc(layer, EDGE)
        assert bound_relative_density_at(shape, ring_radius_pc(last, EDGE), z) >= THRESHOLD
        assert bound_relative_density_at(shape, ring_radius_pc(last + 1, EDGE), z) < THRESHOLD
        assert bounds.describe_miss(last + 1, layer) == \
            f"ring {last + 1} is outside the galaxy (layer {layer} ends at ring {last})"


@pytest.mark.parametrize("plan", ["milky_way", "toy"])
def test_top_and_bottom_layers_are_inside_and_the_next_are_not(plan):
    shape, extents, _outer = _plan(plan)
    bounds = GalaxyBounds(extents, EDGE)
    top = bounds.top_layer_index
    assert top == extents[0][0] == -extents[-1][0] > 0
    assert bounds.contains(0, top) and bounds.contains(0, -top)
    for layer in (top + 1, -top - 1):
        assert not bounds.contains(0, layer)
        assert bounds.describe_miss(0, layer) == \
            f"layer {layer} is outside the galaxy (layers run from {top} down to {-top})"
    assert bound_relative_density_at(shape, ring_radius_pc(0, EDGE), layer_center_z_pc(top + 1, EDGE)) < THRESHOLD


@pytest.mark.parametrize("plan", ["milky_way", "toy"])
def test_outer_face_of_the_outermost_ring_belongs_to_the_ring_outside(plan):
    """Rings are half-open `[i, i+1) * edge`: a point just inside the
    outermost ring's outer face (at the 2*pi seam) is planned, the face
    itself is the next ring and outside the galaxy."""
    _shape_, extents, outer = _plan(plan)
    bounds = GalaxyBounds(extents, EDGE)
    _inner, outer_face = ring_bounds_pc(outer, EDGE)
    n = ring_sector_count(outer)

    inside = (math.nextafter(outer_face, 0.0), -1e-300, 0.0)
    ring, layer, slot = sector_address_at(inside, EDGE)
    assert (ring, layer, slot) == (outer, 0, n - 1)
    assert bounds.contains(ring, layer)
    assert _cell_holds((ring, layer, slot), inside)

    on_face = (outer_face, 0.0, 0.0)
    ring, layer, _slot = sector_address_at(on_face, EDGE)
    assert (ring, layer) == (outer + 1, 0)
    assert not bounds.contains(ring, layer)


@pytest.mark.parametrize("plan", ["milky_way", "toy"])
def test_top_face_of_the_top_layer_is_outside_and_the_bottom_face_inside(plan):
    """Layers are half-open `[j - 1/2, j + 1/2) * edge`: the top layer's
    top face is the layer above (outside), while the bottom layer's bottom
    face is still the bottom layer (inside)."""
    _shape_, extents, _outer = _plan(plan)
    bounds = GalaxyBounds(extents, EDGE)
    top = bounds.top_layer_index
    _bottom_face, top_face = layer_bounds_pc(top, EDGE)
    lowest_face, _ = layer_bounds_pc(-top, EDGE)

    def layer_at(z):
        return sector_address_at((ring_radius_pc(0, EDGE), 0.0, z), EDGE)[1]

    assert layer_at(math.nextafter(top_face, 0.0)) == top
    assert layer_at(top_face) == top + 1 and not bounds.contains(0, top + 1)
    assert layer_at(lowest_face) == -top and bounds.contains(0, -top)
    assert layer_at(math.nextafter(lowest_face, -math.inf)) == -top - 1
    assert not bounds.contains(0, -top - 1)


@pytest.mark.parametrize("plan", ["milky_way", "toy"])
def test_random_addresses_reach_but_never_pass_the_outermost_ring(plan):
    _shape_, extents, outer = _plan(plan)
    bounds = GalaxyBounds(extents, EDGE)
    rng = random.Random(30)
    for max_ring in (outer, outer + 1, outer - 1, 0):
        for _ in range(200):
            ring, layer, slot = bounds.random_address(rng, max_ring=max_ring)
            assert bounds.contains(ring, layer)
            assert ring <= min(outer, max_ring)
            assert 0 <= slot < ring_sector_count(ring)


def test_stored_plan_bounds_match_the_computed_outline(mysql_config):
    """`planetgen plan` stores exactly the outline `build_layer_extents`
    computes, and `store.get_galaxy_bounds` hands it back on the same
    4 pc grid, outermost ring and layer included."""
    _shape_, extents, outer = _plan("toy")
    argv = []
    for key, value in _TOY.items():
        argv += [f"--{key.replace('_', '-')}", str(value)]
    run_cli("plan", argv + ["--no-bright-stars"] + mysql_argv(mysql_config))

    conn = store.get_connection(mysql_config)
    try:
        bounds = store.get_galaxy_bounds(conn)
        columns = {ring: store.get_galaxy_column(conn, ring) for ring in (0, outer, outer + 1)}
    finally:
        conn.close()
    assert bounds.edge_pc == EDGE
    assert bounds.outer_ring == dict(extents)
    assert bounds.outer_ring_index == outer
    expected_columns = {ring: (low, high) for ring, low, high in column_extents(extents)}
    assert columns[0] == expected_columns[0] == (-bounds.top_layer_index, bounds.top_layer_index)
    assert columns[outer] == expected_columns[outer]
    assert columns[outer + 1] is None


# ---------------------------------------------------------------------------
# The nucleus: ring 0, slot 0, layers 0 and -1
# ---------------------------------------------------------------------------

NUCLEUS = [(0, 0, 0), (0, -1, 0)]


def test_nucleus_sectors_sit_one_edge_apart_on_the_same_column():
    upper = sector_position_pc(0, 0, 0, EDGE)
    lower = sector_position_pc(0, -1, 0, EDGE)
    assert upper[:2] == lower[:2]
    assert upper[2] == 0.0 and lower[2] == -EDGE
    assert math.hypot(*upper[:2]) == pytest.approx(EDGE / 2)
    assert ring_sector_count(0) == 3


@pytest.mark.parametrize("z, layer", [
    (0.0, 0), (-0.0, 0), (1.999, 0), (-1.999, 0), (-2.0, 0),
    (math.nextafter(-2.0, -math.inf), -1), (-4.0, -1), (-5.999, -1), (-6.0, -1),
    (math.nextafter(-6.0, -math.inf), -2),
])
def test_axis_points_fall_in_the_nucleus_of_the_right_layer(z, layer):
    """The galactic axis on the +0.0 side is slot 0 of ring 0; the face
    between layers 0 and -1 (z = -2 pc) belongs to layer 0."""
    for x, y in ((0.0, 0.0), (0.0, -0.0)):
        address = sector_address_at((x, y, z), EDGE)
        assert address == (0, layer, 0)
        assert _cell_holds(address, (x, y, z))


@pytest.mark.parametrize("x, y", [(-0.0, 0.0), (-0.0, -0.0)])
@pytest.mark.parametrize("z", [0.0, -2.0, -4.0])
def test_axis_points_with_a_negative_zero_x_still_land_in_a_ring_zero_cell_holding_them(x, y, z):
    """`atan2(+-0.0, -0.0)` is +-pi, so these axis points get ring 0's
    slot 1, not slot 0. Every ring-0 wedge touches the axis, so that slot's
    cell still holds the point."""
    ring, layer, slot = sector_address_at((x, y, z), EDGE)
    assert (ring, layer) == (0, 0 if z > -2.5 else -1)
    assert 0 <= slot < 3
    assert _cell_holds((ring, layer, slot), (x, y, z))


@pytest.mark.parametrize("address", NUCLEUS)
def test_nucleus_cell_holds_the_axis_and_its_own_faces(address):
    ring, layer, slot = address
    bottom, top = layer_bounds_pc(layer, EDGE)
    for z in (bottom, layer_center_z_pc(layer, EDGE), top):
        assert _cell_holds(address, (0.0, 0.0, z))
    for vertex in sector_cell_vertices_pc(ring, layer, slot, EDGE):
        assert _cell_holds(address, vertex)
    assert not _cell_holds(address, (0.0, 0.0, top + 0.01))
    assert not _cell_holds(address, (0.0, 0.0, bottom - 0.01))


def test_nucleus_layers_share_one_face():
    upper = sector_cell_vertices_pc(0, 0, 0, EDGE)
    lower = sector_cell_vertices_pc(0, -1, 0, EDGE)
    # Index 4*r + 2*z + theta: layer 0's bottom corners are layer -1's top.
    for r_bit in (0, 1):
        for t_bit in (0, 1):
            assert upper[4 * r_bit + t_bit] == lower[4 * r_bit + 2 + t_bit]
    assert layer_bounds_pc(0, EDGE)[0] == layer_bounds_pc(-1, EDGE)[1] == -EDGE / 2


def test_nucleus_sectors_are_each_others_neighbors():
    upper = neighbor_addresses(0, 0, 0)
    lower = neighbor_addresses(0, -1, 0)
    assert (0, -1, 0) in upper and (0, 0, 0) in lower
    assert {(0, 0, 1), (0, 0, 2), (0, 1, 0)} <= set(upper)
    assert {(0, -1, 1), (0, -1, 2), (0, -2, 0)} <= set(lower)
    assert all(address[0] in (0, 1) for address in upper + lower)
    assert all(address[1] == -1 for address in lower if address[0] == 1)


def test_nucleus_designations_are_distinct_and_round_trip():
    codes = [provisional_sector_designation(*address) for address in NUCLEUS]
    assert codes[0] != codes[1]
    assert [parse_sector_designation(code) for code in codes] == NUCLEUS


def test_the_central_drill_block_holds_both_nucleus_layers():
    block = parse_drill_key("3.0.0.0")
    assert run_galaxy.block_layers(block) == [-1, 0, 1]
    everything = set(run_galaxy.block_addresses(block))
    assert set(NUCLEUS) <= everything
    assert (0, -1, 0) in set(run_galaxy.block_addresses(block, -1))
    assert (0, 0, 0) in set(run_galaxy.block_addresses(block, 0))


@pytest.mark.parametrize("plan", ["milky_way", "toy"])
def test_every_real_plan_holds_the_nucleus_on_both_layers(plan):
    _shape_, extents, _outer = _plan(plan)
    bounds = GalaxyBounds(extents, EDGE)
    for ring, layer, _slot in NUCLEUS:
        assert bounds.contains(ring, layer)
