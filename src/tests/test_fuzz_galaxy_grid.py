# tests/test_fuzz_galaxy_grid.py

"""
Brute-force and property tests for the galaxy's cylindrical sector grid
(`stellarObjects/galaxyGeometry.py`), its planned outline
(`galaxySkeleton.py`), the density model (`galaxyDensity.py`) and the
Galaxy Map's tile math (`galaxyViewport.py`) -- the geometry every
planning, placing and map request stands on.

Each test states one invariant and lets hypothesis hunt for the address,
point or shape that breaks it (see tests/fuzz_support.py for the
profiles). Where a closed-form shortcut exists (`enumerate_sectors_within_
radius`, `_overlapping_slots`, `GalaxyBounds.random_address`), it is
checked against a slow brute-force enumeration of the same answer.
"""

import math
import random

import pytest
from hypothesis import assume, example, given, settings
from hypothesis import strategies as st

from stellarObjects import galaxyDensity, galaxyGeometry as gg, galaxySkeleton, galaxyViewport as gv
from tests.fuzz_support import finite

EDGES = st.sampled_from([4.0, 1.0, 0.5, 3.26156, 10.0, 7.0])
RINGS = st.integers(min_value=0, max_value=5000)
LAYERS = st.integers(min_value=-4096, max_value=4095)
ANGLE_TOL = 1e-9


@st.composite
def addresses(draw, rings=RINGS, layers=LAYERS):
    ring = draw(rings)
    layer = draw(layers)
    slot = draw(st.integers(min_value=0, max_value=gg.ring_sector_count(ring) - 1))
    return ring, layer, slot


def _small_addresses():
    return addresses(rings=st.integers(0, 60), layers=st.integers(-20, 20))


# --- ring_sector_count --------------------------------------------------------

def test_ring_sector_count_first_rings_match_the_documented_sequence():
    assert [gg.ring_sector_count(i) for i in range(8)] == [3, 9, 15, 21, 27, 36, 42, 48]


@given(st.integers(min_value=0, max_value=10**7))
def test_ring_sector_count_is_a_growing_multiple_of_its_master_count(ring):
    n, n_next = gg.ring_sector_count(ring), gg.ring_sector_count(ring + 1)
    master = gg.ring_master_count(ring)
    assert n >= 3 and n % master == 0
    assert gg.ring_master_count(ring + 1) in (master, 2 * master)
    assert n_next >= n
    assert 0.94 <= 2 * math.pi * (ring + 0.5) / n <= 1.065


@given(st.integers(max_value=-1))
def test_ring_sector_count_rejects_negative_rings(ring):
    with pytest.raises(ValueError):
        gg.ring_sector_count(ring)


@given(st.integers(min_value=0, max_value=10**6), EDGES)
def test_every_ring_cell_is_close_to_a_cube(ring, edge):
    # Module docstring: the arc is within about 6% of an edge.
    cell = gg.SectorCell.for_ring(ring, edge)
    arc = cell.r_center * 2 * cell.half_angle
    if ring >= 1:
        assert 0.94 <= arc / edge <= 1.065


# --- position <-> address ------------------------------------------------------

@given(addresses(), EDGES)
def test_a_sector_center_maps_back_to_its_own_address(address, edge):
    center = gg.sector_position_pc(*address, edge)
    assert gg.sector_address_at(center, edge) == address


@given(st.tuples(finite, finite, finite), EDGES)
@example((0.0, 0.0, 0.0), 4.0)
@example((0.0, 0.0, -2.0), 4.0)  # exactly on a layer boundary
@example((4.0, 0.0, 2.0), 4.0)   # exactly on a ring boundary, +X axis
@example((-4.0, -0.0, 0.0), 4.0)  # theta = pi from below and above
@example((1e-300, -1e-300, 0.0), 4.0)
def test_any_point_lands_in_a_cell_whose_bounds_hold_it(point, edge):
    x, y, z = point
    assume(abs(x) < 1e9 and abs(y) < 1e9 and abs(z) < 1e9)
    ring, layer, slot = gg.sector_address_at(point, edge)
    n = gg.ring_sector_count(ring)
    assert 0 <= slot < n
    r_in, r_out = gg.ring_bounds_pc(ring, edge)
    r = math.hypot(x, y)
    assert r_in - 1e-9 * max(1, r) <= r < r_out + 1e-9 * max(1, r)
    z_lo, z_hi = gg.layer_bounds_pc(layer, edge)
    assert z_lo - 1e-9 * max(1, abs(z)) <= z <= z_hi + 1e-9 * max(1, abs(z))
    if r > 1e-6:
        theta = math.atan2(y, x) % (2 * math.pi)
        t_lo, t_hi = gg.slot_angle_bounds(ring, slot)
        assert t_lo - ANGLE_TOL <= theta <= t_hi + ANGLE_TOL or (slot == n - 1 and theta < ANGLE_TOL)


@given(addresses(rings=st.integers(0, 400), layers=st.integers(-50, 50)), EDGES, st.randoms(use_true_random=False))
def test_points_sampled_inside_a_cell_map_back_to_that_cell(address, edge, rng):
    ring, layer, slot = address
    center = gg.sector_position_pc(ring, layer, slot, edge)
    cell = gg.SectorCell.for_ring(ring, edge)
    for _ in range(20):
        local = cell.sample(rng)
        # Stay a hair inside so a point on a shared face can't be claimed
        # by the neighbor through rounding.
        shrunk = tuple(v * 0.999 for v in local)
        point = gg.local_to_galaxy_pc(center, (shrunk[0], shrunk[1], shrunk[2]))
        if not cell.contains(shrunk):
            continue
        assert gg.sector_address_at(point, edge) == address


@given(addresses(), EDGES)
def test_cell_vertices_are_on_the_cell_bounds(address, edge):
    ring, layer, slot = address
    vertices = gg.sector_cell_vertices_pc(ring, layer, slot, edge)
    assert len(vertices) == 8
    r_bounds = gg.ring_bounds_pc(ring, edge)
    z_bounds = gg.layer_bounds_pc(layer, edge)
    for index, (x, y, z) in enumerate(vertices):
        r_bit, z_bit = index >> 2, (index >> 1) & 1
        assert math.isclose(math.hypot(x, y), r_bounds[r_bit], rel_tol=1e-9, abs_tol=1e-9)
        assert z == z_bounds[z_bit]


@given(st.integers(0, 5000), LAYERS, st.integers())
def test_out_of_range_slots_are_rejected_everywhere(ring, layer, slot):
    n = gg.ring_sector_count(ring)
    assume(not 0 <= slot < n)
    for fn in (gg.sector_position_pc, gg.sector_cell_vertices_pc, gg.describe_sector_cell):
        with pytest.raises(ValueError):
            fn(ring, layer, slot, 4.0)
    with pytest.raises(ValueError):
        gg.neighbor_addresses(ring, layer, slot)


# --- sector-local frame --------------------------------------------------------

@given(st.tuples(finite, finite, finite).filter(lambda p: math.hypot(p[0], p[1]) < 1e8 and abs(p[2]) < 1e8))
def test_sector_orientation_is_a_right_handed_orthonormal_basis(center):
    ax, ay, az = gg.sector_orientation(center)
    dot = lambda a, b: sum(i * j for i, j in zip(a, b))
    for a in (ax, ay, az):
        assert math.isclose(dot(a, a), 1.0, rel_tol=1e-12)
    assert abs(dot(ax, ay)) < 1e-12 and abs(dot(ax, az)) < 1e-12 and abs(dot(ay, az)) < 1e-12
    cross = (ax[1] * ay[2] - ax[2] * ay[1], ax[2] * ay[0] - ax[0] * ay[2], ax[0] * ay[1] - ax[1] * ay[0])
    assert all(math.isclose(c, e, abs_tol=1e-12) for c, e in zip(cross, az))


@given(addresses(rings=st.integers(0, 2000)), EDGES, st.tuples(*[st.floats(-3, 3)] * 3))
def test_local_offsets_keep_their_length_in_the_galaxy_frame(address, edge, offset):
    center = gg.sector_position_pc(*address, edge)
    point = gg.local_to_galaxy_pc(center, offset)
    assert math.isclose(math.dist(point, center), math.hypot(*offset), rel_tol=1e-9, abs_tol=1e-9)


# --- SectorCell ----------------------------------------------------------------

@given(st.integers(0, 5000), EDGES, st.randoms(use_true_random=False))
def test_sampled_points_are_always_inside_the_cell(ring, edge, rng):
    cell = gg.SectorCell.for_ring(ring, edge)
    for _ in range(50):
        assert cell.contains(cell.sample(rng))


@given(st.integers(0, 5000), EDGES)
def test_a_rings_cells_add_up_to_the_whole_annulus(ring, edge):
    cell = gg.SectorCell.for_ring(ring, edge)
    r_in, r_out = gg.ring_bounds_pc(ring, edge)
    annulus = math.pi * (r_out ** 2 - r_in ** 2) * edge
    assert math.isclose(cell.volume * gg.ring_sector_count(ring), annulus, rel_tol=1e-12)


@given(st.integers(0, 500), EDGES)
def test_cell_contains_its_origin_and_rejects_points_just_outside(ring, edge):
    cell = gg.SectorCell.for_ring(ring, edge)
    assert cell.contains((0.0, 0.0, 0.0))
    assert not cell.contains((0.0, 0.0, cell.half_height * 1.0001 + 1e-9))
    assert not cell.contains((0.0, 0.0, -cell.half_height * 1.0001 - 1e-9))
    assert not cell.contains((cell.r_outer - cell.r_center + edge * 1e-6, 0.0, 0.0))
    if ring > 0:
        assert not cell.contains((cell.r_inner - cell.r_center - edge * 1e-6, 0.0, 0.0))
    # Just past the slot's angular edge, on the ring centerline.
    angle = cell.half_angle * 1.001 + 1e-9
    px, py = cell.r_center * math.cos(angle), cell.r_center * math.sin(angle)
    assert not cell.contains((px - cell.r_center, py, 0.0))


# --- neighbors -----------------------------------------------------------------

@given(_small_addresses())
def test_neighbor_relation_is_symmetric(address):
    for neighbor in gg.neighbor_addresses(*address):
        assert address in gg.neighbor_addresses(*neighbor), (address, neighbor)


@given(addresses())
def test_neighbors_are_unique_valid_and_exclude_the_cell(address):
    neighbors = gg.neighbor_addresses(*address)
    assert address not in neighbors
    assert len(neighbors) == len(set(neighbors))
    for ring, layer, slot in neighbors:
        assert ring >= 0
        assert 0 <= slot < gg.ring_sector_count(ring)
    ring = address[0]
    # Two in-ring, two vertical, and 1-3 overlapping slots per adjacent ring
    # (ring 0's three wedges each touch all three... no: ring 0 has no inner
    # neighbor, ring 1 has 9 slots so each ring-0 wedge overlaps 3 of them).
    assert 4 + (0 if ring == 0 else 1) + 1 <= len(neighbors) <= 4 + 3 + 3


def _float_overlaps(from_ring, slot, to_ring):
    a0, a1 = gg.slot_angle_bounds(from_ring, slot)
    result = []
    for other in range(gg.ring_sector_count(to_ring)):
        b0, b1 = gg.slot_angle_bounds(to_ring, other)
        if min(a1, b1) - max(a0, b0) > 1e-12:
            result.append(other)
    return result


@given(st.integers(0, 300), st.data())
def test_overlapping_slots_match_a_brute_force_angle_check(ring, data):
    slot = data.draw(st.integers(0, gg.ring_sector_count(ring) - 1))
    for other in (ring - 1, ring + 1):
        if other < 0:
            continue
        assert gg._overlapping_slots(ring, slot, other) == _float_overlaps(ring, slot, other)


def test_neighbors_of_every_cell_in_the_first_rings_brute_force():
    # Exhaustive, not sampled: every cell in rings 0..40 on one layer.
    for ring in range(41):
        for slot in range(gg.ring_sector_count(ring)):
            neighbors = set(gg.neighbor_addresses(ring, 0, slot))
            expected = {(ring, 0, (slot - 1) % gg.ring_sector_count(ring)),
                        (ring, 0, (slot + 1) % gg.ring_sector_count(ring)),
                        (ring, -1, slot), (ring, 1, slot)}
            for other in (ring - 1, ring + 1):
                if other >= 0:
                    expected.update((other, 0, s) for s in _float_overlaps(ring, slot, other))
            expected.discard((ring, 0, slot))
            assert neighbors == expected, (ring, slot)


# --- designations --------------------------------------------------------------

@given(addresses(rings=st.integers(0, 60000)))
def test_designations_round_trip(address):
    designation = gg.provisional_sector_designation(*address)
    assert designation == designation.upper()
    assert all(c in "0123456789ABCDEF" for c in designation)
    assert gg.parse_sector_designation(designation) == address
    assert gg.parse_sector_designation(f"  {designation.lower()}\n") == address


@given(addresses(rings=st.integers(0, 3000)), addresses(rings=st.integers(0, 3000)))
def test_distinct_addresses_get_distinct_designations(a, b):
    assume(a != b)
    assert gg.provisional_sector_designation(*a) != gg.provisional_sector_designation(*b)


@given(st.integers(0, 100), st.one_of(st.integers(max_value=-4097), st.integers(min_value=4096)))
def test_designation_rejects_layers_outside_its_bit_field(ring, layer):
    with pytest.raises(ValueError):
        gg.provisional_sector_designation(ring, layer, 0)


@given(st.integers(max_value=-1), LAYERS, st.integers(0, 2))
def test_designation_rejects_negative_rings(ring, layer, slot):
    with pytest.raises(ValueError):
        gg.provisional_sector_designation(ring, layer, slot)


@given(st.integers(0, 5000), LAYERS, st.integers(min_value=0, max_value=(1 << 20) - 1))
def test_designation_rejects_slots_the_ring_does_not_have(ring, layer, slot):
    assume(slot >= gg.ring_sector_count(ring))
    with pytest.raises(ValueError):
        gg.provisional_sector_designation(ring, layer, slot)


@given(st.text(max_size=30))
@example("")
@example("-1")
@example("G")
@example("0x")
@example("FFFFF")  # slot 1,048,575 of ring 0
def test_parse_designation_only_ever_raises_value_error(text):
    try:
        ring, layer, slot = gg.parse_sector_designation(text)
    except ValueError:
        return
    assert ring >= 0
    assert -4096 <= layer <= 4095
    assert 0 <= slot < gg.ring_sector_count(ring)


# --- zones and quadrants -------------------------------------------------------

@given(finite, finite)
def test_quadrant_is_always_one_to_four_and_matches_the_angle(x, y):
    quadrant = gg.sector_quadrant(x, y)
    assert quadrant in (1, 2, 3, 4)
    # Within a rounding error of an axis, atan2 may land on either side.
    assume(abs(x) > 1e-9 * abs(y) and abs(y) > 1e-9 * abs(x))
    if x > 0 and y > 0:
        assert quadrant == 1
    elif x < 0 < y:
        assert quadrant == 2
    elif x < 0 and y < 0:
        assert quadrant == 3
    elif y < 0 < x:
        assert quadrant == 4


@given(st.integers(0, 10**6), st.floats(0.1, 1000), st.floats(0.1, 10000))
def test_zones_are_monotone_in_ring(ring, edge_ly, target):
    assert 0 <= gg.sector_zone(ring, edge_ly, target) <= gg.sector_zone(ring + 1, edge_ly, target)


# --- enumerate_sectors_within_radius ------------------------------------------

def _brute_force_within(center, radius, edge):
    cx, cy, cz = center
    r_c = math.hypot(cx, cy)
    eps = 1e-9
    found = set()
    for ring in range(0, int((r_c + radius) / edge) + 2):
        n = gg.ring_sector_count(ring)
        for layer in range(math.floor((cz - radius) / edge) - 1, math.ceil((cz + radius) / edge) + 2):
            for slot in range(n):
                position = gg.sector_position_pc(ring, layer, slot, edge)
                if math.dist(position, center) <= radius + eps:
                    found.add((ring, layer, slot))
    return found


@settings(max_examples=40)
@given(
    st.tuples(st.floats(-60, 60), st.floats(-60, 60), st.floats(-30, 30)),
    st.floats(0, 25),
    st.sampled_from([4.0, 3.0, 7.5]),
)
@example((0.0, 0.0, 0.0), 0.0, 4.0)
@example((0.0, 0.0, 0.0), 2.0, 4.0)  # exactly reaches ring 0's centers
@example((6.0, 0.0, 0.0), 0.0, 4.0)  # exactly on a center
@example((0.0, 0.0, 2.0), 25.0, 4.0)
def test_enumerate_within_radius_matches_brute_force(center, radius, edge):
    fast = {(r, l, s) for r, l, s, *_ in gg.enumerate_sectors_within_radius(center, radius, edge)}
    assert fast == _brute_force_within(center, radius, edge)


@given(st.tuples(st.floats(-200, 200), st.floats(-200, 200), st.floats(-50, 50)), st.floats(0, 40), EDGES)
def test_enumerate_within_radius_yields_consistent_tuples(center, radius, edge):
    seen = set()
    for ring, layer, slot, x, y, z, distance in gg.enumerate_sectors_within_radius(center, radius, edge):
        assert (ring, layer, slot) not in seen
        seen.add((ring, layer, slot))
        assert (x, y, z) == pytest.approx(gg.sector_position_pc(ring, layer, slot, edge))
        assert distance <= radius + 1e-9
        assert distance == pytest.approx(math.dist((x, y, z), center))


@given(st.floats(max_value=-1e-12, allow_nan=False, allow_infinity=False))
def test_enumerate_rejects_a_negative_radius(radius):
    with pytest.raises(ValueError):
        list(gg.enumerate_sectors_within_radius((0.0, 0.0, 0.0), radius, 4.0))


@given(st.floats(max_value=0.0, allow_nan=False, allow_infinity=False))
def test_enumerate_rejects_a_non_positive_edge(edge):
    with pytest.raises(ValueError):
        list(gg.enumerate_sectors_within_radius((0.0, 0.0, 0.0), 1.0, edge))


@given(st.integers(1, 30), st.integers(0, 5000))
def test_slots_near_angle_always_returns_valid_distinct_slots(n_seed, data_seed):
    rng = random.Random(data_seed)
    n = gg.ring_sector_count(n_seed)
    theta = rng.uniform(-10, 10)
    half = rng.uniform(0, 4)
    slots = list(gg._slots_near_angle(n, theta, half))
    assert len(slots) == len(set(slots))
    assert all(0 <= s < n for s in slots)


# --- describe_sector_cell ------------------------------------------------------

@given(addresses(rings=st.integers(0, 3000)), EDGES)
def test_describe_sector_cell_is_self_consistent(address, edge):
    cell = gg.describe_sector_cell(*address, edge)
    ring, layer, slot = address
    x, y, z = cell["cartesian_pc"]
    assert gg.sector_address_at((x, y, z), edge) == address
    assert gg.parse_sector_designation(cell["designation"]) == address
    assert cell["bounds"]["r_inner_pc"] < cell["cylindrical"]["r_pc"] < cell["bounds"]["r_outer_pc"]
    assert cell["bounds"]["z_bottom_pc"] < cell["cylindrical"]["z_pc"] < cell["bounds"]["z_top_pc"]
    assert cell["volume_pc3"] > 0
    assert 0 <= cell["spherical"]["polar_rad"] <= math.pi
    for value in (x, y, z, cell["volume_pc3"], cell["mean_arc_length_pc"]):
        assert math.isfinite(value)


# --- density model -------------------------------------------------------------

SHAPES = st.builds(
    galaxyDensity.build_galaxy_shape,
    disk_scale_length_pc=st.floats(50, 10000),
    disk_scale_height_pc=st.floats(1, 2000),
    bulge_scale_radius_pc=st.floats(10, 5000),
    bulge_amplitude=st.floats(0, 100),
    arm_count=st.integers(1, 8),
    pitch_angle_rad=st.floats(0.05, 1.5),
    arm_amplitude=st.floats(0, 0.99),
)


@given(st.floats(allow_nan=False))
def test_sech_squared_is_bounded_and_never_overflows(x):
    value = galaxyDensity._sech_squared(x)
    assert 0.0 <= value <= 1.0
    assert galaxyDensity._sech_squared(-x) == value


@given(SHAPES, st.tuples(*[st.floats(-1e6, 1e6)] * 3))
def test_relative_density_is_finite_and_never_negative(shape, point):
    value = galaxyDensity.relative_density(point, shape)
    assert math.isfinite(value) and value >= 0


@given(SHAPES)
def test_shape_is_calibrated_to_one_at_its_interarm_point(shape):
    # The calibration point: 2.82 scale lengths out, between two arms.
    r = 2.82 * shape.disk_scale_length_pc
    theta_arm = shape.spiral_reference_angle_rad + math.log(r / shape.spiral_reference_radius_pc) / math.tan(shape.pitch_angle_rad)
    theta = theta_arm + math.pi / shape.arm_count
    value = galaxyDensity.relative_density((r * math.cos(theta), r * math.sin(theta), 0.0), shape)
    assert value == pytest.approx(1.0, rel=1e-9)


@given(SHAPES, st.floats(0, 1e5), st.floats(-1e5, 1e5), st.floats(0, 2 * math.pi))
def test_skeleton_bound_is_an_upper_bound_at_every_angle(shape, r_cyl, z, theta):
    point = (r_cyl * math.cos(theta), r_cyl * math.sin(theta), z)
    bound = galaxySkeleton.bound_relative_density_at(shape, r_cyl, z)
    assert galaxyDensity.relative_density(point, shape) <= bound * (1 + 1e-9) + 1e-300


# --- skeleton ------------------------------------------------------------------

@settings(max_examples=30)
@given(SHAPES, st.sampled_from([4.0, 50.0, 200.0]), st.floats(0.001, 50))
def test_layer_extents_are_symmetric_nested_and_really_qualify(shape, edge, threshold):
    extents, outer, confirmed = galaxySkeleton.build_layer_extents(shape, edge, threshold, max_ring=3000, max_layer=400)
    if not extents:
        assert outer == -1
        assert not galaxySkeleton._qualifies(shape, edge, 0, 0, threshold)
        return
    by_layer = dict(extents)
    assert by_layer[0] == outer
    for layer, ring in extents:
        assert by_layer[-layer] == ring  # symmetric about the plane
        assert 0 <= ring <= outer
        assert by_layer.get(layer - (1 if layer > 0 else -1), outer) >= ring if layer else True
        assert galaxySkeleton._qualifies(shape, edge, ring, layer, threshold)
        if confirmed:  # a capped plane walk caps every layer above it too
            assert not galaxySkeleton._qualifies(shape, edge, ring + 1, layer, threshold)
    top = max(by_layer)
    if top < 400:
        assert not galaxySkeleton._qualifies(shape, edge, 0, top + 1, threshold)


@given(st.lists(st.tuples(st.integers(-30, 30), st.integers(0, 60)), max_size=30, unique_by=lambda t: t[0]))
def test_candidate_sector_count_matches_direct_sum(extents):
    expected = sum(sum(gg.ring_sector_count(r) for r in range(outer + 1)) for _layer, outer in extents)
    assert galaxySkeleton.candidate_sector_count(extents) == expected


@given(st.lists(st.tuples(st.integers(-30, 30), st.integers(0, 60)), max_size=30, unique_by=lambda t: t[0]))
def test_column_extents_are_the_layer_extents_on_their_side(extents):
    columns = galaxySkeleton.column_extents(extents)
    for ring, low, high in columns:
        reaching = [layer for layer, outer in extents if outer >= ring]
        assert low == min(reaching) and high == max(reaching)
    all_rings = {r for _l, outer in extents for r in range(outer + 1)}
    assert [c[0] for c in columns] == sorted(all_rings)


@st.composite
def bounds_and_rng(draw):
    top = draw(st.integers(0, 6))
    plane = draw(st.integers(0, 25))
    rings = [plane]
    for _ in range(top):
        rings.append(draw(st.integers(0, rings[-1])))
    extents = [(j, rings[abs(j)]) for j in range(-top, top + 1)]
    return galaxySkeleton.GalaxyBounds(extents, 4.0), draw(st.integers(0, 2**32))


@given(bounds_and_rng(), st.one_of(st.none(), st.integers(-3, 30)))
def test_random_address_is_inside_and_respects_max_ring(pair, max_ring):
    bounds, seed = pair
    rng = random.Random(seed)
    for _ in range(40):
        address = bounds.random_address(rng, max_ring=max_ring)
        if address is None:
            assert max_ring is not None and max_ring < 0
            return
        ring, layer, slot = address
        assert bounds.contains(ring, layer)
        assert 0 <= slot < gg.ring_sector_count(ring)
        if max_ring is not None:
            assert ring <= max_ring


def test_random_address_is_uniform_over_every_cell():
    # Tiny outline, many draws: every cell hit, none more than ~2x another.
    bounds = galaxySkeleton.GalaxyBounds([(-1, 0), (0, 2), (1, 0)], 4.0)
    rng = random.Random(7)
    counts = {}
    draws = bounds.cell_count() * 400
    for _ in range(draws):
        address = bounds.random_address(rng)
        counts[address] = counts.get(address, 0) + 1
    assert len(counts) == bounds.cell_count() == 3 + 3 + (3 + 9 + 15)
    assert max(counts.values()) < 2 * min(counts.values())


@given(bounds_and_rng(), st.integers(-10, 40), st.integers(-10, 10))
def test_bounds_contains_and_describe_miss_agree(pair, ring, layer):
    bounds, _seed = pair
    if bounds.contains(ring, layer):
        assert 0 <= ring <= bounds.outer_ring[layer]
    else:
        message = bounds.describe_miss(ring, layer)
        assert str(layer) in message or str(ring) in message


def test_empty_bounds_are_falsy_and_hold_nothing():
    bounds = galaxySkeleton.GalaxyBounds([], 4.0)
    assert not bounds
    assert bounds.cell_count() == 0
    assert bounds.random_address(random.Random(1)) is None
    assert bounds.top_layer_index == -1 and bounds.outer_ring_index == -1
    assert not bounds.contains(0, 0)


# --- tiles ---------------------------------------------------------------------

TILE = st.integers(0, gv.TILE_MAX_LEVEL).flatmap(
    lambda level: st.tuples(st.just(level), *[st.integers(0, 2 ** level - 1)] * 3)
)


@given(TILE)
def test_tile_keys_round_trip(tile):
    assert gv.parse_tile_key(gv.tile_key(*tile)) == tile


@given(st.text(max_size=40))
@example("")
@example("0/0/0")
@example("0/0/0/0/0")
@example("13/0/0/0")
@example("-1/0/0/0")
@example("1/2/0/0")
@example("1/-1/0/0")
@example("1e3/0/0/0")
@example("99999999999999999999/0/0/0")
def test_parse_tile_key_only_ever_raises_value_error(text):
    try:
        level, ix, iy, iz = gv.parse_tile_key(text)
    except ValueError:
        return
    assert 0 <= level <= gv.TILE_MAX_LEVEL
    assert all(0 <= i < 2 ** level for i in (ix, iy, iz))


@given(TILE)
def test_tile_boxes_nest_inside_their_parent(tile):
    level, ix, iy, iz = tile
    lo, hi = gv.tile_bounds_pc(level, ix, iy, iz)
    assert all(h - l == pytest.approx(gv.tile_edge_pc(level)) for l, h in zip(lo, hi))
    if level:
        plo, phi = gv.tile_bounds_pc(level - 1, ix // 2, iy // 2, iz // 2)
        assert all(pl <= l and h <= ph for pl, l, h, ph in zip(plo, lo, hi, phi))


HALF_ROOT = gv.TILE_ROOT_EDGE_PC / 2


@given(st.tuples(*[st.floats(-HALF_ROOT * 1.1, HALF_ROOT * 1.1)] * 3))
@example((-HALF_ROOT, -HALF_ROOT, -HALF_ROOT))  # the root cube's closed corner
@example((HALF_ROOT, 0.0, 0.0))  # the open face: outside
@example((0.0, 0.0, 0.0))
def test_tile_keys_containing_really_contain_the_point(point):
    keys = gv.tile_keys_containing(point)
    inside_root = all(-HALF_ROOT <= v < HALF_ROOT for v in point)
    if not inside_root:
        assert keys == []
        return
    assert len(keys) == gv.TILE_MAX_LEVEL + 1
    for level, key in enumerate(keys):
        parsed = gv.parse_tile_key(key)
        assert parsed[0] == level
        lo, hi = gv.tile_bounds_pc(*parsed)
        assert all(l <= v < h for l, v, h in zip(lo, point, hi))


@given(st.floats(allow_nan=False))
def test_tile_level_for_view_radius_is_always_a_valid_level(radius):
    level = gv.tile_level_for_view_radius(radius)
    assert 0 <= level <= gv.TILE_MAX_LEVEL
    if 0 < radius <= gv.TILE_ROOT_EDGE_PC and level < gv.TILE_MAX_LEVEL:
        assert gv.tile_edge_pc(level) >= radius
        assert gv.tile_edge_pc(level + 1) < radius


@settings(max_examples=40)
@given(
    st.integers(0, 5),
    st.tuples(*[st.floats(-HALF_ROOT * 1.2, HALF_ROOT * 1.2)] * 3),
    st.floats(0, HALF_ROOT),
)
def test_tiles_intersecting_sphere_match_brute_force(level, center, radius):
    span = 2 ** level
    expected = set()
    for ix in range(span):
        for iy in range(span):
            for iz in range(span):
                lo, hi = gv.tile_bounds_pc(level, ix, iy, iz)
                d2 = sum((c - min(max(c, l), h)) ** 2 for c, l, h in zip(center, lo, hi))
                if d2 <= radius * radius:
                    expected.add(gv.tile_key(level, ix, iy, iz))
    got = gv.tiles_intersecting_sphere(level, center, radius)
    assert len(got) == len(set(got))
    assert set(got) <= expected
    # Boxes are half-open, so a tile the sphere only touches along its
    # open upper face is rightly left out; nothing else may be.
    for key in expected - set(got):
        lo, hi = gv.tile_bounds_pc(*gv.parse_tile_key(key))
        d2 = sum((c - min(max(c, l), h)) ** 2 for c, l, h in zip(center, lo, hi))
        assert d2 == pytest.approx(radius * radius, abs=1e-6)


@given(TILE.filter(lambda t: t[0] >= 6), st.sampled_from([4.0]))
def test_planned_slots_in_a_tile_lie_in_that_tile(tile, edge):
    shape = galaxyDensity.build_galaxy_shape(3000, 300, 500, 5, 2, 0.2, 0.5)
    e = galaxySkeleton.expected_system_count_at_density_1()
    slots = gv.planned_slots_in_tile(*tile, edge, shape, e, set())
    lo, hi = gv.tile_bounds_pc(*tile)
    for slot in slots:
        position = gg.sector_position_pc(slot["ring_index"], slot["layer_index"], slot["ring_slot_index"], edge)
        assert all(l <= v < h for l, v, h in zip(lo, position, hi))
