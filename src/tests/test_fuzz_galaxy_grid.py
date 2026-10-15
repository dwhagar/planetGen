# tests/test_fuzz_galaxy_grid.py

"""
Brute-force and property tests for the galaxy's spherical sector grid
(`planetgen/galaxy/spherical.py`: shells, bands and wedges), its planned
outline (`skeleton.py`), the density model (`density.py`) and the Galaxy
Map's tile math (`viewport.py`) -- the geometry every planning, placing and
map request stands on.

Each test states one invariant and lets hypothesis hunt for the address,
point or shape that breaks it (see tests/fuzz_support.py for the
profiles). Where a closed-form shortcut exists (`enumerate_sectors_within_
radius`, `_overlapping_wedges`, `neighbor_addresses`,
`GalaxyBounds.random_address`), it is checked against a slow brute-force
enumeration of the same answer.
"""

import math
import random

import pytest
from hypothesis import assume, example, given, settings
from hypothesis import strategies as st

from planetgen import tuning
from planetgen.galaxy import (
    density as galaxyDensity, geometry as gg, skeleton as galaxySkeleton, spherical, viewport as gv,
)
from tests.fuzz_support import finite

EDGES = st.sampled_from([4.0, 1.0, 0.5, 3.26156, 10.0, 7.0])
SHELLS = st.integers(min_value=0, max_value=4095)
ANGLE_TOL = 1e-9


@st.composite
def addresses(draw, shells=SHELLS):
    """A cell of the grid: a shell, one of its bands (caps included) and a wedge of that band."""
    shell = draw(shells)
    low, high = spherical.band_range(shell)
    band = draw(st.integers(min_value=low, max_value=high))
    wedge = draw(st.integers(min_value=0, max_value=spherical.wedge_count(shell, band) - 1))
    return shell, band, wedge


def _small_addresses():
    return addresses(shells=st.integers(0, 40))


def _is_smooth(m):
    for prime in (2, 3, 5, 7):
        while m % prime == 0:
            m //= prime
    return m == 1


# --- the layout of a shell ------------------------------------------------------

def test_the_first_shells_match_the_documented_layout():
    assert spherical.cell_count(0) == 1 and spherical.band_range(0) == (0, 0)
    # Shell 1 is hand cut: one band of 12 wedges between two caps of 37.5 degrees.
    assert spherical.band_range(1) == (-1, 1)
    assert [spherical.wedge_count(1, band) for band in (-1, 0, 1)] == [1, 12, 1]
    assert math.degrees(spherical.band_colatitude_bounds(1, 1)[1]) == pytest.approx(37.5)


@given(st.integers(min_value=1, max_value=10**6))
def test_a_shell_has_an_odd_number_of_bands_between_two_caps(shell):
    layout = spherical.shell_layout(shell)
    low, high = spherical.band_range(shell)
    assert layout.band_count % 2 == 1 and (low, high) == (-(layout.high_band + 1), layout.high_band + 1)
    assert spherical.is_cap(shell, high) and spherical.is_cap(shell, low)
    assert not spherical.is_cap(shell, 0)
    assert 0 < layout.cap_colatitude < math.pi / 2


@given(st.integers(min_value=2, max_value=3000), st.data())
def test_ordinary_bands_hold_six_times_a_smooth_number_of_wedges(shell, data):
    low, high = spherical.band_range(shell)
    band = data.draw(st.integers(low + 1, high - 1))
    n = spherical.wedge_count(shell, band)
    assert n % spherical.CAP_RIM_WEDGES == 0 and _is_smooth(n // spherical.CAP_RIM_WEDGES)
    # North and south halves are mirror images.
    assert spherical.wedge_count(shell, -band) == n
    # A cap is one cell.
    assert spherical.wedge_count(shell, high) == spherical.wedge_count(shell, low) == 1


@given(st.integers(min_value=8, max_value=4000))
def test_the_equatorial_band_has_about_one_edge_of_arc_per_wedge(shell):
    n = spherical.wedge_count(shell, 0)
    arc = 2 * math.pi * shell / n
    assert 0.9 <= arc <= 1.1


@given(st.integers(min_value=-5, max_value=-1))
def test_a_negative_shell_is_rejected(shell):
    with pytest.raises(ValueError):
        spherical.shell_layout(shell)


@given(st.integers(0, 3000))
def test_a_shells_cell_count_is_the_sum_of_its_bands_wedges(shell):
    low, high = spherical.band_range(shell)
    assert spherical.cell_count(shell) == sum(spherical.wedge_count(shell, band) for band in range(low, high + 1))


@given(st.integers(0, 300), st.data())
def test_nth_cell_in_band_range_counts_cells_north_to_south(shell, data):
    high = data.draw(st.integers(0, spherical.band_range(shell)[1]))
    total = spherical.cells_in_band_range(shell, high)
    index = data.draw(st.integers(0, total - 1))
    s, band, wedge = spherical.nth_cell_in_band_range(shell, high, index)
    assert s == shell and -high <= band <= high and 0 <= wedge < spherical.wedge_count(shell, band)
    with pytest.raises(ValueError):
        spherical.nth_cell_in_band_range(shell, high, total)


# --- position <-> address ------------------------------------------------------

@given(addresses(), EDGES)
def test_a_sector_center_maps_back_to_its_own_address(address, edge):
    center = spherical.sector_position_pc(*address, edge)
    assert spherical.sector_address_at(center, edge) == address


@given(st.tuples(finite, finite, finite), EDGES)
@example((0.0, 0.0, 0.0), 4.0)
@example((0.0, 0.0, 2.0), 4.0)    # exactly on the core's surface, on the axis
@example((0.0, 0.0, -2.0), 4.0)
@example((4.0, 0.0, 0.0), 4.0)    # in the middle of shell 1, +X axis
@example((-4.0, -0.0, 0.0), 4.0)  # longitude pi from below and above
@example((1e-300, -1e-300, 0.0), 4.0)
def test_any_point_lands_in_a_cell_whose_bounds_hold_it(point, edge):
    x, y, z = point
    # Past about 1e7 edges the wedge multiplier leaves its table (`spherical.SMOOTH_LIMIT`), far beyond the 262 kpc reach.
    assume(math.sqrt(x * x + y * y + z * z) < 1e6 * edge)
    shell, band, wedge = spherical.sector_address_at(point, edge)
    assert 0 <= wedge < spherical.wedge_count(shell, band)
    r = math.sqrt(x * x + y * y + z * z)
    r_in, r_out = spherical.shell_bounds_pc(shell, edge)
    assert r_in - 1e-9 * max(1, r) <= r < r_out + 1e-9 * max(1, r)
    north, south = spherical.band_colatitude_bounds(shell, band)
    colat = math.atan2(math.hypot(x, y), z)
    assert north - ANGLE_TOL <= colat <= south + ANGLE_TOL
    if math.hypot(x, y) > 1e-6 and spherical.wedge_count(shell, band) > 1:
        theta = math.atan2(y, x) % (2 * math.pi)
        start, end = spherical.wedge_longitude_bounds(shell, band, wedge)
        n = spherical.wedge_count(shell, band)
        assert start - ANGLE_TOL <= theta <= end + ANGLE_TOL or (wedge == n - 1 and theta < ANGLE_TOL)


@given(addresses(shells=st.integers(0, 400)), EDGES, st.randoms(use_true_random=False))
def test_points_sampled_inside_a_cell_map_back_to_that_cell(address, edge, rng):
    center = spherical.sector_position_pc(*address, edge)
    cell = spherical.SectorCell.for_address(*address, edge)
    for _ in range(20):
        local = cell.sample(rng)
        # Stay a hair inside so a point on a shared face can't be claimed
        # by the neighbor through rounding.
        shrunk = tuple(v * 0.999 for v in local)
        point = spherical.local_to_galaxy_pc(center, shrunk)
        if not cell.contains(shrunk):
            continue
        assert spherical.sector_address_at(point, edge) == address


@given(addresses(), EDGES)
def test_cell_vertices_are_on_the_cell_bounds(address, edge):
    shell, band, wedge = address
    vertices = spherical.sector_cell_vertices_pc(shell, band, wedge, edge)
    assert len(vertices) == 8
    r_bounds = spherical.shell_bounds_pc(shell, edge)
    colat_bounds = spherical.band_colatitude_bounds(shell, band)
    for index, (x, y, z) in enumerate(vertices):
        r_bit, colat_bit = index >> 2, (index >> 1) & 1
        assert math.isclose(math.sqrt(x * x + y * y + z * z), r_bounds[r_bit], rel_tol=1e-9, abs_tol=1e-9)
        assert math.isclose(math.atan2(math.hypot(x, y), z), colat_bounds[colat_bit], abs_tol=1e-9) or r_bounds[r_bit] == 0


@given(st.integers(0, 3000), st.integers(), st.integers())
def test_out_of_range_bands_and_wedges_are_rejected_everywhere(shell, band, wedge):
    low, high = spherical.band_range(shell)
    assume(not (low <= band <= high and 0 <= wedge < spherical.wedge_count(shell, band)) if low <= band <= high
           else True)
    for fn in (spherical.sector_position_pc, spherical.sector_cell_vertices_pc, spherical.describe_sector_cell):
        with pytest.raises(ValueError):
            fn(shell, band, wedge, 4.0)
    with pytest.raises(ValueError):
        spherical.neighbor_addresses(shell, band, wedge)


# --- sector-local frame --------------------------------------------------------

@given(st.tuples(finite, finite, finite).filter(lambda p: math.hypot(p[0], p[1]) < 1e8 and abs(p[2]) < 1e8))
def test_sector_orientation_is_a_right_handed_orthonormal_basis(center):
    ax, ay, az = spherical.sector_orientation(center)
    dot = lambda a, b: sum(i * j for i, j in zip(a, b))
    for a in (ax, ay, az):
        assert math.isclose(dot(a, a), 1.0, rel_tol=1e-12)
    assert abs(dot(ax, ay)) < 1e-12 and abs(dot(ax, az)) < 1e-12 and abs(dot(ay, az)) < 1e-12
    cross = (ax[1] * ay[2] - ax[2] * ay[1], ax[2] * ay[0] - ax[0] * ay[2], ax[0] * ay[1] - ax[1] * ay[0])
    assert all(math.isclose(c, e, abs_tol=1e-12) for c, e in zip(cross, az))


@given(addresses(shells=st.integers(0, 2000)), EDGES, st.tuples(*[st.floats(-3, 3)] * 3))
def test_local_offsets_keep_their_length_in_the_galaxy_frame(address, edge, offset):
    center = spherical.sector_position_pc(*address, edge)
    point = spherical.local_to_galaxy_pc(center, offset)
    assert math.isclose(math.dist(point, center), math.hypot(*offset), rel_tol=1e-9, abs_tol=1e-9)


# --- SectorCell ----------------------------------------------------------------

@given(addresses(shells=st.integers(0, 3000)), EDGES, st.randoms(use_true_random=False))
def test_sampled_points_are_always_inside_the_cell(address, edge, rng):
    cell = spherical.SectorCell.for_address(*address, edge)
    for _ in range(50):
        assert cell.contains(cell.sample(rng))


@given(st.integers(0, 40), EDGES)
def test_a_shells_cells_add_up_to_the_whole_shell(shell, edge):
    inner, outer = spherical.shell_bounds_pc(shell, edge)
    total = sum(spherical.cell_volume(*address, edge) for address in spherical.shell_cells(shell))
    assert math.isclose(total, 4 / 3 * math.pi * (outer ** 3 - inner ** 3), rel_tol=1e-9)


@given(addresses(shells=st.integers(1, 500)), EDGES)
def test_a_cells_volume_matches_its_cell_object(address, edge):
    cell = spherical.SectorCell.for_address(*address, edge)
    assert math.isclose(cell.volume, spherical.cell_volume(*address, edge), rel_tol=1e-9)


@given(addresses(shells=st.integers(1, 500)), EDGES)
def test_cell_contains_its_center_and_rejects_points_just_outside(address, edge):
    shell, band, wedge = address
    cell = spherical.SectorCell.for_address(shell, band, wedge, edge)
    assert cell.contains((0.0, 0.0, 0.0))

    def local(radius, colat, lon):
        s = radius * math.sin(colat)
        return spherical.galaxy_to_local_pc(cell.center, (s * math.cos(lon), s * math.sin(lon), radius * math.cos(colat)))

    colat_mid = (cell.colat_north + cell.colat_south) / 2
    lon_mid = (cell.lon_start + cell.lon_end) / 2
    assert cell.contains(local(cell.r_center, colat_mid, lon_mid))
    # Radially past either face.
    assert not cell.contains(local(cell.r_outer + edge * 1e-6, colat_mid, lon_mid))
    assert not cell.contains(local(cell.r_inner - edge * 1e-6, colat_mid, lon_mid))
    # Past a band's latitude face.
    if cell.colat_north > 1e-3:
        assert not cell.contains(local(cell.r_center, cell.colat_north - 1e-6, lon_mid))
    if cell.colat_south < math.pi - 1e-3:
        assert not cell.contains(local(cell.r_center, cell.colat_south + 1e-6, lon_mid))
    # Past a wedge's longitude face.
    if cell.lon_end - cell.lon_start < 2 * math.pi - 1e-9:
        assert not cell.contains(local(cell.r_center, colat_mid, cell.lon_end + 1e-6))
        assert not cell.contains(local(cell.r_center, colat_mid, cell.lon_start - 1e-6))


# --- neighbors -----------------------------------------------------------------

@given(_small_addresses())
def test_neighbor_relation_is_symmetric(address):
    for neighbor in spherical.neighbor_addresses(*address):
        assert address in spherical.neighbor_addresses(*neighbor), (address, neighbor)


@given(addresses())
def test_neighbors_are_unique_valid_and_exclude_the_cell(address):
    neighbors = spherical.neighbor_addresses(*address)
    assert address not in neighbors
    assert len(neighbors) == len(set(neighbors))
    for shell, band, wedge in neighbors:
        assert shell >= 0 and abs(shell - address[0]) <= 1
        low, high = spherical.band_range(shell)
        assert low <= band <= high and 0 <= wedge < spherical.wedge_count(shell, band)
    assert neighbors   # every cell touches something


def _float_overlaps(from_count, wedge, to_count):
    a0, a1 = 2 * math.pi * wedge / from_count, 2 * math.pi * (wedge + 1) / from_count
    result = []
    for other in range(to_count):
        b0, b1 = 2 * math.pi * other / to_count, 2 * math.pi * (other + 1) / to_count
        if min(a1, b1) - max(a0, b0) > 1e-12:
            result.append(other)
    return result


@given(st.integers(1, 400), st.integers(1, 400), st.data())
def test_overlapping_wedges_match_a_brute_force_angle_check(from_count, to_count, data):
    wedge = data.draw(st.integers(0, from_count - 1))
    assert list(spherical._overlapping_wedges(from_count, wedge, to_count)) == _float_overlaps(from_count, wedge, to_count)


def _brute_force_neighbors(shell, band, wedge):
    """Every cell sharing a face of nonzero area, by testing every band of the three nearest shells."""
    n = spherical.wedge_count(shell, band)
    north, south = spherical.band_colatitude_bounds(shell, band)
    found = set()
    for other_shell in (shell - 1, shell, shell + 1):
        if other_shell < 0:
            continue
        low, high = spherical.band_range(other_shell)
        for other in range(low, high + 1):
            o_north, o_south = spherical.band_colatitude_bounds(other_shell, other)
            m = spherical.wedge_count(other_shell, other)
            if other_shell == shell:
                if other == band:
                    if n > 1:
                        found.update({(shell, band, (wedge - 1) % n), (shell, band, (wedge + 1) % n)})
                    continue
                touches = abs(other - band) == 1 and min(south, o_south) - max(north, o_north) > -1e-12
                if not touches:
                    continue
            elif min(south, o_south) - max(north, o_north) <= 1e-12:
                continue
            found.update((other_shell, other, w) for w in _float_overlaps(n, wedge, m))
    found.discard((shell, band, wedge))
    return found


@settings(max_examples=60)
@given(addresses(shells=st.integers(0, 45)))
def test_neighbors_match_a_brute_force_scan_of_the_nearby_bands(address):
    assert set(spherical.neighbor_addresses(*address)) == _brute_force_neighbors(*address)


def test_the_core_touches_all_of_shell_one_and_nothing_else():
    assert set(spherical.neighbor_addresses(0, 0, 0)) == set(spherical.shell_cells(1))


# --- designations --------------------------------------------------------------

def _designatable(address):
    return -spherical.DESIGNATION_BAND_BIAS <= address[1] < spherical.DESIGNATION_BAND_BIAS


@given(addresses().filter(_designatable))
def test_designations_round_trip(address):
    designation = spherical.provisional_sector_designation(*address)
    assert designation == designation.upper()
    assert all(c in "0123456789ABCDEF" for c in designation)
    assert spherical.parse_sector_designation(designation) == address
    assert spherical.parse_sector_designation(f"  {designation.lower()}\n") == address


@given(addresses(shells=st.integers(0, 3000)).filter(_designatable),
       addresses(shells=st.integers(0, 3000)).filter(_designatable))
def test_distinct_addresses_get_distinct_designations(a, b):
    assume(a != b)
    assert spherical.provisional_sector_designation(*a) != spherical.provisional_sector_designation(*b)


@given(addresses(shells=st.integers(1300, 4095)).filter(lambda address: not _designatable(address)))
def test_designation_rejects_bands_outside_its_bit_field(address):
    # The caps and high bands of the outer shells (past about shell 1,300) are beyond the 12-bit band field.
    with pytest.raises(ValueError):
        spherical.provisional_sector_designation(*address)


@given(st.integers(min_value=1 << 12, max_value=10**6))
def test_designation_rejects_shells_outside_its_bit_field(shell):
    with pytest.raises(ValueError):
        spherical.provisional_sector_designation(shell, 0, 0)


@given(st.integers(max_value=-1), st.integers(-3, 3), st.integers(0, 2))
def test_designation_rejects_negative_shells(shell, band, wedge):
    with pytest.raises(ValueError):
        spherical.provisional_sector_designation(shell, band, wedge)


@given(st.integers(0, 3000), st.integers(min_value=0, max_value=(1 << 20) - 1))
def test_designation_rejects_wedges_the_band_does_not_have(shell, wedge):
    assume(wedge >= spherical.wedge_count(shell, 0))
    with pytest.raises(ValueError):
        spherical.provisional_sector_designation(shell, 0, wedge)


@given(st.text(max_size=30))
@example("")
@example("-1")
@example("G")
@example("0x")
@example("FFFFF")  # shell 0, a band and a wedge the core does not have
def test_parse_designation_only_ever_raises_value_error(text):
    try:
        shell, band, wedge = spherical.parse_sector_designation(text)
    except ValueError:
        return
    assert shell >= 0
    low, high = spherical.band_range(shell)
    assert low <= band <= high
    assert 0 <= wedge < spherical.wedge_count(shell, band)


# --- zones and quadrants -------------------------------------------------------

@given(finite, finite)
def test_quadrant_is_always_one_to_four_and_matches_the_angle(x, y):
    quadrant = spherical.sector_quadrant(x, y)
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
def test_zones_are_monotone_in_shell(shell, edge_ly, target):
    assert 0 <= spherical.sector_zone(shell, edge_ly, target) <= spherical.sector_zone(shell + 1, edge_ly, target)


# --- enumerate_sectors_within_radius ------------------------------------------

def _brute_force_within(center, radius, edge):
    reach = math.dist(center, (0.0, 0.0, 0.0))
    eps = 1e-9
    found = set()
    # A cell's center is `shell * edge` from the origin, so only these shells can hold one in range.
    first = max(0, math.floor((reach - radius) / edge) - 1)
    for shell in range(first, int((reach + radius) / edge) + 2):
        for shell_, band, wedge in spherical.shell_cells(shell):
            position = spherical.sector_position_pc(shell, band, wedge, edge)
            if math.dist(position, center) <= radius + eps:
                found.add((shell, band, wedge))
    return found


@settings(max_examples=40)
@given(
    st.tuples(st.floats(-40, 40), st.floats(-40, 40), st.floats(-30, 30)),
    st.floats(0, 25),
    st.sampled_from([4.0, 3.0, 7.5]),
)
@example((0.0, 0.0, 0.0), 0.0, 4.0)
@example((0.0, 0.0, 0.0), 4.0, 4.0)   # exactly reaches shell 1's centers
@example((8.0, 0.0, 0.0), 0.0, 4.0)   # exactly on a center
@example((0.0, 0.0, 4.0), 25.0, 4.0)  # on the axis, on a cap's center
def test_enumerate_within_radius_matches_brute_force(center, radius, edge):
    fast = {(s, b, w) for s, b, w, *_ in spherical.enumerate_sectors_within_radius(center, radius, edge)}
    assert fast == _brute_force_within(center, radius, edge)


@given(st.tuples(st.floats(-200, 200), st.floats(-200, 200), st.floats(-50, 50)), st.floats(0, 40), EDGES)
def test_enumerate_within_radius_yields_consistent_tuples(center, radius, edge):
    seen = set()
    for shell, band, wedge, x, y, z, distance in spherical.enumerate_sectors_within_radius(center, radius, edge):
        assert (shell, band, wedge) not in seen
        seen.add((shell, band, wedge))
        # math.isclose with pytest.approx's default tolerances: this loop
        # can visit a million sectors, and approx() costs ~10x more each.
        expected = spherical.sector_position_pc(shell, band, wedge, edge)
        assert all(_approx(a, b) for a, b in zip((x, y, z), expected)), ((x, y, z), expected)
        assert distance <= radius + 1e-9
        assert _approx(distance, math.dist((x, y, z), center))


def _approx(a, b):
    return math.isclose(a, b, rel_tol=1e-6, abs_tol=1e-12)


@given(st.floats(max_value=-1e-12, allow_nan=False, allow_infinity=False))
def test_enumerate_rejects_a_negative_radius(radius):
    with pytest.raises(ValueError):
        list(spherical.enumerate_sectors_within_radius((0.0, 0.0, 0.0), radius, 4.0))


@given(st.floats(max_value=0.0, allow_nan=False, allow_infinity=False))
def test_enumerate_rejects_a_non_positive_edge(edge):
    with pytest.raises(ValueError):
        list(spherical.enumerate_sectors_within_radius((0.0, 0.0, 0.0), 1.0, edge))


@given(st.integers(1, 400), st.integers(0, 5000))
def test_wedges_near_longitude_returns_valid_distinct_wedges_and_every_close_one(count, data_seed):
    rng = random.Random(data_seed)
    center = rng.uniform(0, 2 * math.pi)
    half = rng.uniform(0, 4)
    wedges = list(spherical._wedges_near_longitude(count, center, half))
    assert len(wedges) == len(set(wedges)) and all(0 <= w < count for w in wedges)
    step = 2 * math.pi / count
    for w in range(count):
        middle = (w + 0.5) * step
        gap = abs((middle - center + math.pi) % (2 * math.pi) - math.pi)
        if gap < half - 1e-6:
            assert w in wedges, (count, center, half, w)


# --- describe_sector_cell ------------------------------------------------------

@given(addresses(shells=st.integers(0, 3000)).filter(_designatable), EDGES)
def test_describe_sector_cell_is_self_consistent(address, edge):
    cell = spherical.describe_sector_cell(*address, edge)
    shell, band, wedge = address
    x, y, z = cell["cartesian_pc"]
    assert spherical.sector_address_at((x, y, z), edge) == address
    assert spherical.parse_sector_designation(cell["designation"]) == address
    assert (cell["shell_index"], cell["band_index"], cell["wedge_index"]) == address
    bounds = cell["bounds"]
    assert bounds["r_inner_pc"] <= cell["spherical"]["r_pc"] < bounds["r_outer_pc"]
    assert bounds["colatitude_north_rad"] <= cell["spherical"]["polar_rad"] <= bounds["colatitude_south_rad"]
    assert math.isclose(math.hypot(x, y), cell["cylindrical"]["r_pc"], abs_tol=1e-9 * max(1.0, abs(z)))
    assert cell["volume_pc3"] > 0
    assert 0 <= cell["spherical"]["polar_rad"] <= math.pi
    assert len(cell["vertices_pc"]) == 8
    for value in (x, y, z, cell["volume_pc3"], cell["mean_band_arc_pc"], cell["mean_wedge_arc_pc"]):
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
    # The calibration point: the solar radius (3.15 scale lengths) out, between two arms.
    r = tuning.GALAXY_SOLAR_RADIUS_TO_SCALE_LENGTH * shape.disk_scale_length_pc
    theta_arm = shape.spiral_reference_angle_rad + math.log(r / shape.spiral_reference_radius_pc) / math.tan(shape.pitch_angle_rad)
    theta = theta_arm + math.pi / shape.arm_count
    value = galaxyDensity.relative_density((r * math.cos(theta), r * math.sin(theta), 0.0), shape)
    assert value == pytest.approx(1.0, rel=1e-9)


@given(SHAPES, st.floats(0, 1e5), st.floats(-1e5, 1e5), st.floats(0, 2 * math.pi))
def test_skeleton_bound_is_an_upper_bound_at_every_angle(shape, r_cyl, z, theta):
    point = (r_cyl * math.cos(theta), r_cyl * math.sin(theta), z)
    bound = galaxySkeleton.bound_relative_density_at(shape, r_cyl, z)
    # The model's own density: the halo floor (GEN.78) is left out of the
    # skeleton's outline on purpose, so it can't widen it.
    model = shape.k_norm * galaxyDensity._raw_density(point, shape)
    assert model <= bound * (1 + 1e-9) + 1e-300
    assert galaxyDensity.relative_density(point, shape) == max(model, tuning.MIN_RELATIVE_DENSITY)


# --- skeleton ------------------------------------------------------------------

@settings(max_examples=30)
@given(SHAPES, st.sampled_from([4.0, 50.0, 200.0]), st.floats(0.001, 50))
def test_shell_extents_run_from_the_core_and_really_qualify(shape, edge, threshold):
    extents, outer, confirmed = galaxySkeleton.build_shell_extents(shape, edge, threshold, max_shell=300)
    if not extents:
        assert outer == -1
        assert not galaxySkeleton._band_qualifies(shape, edge, 0, 0, threshold)
        return
    assert [shell for shell, _high in extents] == list(range(outer + 1))
    for shell, high in extents:
        top = spherical.band_range(shell)[1]
        assert 0 <= high <= top
        assert galaxySkeleton._band_qualifies(shape, edge, shell, 0, threshold)
        assert galaxySkeleton._band_qualifies(shape, edge, shell, high, threshold)
        if high < top:  # the last qualifying band, not just a qualifying one
            assert not galaxySkeleton._band_qualifies(shape, edge, shell, high + 1, threshold)
    if confirmed:
        assert not galaxySkeleton._band_qualifies(shape, edge, outer + 1, 0, threshold)
    else:
        assert outer == 300


@st.composite
def outlines(draw, max_outer=60):
    """`(shell, high_band)` rows for shells `0 .. outer`, each reaching any band the shell has."""
    outer = draw(st.integers(0, max_outer))
    return [(shell, draw(st.integers(0, spherical.band_range(shell)[1]))) for shell in range(outer + 1)]


@given(outlines())
def test_candidate_sector_count_matches_direct_sum(extents):
    expected = sum(sum(spherical.wedge_count(shell, band) for band in range(-high, high + 1))
                   for shell, high in extents)
    assert galaxySkeleton.shell_candidate_sector_count(extents) == expected
    assert [count for _s, _h, count in galaxySkeleton.shell_candidate_counts(extents)] == [
        spherical.cells_in_band_range(shell, high) for shell, high in extents]


@st.composite
def bounds_and_rng(draw):
    return galaxySkeleton.GalaxyBounds(draw(outlines(25)), 4.0), draw(st.integers(0, 2**32))


@given(bounds_and_rng(), st.one_of(st.none(), st.integers(-3, 30)))
def test_random_address_is_inside_and_respects_max_shell(pair, max_shell):
    bounds, seed = pair
    rng = random.Random(seed)
    for _ in range(40):
        address = bounds.random_address(rng, max_shell=max_shell)
        if address is None:
            assert max_shell is not None and max_shell < 0
            return
        shell, band, wedge = address
        assert bounds.contains(shell, band)
        assert 0 <= wedge < spherical.wedge_count(shell, band)
        if max_shell is not None:
            assert shell <= max_shell


def test_random_address_is_uniform_over_every_cell():
    # Tiny outline, many draws: every cell hit, none more than ~2x another.
    bounds = galaxySkeleton.GalaxyBounds([(0, 0), (1, 1), (2, 0)], 4.0)
    rng = random.Random(7)
    counts = {}
    draws = bounds.cell_count() * 400
    for _ in range(draws):
        address = bounds.random_address(rng)
        counts[address] = counts.get(address, 0) + 1
    assert len(counts) == bounds.cell_count() == 1 + (1 + 12 + 1) + spherical.wedge_count(2, 0)
    assert max(counts.values()) < 2 * min(counts.values())


@given(bounds_and_rng(), st.integers(-10, 70), st.integers(-70, 70))
def test_bounds_contains_and_describe_miss_agree(pair, shell, band):
    bounds, _seed = pair
    if bounds.contains(shell, band):
        assert 0 <= shell <= bounds.outer_shell_index and abs(band) <= bounds.high_band[shell]
    else:
        message = bounds.describe_miss(shell, band)
        assert str(shell) in message or str(band) in message


def test_empty_bounds_are_falsy_and_hold_nothing():
    bounds = galaxySkeleton.GalaxyBounds([], 4.0)
    assert not bounds
    assert bounds.cell_count() == 0
    assert bounds.random_address(random.Random(1)) is None
    assert bounds.outer_shell_index == -1
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
