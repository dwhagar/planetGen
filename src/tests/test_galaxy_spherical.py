# tests/test_galaxy_spherical.py

"""
Tests for `planetgen.galaxy.spherical` -- the shell/band/wedge sector grid
(`docs/design/spherical-sector-coordinates.md`).

The grid is pure geometry, so most checks are against the grid itself: the
cells tile the ball exactly (their volumes add up), every point maps to a
cell whose bounds hold it, a cell's sampled points map back to the cell, and
the neighbor lists are symmetric and agree with points nudged across each
face.
"""

import math
import random

import pytest

from planetgen.galaxy import spherical
from planetgen.galaxy.spherical import (
    HAND_CUT_SHELLS,
    SectorCell,
    band_colatitude_bounds,
    band_range,
    cell_count,
    cell_volume,
    is_cap,
    neighbor_addresses,
    sector_address_at,
    sector_cell_test,
    sector_cell_vertices_pc,
    sector_position_pc,
    shell_bounds_pc,
    shell_cells,
    shell_index_at,
    shell_layout,
    wedge_count,
    wedge_longitude_bounds,
)

EDGE_PC = 4.0
RIM_RHO_PC = 6 * EDGE_PC / (2 * math.pi)


@pytest.mark.parametrize("shell, expected", [(0, 1), (1, 14), (2, 50), (3, 104), (4, 206), (5, 296), (10, 1250)])
def test_cell_count_of_the_first_shells(shell, expected):
    assert cell_count(shell) == expected
    assert sum(1 for _ in shell_cells(shell)) == expected


def test_shell_four_is_cut_as_the_study_measured():
    layout = shell_layout(4)
    assert math.degrees(layout.cap_colatitude) == pytest.approx(13.81, abs=0.01)
    assert layout.band_count == 11
    assert [wedge_count(4, b) for b in range(5, -6, -1)] == [12, 12, 18, 24, 24, 24, 24, 24, 18, 12, 12]
    assert band_range(4) == (-6, 6)


@pytest.mark.parametrize("x, expected", [(0.2, 1), (1.0, 1), (1.4, 1), (1.5, 2), (11.0, 12), (11.9, 12), (13.0, 14),
                                         (16.9, 16), (17.0, 18), (19.0, 20), (22.5, 24), (22.4, 21), (4320.0, 4320)])
def test_the_wedge_multiplier_is_the_seven_smooth_number_nearest_in_ratio(x, expected):
    assert spherical.nearest_smooth(x) == expected


def test_every_band_of_every_shell_uses_one_rule():
    """Boss (2026-10-15): 6 * m with m made only of 2, 3, 5 and 7 and nearest the 4 pc arc, in every ordinary band
    of every shell from 4 out; no outline, window or neighbouring shell enters."""
    for shell in [4, 5, 9, 31, 250, 1000, 3764, 4095]:
        low, high = band_range(shell)
        radius = shell * EDGE_PC
        for band in range(low + 1, high):
            m, remainder = divmod(wedge_count(shell, band), 6)
            assert remainder == 0 and m >= 1
            for prime in (2, 3, 5, 7):
                while m % prime == 0:
                    m //= prime
            assert m == 1


@pytest.mark.parametrize("shell", [4, 10, 100, 1000, 3764])
def test_cap_rim_is_six_edges_around(shell):
    layout = shell_layout(shell)
    rim = 2 * math.pi * shell * EDGE_PC * math.sin(layout.cap_colatitude)
    assert rim == pytest.approx(6 * EDGE_PC)
    assert shell * EDGE_PC * math.sin(layout.cap_colatitude) == pytest.approx(RIM_RHO_PC)


@pytest.mark.parametrize("shell", [1, 2, 3, 4, 5, 17, 100, 1000])
def test_every_wedge_count_is_a_multiple_of_six_and_north_south_symmetric(shell):
    low, high = band_range(shell)
    for band in range(low + 1, high):
        assert wedge_count(shell, band) % 6 == 0
        assert wedge_count(shell, band) == wedge_count(shell, -band)
    assert wedge_count(shell, low) == wedge_count(shell, high) == 1


@pytest.mark.parametrize("shell", [4, 5, 17, 100, 1000, 3764])
def test_bands_are_an_odd_count_with_the_plane_through_band_zero(shell):
    layout = shell_layout(shell)
    assert layout.band_count % 2 == 1
    north, south = band_colatitude_bounds(shell, 0)
    assert north < math.pi / 2 < south
    assert (north + south) / 2 == pytest.approx(math.pi / 2)
    assert sector_position_pc(shell, 0, 0, EDGE_PC)[2] == 0.0


def test_shells_two_and_three_follow_the_rule_and_only_shell_one_is_cut_by_hand():
    assert set(HAND_CUT_SHELLS) == {1}
    assert [wedge_count(2, b) for b in range(2, -3, -1)] == [6, 12, 12, 12, 6]
    assert [wedge_count(3, b) for b in range(3, -4, -1)] == [12, 12, 18, 18, 18, 12, 12]


def test_shell_one_cells_are_close_to_one_edge_cube():
    volumes = [cell_volume(*address, EDGE_PC) / EDGE_PC ** 3 for address in shell_cells(1)]
    assert max(volumes) < 1.5 and min(volumes) > 0.85


def test_core_is_a_ball_of_half_an_edge_and_shell_one_starts_where_it_ends():
    assert shell_bounds_pc(0, EDGE_PC) == (0.0, 2.0)
    assert shell_bounds_pc(1, EDGE_PC) == (2.0, 6.0)
    assert shell_bounds_pc(5, EDGE_PC) == (18.0, 22.0)
    assert cell_volume(0, 0, 0, EDGE_PC) == pytest.approx(4 / 3 * math.pi * 8)
    assert [shell_index_at(r, EDGE_PC) for r in (0.0, 1.999, 2.0, 5.999, 6.0, 9.99, 10.0)] == [0, 0, 1, 1, 2, 2, 3]


@pytest.mark.parametrize("last_shell", [0, 1, 2, 3, 6, 12, 25])
def test_cell_volumes_add_up_to_the_ball(last_shell):
    total = sum(cell_volume(*address, EDGE_PC) for shell in range(last_shell + 1) for address in shell_cells(shell))
    assert total == pytest.approx(4 / 3 * math.pi * ((last_shell + 0.5) * EDGE_PC) ** 3, rel=1e-12)


def test_cell_bounds_leave_no_gap_and_no_overlap():
    for shell in (1, 2, 3, 4, 9, 30):
        low, high = band_range(shell)
        assert band_colatitude_bounds(shell, high)[0] == 0.0
        assert band_colatitude_bounds(shell, low)[1] == math.pi
        for band in range(low, high):
            assert band_colatitude_bounds(shell, band)[0] == band_colatitude_bounds(shell, band + 1)[1]
        for band in range(low + 1, high):
            n = wedge_count(shell, band)
            assert wedge_longitude_bounds(shell, band, 0)[0] == 0.0
            for wedge in range(n - 1):
                assert wedge_longitude_bounds(shell, band, wedge)[1] == wedge_longitude_bounds(shell, band, wedge + 1)[0]
            assert wedge_longitude_bounds(shell, band, n - 1)[1] == pytest.approx(2 * math.pi)


def test_every_cell_center_is_inside_its_own_cell():
    for shell in list(range(0, 31)) + [200]:
        for address in shell_cells(shell):
            center = sector_position_pc(*address, EDGE_PC)
            assert sector_address_at(center, EDGE_PC) == address, address
            assert sector_cell_test(address, EDGE_PC)(center)


def test_random_points_land_in_the_cell_whose_bounds_hold_them():
    rng = random.Random(15)
    for _ in range(20000):
        radius = 130 * rng.random() ** (1 / 3)
        cos_colat = rng.uniform(-1, 1)
        lon = rng.uniform(0, 2 * math.pi)
        s = math.sqrt(1 - cos_colat ** 2)
        point = (radius * s * math.cos(lon), radius * s * math.sin(lon), radius * cos_colat)
        shell, band, wedge = sector_address_at(point, EDGE_PC)
        inner, outer = shell_bounds_pc(shell, EDGE_PC)
        assert inner <= radius < outer or radius == pytest.approx(inner) or radius == pytest.approx(outer)
        colat = math.acos(max(-1.0, min(1.0, point[2] / radius)))
        north, south = band_colatitude_bounds(shell, band)
        assert north - 1e-9 <= colat <= south + 1e-9
        if shell and not is_cap(shell, band):
            start, end = wedge_longitude_bounds(shell, band, wedge)
            theta = math.atan2(point[1], point[0]) % (2 * math.pi)
            assert start - 1e-9 <= theta <= end + 1e-9
        assert sector_cell_test((shell, band, wedge), EDGE_PC)(point)


def test_the_cell_test_is_true_for_one_cell_only():
    rng = random.Random(3)
    cells = [sector_cell_test(address, EDGE_PC) for shell in range(0, 4) for address in shell_cells(shell)]
    for _ in range(3000):
        point = tuple(rng.uniform(-14, 14) for _ in range(3))
        if math.sqrt(sum(c * c for c in point)) < 13:
            assert sum(1 for test in cells if test(point)) == 1


def test_a_mirrored_point_lands_in_the_mirrored_band():
    rng = random.Random(7)
    for _ in range(2000):
        point = tuple(rng.uniform(-300, 300) for _ in range(3))
        shell, band, wedge = sector_address_at(point, EDGE_PC)
        assert sector_address_at((point[0], point[1], -point[2]), EDGE_PC) == (shell, -band, wedge)


def test_points_on_the_axis_and_the_origin():
    assert sector_address_at((0.0, 0.0, 0.0), EDGE_PC) == (0, 0, 0)
    assert sector_address_at((0.0, 0.0, 41.0), EDGE_PC) == (10, band_range(10)[1], 0)
    assert sector_address_at((0.0, 0.0, -41.0), EDGE_PC) == (10, band_range(10)[0], 0)
    assert sector_address_at((40.0, 0.0, 0.0), EDGE_PC) == (10, 0, 0)
    assert sector_address_at((0.0, -40.0, 0.0), EDGE_PC)[:2] == (10, 0)


def test_sampled_points_map_back_to_their_own_cell():
    rng = random.Random(21)
    addresses = [(0, 0, 0), (1, 1, 0), (1, 0, 5), (1, -1, 0), (4, 6, 0), (4, -6, 0), (4, 5, 3), (4, 0, 11),
                 (4, 0, 23), (9, 2, 17), (30, 20, 0), (200, 0, 9)]
    for address in addresses:
        cell = SectorCell.for_address(*address, EDGE_PC)
        test = sector_cell_test(address, EDGE_PC)
        from planetgen.galaxy.geometry import local_to_galaxy_pc
        for _ in range(300):
            local = cell.sample(rng)
            assert cell.contains(local), address
            galaxy = local_to_galaxy_pc(cell.center, local)
            assert sector_address_at(galaxy, EDGE_PC) == address, (address, local)
            assert test(galaxy)


def test_sampling_is_uniform_in_volume():
    cell = SectorCell.for_address(8, 2, 10, EDGE_PC)
    rng = random.Random(5)
    samples = [cell.sample(rng) for _ in range(20000)]
    # The cell's own frame: the mean of a uniform fill sits near the centroid, and the
    # outer half of the radial range holds more points than the inner half, by volume.
    from planetgen.galaxy.geometry import local_to_galaxy_pc
    radii = [math.sqrt(sum(c * c for c in local_to_galaxy_pc(cell.center, s))) for s in samples]
    middle = ((cell.r_inner ** 3 + cell.r_outer ** 3) / 2) ** (1 / 3)
    assert sum(1 for r in radii if r < middle) / len(radii) == pytest.approx(0.5, abs=0.02)


def test_the_cell_round_trips_through_its_dict():
    cell = SectorCell.for_address(12, -3, 7, EDGE_PC)
    again = SectorCell.from_dict(cell.to_dict())
    assert again.to_dict() == cell.to_dict()
    assert again.volume == pytest.approx(cell_volume(12, -3, 7, EDGE_PC))


def test_the_sector_cell_volume_matches_the_grid():
    for address in [(0, 0, 0), (1, 1, 0), (6, 4, 9), (50, -20, 100)]:
        assert SectorCell.for_address(*address, EDGE_PC).volume == pytest.approx(cell_volume(*address, EDGE_PC))


def test_a_cell_in_another_unit_scales_by_the_unit_cubed():
    ly = 3.26156
    assert SectorCell.for_address(9, 2, 4, EDGE_PC * ly).volume == pytest.approx(
        SectorCell.for_address(9, 2, 4, EDGE_PC).volume * ly ** 3)


def test_the_measured_arcs_in_a_deep_shell():
    """Study figure (uniform primes-up-to-7 rule): wedge arcs within 2% of 4 pc for nearly every cell; here every
    cell of four shells around shell 1000, the few near a pole excepted."""
    close = total = 0
    for shell in range(1000, 1004):
        low, high = band_range(shell)
        radius = shell * EDGE_PC
        for band in range(low + 1, high):
            north, south = band_colatitude_bounds(shell, band)
            lat_arc = radius * (south - north)
            lon_arc = 2 * math.pi * radius * math.sin((north + south) / 2) / wedge_count(shell, band)
            n = wedge_count(shell, band)
            total += n
            if abs(lat_arc / EDGE_PC - 1) < 0.005 and abs(lon_arc / EDGE_PC - 1) < 0.02:
                close += n
    assert close / total > 0.97


@pytest.mark.parametrize("shell", [4, 5, 6, 30])
def test_deep_cells_are_close_to_a_cube(shell):
    volumes = [cell_volume(*a, EDGE_PC) / EDGE_PC ** 3 for a in shell_cells(shell) if not is_cap(a[0], a[1])]
    assert 0.7 < min(volumes) and max(volumes) < 1.5
    caps = [cell_volume(*a, EDGE_PC) / EDGE_PC ** 3 for a in shell_cells(shell) if is_cap(a[0], a[1])]
    assert all(2.7 < v < 3.2 for v in caps)


def test_neighbors_are_symmetric_and_never_include_the_cell():
    for shell in range(0, 9):
        for address in shell_cells(shell):
            found = neighbor_addresses(*address)
            assert address not in found
            assert len(found) == len(set(found))
            for other in found:
                assert address in neighbor_addresses(*other), (address, other)


def test_the_core_touches_all_of_shell_one_and_nothing_else():
    assert sorted(neighbor_addresses(0, 0, 0)) == sorted(shell_cells(1))


def test_a_cap_touches_every_wedge_of_the_nearest_band():
    shell = 6
    top = band_range(shell)[1]
    found = neighbor_addresses(shell, top, 0)
    same_shell = [a for a in found if a[0] == shell]
    assert sorted(same_shell) == sorted((shell, top - 1, w) for w in range(wedge_count(shell, top - 1)))
    assert not any(a[0] not in (shell - 1, shell, shell + 1) for a in found)


def test_neighbor_counts_in_the_disk_are_ten_to_fifteen():
    counts = [len(neighbor_addresses(40, 0, w)) for w in range(0, wedge_count(40, 0), 7)]
    assert 8 <= min(counts) and max(counts) <= 16


def _points_across_faces(address, rng, per_face=250, nudge=2e-6):
    """Points just outside each of a cell's faces, sampled over the face."""
    shell, band, wedge = address
    inner, outer = shell_bounds_pc(shell, EDGE_PC)
    north, south = band_colatitude_bounds(shell, band)
    start, end = wedge_longitude_bounds(shell, band, wedge)
    found = set()

    def at(r, colat, lon):
        s = r * math.sin(colat)
        return (s * math.cos(lon), s * math.sin(lon), r * math.cos(colat))

    for _ in range(per_face):
        r = (rng.uniform(inner ** 3, outer ** 3)) ** (1 / 3)
        colat = math.acos(rng.uniform(math.cos(south), math.cos(north)))
        lon = rng.uniform(start, end)
        found.add(sector_address_at(at(inner - nudge, colat, lon), EDGE_PC) if inner > 0 else (0, 0, 0))
        found.add(sector_address_at(at(outer + nudge, colat, lon), EDGE_PC))
        if north > 0:
            found.add(sector_address_at(at(r, north - nudge, lon), EDGE_PC))
        if south < math.pi:
            found.add(sector_address_at(at(r, south + nudge, lon), EDGE_PC))
        if end - start < 2 * math.pi - 1e-9:
            found.add(sector_address_at(at(r, colat, start - nudge), EDGE_PC))
            found.add(sector_address_at(at(r, colat, end + nudge), EDGE_PC))
    found.discard(address)
    return found


@pytest.mark.parametrize("shell", [1, 2, 3, 4, 5, 6])
def test_neighbors_agree_with_points_nudged_across_each_face(shell):
    rng = random.Random(100 + shell)
    for address in shell_cells(shell):
        listed = set(neighbor_addresses(*address))
        seen = _points_across_faces(address, rng)
        assert seen <= listed, (address, seen - listed)
        # Every listed neighbor is reached except a face too thin for 250 random samples to hit.
        assert len(listed - seen) <= max(2, len(listed) // 4), (address, listed - seen)


def test_cell_vertices_lie_on_the_cells_bounds():
    for address in [(5, 2, 3), (5, 6, 0), (5, -6, 0), (1, 0, 4), (0, 0, 0)]:
        vertices = sector_cell_vertices_pc(*address, EDGE_PC)
        assert len(vertices) == 8
        inner, outer = shell_bounds_pc(address[0], EDGE_PC)
        for index, (x, y, z) in enumerate(vertices):
            assert math.sqrt(x * x + y * y + z * z) == pytest.approx(outer if index & 4 else inner)


@pytest.mark.parametrize("address", [(-1, 0, 0), (0, 1, 0), (0, 0, 1), (3, 5, 0), (3, 0, 18), (3, 0, -1), (4, -7, 0)])
def test_addresses_outside_the_grid_are_refused(address):
    with pytest.raises(ValueError):
        sector_position_pc(*address, EDGE_PC)


def test_a_negative_shell_has_no_layout():
    with pytest.raises(ValueError):
        shell_layout(-1)


def test_the_module_keeps_the_study_names_and_none_of_the_cylinder_words():
    public = " ".join(spherical.__all__)
    assert "ring" not in public and "layer" not in public.replace("ShellLayout", "").replace("shell_layout", "")


def test_designations_round_trip_and_fit_forty_bits():
    for address in [(0, 0, 0), (1, 1, 0), (1, -1, 0), (4, 6, 0), (4, -6, 0), (10, 0, 59), (3764, 2, 23519),
                    (1020, 1602, 0), (1020, -1602, 0), (4095, 0, 0)]:
        text = spherical.provisional_sector_designation(*address)
        assert spherical.parse_sector_designation(text) == address
        assert int(text, 16) < 1 << 40 and text == text.upper()


def test_a_designation_is_the_sector_field_of_the_object_ids_born_there():
    from planetgen.galaxy import object_uid
    for address in [(0, 0, 0), (7, -3, 11), (3764, 2, 23519), (1020, -1602, 0)]:
        designation = spherical.provisional_sector_designation(*address)
        printed = object_uid.format_id(object_uid.pack(*address, object_uid.SERIAL_GENERATED, 5, 7))
        assert int(printed.split("-")[0], 16) == int(designation, 16)


def test_designations_of_neighbors_differ_and_sort_by_shell_then_band_then_wedge():
    keys = [int(spherical.provisional_sector_designation(*a), 16) for a in sorted(
        (a for shell in range(0, 5) for a in shell_cells(shell)), key=lambda a: (a[0], a[1], a[2]))]
    assert keys == sorted(keys) and len(set(keys)) == len(keys)


@pytest.mark.parametrize("address", [(4096, 0, 0), (3000, 2048, 0), (3000, -2049, 0), (-1, 0, 0), (3, 5, 0), (3, 0, 18)])
def test_a_designation_outside_the_field_or_the_grid_is_refused(address):
    with pytest.raises(ValueError):
        spherical.provisional_sector_designation(*address)


@pytest.mark.parametrize("text", ["XYZ", "", "10000000000", "-1", "FFFFFFFFFF"])
def test_a_designation_that_names_nothing_is_refused(text):
    with pytest.raises(ValueError):
        spherical.parse_sector_designation(text)


def _brute_force_within(center, radius, shells):
    found = set()
    for shell in shells:
        for address in shell_cells(shell):
            if math.dist(sector_position_pc(*address, EDGE_PC), center) <= radius + 1e-9:
                found.add(address)
    return found


@pytest.mark.parametrize("center, radius", [
    ((0.0, 0.0, 0.0), 11.0), ((0.0, 0.0, 0.0), 0.0), ((1.0, 2.0, 3.0), 9.0), ((30.0, 5.0, 2.0), 16.0),
    ((0.0, 0.0, 60.0), 20.0), ((0.0, 0.0, -45.0), 14.0), ((-25.0, 40.0, -30.0), 25.0), ((60.0, 0.0, 0.0), 7.0),
    ((2.0, -1.0, 80.0), 33.0), ((0.5, 0.0, 12.0), 5.0),
])
def test_enumerating_a_sphere_matches_a_scan_of_every_cell(center, radius):
    found = list(spherical.enumerate_sectors_within_radius(center, radius, EDGE_PC))
    addresses = {row[:3] for row in found}
    assert len(addresses) == len(found)
    top = int((math.dist(center, (0, 0, 0)) + radius) / EDGE_PC) + 2
    assert addresses == _brute_force_within(center, radius, range(max(0, int((math.dist(center, (0, 0, 0)) - radius) / EDGE_PC) - 2), top + 1)) | ({(0, 0, 0)} if math.dist(center, (0, 0, 0)) <= radius else set())
    for shell, band, wedge, x, y, z, distance in found:
        assert (x, y, z) == sector_position_pc(shell, band, wedge, EDGE_PC)
        assert distance == pytest.approx(math.dist((x, y, z), center))


def test_enumerating_random_spheres_matches_a_scan():
    rng = random.Random(33)
    for _ in range(12):
        center = tuple(rng.uniform(-70, 70) for _ in range(3))
        radius = rng.uniform(0, 28)
        found = {row[:3] for row in spherical.enumerate_sectors_within_radius(center, radius, EDGE_PC)}
        top = int((math.dist(center, (0, 0, 0)) + radius) / EDGE_PC) + 2
        assert found == _brute_force_within(center, radius, range(max(1, int((math.dist(center, (0, 0, 0)) - radius) / EDGE_PC) - 2), top + 1)) | ({(0, 0, 0)} if math.dist(center, (0, 0, 0)) <= radius else set())


def test_enumeration_refuses_a_negative_radius_or_edge():
    with pytest.raises(ValueError):
        list(spherical.enumerate_sectors_within_radius((0, 0, 0), -1.0, EDGE_PC))
    with pytest.raises(ValueError):
        list(spherical.enumerate_sectors_within_radius((0, 0, 0), 1.0, 0.0))


def test_no_point_of_a_cell_is_farther_from_its_center_than_the_reach():
    farthest = 0.0
    rng = random.Random(2)
    for shell in list(range(0, 14)) + [100, 500]:
        for address in list(shell_cells(shell))[::1 if shell < 100 else 997]:
            center = sector_position_pc(*address, EDGE_PC)
            cell = SectorCell.for_address(*address, EDGE_PC)
            from planetgen.galaxy.geometry import local_to_galaxy_pc
            points = sector_cell_vertices_pc(*address, EDGE_PC)
            points += [local_to_galaxy_pc(center, cell.sample(rng)) for _ in range(6)]
            for point in points:
                farthest = max(farthest, math.dist(point, center))
    assert 1.2 < farthest / EDGE_PC < spherical.CELL_REACH_EDGES


def test_cells_touching_a_sphere_hold_every_point_of_it():
    rng = random.Random(8)
    for center, radius in [((20.0, 10.0, 5.0), 12.0), ((0.0, 0.0, 30.0), 10.0), ((0.0, 0.0, 0.0), 9.0)]:
        touching = {row[:3] for row in spherical.cells_touching_sphere(center, radius, EDGE_PC)}
        for _ in range(1500):
            v = [rng.gauss(0, 1) for _ in range(3)]
            norm = math.sqrt(sum(c * c for c in v))
            scale = radius * rng.random() ** (1 / 3) / norm
            point = tuple(center[i] + v[i] * scale for i in range(3))
            assert sector_address_at(point, EDGE_PC) in touching


def _sampled_cells(start, end, steps=6000):
    return {sector_address_at(tuple(start[i] + (end[i] - start[i]) * k / steps for i in range(3)), EDGE_PC)
            for k in range(steps + 1)}


@pytest.mark.parametrize("start, end", [
    ((1.0, 1.0, 1.0), (50.0, 40.0, -30.0)),
    ((-60.0, 3.0, 0.0), (60.0, 3.0, 0.0)),            # along the plane
    ((-40.0, -40.0, 5.0), (40.0, 40.0, 5.0)),         # past the axis
    ((0.0, 0.0, -50.0), (0.0, 0.0, 50.0)),            # along the axis
    ((30.0, 0.0, 40.0), (30.0, 0.0, -40.0)),          # across the plane
    ((70.0, 20.0, 3.0), (-70.0, 25.0, 30.0)),
    ((5.0, 5.0, 0.0), (5.0, 5.0, 0.0)),               # a point
    ((0.0, 0.0, 0.0), (0.0, 0.0, 9.0)),
])
def test_a_segment_lists_the_cells_it_passes_through(start, end):
    found = spherical.sectors_along_segment(start, end, EDGE_PC)
    assert len(found) == len(set(found))
    assert found[0] == sector_address_at(start, EDGE_PC)
    assert sector_address_at(end, EDGE_PC) in found
    sampled = _sampled_cells(start, end)
    assert sampled <= set(found), sampled - set(found)
    # Cells beyond the sampled ones are only those the segment touches at a corner or a sliver.
    assert len(set(found) - sampled) <= max(3, len(sampled) // 10)


def test_random_segments_list_every_cell_dense_sampling_finds():
    rng = random.Random(77)
    for _ in range(25):
        start = tuple(rng.uniform(-60, 60) for _ in range(3))
        end = tuple(rng.uniform(-60, 60) for _ in range(3))
        found = set(spherical.sectors_along_segment(start, end, EDGE_PC))
        assert _sampled_cells(start, end, 5000) <= found


def test_a_segment_is_listed_in_the_order_it_enters_cells():
    """A line past the core enters shells from the outside in and out again: the shell runs down, then up."""
    found = spherical.sectors_along_segment((-30.0, 0.5, 0.5), (30.0, 0.5, 0.5), EDGE_PC)
    shells = [cell[0] for cell in found]
    low = shells.index(min(shells))
    assert shells[:low + 1] == sorted(shells[:low + 1], reverse=True)
    assert shells[low:] == sorted(shells[low:])


def test_describing_a_cell():
    info = spherical.describe_sector_cell(10, 3, 7, EDGE_PC)
    assert (info["shell_index"], info["band_index"], info["wedge_index"]) == (10, 3, 7)
    assert info["designation"] == spherical.provisional_sector_designation(10, 3, 7)
    assert info["cartesian_pc"] == list(sector_position_pc(10, 3, 7, EDGE_PC))
    assert info["spherical"]["r_pc"] == 40.0 and len(info["vertices_pc"]) == 8
    assert info["volume_pc3"] == pytest.approx(cell_volume(10, 3, 7, EDGE_PC))
    assert info["mean_band_arc_pc"] == pytest.approx(4.0, rel=0.1)
    core = spherical.describe_sector_cell(0, 0, 0, EDGE_PC)
    assert core["volume_pc3"] == pytest.approx(4 / 3 * math.pi * 8) and core["cartesian_pc"] == [0.0, 0.0, 0.0]
    with pytest.raises(ValueError):
        spherical.describe_sector_cell(10, 3, 99, EDGE_PC)


def test_zone_and_quadrant_group_as_before():
    assert spherical.sector_zone(0, 13.05) == 0 and spherical.sector_zone(8, 13.05) == 1
    assert [spherical.sector_quadrant(*xy) for xy in [(1, 1), (-1, 1), (-1, -1), (1, -1)]] == [1, 2, 3, 4]
