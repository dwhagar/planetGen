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
    SMALL_SHELL_TABLE,
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


@pytest.mark.parametrize("shell, expected", [(0, 1), (1, 14), (2, 50), (3, 104), (4, 194), (5, 296), (10, 1250)])
def test_cell_count_of_the_first_shells(shell, expected):
    assert cell_count(shell) == expected
    assert sum(1 for _ in shell_cells(shell)) == expected


def test_shell_four_is_cut_as_the_study_measured():
    layout = shell_layout(4)
    assert math.degrees(layout.cap_colatitude) == pytest.approx(13.81, abs=0.01)
    assert layout.band_count == 11
    assert [wedge_count(4, b) for b in range(5, -6, -1)] == [6, 12, 18, 24, 24, 24, 24, 24, 18, 12, 6][::1]
    assert band_range(4) == (-6, 6)


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


def test_shells_two_and_three_in_the_table_are_the_rim_rule():
    for shell in (2, 3):
        degrees, wedges = SMALL_SHELL_TABLE[shell]
        assert degrees is None
        alpha = math.asin(RIM_RHO_PC / (shell * EDGE_PC))
        n = 2 * int(math.floor(shell * (math.pi - 2 * alpha) / 2)) + 1
        assert n == len(wedges)
        step = (math.pi - 2 * alpha) / n
        for i, table in enumerate(wedges):
            centre = alpha + (i + 0.5) * step
            assert table == 6 * max(1, int(math.floor(2 * math.pi * shell * math.sin(centre) / 6 + 0.5)))


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
    """Study figure: about 99.8% of cells have their band and wedge arcs within 0.5% of 4 pc, in the real galaxy;
    here every cell of four shells around shell 1000."""
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
            if abs(lat_arc / EDGE_PC - 1) < 0.005 and abs(lon_arc / EDGE_PC - 1) < 0.005:
                close += n
    assert close / total > 0.998


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
