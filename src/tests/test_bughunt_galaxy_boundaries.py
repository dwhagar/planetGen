# tests/test_bughunt_galaxy_boundaries.py

"""
Tier 1 bug-hunt coverage: boundary conditions in the spherical sector
grid (`planetgen/galaxy/spherical.py`) -- the core ball and the two polar
caps (whose cells are the whole circle round the galactic axis), the last
wedge in a band (where the longitude wraps back to 0), the shells far from
the core, and points exactly on a cell boundary.

Tier 1 here is crash-freedom and basic finiteness/sanity, plus a
cross-check of every sector center against the shell radius and band
colatitude it claims (once a Tier 2 report, now a hard assert: a systematic
drift there wouldn't crash anything but would be a real placement bug).
"""

import math

import pytest

from planetgen.galaxy import spherical
from planetgen.galaxy.spherical import (
    band_colatitude_bounds, band_range, sector_address_at, sector_cell_vertices_pc, sector_position_pc,
    shell_bounds_pc, shell_radius_pc, wedge_count,
)
from planetgen.tuning import DEFAULT_SECTOR_EDGE_PC

# Every test runs at the real sector edge (4 pc) and at a round 10 pc,
# so a seam that only lines up for one edge length still shows (TEST.5).
EDGE_LENGTHS_PC = (float(DEFAULT_SECTOR_EDGE_PC), 10.0)


@pytest.fixture(params=EDGE_LENGTHS_PC, ids=lambda edge: f"{edge:g}pc")
def edge_pc(request):
    return request.param


def _middle_bands(shell):
    """The bands to test of a shell: its equator, both caps and the bands either side of them."""
    low, high = band_range(shell)
    return sorted({low, low + 1, -1, 0, 1, high - 1, high} & set(range(low, high + 1)))


@pytest.mark.parametrize("shell", [0, 1, 2, 3])
def test_every_wedge_in_every_band_produces_finite_position(shell, edge_pc):
    low, high = band_range(shell)
    for band in range(low, high + 1):
        for wedge in range(wedge_count(shell, band)):
            assert all(math.isfinite(v) for v in sector_position_pc(shell, band, wedge, edge_pc))


@pytest.mark.parametrize("shell", [1, 2, 50, 500])
def test_first_and_last_wedge_cell_vertices_are_finite_and_nondegenerate(shell, edge_pc):
    for band in _middle_bands(shell):
        n = wedge_count(shell, band)
        for wedge in {0, n - 1}:
            vertices = sector_cell_vertices_pc(shell, band, wedge, edge_pc)
            assert len(vertices) == 8
            for vertex in vertices:
                assert all(math.isfinite(v) for v in vertex)
            if n > 1:
                # The outer face's two longitude corners are distinct points.
                assert math.dist(vertices[4], vertices[5]) > 0


def test_the_inner_corners_of_a_polar_cap_sit_on_the_axis(edge_pc):
    for shell in (2, 10, 500):
        for band in (band_range(shell)[1], band_range(shell)[0]):
            vertices = sector_cell_vertices_pc(shell, band, 0, edge_pc)
            # The pole corners (the cap's northern or southern bound is the pole) coincide on the axis.
            pole = [v for v in vertices if math.hypot(v[0], v[1]) < 1e-9]
            assert len(pole) == 4
            assert all(abs(v[2]) == pytest.approx(shell_bounds_pc(shell, edge_pc)[i // 2]) for i, v in enumerate(pole))


def test_sector_address_at_handles_exact_boundaries(edge_pc):
    # The origin, a shell boundary, the +X axis (longitude 0), just below 2*pi and each pole
    # must map to a valid address.
    assert sector_address_at((0.0, 0.0, 0.0), edge_pc) == (0, 0, 0)
    assert sector_address_at((edge_pc / 2, 0.0, 0.0), edge_pc)[0] == 1
    shell, band, wedge = sector_address_at((50.0, -1e-12, 0.0), edge_pc)
    assert 0 <= wedge < wedge_count(shell, band)
    assert wedge == wedge_count(shell, band) - 1
    assert sector_address_at((0.0, 0.0, 50.0), edge_pc)[1] == band_range(sector_address_at((0.0, 0.0, 50.0), edge_pc)[0])[1]
    assert sector_address_at((0.0, 0.0, -50.0), edge_pc)[1] == band_range(sector_address_at((0.0, 0.0, -50.0), edge_pc)[0])[0]
    # The equatorial band is band 0 and holds the plane.
    assert sector_address_at((50.0, 0.0, 0.0), edge_pc)[1] == 0


def test_sector_centers_sit_on_their_shell_and_band(edge_pc):
    for shell in (1, 2, 10, 200):
        nominal = shell_radius_pc(shell, edge_pc)
        lo, hi = shell_bounds_pc(shell, edge_pc)
        for band in _middle_bands(shell):
            north, south = band_colatitude_bounds(shell, band)
            n = wedge_count(shell, band)
            for wedge in range(0, n, max(1, n // 8)):
                position = sector_position_pc(shell, band, wedge, edge_pc)
                radius = math.sqrt(sum(v * v for v in position))
                assert radius == pytest.approx(nominal, rel=1e-9)
                assert lo < radius < hi
                if not spherical.is_cap(shell, band):
                    colatitude = math.acos(position[2] / radius)
                    assert colatitude == pytest.approx((north + south) / 2, abs=1e-9)
                else:
                    assert position[0] == position[1] == 0.0


@pytest.mark.parametrize("shell", [1, 2, 50, 500])
def test_sector_centers_map_back_to_their_own_address(shell, edge_pc):
    """A sector's center lies inside its own cell, first and last wedge
    included, in the caps and on both sides of the plane."""
    for band in _middle_bands(shell):
        n = wedge_count(shell, band)
        for wedge in {0, n // 2, n - 1}:
            position = sector_position_pc(shell, band, wedge, edge_pc)
            assert sector_address_at(position, edge_pc) == (shell, band, wedge)
