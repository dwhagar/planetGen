# tests/test_galaxy_shell_outline.py

"""
Tests for the spherical-shell outline in `planetgen.galaxy.skeleton`
(`build_shell_extents`, `GalaxyBounds`): one band range per shell, built by
bisection from the density bound, checked against an exhaustive scan of every
band of a sample of shells and pinned for the default galaxy.
"""

import math
import random

import pytest

from planetgen import tuning
from planetgen.galaxy import spherical
from planetgen.galaxy.density import build_galaxy_shape, relative_density
from planetgen.physics import mathcheck
from planetgen.galaxy.skeleton import (
    GalaxyBounds,
    _band_qualifies,
    build_shell_extents,
    expected_system_count_at_density_1,
    shell_candidate_counts,
    shell_candidate_sector_count,
)
from tests.toy_shape import toy_shape

EDGE_PC = 4.0
SMALL = toy_shape(
    disk_scale_length_pc=40.0, disk_scale_height_pc=12.0, bulge_scale_radius_pc=10.0, bulge_amplitude=2.0,
    arm_count=2, pitch_angle_rad=math.radians(15), arm_amplitude=0.4)
MILKY_WAY = build_galaxy_shape(
    disk_scale_length_pc=2600.0, disk_scale_height_pc=300.0, bulge_scale_radius_pc=1580.0,
    bulge_amplitude=mathcheck.REFERENCE_BULGE_AMPLITUDE, arm_count=2, pitch_angle_rad=math.radians(15.0),
    arm_amplitude=0.4, arm_level=1.0, nuclear_amplitude=0.0)
THRESHOLD = 1.0 / expected_system_count_at_density_1(tuning.DEFAULT_SECTOR_EDGE_LY)


@pytest.fixture(scope="module")
def milky_way():
    return build_shell_extents(MILKY_WAY, EDGE_PC, THRESHOLD)


def _qualifying_bands(shape, shell, threshold):
    return [b for b in range(0, spherical.band_range(shell)[1] + 1) if _band_qualifies(shape, EDGE_PC, shell, b, threshold)]


def test_the_default_galaxy_reaches_shell_3438_with_the_core_first(milky_way):
    """GEN.404: the galaxy holds the systems its stellar mass makes (it reached shell 3147 against 3762 for the density
    tuned to the real count at the Sun); GEN.544 gave the disk its exponential height, whose midplane is twice as dense
    for the same mass, so the edge sat at shell 3596; GEN.545 derived the mean system mass (0.483 solar masses a system
    where 0.4 had been recalled), so the galaxy holds a sixth fewer systems and the edge sat at shell 3473; the primary mass
    function of GEN.928 makes the mean system mass 0.508 (the Chabrier function and the companions it leaves out of the
    primaries, GEN.927), 5% fewer systems again, and the edge sat at shell 3441; the companions' mass ratios of GEN.951 (Moe and Di Stefano 2017,
    lighter companions for massive primaries and twins for close pairs) move the systems per solar mass by 0.3 percent, and the edge sits at shell 3438."""
    extents, outer, confirmed = milky_way
    assert outer == 3438 and confirmed
    assert [shell for shell, _high in extents] == list(range(3439))
    assert extents[0] == (0, 0)
    assert extents[1:4] == [(1, 1), (2, 3), (3, 4)]      # shells 1 to 3 are the table: whole shells in the bulge


def test_the_outline_reaches_the_axis_in_the_bulge_and_only_there(milky_way):
    extents, _outer, _confirmed = milky_way
    polar = [shell for shell, high in extents if shell and high == spherical.band_range(shell)[1]]
    assert polar[:3] == [1, 2, 3] and max(polar) == 612
    assert len(polar) == 612          # every shell out to 612 has its caps
    for shell, high in extents:
        if shell > 612:
            assert high < spherical.band_range(shell)[1]


def test_the_outline_narrows_with_distance_and_is_a_run_about_the_plane(milky_way):
    extents, _outer, _confirmed = milky_way
    highs = dict(extents)
    assert highs[3000] < highs[2000] < highs[1500]
    assert highs[3438] <= 5


def test_the_default_galaxy_holds_fewer_cells_than_the_study_counted(milky_way):
    """The study counted about 23.93 billion qualifying cells for the density tuned to the Sun; the galaxy that holds
    its stellar mass's systems (GEN.404, with the exponential disk of GEN.544 and the derived mean system mass of GEN.545)
    is 2.1 times sparser at the same shape,
    so fewer cells clear the threshold."""
    extents, _outer, _confirmed = milky_way
    shells = shell_candidate_sector_count(extents)
    assert shells < 23_930_000_000
    # GEN.927: the single-star mass function moves the derived mean system mass, so the pinned count moved 0.005 percent (from 11,271,314,727).
    # GEN.951: the mass ratios of the companions move the systems per solar mass by 0.3 percent, so the pinned count moved 0.26 percent (from 11,270,804,751).
    assert shells == 11_241_821_617


@pytest.mark.parametrize("shell", [1, 2, 3, 4, 5, 6, 7, 9, 13, 21, 55, 144, 377, 610, 799, 800, 801, 900, 1000, 1500, 2500, 3000, 3438])
def test_the_outline_is_tight_and_no_qualifying_band_lies_outside(milky_way, shell):
    extents, _outer, _confirmed = milky_way
    high = dict(extents)[shell]
    qualifying = _qualifying_bands(MILKY_WAY, shell, THRESHOLD)
    assert max(qualifying) == high               # no qualifying band beyond it, and the last one in is real
    assert qualifying[0] == 0


def test_walk_from_a_small_galaxy_matches_a_scan_of_every_shell():
    extents, outer, confirmed = build_shell_extents(SMALL, EDGE_PC, THRESHOLD)
    assert confirmed and outer == extents[-1][0]
    for shell, high in extents:
        assert max(_qualifying_bands(SMALL, shell, THRESHOLD)) == high
    assert not _band_qualifies(SMALL, EDGE_PC, outer + 1, 0, THRESHOLD)


def test_every_sector_that_clears_the_threshold_is_inside_the_outline():
    extents, _outer, _confirmed = build_shell_extents(SMALL, EDGE_PC, THRESHOLD)
    bounds = GalaxyBounds(extents, EDGE_PC)
    rng = random.Random(1)
    inside = 0
    for _ in range(4000):
        shell = rng.randrange(0, 40)
        low, high = spherical.band_range(shell)
        band = rng.randrange(low, high + 1)
        wedge = rng.randrange(spherical.wedge_count(shell, band))
        center = spherical.sector_position_pc(shell, band, wedge, EDGE_PC)
        if relative_density(center, SMALL) >= THRESHOLD:
            inside += 1
            assert bounds.qualifies(shell, band), (shell, band, wedge)
    assert inside > 100


def test_a_galaxy_too_thin_for_even_the_core_has_no_outline():
    assert build_shell_extents(SMALL, EDGE_PC, 1e9) == ([], -1, True)


def test_max_shell_cuts_the_walk_and_says_so():
    extents, outer, confirmed = build_shell_extents(SMALL, EDGE_PC, THRESHOLD, max_shell=5)
    assert outer == 5 and not confirmed and len(extents) == 6


def test_cells_in_a_band_range_add_up():
    for shell in (0, 1, 2, 3, 4, 17, 250):
        top = spherical.band_range(shell)[1]
        for high in {0, min(1, top), top // 2, top}:
            expected = sum(spherical.wedge_count(shell, b) for b in range(-high, high + 1))
            assert spherical.cells_in_band_range(shell, high) == expected
    assert spherical.cells_in_band_range(5, spherical.band_range(5)[1]) == spherical.cell_count(5)


def test_nth_cell_walks_a_range_in_order():
    shell, high = 6, 5
    cells = [spherical.nth_cell_in_band_range(shell, high, i) for i in range(spherical.cells_in_band_range(shell, high))]
    assert cells[0] == (6, 5, 0) and cells[-1] == (6, -5, spherical.wedge_count(6, -5) - 1)
    assert len(set(cells)) == len(cells)
    with pytest.raises(ValueError):
        spherical.nth_cell_in_band_range(shell, high, len(cells))


def test_galaxy_bounds_contains_describes_misses_and_counts():
    extents = [(0, 0), (1, 1), (2, 3), (3, 1)]
    bounds = GalaxyBounds(extents, EDGE_PC)
    assert bounds and bounds.outer_shell_index == 3
    assert bounds.qualifies(2, -3) and bounds.qualifies(2, 0) and not bounds.qualifies(3, 2) and not bounds.qualifies(4, 0)
    assert "shell 4" in bounds.describe_miss(4, 0)
    assert "band 2" in bounds.describe_miss(3, 2) and "-1 to 1" in bounds.describe_miss(3, 2)
    assert bounds.cell_count() == 1 + spherical.cells_in_band_range(1, 1) + spherical.cells_in_band_range(2, 3) + spherical.cells_in_band_range(3, 1)
    assert [row[2] for row in shell_candidate_counts(extents)] == [1, 14, 50, spherical.cells_in_band_range(3, 1)]


def test_galaxy_bounds_takes_stored_counts():
    bounds = GalaxyBounds([(0, 0, 7), (1, 1, 11)], EDGE_PC)
    assert bounds.cell_count() == 18


def test_an_empty_outline_is_falsy_and_offers_no_address():
    bounds = GalaxyBounds([], EDGE_PC)
    assert not bounds and bounds.outer_shell_index == -1 and bounds.cell_count() == 0
    assert bounds.random_address(random.Random(1)) is None


def test_random_addresses_are_inside_the_outline_and_uniform_by_cell():
    extents, _outer, _confirmed = build_shell_extents(SMALL, EDGE_PC, THRESHOLD)
    bounds = GalaxyBounds(extents, EDGE_PC)
    rng = random.Random(9)
    first_shells = 0
    draws = 6000
    for _ in range(draws):
        shell, band, wedge = bounds.random_address(rng)
        assert bounds.qualifies(shell, band)
        assert 0 <= wedge < spherical.wedge_count(shell, band)
        first_shells += shell <= 3
    expected = sum(spherical.cells_in_band_range(shell, high) for shell, high in extents if shell <= 3) / bounds.cell_count()
    assert first_shells / draws == pytest.approx(expected, abs=0.02)


def test_random_address_can_be_limited_to_inner_shells():
    extents, _outer, _confirmed = build_shell_extents(SMALL, EDGE_PC, THRESHOLD)
    bounds = GalaxyBounds(extents, EDGE_PC)
    rng = random.Random(4)
    assert all(bounds.random_address(rng, max_shell=2)[0] <= 2 for _ in range(300))
    assert bounds.random_address(rng, max_shell=0) == (0, 0, 0)


def test_the_default_galaxy_fits_the_object_id_and_its_reach(milky_way):
    extents, _outer, _confirmed = milky_way
    bounds = GalaxyBounds(extents, EDGE_PC)
    assert max(high for _shell, high in extents) == 961
    assert bounds.reach_shell_index == 65535
    assert spherical.wedge_count(3764, 0) < 1 << 19


def test_a_reach_the_object_id_cannot_name_or_that_stops_inside_the_outline_is_refused():
    extents, _outer, _confirmed = build_shell_extents(SMALL, EDGE_PC, THRESHOLD)
    with pytest.raises(ValueError):
        GalaxyBounds(extents, EDGE_PC, reach_shell_index=spherical.MAX_ADDRESSABLE_SHELL + 1)
    with pytest.raises(ValueError):
        GalaxyBounds(extents, EDGE_PC, reach_shell_index=max(row[0] for row in extents) - 1)


def test_contains_answers_for_the_reach_and_qualifies_for_the_outline():
    extents, outer, _confirmed = build_shell_extents(SMALL, EDGE_PC, THRESHOLD)
    bounds = GalaxyBounds(extents, EDGE_PC, reach_shell_index=outer + 50)
    high = bounds.high_band[outer]
    top = spherical.band_range(outer)[1]
    assert bounds.qualifies(outer, high) and bounds.contains(outer, high)
    if top > high:
        assert not bounds.qualifies(outer, top) and bounds.contains(outer, top)
    far = outer + 50
    assert not bounds.qualifies(far, 0) and bounds.contains(far, 0)
    assert bounds.contains(far, spherical.band_range(far)[1]) and not bounds.contains(far, spherical.band_range(far)[1] + 1)
    assert not bounds.contains(far + 1, 0) and not bounds.contains(-1, 0)
    assert bounds.random_address(random.Random(1))[0] <= outer


def test_every_cell_of_the_widest_shell_round_trips_through_its_designation_and_object_id():
    from planetgen.galaxy import object_uid
    shell = 65535
    low, high = spherical.band_range(shell)
    for band in (low, low + 1, -1, 0, 1, high - 1, high):
        for wedge in {0, spherical.wedge_count(shell, band) - 1}:
            text = spherical.provisional_sector_designation(shell, band, wedge)
            assert spherical.parse_sector_designation(text) == (shell, band, wedge)
            value = object_uid.pack(shell, band, wedge, object_uid.SERIAL_GENERATED, 1)
            assert tuple(object_uid.unpack(value))[:3] == (shell, band, wedge)
            assert int(text, 16) == value >> (object_uid.LAYOUT.serial_bits + object_uid.LAYOUT.body_bits)
    with pytest.raises(ValueError):
        spherical.provisional_sector_designation(90000, 131072, 0)
