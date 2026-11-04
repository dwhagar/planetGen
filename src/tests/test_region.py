"""MAP.151: the round region as shell, band and wrapping wedge-arc ranges."""
import math

import pytest

from planetgen.galaxy import spherical
from planetgen.galaxy.region import BALL_EDGE_ALLOWANCE, BandArc, Region
from planetgen.galaxy.skeleton import GalaxyBounds

EDGE = 4.0


def _walked(region):
    """The region by brute force: every sector of the ball the sphere walk finds."""
    return {
        (shell, band, wedge)
        for shell, band, wedge, _x, _y, _z, _d
        in spherical.enumerate_sectors_within_radius(region.center_pc, region.reach_pc, region.edge_pc)
    }


def _bounds(last_shell, limit=None):
    """An outline that holds every band of shells 0 to `last_shell` (or bands -limit..limit)."""
    rows = []
    for shell in range(last_shell + 1):
        high = spherical.band_range(shell)[1]
        rows.append((shell, high if limit is None else min(high, limit)))
    return GalaxyBounds(rows, EDGE)


CENTRES = [
    (0, 0, 0),       # the core
    (1, 0, 0),
    (2, 1, 0),
    (12, 0, 5),
    (12, 3, 7),
    (40, -2, 0),     # wedge 0 on a bigger shell
    (40, 3, 20),
    (150, 0, 400),
]


def _centre_of(address):
    return spherical.sector_position_pc(*address, EDGE)


@pytest.mark.parametrize("centre", CENTRES)
@pytest.mark.parametrize("sectors", [0.5, 1, 2, 5, 13, 30])
def test_a_region_is_the_sectors_a_walk_of_the_sphere_finds(centre, sectors):
    region = Region.around_sector(*centre, sectors, EDGE)
    got = list(region.addresses())
    assert len(got) == len(set(got)) == region.count()
    assert set(got) == _walked(region)


def test_the_centre_sector_is_always_in_and_a_far_sector_is_not():
    region = Region.around_sector(40, 0, 0, 1, EDGE)
    assert region.contains(40, 0, 0)
    assert region.count() >= 1
    assert not region.contains(60, 0, 0)
    assert not region.contains(0, 0, 0)


def test_an_arc_that_crosses_wedge_zero_wraps_and_splits_for_sql():
    region = Region.around_sector(40, 0, 0, 3, EDGE)
    arcs = [arc for arc in region.band_arcs() if arc.shell == 40 and arc.band == 0]
    assert len(arcs) == 1
    arc = arcs[0]
    n = spherical.wedge_count(40, 0)
    assert arc.first > 0 and arc.first + arc.count > n and not arc.whole
    first, second = region.wedge_ranges(arc)
    assert first == (arc.first, n - 1) and second == (0, arc.first + arc.count - 1 - n)
    wedges = [wedge for shell, band, wedge in region.addresses() if (shell, band) == (40, 0)]
    assert wedges[0] == arc.first and set(wedges) == set(range(*first[:1], n)) | set(range(0, second[1] + 1))


def test_a_region_over_the_core_holds_whole_bands_and_the_core():
    region = Region((0, 0, 0), 6 * EDGE, EDGE)
    arcs = region.band_arcs()
    assert arcs[0] == BandArc(0, 0, 0, 1, True)
    assert all(arc.whole for arc in arcs)
    assert {arc.shell for arc in arcs} == set(range(0, 7))


def test_a_region_beside_the_axis_has_whole_caps_and_arcs_in_the_bands():
    region = Region.around_sector(30, 0, 0, 4, EDGE)
    partial = [arc for arc in region.band_arcs() if not arc.whole]
    assert partial and all(arc.count < spherical.wedge_count(arc.shell, arc.band) for arc in partial)
    assert set(region.addresses()) == _walked(region)


def test_a_region_on_a_pole_holds_the_cap():
    north = (30, spherical.band_range(30)[1], 0)
    region = Region.around_sector(*north, 1, EDGE)
    assert BandArc(30, north[1], 0, 1, True) in region.band_arcs()
    assert set(region.addresses()) == _walked(region)


def test_the_count_and_cells_follow_the_outline():
    region = Region.around_sector(12, 0, 5, 6, EDGE)
    bounds = _bounds(40, limit=1)
    listed = list(region.addresses(bounds))
    assert region.count(bounds) == len(listed) < region.count()
    assert all(bounds.contains(shell, band) for shell, band, _wedge in listed)
    assert len(list(region.cells(bounds))) == len(listed)


def test_cells_carry_their_centre_and_distance():
    region = Region.around_sector(12, 3, 7, 3, EDGE)
    for shell, band, wedge, x, y, z, distance in region.cells():
        assert (x, y, z) == pytest.approx(_centre_of((shell, band, wedge)))
        assert distance == pytest.approx(math.dist((x, y, z), region.center_pc))
        assert distance <= region.reach_pc + 1e-9


def test_spans_describe_the_same_sectors():
    region = Region.around_sector(40, 0, 0, 3, EDGE)
    bounds = _bounds(60)
    from_spans = set()
    for span in region.spans():
        from_spans.update(span.addresses(bounds))
    assert from_spans == set(region.addresses(bounds))
    assert sum(span.count(bounds) for span in region.spans()) == region.count(bounds)


def test_a_centre_off_a_sector_centre_uses_its_own_sector():
    position = (60.0, 25.0, 9.0)
    region = Region(position, 2 * EDGE, EDGE)
    assert region.contains(*spherical.sector_address_at(position, EDGE))
    assert set(region.addresses()) == _walked(region)


def test_the_allowance_matches_the_fill_that_uses_it():
    from planetgen.galaxy import region_run
    assert region_run.reach_pc(EDGE, EDGE) == pytest.approx((1 + BALL_EDGE_ALLOWANCE) * EDGE)


def test_a_radius_of_one_sector_holds_the_centre_and_its_face_neighbours():
    for address in [(12, 0, 5), (40, -2, 0), (150, 3, 100)]:
        region = Region.around_sector(*address, 1, EDGE)
        neighbours = set(spherical.neighbor_addresses(*address))
        inside = set(region.addresses())
        assert address in inside
        assert len(neighbours & inside) >= 0.9 * len(neighbours)


@pytest.mark.parametrize("kwargs", [
    {"reach_pc": -1, "edge_pc": 4}, {"reach_pc": 1, "edge_pc": 0},
])
def test_a_region_refuses_nonsense(kwargs):
    with pytest.raises(ValueError):
        Region((0, 0, 0), **kwargs)
