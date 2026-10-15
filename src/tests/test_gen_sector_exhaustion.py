# tests/test_gen_sector_exhaustion.py

"""
`docs/TODO.md` TEST.31 (sector placement exhaustion), for
`planetgen/galaxy/sector.py`:

- "could not place a new object": the exact message, raised after exactly
  `SECTOR_MAX_PLACEMENT_ATTEMPTS` samples (a clear spot on the last
  attempt still wins, one attempt later does not), for systems and
  massive phenomena, cube and galaxy-cell sectors, leaving the sector
  untouched;
- the Poisson cap: `_sample_poisson_count` switches from Knuth's exact
  loop to the normal approximation just above
  `_POISSON_NORMAL_APPROX_MEAN`, and the approximation never goes
  negative; `grow_from_seed`'s default target is drawn from it;
- explicit positions on the real spherical cell's boundary (corners and
  faces, the core's origin and a cap's axis included) are inside and
  kept exactly, and just outside is refused for a pre-placed star;
- `nearest_neighbors` with a bad count (negative, non-integer, None).

`test_fuzz_sector_placement.py` already fuzzes random placement, growth
and the cube's boundary; this file pins the exact edges. Placement draws
from the bound draw stream (`planetgen.util.draw`), and every
test that samples replaces it with a seeded `random.Random` or a scripted
sampler.
"""

import contextlib
import math
import random
import re

import pytest

from planetgen.util import draw
from planetgen.physics import constants as pc
from planetgen import tuning as prog
from planetgen.galaxy import sector as ss
from planetgen.generation.config import SystemConfig
from planetgen.galaxy import spherical
from planetgen.galaxy.spherical import SectorCell
from planetgen.galaxy.sector import SpaceSector

ATTEMPTS = prog.SECTOR_MAX_PLACEMENT_ATTEMPTS


@contextlib.contextmanager
def _seeded_sector_rng(seed):
    with draw.bound(seed) as stream:
        yield stream


class _StubStar:
    def __init__(self, hill_ly, name="stub"):
        self.system_perimeter = hill_ly / pc.AU_TO_LY
        self.name = name


class _StubSystem:
    """Just enough of a `StarSystem` for placement (its Hill radius)."""

    def __init__(self, hill_ly, name="stub"):
        self.star = _StubStar(hill_ly, name)
        self.primary_star = self.star
        self.system_config = SystemConfig()


class _StubRemnant(_StubStar):
    """A massive phenomenon: `system_perimeter` on the object itself."""


class _StubNebula:
    name = "stub nebula"


def _full_message(edge_ly, attempts=ATTEMPTS):
    return (f"Could not place a new object far enough from every existing massive object within a "
            f"{edge_ly} ly sector after {attempts} attempts; the sector may be too full or too small.")


def _scripted(sector, points):
    """Replaces the sector's sampler with `points`, in order, and returns
    the list of points actually drawn."""
    drawn = []
    iterator = iter(points)

    def sample():
        point = next(iterator)
        drawn.append(point)
        return point

    sector._sample_point = sample
    return drawn


# ---------------------------------------------------------------------------
# "could not place a new object"
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("kind", ["system", "remnant"])
def test_full_cube_raises_the_documented_message_after_exactly_the_cap(kind):
    sector = SpaceSector("full", edge_ly=1.0)
    sector.add_system(_StubSystem(10.0), position=(0.0, 0.0, 0.0))
    drawn = _scripted(sector, [(0.1, 0.0, 0.0)] * (ATTEMPTS + 5))
    with pytest.raises(ValueError, match=f"^{re.escape(_full_message(1.0))}$"):
        if kind == "system":
            sector.add_system(_StubSystem(10.0))
        else:
            sector.add_phenomenon(_StubRemnant(10.0, "bh"), "black-hole")
    assert len(drawn) == ATTEMPTS
    assert len(sector.entries) == 1 and sector.phenomena == []


@pytest.mark.parametrize("kind", ["system", "remnant"])
def test_a_clear_spot_on_the_last_attempt_is_still_taken(kind):
    sector = SpaceSector("nearly full", edge_ly=40.0)
    sector.add_system(_StubSystem(5.0), position=(0.0, 0.0, 0.0))
    clear = (19.0, 19.0, 19.0)
    drawn = _scripted(sector, [(1.0, 0.0, 0.0)] * (ATTEMPTS - 1) + [clear])
    if kind == "system":
        entry = sector.add_system(_StubSystem(5.0))
    else:
        entry = sector.add_phenomenon(_StubRemnant(5.0, "ns"), "neutron-star")
    assert entry.position == clear
    assert len(drawn) == ATTEMPTS


def test_a_clear_spot_one_attempt_past_the_cap_is_never_reached():
    sector = SpaceSector("nearly full", edge_ly=40.0)
    sector.add_system(_StubSystem(5.0), position=(0.0, 0.0, 0.0))
    drawn = _scripted(sector, [(1.0, 0.0, 0.0)] * ATTEMPTS + [(19.0, 19.0, 19.0)])
    with pytest.raises(ValueError, match="^Could not place a new object"):
        sector.add_system(_StubSystem(5.0))
    assert len(drawn) == ATTEMPTS and len(sector.entries) == 1


def test_the_message_and_bound_follow_the_attempts_constant(monkeypatch):
    monkeypatch.setattr(prog, "SECTOR_MAX_PLACEMENT_ATTEMPTS", 3)
    sector = SpaceSector("full", edge_ly=2.5)
    sector.add_system(_StubSystem(10.0), position=(0.0, 0.0, 0.0))
    drawn = _scripted(sector, [(0.0, 0.0, 0.0)] * 10)
    with pytest.raises(ValueError, match=f"^{re.escape(_full_message(2.5, attempts=3))}$"):
        sector.add_system(_StubSystem(10.0))
    assert len(drawn) == 3


def test_a_full_galaxy_cell_raises_the_same_message():
    edge = prog.DEFAULT_SECTOR_EDGE_LY
    sector = SpaceSector("nucleus", edge_ly=edge, cell=SectorCell.for_address(0, 0, 0, edge))
    sector.add_system(_StubSystem(100.0), position=(0.0, 0.0, 0.0))
    with _seeded_sector_rng(31):
        with pytest.raises(ValueError, match=f"^{re.escape(_full_message(edge))}$"):
            sector.add_system(_StubSystem(100.0))
        with pytest.raises(ValueError, match="^Could not place a new object"):
            sector.add_phenomenon(_StubRemnant(100.0, "bh"), "black-hole")
    assert len(sector.entries) == 1 and sector.phenomena == []


def test_a_flat_separation_too_wide_for_the_sector_exhausts():
    sector = SpaceSector("flat", edge_ly=10.0)
    sector.add_system(_StubSystem(0.0), position=(0.0, 0.0, 0.0))
    with _seeded_sector_rng(1):
        with pytest.raises(ValueError, match="^Could not place a new object"):
            sector.add_system(_StubSystem(0.0), min_separation_ly=1e9)
        # Zero separation always fits on the first draw.
        sector.add_system(_StubSystem(1e9), min_separation_ly=0.0)
    assert len(sector.entries) == 2


def test_a_full_sector_still_takes_non_massive_phenomena_and_explicit_positions():
    """Only objects with a Hill sphere go through the bounded search."""
    sector = SpaceSector("full", edge_ly=1.0)
    sector.add_system(_StubSystem(10.0), position=(0.0, 0.0, 0.0))
    with _seeded_sector_rng(2):
        nebula = sector.add_phenomenon(_StubNebula(), "nebula")
        explicit = sector.add_system(_StubSystem(10.0), position=(0.1, 0.0, 0.0))
    assert sector.contains(nebula.position)
    assert explicit.position == (0.1, 0.0, 0.0)


def test_a_full_sector_stays_usable_after_exhaustion():
    sector = SpaceSector("full", edge_ly=1.0)
    sector.add_system(_StubSystem(10.0), position=(0.0, 0.0, 0.0))
    with _seeded_sector_rng(3):
        for _ in range(3):
            with pytest.raises(ValueError):
                sector.add_system(_StubSystem(10.0))
        tiny = sector.add_system(_StubSystem(0.0), min_separation_ly=0.0)
    assert sector.entries[-1] is tiny and len(sector.entries) == 2


# ---------------------------------------------------------------------------
# The Poisson cap
# ---------------------------------------------------------------------------

class _CountingRandom(random.Random):
    def __init__(self, seed):
        super().__init__(seed)
        self.uniform_draws = 0
        self.gauss_draws = 0

    def random(self):
        self.uniform_draws += 1
        return super().random()

    def gauss(self, mu=0.0, sigma=1.0):
        self.gauss_draws += 1
        return super().gauss(mu, sigma)


def test_poisson_at_the_cap_still_uses_the_exact_loop():
    cap = ss._POISSON_NORMAL_APPROX_MEAN
    rng = _CountingRandom(28)
    count = ss._sample_poisson_count(cap, rng)
    assert rng.gauss_draws == 0
    assert rng.uniform_draws == count + 1
    assert abs(count - cap) < 6 * math.sqrt(cap)
    # exp(-cap) has not underflowed, which is what capped the loop near 740.
    assert math.exp(-cap) > 0


def test_poisson_just_above_the_cap_uses_one_normal_draw():
    mean = math.nextafter(ss._POISSON_NORMAL_APPROX_MEAN, math.inf)
    rng = _CountingRandom(28)
    count = ss._sample_poisson_count(mean, rng)
    assert isinstance(count, int)
    assert rng.gauss_draws == 1
    assert abs(count - mean) < 6 * math.sqrt(mean)


def test_poisson_normal_branch_never_goes_negative():
    class _Low:
        def gauss(self, mu, sigma):
            return -1e12

        def random(self):  # pragma: no cover - must not be called
            raise AssertionError("the normal branch takes no uniform draws")

    assert ss._sample_poisson_count(1e6, _Low()) == 0


def test_poisson_normal_branch_rounds_to_the_nearest_count():
    class _Fixed:
        def __init__(self, value):
            self.value = value

        def gauss(self, mu, sigma):
            assert sigma == pytest.approx(math.sqrt(mu))
            return self.value

    assert ss._sample_poisson_count(1000.0, _Fixed(1234.4)) == 1234
    assert ss._sample_poisson_count(1000.0, _Fixed(1234.6)) == 1235


def test_grow_from_seed_draws_its_default_target_from_the_expected_count(monkeypatch):
    sector = SpaceSector("default target")
    seen = []

    def fake_poisson(mean, rng=None):
        seen.append(mean)
        return 3

    monkeypatch.setattr(ss, "_sample_poisson_count", fake_poisson)
    with _seeded_sector_rng(4):
        home = sector.add_home_system(_StubSystem(0.1))
        new = sector.grow_from_seed(home, lambda: _StubSystem(0.1))
    assert seen == [sector.expected_system_count()]
    assert len(new) == 2 and len(sector.entries) == 3


def test_grow_from_seed_with_a_zero_draw_adds_nothing(monkeypatch):
    monkeypatch.setattr(ss, "_sample_poisson_count", lambda mean, rng=None: 0)
    sector = SpaceSector("empty draw")
    with _seeded_sector_rng(5):
        home = sector.add_home_system(_StubSystem(0.1))
        assert sector.grow_from_seed(home, lambda: _StubSystem(0.1)) == []
    assert sector.entries == [home]


# ---------------------------------------------------------------------------
# Explicit positions on the cell boundary
# ---------------------------------------------------------------------------

def _cell_address(shell, kind):
    """The `(shell, band, wedge)` of a middle wedge of an ordinary band, the north cap or the core."""
    if shell == 0:
        return (0, 0, 0)
    if kind == "cap":
        return (shell, spherical.band_range(shell)[1], 0)
    return (shell, 0, spherical.wedge_count(shell, 0) // 2)


CELLS = [(0, "core"), (1, "band"), (1, "cap"), (2, "band"), (2, "cap"), (10, "band"), (3855, "band"), (3855, "cap")]
CELL_IDS = [f"shell{shell}-{kind}" for shell, kind in CELLS]


def _local(cell, point):
    """A galaxy-frame point in the cell's own sector-local frame."""
    return spherical.galaxy_to_local_pc(cell.center, point)


def _galaxy_point(radius, colatitude, longitude):
    s = math.sin(colatitude)
    return (radius * s * math.cos(longitude), radius * s * math.sin(longitude), radius * math.cos(colatitude))


def _boundary_points(cell):
    """Every corner of the cell plus the middle of each face, in the
    cell's own local frame."""
    points = []
    for r in (cell.r_inner, cell.r_outer):
        for colat in (cell.colat_north, cell.colat_south):
            for lon in (cell.lon_start, cell.lon_end):
                points.append(_galaxy_point(r, colat, lon))
    mid_colat = (cell.colat_north + cell.colat_south) / 2
    mid_lon = (cell.lon_start + cell.lon_end) / 2
    mid_r = (cell.r_inner + cell.r_outer) / 2
    points += [
        _galaxy_point(cell.r_inner, mid_colat, mid_lon),
        _galaxy_point(cell.r_outer, mid_colat, mid_lon),
        _galaxy_point(mid_r, cell.colat_north, mid_lon),
        _galaxy_point(mid_r, cell.colat_south, mid_lon),
        _galaxy_point(mid_r, mid_colat, cell.lon_start),
        _galaxy_point(mid_r, mid_colat, cell.lon_end),
    ]
    return [_local(cell, point) for point in points]


def _outside_points(cell):
    """Just past each face the cell has: a millionth of the edge outward."""
    step = 1e-6 * (cell.r_outer - cell.r_inner)
    mid_colat = (cell.colat_north + cell.colat_south) / 2
    mid_lon = (cell.lon_start + cell.lon_end) / 2
    mid_r = (cell.r_inner + cell.r_outer) / 2
    points = [_galaxy_point(cell.r_outer + step, mid_colat, mid_lon)]
    if cell.r_inner > 0:
        points.append(_galaxy_point(cell.r_inner - step, mid_colat, mid_lon))
    if cell.colat_north > 0:
        points.append(_galaxy_point(mid_r, cell.colat_north - 1e-6, mid_lon))
    if cell.colat_south < math.pi:
        points.append(_galaxy_point(mid_r, cell.colat_south + 1e-6, mid_lon))
    if cell.lon_end - cell.lon_start < 2 * math.pi - 1e-9:
        points.append(_galaxy_point(mid_r, mid_colat, cell.lon_start - 1e-6))
        points.append(_galaxy_point(mid_r, mid_colat, cell.lon_end + 1e-6))
    return [_local(cell, point) for point in points]


def _cell_sector(shell, kind):
    edge = prog.DEFAULT_SECTOR_EDGE_LY
    cell = SectorCell.for_address(*_cell_address(shell, kind), edge)
    return SpaceSector(f"shell {shell} {kind}", edge_ly=edge, cell=cell)


@pytest.mark.parametrize("shell,kind", CELLS, ids=CELL_IDS)
def test_explicit_positions_on_the_cell_boundary_are_inside_and_kept(shell, kind):
    sector = _cell_sector(shell, kind)
    for point in _boundary_points(sector.cell):
        assert sector.contains(point), point
        system = sector.add_system(_StubSystem(0.0), position=point)
        assert system.position == point
        remnant = sector.add_phenomenon(_StubRemnant(0.0, "bh"), "black-hole", position=point)
        assert remnant.position == point


@pytest.mark.parametrize("shell,kind", CELLS, ids=CELL_IDS)
def test_preplaced_stars_on_the_cell_boundary_are_accepted(shell, kind):
    sector = _cell_sector(shell, kind)
    for point in _boundary_points(sector.cell):
        entry = sector.add_preplaced_system(_StubSystem(0.0), point)
        assert entry.preplaced and entry.position == point


@pytest.mark.parametrize("shell,kind", CELLS, ids=CELL_IDS)
def test_preplaced_stars_just_outside_the_cell_are_refused(shell, kind):
    sector = _cell_sector(shell, kind)
    for point in _outside_points(sector.cell):
        assert not sector.contains(point), point
        with pytest.raises(ValueError, match="pre-placed position .* is outside the sector"):
            sector.add_preplaced_system(_StubSystem(0.0), point)
        # add_system deliberately keeps an explicit position as given.
        assert sector.add_system(_StubSystem(0.0), position=point).position == point
    assert not any(entry.preplaced for entry in sector.entries)


def test_the_cap_tip_on_the_galactic_axis_is_inside_and_the_other_pole_is_not():
    """A cap is the whole circle round the axis: the pole point of the shell, at any radius inside it, is
    in the cap, and the point across the galaxy is in the south cap, not this one."""
    sector = _cell_sector(10, "cap")
    cell = sector.cell
    assert cell.colat_north == 0.0 and cell.lon_end - cell.lon_start == pytest.approx(2 * math.pi)
    for r in (cell.r_inner, cell.r_center, cell.r_outer):
        tip = _local(cell, (0.0, 0.0, r))
        assert sector.contains(tip)
        assert sector.add_preplaced_system(_StubSystem(0.0), tip).position == tip
    assert not sector.contains(_local(cell, (0.0, 0.0, -cell.r_center)))


def test_the_core_ball_holds_its_origin_and_its_radius_but_not_a_point_past_it():
    sector = _cell_sector(0, "core")
    cell = sector.cell
    assert cell.r_inner == 0.0
    assert sector.contains(_local(cell, (0.0, 0.0, 0.0)))
    assert sector.contains(_local(cell, (cell.r_outer, 0.0, 0.0)))
    assert not sector.contains(_local(cell, (cell.r_outer * (1 + 1e-6), 0.0, 0.0)))


def test_boundary_positions_block_random_placement_like_any_other():
    """A system pinned to a cell corner still keeps later random
    placements out of its Hill sphere."""
    sector = _cell_sector(10, "band")
    corner = _boundary_points(sector.cell)[0]
    sector.add_system(_StubSystem(2.0), position=corner)
    with _seeded_sector_rng(6):
        for _ in range(5):
            entry = sector.add_system(_StubSystem(0.5))
            assert ss.distance_between(entry.position, corner) >= 2.5
            assert sector.contains(entry.position)


# ---------------------------------------------------------------------------
# nearest_neighbors with a bad count
# ---------------------------------------------------------------------------

def _line_sector(n=6):
    sector = SpaceSector("line", edge_ly=20.0)
    entries = [sector.add_system(_StubSystem(0.0), position=(float(i), 0.0, 0.0)) for i in range(n)]
    return sector, entries


@pytest.mark.parametrize("count", [-1, -2, -10**9])
def test_nearest_neighbors_negative_count_message(count):
    sector, entries = _line_sector()
    with pytest.raises(ValueError, match=f"^nearest_neighbors: count must be >= 0, got {count}$"):
        sector.nearest_neighbors(entries[0], count)


def test_nearest_neighbors_negative_count_is_refused_even_in_an_empty_sector():
    sector = SpaceSector("empty")
    lone = SpaceSector("other").add_system(_StubSystem(0.0), position=(0.0, 0.0, 0.0))
    with pytest.raises(ValueError, match="count must be >= 0"):
        sector.nearest_neighbors(lone, -1)
    assert sector.nearest_neighbors(lone, 5) == []


@pytest.mark.parametrize("count", [1.5, 2.0, math.nan, math.inf, None, "2"])
def test_nearest_neighbors_non_integer_count_is_an_error_not_a_guess(count):
    """A float, NaN, infinite, None or string count is refused with an
    exception rather than silently truncated or ignored."""
    sector, entries = _line_sector()
    with pytest.raises((TypeError, ValueError)):
        sector.nearest_neighbors(entries[0], count)
    assert len(sector.entries) == 6


def test_nearest_neighbors_zero_and_oversized_counts():
    sector, entries = _line_sector()
    assert sector.nearest_neighbors(entries[0], 0) == []
    assert sector.nearest_neighbors(entries[0], 10**9) == entries[1:]
    assert sector.nearest_neighbors(entries[0], len(entries) - 1) == entries[1:]
    assert sector.nearest_neighbors(entries[3], 2) in ([entries[2], entries[4]], [entries[4], entries[2]])


def test_nearest_neighbors_of_a_lone_system_is_empty():
    sector = SpaceSector("lone")
    entry = sector.add_system(_StubSystem(0.0), position=(0.0, 0.0, 0.0))
    assert sector.nearest_neighbors(entry) == []
    assert sector.nearest_neighbors(entry, 0) == []


def _scattered_sector(seed, count, snap=None):
    rng = random.Random(seed)
    sector = SpaceSector("scatter")
    for _ in range(count):
        position = tuple(rng.uniform(-5.0, 5.0) for _ in range(3))
        if snap:
            position = tuple(round(c / snap) * snap for c in position)   # many equal distances
        sector.entries.append(ss.SectorSystemEntry(_StubSystem(0.0), position, SystemConfig()))
    return sector


@pytest.mark.parametrize("seed,count,snap", [(1, 40, None), (2, 300, None), (3, 120, 2.5), (4, 2, None), (5, 1, None)])
@pytest.mark.parametrize("neighbors", [0, 1, 3, 7, 500])
def test_all_nearest_neighbors_match_nearest_neighbors_exactly(seed, count, snap, neighbors):
    """PERF.81: the one-tree answer is the per-entry one, ties included, for every entry."""
    sector = _scattered_sector(seed, count, snap)
    expected = [sector.nearest_neighbors(entry, neighbors) for entry in sector.entries]
    assert sector.all_nearest_neighbors(neighbors) == expected


def test_all_nearest_neighbors_refuses_a_negative_count_and_handles_an_empty_sector():
    assert SpaceSector("empty").all_nearest_neighbors(3) == []
    with pytest.raises(ValueError, match="count must be >= 0"):
        _scattered_sector(1, 5).all_nearest_neighbors(-1)
