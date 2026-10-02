# tests/test_fuzz_sector_placement.py

"""
Property-based fuzzing of `stellarObjects/spaceSector.py`: random and
explicit placement (`add_system`/`add_phenomenon`/`add_home_system`),
Poisson-disk growth (`grow_from_seed`), `nearest_neighbors`, the
Poisson count sampler, the galaxy-cell geometry sectors place into, and
the `to_dict`/`from_dict`/`save`/`load` round trip.

`test_bughunt_sector_placement.py` crams a handful of real large systems
into one small cube. This file instead drives the placement code with
hypothesis-chosen *geometry*: sector edges from 0 to 10^4 ly, real
cylindrical galaxy cells from ring 0 (a wedge touching the galactic axis)
to ring 10^6 (a thin slab), Hill radii from 0 to far larger than the
sector, flat `min_separation_ly` overrides including 0, negative and NaN,
target counts from 0 to 10^6, `k` from 0 up, and mixes of star systems
with massive (black hole / neutron star) and non-massive phenomena.

Placement geometry only needs an object's Hill radius, so most tests use
tiny stand-ins (`_StubSystem`/`_StubRemnant`) whose radius hypothesis
picks directly -- thousands of placements per second. The round-trip
tests use real `StarSystem`s and real phenomena.

Reproducibility: `spaceSector._rng` reads the module-level `random`
stream (GEN.39). `_seeded_sector_rng` shadows its `random`/
`getrandbits` methods on the instance with a seeded `random.Random`'s, so
`uniform`/`choice`/`random` -- and every helper that captured `_rng` as a
default argument at import time -- become a pure function of the
hypothesis-drawn seed.

Invariants: every placed entry lies inside the sector (`contains`) at a
finite position; every pair of massive objects is at least its required
separation apart; placement and growth always terminate quickly; the
round trip is lossless; `nearest_neighbors` agrees with brute force.
"""

import contextlib
import json
import math
import random
import signal
import time
from unittest import mock

import pytest
from hypothesis import HealthCheck, assume, example, given, note, settings
from hypothesis import strategies as st

from stellarObjects import physical_constants as pc, program_constants as prog
from stellarObjects import spaceSector as ss
from stellarObjects.config import SystemConfig
from stellarObjects.galaxyGeometry import SectorCell, ring_sector_count
from stellarObjects.spaceSector import SpaceSector, distance_between, required_separation_ly
from stellarObjects.systemData import StarSystem

from tests.fuzz_support import any_float, deterministic_entropy as _seeded_generation, hostile_text, scaled

DEFAULT_EDGE = prog.DEFAULT_SECTOR_EDGE_LY
SEEDS = st.integers(min_value=0, max_value=2**64 - 1)
TOL = prog.SECTOR_GROWTH_FLOATING_POINT_TOLERANCE_LY


# ---------------------------------------------------------------------------
# Harness: determinism, time limits, stand-in objects
# ---------------------------------------------------------------------------

@contextlib.contextmanager
def _seeded_sector_rng(seed):
    rng = random.Random(seed)
    with mock.patch.object(ss._rng, "random", rng.random), \
            mock.patch.object(ss._rng, "getrandbits", rng.getrandbits):
        yield rng


class _Timeout(Exception):
    pass


@contextlib.contextmanager
def _time_limit(seconds):
    """Raises `_Timeout` inside the block after `seconds` of wall time --
    the placement loops are pure Python, so SIGALRM interrupts them."""
    def _raise(signum, frame):
        raise _Timeout(f"did not finish within {seconds}s")

    previous = signal.signal(signal.SIGALRM, _raise)
    signal.setitimer(signal.ITIMER_REAL, seconds)
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous)


class _StubStar:
    def __init__(self, hill_ly, name):
        self.system_perimeter = hill_ly / pc.AU_TO_LY
        self.name = name


class _StubSystem:
    """Just enough of a `StarSystem` for placement: `hill_radius_ly` reads
    `.star.system_perimeter`, entries read `.system_config`."""

    def __init__(self, hill_ly, name="stub"):
        self.star = _StubStar(hill_ly, name)
        self.primary_star = self.star
        self.system_config = SystemConfig()


class _StubRemnant(_StubStar):
    """A massive phenomenon: exposes `system_perimeter` directly, exactly
    like `compactRemnant.BlackHole`/`NeutronStar`."""


class _StubNebula:
    name = "stub nebula"


def _massive(sector):
    """Every placed object that has a Hill sphere, with its position."""
    return list(sector._massive_neighbors())


def assert_sector_sane(sector, *, min_separation_ly=None):
    problems = []
    for entry in list(sector.entries) + list(sector.phenomena):
        if not all(math.isfinite(c) for c in entry.position):
            problems.append(f"non-finite position {entry.position}")
        elif not sector.contains(entry.position):
            problems.append(f"position {entry.position} outside the sector")
    massive = _massive(sector)
    for i in range(len(massive)):
        for j in range(i + 1, len(massive)):
            (a, pa), (b, pb) = massive[i], massive[j]
            required = required_separation_ly(a, b) if min_separation_ly is None else min_separation_ly
            dist = distance_between(pa, pb)
            if dist < required - max(TOL, 1e-12 * required):
                problems.append(f"objects {i} and {j} are {dist} ly apart, need {required}")
    assert not problems, "\n".join(problems[:10])


# ---------------------------------------------------------------------------
# Geometry strategies
# ---------------------------------------------------------------------------

EDGES = st.one_of(
    st.sampled_from([1e-6, 0.5, 1.0, DEFAULT_EDGE, 50.0, 1e4]),  # edge <= 0 is rejected (see below)
    st.floats(min_value=1e-3, max_value=200.0),
)
HILL_LY = st.one_of(
    st.sampled_from([0.0, 1e-12, 0.05, 0.5, 1.71, 7.5, 1e3, 1e6]),
    st.floats(min_value=0.0, max_value=20.0),
)
RINGS = st.one_of(st.sampled_from([0, 1, 2, 3, 100, 10**6]), st.integers(0, 5000))


@st.composite
def sectors(draw):
    edge = draw(EDGES)
    cell = None
    if edge > 0 and draw(st.booleans()):
        cell = SectorCell.for_ring(draw(RINGS), edge)
    name = draw(hostile_text)
    return SpaceSector(name, edge_ly=edge, cell=cell)


# ---------------------------------------------------------------------------
# Random placement
# ---------------------------------------------------------------------------

PLACEMENTS = st.lists(
    st.tuples(st.sampled_from(["system", "remnant", "nebula"]), HILL_LY), max_size=14,
)
MIN_SEPARATION = st.one_of(st.none(), st.sampled_from([0.0, -1.0, 1e-9, 0.5, 5.0, 1e9, math.nan]))


@settings(max_examples=scaled(120))
@given(seed=SEEDS, sector=sectors(), placements=PLACEMENTS, min_sep=MIN_SEPARATION)
def test_random_placement_stays_inside_and_separated(seed, sector, placements, min_sep):
    note(f"edge={sector.edge_ly} cell={vars(sector.cell) if sector.cell else None}")
    placed_massive = 0
    with _seeded_sector_rng(seed), _time_limit(10):
        for kind, hill in placements:
            try:
                if kind == "system":
                    sector.add_system(_StubSystem(hill), min_separation_ly=min_sep)
                    placed_massive += 1
                elif kind == "remnant":
                    sector.add_phenomenon(_StubRemnant(hill, "bh"), "black-hole", min_separation_ly=min_sep)
                    placed_massive += 1
                else:
                    sector.add_phenomenon(_StubNebula(), "nebula")
            except ValueError:
                pass  # documented: no room left for this massive object
    flat = min_sep if (min_sep is not None and math.isfinite(min_sep)) else None
    assert_sector_sane(sector, min_separation_ly=flat)
    if min_sep is not None and math.isnan(min_sep):
        # Every distance comparison against NaN is False: only the very
        # first massive object (no neighbors to compare against) fits.
        assert len(_massive(sector)) <= 1


def test_full_sector_raises_after_bounded_attempts_not_hang():
    sector = SpaceSector("full", edge_ly=1.0)
    with _seeded_sector_rng(7), _time_limit(5):
        sector.add_system(_StubSystem(10.0))
        started = time.monotonic()
        with pytest.raises(ValueError):
            sector.add_system(_StubSystem(10.0))
        with pytest.raises(ValueError):
            sector.add_phenomenon(_StubRemnant(10.0, "bh"), "black-hole")
    assert time.monotonic() - started < 5.0
    assert len(sector.entries) == 1 and sector.phenomena == []


@settings(max_examples=scaled(60))
@given(seed=SEEDS, sector=sectors(), jitter=st.one_of(
    st.sampled_from([0.0, -1.0, 1e-9, 1e9, math.nan, math.inf, -math.inf]), st.floats(0.0, 20.0)))
def test_home_system_is_always_inside(seed, sector, jitter):
    with _seeded_sector_rng(seed):
        entry = sector.add_home_system(_StubSystem(1.0), jitter_ly=jitter)
    assert all(math.isfinite(c) for c in entry.position)
    assert sector.contains(entry.position)


@settings(max_examples=scaled(40))
@given(sector=sectors(), signs=st.tuples(*[st.sampled_from([-1.0, 0.0, 1.0])] * 3))
def test_explicit_positions_on_the_cube_boundary_are_inside(sector, signs):
    assume(sector.cell is None)
    half = sector.edge_ly / 2
    position = tuple(s * half for s in signs)
    entry = sector.add_system(_StubSystem(1.0), position=position)
    assert entry.position == position and sector.contains(position)


@settings(max_examples=scaled(30))
@given(bad=st.tuples(any_float, any_float, any_float).filter(lambda p: not all(map(math.isfinite, p))),
       method=st.sampled_from(["add_system", "add_phenomenon"]))
@example(bad=(math.nan, 0.0, 0.0), method="add_system")
@example(bad=(0.0, math.inf, 0.0), method="add_system")
@example(bad=(0.0, 0.0, -math.inf), method="add_phenomenon")
def test_explicit_non_finite_position_is_rejected(bad, method):
    """Regression: a NaN position was accepted, after which every
    `distance >= required` check was False and the sector was poisoned."""
    sector = SpaceSector("poison")
    with pytest.raises(ValueError):
        if method == "add_system":
            sector.add_system(_StubSystem(1.0), position=bad)
        else:
            sector.add_phenomenon(_StubRemnant(1.0, "bh"), "black-hole", position=bad)
    assert not sector.entries and not sector.phenomena
    sector.add_system(_StubSystem(1.0))  # still usable


@pytest.mark.parametrize("position", [(1e9, 0.0, 0.0), (5.000001, 0.0, 0.0), (-6.0, -6.0, -6.0)])
def test_explicit_position_outside_the_sector_is_kept_as_given(position):
    """Containment is deliberately not enforced for an explicit position
    (generate.py anchors the quasar at the galactic centre's sector-local
    offset, and older tests place neighbours past the cube), so a finite
    one outside the cube is stored exactly as given, never moved."""
    sector = SpaceSector("outside", edge_ly=10.0)
    entry = sector.add_system(_StubSystem(1.0), position=position)
    assert entry.position == position


@given(edge=st.one_of(st.sampled_from([math.inf, -math.inf, math.nan, 0.0, -0.0]), st.floats(max_value=0.0)))
@example(edge=math.inf)
@example(edge=math.nan)
@example(edge=-10.0)
def test_non_finite_or_non_positive_edge_is_rejected(edge):
    """Regression: edge inf/nan placed systems at NaN, a negative edge put
    them outside the sector's own cube."""
    with pytest.raises(ValueError):
        SpaceSector("bad edge", edge_ly=edge)


# ---------------------------------------------------------------------------
# Galaxy cells: sampling at the faces
# ---------------------------------------------------------------------------

class _EndpointRng:
    """`uniform(a, b)` returns exactly `a` or `b` per a drawn bit pattern --
    the faces a real RNG reaches only by rounding."""

    def __init__(self, bits):
        self.bits = list(bits)

    def uniform(self, a, b):
        return b if self.bits.pop(0) else a


@given(ring=RINGS, edge=st.floats(1e-3, 1e4), bits=st.lists(st.booleans(), min_size=3, max_size=3))
def test_cell_samples_on_every_face_are_inside(ring, edge, bits):
    cell = SectorCell.for_ring(ring, edge)
    point = cell.sample(_EndpointRng(bits))
    assert all(math.isfinite(c) for c in point)
    assert cell.contains(point), (ring, edge, bits, point)


@given(ring=RINGS, edge=st.floats(1e-3, 1e4), seed=SEEDS)
def test_cell_samples_are_inside_and_volume_matches(ring, edge, seed):
    cell = SectorCell.for_ring(ring, edge)
    rng = random.Random(seed)
    for _ in range(20):
        assert cell.contains(cell.sample(rng))
    # Every ring's slots tile the annulus: n * slot volume == annulus volume.
    annulus = math.pi * (cell.r_outer ** 2 - cell.r_inner ** 2) * 2 * cell.half_height
    assert cell.volume * ring_sector_count(ring) == pytest.approx(annulus, rel=1e-9)


# ---------------------------------------------------------------------------
# Growth
# ---------------------------------------------------------------------------

@settings(max_examples=scaled(80))
@given(
    seed=SEEDS,
    sector=sectors(),
    home_hill=HILL_LY.filter(lambda h: h <= 20.0),
    hills=st.lists(st.one_of(st.sampled_from([0.3, 0.5, 1.71, 7.5, 1e3]), st.floats(0.3, 10.0)),
                   min_size=1, max_size=6),
    target=st.one_of(st.none(), st.sampled_from([0, 1, 2, 5]), st.integers(0, 60)),
    k=st.sampled_from([0, 1, 2, 5, prog.SECTOR_GROWTH_POISSON_DISK_K]),
)
def test_grow_from_seed_terminates_inside_and_separated(seed, sector, home_hill, hills, target, k):
    # grow_from_seed is O(n^2 * k) in the sector's final size; keep the
    # Poisson-drawn default target (target=None) to sectors that expect a
    # realistic handful (a huge sector's own default is ~740, see the
    # capped-Poisson xfail below) so one example stays well under a second.
    assume(target is not None or sector.expected_system_count() <= 30)
    note(f"edge={sector.edge_ly} cell={vars(sector.cell) if sector.cell else None}")
    factory_calls = []

    def factory():
        factory_calls.append(1)
        return _StubSystem(hills[len(factory_calls) % len(hills)])

    with _seeded_sector_rng(seed), _time_limit(20):
        home = sector.add_home_system(_StubSystem(home_hill))
        before = list(sector.entries)
        new = sector.grow_from_seed(home, factory, target_count=target, k=k)
    assert sector.entries == before + new
    assert_sector_sane(sector)
    if target is not None:
        assert len(sector.entries) <= max(target, 1)
    if k <= 0 or (target is not None and target <= 1):
        assert new == []
    # Every accepted candidate cost one factory call; retired parents cost k.
    assert len(new) <= len(factory_calls)


@settings(max_examples=scaled(10))
@given(seed=SEEDS, hill=st.floats(2.0, 8.0))
def test_grow_toward_an_unreachable_target_stops_when_full(seed, hill):
    sector = SpaceSector("greedy")
    with _seeded_sector_rng(seed), _time_limit(30):
        home = sector.add_home_system(_StubSystem(hill))
        sector.grow_from_seed(home, lambda: _StubSystem(hill), target_count=10**6)
    assert 1 <= len(sector.entries) < 200
    assert_sector_sane(sector)


def test_grow_from_seed_rejects_a_foreign_seed():
    sector = SpaceSector("a")
    other = SpaceSector("b")
    foreign = other.add_home_system(_StubSystem(1.0))
    with pytest.raises(ValueError):
        sector.grow_from_seed(foreign, lambda: _StubSystem(1.0), target_count=5)


def test_growth_respects_massive_phenomena():
    sector = SpaceSector("black hole nearby")
    with _seeded_sector_rng(1), _time_limit(20):
        home = sector.add_system(_StubSystem(0.5), position=(0.0, 0.0, 0.0))
        sector.add_phenomenon(_StubRemnant(3.0, "bh"), "black-hole", position=(4.0, 0.0, 0.0))
        sector.grow_from_seed(home, lambda: _StubSystem(0.5), target_count=40)
    assert_sector_sane(sector)


# ---------------------------------------------------------------------------
# Poisson count sampler
# ---------------------------------------------------------------------------

@given(mean=st.one_of(st.floats(max_value=0.0, allow_nan=False, allow_infinity=False), st.sampled_from([0.0, -0.0])))
def test_poisson_non_positive_mean_is_zero(mean):
    assert ss._sample_poisson_count(mean, random.Random(0)) == 0


@settings(max_examples=scaled(30))
@given(mean=st.floats(min_value=1e-9, max_value=50.0), seed=SEEDS)
def test_poisson_small_means_are_unbiased(mean, seed):
    rng = random.Random(seed)
    n = 300
    with _time_limit(10):
        samples = [ss._sample_poisson_count(mean, rng) for _ in range(n)]
    assert all(isinstance(x, int) and x >= 0 for x in samples)
    average = sum(samples) / n
    assert abs(average - mean) <= 6 * math.sqrt(mean / n) + 1e-9


@given(mean=st.sampled_from([math.nan, math.inf, -math.inf]))
@example(mean=math.nan)
def test_poisson_non_finite_mean_is_a_value_error(mean):
    """Regression: a NaN mean never returned (reachable from
    'generate.py sector --density nan')."""
    with _time_limit(3), pytest.raises(ValueError):
        ss._sample_poisson_count(mean, random.Random(0))


@settings(max_examples=scaled(30))
@given(mean=st.floats(min_value=100.0, max_value=1e15), seed=SEEDS)
@example(mean=5000.0, seed=1)  # was capped at 740
@example(mean=1000.0, seed=0)
@example(mean=1e9, seed=0)
def test_poisson_large_means_are_not_capped(mean, seed):
    rng = random.Random(seed)
    n = 20
    with _time_limit(20):
        samples = [ss._sample_poisson_count(mean, rng) for _ in range(n)]
    assert all(isinstance(x, int) and x >= 0 for x in samples)
    assert abs(sum(samples) / n - mean) <= 6 * math.sqrt(mean / n) + 1


# ---------------------------------------------------------------------------
# nearest_neighbors
# ---------------------------------------------------------------------------

POSITIONS = st.lists(
    st.tuples(*[st.one_of(st.sampled_from([0.0, -0.0, 5.0, -5.0]), st.floats(-5.0, 5.0))] * 3),
    min_size=1, max_size=25,
)


@settings(max_examples=scaled(100))
@given(positions=POSITIONS, pick=st.integers(0, 24), count=st.integers(0, 30))
def test_nearest_neighbors_matches_brute_force(positions, pick, count):
    sector = SpaceSector("nn", edge_ly=10.0)
    entries = [sector.add_system(_StubSystem(0.0), position=p) for p in positions]
    entry = entries[pick % len(entries)]
    result = sector.nearest_neighbors(entry, count)
    others = [e for e in entries if e is not entry]
    assert len(result) == min(count, len(others))
    assert all(r is not entry for r in result)
    assert len({id(r) for r in result}) == len(result)
    expected = sorted(math.dist(entry.position, e.position) for e in others)[:count]
    assert [math.dist(entry.position, r.position) for r in result] == expected
    # Nothing left out is closer than the farthest thing returned.
    if result:
        cutoff = math.dist(entry.position, result[-1].position)
        left_out = [e for e in others if all(e is not r for r in result)]
        assert all(math.dist(entry.position, e.position) >= cutoff for e in left_out)


@given(count=st.integers(max_value=-1))
@example(count=-1)
@example(count=-5)
def test_nearest_neighbors_negative_count_is_a_value_error(count):
    """Regression: count=-1 returned others[:-1] (all but the farthest)."""
    sector = SpaceSector("nn", edge_ly=10.0)
    entries = [sector.add_system(_StubSystem(0.0), position=(float(i) / 10, 0.0, 0.0)) for i in range(10)]
    with pytest.raises(ValueError):
        sector.nearest_neighbors(entries[0], count)


# ---------------------------------------------------------------------------
# Round trip with real systems and phenomena
# ---------------------------------------------------------------------------

_POOL = {}


def _real_system(index):
    """A small, deterministic pool of real systems (single, close binary,
    wide binary, with and without planets), generated once per process."""
    if index not in _POOL:
        cfg = SystemConfig()
        cfg.BINARY_SYSTEM, cfg.WIDE_BINARY = [(False, None), (True, False), (True, True)][index % 3]
        cfg.PLANETS = index % 2 == 0 or None
        cfg.MOONS = True if index % 4 == 0 else None
        cfg.COMETS = index % 5 == 0
        with _seeded_generation(index):
            _POOL[index] = (StarSystem(system_config=cfg), cfg)
    return _POOL[index]


def _real_phenomenon(kind, seed):
    with _seeded_generation(seed):
        return ss._PHENOMENON_CLASSES_BY_TYPE[kind](SystemConfig())


def _canonical(data):
    return json.dumps(data, sort_keys=True)


@settings(max_examples=scaled(25), suppress_health_check=[HealthCheck.too_slow])
@given(
    seed=SEEDS,
    name=hostile_text,
    system_indices=st.lists(st.integers(0, 11), max_size=5),
    phenomena=st.lists(st.sampled_from(sorted(ss._PHENOMENON_CLASSES_BY_TYPE)), max_size=4),
    edge=st.sampled_from([DEFAULT_EDGE, 30.0, 200.0]),
)
def test_sector_round_trip_is_lossless(tmp_path_factory, seed, name, system_indices, phenomena, edge):
    sector = SpaceSector(name, edge_ly=edge)
    with _seeded_sector_rng(seed):
        for index in system_indices:
            system, cfg = _real_system(index)
            try:
                sector.add_system(system, system_config=cfg)
            except ValueError:
                pass
        for i, kind in enumerate(phenomena):
            try:
                sector.add_phenomenon(_real_phenomenon(kind, seed + i), kind)
            except ValueError:
                pass
    assert_sector_sane(sector)
    before = sector.to_dict()
    strict = json.dumps(before, allow_nan=False)  # also: every number finite
    rebuilt = SpaceSector.from_dict(json.loads(strict))
    after = rebuilt.to_dict()
    assert _canonical(after) == _canonical(before)
    assert SpaceSector.from_dict(after).to_dict() == after
    assert [e.position for e in rebuilt.entries] == [e.position for e in sector.entries]
    assert [e.phenomenon_type for e in rebuilt.phenomena] == [e.phenomenon_type for e in sector.phenomena]
    for a, b in zip(sector.entries, rebuilt.entries):
        assert a.named_location() == b.named_location()
        assert [n.position for n in sector.nearest_neighbors(a, 3)] == \
               [n.position for n in rebuilt.nearest_neighbors(b, 3)]
    path = tmp_path_factory.mktemp("sector") / "sector.json"
    sector.save(str(path))
    assert _canonical(SpaceSector.load(str(path)).to_dict()) == _canonical(before)


def test_from_dict_without_generated_key_rebuilds_from_the_recipe():
    system, cfg = _real_system(1)
    sector = SpaceSector("recipe")
    sector.add_system(system, system_config=cfg, position=(1.0, 2.0, 3.0))
    data = sector.to_dict()
    del data["systems"][0]["generated"]
    with _seeded_generation(3):
        rebuilt = SpaceSector.from_dict(data)
    assert rebuilt.entries[0].position == (1.0, 2.0, 3.0)
    assert rebuilt.entries[0].star_system.star.name == system.star.name  # name pinned in the recipe


@given(ring=st.integers(0, 10**6), edge=st.floats(1e-3, 1e4))
@example(ring=5, edge=DEFAULT_EDGE)
def test_round_trip_keeps_the_galaxy_cell(ring, edge):
    cell = SectorCell.for_ring(ring, edge)
    sector = SpaceSector("celled", edge_ly=edge, cell=cell)
    rebuilt = SpaceSector.from_dict(json.loads(json.dumps(sector.to_dict())))
    assert rebuilt.volume_ly3 == pytest.approx(sector.volume_ly3)
    assert rebuilt.expected_system_count() == pytest.approx(sector.expected_system_count())
    for point in [(0.0, 0.0, 0.0), (edge * 0.49, 0.0, 0.0), (0.0, 0.0, edge * 0.51)]:
        assert rebuilt.contains(point) == sector.contains(point)


@pytest.mark.parametrize("ring", [0, 5, 40])
def test_db_load_rebuilds_the_galaxy_cell_from_the_ring_address(mysql_config, ring):
    """The `sectors` table has no cell column; `_db.load_sector` must
    rebuild it from the stored ring address and edge (a galaxy-placed
    sector used to come back as a cube)."""
    from stellarObjects import _db
    edge = DEFAULT_EDGE
    sector = SpaceSector("celled db", edge_ly=edge, cell=SectorCell.for_ring(ring, edge))
    placement = {"center_x_pc": 1.0, "center_y_pc": 0.0, "center_z_pc": 0.0, "galactic_radius_pc": 1.0,
                 "ring_index": ring, "layer_index": 0, "ring_slot_index": 0}
    sector_id = _db.save_sector(sector, config=mysql_config, galaxy_position=placement)
    conn = _db.get_connection(mysql_config)
    try:
        loaded = _db.load_sector(conn, sector_id)
    finally:
        conn.close()
    assert loaded.cell is not None
    assert loaded.volume_ly3 == pytest.approx(sector.volume_ly3, rel=1e-9)
    assert loaded.expected_system_count() == pytest.approx(sector.expected_system_count(), rel=1e-9)
