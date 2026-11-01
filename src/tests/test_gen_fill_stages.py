"""
GEN.251: a sector fill as stages on their own random streams (`planetgen.generation.fill_stages`), and the system pass
(`planetgen.generation.system_pass`) that gives a scattered or below-cut compact remnant its planets at fill, on its
own seed, with the odds of `tuning.SYSTEM_PASS_ODDS`.
"""

import math

import pytest

from planetgen import tuning
from planetgen.galaxy.profile import profile_at
from planetgen.galaxy.sector import SpaceSector
from planetgen.galaxy.spherical import sector_position_pc
from planetgen.generation import bright_stars as brightStars
from planetgen.generation import fill_stages, run_galaxy, run_plan, run_sector, system_pass
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena import hosts, registry
from planetgen.generation.phenomena.compact_remnant import BlackHole, NeutronStar
from planetgen.generation.phenomena.rogue import InterstellarComet, RoguePlanet
from planetgen.physics.units import ly_to_pc
from planetgen.util import draw
from tests.test_host_tables import reweighted
from tests.toy_shape import toy_shape


def _sector(edge_ly=100.0):
    return SpaceSector("Stage Test", edge_ly=edge_ly)


def _sample(sector, name, sub=None, n=3):
    with fill_stages.stage(sector, name, sub):
        return [draw.random() for _ in range(n)]


def _begun(seed):
    with draw.bound(seed):
        sector = _sector()
        fill_stages.begin(sector)
    return sector


# --- The streams --------------------------------------------------------------------------------------------------


def test_a_fill_draws_its_base_from_the_stream_it_is_bound_to():
    assert _begun(5).stage_base == _begun(5).stage_base
    assert _begun(5).stage_base != _begun(6).stage_base
    assert 0 <= _begun(5).stage_base < 1 << fill_stages.BASE_BITS


def test_each_stage_and_each_stream_within_a_stage_draws_from_a_stream_of_its_own():
    a, b = _begun(5), _begun(5)
    assert _sample(a, fill_stages.SYSTEMS) == _sample(b, fill_stages.SYSTEMS)
    draws = {_sample(a, name)[0] for name in fill_stages.STAGES}
    assert len(draws) == len(fill_stages.STAGES)
    assert _sample(a, fill_stages.PER_STAR_KINDS, "comet") != _sample(a, fill_stages.PER_STAR_KINDS, "rogue-planet")
    assert _sample(a, fill_stages.FLAGS, "nucleus") != _sample(a, fill_stages.FLAGS, "fast-stars")


def test_drawing_in_one_stage_moves_no_draw_of_another():
    a, b = _begun(11), _begun(11)
    with fill_stages.stage(a, fill_stages.SCATTERED):
        for _ in range(500):
            draw.random()
    assert _sample(a, fill_stages.SYSTEMS, n=5) == _sample(b, fill_stages.SYSTEMS, n=5)
    assert _sample(a, fill_stages.FLAGS) == _sample(b, fill_stages.FLAGS)


def test_a_sector_no_fill_began_on_keeps_the_stream_it_has():
    sector = _sector()
    assert fill_stages.stream_seed(sector, fill_stages.SYSTEMS) is None
    with draw.bound(3):
        before = draw.random()
    with draw.bound(3):
        assert _sample(sector, fill_stages.SYSTEMS, n=1) == [before]


def test_the_counts_are_kept_by_stage_and_kind_and_described_in_fill_order():
    sector = _begun(1)
    fill_stages.count(sector, fill_stages.SYSTEMS, "star-system", 3)
    fill_stages.count(sector, fill_stages.SCATTERED, "black-hole")
    fill_stages.count(sector, fill_stages.SCATTERED, "black-hole")
    fill_stages.count(sector, fill_stages.SCATTERED, "neutron-star")
    fill_stages.count(sector, fill_stages.FLAGS, "fast-star", 0)
    assert fill_stages.counts(sector) == {fill_stages.SCATTERED: {"black-hole": 2, "neutron-star": 1},
                                          fill_stages.SYSTEMS: {"star-system": 3}}
    assert fill_stages.describe(sector) == ["scattered: 2 black-hole, 1 neutron-star", "systems: 3 star-system"]
    fill_stages.begin(sector)
    assert fill_stages.counts(sector) == {}


# --- A generated sector --------------------------------------------------------------------------------------------


def _args():
    args = run_galaxy._default_generation_args()
    args.num_systems = 3
    args.density = None
    return args


def _everything(sector):
    systems = {(entry.star_system.name, tuple(round(c, 9) for c in entry.position)) for entry in sector.entries}
    phenomena = {(entry.phenomenon_type, getattr(entry.phenomenon, "name", None),
                  tuple(round(c, 9) for c in entry.position)) for entry in sector.phenomena}
    return systems, phenomena


def _generate(seed):
    with draw.bound(seed):
        return run_sector.generate_sector(_args())[1]


def test_a_sector_is_the_same_from_the_same_seed():
    assert _everything(_generate(21)) == _everything(_generate(21))
    assert _everything(_generate(21)) != _everything(_generate(22))


def test_a_kind_added_to_the_per_star_rolls_moves_no_other_draw(monkeypatch):
    real_rate = tuning.phenomenon_rate_per_star

    def extra_factory(config, _distance):
        comet = InterstellarComet(config)
        comet.name = "Extra Kind"
        return comet

    baseline = _generate(33)
    monkeypatch.setattr(tuning, "phenomenon_rate_per_star", lambda kind: 4.0 if kind == "extra-kind" else real_rate(kind))
    monkeypatch.setattr(run_sector, "SECTOR_PHENOMENON_KINDS",
                        (("extra-kind", "comet", extra_factory),) + run_sector.SECTOR_PHENOMENON_KINDS)
    changed = _generate(33)
    extras = [entry for entry in changed.phenomena if getattr(entry.phenomenon, "name", None) == "Extra Kind"]
    assert len(extras) > 1
    systems, phenomena = _everything(changed)
    assert systems == _everything(baseline)[0]
    assert {item for item in phenomena if item[1] != "Extra Kind"} == _everything(baseline)[1]
    assert changed.stage_counts[(fill_stages.PER_STAR_KINDS, "extra-kind")] == len(extras)


def test_the_stages_count_what_they_made():
    sector = _generate(33)
    made = fill_stages.counts(sector)
    assert made[fill_stages.SYSTEMS]["star-system"] == len(sector.entries)
    assert sum(n for kind, n in made.get(fill_stages.PER_STAR_KINDS, {}).items() if kind != "star-system") \
        + sum(made.get(fill_stages.NEBULAE, {}).values()) == len(sector.phenomena)
    assert all(stage in fill_stages.STAGES for stage in made)


# --- The system pass -----------------------------------------------------------------------------------------------


def _millisecond_pulsar(seed):
    """A neutron star built from `seed` that came out a millisecond pulsar (the odds of that are 0.2%, so the recycled
    chance is forced), or `None`."""
    config = SystemConfig()
    with draw.bound(seed):
        star = NeutronStar(config, age_gy=5.0)
    return star if star.pulsar_type == "millisecond" else None


@pytest.fixture
def recycled(monkeypatch):
    monkeypatch.setattr(tuning, "PULSAR_RECYCLED_CHANCE", 1.0)


def _reweigh(monkeypatch, table, row, weight):
    """Gives the row of a host table another weight for one test."""
    monkeypatch.setitem(hosts.REGISTRY_TABLES, table, reweighted(hosts.REGISTRY_TABLES[table], row, weight))


@pytest.fixture
def certain(monkeypatch, recycled):
    _reweigh(monkeypatch, "remnant-satellites", "planets", 1.0)


def _with_pulsar(seed, system_seed=None, scatter_id=None, sector=None):
    sector = sector or _begun(seed)
    entry = sector.add_phenomenon(_millisecond_pulsar(seed), "neutron-star", position=(0.0, 0.0, 0.0))
    entry.system_seed = seed if system_seed is None else system_seed
    if scatter_id is not None:
        entry.scatter_id = scatter_id
    return sector, entry


def test_an_object_key_names_the_kind_and_for_a_remnant_its_type(recycled):
    config = SystemConfig()
    sector = _sector()
    neutron_star = sector.add_phenomenon(_millisecond_pulsar(2), "neutron-star", position=(0.0, 0.0, 0.0))
    assert system_pass.object_key(neutron_star) == "neutron-star:millisecond"
    with draw.bound(4):
        black_hole = BlackHole(config, mass_class="stellar")
    assert system_pass.object_key(sector.add_phenomenon(black_hole, "black-hole", position=(30.0, 0.0, 0.0))) \
        == "black-hole:stellar"
    assert system_pass.object_key(sector.add_phenomenon(InterstellarComet(config), "comet")) == "comet"


def test_the_odds_are_the_rows_of_the_host_tables():
    for table_name in ("remnant-satellites", "planet-satellites"):
        table = hosts.REGISTRY_TABLES[table_name]
        assert table.row(hosts.NO_SATELLITES)
        for row in table.rows:
            assert row.state in ("partly verified", "stand-in") and row.source
    assert hosts.chance("remnant-satellites", "neutron-star:millisecond") == pytest.approx(0.007)
    assert hosts.chance("planet-satellites", "rogue-planet") == pytest.approx(0.2)
    assert registry.BY_KIND["neutron-star"].host == "remnant" and registry.BY_KIND["rogue-planet"].host == "planet"
    for host in hosts.HOSTS:
        assert all(name is None or name in hosts.REGISTRY_TABLES for name in (host.companions, host.satellites))


def test_the_share_of_millisecond_pulsars_that_roll_planets_is_the_table_within_poisson_bounds():
    key, seeds = "neutron-star:millisecond", 100_000
    expected = seeds * hosts.chance("remnant-satellites", key)
    hits = sum(system_pass.hit(key, seed) for seed in range(seeds))
    assert abs(hits - expected) < 4 * math.sqrt(expected), (hits, expected)


def test_a_kind_without_odds_never_rolls_planets():
    for key in ("neutron-star:magnetar", "neutron-star:young", "black-hole:stellar", "black-hole:intermediate",
                "black-hole:supermassive", "comet", "asteroid-field", "brown-dwarf"):
        assert hosts.chance("remnant-satellites", key) == 0.0
        assert not any(system_pass.hit(key, seed) for seed in range(2000))


def test_a_pulsar_that_rolls_planets_becomes_a_system_at_its_own_point(certain):
    made = None
    for seed in range(1, 16):
        sector, entry = _with_pulsar(seed, scatter_id=77)
        star = entry.phenomenon
        position = entry.position
        systems = system_pass.run(sector, _args(), [entry])
        key = "neutron-star:millisecond"
        assert sector.stage_counts[(fill_stages.SYSTEM_PASS, key)] == 1
        if not systems:
            assert sector.stage_counts[(fill_stages.SYSTEM_PASS, key + ":empty")] == 1
            assert sector.phenomena == [entry] and sector.entries == []
            continue
        (made,) = systems
        assert sector.phenomena == [] and sector.entries == [made]
        assert made.star_system.star is star and made.star_system.planets
        assert made.position == position and made.preplaced and made.scatter_id == 77
        assert sector.stage_counts[(fill_stages.SYSTEM_PASS, key + ":system")] == 1
        break
    assert made is not None


def test_the_system_is_the_same_from_the_same_seed_and_in_any_order(certain):
    def planets(order):
        sector = _begun(1)
        entries = [sector.add_phenomenon(_millisecond_pulsar(s), "neutron-star", position=(40.0 * i, 0.0, 0.0))
                   for i, s in enumerate(order)]
        for entry, s in zip(entries, order):
            entry.system_seed = s
        system_pass.run(sector, _args(), entries)
        return {entry.star_system.star.mass: [(type(planet).__name__, round(planet.distance, 9))
                                              for planet in entry.star_system.planets]
                for entry in sector.entries}

    # Only a disc origin (a quarter of the pulsars that roll planets, GEN.317) leaves a system: take two seeds that do.
    first, second = [seed for seed in range(1, 60) if planets([seed])][:2]
    forward, backward = planets([first, second]), planets([second, first])
    assert forward == backward and len(forward) == 2


def test_a_roll_that_misses_leaves_the_bare_object(recycled, monkeypatch):
    _reweigh(monkeypatch, "remnant-satellites", "planets", 0.0)
    sector, entry = _with_pulsar(9)
    assert system_pass.run(sector, _args(), [entry]) == []
    assert sector.phenomena == [entry] and sector.entries == []
    assert sector.stage_counts == {(fill_stages.SYSTEM_PASS, "neutron-star:millisecond"): 1}


def test_a_black_hole_or_a_magnetar_is_not_rolled_at_all(monkeypatch):
    sector = _begun(2)
    config = SystemConfig()
    with draw.bound(8):
        hole = BlackHole(config, mass_class="stellar")
    entry = sector.add_phenomenon(hole, "black-hole", position=(0.0, 0.0, 0.0))
    entry.system_seed = 8
    assert system_pass.run(sector, _args(), [entry]) == []
    assert sector.stage_counts == {}


def test_a_remnant_without_a_seed_is_left_alone(certain):
    sector, entry = _with_pulsar(3)
    del entry.system_seed
    assert system_pass.run(sector, _args(), [entry]) == []
    assert sector.phenomena == [entry]


def test_compact_remnants_are_rolled_before_the_rest(certain, monkeypatch):
    sector, pulsar = _with_pulsar(3)
    comet = sector.add_phenomenon(InterstellarComet(SystemConfig()), "comet")
    comet.system_seed = 4
    rolled = []
    monkeypatch.setattr(system_pass, "object_key", lambda entry: rolled.append(entry.phenomenon_type) or
                        ("neutron-star:millisecond" if entry is pulsar else "comet"))
    system_pass.run(sector, _args(), [comet, pulsar])
    assert rolled == ["neutron-star", "comet"]


def test_replacing_a_phenomenon_keeps_the_neighbour_check_current():
    with draw.bound(1):
        sector = _sector(edge_ly=200.0)
        first = sector.add_phenomenon(NeutronStar(SystemConfig()), "neutron-star")
        second = sector.add_phenomenon(NeutronStar(SystemConfig()), "neutron-star")
        sector._neighbor_cache._refresh(sector)
        assert sector._neighbor_cache.count == 2
        system = system_pass.StarSystem(system_config=SystemConfig(), compact_remnant=first.phenomenon)
        replaced = sector.replace_phenomenon_with_system(first, system)
        assert sector.phenomena == [second] and sector.entries == [replaced]
        third = sector.add_phenomenon(NeutronStar(SystemConfig()), "neutron-star")
        sector._neighbor_cache._refresh(sector)
        assert sector._neighbor_cache.count == 3
        assert sector.phenomena == [second, third]
        with pytest.raises(ValueError):
            sector.replace_phenomenon_with_system(first, system)


def test_a_rogue_planet_keeps_its_moon_at_the_chance_of_the_table(monkeypatch):
    _reweigh(monkeypatch, "planet-satellites", "captured moon", 1.0)
    with draw.bound(1):
        assert all(RoguePlanet(SystemConfig()).has_moons for _ in range(5))
    _reweigh(monkeypatch, "planet-satellites", "captured moon", 0.0)
    with draw.bound(1):
        assert not any(RoguePlanet(SystemConfig()).has_moons for _ in range(5))


# --- A scattered fill ----------------------------------------------------------------------------------------------

EDGE_PC = ly_to_pc(tuning.DEFAULT_SECTOR_EDGE_LY)
SHAPE = toy_shape(disk_scale_length_pc=40.0, disk_scale_height_pc=12.0, bulge_scale_radius_pc=10.0,
                           bulge_amplitude=2.0, arm_count=2, pitch_angle_rad=math.radians(15), arm_amplitude=0.4)


def _row(row_id, seed, center, kind="neutron-star"):
    return {"id": row_id, "kind": kind, "subtype": None, "seed": seed, "mass_solar": 1.4,
            "position_x_mpc": round(center[0] * 1000), "position_y_mpc": round(center[1] * 1000),
            "position_z_mpc": round(center[2] * 1000)}


def test_a_filled_sector_gives_its_scattered_pulsar_planets_at_fill_from_the_row_seed(certain):
    center = sector_position_pc(2, 0, 3, EDGE_PC)
    profile = profile_at(center, SHAPE, EDGE_PC, 2000.0)
    seeds = []
    for seed in range(100, 140):
        with draw.bound(seed):
            probe = NeutronStar(SystemConfig(), mass_solar=1.4)
        if probe.pulsar_type == "millisecond":
            seeds.append(seed)
    assert seeds, "no millisecond pulsar among the probe seeds"

    def fill(order):
        rows = [_row(i + 1, s, center) for i, s in enumerate(order)]
        context = brightStars.FillContext(profile, phenomenon_rows=rows)
        args = _args()
        args.num_systems = 0
        with draw.bound(5):
            return run_sector.generate_sector(args, galactic_center_dist_ly=26000.0, fill=context)[1]

    sector = fill(seeds[:2])
    made = fill_stages.counts(sector)
    assert made[fill_stages.SCATTERED]["neutron-star"] == 2
    rolled = made[fill_stages.SYSTEM_PASS]
    assert rolled["neutron-star:millisecond"] == 2
    assert rolled.get("neutron-star:millisecond:system", 0) + rolled.get("neutron-star:millisecond:empty", 0) == 2
    systems = [entry for entry in sector.entries if isinstance(entry.star_system.star, NeutronStar)]
    bare = [entry for entry in sector.phenomena if entry.phenomenon_type == "neutron-star"]
    assert len(systems) == rolled.get("neutron-star:millisecond:system", 0)
    assert len(bare) == rolled.get("neutron-star:millisecond:empty", 0)
    assert {entry.scatter_id for entry in systems + bare} == {1, 2}
    assert all(entry.star_system.planets for entry in systems)
    if len(seeds) >= 2:
        again = fill(list(reversed(seeds[:2])))
        key = lambda s: sorted((entry.star_system.star.mass, len(entry.star_system.planets)) for entry in s.entries
                               if isinstance(entry.star_system.star, NeutronStar))
        assert key(again) == key(sector)


def test_a_saved_fill_stores_the_pulsar_systems_and_marks_every_scatter_row_built(mysql_config, certain, monkeypatch):
    from planetgen.db import store
    from tests import first_fill
    from tests.test_phenomenon_scatter import _busiest_cell, _scatter_args, _seed_galaxy
    first_fill.without_the_deepest_floor(monkeypatch)   # the toy galaxy's first fill would draw its thousands of stars (GEN.633, TEST.422)
    _seed_galaxy(mysql_config)
    run_plan.scatter_phenomena(_scatter_args(mysql_config))
    address = _busiest_cell(mysql_config, "kind = 'neutron-star'")
    conn = store.get_connection(mysql_config)
    try:
        rows = [row for row in store.phenomena_for_sector(conn, *address) if row["kind"] == "neutron-star"]
    finally:
        conn.close()
    args = run_galaxy._default_generation_args(config=mysql_config)
    args.num_systems = 0
    _id, _name, sector = run_galaxy.generate_and_save_sector_at(
        args, address, sector_position_pc(*address, ly_to_pc(tuning.DEFAULT_SECTOR_EDGE_LY)),
        ly_to_pc(tuning.DEFAULT_SECTOR_EDGE_LY))
    systems = [entry for entry in sector.entries if isinstance(entry.star_system.star, NeutronStar)]
    bare = [entry for entry in sector.phenomena if entry.phenomenon_type == "neutron-star"]
    assert len(systems) + len(bare) == len(rows)
    made = fill_stages.counts(sector).get(fill_stages.SYSTEM_PASS, {})
    assert made.get("neutron-star:millisecond", 0) <= len(rows)
    assert made.get("neutron-star:millisecond:system", 0) == len(systems)
    conn = store.get_connection(mysql_config)
    try:
        assert store.phenomena_for_sector(conn, *address) == []
        ids = [row["id"] for row in conn.execute("SELECT id FROM star_systems ORDER BY id").fetchall()]
        assert len(ids) == len(sector.entries)
        for system_id in ids:
            assert store.load_star_system(conn, system_id).stars
    finally:
        conn.close()
