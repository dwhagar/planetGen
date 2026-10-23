"""
GEN.329: runaway bookkeeping. An O or a B star is a runaway at the observed rates (20% and 7.5%, within the catalogue's
10% to 30% and 5% to 10%), by the supernova channel for two thirds of them and the dynamical one for the rest; a supernova
runaway keeps a compact companion at the share the collapse leaves bound; a hypervelocity star stays single; every
runaway carries its origin, in the tree of its scatter row, on the system and in the flags of its `stars` row.
"""

import math
import random

import pytest

from planetgen.db import store
from planetgen.generation import multistar, run_galaxy, run_plan, run_sector, scatter_companions
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena import system_types
from planetgen.generation.system import StarSystem
from planetgen.physics import constants
from planetgen.util import draw
from tests.test_bright_star_scatter import EDGE_PC, _scatter_args, _seed_galaxy, light_steps  # noqa: F401
from tests.test_gen_scatter_companions import STAR, _star_row

RATES = system_types.SYSTEM_TYPE_RATES


def _sigma(share, count):
    return math.sqrt(share * (1.0 - share) / count)


def _units(count, seed=1):
    rng = random.Random(seed)
    return [(rng.random(), rng.random(), rng.random()) for _ in range(count)]


# --- The draw -----------------------------------------------------------------------------------------------------


def test_only_o_and_b_stars_have_an_observed_rate():
    assert multistar.runaway_rate(30.0) == RATES["runaway-star"].default == 0.20
    assert multistar.runaway_rate(5.0) == RATES["runaway-star"].by_class["B"] == 0.075
    assert RATES["runaway-star"].low <= multistar.runaway_rate(30.0) <= RATES["runaway-star"].high
    for mass in (0.4, 1.0, 1.5, 2.0):
        assert multistar.runaway_rate(mass) is None
        assert multistar.draw_runaway_origin(mass, 0.0, 0.0, 0.0) is None


def test_the_channel_split_is_two_thirds_supernova_and_the_companion_share_is_what_the_collapse_leaves_bound():
    assert multistar.runaway_origin_at(0.0) == "supernova" and multistar.runaway_origin_at(0.66) == "supernova"
    assert multistar.runaway_origin_at(0.67) == "dynamical" and multistar.runaway_origin_at(0.999) == "dynamical"
    assert multistar.bound_companion_share() == pytest.approx(1.0 - 0.83)


@pytest.mark.parametrize("mass, rate", [(30.0, 0.20), (5.0, 0.075)])
def test_runaway_counts_are_the_observed_rates_within_four_sigma_and_the_catalogues_ranges(mass, rate):
    count = 40000
    draws = [multistar.draw_runaway_origin(mass, *units) for units in _units(count)]
    found = [d for d in draws if d is not None]
    assert abs(len(found) / count - rate) < 4.0 * _sigma(rate, count)
    low, high = (0.10, 0.30) if mass > 16.0 else (0.05, 0.10)
    assert low <= len(found) / count <= high
    supernova = [d for d in found if d.origin == "supernova"]
    assert abs(len(supernova) / len(found) - 2 / 3) < 4.0 * _sigma(2 / 3, len(found))
    assert all(d.origin in ("supernova", "dynamical") for d in found)     # no runaway has a missing origin
    bound = [d for d in found if d.bound_companion]
    assert all(d.origin == "supernova" for d in bound)
    share = multistar.bound_companion_share()
    assert abs(len(bound) / len(supernova) - share) < 4.0 * _sigma(share, len(supernova))


def test_the_scatter_derives_the_draw_from_the_rows_seed():
    seeds = random.Random(5).sample(range(1, 10 ** 9), 4000)
    draws = [scatter_companions.runaway_draw(seed, 25.0) for seed in seeds]
    assert draws == [scatter_companions.runaway_draw(seed, 25.0) for seed in seeds]
    share = sum(d is not None for d in draws) / len(seeds)
    assert abs(share - 0.20) < 4.0 * _sigma(0.20, len(seeds))
    assert scatter_companions.runaway_draw(seeds[0], 1.0) is None


# --- The tree -----------------------------------------------------------------------------------------------------


def _hit(origin, bound, start=1, mass=20.0):
    """A seed whose draw for a star of `mass` solar masses is `origin` with `bound` compact companion."""
    rng = random.Random(start)
    while True:
        seed = rng.getrandbits(63)
        draw_ = scatter_companions.runaway_draw(seed, mass)
        if draw_ is not None and draw_.origin == origin and draw_.bound_companion == bound:
            return seed


def _row(seed, mass=20.0, age_gy=0.005):
    return _star_row(seed, initial_mass=mass, age_gy=age_gy)


@pytest.mark.parametrize("origin", ["supernova", "dynamical"])
def test_a_runaways_tree_carries_the_flag_of_its_channel_on_member_zero(origin):
    seeds = []
    start = 1
    while len(seeds) < 6:
        seed = _hit(origin, False, start)
        start += 1
        if scatter_companions.star_tree(_row(seed)) is not None:
            seeds.append(seed)
    for seed in seeds:
        tree = scatter_companions.star_tree(_row(seed))
        assert tree.members[0].flags & multistar.RUNAWAY_ORIGIN_FLAGS[origin]
        other = "dynamical" if origin == "supernova" else "supernova"
        assert not tree.members[0].flags & multistar.RUNAWAY_ORIGIN_FLAGS[other]


def test_a_supernova_runaway_that_kept_its_companion_has_a_neutron_star_or_black_hole_in_its_tree():
    # a 12 solar mass star at 15 Myr: heavier stars have had time to explode before it. (At 5 Myr a 20 solar mass
    # runaway finds a progenitor in about a quarter of the draws: only stars over 50 solar masses have exploded by then.)
    found = kept = 0
    rng = random.Random(9)
    while found < 20:
        seed = _hit("supernova", True, rng.getrandbits(32), mass=12.0)
        found += 1
        tree = scatter_companions.star_tree(_row(seed, mass=12.0, age_gy=0.015))
        assert tree is not None
        compact = [m for m in tree.members[1:] if m.kind in multistar.COMPACT_KINDS]
        kept += bool(compact)
        assert tree.members[0].flags & multistar.MEMBER_RUNAWAY_SUPERNOVA
    assert kept == found


def test_a_star_that_is_not_a_runaway_has_no_origin_bit_and_its_tree_is_the_one_it_always_drew():
    rng = random.Random(3)
    checked = 0
    while checked < 30:
        seed = rng.getrandbits(63)
        if scatter_companions.runaway_draw(seed, 20.0) is not None:
            continue
        tree = scatter_companions.star_tree(_row(seed))
        checked += 1
        if tree is not None:
            assert not any(m.flags & (multistar.MEMBER_RUNAWAY_SUPERNOVA | multistar.MEMBER_RUNAWAY_DYNAMICAL)
                           for m in tree.members)


def test_a_white_dwarf_row_is_never_a_runaway_by_the_catalogue():
    row = list(_row(123, mass=20.0))
    row[STAR["yerkes_class"]] = "VII"
    tree = scatter_companions.star_tree(tuple(row))
    assert tree is None or not tree.members[0].flags & (multistar.MEMBER_RUNAWAY_SUPERNOVA
                                                        | multistar.MEMBER_RUNAWAY_DYNAMICAL)


def test_hypervelocity_stars_stay_single():
    for subtype in (None, "star", "hypervelocity"):
        assert scatter_companions.phenomenon_companion("hypervelocity-star", subtype, 77) is None


# --- The fill -----------------------------------------------------------------------------------------------------


class _Entry:
    def __init__(self, system):
        self.star_system = system


class _Sector:
    def __init__(self, systems):
        self.entries = [_Entry(system) for system in systems]


def _system(seed):
    with draw.bound(seed):
        return StarSystem(SystemConfig())


def test_the_field_roll_gives_every_runaway_an_origin(monkeypatch):
    monkeypatch.setattr(run_sector.program_constants, "phenomenon_rate_per_star", lambda kind: 0.5)
    systems = [_system(seed) for seed in range(12)]
    with draw.bound(4):
        flagged = run_sector.flag_fast_stars(_Sector(systems), galactic_center_dist_ly=26000.0, scattered=True)
    runaways = [s for s in systems if s.runaway_class == "runaway"]
    assert flagged == len(runaways) > 3
    assert all(s.runaway_origin in ("supernova", "dynamical") for s in runaways)
    assert all(s.runaway_origin is None for s in systems if s.runaway_class != "runaway")


def test_a_preplaced_o_or_b_star_is_a_runaway_by_its_seed_and_the_field_roll_leaves_it_alone():
    seed = _hit("supernova", False)
    row = {"seed": seed, "initial_mass_sol": 20.0, "yerkes_class": "V"}
    first, second = _system(1), _system(2)
    run_sector._flag_scattered_runaway(first, row)
    run_sector._flag_scattered_runaway(second, row)
    assert first.runaway_class == "runaway" and first.runaway_origin == "supernova" and first.runaway_decided
    assert (first.runaway_speed_kms, first.runaway_direction) == (second.runaway_speed_kms, second.runaway_direction)
    low, high = run_sector.program_constants.RUNAWAY_STAR_SPEED_RANGE_KMS
    assert low <= first.runaway_speed_kms <= high
    before = (first.runaway_class, first.runaway_speed_kms, first.runaway_origin)
    with draw.bound(8):
        run_sector.flag_fast_stars(_Sector([first]), galactic_center_dist_ly=26000.0, scattered=False)
    assert (first.runaway_class, first.runaway_speed_kms, first.runaway_origin) == before


def test_a_preplaced_star_that_is_not_a_runaway_is_decided_too_and_a_white_dwarf_is_left_to_the_field():
    rng = random.Random(2)
    seed = next(s for s in (rng.getrandbits(63) for _ in range(100)) if scatter_companions.runaway_draw(s, 20.0) is None)
    plain = _system(3)
    run_sector._flag_scattered_runaway(plain, {"seed": seed, "initial_mass_sol": 20.0, "yerkes_class": "V"})
    assert plain.runaway_class is None and plain.runaway_decided
    dwarf = _system(4)
    run_sector._flag_scattered_runaway(dwarf, {"seed": seed, "initial_mass_sol": 20.0, "yerkes_class": "VII"})
    assert not dwarf.runaway_decided and dwarf.runaway_class is None


def test_the_origin_survives_the_json_form():
    system = _system(6)
    system.runaway_class, system.runaway_speed_kms, system.runaway_origin = "runaway", 60.0, "dynamical"
    again = StarSystem.from_dict(system.to_dict())
    assert again.runaway_origin == "dynamical"
    legacy = system.to_dict()
    legacy.pop("runaway_origin")
    assert StarSystem.from_dict(legacy).runaway_origin is None


# --- The database -------------------------------------------------------------------------------------------------


def test_a_saved_runaway_keeps_its_origin_in_its_star_rows_flags(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        ids = {}
        for origin in ("supernova", "dynamical", None):
            system = _system(31)
            if origin is not None:
                system.runaway_class, system.runaway_speed_kms, system.runaway_origin = "runaway", 55.0, origin
            with conn:
                system.star.name = f"Runaway {origin}"
                ids[origin] = store.insert_star_system(conn, system, SystemConfig())
        for origin, system_id in ids.items():
            flags = conn.execute("SELECT flags FROM stars WHERE star_system_id = ?", (system_id,)).fetchone()["flags"]
            assert flags == (0 if origin is None else multistar.RUNAWAY_ORIGIN_FLAGS[origin])
            assert store.load_star_system(conn, system_id).runaway_origin == origin
    finally:
        conn.close()


def test_every_runaway_a_fill_builds_carries_an_origin(mysql_config):
    _seed_galaxy(mysql_config)
    run_plan.scatter_bright_stars(_scatter_args(mysql_config, "--mass"))
    conn = store.get_connection(mysql_config)
    try:
        rows = conn.execute("SELECT shell_index, band_index, wedge_index, seed, initial_mass_sol, yerkes_class"
                            " FROM bright_stars WHERE star_system_id IS NULL").fetchall()
    finally:
        conn.close()
    hits = [r for r in rows if r["yerkes_class"] != "VII" and scatter_companions.runaway_draw(
        r["seed"], r["initial_mass_sol"]) is not None]
    assert hits, "the toy galaxy scattered no O or B runaway"
    address = (hits[0]["shell_index"], hits[0]["band_index"], hits[0]["wedge_index"])
    args = run_galaxy._default_generation_args(config=mysql_config)
    args.num_systems = 0
    from planetgen.galaxy.spherical import sector_position_pc
    sector_id, _name, _sector = run_galaxy.generate_and_save_sector_at(
        args, address, sector_position_pc(*address, EDGE_PC), EDGE_PC)
    conn = store.get_connection(mysql_config)
    try:
        saved = conn.execute("SELECT id, runaway_class FROM star_systems WHERE sector_id = ?", (sector_id,)).fetchall()
        runaways = [r["id"] for r in saved if r["runaway_class"] == "runaway"]
        assert len(runaways) >= 1
        for system_id in runaways:
            assert store.load_star_system(conn, system_id).runaway_origin in ("supernova", "dynamical")
    finally:
        conn.close()


def test_the_generator_the_catalogue_names_for_runaways_exists():
    for key in ("runaway-star", "unbound-by-collapse"):
        assert RATES[key].generator == "draw_runaway_origin" and hasattr(multistar, RATES[key].generator)
