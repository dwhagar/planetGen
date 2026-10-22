# tests/test_census_nebula_wholesphere.py

"""
TEST.184: the whole-sphere tests of the nebula and Wolf-Rayet side, on the census harness (`tests.audit`): a planetary
nebula always has its white dwarf; a core-collapse remnant has a compact core at the stated share; every hosted nebula
holds its host star; a Wolf-Rayet nebula holds a W star; tiling and the object ID round trip at the widest shell of the
reach; a profile and an empty sector are deterministic. TEST.185 holds the component side; TEST.159 is done when both are.

Three of the items were already tests of their own and are listed here so the umbrella can be read in one place; they run
in the normal suite:

- the Milky Way holds 1,200 +- 300 Wolf-Rayet stars and none is born in an old population:
  `test_wolf_rayet.py::test_the_milky_way_holds_about_twelve_hundred_wolf_rayet_stars` and
  `::test_none_is_born_old_enough_to_be_one_in_an_old_population`;
- the King solver against the textbook concentrations: `test_king.py::test_the_solver_reproduces_the_textbook_concentration`;
- a Wolf-Rayet star lights a ring in 40%, and the ring names it: `test_wolf_rayet.py`.
"""

import collections
import math
import random
import types

import pytest

from planetgen import tuning
from planetgen.galaxy import object_uid, profile, spherical
from planetgen.generation import nebula_contents, run_sector
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.supernova_remnant import SupernovaRemnant
from planetgen.galaxy.sector import SpaceSector
from planetgen.util import draw
from tests.audit import census as audit
from tests.audit import representation as rep
from tests.test_sector_gen import _make_cheap_system
from tests.test_wolf_rayet import _wolf_rayet_star

DISK_REGIONS = ("thin-disk-arm", "bulge", "satellite")
"""The regions of the audit that make planetary nebulae in a boosted fill."""


# --- A planetary nebula always has its white dwarf -----------------------------------------------------------------------


@pytest.fixture(scope="module")
def planetary_nebulae():
    """Every planetary nebula the boosted fills of three regions made, with its sector."""
    patch = pytest.MonkeyPatch()
    try:
        rep.boosted_rates(patch)
        found = []
        for case in rep.region_cases(DISK_REGIONS, sectors=60):
            for index in range(case.sectors):
                sector = audit.fill_region_sector(case, audit.seed_for("pn-" + case.name, index))
                found.extend((sector, entry) for entry in sector.phenomena
                             if getattr(entry.phenomenon, "nebula_type", None) == "planetary")
        yield found
    finally:
        patch.undo()


def test_every_planetary_nebula_sits_on_a_hot_white_dwarf_system_of_its_own(planetary_nebulae):
    assert len(planetary_nebulae) >= 10
    for sector, entry in planetary_nebulae:
        hosts = [e for e in sector.entries if e.position == entry.position]
        assert len(hosts) == 1, entry.phenomenon.name
        primary = hosts[0].star_system.stars[0]
        assert run_sector._spectral_code(primary) in tuning.PLANETARY_NEBULA_CENTRAL_STAR_TYPES


def test_a_planetary_nebula_with_a_companion_still_has_the_white_dwarf_as_its_primary(planetary_nebulae):
    paired = 0
    for sector, entry in planetary_nebulae:
        system = next(e for e in sector.entries if e.position == entry.position).star_system
        paired += len(system.stars) > 1
        assert run_sector._spectral_code(system.stars[0]) in tuning.PLANETARY_NEBULA_CENTRAL_STAR_TYPES
    assert paired < len(planetary_nebulae)       # the companion is a minority (a fifth)


# --- A core-collapse remnant has a compact core at the stated share ---------------------------------------------------------


def _remnants(n, type_ia_share, seed):
    draw.set_run_seed(seed)
    config = SystemConfig()
    config.TYPE_IA_SHARE = type_ia_share
    return [SupernovaRemnant(config) for _ in range(n)]


def test_a_core_collapse_remnant_holds_a_compact_core_in_the_stated_share_and_a_type_ia_one_holds_none():
    found = _remnants(3000, 0.0, 184)
    assert all(r.progenitor_type == "core-collapse" for r in found)
    cores = sum(r.compact_remnant is not None for r in found)
    wanted = tuning.SUPERNOVA_CORE_COLLAPSE_REMNANT_VISIBLE_CHANCE * len(found)
    assert cores == pytest.approx(wanted, abs=4.0 * math.sqrt(wanted * (1.0 - tuning.SUPERNOVA_CORE_COLLAPSE_REMNANT_VISIBLE_CHANCE)))
    holes = sum(type(r.compact_remnant).__name__ == "BlackHole" for r in found if r.compact_remnant is not None)
    wanted_holes = tuning.SUPERNOVA_CORE_COLLAPSE_BLACK_HOLE_CHANCE * cores
    assert holes == pytest.approx(wanted_holes, abs=4.0 * math.sqrt(wanted_holes))
    assert all((r.compact_offset_ly is None) == (r.compact_remnant is None) for r in found)
    assert all(r.compact_remnant is None for r in _remnants(500, 1.0, 185))


# --- Every hosted nebula holds its host star; a Wolf-Rayet ring holds a W star ------------------------------------------------


HOT_TYPES = ("O3V", "O5V", "O9V", "B0V", "B2V", "B5V", "B9V", "A0V", "A5V")


def test_every_hosted_nebula_is_on_a_hot_star_of_its_own_and_names_it(monkeypatch):
    monkeypatch.setattr(run_sector, "_sample_poisson_count", lambda mean: 0)
    monkeypatch.setattr(tuning, "NEBULA_HOST_RULES", tuple((*rule[:4], 1.0) for rule in tuning.NEBULA_HOST_RULES))
    sector = SpaceSector("Hosted Nebulae Audit", edge_ly=200.0)
    for star_type in HOT_TYPES * 2 + ("G2V", "K5V", "M3V"):
        system, config = _make_cheap_system(star_type)
        config.BINARY_SYSTEM = False
        sector.add_system(system, system_config=config)
    for entry in sector.entries:
        entry.star_system.stars[0].age = 0.003     # 3 Myr, still in its cloud
    made = run_sector.add_star_hosted_nebulae(sector, types.SimpleNamespace())
    assert len(made) == len(HOT_TYPES) * 2          # every hot star hosts one, no cool star does
    positions = collections.Counter(entry.position for entry in made)
    assert max(positions.values()) == 1
    for entry in made:
        host = next(e for e in sector.entries if e.position == entry.position)
        code = run_sector._spectral_code(host.star_system.stars[0])
        assert code in HOT_TYPES
        assert entry.phenomenon.hosts.split(" x")[0] == code or code in entry.phenomenon.hosts
        assert entry.phenomenon.age_myr == pytest.approx(3.0)


def test_every_wolf_rayet_ring_holds_the_w_star_that_lit_it():
    rings = 0
    for seed in range(200):
        star = _wolf_rayet_star(seed)
        nebula = nebula_contents.hosted_nebula_of(star)
        if nebula is None:
            continue
        rings += 1
        assert nebula.nebula_class == "X"
        assert star.type.startswith("W")
        assert nebula.hosts == f"{star.type.split(' ', 1)[0]} x1"
    assert rings > 40


# --- Tiling and the object ID round trip at the widest shell of the reach ----------------------------------------------------


REACH_SHELL = 65535
EDGE_PC = tuning.DEFAULT_SECTOR_EDGE_PC


def _cells_of_the_widest_shell():
    lowest, highest = spherical.band_range(REACH_SHELL)
    middle = (lowest + highest) // 2
    bands = sorted({lowest, lowest + 1, -1000, -1, 0, 1, 1000, middle, highest - 1, highest})
    cells = []
    for band in bands:
        count = spherical.wedge_count(REACH_SHELL, band)
        for wedge in sorted({w for w in (0, 1, count // 3, count // 2, count - 1) if w < count}):
            cells.append((REACH_SHELL, band, wedge))
    return cells


def test_the_reach_is_nameable_and_the_next_shells_past_it_are_refused_in_time():
    spherical.check_reach(REACH_SHELL)
    spherical.check_reach(spherical.MAX_ADDRESSABLE_SHELL)
    with pytest.raises(ValueError):
        spherical.check_reach(spherical.MAX_ADDRESSABLE_SHELL + 1)


@pytest.mark.parametrize("address", _cells_of_the_widest_shell())
def test_a_cell_of_the_widest_shell_holds_its_own_centre_and_names_itself_by_its_designation(address):
    centre = spherical.sector_position_pc(*address, EDGE_PC)
    assert spherical.sector_address_at(centre, EDGE_PC) == address
    assert spherical.sector_cell_test(address, EDGE_PC)(centre)
    assert spherical.parse_sector_designation(spherical.provisional_sector_designation(*address)) == address
    for serial_kind in object_uid.SERIAL_KINDS:
        value = object_uid.pack(*address, serial_kind, 12345, body=7)
        assert object_uid.unpack(value) == object_uid.ObjectId(*address, serial_kind, 12345, 7)
        assert object_uid.from_bytes(object_uid.to_bytes(value)) == value
        assert object_uid.parse_id(object_uid.format_id(value)) == value
        assert object_uid.sector_bytes(*address) == object_uid.to_bytes(value)[:7]


@pytest.mark.parametrize("address", _cells_of_the_widest_shell())
def test_the_neighbours_of_a_cell_of_the_widest_shell_are_symmetric_and_share_no_point(address):
    centre = spherical.sector_position_pc(*address, EDGE_PC)
    neighbours = spherical.neighbor_addresses(*address)
    assert address not in neighbours and len(set(neighbours)) == len(neighbours) >= 4
    for other in neighbours:
        assert address in spherical.neighbor_addresses(*other), (address, other)
        assert not spherical.sector_cell_test(other, EDGE_PC)(centre)


def test_random_points_of_the_widest_shell_each_land_in_exactly_one_cell():
    rng = random.Random(184)
    inner, outer = spherical.shell_bounds_pc(REACH_SHELL, EDGE_PC)
    landed = collections.Counter()
    for _ in range(300):
        radius = rng.uniform(inner, outer)
        cosine = rng.uniform(-1.0, 1.0)
        sine = math.sqrt(1.0 - cosine * cosine)
        angle = rng.uniform(0.0, 2.0 * math.pi)
        point = (radius * sine * math.cos(angle), radius * sine * math.sin(angle), radius * cosine)
        address = spherical.sector_address_at(point, EDGE_PC)
        landed[address[0]] += 1
        assert spherical.sector_cell_test(address, EDGE_PC)(point)
        for other in spherical.neighbor_addresses(*address):
            assert not spherical.sector_cell_test(other, EDGE_PC)(point)
        assert spherical.parse_sector_designation(spherical.provisional_sector_designation(*address)) == address
    assert landed[REACH_SHELL] >= 290          # a point on a shell edge may round into its neighbour, a few at most


# --- A profile and an empty sector are deterministic -------------------------------------------------------------------------


@pytest.mark.parametrize("name", list(rep.REGION_POINTS))
def test_a_profile_is_the_same_every_time_it_is_asked_for_in_every_region(name):
    point = rep.REGION_POINTS[name][0]
    first = profile.profile_at(point, rep.WITH_COMPONENTS, rep.EDGE, rep.E1)
    again = profile.profile_at(point, rep.WITH_COMPONENTS, rep.EDGE, rep.E1)
    assert first == again and first.mix == again.mix and first.kind_factors == again.kind_factors
    assert profile.profile_at((point[0] + 1.0e-6, point[1], point[2]), rep.WITH_COMPONENTS, rep.EDGE, rep.E1).mix.keys() \
        == first.mix.keys()


@pytest.mark.parametrize("name", list(rep.REGION_POINTS))
def test_an_empty_sector_is_the_same_on_every_fill_in_every_region(name):
    case = audit.RegionCase(name, rep.REGION_POINTS[name][0], rep.WITH_COMPONENTS, rep.EDGE, rep.E1, sectors=1,
                            num_systems=0)
    first = audit.fill_region_sector(case, "empty-sector")
    again = audit.fill_region_sector(case, "empty-sector")
    assert first.stage_counts == again.stage_counts == collections.Counter()
    assert len(first.entries) == len(again.entries) == 0
    assert len(first.phenomena) == len(again.phenomena) == 0
    assert first.name == again.name
