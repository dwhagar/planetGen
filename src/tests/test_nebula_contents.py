# tests/test_nebula_contents.py

"""
GEN.320: the stars a nebula needs, made with it (`generation/nebula_contents.py`): the hosted roll from the star's own
seed, the planetary nebula's close companion, and the protostars of Bok globules and star-forming cores, built by the
sector that holds each one's point.
"""

import collections
import types


from planetgen import tuning
from planetgen.generation import nebula_contents, run_sector
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.nebula import Nebula
from planetgen.generation.system import StarSystem
from planetgen.physics.units import ly_to_pc
from planetgen.util import draw
from tests.test_sector_gen import _real_sector_args

SEED = bytes(range(16))


def _star(number, type_="O5V"):
    return types.SimpleNamespace(type=f"{type_} Blue", age=0.0005, mass=2.0e31 * (1.0 + number * 1e-6),
                                 temperature=40000.0)


def test_a_star_rolls_the_same_nebula_wherever_it_is_made():
    """The roll is on a stream of the star's own, so a pre-placed star and the same star in a fill agree."""
    star = _star(1)
    first = nebula_contents.hosted_nebula_of(star, 1.0)
    for _ in range(3):
        again = nebula_contents.hosted_nebula_of(star, 1.0)
        assert (first is None) == (again is None)
        if first is not None:
            assert (first.nebula_class, first.radius_ly, first.hosts, first.density_cm3) == (
                again.nebula_class, again.radius_ly, again.hosts, again.density_cm3)
    draw.set_run_seed(1)
    after_other_draws = nebula_contents.hosted_nebula_of(star, 1.0)
    assert (first is None) == (after_other_draws is None)


def test_a_whole_kelvin_temperature_rolls_the_same_as_its_float():
    """A star made in a fill holds 36900 and its stored row gives back 36900.0: the roll must not tell them apart."""
    whole, stored = _star(5), _star(5)
    whole.temperature, stored.temperature = 36900, 36900.0
    assert nebula_contents.hosted_seed(whole) == nebula_contents.hosted_seed(stored)
    first, second = nebula_contents.hosted_nebula_of(whole), nebula_contents.hosted_nebula_of(stored)
    assert (first is None) == (second is None)


def test_the_hosted_chance_is_the_rules_and_the_profile_scales_it():
    stars = [_star(n) for n in range(1500)]
    hits = sum(nebula_contents.hosted_nebula_of(star, 1.0) is not None for star in stars) / len(stars)
    assert abs(hits - 0.25) < 0.04
    half = sum(nebula_contents.hosted_nebula_of(star, 0.4) is not None for star in stars) / len(stars)
    assert abs(half - 0.10) < 0.035
    assert all(nebula_contents.hosted_nebula_of(star, 0.0) is None for star in stars[:100])


def test_a_star_too_old_or_of_another_type_lights_nothing():
    old = _star(2)
    old.age = 0.5
    assert nebula_contents.hosted_nebula_of(old, 1.0) is None
    assert nebula_contents.hosted_nebula_of(_star(3, "G2V"), 1.0) is None
    assert nebula_contents.hosted_nebula_of(_star(4, "O5III"), 1.0) is None


def test_a_hosted_nebula_holds_its_host_in_the_star_type_it_was_rolled_for():
    found = [n for n in (nebula_contents.hosted_nebula_of(_star(k), 1.0) for k in range(200)) if n is not None]
    assert found and all(n.nebula_class in "CDE" and n.hosts for n in found)


def test_a_planetary_nebulas_central_star_has_a_close_companion_in_a_fifth():
    draw.set_run_seed(11)
    config_hits = 0
    for _ in range(3000):
        config = nebula_contents.planetary_central_config(SystemConfig())
        assert config.WIDE_BINARY is False
        config_hits += config.BINARY_SYSTEM
    assert abs(config_hits / 3000 - tuning.PLANETARY_NEBULA_COMPANION_CHANCE) < 0.03


def test_a_forced_companion_makes_a_binary_around_the_white_dwarf_and_none_makes_a_single():
    for seed in range(8):
        draw.set_run_seed(seed)
        config = SystemConfig()
        config.STAR_TYPE = "O5VII"
        config.BINARY_SYSTEM = True
        config.WIDE_BINARY = False
        config.PLANETS = False
        system = StarSystem(system_config=config)
        assert system.binary_type is not None and system.primary_star.yerkes_class == "VII"
        config.BINARY_SYSTEM = False
        assert StarSystem(system_config=config).binary_type is None


def _cloud(letter):
    draw.set_run_seed(5)
    return Nebula(SystemConfig(), nebula_class=letter)


def test_a_globule_has_a_protostar_in_three_of_ten_and_a_core_one_to_four():
    globule = _cloud("P")
    held = [len(nebula_contents.cloud_members(SEED, globule, (1000.0 + k, 0.0, 0.0))) for k in range(1500)]
    assert set(held) <= {0, 1} and abs(sum(held) / len(held) - 0.30) < 0.04
    core = _cloud("Q")
    counts = collections.Counter(len(nebula_contents.cloud_members(SEED, core, (1000.0 + k, 0.0, 0.0)))
                                 for k in range(1500))
    assert set(counts) == {1, 2, 3, 4} and counts[0] == 0
    assert counts[1] > counts[4]
    assert nebula_contents.cloud_members(SEED, _cloud("M"), (0.0, 0.0, 0.0)) == []


def test_the_members_are_young_light_stars_inside_the_cloud_and_the_same_every_time():
    core = _cloud("Q")
    radius = ly_to_pc(core.radius_ly)
    for k in range(60):
        centre = (500.0 + 3.7 * k, 20.0, -4.0)
        members = nebula_contents.cloud_members(SEED, core, centre)
        again = nebula_contents.cloud_members(SEED, core, centre)
        assert [(m.position_pc, m.mass_sol, m.age_yr, m.seed) for m in members] == [
            (m.position_pc, m.mass_sol, m.age_yr, m.seed) for m in again]
        for m in members:
            assert sum((m.position_pc[i] - centre[i]) ** 2 for i in range(3)) <= radius ** 2 * (1 + 1e-9)
            assert 0.08 <= m.mass_sol < tuning.PRE_MAIN_SEQUENCE_MAX_MASS_SOL
            assert 0.0 <= m.age_yr < tuning.PROTOSTAR_MAX_AGE_YR
    other = nebula_contents.cloud_members(bytes(16), core, (500.0, 20.0, -4.0))
    assert [m.seed for m in other] != [m.seed for m in nebula_contents.cloud_members(SEED, core, (500.0, 20.0, -4.0))]


def _args():
    return _real_sector_args(1, "t", planets=False)


class _Sector:
    """A stand-in sector that holds one side of a plane and records the systems it is given."""

    def __init__(self, clouds, side):
        self.field_nebulae = clouds
        self.side = side
        self.systems = []

    def contains(self, local_ly):
        return (local_ly[0] >= 0.0) == self.side

    def add_preplaced_system(self, system, position, system_config=None):
        entry = types.SimpleNamespace(star_system=system, position=position)
        self.systems.append(entry)
        return entry

    stage_base = None
    stage_counts = collections.Counter()


def test_each_sector_builds_only_the_members_inside_it_and_together_they_make_every_one():
    core = _cloud("Q")
    centre_pc = (8000.0, 0.0, 0.0)
    cloud_centre = (8000.0, 0.0, 0.0)
    cloud_field = (SEED, None, centre_pc, 3.0)
    args = _args()
    expected = 0
    built = []
    for k in range(10):
        cloud = (core, (cloud_centre[0] + 0.01 * k, 0.0, 0.0))
        expected += len(nebula_contents.cloud_members(SEED, core, cloud[1]))
        halves = [_Sector([cloud], side) for side in (True, False)]
        for sector in halves:
            run_sector.add_cloud_protostars(sector, args, cloud_field, 26000.0)
            built += sector.systems
    assert expected >= 10 and len(built) == expected  # none twice, none missing
    for entry in built:
        star = entry.star_system.primary_star
        assert star.yerkes_class == "protostar" and star.disk_kind in ("class-0", "class-i")
        assert [p for p in entry.star_system.planets if p.body_type != "a"] == []


def test_a_cloud_without_protostars_adds_no_system():
    sector = _Sector([(_cloud("M"), (8000.0, 0.0, 0.0))], True)
    assert run_sector.add_cloud_protostars(sector, _args(),
                                           (SEED, None, (8000.0, 0.0, 0.0), 3.0)) == []
