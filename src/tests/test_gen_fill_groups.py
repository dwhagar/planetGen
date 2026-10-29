# tests/test_gen_fill_groups.py

"""
GEN.328: young groups at fill. A massive member of an emission region (D or E) or a field cloud rolls the `young-group`
share on a stream of its own seed (`nebula_contents.member_group`); a group that has not decayed is built as its member and
the mates (`run_sector.add_group_stars`), a group older than its decay age as the hierarchy it left, and a group across a
sector boundary is made once, each star by the sector that holds it.
"""

import collections
import types

from planetgen import tuning
from planetgen.generation import nebula_contents, run_sector
from planetgen.generation.phenomena.system_types import SYSTEM_TYPE_RATES
from planetgen.physics import constants as physical_constants
from tests.test_sector_gen import _real_sector_args

SEED = bytes(range(16))
LIVE_AGE_YR = 5.0e3      # under every decay age (1e4 to 4e4 years)
DECAYED_AGE_YR = 2.0e5   # over every decay age, under the 0.3 Myr limit


def _member(number, mass=25.0, age_yr=LIVE_AGE_YR, position=(8000.0, 0.0, 0.0)):
    return nebula_contents.CloudMember(position, mass, age_yr, f"gen328:{number}")


def _groups(count, **kwargs):
    found = []
    for number in range(count):
        member = _member(number, **kwargs)
        tree = nebula_contents.member_group(member)
        if tree is not None:
            found.append((member, tree))
    return found


def test_a_member_draws_the_same_group_wherever_it_is_asked():
    for member, tree in _groups(400):
        again = nebula_contents.member_group(member)
        assert again == tree


def test_a_massive_young_member_is_in_a_group_at_the_default_share():
    n = 3000
    share = len(_groups(n)) / n
    default = SYSTEM_TYPE_RATES["young-group"].default
    assert 0.6 * default < share < 1.5 * default


def test_only_a_massive_member_younger_than_the_limit_is_in_a_group():
    assert _groups(600, mass=1.0) == []                    # a protostar of a core
    assert _groups(600, mass=tuning.PRE_MAIN_SEQUENCE_MAX_MASS_SOL - 0.5) == []
    assert _groups(600, age_yr=tuning.YOUNG_GROUP_AGE_MAX_YR * 1.01) == []   # an emission region's usual age


def test_a_group_that_has_not_decayed_holds_the_member_and_its_mates_younger_than_the_decay_age():
    found = _groups(400)
    assert found
    for member, tree in found:
        assert nebula_contents.is_live_group(tree)
        assert 3 <= tree.member_count <= 5
        assert member.mass_sol == tree.members[0].mass_solar
        for mate in tree.members:
            assert mate.offset_au is not None
            assert member.age_yr < mate.decay_age_yr


def test_a_group_past_its_decay_age_is_the_hierarchy_it_left():
    found = _groups(600, age_yr=DECAYED_AGE_YR)
    assert found
    for member, tree in found:
        assert not nebula_contents.is_live_group(tree)
        assert tree.member_count >= 2 and tree.members[0].mass_solar == member.mass_sol
        assert all(mate.offset_au is None for mate in tree.members)


def test_an_emission_regions_reach_holds_the_widest_groups_mates():
    widest_pc = 2.0 * tuning.YOUNG_GROUP_RADIUS_AU_RANGE[1] / physical_constants.AU_PER_PARSEC
    assert nebula_contents.group_reach_pc() >= widest_pc
    assert nebula_contents.emission_reach_pc() > nebula_contents.group_reach_pc()


class _Sector:
    """A stand-in sector that holds one side of a plane and records the systems it is given."""

    def __init__(self, side):
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


def _args():
    return _real_sector_args(1, "t", planets=False)


def test_a_group_across_a_sector_boundary_is_made_once_each_star_by_the_sector_that_holds_it():
    centre_pc = (8000.0, 0.0, 0.0)
    args = _args()
    built = []
    expected = 0
    for member, tree in _groups(400)[:4]:
        expected += sum(1 for mate in tree.members if mate.kind in ("main-sequence", "pre-main-sequence", "protostar",
                                                                   "giant", "white-dwarf"))
        for side in (True, False):
            sector = _Sector(side)
            run_sector._add_members(sector, args, [member], centre_pc, 26000.0, "emission-member")
            built += sector.systems
    assert built and len(built) == expected          # none twice, none missing
    by_group = collections.Counter(entry.group_seed for entry in built)
    assert all(3 <= count <= 5 for count in by_group.values())
    assert any(entry.position[0] < 0.0 for entry in built) and any(entry.position[0] >= 0.0 for entry in built)
    for entry in built:
        assert entry.group_decay_age_yr > LIVE_AGE_YR
        assert entry.group_offset_au is not None


def test_a_member_in_a_decayed_group_is_built_from_the_hierarchy_and_marks_no_group():
    centre_pc = (8000.0, 0.0, 0.0)
    args = _args()
    for member, tree in _groups(600, age_yr=DECAYED_AGE_YR)[:3]:
        sector = _Sector(True)
        run_sector._add_members(sector, args, [member], centre_pc, 26000.0, "emission-member")
        assert len(sector.systems) == 1
        entry = sector.systems[0]
        assert not hasattr(entry, "group_seed")
        assert len(entry.star_system.stars) >= 2


def test_a_member_in_no_group_is_built_as_before():
    centre_pc = (8000.0, 0.0, 0.0)
    member = _member(0, mass=1.0)
    sector = _Sector(True)
    run_sector._add_members(sector, _args(), [member], centre_pc, 26000.0, "protostar")
    assert len(sector.systems) == 1 and not hasattr(sector.systems[0], "group_seed")
