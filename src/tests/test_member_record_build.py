# tests/test_member_record_build.py

"""
GEN.1009: the tree and the build read one member record. The radius the tree gave a member is the radius of the star the build
makes for it, for every member of a drawn tree (the guard that replaces `stored_trees.MASS_TOLERANCE`), and a system built from
a drawn tree keeps all its buildable members.
"""

import pytest

from planetgen.generation import member_record, multistar, stored_trees
from planetgen.generation.config import SystemConfig
from planetgen.generation.star import Star
from planetgen.generation.system import StarSystem
from planetgen.physics import constants
from planetgen.util import draw

GY = 1.0e9


def _radius_au(star):
    return star.radius / constants.AU_TO_KM


@pytest.mark.parametrize("count", [2, 3, 4])
def test_the_built_star_of_every_member_has_the_radius_the_tree_drew(count):
    checked = 0
    for seed in range(60):
        with draw.bound(seed):
            primary = Star(SystemConfig(), name="Placed", initial_mass_sol=1.0 + 0.02 * seed, age_gy=2.0 + 0.1 * seed,
                           metallicity_fe_h=0.1)
        placed = multistar.Placed("star", primary.mass / constants.SOLAR_MASS_TO_KG, primary.initial_mass_sol,
                                  fe_h=primary.metallicity_fe_h, record=member_record.of_star(primary))
        plan = stored_trees.drawn(
            multistar.draw_member_tree(placed, primary.age * GY, draw.Stream(seed), n_stars=count, pair_placed=True))
        for member in plan.lit_members[1:]:
            record = member_record.make(member.seed, member.initial_mass_solar, member.age_yr or primary.age * GY,
                                        primary.metallicity_fe_h,
                                        white_dwarf_mass_sol=stored_trees.helium_core_mass_sol(member))
            built = Star.from_record(record, SystemConfig(), name="Member")
            assert stored_trees.matches(built, member)
            assert _radius_au(built) == pytest.approx(member.radius_au, rel=1.0e-9)
            checked += 1
    assert checked > 20


def test_a_forced_triple_keeps_all_three_stars():
    built = 0
    for seed in range(40):
        config = SystemConfig()
        config.STARS = 3
        with draw.bound(seed):
            system = StarSystem(config)
        built += len(system.stars) == 3
    assert built >= 30


def test_the_primarys_record_is_the_primarys_own_radius():
    with draw.bound(3):
        primary = Star(SystemConfig(), name="Placed", initial_mass_sol=1.2, age_gy=6.0, metallicity_fe_h=0.0)
    record = member_record.of_star(primary)
    assert record.radius_au == pytest.approx(_radius_au(primary), rel=1.0e-9)
    assert record.radius_max_au >= record.radius_au
    assert record.mass_max_solar >= record.mass_solar
