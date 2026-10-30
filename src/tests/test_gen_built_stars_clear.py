# tests/test_gen_built_stars_clear.py

"""
A system of three or more stars never holds one star inside another. A drawn tree is built again from each member's birth
mass and the system's age, so a member the tree had evolved (a white dwarf of a post-common-envelope pair) can come out
a giant wider than its orbit; `StarSystem` draws the tree again when the stars as built do not clear its orbits
(`multistar.built_tree_clears`).
"""

import pytest

from planetgen.generation import multistar
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import constants
from planetgen.physics.system_tree import OrbitNode
from planetgen.util import draw

AU_KM = constants.AU_TO_KM


def _node(node_no, parent, axis_au, e=0.0):
    return OrbitNode(node_no=node_no, parent_node_no=parent, slot=0 if parent is not None else 0, mass_solar=1.0,
                     separation_km=axis_au * AU_KM, eccentricity=e, inclination_deg=0.0, ascending_node_deg=0.0,
                     arg_periapsis_deg=0.0, mean_anomaly_deg=0.0, epoch_unix=0.0, period_years=1.0, kind="close",
                     interaction="detached", flags=0)


def test_a_pair_with_a_giant_wider_than_its_orbit_does_not_clear():
    nodes = (_node(0, None, 6.8), _node(1, 0, 0.0098), _node(2, 0, 0.0088))
    giant, dwarf = (2, 0.0335, 1.23), (2, 0.00226, 0.51)
    assert not multistar.built_tree_clears(nodes, [(1, 0.0029, 0.64), (1, 0.00006, 0.63), dwarf, giant])


def test_the_same_pair_with_the_white_dwarf_the_tree_drew_clears():
    nodes = (_node(0, None, 6.8), _node(1, 0, 0.0098), _node(2, 0, 0.0088))
    assert multistar.built_tree_clears(nodes, [(1, 0.0029, 0.64), (1, 0.00006, 0.63), (2, 0.00226, 0.51),
                                               (2, 0.00006, 0.60)])


def test_a_third_star_whose_periapsis_is_inside_the_inner_pair_does_not_clear():
    nodes = (_node(0, None, 0.4, e=0.5), _node(1, 0, 0.3))
    # the inner pair reaches 0.3 * 1.0 = 0.3 AU out (e = 0); the outer star comes in to 0.2 AU
    assert not multistar.built_tree_clears(nodes, [(0, 0.005, 0.5), (1, 0.005, 0.5), (1, 0.005, 0.5)])


@pytest.mark.parametrize("seed", [241, 297])
def test_a_forced_triple_the_first_draw_built_a_star_inside_another_clears_now(seed):
    config = SystemConfig()
    config.PLANETS = False
    config.STARS = 3
    with draw.bound(seed):
        system = StarSystem(system_config=config)
    assert len(system.stars) == 3
    assert system._built_stars_clear()
