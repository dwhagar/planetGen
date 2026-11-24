# tests/test_pair_orbit.py

"""
GEN.640 slice b: `generation.pair_orbit.PairOrbit`, the pair of a system read from its orbit node (it replaced the `WideBinaryPair` class, which kept a copy).
The numbers are the node's: the axis, the eccentricity, the period, the angles; the position is the Kepler ellipse of `system_tree.node_relative_au` and
the two stars' places about the center of mass follow from the mass fraction.
"""

import dataclasses
import math

import pytest

from planetgen.generation import pair_orbit
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import constants
from planetgen.util import draw

AU_KM = constants.AU_TO_KM


def _wide(seed):
    config = SystemConfig()
    config.BINARY_SYSTEM, config.WIDE_BINARY, config.PLANETS = True, True, False
    with draw.bound(seed):
        return StarSystem(config)


@pytest.mark.parametrize("seed", range(6))
def test_a_wide_pair_is_its_node(seed):
    system = _wide(seed)
    pair = system.pair_orbit
    (node,) = system.orbit_nodes
    assert pair.separation_au * AU_KM == pytest.approx(node.separation_km)
    assert pair.eccentricity == node.eccentricity
    assert pair.period_years == node.period_years and pair.phase_deg == node.mean_anomaly_deg
    assert pair.periapsis_au == pytest.approx(pair.separation_au * (1.0 - pair.eccentricity))
    assert pair.apoapsis_au == pytest.approx(pair.separation_au * (1.0 + pair.eccentricity))
    assert pair.secondary_mass_fraction == pytest.approx(
        system.secondary_star.mass / (system.primary_star.mass + system.secondary_star.mass))


@pytest.mark.parametrize("seed", range(6))
def test_the_two_stars_sit_either_side_of_the_center_of_mass(seed):
    pair = _wide(seed).pair_orbit
    for axis in "xyz":
        mutual = getattr(pair, f"position_{axis}_au")
        primary, secondary = getattr(pair, f"primary_position_{axis}_au"), getattr(pair, f"secondary_position_{axis}_au")
        assert secondary - primary == pytest.approx(mutual)
        assert primary == pytest.approx(-pair.secondary_mass_fraction * mutual)
    distance = math.sqrt(sum(getattr(pair, f"position_{axis}_au") ** 2 for axis in "xyz"))
    assert pair.periapsis_au * (1 - 1e-9) <= distance <= pair.apoapsis_au * (1 + 1e-9)


def test_a_pair_the_other_way_round_turns_the_orbit_half_way():
    node = _wide(1).orbit_nodes[0]
    straight = pair_orbit.pair_orbit(node, 1.0, 0.5)
    turned = pair_orbit.pair_orbit(node, 1.0, 0.5, turned=True)
    assert turned.separation_au == straight.separation_au and turned.period_years == straight.period_years
    assert turned.position_x_au == pytest.approx(-straight.position_x_au) and turned.position_y_au == pytest.approx(-straight.position_y_au)
    assert turned.position_z_au == pytest.approx(-straight.position_z_au)


def test_the_fraction_is_the_lighter_stars_whichever_way_the_masses_are_given():
    node = _wide(2).orbit_nodes[0]
    assert pair_orbit.pair_orbit(node, 1.0, 0.25).secondary_mass_fraction == pytest.approx(0.2)
    assert pair_orbit.pair_orbit(node, 0.25, 1.0).secondary_mass_fraction == pytest.approx(0.2)


def test_a_single_star_has_no_pair_and_the_pair_follows_its_node():
    config = SystemConfig()
    config.BINARY_SYSTEM, config.PLANETS = False, False
    assert StarSystem(config).pair_orbit is None
    system = _wide(3)
    node = system.orbit_nodes[0]
    system.orbit_nodes = (dataclasses.replace(node, separation_km=300.0 * AU_KM, eccentricity=0.5),)
    assert system.pair_orbit.separation_au == pytest.approx(300.0) and system.pair_orbit.periapsis_au == pytest.approx(150.0)
