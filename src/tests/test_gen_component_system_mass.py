# tests/test_gen_component_system_mass.py -- GEN.598: a cluster's or satellite's systems are its mass over the mean mass of
# its own systems, worked out from its age with the same mass function, remnants and multiplicity the galaxy's are (GEN.545).

import math
import statistics

import pytest

from planetgen.galaxy import components
from planetgen.generation import multistar, system_mass
from planetgen.physics import remnant_mass, remnant_progenitors, stellar_evolution
from planetgen.util import draw


def _component(age_gy, kind="cluster", mass=2.0e5):
    return components.Component(
        kind=kind, name="Test", source="generated", center_pc=(0.0, 0.0, 0.0), velocity_kms=(0.0, 0.0, 0.0), epoch_unix=0.0,
        mass_solar=mass, half_mass_radius_pc=5.0, profile="plummer", age_gy=age_gy, metallicity_fe_h=-0.5, seed=1,
        tidal_radius_pc=50.0)


@pytest.mark.parametrize("age_gy", [0.03, 0.3, 3.0, 12.0])
def test_a_components_mean_is_a_seeded_draw_of_systems_of_its_age(age_gy):
    rng = draw.Stream(f"gen598-systems-{age_gy}")
    found = []
    for _ in range(4000):
        mass = stellar_evolution.sample_imf_mass_sol(rng=rng)
        if stellar_evolution.evolve_star(mass, age_gy, rng) is None:
            kind = "black-hole" if remnant_progenitors.fails(mass) else "neutron-star"
            found.append(remnant_mass.draw_mass_sol(kind, rng=rng))
            continue
        found.append(multistar.draw_member_tree(multistar.Placed("star", mass), age_gy * 1.0e9, rng=rng).system_mass_solar)
    error = statistics.stdev(found) / math.sqrt(len(found))
    assert abs(statistics.fmean(found) - _component(age_gy).mean_system_mass_solar) < 4.0 * error + 0.02 * statistics.fmean(found)


def test_the_systems_are_the_mass_over_the_mean_and_a_younger_component_holds_fewer():
    old, young = _component(12.0), _component(0.03)
    assert old.mean_system_mass_solar == pytest.approx(system_mass.mean_system_mass_at_age_solar(12.0))
    assert old.expected_systems == pytest.approx(old.mass_solar / old.mean_system_mass_solar)
    assert young.mean_system_mass_solar > old.mean_system_mass_solar
    assert young.expected_systems < old.expected_systems
