# tests/test_gen367_cloud_gas.py

"""
GEN.367: a Magellanic Cloud holds the molecular clouds its gas would, not the clouds its whole young stellar ball would.

The nebula field draws clouds from a tracer that follows a satellite's young stars (GEN.528). Those stars fill a ball out to six
half-mass radii, so the Large Cloud came out with about 2.5 times the clouds of its gas. `tuning.SATELLITE_NEBULA_GAS_SCALE` scales a
named satellite's tracer, and these tests work the target again from the sources [R3102, R5123, R5124] and hold the scale to it.
"""

import itertools
import math

import pytest

from planetgen import tuning
from planetgen.galaxy import nebula_field

from tests.test_gen_magellanic_classes import FORNAX, LMC, SHAPE, SMC, SUN, WITH_SATELLITES

GRID_PC = 400.0
"""float: The step of the grid the clouds are summed on; 300 and 400 pc agree to within 1 percent, 300 and 200 likewise."""


def _class_m_kept_per_drawn():
    dark = {letter: tuning.NEBULA_CLASSES[letter]["frequency"] for letter in nebula_field.FIELD_CLASSES}
    return dark["M"] / sum(dark.values()) * tuning.NEBULA_FIELD_KEEP["M"]


def _expected_class_m_clouds(cloud, shape):
    """The class M clouds the field keeps in a Cloud, on average: its drawn mean over the tidal ball, times the share kept."""
    density = tuning.PHENOMENON_DENSITY_PC3["molecular-cloud"] * tuning.PHENOMENON_RATE_SCALE.get("molecular-cloud", 1.0)
    reach = cloud.half_mass_radius_pc * tuning.SATELLITE_TIDAL_OVER_HALF_MASS
    n = int(reach // GRID_PC)
    total = 0.0
    for i, j, k in itertools.product(range(-n, n + 1), repeat=3):
        point = (cloud.center_pc[0] + i * GRID_PC, cloud.center_pc[1] + j * GRID_PC, cloud.center_pc[2] + k * GRID_PC)
        if math.dist(point, cloud.center_pc) <= reach:
            total += density * GRID_PC ** 3 * nebula_field.gas_factor(point, shape)
    return total * _class_m_kept_per_drawn()


def _target(name):
    return (tuning.MILKY_WAY_CATALOGUED_CLOUDS * tuning.SATELLITE_MOLECULAR_GAS_MASS_SOLAR[name]
            / tuning.MILKY_WAY_MOLECULAR_GAS_MASS_SOLAR)


def test_the_sourced_numbers_are_the_papers():
    assert tuning.MILKY_WAY_CATALOGUED_CLOUDS == 8107 and tuning.MILKY_WAY_MOLECULAR_GAS_MASS_SOLAR == 1.2e9
    assert tuning.SATELLITE_MOLECULAR_GAS_MASS_SOLAR == {"Large Magellanic Cloud": 5.0e7, "Small Magellanic Cloud": 3.2e7}
    assert _target("Large Magellanic Cloud") == pytest.approx(337.8, abs=0.1)
    assert _target("Small Magellanic Cloud") == pytest.approx(216.2, abs=0.1)


def test_the_large_clouds_class_m_clouds_are_its_share_of_the_milky_ways_by_molecular_mass():
    found = _expected_class_m_clouds(LMC, WITH_SATELLITES)
    assert found == pytest.approx(_target(LMC.name), rel=0.15)
    assert 0.8 < found / 272 < 1.6      # the 272 clouds the CO survey catalogued above its limit [R5123]


def test_the_small_clouds_class_m_clouds_are_its_share_of_the_milky_ways_by_molecular_mass():
    assert _expected_class_m_clouds(SMC, WITH_SATELLITES) == pytest.approx(_target(SMC.name), rel=0.15)


def test_a_satellite_that_is_not_named_keeps_the_young_star_tracer(monkeypatch):
    """The Large Cloud, renamed, has no scale: its gas tracer is its young stars' density as before GEN.367."""
    point = tuple(LMC.center_pc)
    components = WITH_SATELLITES.components
    scaled = components.gas_density(point)
    young = components.densities(point)["satellite_young"]
    assert scaled == pytest.approx(young * tuning.SATELLITE_NEBULA_GAS_SCALE[LMC.name], rel=1e-9)
    monkeypatch.setattr(tuning, "SATELLITE_NEBULA_GAS_SCALE", {})
    assert components.gas_density(point) == pytest.approx(young, rel=1e-12)


def test_away_from_the_clouds_and_for_a_dwarf_without_young_stars_the_gas_is_what_it_was():
    components = WITH_SATELLITES.components
    assert components.gas_density(SUN) == components.densities(SUN)["satellite_young"] == 0.0
    assert components.gas_density(tuple(FORNAX.center_pc)) == 0.0
    assert nebula_field.gas_factor(SUN, WITH_SATELLITES) == nebula_field.gas_factor(SUN, SHAPE)


def test_the_large_clouds_scale_is_the_target_over_the_unscaled_count_to_the_power_one_over_the_exponent(monkeypatch):
    """Unscaled the field keeps about 847 class M clouds in the Large Cloud; the count goes as the tracer to the 1.4 there."""
    monkeypatch.setattr(tuning, "SATELLITE_NEBULA_GAS_SCALE", {})
    unscaled = _expected_class_m_clouds(LMC, WITH_SATELLITES)
    assert unscaled == pytest.approx(847.0, rel=0.03)
    solved = (_target(LMC.name) / unscaled) ** (1.0 / tuning.GMC_GAS_DENSITY_EXPONENT)
    assert solved == pytest.approx(0.52, abs=0.02)
