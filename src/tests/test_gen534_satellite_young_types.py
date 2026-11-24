# tests/test_gen534_satellite_young_types.py

"""
GEN.534: a satellite galaxy's young stars reach the system-type tables. The types the age gates name a disk population as host; a satellite's stars of the same
ages (`satellite_young`, `satellite_intermediate`, `satellite`, GEN.526) now name them too, so a system asked inside the Large Magellanic Cloud can be a young-star
type with the weight of the Cloud's own young share, and one inside a dwarf that holds no young stars cannot.
"""

import random

import pytest

from planetgen import tuning
from planetgen.galaxy import component_field, profile
from planetgen.generation import multistar, star_population
from planetgen.generation.phenomena import system_types as types
from tests import test_satellites as satellite_tests
from tests.test_satellites import E1, EDGE, SHAPE

PLACED = satellite_tests._placed()
WITH_SATELLITES = SHAPE._replace(components=component_field.ComponentField(PLACED, EDGE, E1))
LMC, FORNAX = (next(c for c in PLACED if c.name == name) for name in ("Large Magellanic Cloud", "Fornax dwarf"))
YOUNG_ONLY = ("young-group", "multiple-with-disk", "be-star", "wolf-rayet-binary", "high-mass-xray-binary")
"""tuple: Types of the young population alone (the disk's list before GEN.534)."""


def _host_share(key, component):
    """The share of the stars in the cell at `component`'s center whose population may host the type: the profile's own mix, with no weight of the type's."""
    mix = profile.profile_at(component.center_pc, WITH_SATELLITES, EDGE, E1).mix
    return sum(mix[name] for name in types.SYSTEM_TYPE_RATES[key].populations) / sum(mix.values())


def test_every_type_that_names_a_disk_population_names_its_satellite_twin():
    for key, entry in types.SYSTEM_TYPE_RATES.items():
        for disk, twin in types.SATELLITE_OF.items():
            assert (disk in entry.populations) <= (twin in entry.populations), key
        assert set(entry.populations) <= set(types.POPULATIONS), key
    assert types.SYSTEM_TYPE_RATES["algol"].populations == ("young", "intermediate", "satellite_young", "satellite_intermediate")
    assert types.SYSTEM_TYPE_RATES["be-star"].populations == ("young", "satellite_young")


@pytest.mark.parametrize("key", YOUNG_ONLY)
def test_a_young_star_type_is_hosted_in_the_large_cloud_by_exactly_its_young_share_and_in_a_dwarf_with_no_young_stars_not_at_all(key):
    found = component_field.population_shares(LMC)
    assert found["satellite_young"] > 0.005
    assert _host_share(key, LMC) == pytest.approx(found["satellite_young"], rel=1e-3)
    assert component_field.population_shares(FORNAX)["satellite_young"] == 0.0
    assert _host_share(key, FORNAX) == 0.0


def test_a_type_of_any_population_is_hosted_everywhere():
    assert _host_share("close-binary", LMC) == pytest.approx(1.0) and _host_share("close-binary", FORNAX) == pytest.approx(1.0)


def test_the_stars_a_cloud_draws_include_young_ones_that_the_young_star_types_gate_on_and_a_dwarfs_do_not():
    low, high = tuning.STELLAR_POPULATION_AGE_RANGES_GY["satellite_young"]
    assert high <= tuning.YOUNG_POPULATION_MAX_AGE_GY, "a satellite's young star is young to the gate (`multistar._young`)"
    for component, wanted in ((LMC, True), (FORNAX, False)):
        mix = profile.profile_at(component.center_pc, WITH_SATELLITES, EDGE, E1).mix
        rng = random.Random(7)
        picks = [star_population.pick_population(mix, rng) for _ in range(6000)]
        assert ("satellite_young" in picks) is wanted
    be, _stripped = multistar.draw_be_star(4.0, high * 1.0e9 / 2.0, 0.0, 1.0)
    assert be, "a B star of a satellite's young window can be a Be star"
