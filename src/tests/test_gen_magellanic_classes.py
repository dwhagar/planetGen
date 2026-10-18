# tests/test_gen_magellanic_classes.py

"""
GEN.528: a cell inside a Magellanic Cloud can show every class its profile allows, and the planetary nebulae of the
scatters follow the old stars' share as a sector's do. A satellite's young stars (GEN.526) count in the gas tracer, so
its cells hold molecular clouds, nebulae and supernova remnants; its component pass places supernova remnants where the
gas is and magnetars in proportion to its young stars; and the scatters weigh a planetary nebula by the share of the
stars that are not young against the solar circle's (`profile.planetary_nebula_factor_array`).
"""

import numpy as np
import pytest

from planetgen import tuning
from planetgen.galaxy import component_field, density, nebula_field, profile
from planetgen.generation import component_pass, phenomenon_reach, phenomenon_scatter
from tests import test_satellites as satellite_tests
from tests.test_satellites import E1, EDGE, SHAPE

KIND = phenomenon_scatter.OLD_STAR_KIND
PLACED = satellite_tests._placed()
WITH_SATELLITES = SHAPE._replace(components=component_field.ComponentField(PLACED, EDGE, E1))
SMC, LMC, FORNAX = (next(c for c in PLACED if c.name == name) for name in
                    ("Small Magellanic Cloud", "Large Magellanic Cloud", "Fornax dwarf"))
SUN = (8128.0, 0.0, 0.0)
GALAXY_SEED = bytes(range(16))


def _young_rich_point():
    """A point of the inner disk's plane where the young stars are most of the continuous populations' mix."""
    best = None
    for radius in np.linspace(1500.0, 6000.0, 19):
        for angle in np.linspace(0.0, 2.0 * np.pi, 72, endpoint=False):
            point = (radius * np.cos(angle), radius * np.sin(angle), 0.0)
            value = float(profile.planetary_nebula_factor_array(point[0], point[1], point[2], SHAPE))
            if best is None or value < best[0]:
                best = (value, point)
    return best[1]


# --- The planetary nebulae of the scatters ---------------------------------------------------------------------------


def test_the_array_factor_is_the_profiles_planetary_nebula_factor_at_the_same_point():
    for position in (SUN, (3000.0, 0.0, 0.0), (20000.0, 3000.0, 15000.0), _young_rich_point()):
        here = profile.profile_at(position, SHAPE, EDGE, E1).kind_factors[KIND]
        assert float(profile.planetary_nebula_factor_array(*position, SHAPE)) == pytest.approx(here, rel=1e-9)


def test_a_young_rich_place_makes_fewer_planetary_nebulae_than_an_old_one_in_every_scatter_weight():
    """The young are a few percent of the stars even on an arm's crest, so the factor is a few percent under 1 there
    (the profile's own test has the disk within 3%), and exactly 1 where every star is old."""
    young_rich = _young_rich_point()
    old = (young_rich[0], young_rich[1], 3000.0)
    factor = float(profile.planetary_nebula_factor_array(*young_rich, SHAPE))
    assert 0.9 < factor < 0.99
    assert float(profile.planetary_nebula_factor_array(*old, SHAPE)) == 1.0
    for place, wanted in ((young_rich, factor), (old, 1.0)):
        # The reach weighs the class by the stars' density and the factor alone: no remnant or gas term.
        assert float(phenomenon_reach._weights(KIND, None, *place, SHAPE)) == pytest.approx(
            wanted * float(density.total_density_array(*place, SHAPE)), rel=1e-9)
    # The progress bar's estimate follows the same factor.
    plane = float(np.hypot(*young_rich[:2]))
    assert phenomenon_scatter._expectation_factor(SHAPE, KIND, None, plane, 0.0) == pytest.approx(
        float(profile.planetary_nebula_factor_array(plane, 0.0, 0.0, SHAPE)), rel=1e-9)
    assert phenomenon_scatter._expectation_factor(SHAPE, KIND, None, plane, 3000.0) == 1.0


# --- The gas of a Cloud ----------------------------------------------------------------------------------------------


def test_a_clouds_young_stars_give_it_gas_and_a_dwarf_without_young_stars_has_none():
    assert nebula_field.gas_factor(LMC.center_pc, WITH_SATELLITES) > 0.1
    # (both Clouds sit at the cap now that the disk holds its galaxy's mass, so the tracer is compared, not the factor)
    assert nebula_field._gas_tracer(SMC.center_pc, WITH_SATELLITES) > nebula_field._gas_tracer(LMC.center_pc, WITH_SATELLITES)
    assert nebula_field.gas_factor(FORNAX.center_pc, WITH_SATELLITES) == 0.0
    # Away from every satellite nothing changes: the disk's tracer is the same with the components in or out.
    assert nebula_field.gas_factor(SUN, WITH_SATELLITES) == nebula_field.gas_factor(SUN, SHAPE)


def test_a_cell_inside_a_cloud_shows_every_class_the_profile_allows_and_a_dwarfs_shows_none_of_the_gas_ones():
    cloud = profile.profile_at(LMC.center_pc, WITH_SATELLITES, EDGE, E1)
    for kind in ("supernova-remnant", "molecular-cloud", "nebula"):
        assert 0.0 < cloud.kind_factors[kind] <= 1.0, kind
    assert 0.0 < cloud.kind_factors[KIND] <= 1.0
    assert cloud.gas.cold_fraction > 0.0
    dwarf = profile.profile_at(FORNAX.center_pc, WITH_SATELLITES, EDGE, E1)
    for kind in ("supernova-remnant", "molecular-cloud", "nebula"):
        assert dwarf.kind_factors[kind] == 0.0, kind


def test_the_field_cells_of_a_cloud_hold_molecular_clouds():
    index = nebula_field.cell_index(SMC.center_pc)
    middle = tuple((i + 0.5) * tuning.NEBULA_FIELD_CELL_PC for i in index)
    assert nebula_field.gas_factor(middle, WITH_SATELLITES) > 0.0
    assert sum(len(nebula_field.cell_draw(GALAXY_SEED, WITH_SATELLITES, (index[0] + i, index[1], index[2])))
               for i in range(-2, 3)) > 0
    assert nebula_field.cell_draw(GALAXY_SEED, SHAPE, index) == []


# --- The component pass ----------------------------------------------------------------------------------------------


def test_a_satellites_planetary_nebulae_and_magnetars_follow_its_history_and_a_clusters_are_unchanged():
    rates = {(kind, subtype): rate for kind, subtype, rate in
             component_pass.phenomenon_rates(1.0) + component_pass.satellite_rates(1.0)}
    shares = component_field.population_shares(SMC)
    pn = rates[(KIND, None)]
    factor = component_pass.planetary_nebula_factor(SMC, WITH_SATELLITES)
    assert 0.0 < factor <= 1.0
    assert component_pass._class_mean(SMC, KIND, None, pn, WITH_SATELLITES) == pytest.approx(
        component_pass.systems_of(SMC) * pn * factor)
    magnetar = rates[("neutron-star", phenomenon_scatter.MAGNETAR_SUBTYPE)]
    assert component_pass._class_mean(SMC, "neutron-star", phenomenon_scatter.MAGNETAR_SUBTYPE, magnetar,
                                      WITH_SATELLITES) == pytest.approx(
        component_pass.systems_of(SMC) * magnetar * shares["satellite_young"])
    assert component_pass._class_mean(FORNAX, "neutron-star", phenomenon_scatter.MAGNETAR_SUBTYPE, magnetar,
                                      WITH_SATELLITES) == 0.0
    assert component_pass.planetary_nebula_factor(FORNAX, WITH_SATELLITES) == 1.0     # no young stars: every star is old


def _only_the_satellite_classes(monkeypatch):
    monkeypatch.setattr(component_pass, "phenomenon_rates", lambda _cut: [])


def _remnants(component, shape, seeds, monkeypatch):
    _only_the_satellite_classes(monkeypatch)
    found = []
    for seed in seeds:
        found.extend(row for row in component_pass.scatter_phenomena(
            [component], EDGE, seed, min_mass_solar=tuning.PHENOMENON_MIN_MASS_SOLAR, shape=shape)
            if row[3] == "supernova-remnant")
    return found


def test_a_satellites_remnants_are_kept_where_its_gas_is(monkeypatch):
    rate = dict(((k, s), r) for k, s, r in component_pass.satellite_rates(tuning.PHENOMENON_MIN_MASS_SOLAR))[
        ("supernova-remnant", None)]
    full = component_pass.systems_of(SMC) * rate
    seeds = range(40)
    rows = _remnants(SMC, WITH_SATELLITES, seeds, monkeypatch)
    radii, gas = component_pass.gas_by_radius(SMC, WITH_SATELLITES)
    assert 0.0 < gas[0] <= 1.0 and gas[-1] < gas[0]
    assert rows, "the Small Cloud holds some remnants"
    centre = np.asarray(SMC.center_pc)
    distance = [float(np.linalg.norm(np.asarray(row[5:8]) / 1000.0 - centre)) for row in rows]
    # Kept ones crowd the centre, where the gas is: the median kept radius is under the median of all the satellite's.
    assert np.median(distance) < SMC.half_mass_radius_pc * 1.3
    # The kept count is the candidates (full * seeds) times the mean gas fraction over the profile, well under the whole.
    assert len(rows) < 0.75 * full * len(seeds)


def test_a_satellite_without_gas_has_no_remnants_or_magnetars(monkeypatch):
    assert _remnants(FORNAX, WITH_SATELLITES, range(10), monkeypatch) == []
    _only_the_satellite_classes(monkeypatch)
    rows = list(component_pass.scatter_phenomena(
        [FORNAX], EDGE, 3, min_mass_solar=tuning.PHENOMENON_MIN_MASS_SOLAR, shape=WITH_SATELLITES))
    assert not [row for row in rows if row[3] == "neutron-star"]


def test_the_classes_every_component_has_keep_their_draws(monkeypatch):
    """The satellite's remnants and magnetars come after the others in its stream, so adding them moves none of the
    black holes, neutron stars or planetary nebulae."""
    cut = tuning.PHENOMENON_MIN_MASS_SOLAR
    everything = list(component_pass.scatter_phenomena([SMC], EDGE, 9, min_mass_solar=cut, shape=WITH_SATELLITES))
    generic = [row for row in everything if row[3] not in ("supernova-remnant",)
               and not (row[3] == "neutron-star" and row[4] == phenomenon_scatter.MAGNETAR_SUBTYPE)]
    monkeypatch.setattr(component_pass, "satellite_rates", lambda _cut: [])
    only = list(component_pass.scatter_phenomena([SMC], EDGE, 9, min_mass_solar=cut, shape=WITH_SATELLITES))
    assert generic == only
