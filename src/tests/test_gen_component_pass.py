# tests/test_gen_component_pass.py

"""
GEN.324: the clusters and the satellite galaxies place their own bright and massive stars and phenomena
(`generation/component_pass.py`): a count from the component's systems and the share its population keeps, positions
from its profile, the tier and the cell from where each star lands, a stream of its own per component, and the same
share the sector fill already takes off the component's cells, so no pre-placed star is drawn twice.
"""

import collections
import math

import numpy as np
import pytest

from planetgen import tuning
from planetgen.db import store
from planetgen.galaxy import component_field, components, scatter_limits, spherical
from planetgen.galaxy.scatter_limits import ScatterLimits
from planetgen.generation import bright_stars as brightStars
from planetgen.generation import component_pass, phenomenon_scatter, run_plan
from planetgen.physics import component_profiles, king
from tests import test_gen_scatter_tiers as tiers
from tests.test_bright_star_scatter import EXTENTS
from tests.test_phenomenon_scatter import _plan_args as _phenomenon_args

EDGE = tiers.EDGE_PC
E_VALUE = 6.3
FLOOR = 300.0
"""A floor the cluster's giants reach: a few percent of an old population's stars are this bright."""
SEED = 17


def _cluster(name="NGC Test", mass=3.0e6, half=5.0, w0=7, center=(2000.0, 500.0, 300.0)):
    return components.Component(
        kind="cluster", name=name, source="generated", center_pc=center, velocity_kms=(0.0, 0.0, 0.0), epoch_unix=0.0,
        mass_solar=mass, half_mass_radius_pc=half, profile="king", age_gy=11.0, metallicity_fe_h=-1.5, seed=1,
        king_w0=w0, tidal_radius_pc=king.tidal_radius_pc(w0, half))


def _satellite(name="Test Cloud", mass=2.0e8, half=900.0, center=(6000.0, -4000.0, 2500.0)):
    return components.Component(
        kind="satellite", name=name, source="generated", center_pc=center, velocity_kms=(0.0, 0.0, 0.0), epoch_unix=0.0,
        mass_solar=mass, half_mass_radius_pc=half, profile="exponential", age_gy=6.0, metallicity_fe_h=-0.8, seed=2,
        tidal_radius_pc=component_profiles.cut_radius_pc("exponential", half))


def _pass(floor=FLOOR):
    fractions = brightStars.band_fractions(floor)
    return fractions, [scatter_limits.PassBand(floor, None, None)]


def _stars(found, seed=SEED, floor=FLOOR, **kwargs):
    fractions, bands = _pass(floor)
    return list(component_pass.scatter_stars(found, EDGE, seed, fractions, bands, **kwargs))


def _uniform(rng_seed, count):
    return np.random.default_rng(rng_seed).random((3, count))


@pytest.mark.parametrize("make", [_cluster, _satellite])
def test_the_positions_follow_the_profile_about_the_components_centre(make):
    component = make()
    points = component_pass.positions_mpc(component, _uniform(1, 60000)) / 1000.0
    offsets = points - np.asarray(component.center_pc)[:, None]
    radius = np.sqrt((offsets ** 2).sum(axis=0))
    table = component_profiles.table_for(component)
    assert float(np.mean(radius < component.half_mass_radius_pc)) == pytest.approx(0.5, abs=0.01)
    assert radius.max() <= component.tidal_radius_pc * 1.001
    scale = component_profiles.scale_radius_pc(component)
    for fraction in (0.1, 0.9):
        assert float(np.mean(radius < scale * table.radius_x(fraction))) == pytest.approx(fraction, abs=0.01)
    assert np.abs(offsets.mean(axis=1)).max() < 0.03 * radius.mean()          # isotropic: no direction preferred


def test_the_count_is_the_systems_times_the_share_the_population_keeps():
    fractions, _bands = _pass()
    component = _cluster()
    wanted = component_pass.systems_of(component) * fractions["cluster"]
    assert wanted > 200
    counts = [len(_stars([component], seed=seed)) for seed in range(8)]
    assert np.mean(counts) == pytest.approx(wanted, abs=5.0 * math.sqrt(wanted / len(counts)))
    assert component_pass.expected_stars([component], fractions) == pytest.approx(wanted)
    shares = component_field.population_shares(_satellite())
    assert component_pass.expected_stars([_satellite()], fractions) == pytest.approx(
        component_pass.systems_of(_satellite()) * sum(share * fractions[name] for name, share in shares.items()))


def test_a_stars_row_sits_in_its_cell_and_is_a_star_of_its_components_population():
    rows = _stars([_cluster(), _satellite()])
    assert len(rows) > 100
    columns = dict(zip(store.BRIGHT_STAR_COLUMNS, rows[0]))
    assert columns["population"] in component_field.POPULATIONS_IN_THE_GATE
    populations = {row[6] for row in rows}
    assert "cluster" in populations and populations <= set(component_field.POPULATIONS_IN_THE_GATE)
    assert "satellite_young" in populations       # a Cloud with gas places its massive young stars (GEN.526)
    for row in rows:
        address = spherical.sector_address_at(tuple(value / 1000.0 for value in row[3:6]), EDGE)
        assert address == tuple(row[:3])
        assert row[12] >= FLOOR * brightStars.constants.SOLAR_LUMINOSITY * 0.999


def test_the_same_seed_gives_the_same_stars_and_a_component_does_not_depend_on_the_others():
    first, second = _cluster(), _satellite()
    both = _stars([first, second])
    assert both == _stars([first, second]) and both != _stars([first, second], seed=SEED + 1)
    by_population = lambda rows, population: [row for row in rows if row[6] == population]
    assert by_population(both, "cluster") == by_population(_stars([first]), "cluster")
    for population in component_field.SATELLITE_POPULATIONS:
        assert by_population(both, population) == by_population(
            _stars([second, _cluster("Another", mass=3.0e5)]), population)
    assert _stars([_cluster("Renamed")]) != _stars([first])         # the stream is the component's own


def test_sectors_already_filled_get_no_star():
    rows = _stars([_cluster()])
    skipped = {tuple(row[:3]) for row in rows[:40]}
    kept = _stars([_cluster()], skip_addresses=skipped)
    # The stars are drawn after the candidates are placed, so a skipped cell moves the draws of the rows after it.
    assert [tuple(row[:3]) for row in rows if tuple(row[:3]) not in skipped] == [tuple(row[:3]) for row in kept]


def test_each_stars_tier_is_the_one_of_the_place_it_lands_in():
    skeleton = tiers._skeleton()
    limits = ScatterLimits.of(skeleton, tiers.FLOOR, tiers.MASS)
    tiered = scatter_limits.luminosity_pass(limits)
    by_tier = brightStars.tier_fractions_of(tiered)
    fractions, bands = brightStars.pass_of(tiered, tiers.FLOOR, None, None)
    component = _cluster(mass=6.0e6, center=(1000.0, 0.0, 6000.0))
    rows = list(component_pass.scatter_stars([component], EDGE, SEED, fractions, bands, tiered=tiered))
    assert rows
    assert {limits.tier_at(row[0], row[1]) for row in rows} == {scatter_limits.HALO}
    for row in rows:
        tier = limits.tier_at(row[0], row[1])
        assert row[12] >= bands[tier].min_luminosity_sol * brightStars.constants.SOLAR_LUMINOSITY * 0.999
    wanted = component_pass.systems_of(component) * max(item["cluster"] for item in by_tier)
    assert len(rows) <= wanted * 1.3 + 20


def test_the_stars_placed_in_a_cell_are_the_share_the_fill_takes_off_it():
    """No double count: a cell's expected pre-placed stars are its systems times the share the pass keeps, which is what
    `FillContext.bright_share` removes from the cell's own draw."""
    component = _cluster(mass=1.0e6, half=4.0, w0=8, center=(0.0, 0.0, 0.0))
    field = component_field.ComponentField([component], EDGE, E_VALUE)
    heavy = collections.Counter(tuple(row[:3]) for seed in range(12) for row in _stars([component], seed=seed))
    address = heavy.most_common(1)[0][0]
    centre = spherical.sector_position_pc(*address, EDGE)
    systems = field.relative_density(centre) * E_VALUE
    share = brightStars.placed_star_fraction(FLOOR, None, "cluster")
    wanted = systems * share
    assert wanted > 20
    found = heavy[address] / 12
    assert found == pytest.approx(wanted, rel=0.35)


def test_the_phenomena_of_a_component_are_its_systems_times_each_classs_rate():
    cluster = _cluster(mass=4.0e6)
    cut = 1.0
    rates = {(kind, subtype): rate for kind, subtype, rate in component_pass.phenomenon_rates(cut)}
    assert ("supernova-remnant", None) not in rates and ("neutron-star", phenomenon_scatter.MAGNETAR_SUBTYPE) not in rates
    assert ("black-hole", "stellar") in rates and ("neutron-star", None) in rates
    rows = [list(component_pass.scatter_phenomena([cluster], EDGE, seed, min_mass_solar=cut)) for seed in range(6)]
    counts = collections.Counter(phenomenon_scatter.class_label(row[3], row[4]) for chunk in rows for row in chunk)
    for (kind, subtype), rate in rates.items():
        wanted = component_pass.systems_of(cluster) * rate * 6
        if wanted > 30:
            assert counts[phenomenon_scatter.class_label(kind, subtype)] == pytest.approx(wanted, abs=5.0 * math.sqrt(wanted))
    for row in rows[0]:
        assert spherical.sector_address_at(tuple(value / 1000.0 for value in row[5:8]), EDGE) == tuple(row[:3])
        assert row[12] is None or row[12] >= cut
    assert len({row[11] for row in rows[0]}) == len(rows[0])
    assert rows[0] == list(component_pass.scatter_phenomena([cluster], EDGE, 0, min_mass_solar=cut))
    skipped = {tuple(row[:3]) for row in rows[0][:30]}
    kept = list(component_pass.scatter_phenomena([cluster], EDGE, 0, skip_addresses=skipped, min_mass_solar=cut))
    assert kept == [row for row in rows[0] if tuple(row[:3]) not in skipped]


@pytest.mark.db
def test_a_plan_scatter_places_the_components_stars_and_phenomena(mysql_config):
    from tests.test_bright_star_scatter import E_VALUE as e_value, EDGE_PC as edge_pc, SHAPE as shape
    found = [_cluster(mass=1.0e8, half=30.0, w0=5, center=(150.0, -90.0, 60.0)),
             _satellite(mass=1.0e8, half=40.0, center=(150.0, -90.0, 60.0))]
    store.save_galaxy_shape(shape, edge_pc=edge_pc, outer_shell_index=8, expected_system_count_at_density_1=e_value,
                            config=mysql_config, cell_star_ceiling=1e12, galaxy_seed=bytes(range(16)))
    store.replace_galaxy_shells(EXTENTS, config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        store.replace_galaxy_components(conn, found)
        conn.commit()
    finally:
        conn.close()
    stars = run_plan.scatter_bright_stars(_phenomenon_args(mysql_config, "--bright-stars-only", "--workers", "1",
                                                           "--bright-star-min-luminosity", "2500"))
    phenomena = run_plan.scatter_phenomena(_phenomenon_args(mysql_config, "--workers", "1"))
    conn = store.get_connection(mysql_config)
    try:
        populations = {row["population"]: row["n"] for row in conn.execute(
            "SELECT population, COUNT(*) AS n FROM bright_stars GROUP BY population").fetchall()}
        placed = conn.execute("SELECT COUNT(*) AS n FROM phenomenon_scatter").fetchone()["n"]
    finally:
        conn.close()
    assert populations.get("cluster", 0) > 50 and populations.get("satellite", 0) > 50
    assert stars["total"] >= sum(populations.values())
    assert phenomena["total"] == placed > 0
