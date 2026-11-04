# tests/test_gen_component_pass.py

"""
GEN.324: the clusters and the satellite galaxies place their own bright and massive stars and phenomena
(`generation/component_pass.py`): a count from the component's systems and the share its population keeps, positions
from its profile, the tier and the cell from where each star lands, a stream of its own per component, and the same
share the sector fill already takes off the component's cells, so no pre-placed star is drawn twice.
"""

import collections
import datetime
import math

import numpy as np
import pytest

from planetgen.db import store
from planetgen.galaxy import remnant_density, component_field, components, scatter_limits, spherical
from planetgen.galaxy.scatter_limits import ScatterLimits
from planetgen.generation import bright_stars as brightStars
from planetgen.generation import component_pass, phenomenon_scatter, run_plan
from planetgen.physics import component_profiles, king
from planetgen.util import draw
from tests import test_gen_scatter_tiers as tiers
from tests.test_bright_star_scatter import EXTENTS
from tests.test_phenomenon_scatter import _scatter_args as _phenomenon_args

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


SMALL_SATELLITE_MASS = 2.0e7
"""A Cloud a tenth the default's mass: it still places ten thousand stars, and a test that draws it five times takes a minute and not twelve."""

YOUNG_SATELLITE_MASS = 5.0e7
"""The least Cloud mass at which the gas-rich Cloud places young massive stars at this floor (2e7 draws none of them)."""


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
    rows = _stars([_cluster(), _satellite(mass=YOUNG_SATELLITE_MASS)])
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
    first, second = _cluster(), _satellite(mass=SMALL_SATELLITE_MASS)
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


def _tiered_pass():
    skeleton = tiers._skeleton()
    limits = ScatterLimits.of(skeleton, tiers.FLOOR, tiers.MASS)
    tiered = scatter_limits.luminosity_pass(limits)
    fractions, bands = brightStars.pass_of(tiered, tiers.FLOOR, None, None)
    return tiered, fractions, bands


@pytest.mark.parametrize("component", [
    _cluster(mass=6.0e6, center=(1000.0, 0.0, 6000.0)), _cluster(mass=6.0e6, center=(10.0, 0.0, 5.0)),
    _satellite(mass=5.0e8)])
def test_gen678_the_stars_a_component_is_expected_to_place_are_the_stars_its_tiers_keep(component):
    """GEN.678: a unit was expected to place about 23,000 stars and placed 1,000 to 3,000: the expectation counted every
    candidate at the share of the largest tier, but a component lies in the tier of its place, which keeps less."""
    tiered, fractions, bands = _tiered_pass()
    rows = list(component_pass.scatter_stars([component], EDGE, SEED, fractions, bands, tiered=tiered))
    kept = component_pass.kept_stars(component, fractions, tiered, EDGE)
    candidates = dict(component_pass.star_weights([component], fractions, tiered))[component]
    assert kept <= candidates * 1.0000001
    assert len(rows) == pytest.approx(kept, abs=4.0 * math.sqrt(kept) + 3.0)
    assert dict(component_pass.kept_weights([component], fractions, tiered, EDGE)) == {component: kept}


def test_gen678_a_component_in_a_tier_that_keeps_none_is_expected_to_place_none_though_its_candidates_are_many():
    tiered, fractions, bands = _tiered_pass()
    inner = _cluster(mass=6.0e6, center=(10.0, 0.0, 5.0))
    assert component_pass.kept_stars(inner, fractions, tiered, EDGE) == 0.0
    assert dict(component_pass.star_weights([inner], fractions, tiered))[inner] > 100.0
    assert component_pass.kept_stars(inner, brightStars.band_fractions(FLOOR), None, EDGE) == pytest.approx(
        dict(component_pass.star_weights([inner], brightStars.band_fractions(FLOOR)))[inner]), \
        "an untiered pass keeps the share it was given everywhere"


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
        if (kind, subtype) == ("black-hole", "stellar"):      # the failed share of a globular cluster's old, metal-poor stars
            wanted *= remnant_density.component_factor({"cluster": 1.0})
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
    from tests.test_bright_star_scatter import E_VALUE as e_value, EDGE_PC as edge_pc, SEED, SHAPE as shape
    found = [_cluster(mass=1.0e8, half=30.0, w0=5, center=(150.0, -90.0, 60.0)),
             _satellite(mass=1.0e8, half=40.0, center=(150.0, -90.0, 60.0))]
    store.save_galaxy_shape(shape, edge_pc=edge_pc, outer_shell_index=8, expected_system_count_at_density_1=e_value,
                            config=mysql_config, cell_star_ceiling=1e12, galaxy_seed=SEED)
    store.replace_galaxy_shells(EXTENTS, config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        store.replace_galaxy_components(conn, found)
        conn.commit()
    finally:
        conn.close()
    stars = run_plan.scatter_bright_stars(_phenomenon_args(mysql_config, "--mass", "--workers", "1",
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


def test_components_are_cut_into_tasks_by_their_expected_rows_alone():
    huge, big, small = _satellite("Huge", mass=1.0e9), _satellite("Big", mass=2.0e8), _cluster("Small", mass=1.0e5)
    weighted = [(small, 10.0), (huge, 250.0), (_cluster("Dead"), 0.0), (big, 100.0), (_cluster("Other"), 30.0)]
    tasks = component_pass.parts(weighted, 100.0)
    assert [(round(expected), [part[1:] for part in task]) for expected, task in tasks] == [
        (83, [["Huge", 0, 3]]), (83, [["Huge", 1, 3]]), (83, [["Huge", 2, 3]]),
        (10, [["Small", 0, 1]]), (100, [["Big", 0, 1]]), (30, [["Other", 0, 1]])], "a large component in slices, the rest bundled"
    assert sum(expected for expected, _task in tasks) == pytest.approx(10.0 + 250.0 + 100.0 + 30.0)
    assert max(expected for expected, _task in tasks) <= 100.0, "no task holds more than the size"
    bundled = component_pass.parts([(_cluster(f"C{index}"), 30.0) for index in range(7)], 100.0)
    assert [round(expected) for expected, _task in bundled] == [90, 90, 30]
    assert component_pass.parts(weighted, 100.0) == tasks, "the same rows are the same tasks"


def test_a_task_names_its_components_and_their_slices_and_a_missing_component_is_said():
    found = [_cluster("One"), _satellite("Two"), _cluster("Three")]
    chosen, slices = component_pass.chosen(found, [["satellite", "Two", 1, 4], ["cluster", "One", 0, 1]])
    assert [component.name for component in chosen] == ["Two", "One"]
    assert slices == {("satellite", "Two"): (1, 4), ("cluster", "One"): (0, 1)}
    with pytest.raises(ValueError, match="Gone"):
        component_pass.chosen(found, [["cluster", "Gone", 0, 1]])


def test_a_component_in_one_slice_keeps_its_stars_and_its_slices_are_draws_of_its_own():
    cloud = _satellite(mass=5.0e5)
    whole = _stars([cloud])
    assert _stars([cloud], slices={("satellite", cloud.name): (0, 1)}) == whole
    first, second = (_stars([cloud], slices={("satellite", cloud.name): (index, 2)}) for index in range(2))
    assert first and second and first != second and first != whole
    # Together the slices are a draw of the whole component: their Poisson means add to its mean.
    assert len(first) + len(second) == pytest.approx(len(whole), abs=6.0 * math.sqrt(len(whole)))
    again = _stars([cloud], slices={("satellite", cloud.name): (1, 2)})
    assert again == second, "a slice does not depend on the other slices or on the worker that draws it"


def test_the_phenomena_of_a_component_in_slices_are_draws_of_its_own_that_add_up_to_its_whole():
    cut = 1.0
    cluster = _cluster(mass=4.0e6)
    whole = list(component_pass.scatter_phenomena([cluster], EDGE, SEED, min_mass_solar=cut))
    assert list(component_pass.scatter_phenomena(
        [cluster], EDGE, SEED, min_mass_solar=cut, slices={("cluster", cluster.name): (0, 1)})) == whole
    slices = [list(component_pass.scatter_phenomena(
        [cluster], EDGE, SEED, min_mass_solar=cut, slices={("cluster", cluster.name): (index, 3)}))
        for index in range(3)]
    assert all(slices) and slices[0] != slices[1]
    assert sum(len(rows) for rows in slices) == pytest.approx(len(whole), abs=6.0 * math.sqrt(len(whole)))
    assert slices[2] == list(component_pass.scatter_phenomena(
        [cluster], EDGE, SEED, min_mass_solar=cut, slices={("cluster", cluster.name): (2, 3)}))


def _seed_components(mysql_config, found):
    from tests.test_bright_star_scatter import E_VALUE as e_value, EDGE_PC as edge_pc, SEED as galaxy_seed, SHAPE as shape
    store.save_galaxy_shape(shape, edge_pc=edge_pc, outer_shell_index=8, expected_system_count_at_density_1=e_value,
                            config=mysql_config, cell_star_ceiling=1e12, galaxy_seed=galaxy_seed)
    store.replace_galaxy_shells(EXTENTS, config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        store.replace_galaxy_components(conn, found)
        conn.commit()
    finally:
        conn.close()


def _stored(mysql_config, table):
    conn = store.get_connection(mysql_config)
    try:
        return sorted(tuple(str(value) for key, value in row.items()
                            if key != "id" and not isinstance(value, datetime.datetime))
                      for row in conn.execute(f"SELECT * FROM {table}").fetchall())
    finally:
        conn.close()


@pytest.mark.fresh_database
def test_a_large_component_is_a_unit_per_slice_that_reports_its_progress_and_adds_up_at_any_order(
        mysql_config, monkeypatch):
    import json

    from planetgen.queue import wrapper
    monkeypatch.setattr(component_pass, "STAR_SLICE_STARS", 100.0)
    monkeypatch.setattr(run_plan, "COMPONENT_REPORT_ROWS", 4)
    monkeypatch.setattr(run_plan, "CHANNEL_INTERVAL_SECONDS", 0.0)
    _seed_components(mysql_config, [_satellite("Big Cloud", mass=1.0e7, half=40.0, center=(150.0, -90.0, 60.0)),
                                    _cluster("Small One", mass=1.0e6, half=30.0, w0=5, center=(160.0, -90.0, 60.0)),
                                    _cluster("Small Two", mass=1.0e6, half=30.0, w0=5, center=(170.0, -90.0, 60.0))])
    params = {"database": wrapper.database_args(mysql_config), "min_luminosity_sol": 2500.0, "mass_limit_sol": 8.0,
              "seed": 77}
    expansion = run_plan.expand_scatter_luminosity(params)
    shape = run_plan._ScatterShape.from_json(expansion.plan)
    assert expansion.plan == json.loads(json.dumps(expansion.plan))
    assert {key: expansion.plan[key] for key in shape.to_json()} == shape.to_json()
    named = [parts for _key, parts in sorted(shape.component_parts.items())]
    assert len(named) > 2 and sum(1 for task in named if task[0][1] == "Big Cloud") == max(
        part[3] for task in named for part in task if part[1] == "Big Cloud") > 1
    assert all(len(task) == 1 for task in named if task[0][1] == "Big Cloud")
    assert sum(shape.tasks[key][0][2] for key in shape.component_parts) == pytest.approx(
        shape.expected[run_plan.COMPONENT_TASK_SHELL])
    assert all(shape.tasks[key][0][2] <= 100.0 * 1.0001 for key in shape.component_parts)
    mine = [(args, cost) for args, cost in expansion.units if args["component_parts"]]
    assert len(mine) == len(named)

    def run_all(units):
        conn = store.get_connection(mysql_config)
        try:
            conn.execute("DELETE FROM bright_stars")
            conn.commit()
        finally:
            conn.close()
        reports = []
        for args, _cost in units:
            run_plan.scatter_shells_unit(json.loads(json.dumps(args)), lambda done, estimate: reports.append(
                (args["key"], done, estimate)))
        return _stored(mysql_config, "bright_stars"), reports

    stars, reports = run_all(mine)
    assert len(stars) > 60
    heaviest = max(mine, key=lambda unit: unit[1])[0]["key"]
    seen = [(done, estimate) for key, done, estimate in reports if key == heaviest]
    assert len(seen) >= 3 and seen == sorted(seen), "a component's stars report their progress as they are drawn"
    assert all(estimate >= done for done, estimate in seen)
    assert run_all(list(reversed(mine)))[0] == stars, "the stars do not depend on the order the units run in"


@pytest.mark.fresh_database
def test_a_large_components_phenomena_are_a_unit_per_slice_that_add_up_at_any_order(mysql_config, monkeypatch):
    import json

    from planetgen.queue import wrapper
    monkeypatch.setattr(component_pass, "PHENOMENON_SLICE_ROWS", 40000.0)
    _seed_components(mysql_config, [_cluster("Rich", mass=1.0e7, half=30.0, w0=5, center=(150.0, -90.0, 60.0)),
                                    _satellite("Cloud", mass=1.0e7, half=40.0, center=(150.0, -90.0, 60.0))])
    params = {"database": wrapper.database_args(mysql_config), "min_mass_solar": 8.0, "seed": 77}
    expansion = run_plan.expand_scatter_phenomena(params)
    shape = run_plan._PhenomenaShape.from_json(expansion.plan)
    # The plan also says how many outline shells and what mass limit the stage's clock reads (PERF.212, `outline_shells`, `min_mass_solar`): not the shape's.
    plan = {key: value for key, value in expansion.plan.items() if key not in ("filled", "outline_shells", "min_mass_solar")}
    assert shape.to_json() == plan == json.loads(json.dumps(plan))
    assert len(shape.component_parts) > 2
    assert all(shape.tasks[key][0][2] <= 40000.0 * 1.0001 for key in shape.component_parts)
    mine = [(args, cost) for args, cost in expansion.units if args["component_parts"]]
    assert len(mine) == len(shape.component_parts)

    def run_all(units):
        conn = store.get_connection(mysql_config)
        try:
            conn.execute("DELETE FROM phenomenon_scatter")
            conn.commit()
        finally:
            conn.close()
        for args, _cost in units:
            run_plan.scatter_phenomena_unit(json.loads(json.dumps(args)), lambda done, estimate: None)
        return _stored(mysql_config, "phenomenon_scatter")

    rows = run_all(mine)
    assert len(rows) > 100 and run_all(list(reversed(mine))) == rows


def test_a_spreads_tiers_are_looked_up_together_and_are_each_bands_own():
    """PERF.146: the pass asks for the tiers of all its places at once; each star still gets the tier `tier_at` gives
    its band."""
    limits = ScatterLimits.of(tiers._skeleton(), tiers.FLOOR, tiers.MASS)
    tiered = scatter_limits.luminosity_pass(limits)
    fractions, bands = brightStars.pass_of(tiered, tiers.FLOOR, None, None)
    component = _satellite(mass=1.0e8, half=3000.0, center=(1500.0, 0.0, 500.0))
    asked = []

    def counted(places):
        asked.append(list(places))
        return limits.tiers_of_places(places)

    by_tier = brightStars.tier_fractions_of(tiered)
    spots = list(component_pass.component_spots(component, EDGE, by_tier, counted, draw.Stream("x"), set()))
    assert len(spots) > 200
    assert len(asked) <= 3                                # one lookup per population
    assert len({spot[-1] for spot in spots}) > 1                # the baseline's and the halo's
    for shell, band, _wedge, _point, _population, tier in spots:
        assert tier == limits.tier_at(shell, band)
