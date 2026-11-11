# tests/test_sky_charting.py

"""
VIEW.25: charting builds the sky tiers (`generation/sky_charting.py`, `db/sky_tier_ledger.py`).

* the ledger: how deep a cell is complete under the spheres drawn, the chain of touching bands, the spheres a new one is inside;
* a tier drawn over a sphere gives the stars the cell backfill gives for the cells of the sphere, writes no `sector_stats` row, records the sphere, and is not drawn twice;
* a stopped draw goes on from its last batch and ends with the stars of an unbroken one;
* the rings follow the luminosity bands, each cell as it stands.
"""

import json
import math

import pytest

from planetgen import tuning
from planetgen.db import sky_tier_ledger as ledger_module
from planetgen.db import store
from planetgen.db.sky_tier_ledger import Draw, Ledger
from planetgen.galaxy import sky_band, spherical
from planetgen.generation import bright_stars as brightStars
from planetgen.generation import run_plan, sky_charting, sky_tiers
from tests.test_bright_star_scatter import EDGE_PC, SMALL_E_VALUE, THRESHOLD, _scatter_args, _seed_galaxy, _stored_stars

CENTER = (0.0, 0.0, 0.0)


def lum(low, high, center=CENTER, radius=50.0):
    return Draw("luminosity", low, high, center, radius)


def mass(low, high, center=CENTER, radius=50.0):
    return Draw("mass", low, high, center, radius)


# --- the ledger -------------------------------------------------------------------------------------------------------------------------------------

def test_a_cell_is_complete_down_through_the_bands_whose_ranges_touch():
    ledger = Ledger([lum(50.0, 100.0), lum(100.0, 200.0), lum(200.0, None), lum(1000.0, 2000.0)])
    cell = (1, 0, 0)
    assert ledger.luminosity_reached([cell], EDGE_PC, {cell: 250.0}) == {cell: 50.0}          # 250 down through 200, 100, 50
    assert ledger.luminosity_reached([cell], EDGE_PC, {cell: 150.0}) == {cell: 50.0}
    assert ledger.luminosity_reached([cell], EDGE_PC, {cell: 40.0}) == {}                       # already deeper than any band
    assert ledger.luminosity_reached([cell], EDGE_PC, {cell: 3000.0}) == {cell: 50.0}          # the 200 band has no top, so the 1000 to 2000 one inside it adds nothing


def test_a_gap_in_the_bands_stops_the_chain():
    ledger = Ledger([lum(50.0, 100.0), lum(200.0, 300.0)])
    cell = (1, 0, 0)
    assert ledger.luminosity_reached([cell], EDGE_PC, {cell: 300.0}) == {cell: 200.0}          # 100 to 200 was not drawn, so 50 does not count
    assert ledger.luminosity_reached([cell], EDGE_PC, {cell: 100.0}) == {cell: 50.0}


def test_a_cell_outside_a_sphere_is_not_reached_and_the_mass_rings_chain_the_same_way():
    ledger = Ledger([lum(50.0, 100.0, radius=5.0), mass(1.0, 2.0), mass(2.0, 5.0)])
    near = (1, 0, 0)
    far = (6, 0, 0)
    assert spherical.cell_in_sphere(near, CENTER, 50.0, EDGE_PC) and not spherical.cell_in_sphere(far, CENTER, 5.0, EDGE_PC)
    assert ledger.luminosity_reached([far], EDGE_PC, {far: 100.0}) == {}
    assert ledger.mass_reached([near], EDGE_PC, {near: 5.0}) == {near: 1.0}
    assert ledger.mass_reached([near], EDGE_PC, {near: 8.0}) == {}                              # 5 to 8 was not drawn


def test_a_sphere_is_drawn_when_a_recorded_one_of_the_same_band_holds_it():
    ledger = Ledger([lum(50.0, 100.0, center=(0.0, 0.0, 0.0), radius=100.0)])
    assert ledger.drawn("luminosity", 50.0, 100.0, (10.0, 0.0, 0.0), 80.0)
    assert not ledger.drawn("luminosity", 50.0, 100.0, (10.0, 0.0, 0.0), 95.0)
    assert not ledger.drawn("luminosity", 100.0, 200.0, (10.0, 0.0, 0.0), 10.0)
    assert not ledger.drawn("mass", 50.0, 100.0, (10.0, 0.0, 0.0), 10.0)


# --- the limit a run charts to ----------------------------------------------------------------------------------------------------------------------

def test_the_setting_gives_the_sky_limit_a_run_charts_to_and_the_command_line_overrides_it(monkeypatch):
    import argparse
    from planetgen.util import settings
    monkeypatch.delenv("PLANETGEN_SKY_LIMIT", raising=False)
    assert settings.Queue().sky_limit == 7.5                                          # on by default: charting builds the tiers when a system is charted
    assert sky_charting.limit_of(argparse.Namespace()) == 7.5
    monkeypatch.setenv("PLANETGEN_SKY_LIMIT", "6.0")
    assert sky_charting.limit_of(argparse.Namespace(sky_limit=None)) == 6.0
    assert sky_charting.limit_of(argparse.Namespace(sky_limit=8.0)) == 8.0
    assert sky_charting.limit_of(argparse.Namespace(sky_limit=0.0)) is None           # the command line turns it off for a run
    monkeypatch.setenv("PLANETGEN_SKY_LIMIT", "0")
    assert sky_charting.limit_of(argparse.Namespace(sky_limit=None)) is None
    assert sky_charting.limit_of(argparse.Namespace(sky_limit=7.5)) == 7.5


def test_the_test_suite_charts_without_the_sky_tiers():
    from planetgen.util.settings import get_settings
    assert get_settings().queue.sky_limit == 0                                        # tests/conftest.py: a toy galaxy would take minutes to walk a real table


# --- a tier drawn over a sphere ---------------------------------------------------------------------------------------------------------------------

MASS_PASS = ("--mass", "--phenomenon-min-mass", "20")
SPHERE_CENTER = spherical.sector_position_pc(4, 0, 5, EDGE_PC)


@pytest.fixture
def galaxy(mysql_config, monkeypatch):
    """The toy galaxy scattered down to 2,500 Lsun and 20 solar masses, with the observer's cell padding cut to a few sectors so a sphere stays small."""
    _seed_galaxy(mysql_config, e_value=SMALL_E_VALUE * 5)
    run_plan.scatter_bright_stars(_scatter_args(mysql_config, *MASS_PASS))
    monkeypatch.setattr(sky_charting, "CELL_PAD_PC", 2.0 * EDGE_PC)
    conn = store.get_connection(mysql_config)
    try:
        yield conn
    finally:
        conn.close()


def a_tier(ledger="luminosity", low=500.0, high=1000.0, reach=2.0 * EDGE_PC):
    return sky_tiers.Tier(ledger, low, low, high, reach, 0, 0)


def new_stars(conn, before):
    """The bright stars written since `before` (a `_stored_stars` list), sorted."""
    remaining = list(_stored_stars(conn))
    for row in before:
        remaining.remove(row)
    return remaining


def expected_by_cells(conn, tier, center, radius, level, limit):
    context = sky_charting._Context(conn)
    cells = [row[:3] for row in spherical.enumerate_sectors_within_radius(center, radius, EDGE_PC)]
    top = level if tier.high is None else min(tier.high, level)
    return sorted(tuple(row) for row in brightStars.backfill_cells(context.skeleton.shape, cells, EDGE_PC, context.expected, tier.low, top, context.seed,
                                                                   mass_range=(None, limit)))


def test_a_luminosity_tier_writes_the_stars_of_its_sphere_records_it_and_is_not_drawn_twice(galaxy):
    conn = galaxy
    tier = a_tier()
    radius = sky_charting.padded_radius(tier)
    before = _stored_stars(conn)
    done = sky_charting.draw_tier(sky_charting._Context(conn), tier, SPHERE_CENTER)
    assert not done.skipped and done.stars > 0 and done.batches >= 1 and done.radius_pc == radius
    added = new_stars(conn, before)
    limit = store.bright_star_mass_limit(conn)
    assert len(added) == done.stars
    assert added == expected_by_cells(conn, tier, SPHERE_CENTER, radius, THRESHOLD, limit)
    levels = {row["bright_level_sol"] for row in conn.execute("SELECT bright_level_sol FROM sector_stats").fetchall()}
    assert levels == {store.UNTOUCHED_LEVEL}          # the band stops short of the cells' level (2,500 and over), so it is the ledger's, not a cell row's, that holds the depth
    ledger = ledger_module.read(conn)
    (draw,) = ledger.draws
    assert (draw.kind, draw.low, draw.high, draw.radius_pc) == ("luminosity", 500.0, 1000.0, radius) and draw.stars == done.stars
    again = sky_charting.draw_tier(sky_charting._Context(conn), tier, SPHERE_CENTER)
    assert again.skipped and new_stars(conn, before) == added


def test_a_stopped_draw_goes_on_and_ends_with_the_stars_of_an_unbroken_one(galaxy, monkeypatch):
    conn = galaxy
    monkeypatch.setattr(brightStars, "SPHERE_BATCH_CELLS", 4)           # the toy sphere holds few cells with a star: cut it into several batches
    monkeypatch.setattr(brightStars, "SPHERE_BATCH_ROWS", 3)
    tier = a_tier(low=100.0, high=1000.0)
    unbroken_before = _stored_stars(conn)
    sky_charting.draw_tier(sky_charting._Context(conn), tier, SPHERE_CENTER)
    unbroken = new_stars(conn, unbroken_before)
    batches = conn.execute("SELECT batches_done FROM sky_tier_draws").fetchone()["batches_done"]
    assert unbroken and batches >= 2, "the sphere must take two batches for this test"
    conn.execute("DELETE FROM bright_stars")
    conn.execute("DELETE FROM sky_tier_draws")
    conn.commit()
    # the same galaxy's stars again, stopped after the first batch
    class Stop(Exception):
        pass

    def stop(batch, cells, stars):
        raise Stop

    with pytest.raises(Stop):
        sky_charting.draw_tier(sky_charting._Context(conn), tier, SPHERE_CENTER, on_batch=stop)
    conn.rollback()
    assert conn.execute("SELECT complete, batches_done FROM sky_tier_draws").fetchone() == {"complete": 0, "batches_done": 1}
    assert not ledger_module.read(conn).draws                                                      # a half-drawn sphere does not count
    sky_charting.draw_tier(sky_charting._Context(conn), tier, SPHERE_CENTER)
    assert sorted(_stored_stars(conn)) == sorted(unbroken)
    assert conn.execute("SELECT complete FROM sky_tier_draws").fetchone()["complete"] == 1


def test_a_mass_ring_draws_the_cells_as_the_mass_backfill_does_and_records_the_mass(galaxy):
    conn = galaxy
    tier = a_tier("mass", low=5.0, high=8.0)
    radius = sky_charting.padded_radius(tier)
    context = sky_charting._Context(conn)
    before = _stored_stars(conn)
    done = sky_charting.draw_tier(context, tier, SPHERE_CENTER)
    cells = [row[:3] for row in spherical.enumerate_sectors_within_radius(SPHERE_CENTER, radius, EDGE_PC)]
    expected = sorted(tuple(row) for row in brightStars.backfill_mass_cells(context.skeleton.shape, cells, EDGE_PC, context.expected, 5.0, 8.0, THRESHOLD,
                                                                           context.seed, context.edges))
    assert done.stars > 0 and new_stars(conn, before) == expected
    masses = {row["bright_mass_sol"] for row in conn.execute("SELECT bright_mass_sol FROM sector_stats").fetchall()}
    assert masses == {5.0}


def delete_since(conn, first_id):
    conn.execute("DELETE FROM bright_stars WHERE id > ?", (first_id,))
    conn.execute("DELETE FROM sector_stats")
    conn.execute("DELETE FROM sky_tier_draws")
    conn.commit()


def test_a_sphere_then_the_sector_backfill_gives_the_stars_of_the_sector_backfill_alone(galaxy):
    # The consistency the ledger is for: cells inside a drawn sphere are complete down to its floor, so the per-sector backfill that follows draws nothing twice and nothing is
    # missing, whether or not the cell held a star the sphere wrote.
    conn = galaxy
    tier = a_tier(low=1000.0, high=3000.0)                 # reaches up to the scatter's 2,500, so the chain touches
    radius = sky_charting.padded_radius(tier)
    cells = sorted(row[:3] for row in spherical.enumerate_sectors_within_radius(SPHERE_CENTER, radius, EDGE_PC))
    first_id = conn.execute("SELECT MAX(id) AS id FROM bright_stars").fetchone()["id"]
    context = sky_charting._Context(conn)

    sky_charting.draw_tier(context, tier, SPHERE_CENTER)
    run_plan._draw_sector_bands(conn, context.skeleton, cells, 1000.0, context.limits, context.seed)
    both = sorted(_stored_stars(conn))
    levels = {row["bright_level_sol"] for row in conn.execute("SELECT bright_level_sol FROM sector_stats").fetchall()}
    assert levels == {1000.0}                              # every cell of the sphere, with a star or not, now carries the depth

    delete_since(conn, first_id)
    run_plan._draw_sector_bands(conn, context.skeleton, cells, 1000.0, context.limits, context.seed)
    assert sorted(_stored_stars(conn)) == both


# --- the queue's unit -------------------------------------------------------------------------------------------------------------------------------

def test_the_progress_of_a_tier_rises_with_every_cell_even_past_what_was_predicted():
    """OPS.134: the counter of a sky unit never stands still while it works; a tier that runs over its prediction creeps toward the tier's end and does not reach it."""
    steps = [sky_charting.tier_progress(100.0, seen) for seen in range(0, 5000, 7)]
    assert steps == sorted(set(steps)), "strictly rising"
    assert sky_charting.tier_progress(100.0, 45) == 45.0 and sky_charting.tier_progress(100.0, 90) == 90.0
    assert 90.0 < sky_charting.tier_progress(100.0, 91) < sky_charting.tier_progress(100.0, 10 ** 6) < 100.0
    assert sky_charting.tier_progress(0.0, 0) == 0.0 and sky_charting.tier_progress(0.0, 3) < 1.0


def test_the_queue_draws_the_sky_tiers_of_each_observer_cell_once(galaxy, mysql_config, monkeypatch):
    """`_sky_units` makes one unit for each observer cell the ledger does not hold all of; the unit draws its tiers, and a second look finds nothing owed."""
    from planetgen.generation import run_galaxy
    conn = galaxy
    real_plan = sky_charting.plan

    def small_plan(observer_pc, limit_mag, edge_pc, **options):         # a toy table: one band and one ring, a few cells across
        center, tiers = real_plan(observer_pc, limit_mag, edge_pc, **options)
        return center, tiers._replace(luminosity=(a_tier(low=500.0, high=1000.0),), rings=(a_tier("mass", low=5.0, high=8.0),), dropped=())

    monkeypatch.setattr(sky_charting, "plan", small_plan)
    named = {"host": mysql_config.host, "port": mysql_config.port, "user": mysql_config.user, "database": mysql_config.database}
    params = {"database": named, "sky_limit": 7.5}
    units = run_galaxy._sky_units(params, [SPHERE_CENTER, SPHERE_CENTER])       # two sectors of one cell: one unit
    assert len(units) == 1 and units[0][1] >= 1.0 and units[0][0]["sky"]["limit"] == 7.5
    assert units[0][0]["expected_seconds"] >= 0.0          # OPS.134: the unit's own time, for its own limit (the toy tiers predict no visits)
    before = _stored_stars(conn)
    reports = []
    result = run_galaxy.backfill_bright_stars_unit(units[0][0], lambda done, total: reports.append((done, total)))
    conn.commit()
    # OPS.134: the counter moved with the batches and the tiers (it read 0 of 1 for the whole unit), up to every owed tier's visits
    assert len(reports) >= 2 and len({total for _done, total in reports}) == 1
    assert [done for done, _total in reports] == sorted(done for done, _total in reports) and reports[-1][0] == reports[-1][1] > reports[0][0] > 0.0                                                       # a new snapshot: the unit wrote through its own connection
    assert result["sky"]["spheres"] == 2 and result["stars"] == result["sky"]["stars"] == len(new_stars(conn, before))
    assert len(ledger_module.read(conn).draws) == 2
    floors = result["sky"]["floors"]                                  # VIEW.47: the unit reports the floor against distance, whole enough to cross the queue as JSON
    assert [row["distance_pc"] for row in floors] == list(sky_tiers.FLOOR_DISTANCES_PC)
    assert floors == json.loads(json.dumps(floors)) and floors[0]["floor_sol"] < floors[-1]["floor_sol"]
    assert run_galaxy._sky_units(params, [SPHERE_CENTER]) == []


def test_a_charted_cell_is_built_out_to_its_spheres_and_a_half_chart_is_not(galaxy):
    """VIEW.42: `built_to_pc` (the sky route's) reads the ledger: 0 before a chart and with only the luminosity chain, the sphere's reach once the luminosity bands reach the rings' 50 Lsun and the rings reach the lightest mass."""
    conn = galaxy
    top = store.bright_star_mass_limit(conn)
    center = tuple(sky_band.cell_center(sky_band.cell_of(SPHERE_CENTER)))
    tiers = sky_tiers.SkyTiers(7.5, SPHERE_CENTER, EDGE_PC, 0.0, (), (), (), None)
    bands = tiers._replace(luminosity=(a_tier(low=sky_tiers.RING_TOP_SOL, high=THRESHOLD),), rings=(), dropped=())
    toy = tiers._replace(luminosity=bands.luminosity, rings=(a_tier("mass", low=min(sky_tiers.ring_edges()), high=top),), dropped=())
    assert sky_charting.built_to_pc(conn, center) == 0.0
    sky_charting.chart(conn, center, 7.5, tiers=bands, center=center)
    conn.commit()
    assert sky_charting.built_to_pc(conn, center) == 0.0              # the luminosity chain is whole, the rings are not drawn
    sky_charting.chart(conn, center, 7.5, tiers=toy, center=center)
    conn.commit()
    built = sky_charting.built_to_pc(conn, center)
    assert built == pytest.approx(sky_charting.padded_radius(toy.luminosity[0]) - 0.0) and built > 0.0
    far = tuple(c + 10.0 * EDGE_PC for c in center)
    assert sky_charting.built_to_pc(conn, far) == 0.0                         # no drawn sphere holds an observer out there
