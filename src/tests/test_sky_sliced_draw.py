# tests/test_sky_sliced_draw.py

"""
OPS.253 part B: the sky draw of a backfill stage runs as several units side by side (`run_galaxy._sky_units`, `sky_charting.chart_part`, `sky_tier_ledger.conclude_parts`).

* the rows of a sphere's walk are dealt out to the parts: disjoint, together the whole walk, and the cells they yield are the whole walk's cells;
* the parts of a stage draw the stars the draw one cell at a time draws, and leave every cell's level, its mass and the ledger as that draw leaves them;
* the order the parts run in changes nothing, and a part stopped part way goes on from its last batch;
* nothing the parts do counts in the ledger until the stage ends, and a tier a cell's sphere already holds is not owed again.
"""

import pytest

from planetgen.db import sky_tier_ledger as ledger_module
from planetgen.db import store
from planetgen.generation import bright_stars as brightStars
from planetgen.generation import run_galaxy, sky_charting
from planetgen.queue import commands
from tests.test_bright_star_scatter import EDGE_PC, _stored_stars
from tests.test_sky_charting import SPHERE_CENTER, a_tier, delete_since, galaxy  # noqa: F401 (the fixture)

PARTS = 3


@pytest.fixture
def toy_tiers(monkeypatch):
    """A toy table of three tiers (two luminosity bands that reach the cells' level, one ring), a few sectors across, in place of the real one."""
    real_plan = sky_charting.plan

    def small_plan(observer_pc, limit_mag, edge_pc, **options):
        center, tiers = real_plan(observer_pc, limit_mag, edge_pc, **options)
        return center, tiers._replace(luminosity=(a_tier(low=1000.0, high=3000.0), a_tier(low=500.0, high=1000.0)),
                                      rings=(a_tier("mass", low=5.0, high=8.0),), dropped=())

    monkeypatch.setattr(sky_charting, "plan", small_plan)
    monkeypatch.setattr(run_galaxy, "SKY_PARTS", PARTS)
    monkeypatch.setattr(brightStars, "SPHERE_BATCH_CELLS", 4)           # several batches to a part
    monkeypatch.setattr(brightStars, "SPHERE_BATCH_ROWS", 3)


def named(mysql_config):
    return {"host": mysql_config.host, "port": mysql_config.port, "user": mysql_config.user, "database": mysql_config.database}


def observers_of():
    """Two observer cells side by side."""
    other = (SPHERE_CENTER[0] + 50.0, SPHERE_CENTER[1], SPHERE_CENTER[2])
    return [SPHERE_CENTER, other]


def cell_state(conn):
    """What every cell's row holds: its luminosity level and its mass."""
    return sorted((row["shell_index"], row["band_index"], row["wedge_index"], row["bright_level_sol"], row["bright_mass_sol"])
                  for row in conn.execute("SELECT shell_index, band_index, wedge_index, bright_level_sol, bright_mass_sol FROM sector_stats").fetchall())


def spheres(conn):
    return sorted((d.kind, d.low, d.high, d.center, d.radius_pc) for d in ledger_module.read(conn).draws)


def run_stage(conn, mysql_config, order=None, stop=None):
    """What the backfill stage does for the sky: cut the units, run them in `order`, conclude. Returns `(units, detail)`."""
    params = {"database": named(mysql_config), "sky_limit": 7.5}
    units, detail = run_galaxy._sky_units(params, observers_of())
    for index in order if order is not None else range(len(units)):
        commands.run_alone("scatter-backfill", units[index][0])
        conn.commit()
    run_galaxy.conclude_backfill(params, detail, [], 0.0)
    conn.commit()
    return units, detail


def test_the_rows_of_a_walk_are_dealt_out_to_the_parts_and_the_parts_are_the_whole_walk(galaxy, toy_tiers):
    context = sky_charting._Context(galaxy)
    tier = a_tier(low=500.0, high=1000.0)
    radius = sky_charting.padded_radius(tier)
    args = (context.skeleton.shape, SPHERE_CENTER, radius, EDGE_PC, context.expected, 500.0, 1000.0, context.seed)
    whole = [a for addresses, _memo in brightStars.luminosity_sphere_batches(*args, mass_range=context.walk_mass_range) for a in addresses]
    found = [[a for addresses, _memo in brightStars.luminosity_sphere_batches(*args, mass_range=context.walk_mass_range, part=(index, PARTS)) for a in addresses]
             for index in range(PARTS)]
    assert len(whole) > PARTS and sorted(whole) == sorted(a for part in found for a in part)
    assert len(set(whole)) == len(whole) and not any(set(found[i]) & set(found[j]) for i in range(PARTS) for j in range(i))
    assert all(found), "every part has cells to draw in this sphere"


def test_the_parts_of_a_stage_draw_what_one_cell_at_a_time_draws(galaxy, mysql_config, toy_tiers):
    """The census: the stars (every column but the ids), every cell's level and mass, and the spheres the ledger holds are the same."""
    conn = galaxy
    first_id = conn.execute("SELECT MAX(id) AS id FROM bright_stars").fetchone()["id"]
    before = _stored_stars(conn)
    for observer in observers_of():
        sky_charting.chart(conn, observer, 7.5)
    conn.commit()
    serial_stars, serial_cells, serial_spheres = _stored_stars(conn), cell_state(conn), spheres(conn)
    assert len(serial_stars) > len(before) and len(serial_spheres) == 6 and serial_cells
    delete_since(conn, first_id)
    assert _stored_stars(conn) == before and not spheres(conn)

    units, detail = run_stage(conn, mysql_config, order=[2, 0, 1])          # the order the parts run in changes nothing
    assert len(units) == PARTS and detail["sky_cells"] == 2 and len(detail["sky_spheres"]) == 6
    assert _stored_stars(conn) == serial_stars
    assert cell_state(conn) == serial_cells
    assert spheres(conn) == serial_spheres
    assert conn.execute("SELECT COUNT(*) AS n FROM sky_tier_draws WHERE part_count > 0").fetchone()["n"] == 0, "the parts' rows are gone with the stage"


def test_nothing_counts_in_the_ledger_until_the_stage_ends(galaxy, mysql_config, toy_tiers):
    conn = galaxy
    params = {"database": named(mysql_config), "sky_limit": 7.5}
    units, detail = run_galaxy._sky_units(params, observers_of())
    for args, _cost in units:
        commands.run_alone("scatter-backfill", args)
        conn.commit()
    assert not ledger_module.read(conn).draws
    assert conn.execute("SELECT COUNT(*) AS n FROM sky_tier_draws WHERE part_count > 0 AND complete = 1").fetchone()["n"] == PARTS * 6
    run_galaxy.conclude_backfill(params, detail, [], 0.0)
    conn.commit()
    assert len(ledger_module.read(conn).draws) == 6
    run_galaxy.conclude_backfill(params, detail, [], 0.0)           # a second conclusion (a retry) records nothing more
    conn.commit()
    assert len(ledger_module.read(conn).draws) == 6


def test_a_part_stopped_part_way_goes_on_and_the_stage_ends_with_the_same_stars(galaxy, mysql_config, toy_tiers):
    conn = galaxy
    first_id = conn.execute("SELECT MAX(id) AS id FROM bright_stars").fetchone()["id"]
    run_stage(conn, mysql_config)
    unbroken, cells = _stored_stars(conn), cell_state(conn)
    delete_since(conn, first_id)

    class Stop(Exception):
        pass

    params = {"database": named(mysql_config), "sky_limit": 7.5}
    units, detail = run_galaxy._sky_units(params, observers_of())
    real = sky_charting.draw_tier
    stopped = []

    def draw_and_stop(context, tier, center, **options):
        def stop(batch, cells_, stars):
            if not stopped and batch >= 1 and options.get("part") == (1, PARTS):
                stopped.append(batch)
                raise Stop
        return real(context, tier, center, **{**options, "on_batch": stop})

    sky_charting.draw_tier = draw_and_stop
    try:
        with pytest.raises(Stop):
            commands.run_alone("scatter-backfill", units[1][0])
    finally:
        sky_charting.draw_tier = real
    conn.rollback()
    assert stopped, "the part must have been stopped after a first batch"
    for args, _cost in units:
        commands.run_alone("scatter-backfill", args)
        conn.commit()
    run_galaxy.conclude_backfill(params, detail, [], 0.0)
    conn.commit()
    assert _stored_stars(conn) == unbroken and cell_state(conn) == cells


def test_a_tier_the_ledger_holds_is_not_owed_again_and_a_held_cell_makes_no_unit(galaxy, mysql_config, toy_tiers):
    conn = galaxy
    params = {"database": named(mysql_config), "sky_limit": 7.5}
    units, detail = run_galaxy._sky_units(params, [SPHERE_CENTER])
    assert len(detail["sky_spheres"]) == 3
    held = detail["sky_spheres"][0]
    draw = ledger_module.Draw(held["kind"], held["low"], held["high"], tuple(held["center"]), held["radius"])
    draw_id, _done = ledger_module.begin(conn, draw)
    ledger_module.finish(conn, draw_id, 0, 0.0)
    conn.commit()
    _units, again = run_galaxy._sky_units(params, [SPHERE_CENTER])
    assert [(s["kind"], s["low"]) for s in again["sky_spheres"]] == [(s["kind"], s["low"]) for s in detail["sky_spheres"][1:]]
    for sphere in detail["sky_spheres"][1:]:
        item = ledger_module.Draw(sphere["kind"], sphere["low"], sphere["high"], tuple(sphere["center"]), sphere["radius"])
        draw_id, _done = ledger_module.begin(conn, item)
        ledger_module.finish(conn, draw_id, 0, 0.0)
    conn.commit()
    assert run_galaxy._sky_units(params, [SPHERE_CENTER]) == ([], {})
    assert units and all(unit["sky"]["part"] == [index, PARTS] for index, (unit, _cost) in enumerate(units))


def test_the_stage_command_concludes_the_spheres(galaxy):
    from planetgen.queue import commands
    assert commands.conclude_function("scatter-backfill") is run_galaxy.conclude_backfill
    assert run_galaxy.conclude_backfill({"database": {}}, {"sectors": 5}, [], 0.0) == {}          # no sky in the stage: nothing recorded, nothing connected to
