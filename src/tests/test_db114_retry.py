# tests/test_db114_retry.py

"""
DB.114: the backfills, the special phenomena and the sky charting write their stars and Object IDs through `db.contention.run_with_retries`, the way the
scatter units do since DB.113. Each batch is one whole transaction (locks, draw, stars with their IDs, levels, ledger, commit) at READ COMMITTED that a
deadlock (1213) rolls back and writes again with nothing carried over.

A deadlock is made the way it happens: the batch's rows are written (and given their IDs), then MySQL kills the transaction. The tests check that the
rollback takes the rows back, that the retry writes exactly the stars an undisturbed run writes, and that a progress tick is not given twice.
"""

import pymysql
import pytest

from planetgen.db import contention, store
from planetgen.generation import run_galaxy, run_plan, scatter_companions, sky_charting
from tests.test_bright_star_scatter import EDGE_PC, SMALL_E_VALUE, THRESHOLD, _seed_galaxy, _stored_stars  # noqa: F401
from tests.test_sky_charting import SPHERE_CENTER, a_tier, galaxy  # noqa: F401


class Deadlock:
    """`insert_bright_stars` that writes its rows and then fails with a deadlock on its first `times` calls."""

    def __init__(self, real, times=1):
        self.real, self.times, self.calls = real, times, 0

    def __call__(self, conn, rows, *args, **kwargs):
        self.calls += 1
        written = self.real(conn, rows, *args, **kwargs)
        if self.calls <= self.times and rows:
            raise pymysql.err.OperationalError(1213, "Deadlock found when trying to get lock; try restarting transaction")
        return written


@pytest.fixture
def deadlock(monkeypatch):
    """Installs a `Deadlock` on the scatter writer all the DB.114 callers use; no waiting between tries."""
    monkeypatch.setattr(contention.time, "sleep", lambda seconds: None)

    def install(times=1):
        wrapper = Deadlock(scatter_companions.insert_bright_stars, times)
        monkeypatch.setattr(scatter_companions, "insert_bright_stars", wrapper)
        return wrapper

    return install


def _unit(conn):
    skeleton = store.get_galaxy_shape(conn)
    limits, seed = run_plan._scatter_limits_and_seed(conn, skeleton)
    return skeleton, limits, seed


@pytest.mark.fresh_database
def test_the_mass_backfill_after_a_deadlock_writes_the_stars_of_an_undisturbed_run(mysql_config, deadlock, monkeypatch):
    from planetgen.galaxy.spherical import sector_position_pc
    from tests import first_fill
    from tests.test_bright_star_scatter import EXTENTS, SEED, SHAPE
    first_fill.without_the_deepest_floor(monkeypatch)    # the rows below are the backfill's alone (GEN.633)
    store.save_galaxy_shape(SHAPE, edge_pc=EDGE_PC, outer_shell_index=8, expected_system_count_at_density_1=100.0,
                            config=mysql_config, cell_star_ceiling=1e12, galaxy_seed=SEED)
    store.replace_galaxy_shells(EXTENTS, config=mysql_config)
    args = run_galaxy._default_generation_args(config=mysql_config)
    args.num_systems = 0
    address = (3, 0, 4)
    center = sector_position_pc(*address, EDGE_PC)
    run_galaxy.generate_and_save_sector_at(args, address, center, EDGE_PC)
    conn = store.get_connection(mysql_config)
    try:
        skeleton, limits, seed = _unit(conn)
        todo = run_galaxy.backfill_targets(conn, skeleton, limits, [center])
        assert todo
        ticks = []
        wrapper = deadlock(times=2)
        sectors, stars = run_plan._draw_sector_masses(conn, skeleton, todo, limits, seed, on_sector=lambda: ticks.append(1))
        assert wrapper.calls >= 3 and stars > 0, "the writer should have been killed twice and then landed"
        after_deadlocks = _stored_stars(conn)
        assert len(after_deadlocks) == stars                                  # nothing from a rolled-back try is left
        assert len(ticks) == len(todo)                                        # a retried chunk does not tick its sectors twice
        store.clear_bright_stars(conn)
        wrapper.times = 0
        assert run_plan._draw_sector_masses(conn, skeleton, todo, limits, seed) == (sectors, stars)
        assert _stored_stars(conn) == after_deadlocks
    finally:
        conn.close()


def test_the_luminosity_backfill_after_a_deadlock_lands_each_star_once(galaxy, deadlock):
    conn = galaxy
    skeleton, limits, seed = _unit(conn)
    addresses = sorted({row[:3] for row in _stored_stars(conn)})[:6]
    assert addresses
    before = _stored_stars(conn)
    wrapper = deadlock(times=1)
    ticks = []
    _sectors, stars = run_plan._draw_sector_bands(conn, skeleton, addresses, 500.0, limits, seed, on_sector=lambda: ticks.append(1))
    assert wrapper.calls == 2 and stars > 0, "the writer should have been killed once and then landed"
    after = _stored_stars(conn)
    assert len(after) == len(before) + stars and len(set(after)) == len(after)       # nothing doubled by the retry
    assert all(row in after for row in before)
    assert len(ticks) == len(addresses)
    assert conn.execute("SELECT COUNT(*) AS n FROM bright_stars WHERE uid IS NULL").fetchone()["n"] == 0
    assert conn.execute("SELECT COUNT(*) AS n, COUNT(DISTINCT uid) AS d FROM bright_stars").fetchone()["n"] == len(after)


def test_a_sky_tier_batch_after_a_deadlock_writes_the_stars_of_an_undisturbed_draw(galaxy, deadlock):
    conn = galaxy
    tier = a_tier()
    base = _stored_stars(conn)
    seeds = [row[store.BRIGHT_STAR_COLUMNS.index("seed")] for row in base]
    sky_charting.draw_tier(sky_charting._Context(conn), tier, SPHERE_CENTER)
    clean = _stored_stars(conn)
    assert len(clean) > len(base)
    conn.execute("DELETE FROM bright_stars WHERE seed NOT IN (%s)" % ", ".join(["%s"] * len(seeds)), seeds) if seeds else None
    conn.execute("DELETE FROM sky_tier_draws")
    conn.commit()
    assert _stored_stars(conn) == base
    wrapper = deadlock(times=2)
    done = sky_charting.draw_tier(sky_charting._Context(conn), tier, SPHERE_CENTER)
    assert wrapper.calls >= 3 and done.stars > 0
    assert _stored_stars(conn) == clean
    assert conn.execute("SELECT complete, stars FROM sky_tier_draws").fetchone() == {"complete": 1, "stars": done.stars}


def test_the_tick_helper_gives_each_sector_one_tick_whatever_the_number_of_tries():
    ticks, ticked = [], [0]
    for _try in range(3):
        for index in range(4):
            run_plan._tick(lambda: ticks.append(1), ticked, index)
    assert len(ticks) == 4
