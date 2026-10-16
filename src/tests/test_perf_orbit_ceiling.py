# tests/test_perf_orbit_ceiling.py

"""
PERF.103: a ceiling on the rows one orbital-update run moves. A run over the ceiling moves at most the ceiling (the
oldest due of each table first), warns with the backlog, and leaves the rest due for the next run, which takes
them first; nothing is skipped silently.
"""

import pytest

from planetgen.cli import orbits as updateOrbits
from planetgen.db import store
from planetgen.db.store import RowBudget
from tests.test_galactic_motion import SCATTER_RADIUS_PC, _insert_scatter_point, _seed_scatter_grid
from tests.test_update_due import _planet_system

LONG = 1e6   # years: every planet, moon and scatter point is due


def _budget(candidates, ceiling):
    budget = RowBudget(ceiling)
    budget.candidates = dict(candidates)
    budget.allocate()
    return budget


def test_everything_fits_when_the_due_rows_are_under_the_ceiling_or_there_is_no_ceiling():
    assert _budget({"planets": 60, "moons": 40}, 100).limits == {}
    assert _budget({"planets": 600, "moons": 400}, None).limits == {}
    assert _budget({"planets": 600, "moons": 400}, 0).limits == {}


@pytest.mark.parametrize("candidates, ceiling", [
    ({"planets": 1_000_000, "moons": 3_000_000, "comets": 5, "bright_stars": 40_000_000}, 10_000_000),
    ({"planets": 10, "moons": 10, "comets": 10}, 7),
    ({"planets": 1_000_000, "moons": 1, "bright_stars": 0}, 100),
    ({"planets": 9, "moons": 9, "bright_stars": 9, "black_holes": 9}, 3),
])
def test_the_limits_never_add_up_to_more_than_the_ceiling_and_never_starve_a_small_table(candidates, ceiling):
    budget = _budget(candidates, ceiling)
    assert set(budget.limits) == set(candidates)
    assert sum(budget.limits.values()) <= ceiling
    assert sum(budget.limits.values()) == min(ceiling, sum(candidates.values()))
    floor = ceiling // (4 * len(candidates))
    for key, count in candidates.items():
        assert budget.limits[key] <= count
        assert budget.limits[key] >= min(count, floor)


def test_the_rest_of_the_ceiling_is_shared_in_proportion_to_what_is_left_due():
    # a quarter of the ceiling (100) is the floor, 50 each; the other 300 goes 950 : 2950 to what is left due
    budget = _budget({"planets": 1000, "moons": 3000}, 400)
    assert budget.limits == {"planets": 50 + 73, "moons": 50 + 227}


def test_take_counts_a_refused_row_as_deferred_and_a_table_without_a_limit_takes_all():
    budget = _budget({"comets": 5, "planets": 100}, 4)
    taken = [budget.take("comets") for _ in range(budget.limits["comets"] + 3)]
    assert taken.count(True) == budget.limits["comets"] and budget.deferred["comets"] == 3
    free = RowBudget(None)
    assert all(free.take("comets") for _ in range(50)) and free.backlog == 0 and free.warning() is None


def test_the_warning_names_the_backlog_and_where_it_sits():
    budget = _budget({"planets": 100, "moons": 10}, 20)
    budget.defer("planets", 30)
    budget.defer("moons", 5)
    text = budget.warning()
    assert "35 due row(s)" in text and "planets 30" in text and "moons 5" in text and "ceiling of 20" in text


def test_a_run_over_the_ceiling_moves_the_oldest_planets_and_moons_and_the_next_run_takes_the_rest(mysql_config):
    system, cfg = _planet_system()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            store.insert_star_system(conn, system, cfg)
        tick = store.orbit_clock(conn, 1.0 / store.physical_constants.SECONDS_PER_YEAR)
        store.advance_orbital_phases(conn, tick)       # schedules every planet and moon
        before = {name: {r["id"]: r for r in conn.execute(f"SELECT id, next_update_due, epoch_unix FROM {name}").fetchall()}
                  for name in ("planets", "moons")}
        total = sum(len(rows) for rows in before.values())
        assert total >= 3
        ceiling = total - 2
        clock = store.orbit_clock(conn, LONG)
        budget = RowBudget(ceiling)
        budget.plan(conn, clock)
        counts = store.advance_orbital_phases(conn, clock, budget=budget)
        moved = counts["planets"] + counts["moons"]
        assert 0 < moved <= ceiling
        assert budget.deferred.get("planets", 0) + budget.deferred.get("moons", 0) == total - moved
        assert budget.backlog == total - moved
        for name in ("planets", "moons"):
            now = {r["id"]: r for r in conn.execute(f"SELECT id, next_update_due, epoch_unix FROM {name}").fetchall()}
            done = [i for i, r in now.items() if r["epoch_unix"] == pytest.approx(clock.end_unix)]
            left = [i for i in now if i not in done]
            # the oldest due went first, and what is left is still due with its old epoch
            assert all(before[name][i]["next_update_due"] <= before[name][j]["next_update_due"] for i in done for j in left)
            assert all(now[i]["next_update_due"] <= clock.end_unix and now[i]["epoch_unix"] == before[name][i]["epoch_unix"]
                       for i in left)
        again = store.advance_orbital_phases(conn, clock, budget=RowBudget(ceiling))
        assert again["planets"] + again["moons"] == total - moved
    finally:
        conn.close()


def test_scatter_points_over_the_ceiling_wait_with_their_epoch_and_move_from_it_next_time(mysql_config):
    _seed_scatter_grid(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            ids = [_insert_scatter_point(conn, "bright_stars", (SCATTER_RADIUS_PC, float(k), 0.0)) for k in range(5)]
        clock = store.orbit_clock(conn, LONG)
        budget = RowBudget(2)
        budget.plan(conn, clock)
        assert budget.candidates["bright_stars"] == 5 and budget.limits["bright_stars"] == 2
        motion = store.advance_scatter_points(conn, clock, budget=budget)
        assert motion["counts"]["bright_stars"] == 2 and budget.deferred == {"bright_stars": 3}
        rows = {r["id"]: r for r in conn.execute("SELECT id, epoch_unix, next_update_due FROM bright_stars").fetchall()}
        moved = [i for i in ids if rows[i]["epoch_unix"] == pytest.approx(clock.end_unix)]
        waiting = [i for i in ids if i not in moved]
        assert len(moved) == 2 and len(waiting) == 3
        # a point that waits has its own epoch (the run's start), so it is due and moves from there next time
        assert all(rows[i]["epoch_unix"] == pytest.approx(clock.start_unix) and rows[i]["next_update_due"] <= clock.end_unix
                   for i in waiting)
        # the same span again: the two that moved are not due until their next interval, so the three that waited
        # are the candidates, two of them move (oldest first) and one stays
        for moved_now, left in ((2, 1), (1, 0)):
            budget = RowBudget(2)
            budget.plan(conn, clock)
            assert budget.candidates["bright_stars"] == moved_now + left
            assert store.advance_scatter_points(conn, clock, budget=budget)["counts"]["bright_stars"] == moved_now
            assert budget.backlog == left
    finally:
        conn.close()


def test_a_run_with_no_budget_or_room_enough_moves_everything_as_before(mysql_config):
    _seed_scatter_grid(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            for k in range(4):
                _insert_scatter_point(conn, "bright_stars", (SCATTER_RADIUS_PC, float(k), 0.0))
        clock = store.orbit_clock(conn, LONG)
        roomy = RowBudget(1000)
        roomy.plan(conn, clock)
        assert roomy.limits == {}
        assert store.advance_scatter_points(conn, clock, budget=roomy)["counts"]["bright_stars"] == 4
        assert roomy.backlog == 0
    finally:
        conn.close()


def _argv(mysql_config, *extra):
    return ["planetgen.cli.orbits", "--mysql-host", mysql_config.host, "--mysql-port", str(mysql_config.port),
            "--mysql-user", mysql_config.user, "--mysql-password", mysql_config.password or "",
            "--mysql-database", mysql_config.database, *extra]


def test_the_command_warns_with_the_backlog_when_the_ceiling_is_low(mysql_config, monkeypatch, capsys):
    _seed_scatter_grid(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            for k in range(6):
                _insert_scatter_point(conn, "bright_stars", (SCATTER_RADIUS_PC, float(k), 0.0))
    finally:
        conn.close()
    monkeypatch.setattr("sys.argv", _argv(mysql_config, "--row-ceiling", "2"))
    updateOrbits.main()                      # the first run only sets the clock
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            conn.execute("UPDATE orbit_simulation_state SET last_updated_at = NOW() - INTERVAL 400 DAY")
            conn.execute("UPDATE bright_stars SET epoch_unix = NULL, next_update_due = NULL")
    finally:
        conn.close()
    capsys.readouterr()
    updateOrbits.main()
    captured = capsys.readouterr()
    assert "warning:" in captured.err and "due row(s) are over the ceiling of 2" in captured.err
    assert "bright_stars" in captured.err


def test_a_negative_ceiling_is_refused(mysql_config, monkeypatch, capsys):
    monkeypatch.setattr("sys.argv", _argv(mysql_config, "--row-ceiling", "-1"))
    with pytest.raises(SystemExit):
        updateOrbits.main()
    assert "--row-ceiling must be 0" in capsys.readouterr().err
