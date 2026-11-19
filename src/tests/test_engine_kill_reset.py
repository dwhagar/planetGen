"""OPS.225, slice 4 (the real-MySQL kill tests, the reset and the site reopening): a galaxy reset with units killed at drawn commits leaves the empty galaxy an undisturbed reset leaves, with the Object ID counters kept (DB.20); the reopening of the site after a map rebuild, killed, ends with the site open (`tests/kill_harness.py`).

Commands: reset-galaxy, end-map-rebuild-maintenance. Replay a seed with `pytest -k seed_3`."""

import random

import pytest

from planetgen.db import store
from planetgen.generation import run_plan
from tests import inline_queue, kill_harness
from tests.test_engine_kill_scatter import SEEDS, _scatter_all

pytestmark = [pytest.mark.fresh_database, pytest.mark.slow]


def _rows_and_counters(config):
    """The number of rows in every table of the database but the engine's own ledgers (they hold the runs, which the reset does not take back), and the Object ID counters (a reset keeps these; the `sdel` stamp it leaves for the page caches is not one)."""
    conn = store.get_connection(config)
    try:
        tables = [next(iter(row.values())) for row in conn.execute("SHOW TABLES").fetchall()]
        counts = {table: conn.execute(f"SELECT COUNT(*) AS n FROM `{table}`").fetchone()["n"] for table in tables if not table.startswith("engine_")}
        counters = sorted(repr(tuple(row.values())) for row in conn.execute("SELECT * FROM id_counters").fetchall() if "sdel" not in repr(tuple(row.values())))
        return counts, counters
    finally:
        conn.close()


def _kill_at_the_steps(monkeypatch, plan):
    """The wipe is a table at a time with no `commit` (a TRUNCATE ends its own transaction), so a step of it is the point a kill lands on: before or after a table is wiped, or the tidying."""
    from planetgen.cli import reset as reset_galaxy
    kill_harness.install(monkeypatch, plan)
    for name in ("_wipe_table", "_tidy_after_wipe"):
        real = getattr(reset_galaxy, name)

        def stepped(*args, _real=real, **kwargs):
            number, plan.seen = plan.seen, plan.seen + 1
            mode = plan.kills.get(number)
            if mode is None or "reset" in plan.units:
                return _real(*args, **kwargs)
            plan.units.add("reset")
            plan.fired.append((number, mode, "reset"))
            if mode == "after":
                _real(*args, **kwargs)
            plan.dead = True
            raise kill_harness.Killed(f"step {number} {mode}")

        monkeypatch.setattr(reset_galaxy, name, stepped)


def _reset(config):
    from planetgen.cli import reset as reset_galaxy
    with inline_queue.inline_queue():
        assert reset_galaxy.reset_database(config, assume_yes=True)


@pytest.mark.parametrize("seed", SEEDS, ids=lambda seed: f"seed_{seed}")
def test_a_reset_with_killed_units_leaves_the_empty_galaxy_an_undisturbed_reset_leaves(mysql_config, monkeypatch, seed):
    import _pytest.monkeypatch
    _scatter_all(mysql_config)
    kept = _rows_and_counters(mysql_config)[1]
    counter = kill_harness.Plan()
    with _pytest.monkeypatch.MonkeyPatch.context() as patch:
        _kill_at_the_steps(patch, counter)
        _reset(mysql_config)
    expected, commits = _rows_and_counters(mysql_config), counter.seen
    assert commits > 0 and not any(expected[0].get(table) for table in ("bright_stars", "phenomenon_scatter", "sectors")), "the undisturbed reset emptied the galaxy"
    assert expected[1] == kept, "the reset kept the Object ID counters"
    _scatter_all(mysql_config)
    rng = random.Random(str(seed))
    kills = {} if seed == "control" else {number: rng.choice(("before", "after")) for number in rng.sample(range(commits), min(commits, rng.randint(1, 4)))}
    plan = kill_harness.Plan(kills)
    _kill_at_the_steps(monkeypatch, plan)
    _reset(mysql_config)
    assert plan.fired or seed == "control", f"seed {seed}: no kill happened (commits {commits}, plan {kills})"
    found = _rows_and_counters(mysql_config)
    assert found[0] == expected[0], f"seed {seed}: kills {plan.fired} left rows behind"
