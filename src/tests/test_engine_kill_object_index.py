"""OPS.225, slice 5 (the real-MySQL kill tests, the warm-ups that write): `warm-object-index` killed in the middle of a kind (before or after any write or sweep of the index) leaves, once it runs again, an index that holds every object as the tables do (`objectindex.check`: nothing missing, different or extra) with each row once.

The read-only commands of the slice (the `check-*` family, `check-math`, `check-datasets`) write nothing a kill could leave half done; the page warm-ups write only caches that a run writes again. Replay a seed with `pytest -k seed_3`."""

import random

import pytest

from planetgen.db import objectindex, store
from planetgen.queue import commands, wrapper
from tests import inline_queue, kill_harness
from tests.test_engine_kill_population import SEEDS, _stored_galaxy

pytestmark = [pytest.mark.fresh_database, pytest.mark.slow]


def _warm(config, mode="build"):
    step = commands.Step("a", "warm-object-index", "A", {"database": config.database, "index_mode": mode})
    with inline_queue.inline_queue():
        wrapper.run("Warm", [step], lambda index, plan: None, lambda index, unit, result: None)


def _counts(config):
    conn = store.get_connection(config)
    try:
        report = objectindex.check(conn)
        return report, conn.execute("SELECT COUNT(*) AS n, COUNT(DISTINCT kind, uid) AS d FROM object_index").fetchone()
    finally:
        conn.close()


def _kill_at_the_writes(monkeypatch, plan):
    """The index writes with its own connection (a commit it makes itself), so the points a kill lands on are its writes and sweeps."""
    kill_harness.install(monkeypatch, plan)
    for name in ("_write", "_sweep"):
        real = getattr(objectindex, name)

        def stepped(*args, _real=real, **kwargs):
            number, plan.seen = plan.seen, plan.seen + 1
            mode = plan.kills.get(number)
            if mode is None or "index" in plan.units:
                return _real(*args, **kwargs)
            plan.units.add("index")
            plan.fired.append((number, mode, "index"))
            if mode == "after":
                _real(*args, **kwargs)
            plan.dead = True
            raise kill_harness.Killed(f"step {number} {mode}")

        monkeypatch.setattr(objectindex, name, stepped)


@pytest.mark.parametrize("seed", SEEDS, ids=lambda seed: f"seed_{seed}")
def test_an_object_index_warmed_with_a_killed_unit_holds_every_object_once_when_it_runs_again(mysql_config, monkeypatch, seed):
    _stored_galaxy(mysql_config)
    counter = kill_harness.Plan()
    with monkeypatch.context() as patch:
        _kill_at_the_writes(patch, counter)
        _warm(mysql_config)
    report, rows = _counts(mysql_config)
    steps = counter.seen
    assert steps > 0 and rows["n"] > 0 and all(not (found["missing"] or found["extra"] or found["different"]) for found in report.values()), "the undisturbed warm-up indexed every object"
    conn = store.get_connection(mysql_config)
    try:
        conn.execute("DELETE FROM object_index")
        conn.commit()
    finally:
        conn.close()
    rng = random.Random(str(seed))
    kills = {} if seed == "control" else {number: rng.choice(("before", "after")) for number in rng.sample(range(steps), min(steps, rng.randint(1, 4)))}
    plan = kill_harness.Plan(kills)
    _kill_at_the_writes(monkeypatch, plan)
    _warm(mysql_config)
    assert plan.fired or seed == "control", f"seed {seed}: no kill happened (steps {steps}, plan {kills})"
    found, again = _counts(mysql_config)
    assert all(not (kind["missing"] or kind["extra"] or kind["different"]) for kind in found.values()), f"seed {seed}: kills {plan.fired} left the index different from the tables: {found}"
    assert (again["n"], again["d"]) == (rows["n"], rows["n"]) and rows["d"] == rows["n"], f"seed {seed}: each object has one row"
