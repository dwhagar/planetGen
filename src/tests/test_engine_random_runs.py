"""OPS.245: the proof of the job engine's stage 5. A seed draws a whole run (two to four stages of units made of pieces that commit, one stage written in windows by the spawner, units that fail for good so the run
is held, kills at commits, a spawner that stops in the middle of a window, an end of a stage that stops) and the run is carried through all of it to its end. `pytest -k seed_7` replays one; more seeds with
`PLANETGEN_ENGINE_SOAK_SEEDS=200`. Real MySQL and Redis, the inline launcher and `kill_harness` of OPS.225 (a unit killed at a commit ends as a job process that died of SIGKILL).

What must hold at the end, for every seed:

* the run is `done`, and what its units wrote is exactly what an undisturbed run writes: every piece of every unit once (a piece done twice, or never, is a row twice, or none);
* nothing ran before its needs: no unit of stage n + 1 had an attempt begin before stage n's `release` row was done, and no attempt began on a row that was not released;
* nothing ran twice: each unit row has one attempt that is `done`; a windowed stage has a row for each of its items, once, with no gap; each end row of each stage is done.
"""

import os
import random

import pytest

from planetgen.db import store
from planetgen.queue import commands, held, launcher, ledger, planrows, redisqueue, runrecord, units, wrapper
from tests import inline_queue, kill_harness

SEEDS = range(int(os.environ.get("PLANETGEN_ENGINE_SOAK_SEEDS", "6")))
BROKEN = set()
"""set: The `(stage, unit)` whose piece fails for good while it is in this set (the fault a person fixes before `Retry held`)."""


def _cut(params):
    level, count, window = params["level"], params["count"], int(params.get("window") or params["count"])
    shared = {"level": level, "pieces": params["pieces"], "database": params["database"]}
    first = [({**shared, "n": n}, 1.0) for n in range(min(count, window))]
    later = shared if count > window else None
    return commands.Expansion(first, {"count": count}, total=count, cost=float(count), later=later)


def _window(later, start, count):
    return [({**later, "n": n}, 1.0) for n in range(start, start + count)]


def _pieces(args, cursor):
    return list(range(args["pieces"]))[cursor:]


def _do(args, piece, ctx):
    if (args["level"], args["n"]) in BROKEN:
        raise ValueError(f"unit {args['n']} of stage {args['level']} is broken")
    ctx.db.execute("INSERT INTO proof_rows (stage, unit, piece) VALUES (?, ?, ?)", (args["level"], args["n"], piece))
    ctx.commit_piece({"piece": piece})        # a piece that would write twice if it ran twice commits as finished (OPS.222)
    return {"piece": piece}


class Stop(BaseException):
    """What a stop of the spawner looks like to the code it interrupts."""


def _shape(rng):
    """The run of a seed: the stages' unit counts and pieces, which stage is written in windows and how big a window is."""
    stages = rng.randint(2, 4)
    windowed = rng.randrange(stages)
    shape = [{"count": rng.randint(8, 22) if index == windowed else rng.randint(1, 10), "pieces": rng.randint(1, 3)} for index in range(stages)]
    return shape, windowed, rng.randint(3, 6)


def _expected(shape):
    return {(stage, unit, piece) for stage, one in enumerate(shape) for unit in range(one["count"]) for piece in range(one["pieces"])}


def _pump(conn, run_id, turns=300):
    """What the supervisor does once the run's follower is gone: hands out the ready rows and ends the stages whose end was stopped, until the run is over."""
    for _turn in range(turns):
        inline_queue.run_ready_rows()
        for index in range(len(runrecord.declared_stages(conn, run_id))):
            if units.stage_unended(conn, run_id, index):
                conn.delete(runrecord.ending_key(run_id, index))        # the lease of the attempt that stopped has run out
                units._end_stage(conn, run_id, index, wrapper, sleep=lambda seconds: None)
        if runrecord.run_state(conn, run_id) not in runrecord.LIVE_STATES:
            return


def _ledger_rows(control, conn, run_id):
    engine_run = int(runrecord._s(conn.hget(runrecord.run_key(run_id), planrows.ENGINE_RUN_FIELD)))
    return engine_run, control.execute("SELECT * FROM engine_jobs WHERE run_id = ? ORDER BY id", (engine_run,)).fetchall()


@pytest.fixture
def conn(redis_server):
    return redisqueue.connect(redis_server)


@pytest.mark.timeout(900)        # a seed takes under a minute; a run that never ends (a gap in the windows, a row never released) is a failure, not a hang
@pytest.mark.fresh_database
@pytest.mark.parametrize("seed", SEEDS, ids=lambda seed: f"seed_{seed}")
def test_a_run_drawn_by_a_seed_ends_done_with_every_piece_once_and_every_stage_in_order(mysql_config, conn, monkeypatch, seed):
    rng = random.Random(seed)
    shape, windowed, window = _shape(rng)
    keys = [f"rprobe{index}" for index in range(len(shape))]
    for key in keys:
        monkeypatch.setitem(commands.REGISTRY, key, commands.CommandSpec(key, "A test", pieces=(__name__, "_pieces"), do=(__name__, "_do"), resumes=True, expand=(__name__, "_cut"),
                                                                         window=(__name__, "_window")))
    monkeypatch.setattr(planrows, "window_size", lambda workers: window)
    database = wrapper.database_args(mysql_config)
    steps = [commands.Step(f"s{index}", key, f"Stage {index}", {"level": index, "count": one["count"], "pieces": one["pieces"], "database": database})
             for index, (key, one) in enumerate(zip(keys, shape))]
    galaxy = store.get_connection(mysql_config)
    galaxy.execute("CREATE TABLE IF NOT EXISTS proof_rows (stage INT, unit INT, piece INT)")
    galaxy.execute("DELETE FROM proof_rows")
    galaxy.commit()

    total_pieces = sum(one["count"] * one["pieces"] for one in shape)
    kills = {number: rng.choice(("before", "after")) for number in rng.sample(range(total_pieces), min(total_pieces, rng.randint(0, 5)))}
    plan = kill_harness.Plan(kills)
    kill_harness.install(monkeypatch, plan)
    BROKEN.clear()
    if rng.random() < 0.5:                                            # units that fail for good: the run is held, a person fixes the fault and retries
        BROKEN.update((stage, unit) for stage in rng.sample(range(len(shape)), rng.randint(1, len(shape))) for unit in rng.sample(range(shape[stage]["count"]), 1))
    stops = {"after rows": rng.random() < 0.7, "before release": rng.random() < 0.7}        # the spawner stops once after its rows went in, once after the cursor moved
    real_add, real_release = ledger.add_jobs, planrows._release_window

    def add_then_stop(control, run, jobs):
        real_add(control, run, jobs)
        if stops["after rows"] and any(job.get("item_key") not in (None, "", "0") and int(job["item_key"]) >= window for job in jobs):
            stops["after rows"] = False
            raise Stop("the spawner stopped after its rows went in and before the cursor moved")

    def release_then_stop(conn_, control, run, stage, row_id, arguments):
        if stops["before release"] and arguments["released"] < arguments["next"]:
            stops["before release"] = False
            raise Stop("the spawner stopped after the cursor moved and before the rows were released")
        return real_release(conn_, control, run, stage, row_id, arguments)

    real_spawn = launcher.Launcher.spawn

    def spawn_that_may_die(self):
        try:
            return real_spawn(self)
        except Stop:                    # the supervisor died here: what it had written stands, what it had not committed is rolled back with its connection, and the next turn is a new supervisor's
            return 0

    monkeypatch.setattr(ledger, "add_jobs", add_then_stop)
    monkeypatch.setattr(planrows, "_release_window", release_then_stop)
    monkeypatch.setattr(launcher.Launcher, "spawn", spawn_that_may_die)

    before = set(runrecord.recent_runs(conn, limit=1000))
    try:
        wrapper.run("Random", steps, lambda index, plan: None, lambda index, unit, result: None)
    except (units.RunStopped, Stop):
        pass
    run_id = max(set(runrecord.recent_runs(conn, limit=1000)) - before)         # the ids begin with the second the run started in
    for _round in range(6):
        _pump(conn, run_id)
        state = runrecord.run_state(conn, run_id)
        if state == runrecord.HELD:
            BROKEN.clear()                                              # the fault is fixed
            assert held.retry_held(conn, run_id, who="test") is not None
        elif state == "done":
            break
    assert runrecord.run_state(conn, run_id) == "done", f"seed {seed}: the run ended {runrecord.run_state(conn, run_id)} (shape {shape}, window {window}, kills {kills}, stops {stops}): {runrecord.read_run(conn, run_id).get('error')}"

    assert plan.fired or not kills, f"seed {seed}: the kills {kills} never happened"
    rows = galaxy.execute("SELECT stage, unit, piece, COUNT(*) AS n FROM proof_rows GROUP BY stage, unit, piece").fetchall()
    assert {(row["stage"], row["unit"], row["piece"]) for row in rows} == _expected(shape), f"seed {seed}: what the units wrote is not what an undisturbed run writes"
    assert {row["n"] for row in rows} == {1}, f"seed {seed}: a piece was done twice (kills {kills})"

    control = store.get_control_connection(store.control_mysql_config(mysql_config))
    try:
        engine_run, jobs = _ledger_rows(control, conn, run_id)
        released = {job["stage"]: job for job in jobs if job["command"] == "release" and not job["item_key"]}
        for index, one in enumerate(shape):
            stage = steps[index].key
            unit_rows = [job for job in jobs if job["stage"] == stage and job["item_key"]]
            assert sorted(int(job["item_key"]) for job in unit_rows) == list(range(one["count"])), f"seed {seed}: stage {index} has a gap or a row twice"
            assert {job["state"] for job in unit_rows} == {"done"}, f"seed {seed}: stage {index} has a unit that is not done"
            done_attempts = control.execute("SELECT j.id AS id, COUNT(*) AS n FROM engine_attempts a JOIN engine_jobs j ON j.id = a.job_id WHERE j.run_id = ? AND j.stage = ? AND j.item_key <> '' "
                                            "AND a.outcome = 'done' GROUP BY j.id", (engine_run, stage)).fetchall()
            assert len(done_attempts) == one["count"] and {row["n"] for row in done_attempts} == {1}, f"seed {seed}: a unit of stage {index} did not end done exactly once"
            assert released[stage]["state"] == "done", f"seed {seed}: stage {index} was never released"
            if index:
                earlier = released[steps[index - 1].key]["updated_at"]
                first = control.execute("SELECT MIN(a.started_at) AS at FROM engine_attempts a JOIN engine_jobs j ON j.id = a.job_id WHERE j.run_id = ? AND j.stage = ?", (engine_run, stage)).fetchone()["at"]
                assert first >= earlier, f"seed {seed}: a unit of stage {index} began at {first}, before stage {index - 1} was released at {earlier}"
        spawn = [job for job in jobs if job["command"] == planrows.SPAWN]
        assert all(job["state"] == "done" for job in spawn), f"seed {seed}: a spawn row did not end done"
        assert len(spawn) == sum(one["count"] > window for one in shape), f"seed {seed}: a stage with more units than a window has one spawn row, the others none"
    finally:
        control.close()
        galaxy.close()
