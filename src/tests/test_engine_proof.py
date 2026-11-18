"""OPS.225, slice 1: the proof of the job engine's stage 3 for what is built (OPS.215 to OPS.218), as seeded scenarios that replay exactly (`pytest -k seed_7`).

* Two attempts never both commit: attempts of several jobs begin and commit in a drawn order; a commit lands if and only if its attempt is the job's newest, the data of every other
  commit is rolled back, and a job is done in the galaxy only by its newest attempt.
* A finished stage advances once: the end of a stage (conclusion, finish, the opening of the next stage) is killed at drawn points and retried; each step is done exactly once and the
  stage is released exactly once, however many retries it took.
* A held run resumes to the same result (slice 2): units of a two-stage run fail for good, the run is held and the second stage does not open; after the fault is fixed and Retry held the run ends done
  with every unit's data equal to an undisturbed run's, and a unit that was done before the hold ran only once.
"""

import random
import types

import pytest

from planetgen.db import store
from planetgen.queue import attempts, commands, ledger, planrows, redisqueue, runrecord, units, wrapper

SEEDS = range(10)


class Killed(BaseException):
    """What a kill looks like to the code it interrupts: nothing in the engine may catch it by accident."""


def _control(mysql_config):
    return store.get_control_connection(store.control_mysql_config(mysql_config))


@pytest.mark.fresh_database
@pytest.mark.parametrize("seed", SEEDS, ids=lambda seed: f"seed_{seed}")
def test_two_attempts_never_both_commit(mysql_config, seed):
    rng = random.Random(seed)
    control, galaxy = _control(mysql_config), store.get_connection(mysql_config)
    try:
        galaxy.execute("CREATE TABLE IF NOT EXISTS proof_probe (job BIGINT, attempt INT, step INT)")
        galaxy.execute("DELETE FROM proof_probe")
        galaxy.commit()
        run = ledger.add_run(control, "generate", f"proof-{seed}")
        ledger.add_jobs(control, run, [{"stage": "s", "command": "c", "item_key": str(n)} for n in range(3)])
        control.commit()
        rows = [row["id"] for row in control.execute("SELECT id FROM engine_jobs WHERE run_id = ? ORDER BY id", (run,)).fetchall()]
        held = {row: [] for row in rows}                  # the attempts each job has been handed, oldest first
        newest, done_by = {}, {}
        landed = set()
        for step in range(60):
            row = rng.choice(rows)
            if not held[row] or rng.random() < 0.25:     # a retry starts while the attempts before it may still be alive
                held[row].append(attempts.Attempt(control, row, ledger.begin_attempt(control, row, galaxy), mysql_config))
                newest[row] = held[row][-1].number
                done_by.pop(row, None)                    # a new attempt takes the job back: what an older attempt finished is no longer done
                continue
            attempt = rng.choice(held[row])
            finished = rng.random() < 0.3
            galaxy.execute("INSERT INTO proof_probe VALUES (?, ?, ?)", (row, attempt.number, step))
            try:
                attempt.commit(galaxy, cursor_key=str(step), done=finished)
            except attempts.Replaced:
                assert attempt.number != newest[row], f"seed {seed}: the newest attempt of job {row} was refused"
                continue
            assert attempt.number == newest[row], f"seed {seed}: attempt {attempt.number} of job {row} committed after attempt {newest[row]} took it over"
            landed.add((row, attempt.number, step))
            if finished:
                done_by[row] = attempt.number
        stored = {(r["job"], r["attempt"], r["step"]) for r in galaxy.execute("SELECT job, attempt, step FROM proof_probe").fetchall()}
        assert stored == landed, f"seed {seed}: the galaxy holds data of a commit that was refused, or lost data of one that landed"
        for row in rows:
            item = ledger.item_state(galaxy, row)
            assert item["done"] == (row in done_by), f"seed {seed}: job {row} is done in the galaxy only if its newest attempt finished it"
    finally:
        galaxy.execute("DROP TABLE IF EXISTS proof_probe")
        control.close()
        galaxy.close()


@pytest.fixture
def conn(redis_server):
    return redisqueue.connect(redis_server)


@pytest.mark.parametrize("seed", SEEDS, ids=lambda seed: f"seed_{seed}")
def test_a_finished_stage_advances_exactly_once_however_often_its_end_is_killed(conn, monkeypatch, seed):
    rng = random.Random(seed)
    run_id = runrecord.start_run(conn, "T", [runrecord.Stage("a", "A", units=1), runrecord.Stage("b", "B", units=1)])
    runrecord.begin_stage(conn, run_id, 0)
    assert runrecord.unit_done(conn, run_id, 0, 0).advance
    done = {"conclude": 0, "finish": 0, "open": 0}
    kills = {"conclude": rng.randint(0, 2), "finish": rng.randint(0, 2), "open": rng.randint(0, 3)}        # kills each step suffers before it completes
    asked = {"conclude": 0, "finish": 0, "open": 0}

    def fault(step):
        asked[step] += 1
        if asked[step] <= kills[step]:
            raise Killed(step)

    def conclude(c, r, i, label=None):
        fault("conclude")
        done["conclude"] += 1
        planrows.mark_ended(c, r, i, planrows.CONCLUDE)
        return 0.0

    def finish(c, r, i, took=0.0):
        fault("finish")
        done["finish"] += 1
        planrows.mark_ended(c, r, i, planrows.FINISH)

    def open_next(c, r, i):
        fault("open")
        done["open"] += 1

    monkeypatch.setattr(wrapper, "conclude_stage_job", conclude)
    monkeypatch.setattr(wrapper, "finish_stage_job", finish)
    monkeypatch.setattr(units, "enqueue_open", open_next)
    for retry in range(40):
        try:
            units._end_stage(conn, run_id, 0, wrapper, sleep=lambda seconds: None)
        except Killed:
            conn.delete(runrecord.ending_key(run_id, 0))      # the killed attempt's lease runs out
            continue
        except ConnectionError:
            continue
        if not units.stage_unended(conn, run_id, 0):
            break
    assert not units.stage_unended(conn, run_id, 0), f"seed {seed}: the end was never finished"
    assert done == {"conclude": 1, "finish": 1, "open": 1}, f"seed {seed}: {done} with kills {kills}"
    units._end_stage(conn, run_id, 0, wrapper, sleep=lambda seconds: None)          # a late retry finds the stage released and does nothing
    assert done == {"conclude": 1, "finish": 1, "open": 1}
    assert set(planrows.ended(conn, run_id, 0)) == {planrows.CONCLUDE, planrows.FINISH, planrows.RELEASE}


# --- a held run resumes to the same result (OPS.225, slice 2) ------------------------------------------------------------------------------------------------------------------

def _two_stage_run(monkeypatch, broken, stored, counts):
    """A run of two stages of test units; a unit whose number is in `broken` fails for good while `broken` holds it. What a unit makes is stored under (stage, unit), as a unit's data is in the galaxy."""
    from tests import test_work_units as work

    def make(stage):
        def handler(args, report):
            counts[(stage, args["n"])] = counts.get((stage, args["n"]), 0) + 1
            if args["n"] in broken:
                raise ValueError(f"unit {args['n']} of {stage} is broken")
            stored[(stage, args["n"])] = args["n"] ** 2 + len(stage)
            return {"value": stored[(stage, args["n"])]}
        return handler

    handlers = {"alpha": make("alpha"), "beta": make("beta")}
    work._handlers(monkeypatch, **handlers)
    monkeypatch.setattr(commands, "unit_function", lambda kind: handlers[kind])
    return [work._step(5, key="alpha", kind="alpha"), work._step(4, key="beta", kind="beta")]


@pytest.mark.parametrize("seed", SEEDS, ids=lambda seed: f"seed_{seed}")
def test_a_held_run_that_is_retried_makes_what_an_undisturbed_run_makes(conn, monkeypatch, seed):
    from planetgen.queue import held
    from tests import inline_queue
    rng = random.Random(seed)
    undisturbed, ran = {}, {}
    steps = _two_stage_run(monkeypatch, set(), undisturbed, ran)
    wrapper.run("Undisturbed", steps, lambda index, plan: None, lambda index, unit, result: None)
    assert len(undisturbed) == 9 and set(ran.values()) == {1}

    broken = set(rng.sample(range(4), rng.randint(1, 3)))           # units of the first stage that fail for good (the second has four, the first five)
    disturbed, ran = {}, {}
    steps = _two_stage_run(monkeypatch, broken, disturbed, ran)
    with pytest.raises(units.RunStopped):
        wrapper.run("Disturbed", steps, lambda index, plan: None, lambda index, unit, result: None)
    run_id = next(found for found in runrecord.recent_runs(conn) if runrecord.read_run(conn, found)["title"] == "Disturbed")
    assert runrecord.run_state(conn, run_id) == runrecord.HELD, f"seed {seed}: the run is held for a person"
    assert ("beta", 0) not in disturbed, f"seed {seed}: the second stage does not open while the first holds units"
    done_before = dict(disturbed)
    broken.clear()                                                 # the fault is fixed
    assert held.retry_held(conn, run_id, who="test") is not None
    for _turn in range(60):
        inline_queue.run_ready_rows()
        if runrecord.run_state(conn, run_id) == "done":
            break
    assert runrecord.run_state(conn, run_id) == "done", f"seed {seed}: the retried run did not end done"
    assert disturbed == undisturbed, f"seed {seed}: the retried run made another result than an undisturbed one"
    for key, value in done_before.items():
        assert ran[key] == 1, f"seed {seed}: unit {key} was done and ran again after the retry"
