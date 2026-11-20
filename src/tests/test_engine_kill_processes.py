"""OPS.236, the kill tests of stage 4 with real processes: a real template, units forked from it that take their attempts on a real control database, and a real SIGKILL of the template, of a job, of the
process that is the supervisor and of a job as the kernel would kill it. Each ends the same way: the unit is done once, each of its pieces was finished once whatever was killed, and its attempts say what
happened to each process.

The units are `tests/engine_kill_runners.unit`: pieces written to a file, so what a kill cut off is seen. The sector and scatter kill tests (`test_engine_kill_sectors.py`, `test_engine_kill_scatter.py`, OPS.225)
kill a job at a commit of a real command with the inline launcher; these kill the processes around it."""

import json
import os
import signal
import subprocess
import sys
import time
from types import SimpleNamespace

import pytest

from planetgen.db import store
from planetgen.queue import allowlist, claims, jobhost, launcher, ledger, memorygate, redisqueue

R = "tests.engine_kill_runners"
SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
HOST = memorygate.Host(16000.0, 12000.0, 8, 0.2)


def _named(config):
    return {"host": config.host, "port": config.port, "user": config.user, "password": config.password, "database": config.database}


@pytest.fixture
def engine(mysql_config, monkeypatch, tmp_path):
    monkeypatch.setenv(allowlist.EXTRA_ENV, f"{R}:*")
    control = store.get_control_connection(mysql_config, ensure_schema=True)
    control.execute("DELETE FROM engine_runs")
    control.commit()
    forker = jobhost.TemplateSet(warm_sets={"plain": []}, ready_seconds=60)
    made = SimpleNamespace(control=control, forker=forker, config=mysql_config, tmp=tmp_path, children=[])
    try:
        yield made
    finally:
        for child in made.children:
            try:
                os.killpg(child.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            child.wait()
        forker.stop()
        control.execute("DELETE FROM engine_runs")
        control.commit()
        control.close()


def _units(engine, count=1, pieces=6, piece_seconds=0.3, **faults):
    """A run whose ready rows are `count` units of `pieces` pieces each; returns `(run, [row ids], [marks files])`."""
    run = ledger.add_run(engine.control, "generate", "g1")
    marks = [str(engine.tmp / f"marks{number}.txt") for number in range(count)]
    ledger.add_jobs(engine.control, run, [{"stage": "kill", "command": "open", "arguments": {"index": 0, "plan": {}, "run_id": "20261118-kill"}}]
                    + [{"stage": "kill", "command": "kill-unit", "item_key": str(number), "claim_mb": 500,
                        "arguments": {"args": {"database": _named(engine.config), "marks": marks[number], "pieces": pieces, "piece_seconds": piece_seconds, **faults}, "cost": 1.0}}
                       for number in range(count)])
    engine.control.execute("UPDATE engine_jobs SET state = 'done' WHERE run_id = ? AND item_key = ''", (run,))
    engine.control.execute("UPDATE engine_jobs SET state = 'ready' WHERE run_id = ? AND item_key <> ''", (run,))
    engine.control.commit()
    ids = [row["id"] for row in engine.control.execute("SELECT id FROM engine_jobs WHERE run_id = ? AND item_key <> '' ORDER BY item_key", (run,)).fetchall()]
    return run, ids, marks


def _launcher(engine, **options):
    return launcher.Launcher(engine.control, engine.forker, runner=f"{R}.unit", read_host=lambda: HOST, **options)


def _until(check, seconds, what, step=None):
    deadline = time.monotonic() + seconds
    while time.monotonic() < deadline:
        found = check()
        if found:
            return found
        if step:
            step()
        time.sleep(0.1)
    raise AssertionError(f"{what} did not happen within {seconds} seconds")


def _state(engine, row):
    engine.control.commit()
    return engine.control.execute("SELECT state, attempt FROM engine_jobs WHERE id = ?", (row,)).fetchone()


def _attempts(engine, row):
    engine.control.commit()
    return [(found["outcome"], found["error_kind"]) for found in engine.control.execute("SELECT outcome, error_kind FROM engine_attempts WHERE job_id = ? ORDER BY attempt", (row,)).fetchall()]


def _pieces(path):
    """`(starts, dones)`: how many times each piece was started and finished."""
    starts, dones = {}, {}
    with open(path, encoding="ascii") as handle:
        for line in handle:
            word, piece, _attempt = line.split()
            (starts if word == "start" else dones)[int(piece)] = (starts if word == "start" else dones).get(int(piece), 0) + 1
    return starts, dones


def _each_piece_finished_once(path, pieces):
    _starts, dones = _pieces(path)
    assert dones == {piece: 1 for piece in range(pieces)}, f"every piece is finished once, whatever was killed: {dones}"


def test_a_template_killed_while_its_job_runs_leaves_the_job_to_finish_and_the_next_job_gets_a_new_template(engine):
    _run, (first, second), marks = _units(engine, count=2, pieces=5)
    one = _launcher(engine)
    one.start({**engine.control.execute("SELECT id, run_id, stage, command, item_key, arguments, claim_mb, job_class, attempt FROM engine_jobs WHERE id = ?", (first,)).fetchone()})
    template_pid = engine.forker.hosts["plain"].template_pid
    _until(lambda: os.path.exists(marks[0]), 20, "the first piece")
    os.kill(template_pid, signal.SIGKILL)
    _until(lambda: _state(engine, first)["state"] == "done", 30, "the job whose template was killed to finish", step=one.turn)
    _until(lambda: first not in one.running, 20, "the launcher to settle the job's end", step=one.turn)
    assert _attempts(engine, first) == [("done", None)] and _state(engine, first)["attempt"] == 1
    _each_piece_finished_once(marks[0], 5)
    _until(lambda: _state(engine, second)["state"] == "done", 60, "the next job, on a new template", step=one.turn)
    assert engine.forker.hosts["plain"].template_pid != template_pid and _attempts(engine, second) == [("done", None)]
    _each_piece_finished_once(marks[1], 5)


def test_a_template_and_its_job_both_killed_have_the_job_handed_out_again_once_and_finished(engine):
    _run, (job,), (marks,) = _units(engine, pieces=6, piece_seconds=0.5)
    one = _launcher(engine)
    one.turn()
    host = engine.forker.hosts["plain"]
    _until(lambda: os.path.exists(marks), 20, "the first piece")
    process = next(iter(host.jobs.values()))
    os.kill(host.template_pid, signal.SIGKILL)
    os.kill(process.pid, signal.SIGKILL)
    _until(lambda: _state(engine, job)["state"] == "done", 60, "the job to be handed out again and finish", step=one.turn)
    assert _attempts(engine, job) == [("lost", "lost"), ("done", None)]
    _each_piece_finished_once(marks, 6)


def test_a_job_the_kernel_kills_runs_again_alone_on_a_claim_grown_by_one_and_a_half_and_finishes_once(engine, redis_server, monkeypatch):
    monkeypatch.setattr(claims, "unit_table_mb", lambda kind: 2000.0)
    monkeypatch.setattr(claims, "unit_tokens", lambda kind: ())
    _run, (job,), (marks,) = _units(engine, pieces=5, piece_seconds=0.2, kill_at=[1, 2])
    one = _launcher(engine, redis=redisqueue.connect(redis_server))
    one.turn()
    first = one.running[job].claim
    assert first.mb == 2000.0 and not first.solo
    _until(lambda: _state(engine, job)["attempt"] == 2 and job in one.running, 60, "the retry after the kill", step=one.turn)
    retry = one.running[job].claim
    assert retry.mb == pytest.approx(3000.0, rel=1e-6) and retry.solo
    _until(lambda: _state(engine, job)["state"] == "done", 60, "the retry to finish", step=one.turn)
    assert _attempts(engine, job) == [("killed_memory", "killed"), ("done", None)]
    _each_piece_finished_once(marks, 5)


def test_a_job_killed_before_it_took_its_attempt_is_handed_out_again_and_finished_once(engine):
    _run, (job,), (marks,) = _units(engine, pieces=4, die_before_attempt=True)
    one = _launcher(engine)
    one.turn()
    _until(lambda: _state(engine, job)["state"] == "done", 60, "the job to be handed out again and finish", step=one.turn)
    assert _attempts(engine, job) == [("done", None)], "the attempt the first process never took is not counted: the second process took the first attempt"
    _each_piece_finished_once(marks, 4)


def _supervisor(engine):
    spec = json.dumps({"database": _named(engine.config), "runner": f"{R}.unit"})
    env = {**os.environ, "PYTHONPATH": os.pathsep.join([SRC, *sys.path]), allowlist.EXTRA_ENV: f"{R}:*"}
    child = subprocess.Popen([sys.executable, "-m", "tests.engine_kill_supervisor", spec], cwd=SRC, env=env, start_new_session=True)
    engine.children.append(child)
    return child


def _alive(pid):
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    return True


def test_a_supervisor_killed_with_a_job_running_leaves_a_job_that_stops_at_its_piece_and_the_next_supervisor_finishes_it(engine):
    _run, (job,), (marks,) = _units(engine, pieces=8, piece_seconds=0.5)
    child = _supervisor(engine)
    _until(lambda: os.path.exists(marks), 60, "the job the supervisor handed out to begin")
    engine.control.commit()
    pid = engine.control.execute("SELECT pid FROM engine_attempts WHERE job_id = ?", (job,)).fetchone()["pid"]
    os.kill(child.pid, signal.SIGKILL)
    child.wait()
    _until(lambda: not _alive(pid), 30, "the job to see its supervisor gone and exit")
    starts, dones = _pieces(marks)
    assert set(starts) == set(dones) and len(dones) < 8, "it finished the piece it was in and no more: it did not run to the end of its eight pieces"
    assert _state(engine, job)["state"] == "running", "nothing wrote its end: the row is left running by a supervisor that is gone"
    one = _launcher(engine)
    _until(lambda: _state(engine, job)["state"] == "done", 90, "the next supervisor to hand the unit out again and finish it", step=one.turn)
    assert _attempts(engine, job) == [("lost", "lost"), ("done", None)]
    _each_piece_finished_once(marks, 8)


def test_a_supervisor_killed_after_its_job_wrote_its_end_has_nothing_to_recover(engine):
    _run, (job,), (marks,) = _units(engine, pieces=2, piece_seconds=0.1)
    child = _supervisor(engine)
    _until(lambda: _state(engine, job)["state"] == "done", 60, "the job to finish")
    os.kill(child.pid, signal.SIGKILL)
    child.wait()
    one = _launcher(engine)
    assert one.recover() == [] and one.turn() == []
    assert _state(engine, job)["state"] == "done" and _attempts(engine, job) == [("done", None)]
    _each_piece_finished_once(marks, 2)
