"""OPS.146: a unit queued again after being found stuck runs in debug mode and writes a debug file of its own, bounded and rotated, with the first attempt's stack at the top."""

import glob
import json
import os
import time

import pytest
import rq

from planetgen.queue import commands, redisqueue, runrecord, unitbeat, unitdebug, units, unitstuck
from planetgen.util import log


@pytest.fixture
def folder(tmp_path, monkeypatch):
    monkeypatch.setenv("PLANETGEN_LOG_DIR", str(tmp_path))
    return tmp_path


@pytest.fixture
def conn(redis_server):
    return redisqueue.connect(redis_server)


@pytest.fixture
def on(monkeypatch):
    monkeypatch.setattr(unitdebug, "enabled", lambda: True)


@pytest.fixture
def off(monkeypatch):
    monkeypatch.setattr(unitdebug, "enabled", lambda: False)


def test_the_setting_is_on_by_default_and_named_with_the_other_queue_keys():
    from planetgen.util import settings
    assert settings.Queue().debug_restarted_units is True
    assert "debug_restarted_units" in settings.Queue.model_fields


def test_the_file_is_named_with_the_job_the_unit_and_the_attempt(folder):
    job_id = units.job_id("20261109-120000-abcd", 2, 7)
    assert unitdebug.file_name(job_id, 2) == "unit-debug-20261109-120000-abcd-s2-u7-attempt2.log"
    assert unitdebug.file_path(job_id, 2) == str(folder / "unit-debug-20261109-120000-abcd-s2-u7-attempt2.log")


def test_begin_writes_why_and_the_first_stack_and_nothing_when_off(folder, monkeypatch):
    monkeypatch.setattr(unitdebug, "enabled", lambda: True)
    path = unitdebug.begin("r-s0-u1", 2, "unit 1 of stage 0 of run r", "unit 1 is quiet: the count stood still", "File x.py, line 3, in spin\n  while True: pass")
    text = open(path).read()
    assert os.path.dirname(path) == str(folder) and oct(os.stat(path).st_mode & 0o777) == "0o660"
    assert "unit 1 is quiet: the count stood still" in text and "in spin" in text and "attempt 2 begins" in text
    monkeypatch.setattr(unitdebug, "enabled", lambda: False)
    assert unitdebug.begin("r-s0-u2", 2, "unit 2", "why", "stack") is None
    assert len(os.listdir(folder)) == 1, "with the setting off nothing is written"


def test_a_folder_that_cannot_be_written_costs_the_unit_nothing(tmp_path, monkeypatch, on):
    blocked = tmp_path / "file"
    blocked.write_text("not a folder")
    monkeypatch.setenv("PLANETGEN_LOG_DIR", str(blocked / "logs"))
    assert unitdebug.begin("r-s0-u1", 2, "unit 1", "why", "stack") is None
    with unitdebug.session(str(blocked / "logs" / "x.log")) as path:
        assert path is None


def test_only_the_newest_files_are_kept(folder, monkeypatch, on):
    for number in range(5):
        path = folder / unitdebug.file_name(f"r-s0-u{number}", 2)
        path.write_text("x")
        (folder / (path.name + ".1")).write_text("older copy")
        os.utime(path, (1000 + number, 1000 + number))
    (folder / "error.log").write_text("not ours")
    unitdebug.prune(str(folder), keep=2)
    assert sorted(name for name in os.listdir(folder) if name.startswith("unit-debug-") and not name.endswith(".1")) == [
        unitdebug.file_name("r-s0-u3", 2), unitdebug.file_name("r-s0-u4", 2)]
    assert (folder / "error.log").exists() and not (folder / (unitdebug.file_name("r-s0-u0", 2) + ".1")).exists()


def test_a_session_records_debug_and_sql_lines_but_not_random_draws_and_leaves_the_log_as_it_found_it(folder):
    path = str(folder / "one.log")
    was = (log._logger.level, log._tracing, list(log._logger.handlers))
    with unitdebug.session(path):
        log.debug("a decision", area="decision")
        log.trace("SELECT 1 in 0.4 ms", area="sql")
        log.trace("draw 0.5", area="draw")
        log.debug("a queue line", area="queue")
        assert log._tracing
    text = open(path).read()
    assert "a decision" in text and "SELECT 1 in 0.4 ms" in text and "a queue line" in text and "draw 0.5" not in text
    assert (log._logger.level, log._tracing, list(log._logger.handlers)) == was
    log.debug("after the session", area="decision")
    assert "after the session" not in open(path).read()


def test_the_file_is_bounded_by_the_size_cap_and_the_copies_kept(folder, monkeypatch):
    monkeypatch.setattr(unitdebug, "ROTATE_BYTES", 2000)
    monkeypatch.setattr(unitdebug, "ROTATE_COPIES", 2)
    path = str(folder / "big.log")
    with unitdebug.session(path):
        for number in range(400):
            log.debug(f"line {number} " + "x" * 80, area="decision")
    files = sorted(glob.glob(path + "*"))
    assert len(files) == 3 and all(os.path.getsize(name) <= 2000 + 400 for name in files)


def test_a_beat_records_phases_waits_counter_and_memory(folder):
    clock = [1000.0]

    class Conn:
        def hset(self, *args):
            pass

        def time(self):
            return (int(clock[0]), 0)

    beat = unitbeat.UnitBeat(Conn(), "r", 0, 3, clock=lambda: clock[0], resident=lambda: 812.5, hash_key="k")
    path = str(folder / "beat.log")
    with unitdebug.session(path):
        beat.phase("weigh shells")
        beat(5, 20)
        with beat.waiting("gate"):
            clock[0] += 4.0
        clock[0] += 1.0
        beat.write()
    text = open(path).read()
    assert "unit 3: phase 'weigh shells' begins" in text and "unit 3: waits for gate" in text and "waited 4.00 s for gate" in text
    assert "counter 5 of 20" in text and "memory 812.5 MB" in text and "phase 'weigh shells'" in text


# --- the follower and the horse ---------------------------------------------------------------------------------

@pytest.fixture
def run_id(conn, monkeypatch):
    monkeypatch.setitem(commands.REGISTRY, "dbgkind", commands.CommandSpec("dbgkind", "A test", unit=(__name__, "_unit"), expand=(__name__, "_unit")))
    run = runrecord.start_run(conn, "T", [runrecord.Stage("dbgkind", "Dbg", units=2)])
    runrecord.begin_stage(conn, run, 0)
    units.enqueue_units(conn, run, 0, "dbgkind", [({"n": 0}, 1.0), ({"n": 1}, 1.0)])
    yield run
    redisqueue.queue(redisqueue.BULK_QUEUE, conn).empty()


def _unit(args, report):
    log.debug(f"the unit works on {args['n']}", area="decision")
    report.phase("working")
    report(1, 2)
    return {}


def _requeue_after_a_stop(conn, run_id, stack="FIRST STACK"):
    from planetgen.queue import beatjudge
    action = unitstuck.Stuck(conn, run_id, 0, dump=lambda raw: stack, stop=lambda job_id: None, record_error=lambda message, **fields: None, solo=False)
    now = unitbeat.server_time(conn)
    conn.hset(runrecord.progress_key(run_id, 0), 1, json.dumps({"done": 3.0, "estimate": 9.0, "beat": now, "started": now - 3000, "moved": now - 2900, "phase": "p",
                                                                "waiting_on": "", "restarts": 0, "pid": 999_999_999, "host": "elsewhere"}))
    action.step([1], {1: beatjudge.Verdict(beatjudge.QUIET, reason="the count has stood at 3 for 2900 s")})
    job = redisqueue.fetch_job(units.job_id(run_id, 0, 1), conn)
    redisqueue.queue(redisqueue.BULK_QUEUE, conn).remove(job.id)
    job.set_status(rq.job.JobStatus.STOPPED)
    action.step([], {})
    return action, job


def test_a_unit_queued_again_gets_a_debug_file_with_the_first_stack_and_the_run_log_points_at_it(conn, run_id, folder, on):
    action, job = _requeue_after_a_stop(conn, run_id)
    record = action.records()[1]
    path = record["debug_file"]
    assert path == str(folder / f"unit-debug-{job.id}-attempt2.log") and record["phase"] == "requeued"
    text = open(path).read()
    assert "FIRST STACK" in text and "unit 1 is quiet" in text
    lines = [line for line in runrecord.log_text(conn, run_id).splitlines() if path in line]
    assert len(lines) == 1 and "unit 1" in lines[0]
    assert unitstuck.debug_file_of(conn, run_id, 0, 1) == path and unitstuck.debug_file_of(conn, run_id, 0, 0) == ""


def test_the_restarted_unit_writes_into_that_file_in_its_horse_and_the_other_unit_does_not(conn, run_id, folder, on):
    _requeue_after_a_stop(conn, run_id)
    path = unitstuck.debug_file_of(conn, run_id, 0, 1)
    units._do_unit(conn, run_id, 0, 1, "dbgkind", {"n": 1}, 1.0)
    text = open(path).read()
    assert "the unit works on 1" in text and "unit 1: phase 'working' begins" in text and "FIRST STACK" in text
    units._do_unit(conn, run_id, 0, 0, "dbgkind", {"n": 0}, 1.0)
    assert "the unit works on 0" not in open(path).read()
    assert len(glob.glob(str(folder / "unit-debug-*"))) == 1


def test_with_the_setting_off_nothing_is_written_and_the_run_log_has_no_line(conn, run_id, folder, off):
    action, job = _requeue_after_a_stop(conn, run_id)
    assert "debug_file" not in action.records()[1] and os.listdir(folder) == []
    assert "debug mode" not in runrecord.log_text(conn, run_id)
    units._do_unit(conn, run_id, 0, 1, "dbgkind", {"n": 1}, 1.0)
    assert os.listdir(folder) == []
