# tests/test_pagerun.py

"""
A Generate page job as a run of the ledger (OPS.244, slice B; `planetgen/queue/pagerun.py`): the rows `start` writes, what `read` says of them
at each point of a job's life, how a step that ends moves the job on, and how the rows are deleted.
"""

import os
import time

import pytest

from planetgen.queue import ledger, pagerun
from tests import page_jobs


@pytest.fixture
def control():
    page_jobs.forget_all_jobs()
    with page_jobs.ledger_rows() as connection:
        yield connection
    page_jobs.forget_all_jobs()


def _job(job_id="jobone", steps=2, created=None):
    return {"id": job_id, "kind": "reset", "title": "A job", "admin": "boss", "database": "planetgen_x", "draft": None,
            "created_at": time.time() if created is None else created, "cwd": "/tmp", "root": "/tmp/jobs"}, \
        [{"label": f"Step {number}", "argv": ["true", str(number)], "stages": None} for number in range(1, steps + 1)]


def _start(control, job_id="jobone", steps=2, created=None):
    job, rows = _job(job_id, steps, created)
    return pagerun.start(control, job, rows)


def _states(control, run):
    return [row["state"] for row in control.execute("SELECT state FROM engine_jobs WHERE run_id = ? AND command = ? ORDER BY item_key", (run, pagerun.STEP)).fetchall()]


def _step_ids(control, run):
    return [row["id"] for row in control.execute("SELECT id FROM engine_jobs WHERE run_id = ? AND command = ? ORDER BY item_key", (run, pagerun.STEP)).fetchall()]


def test_a_job_is_a_row_for_the_job_and_a_row_for_each_step_with_the_first_ready(control):
    run = _start(control, steps=3)
    assert _states(control, run) == ["ready", "waiting", "waiting"]
    kept = pagerun.job_of(control, run)
    assert kept["id"] == "jobone" and kept["title"] == "A job" and kept["root"] == "/tmp/jobs"
    assert [step["label"] for step in kept["steps"]] == ["Step 1", "Step 2", "Step 3"]
    assert kept["steps"][1]["argv"] == ["true", "2"]
    assert pagerun.run_of(control, "jobone") == run
    assert pagerun.run_of(control, "nosuchjob") is None
    assert pagerun.job_of(control, 987654321) is None
    assert control.execute("SELECT kind FROM engine_runs WHERE id = ?", (run,)).fetchone()["kind"] == pagerun.KIND


def test_a_job_that_nothing_has_begun_is_starting(control):
    _start(control)
    job = pagerun.read(control, "jobone")
    assert job["status"] == "starting" and job["step"] == 0 and job["runner_pid"] is None and job["started_at"] is None
    assert pagerun.read(control, "nosuchjob") is None


def test_a_step_that_runs_makes_the_job_running_and_names_its_process(control):
    run = _start(control)
    page_jobs.running("jobone")
    job = pagerun.read(control, "jobone")
    assert job["status"] == "running" and job["step"] == 1 and job["runner_pid"] == os.getpid()
    assert job["started_at"] and job["finished_at"] is None and job["run"] == run


def test_a_step_that_ends_well_makes_the_next_one_ready_and_the_last_ends_the_job(control):
    run = _start(control)
    first, second = _step_ids(control, run)
    ledger.begin_attempt(control, first, pid=os.getpid())
    assert pagerun.step_done(control, first, run) is False
    assert _states(control, run) == ["done", "ready"]
    ledger.begin_attempt(control, second, pid=os.getpid())
    assert pagerun.step_done(control, second, run) is True
    job = pagerun.read(control, "jobone")
    assert job["status"] == "succeeded" and job["step"] == 2 and job["finished_at"] and job["error"] is None
    assert control.execute("SELECT state FROM engine_runs WHERE id = ?", (run,)).fetchone()["state"] == "done"


def test_a_step_that_fails_stops_the_job_and_cancels_the_steps_after_it(control):
    run = _start(control, steps=3)
    first, second, _third = _step_ids(control, run)
    ledger.begin_attempt(control, first, pid=os.getpid())
    pagerun.step_done(control, first, run)
    ledger.begin_attempt(control, second, pid=os.getpid())
    pagerun.step_failed(control, second, run, "Step 2 exited with status 3.")
    assert _states(control, run) == ["done", "failed", "canceled"]
    job = pagerun.read(control, "jobone")
    assert job["status"] == "failed" and job["step"] == 2 and job["error"] == "Step 2 exited with status 3."
    assert control.execute("SELECT state, held_reason FROM engine_runs WHERE id = ?", (run,)).fetchone()["state"] == "held"


def test_a_step_whose_process_is_gone_is_an_interrupted_job(control):
    _start(control)
    page_jobs.running("jobone", pid=page_jobs.DEAD_PID)
    assert pagerun.read(control, "jobone")["status"] == "interrupted"


def test_a_step_the_launcher_found_lost_is_an_interrupted_job_and_one_that_failed_is_failed(control):
    run = _start(control)
    first = _step_ids(control, run)[0]
    ledger.begin_attempt(control, first, pid=page_jobs.DEAD_PID)
    ledger.end_running(control, first, "lost", error_kind="lost", error_text="the process was gone")
    ledger.set_job_state(control, first, "failed")
    assert pagerun.read(control, "jobone")["status"] == "interrupted"


def test_a_cancel_that_came_first_is_never_undone_by_the_step_that_then_ends(control):
    from planetgen.queue import controls

    run = _start(control)
    first = _step_ids(control, run)[0]
    ledger.begin_attempt(control, first, pid=os.getpid())
    controls.end_rows(control, run)
    control.commit()
    assert pagerun.step_done(control, first, run) is True
    pagerun.step_failed(control, first, run, "late")
    job = pagerun.read(control, "jobone")
    assert job["status"] == "canceled"
    assert control.execute("SELECT state FROM engine_runs WHERE id = ?", (run,)).fetchone()["state"] == "canceled"


def test_a_job_nothing_started_for_a_long_time_is_interrupted_unless_the_launcher_explained_the_wait(control):
    run = _start(control, created=time.time() - pagerun.STARTING_LIMIT_SECONDS - 5)
    assert pagerun.read(control, "jobone")["status"] == "interrupted"
    ready = _step_ids(control, run)[0]
    ledger.set_job_state(control, ready, "ready", waited="the memory gate holds it")
    control.commit()
    assert pagerun.read(control, "jobone")["status"] == "starting"


def test_the_estimates_are_kept_on_the_jobs_row_and_read_back(control):
    run = _start(control)
    pagerun.set_estimates(control, run, [12.5, 30.0], [[1.0], [2.0]])
    job = pagerun.read(control, "jobone")
    assert job["step_estimates"] == [12.5, 30.0] and job["stage_estimates"] == [[1.0], [2.0]]
    assert pagerun.read(control, "jobone")["title"] == "A job"


def test_only_a_job_whose_run_is_not_over_is_live(control):
    run = _start(control, "jobone")
    _start(control, "jobtwo")
    pagerun.fail_rest(control, pagerun.run_of(control, "jobtwo"), "stopped")
    ledger.set_run_state(control, pagerun.run_of(control, "jobtwo"), "canceled")
    control.commit()
    assert pagerun.live_ids(control, ["jobone", "jobtwo", "nosuchjob"]) == {"jobone"}
    assert run


def test_deleting_a_job_deletes_its_run_and_its_rows(control):
    run = _start(control)
    pagerun.delete(control, "jobone")
    assert pagerun.run_of(control, "jobone") is None
    assert control.execute("SELECT COUNT(*) AS n FROM engine_jobs WHERE run_id = ?", (run,)).fetchone()["n"] == 0
    pagerun.delete(control, "jobone")      # a job that is gone is not an error


def test_no_control_database_is_an_oserror(monkeypatch):
    from planetgen.db import store

    def refuse(*args, **kwargs):
        raise RuntimeError("down")

    monkeypatch.setattr(store, "get_control_connection", refuse)
    with pytest.raises(OSError, match="no control database"):
        pagerun.connect()


def test_a_job_of_no_steps_is_over_at_once_and_succeeded(control):
    run = _start(control, steps=0)
    assert control.execute("SELECT state FROM engine_runs WHERE id = ?", (run,)).fetchone()["state"] == "done"
    job = pagerun.read(control, "jobone")
    assert job["status"] == "succeeded" and job["step"] == 0 and job["steps"] == []


def _queue(control, ticket, steps=2, after=None, root="/tmp/jobs", now=None):
    job, rows = _job("unused", steps)
    job["root"] = root
    job["title"] = f"Job {ticket}"
    return pagerun.queue(control, job, rows, ticket, after_job=after, now=now)


def test_a_queued_job_is_a_waiting_run_whose_rows_all_wait_and_is_no_job_yet(control):
    run = _queue(control, "ticketone", steps=2, after="jobbefore", now=1000.0)
    assert _states(control, run) == ["waiting", "waiting"]
    assert control.execute("SELECT state FROM engine_runs WHERE id = ?", (run,)).fetchone()["state"] == "waiting"
    (entry,) = pagerun.waiting(control)
    assert (entry["ticket"], entry["title"], entry["after_job"], entry["queued_at"], entry["held"], entry["run"], entry["position"]) == ("ticketone", "Job ticketone", "jobbefore", 1000.0, None, run, 1)
    assert pagerun.live_ids(control, ["ticketone"]) == {"ticketone"}
    assert pagerun.job_of(control, run)["steps"][1]["argv"] == ["true", "2"]


def test_the_waiting_jobs_come_in_the_order_they_queued_and_can_be_told_by_folder(control):
    _queue(control, "first", root="/tmp/a")
    _queue(control, "second", root="/tmp/b")
    _queue(control, "third", root="/tmp/a")
    assert [(e["ticket"], e["position"]) for e in pagerun.waiting(control)] == [("first", 1), ("second", 2), ("third", 3)]
    assert [(e["ticket"], e["position"]) for e in pagerun.waiting(control, root="/tmp/a")] == [("first", 1), ("third", 2)]


def test_a_waiting_job_can_be_held_and_told_to_wait_for_another_job(control):
    _queue(control, "ticketone")
    (entry,) = pagerun.waiting(control)
    pagerun.hold_queued(control, entry["row"], "The job before it did not succeed (failed).")
    pagerun.set_after(control, entry["row"], "jobafter")
    (held,) = pagerun.waiting(control)
    assert held["held"] == "The job before it did not succeed (failed)." and held["after_job"] == "jobafter"


def test_only_one_process_claims_a_waiting_job_and_a_claimed_job_leaves_the_line(control):
    _queue(control, "ticketone")
    (entry,) = pagerun.waiting(control)
    assert pagerun.claim_queued(control, entry["row"]) is True
    assert pagerun.claim_queued(control, entry["row"]) is False
    assert pagerun.waiting(control) == []


def test_a_claimed_job_that_could_not_start_goes_back_to_its_place_with_the_reason(control):
    _queue(control, "first")
    _queue(control, "second")
    first, _second = pagerun.waiting(control)
    pagerun.claim_queued(control, first["row"])
    pagerun.restore_queued(control, first["row"], "It could not be started: no room")
    assert [(e["ticket"], e["held"]) for e in pagerun.waiting(control)] == [("first", "It could not be started: no room"), ("second", None)]
    pagerun.restore_queued(control, first["row"])
    assert pagerun.waiting(control)[0]["held"] is None


def test_a_released_job_takes_its_real_id_and_its_first_step_is_ready(control):
    run = _queue(control, "ticketone", steps=2)
    (entry,) = pagerun.waiting(control)
    pagerun.claim_queued(control, entry["row"])
    job, _rows = _job("jobreal", 2)
    pagerun.release_queued(control, entry["run"], entry["row"], job)
    assert pagerun.run_of(control, "jobreal") == run and pagerun.run_of(control, "ticketone") is None
    assert _states(control, run) == ["ready", "waiting"]
    assert pagerun.read(control, "jobreal")["status"] == "starting"
    assert pagerun.waiting(control) == []


def test_a_waiting_job_that_is_discarded_is_gone_with_its_rows(control):
    run = _queue(control, "ticketone")
    control.execute("DELETE FROM engine_runs WHERE id = ?", (run,))
    control.commit()
    assert pagerun.waiting(control) == [] and pagerun.live_ids(control, ["ticketone"]) == set()
