# planetgen/queue/pagerun.py

"""
A Generate page job as a run of the ledger (OPS.244, slice B; `docs/design/operations/job-engine.md`).

A job the admin Generate page starts (`web/jobs.py`) is a ledger run of kind `page` with one row for the job and one row for each of its steps, in the control database:

* **The job's row** (stage `page`, command `page`, item key the job's id, state `done` from the start) holds what the job is: `{"id", "kind", "title", "admin", "database", "draft", "created_at", "cwd", "root"}` (`root` the jobs directory whose folder holds the job's logs), and
  once the first step has begun, the estimates of its steps and stages (`step_estimates`, `stage_estimates`).
* **A step's row** (command `page-step`, item key its number, three digits) holds `{"label", "argv", "stages"}`. Only the next step is `ready`; the step that ends well makes the one after it `ready`,
  so the steps run in order and the first that fails stops the job. A step is one job process forked by the launcher (`launcher.Launcher`, `forkedjob.run`) that runs the step's command line as a process
  of its own and follows it (`web.job_runner.run_step`).

The page asks the rows what the job is doing (`read`: its status, the step it is at, when it started and ended, why it failed) and stops it with a request to the supervisor (`controls.stop_engine_run`).
What stays in the job's folder under the jobs directory is what a process writes as it goes: `output.log`, `records.log`, `progress.json`, `memory.json`.

A step is never run again by the launcher when its process dies (`launcher.Launcher.ended`): a reset or a run that was cut off is the admin's to start again (the page's Retry).
"""

import json
import time

from planetgen.queue import ledger
from planetgen.util import log, process

KIND = "page"
"""str: The kind of the ledger run of a Generate page job."""

STAGE = "page"
"""str: The stage of every row of the run."""

JOB = "page"
"""str: The command of the row that holds what the job is."""

STEP = "page-step"
"""str: The command of a step's row."""

CLASS = "page"
"""str: The job class of the rows of the run (the launcher hands them out on the lane for classes it does not name)."""

STEP_MB = 700
"""int: What a step claims: the runner that follows the step's process is small, and the process it starts is counted under it (`memorygate.held_mb` takes every process under a job)."""

STARTING_LIMIT_SECONDS = 300
"""int: How long a job none of whose steps has begun is `starting` (a supervisor may have to start, and its template to be built, for the first step); older than that it is `interrupted`, unless the
launcher has said why its step waits (the row's `waited`: the dispatch rule held it, so a supervisor is looking at it)."""

FINISHED = ("succeeded", "failed", "canceled", "interrupted")
"""tuple: The statuses a job that is over has."""

ENDED_BY_PROCESS = ("lost", "killed_memory", "timed_out", "stuck")
"""tuple: The outcomes of an attempt whose process was gone when it was found, which the page calls `interrupted`."""


def connect():
    """A connection to the control database.

    Raises:
        OSError: There is none (the ledger is where a job's state is kept, so no job runs without it).
    """
    from planetgen.db import store
    try:
        return store.get_control_connection(store.control_mysql_config())
    except Exception as exc:  # noqa: BLE001 -- the callers handle the OSError of a missing server
        raise OSError(f"no control database answers (config.json's mysql section); it holds the state of the jobs: {exc}") from exc


def _step_key(number):
    return f"{number:03d}"


def start(control, job, steps, now=None):
    """
    Writes the run of `job` (`{"id", "kind", "title", "admin", "database", "draft", "created_at", "cwd"}`) with a row for each of its `steps` (`{"label", "argv", "stages"}`), the first one `ready`, and
    commits. Returns the ledger run's id.
    """
    run = ledger.add_run(control, KIND, job.get("database"))
    named = {name: job.get(name) for name in ("id", "kind", "title", "admin", "database", "draft", "created_at", "cwd", "root")}
    rows = [{"stage": STAGE, "command": JOB, "item_key": job["id"], "arguments": named, "job_class": CLASS}]
    for number, step in enumerate(steps, start=1):
        rows.append({"stage": STAGE, "command": STEP, "item_key": _step_key(number), "job_class": CLASS, "claim_mb": STEP_MB,
                     "arguments": {"label": step["label"], "argv": list(step["argv"]), "stages": step.get("stages")}})
    ledger.add_jobs(control, run, rows)
    ids = [row["id"] for row in control.execute("SELECT id FROM engine_jobs WHERE run_id = ? ORDER BY id", (run,)).fetchall()]
    ledger.set_job_state(control, ids[0], "done")
    if len(ids) > 1:
        ledger.set_job_state(control, ids[1], "ready")
    else:
        ledger.set_run_state(control, run, "done")                  # a job of no steps has nothing to run
    control.commit()
    return run


def _f(value):
    return None if value is None else float(value)


def _loaded(text):
    try:
        loaded = json.loads(text) if text else {}
    except ValueError:
        return {}
    return loaded if isinstance(loaded, dict) else {}


def _found(control, job_id):
    return control.execute(
        "SELECT j.id AS job_row, j.run_id AS run, j.arguments AS arguments, r.state AS state, r.held_reason AS held, UNIX_TIMESTAMP(r.started_at) AS started, UNIX_TIMESTAMP(r.finished_at) AS finished"
        " FROM engine_jobs j JOIN engine_runs r ON r.id = j.run_id WHERE r.kind = ? AND j.stage = ? AND j.command = ? AND j.item_key = ?", (KIND, STAGE, JOB, job_id)).fetchone()


def job_of(control, run):
    """The job of the ledger run `run` as the runner needs it: the fields of the job's row and its `steps` (`{"label", "argv", "stages"}`), or `None` for a run that is not a page job's."""
    row = control.execute("SELECT arguments FROM engine_jobs WHERE run_id = ? AND stage = ? AND command = ?", (run, STAGE, JOB)).fetchone()
    if row is None:
        return None
    steps = [_loaded(found["arguments"]) for found in control.execute("SELECT arguments FROM engine_jobs WHERE run_id = ? AND command = ? ORDER BY item_key", (run, STEP)).fetchall()]
    return {**_loaded(row["arguments"]), "steps": steps}


def run_of(control, job_id):
    """The ledger run of the job `job_id`, or `None` for a job the ledger does not know."""
    found = _found(control, job_id)
    return None if found is None else int(found["run"])


def _step_rows(control, run):
    return control.execute(
        "SELECT j.id, j.item_key, j.state, j.waited, j.arguments, j.attempt, a.pid AS pid, a.outcome AS outcome, a.error_kind AS error_kind, a.error_text AS error_text,"
        " UNIX_TIMESTAMP(a.started_at) AS began, UNIX_TIMESTAMP(a.ended_at) AS ended"
        " FROM engine_jobs j LEFT JOIN engine_attempts a ON a.job_id = j.id AND a.attempt = j.attempt WHERE j.run_id = ? AND j.command = ? ORDER BY j.item_key", (run, STEP)).fetchall()


def read(control, job_id, now=None):
    """
    What the ledger says of the job `job_id`, or `None` for an unknown job: the job's own fields (`id`, `kind`, `title`, `admin`, `database`, `draft`, `created_at`), the `steps` (`[{"label", "argv", "stages"}]`)
    and the state: `status` (`starting`, `running`, `succeeded`, `failed`, `canceled`, `interrupted`), `step` (1-based, 0 before the first), `step_started_at`, `started_at`, `finished_at`, `error`,
    `step_estimates`, `stage_estimates`, `run` (the ledger run) and `runner_pid` (the process of the step that runs, else `None`).

    A job whose step's process is gone and not yet settled by the launcher is `interrupted`, as is a job nothing started in `STARTING_LIMIT_SECONDS`.
    """
    found = _found(control, job_id)
    if found is None:
        return None
    now = time.time() if now is None else now
    job = _loaded(found["arguments"])
    rows = _step_rows(control, found["run"])
    steps = [{"label": _loaded(row["arguments"]).get("label"), "argv": _loaded(row["arguments"]).get("argv"), "stages": _loaded(row["arguments"]).get("stages")} for row in rows]
    began = [number for number, row in enumerate(rows, start=1) if row["state"] in ("running", "done", "failed") or (row["state"] == "canceled" and row["attempt"])]
    step = max(began) if began else 0
    current = rows[step - 1] if step else None
    failed = next((row for row in rows if row["state"] == "failed"), None)
    error = None
    if found["state"] == "canceled":
        status = "canceled"
    elif failed is not None:
        status = "interrupted" if failed["outcome"] in ENDED_BY_PROCESS and failed["error_kind"] in ("lost", "killed") else "failed"
        error = failed["error_text"] or found["held"]
    elif found["state"] == "done":
        status = "succeeded"
    elif current is not None and current["state"] == "running":
        status = "running" if process.pid_alive(current["pid"]) else "interrupted"
    elif now - float(job.get("created_at") or 0) < STARTING_LIMIT_SECONDS or any(row["state"] == "ready" and row["waited"] for row in rows):
        status = "starting" if current is None else "running"
    else:
        status = "interrupted"
    ended = None
    if status in FINISHED:
        ended = _f(found["finished"]) or (_f(current["ended"]) if current is not None and current["ended"] else None) or now
    return {
        **job, "steps": steps, "status": status, "step": step, "run": int(found["run"]), "error": error,
        "runner_pid": current["pid"] if current is not None and current["state"] == "running" else None,
        "started_at": _f(found["started"]) or (_f(rows[0]["began"]) if rows and rows[0]["began"] else None),
        "finished_at": ended, "step_started_at": _f(current["began"]) if current is not None else None,
        "step_estimates": job.get("step_estimates") or [], "stage_estimates": job.get("stage_estimates") or [],
    }


def set_estimates(control, run, step_estimates, stage_estimates):
    """Keeps what the job's steps and stages took in earlier jobs on the job's row (`read` hands them to the page). Commits."""
    row = control.execute("SELECT id, arguments FROM engine_jobs WHERE run_id = ? AND stage = ? AND command = ?", (run, STAGE, JOB)).fetchone()
    if row is None:
        return
    named = _loaded(row["arguments"])
    named.update(step_estimates=step_estimates, stage_estimates=stage_estimates)
    control.execute("UPDATE engine_jobs SET arguments = ?, updated_at = NOW(6) WHERE id = ?", (json.dumps(named), row["id"]))
    control.commit()


def begin(control, run):
    """The first step has begun: the run is `running` (its start time is stamped once). Commits."""
    ledger.set_run_state(control, run, "running")
    control.commit()


def step_done(control, row_id, run, peak_mb=None):
    """
    A step ended well: its attempt and its row are `done`, and the step after it is `ready`, or the run is `done` when it was the last. Only a step still `running` is changed (a Cancel that came first
    stays). Commits. Returns whether the job is over (the step was the last, or the job was ended before).
    """
    if control.execute("SELECT state FROM engine_jobs WHERE id = ?", (row_id,)).fetchone()["state"] != "running":
        control.rollback()          # a Cancel came first: the job is over
        return True
    ledger.end_running(control, row_id, "done", peak_mb=peak_mb)
    ledger.set_job_state(control, row_id, "done")
    waiting = control.execute("SELECT id FROM engine_jobs WHERE run_id = ? AND command = ? AND state = 'waiting' ORDER BY item_key LIMIT 1", (run, STEP)).fetchone()
    if waiting is not None:
        ledger.set_job_state(control, waiting["id"], "ready")
    else:
        ledger.set_run_state(control, run, "done")
    control.commit()
    return waiting is None


def step_failed(control, row_id, run, text, peak_mb=None, kind="failed"):
    """A step failed (`text` says how): its attempt and row are `failed`, the steps after it are `canceled` and the run is `held` with the reason. Only a step still `running` is changed. Commits."""
    if control.execute("SELECT state FROM engine_jobs WHERE id = ?", (row_id,)).fetchone()["state"] == "running":
        ledger.end_running(control, row_id, "failed", peak_mb=peak_mb, error_kind=kind, error_text=text)
        ledger.set_job_state(control, row_id, "failed", waited=text[:ledger.WAITED_CHARS])
    fail_rest(control, run, text)


def fail_rest(control, run, text):
    """The steps of the run that have not begun are `canceled` and the run is `held` for `text`, unless it was canceled. Commits."""
    control.execute("UPDATE engine_jobs SET state = 'canceled', waited = NULL, updated_at = NOW(6) WHERE run_id = ? AND command = ? AND state IN ('waiting', 'ready')", (run, STEP))
    if control.execute("SELECT state FROM engine_runs WHERE id = ?", (run,)).fetchone()["state"] not in ("canceled", "done"):
        ledger.set_run_state(control, run, "held", text[:ledger.WAITED_CHARS])
    control.commit()


def delete(control, job_id):
    """Deletes the run of the job `job_id` with its rows (the job's folder is pruned with it). Commits."""
    run = run_of(control, job_id)
    if run is not None:
        control.execute("DELETE FROM engine_runs WHERE id = ?", (run,))
    control.commit()


def live_ids(control, job_ids):
    """Which of the page jobs `job_ids` have a run that is not over (waiting or running)."""
    live = set()
    for job_id in job_ids:
        found = _found(control, job_id)
        if found is not None and found["state"] in ("waiting", "running"):
            live.add(job_id)
    return live


def tidy(control, connection):
    """Forgets the jobs of the supervisor's keep-alive list (`supervisor.pages_key`) whose run is over or gone, for a job that ended without anyone to say so. Best effort."""
    from planetgen.queue import runrecord, supervisor
    try:
        kept = [runrecord._s(name) for name in connection.hkeys(supervisor.pages_key())]
        for name in set(kept) - live_ids(control, kept):
            connection.hdel(supervisor.pages_key(), name)
    except Exception as exc:  # noqa: BLE001 -- a job left on the list only keeps the supervisor awake a little longer
        log.swallowed("queue/pagerun", "tidy", exc)
