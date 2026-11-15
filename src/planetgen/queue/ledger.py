# planetgen/queue/ledger.py

"""
The job engine's ledger (OPS.215, L1 of `docs/design/operations/job-engine.md`): the rows that say what ran.

Four tables. Three are in the control database (`control_schema.sql`, v23): `engine_runs` (a run), `engine_jobs` (the jobs of a run, with the number of
the newest attempt) and `engine_attempts` (one row per attempt, with how it ended). The fourth, `engine_items`, is in each galaxy's own database (migration
0261): one row per job that writes to that galaxy, so a job's commits and its done mark are one transaction.

**Nothing reads or writes these rows yet.** The supervisor (OPS.207) is the writer and the plan-as-rows item (OPS.221) the first reader; this module is what
both build on, so the column names live in one place.

* **A plan.** `add_run` and `add_jobs` insert a run and its jobs (the cut of a stage, in one transaction with its done mark: the caller commits). `add_jobs`
  is safe to repeat: a job of the same stage, command and item key is the same job, and is not added twice.
* **An attempt.** `begin_attempt` is the supervisor's dispatch write: it raises the job's attempt number, opens an attempt row and, given the galaxy
  connection, writes the number to `engine_items` before the job starts. It is the only write the supervisor waits for.
* **The fence.** `commit_item` is what a job calls in the same transaction as its data: it updates its `engine_items` row `WHERE attempt = mine`. No row
  touched means another attempt took the job: the caller rolls back and stops (`Replaced`). With `done=True` it sets the done time, so the data and the mark
  commit together.
* **The writer.** `Writer` runs ledger calls on its own thread and connection, so the supervisor's loop never waits on MySQL.
"""

import json
import queue
import threading

from planetgen.util import log

RUN_STATES = ("waiting", "running", "held", "done", "canceled")
"""tuple: The states of a run (the column's CHECK)."""

JOB_STATES = ("waiting", "ready", "running", "done", "failed", "canceled")
"""tuple: The states of a job (the column's CHECK)."""

ATTEMPT_OUTCOMES = ("running", "done", "failed", "killed_memory", "timed_out", "stuck", "replaced", "lost")
"""tuple: How an attempt ended, `running` while it has not (the column's CHECK)."""

JOB_FIELDS = ("stage", "command", "item_key", "arguments", "needs", "makes", "job_class", "claim_mb", "tokens")
"""tuple: What a job of a plan says; `stage` and `command` are needed, the rest have defaults (`item_key` '', `job_class` 'bulk', the others empty)."""


class Replaced(Exception):
    """The job's row in `engine_items` belongs to another attempt: this one was replaced and must roll back and stop."""


def _text(value):
    """`value` as the JSON text a column holds, or `None` for nothing; text is stored as it is."""
    if value is None or isinstance(value, str):
        return value
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def add_run(conn, kind, galaxy=None):
    """Inserts a run in state `waiting` and returns its id. The caller commits."""
    cursor = conn.execute("INSERT INTO engine_runs (kind, galaxy, state, created_at) VALUES (?, ?, 'waiting', NOW(6))", (kind, galaxy))
    return cursor.lastrowid


def set_run_state(conn, run_id, state, reason=None):
    """Moves a run to `state` (`reason` is why it is held, else cleared) and stamps the time it started or finished. The caller commits."""
    if state not in RUN_STATES:
        raise ValueError(f"A run is never {state!r}: {RUN_STATES}.")
    conn.execute("UPDATE engine_runs SET state = ?, held_reason = ?,"
                 " started_at = IF(? = 'running' AND started_at IS NULL, NOW(6), started_at),"
                 " finished_at = IF(? IN ('done', 'canceled'), NOW(6), NULL) WHERE id = ?", (state, reason if state == "held" else None, state, state, run_id))


def add_jobs(conn, run_id, jobs):
    """
    Inserts the jobs of a plan (dicts with `JOB_FIELDS`; `arguments`, `needs` and `makes` may be any JSON value) in one statement. A job already in the run
    (same stage, command and item key) is left as it is. Returns the number added. The caller commits.
    """
    rows = []
    for job in jobs:
        unknown = set(job) - set(JOB_FIELDS)
        if unknown or "stage" not in job or "command" not in job:
            raise ValueError(f"A job of a plan has a stage and a command and only {JOB_FIELDS}: {sorted(unknown) or sorted(job)}.")
        rows.append((run_id, job["stage"], job["command"], job.get("item_key", ""), _text(job.get("arguments")), _text(job.get("needs")), _text(job.get("makes")),
                     job.get("job_class", "bulk"), job.get("claim_mb"), job.get("tokens")))
    if not rows:
        return 0
    marks = "(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, 'waiting', NOW(6), NOW(6))"
    cursor = conn.execute("INSERT IGNORE INTO engine_jobs (run_id, stage, command, item_key, arguments, needs, makes, job_class, claim_mb, tokens, state, created_at,"
                          f" updated_at) VALUES {', '.join([marks] * len(rows))}", tuple(value for row in rows for value in row))
    return cursor.rowcount


def set_job_state(conn, job_id, state, waited=None, not_before=None):
    """Moves a job to `state`, with the last reason it was not started (`waited`) and the time a retry may start (`not_before`, a datetime or None). The caller commits."""
    if state not in JOB_STATES:
        raise ValueError(f"A job is never {state!r}: {JOB_STATES}.")
    conn.execute("UPDATE engine_jobs SET state = ?, waited = ?, not_before = ?, updated_at = NOW(6) WHERE id = ?", (state, waited, not_before, job_id))


WAITED_CHARS = 255
"""int: How much of a reason a job's `waited` column keeps (`VARCHAR(255)`)."""


def set_wait(conn, job_id, reason):
    """
    Notes why a job that has not started is not starting (`reason`, the words the memory gate gave a worker that held it, OPS.224), on the row of a job that is still `waiting`; a job that has begun
    (or is over) keeps its state. Returns whether the row took it. The caller commits.
    """
    reason = None if reason is None else (reason if len(reason) <= WAITED_CHARS else reason[:WAITED_CHARS - 3] + "...")
    return conn.execute("UPDATE engine_jobs SET waited = ?, updated_at = NOW(6) WHERE id = ? AND state = 'waiting'", (reason, job_id)).rowcount > 0


def begin_attempt(conn, job_id, galaxy_conn=None, pid=None):
    """
    The dispatch write: raises the job's attempt number, opens its attempt row (`running`) and sets the job running, and returns the new number. Given the galaxy's
    connection, it also writes the number to the job's `engine_items` row, before the job starts, so the fence has a number to hold. Commits both.
    """
    end_running(conn, job_id, "lost", error_kind="lost", error_text="the attempt wrote no end: its process was gone when the next attempt began")       # OPS.223: an attempt nobody closed
    conn.execute("UPDATE engine_jobs SET attempt = attempt + 1, state = 'running', waited = NULL, updated_at = NOW(6) WHERE id = ?", (job_id,))
    attempt = conn.execute("SELECT attempt FROM engine_jobs WHERE id = ?", (job_id,)).fetchone()["attempt"]
    conn.execute("INSERT INTO engine_attempts (job_id, attempt, pid, started_at, outcome) VALUES (?, ?, ?, NOW(6), 'running')", (job_id, attempt, pid))
    conn.commit()
    if galaxy_conn is not None:
        galaxy_conn.execute("INSERT INTO engine_items (job_id, attempt, updated_at) VALUES (?, ?, NOW(6))"
                            " ON DUPLICATE KEY UPDATE attempt = VALUES(attempt), cursor_key = NULL, done_at = NULL, updated_at = NOW(6)", (job_id, attempt))
        galaxy_conn.commit()
    return attempt


def end_attempt(conn, job_id, attempt, outcome, peak_mb=None, cpu_seconds=None, error_kind=None, error_text=None, stack_file=None):
    """Closes the attempt row with how it ended and what it cost; the caller sets the job's state and commits."""
    if outcome not in ATTEMPT_OUTCOMES or outcome == "running":
        raise ValueError(f"An attempt ends as one of {ATTEMPT_OUTCOMES[1:]}, not {outcome!r}.")
    conn.execute("UPDATE engine_attempts SET ended_at = NOW(6), outcome = ?, peak_mb = ?, cpu_seconds = ?, error_kind = ?, error_text = ?, stack_file = ?"
                 " WHERE job_id = ? AND attempt = ?", (outcome, peak_mb, cpu_seconds, error_kind, _clip(error_text), stack_file, job_id, attempt))


def end_running(conn, job_id, outcome, peak_mb=None, cpu_seconds=None, error_kind=None, error_text=None, stack_file=None):
    """
    Closes the attempts of the job that are still `running` with `outcome` (OPS.223): the end of an attempt whose process could not write it itself (killed for memory, stopped as stuck,
    gone), written by whoever saw it go, or by the next `begin_attempt` as `lost`. Returns the numbers of the attempts closed, oldest first. The caller commits.
    """
    if outcome not in ATTEMPT_OUTCOMES or outcome == "running":
        raise ValueError(f"An attempt ends as one of {ATTEMPT_OUTCOMES[1:]}, not {outcome!r}.")
    numbers = [row["attempt"] for row in conn.execute("SELECT attempt FROM engine_attempts WHERE job_id = ? AND outcome = 'running' ORDER BY attempt", (job_id,)).fetchall()]
    if numbers:
        conn.execute("UPDATE engine_attempts SET ended_at = NOW(6), outcome = ?, peak_mb = ?, cpu_seconds = ?, error_kind = ?, error_text = ?, stack_file = ?"
                     " WHERE job_id = ? AND outcome = 'running'", (outcome, peak_mb, cpu_seconds, error_kind, _clip(error_text), stack_file, job_id))
    return numbers


def attempt_row(conn, job_id, number):
    """The attempt row (`outcome`, `ended_at`, `error_kind`, `error_text`, `stack_file`, `peak_mb`, `cpu_seconds`, `started_at`, `pid`) of the job's attempt `number`, or `None`."""
    return conn.execute("SELECT outcome, started_at, ended_at, pid, peak_mb, cpu_seconds, error_kind, error_text, stack_file FROM engine_attempts WHERE job_id = ? AND attempt = ?",
                        (job_id, number)).fetchone()


ERROR_TEXT_CHARS = 8000
"""int: How much of an error's text an attempt row keeps (the column is TEXT, 64 KB; the whole stack is in the file the row names)."""


def _clip(text):
    return text if text is None or len(text) <= ERROR_TEXT_CHARS else text[:ERROR_TEXT_CHARS] + f"... ({len(text)} characters)"


def commit_item(galaxy_conn, job_id, attempt, cursor_key=None, done=False):
    """
    The fence. Called in the same transaction as the job's data: updates the job's `engine_items` row `WHERE attempt = mine`, with the place reached
    (`cursor_key`) and, for the last commit, the done time. Raises `Replaced` when the row is not this attempt's; the caller rolls back and stops. Does not commit.
    """
    touched = galaxy_conn.execute("UPDATE engine_items SET cursor_key = COALESCE(?, cursor_key), done_at = IF(?, NOW(6), done_at), updated_at = NOW(6)"
                                  " WHERE job_id = ? AND attempt = ?", (cursor_key, 1 if done else 0, job_id, attempt)).rowcount
    if touched:
        return
    row = galaxy_conn.execute("SELECT attempt FROM engine_items WHERE job_id = ? FOR UPDATE", (job_id,)).fetchone()     # an update that changed nothing still held the attempt
    if row is None or row["attempt"] != attempt:
        raise Replaced(f"Job {job_id} attempt {attempt} was replaced" + (f" by attempt {row['attempt']}." if row else " (it has no row)."))


def item_state(galaxy_conn, job_id):
    """`{"attempt", "cursor_key", "done"}` of the job's `engine_items` row, or `None` for a job with no row: what the supervisor reads to learn a job is done."""
    row = galaxy_conn.execute("SELECT attempt, cursor_key, done_at FROM engine_items WHERE job_id = ?", (job_id,)).fetchone()
    return None if row is None else {"attempt": row["attempt"], "cursor_key": row["cursor_key"], "done": row["done_at"] is not None}


class Writer:
    """
    Runs ledger calls on a thread of its own, on one connection (`connect()` makes it), so the caller never waits on MySQL. `submit(call)` queues `call(conn)`;
    the writer runs it, commits, and logs (never raises) a call that failed, then goes on. `flush(timeout)` waits for what was queued; `close()` flushes and stops.
    """

    def __init__(self, connect):
        self._connect = connect
        self._queue = queue.Queue()
        self._thread = threading.Thread(target=self._run, name="ledger-writer", daemon=True)
        self._thread.start()

    def submit(self, call):
        self._queue.put(("call", call))

    def flush(self, timeout=30.0):
        """Waits until everything queued so far has been written; returns whether it was within `timeout` seconds."""
        done = threading.Event()
        self._queue.put(("flag", done))
        return done.wait(timeout)

    def close(self, timeout=30.0):
        self.flush(timeout)
        self._queue.put(("stop", None))
        self._thread.join(timeout)

    def _run(self):
        conn = None
        while True:
            what, call = self._queue.get()
            if what == "stop":
                break
            if what == "flag":
                call.set()
                continue
            try:
                if conn is None:
                    conn = self._connect()
                call(conn)
                conn.commit()
            except Exception as exc:        # noqa: BLE001 -- a ledger write never stops the loop that asked for it
                log.swallowed("queue", "writing the job ledger", exc)
                try:
                    if conn is not None:
                        conn.rollback()
                except Exception as bad:        # noqa: BLE001
                    log.swallowed("queue", "rolling back a failed ledger write", bad)
                    conn = None
        if conn is not None:
            conn.close()
