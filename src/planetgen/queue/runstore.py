# planetgen/queue/runstore.py

"""
The durable copy of a run's record (OPS.284, docs/design/operations/job-engine.md).

A run's record (`runrecord.py`) is Redis hashes: the run's own (its title, kind, database, steps, galaxy seed and build, state) and one for each stage (its counters, clocks and conclusion). The rows of
the ledger (`ledger.py`, MySQL) say what ran and what did not, and everything `finish-run` carries over comes from them; the record is what names the run, its steps and the galaxy and build it began on.
Redis does not keep a record through a restart of the machine, so a run that lost its Redis could not be carried on, and its ledger rows stayed `running` for ever.

`persist` writes the record to the control database (`engine_run_records`) at each change that matters (a run's rows are linked, a stage opens, ends or concludes, the run is held, ends or goes on). Redis stays the
cache the pages read: `restore_lost` puts a run back into Redis when Redis does not have it, and `settle_lost` then ends the runs that were live when the machine went (nobody follows them any more: failed,
killed, the ledger's rows ended as for any clear). A run with no ledger rows (begun by hand, a test) has no durable copy: there is nothing to carry on from.
"""

import json
import time

from planetgen.queue import runrecord
from planetgen.util import log

KEEP_DAYS = runrecord.KEEP_SECONDS // 86400
"""int: How long a record is kept after its last change, as long as Redis keeps a finished run; older ones are dropped by `restore_lost`."""

LOST_TEXT = "the machine (or its Redis) restarted and took the run's process and record away before the run ended; the record is back from the control database"
"""str: What a run found live after its record was restored says as its error."""


def _linked(conn, run_id):
    """Whether the run has rows in the ledger (`planrows.engine_run`): only such a run has a durable copy."""
    from planetgen.queue import planrows
    return bool(conn.hget(runrecord.run_key(run_id), planrows.ENGINE_RUN_FIELD))


def _record(conn, run_id):
    """The run's hashes as JSON text: `{"run": {...}, "stages": [{...}]}`, every field and value text."""
    run = runrecord.text_hash(conn.hgetall(runrecord.run_key(run_id)))
    count = len(json.loads(run.get("stages") or "[]"))
    return {"run": run, "stages": [runrecord.text_hash(conn.hgetall(runrecord.stage_key(run_id, index))) for index in range(count)]}


def persist(conn, run_id):
    """
    Writes the run's record to the control database. Best effort: a record that could not be written is said (`log.swallowed`) and the run goes on, with the copy it had; nothing that waits on a run waits on
    this. A run with no ledger rows is left alone.

    Returns:
        bool: Whether the copy was written.
    """
    try:
        if not _linked(conn, run_id):
            return False
        from planetgen.queue import planrows
        record = _record(conn, run_id)
        control = planrows.control_connection(conn, run_id)
        try:
            control.execute(
                "INSERT INTO engine_run_records (run_id, ledger_run, state, record, updated_at) VALUES (?, ?, ?, ?, NOW(6))"
                " ON DUPLICATE KEY UPDATE ledger_run = VALUES(ledger_run), state = VALUES(state), record = VALUES(record), updated_at = NOW(6)",
                (run_id, int(record["run"][planrows.ENGINE_RUN_FIELD]), record["run"].get("state", ""), json.dumps(record, separators=(",", ":"))))
            control.commit()
        finally:
            control.close()
        return True
    except Exception as exc:  # noqa: BLE001 -- see above
        log.swallowed("queue/runstore", "persist", exc)
        return False


def stored(control, run_id):
    """The durable record of `run_id` (`{"run": {...}, "stages": [...]}`), or `None`."""
    row = control.execute("SELECT record FROM engine_run_records WHERE run_id = ?", (run_id,)).fetchone()
    return None if row is None else json.loads(row["record"])


def restore_lost(conn, control=None):
    """
    Puts back into Redis every run whose Redis record is gone (the run hash is the sign) and that Redis would still hold: a live or held run, whatever its age, and any other run for as long as Redis keeps a finished one
    (`runrecord.KEEP_SECONDS`); a run the Queue page deleted has no record (`forget`). Records older than `KEEP_DAYS` are dropped. The fields come back as they
    were written; a run that was live comes back live (`settle_lost` ends it). The ledger's own `engine_run` link is in the record, so the run finds its rows again.

    Args:
        conn: The Redis connection.
        control: The control database connection; the configured one where left out.

    Returns:
        list[str]: The ids of the runs restored.
    """
    own = control is None
    if own:
        from planetgen.db import store
        control = store.get_control_connection(store.control_mysql_config())
    restored = []
    try:
        control.execute("DELETE FROM engine_run_records WHERE updated_at < NOW(6) - INTERVAL ? DAY", (KEEP_DAYS,))
        control.commit()
        for row in control.execute("SELECT run_id, state, record FROM engine_run_records WHERE state IN (?, ?, ?) OR updated_at >= NOW(6) - INTERVAL ? SECOND ORDER BY updated_at",
                                   (*runrecord.LIVE_STATES, runrecord.HELD, runrecord.KEEP_SECONDS)).fetchall():
            if conn.exists(runrecord.run_key(row["run_id"])):
                continue
            record = json.loads(row["record"])
            pipe = conn.pipeline()
            pipe.hset(runrecord.run_key(row["run_id"]), mapping=record["run"])
            for index, stage in enumerate(record["stages"]):
                if stage:
                    pipe.hset(runrecord.stage_key(row["run_id"], index), mapping=stage)
            pipe.execute()
            restored.append(row["run_id"])
    finally:
        if own:
            control.close()
    return restored


def settle_lost(conn, control=None):
    """
    Restores the runs Redis lost (`restore_lost`) and ends each that was live when it went, as a run whose process is gone: `failed` and killed with `LOST_TEXT`, its ledger rows ended (`runcontrol.clear_run`).
    A run that was held, done, failed or canceled comes back as it was and is left alone. Then `settle_ended` ends the ledger rows of the runs that ended in Redis with theirs still open, and `settle_stuck` ends a live run nothing can run in (OPS.314).

    Returns:
        list[str]: The ids of the runs it ended.
    """
    from planetgen.queue import runcontrol
    ended = []
    for run_id in restore_lost(conn, control):
        if runrecord.run_state(conn, run_id) in runrecord.LIVE_STATES and runcontrol.clear_run(conn, run_id, who="the supervisor", reason=LOST_TEXT) is not None:
            ended.append(run_id)
    settle_ended(conn, control)
    ended.extend(settle_stuck(conn, control))
    return ended


def settle_ended(conn, control=None):
    """
    Ends the ledger rows of every run whose Redis record says it is over (`done`, `failed`, `canceled`) while its ledger run is still `running` or `waiting`: a run `finish-run` found with no process (its process was
    killed, or could not even write its own end: a full disk), whose rows nothing else ends (found by the OPS.247 trial, run 13). The unfinished jobs become `canceled` and the run row says so (`controls.end_rows`); a
    finished job is left as it is, and a run that is live or held is not touched.

    Returns:
        list[str]: The ids of the runs whose rows it ended.
    """
    from planetgen.queue import controls
    own = control is None
    if own:
        from planetgen.db import store
        control = store.get_control_connection(store.control_mysql_config())
    ended = []
    try:
        for row in control.execute("SELECT r.run_id, r.ledger_run FROM engine_run_records r JOIN engine_runs e ON e.id = r.ledger_run WHERE e.state IN ('running', 'waiting')").fetchall():
            if runrecord.run_state(conn, row["run_id"]) in runrecord.FINAL_STATES:
                controls.end_rows(control, int(row["ledger_run"]))
                ended.append(row["run_id"])
        control.commit()
    finally:
        if own:
            control.close()
    return ended


STUCK_FAILED_SECONDS = 600.0
"""float: How long a live run may stand with a failed job and nothing running or ready before it is ended as failed (OPS.314)."""

STUCK_WAITING_SECONDS = 1800.0
"""float: How long a live run may stand with jobs waiting and nothing running, ready or failed before it is ended as failed (OPS.314)."""

_idle_since = {}
"""dict: `{run id: first time it was seen with nothing able to run}`, in this process; a run that moves again, or ends, is dropped."""


def stuck_reason(control, ledger_run, counts, seconds):
    """What to say of a run nothing can run in: the first failed job with its error, or the jobs that wait on something that never finishes, with how long it stood (HMS)."""
    from planetgen.util import logrecord
    stood = logrecord.duration_text(seconds)
    failed = control.execute("SELECT j.stage, j.command, j.item_key, a.error_text FROM engine_jobs j LEFT JOIN engine_attempts a ON a.job_id = j.id AND a.attempt = j.attempt "
                             "WHERE j.run_id = ? AND j.state = 'failed' ORDER BY j.id LIMIT 1", (ledger_run,)).fetchone()
    total = sum(counts.values())
    done = counts.get("done", 0)
    if failed is not None:
        where = f"{failed['command']}" + (f" {failed['item_key']}" if failed["item_key"] else "") + f" of stage {failed['stage']}"
        text = (failed["error_text"] or "no error was kept").strip().splitlines()[0][:300] if (failed["error_text"] or "").strip() else "no error was kept"
        return (f"nothing could run for {stood}: job {where} failed ({text}), no job was running, and {counts.get('waiting', 0)} job(s) wait on it ({done} of {total} jobs done); "
                f"nothing was following the run, so it is ended here")
    return (f"nothing could run for {stood}: {counts.get('waiting', 0)} job(s) wait on a job that never finishes, none is running or ready ({done} of {total} jobs done); "
            f"nothing was following the run, so it is ended here")


def settle_stuck(conn, control=None, now=None):
    """
    OPS.314: ends a live run in which nothing can run any more. A run whose ledger rows have no job `running` or `ready`, and have a job `failed` (or only jobs `waiting`), stands for
    `STUCK_FAILED_SECONDS` (`STUCK_WAITING_SECONDS` with no failed job) and then ends `failed` with the cause said in plain words (`stuck_reason`), its rows ended (`controls.end_rows`), so
    the Queue page does not show it Running for days with no worker alive. A run that has a job running or ready, is held or is over is left alone.

    Returns:
        list[str]: The ids of the runs it ended.
    """
    import time
    from planetgen.queue import controls
    now = time.time() if now is None else now
    own = control is None
    if own:
        from planetgen.db import store
        control = store.get_control_connection(store.control_mysql_config())
    ended = []
    try:
        for row in control.execute("SELECT r.run_id, r.ledger_run FROM engine_run_records r JOIN engine_runs e ON e.id = r.ledger_run WHERE e.state IN ('running', 'waiting')").fetchall():
            run_id, ledger_run = row["run_id"], int(row["ledger_run"])
            if runrecord.run_state(conn, run_id) not in runrecord.LIVE_STATES:
                _idle_since.pop(run_id, None)
                continue
            counts = {r["state"]: int(r["n"]) for r in control.execute("SELECT state, COUNT(*) AS n FROM engine_jobs WHERE run_id = ? GROUP BY state", (ledger_run,)).fetchall()}
            if counts.get("running") or counts.get("ready") or not (counts.get("failed") or counts.get("waiting")):
                _idle_since.pop(run_id, None)
                continue
            first = _idle_since.setdefault(run_id, now)
            limit = STUCK_FAILED_SECONDS if counts.get("failed") else STUCK_WAITING_SECONDS
            if now - first < limit:
                continue
            reason = stuck_reason(control, ledger_run, counts, now - first)
            runrecord.finish_run(conn, run_id, "failed", error=reason)
            controls.end_rows(control, ledger_run)
            _idle_since.pop(run_id, None)
            ended.append(run_id)
        control.commit()
    finally:
        if own:
            control.close()
    return ended


def forget(conn, run_id):
    """Drops the durable copy of a run the Queue page deleted (`runrecord.delete_run`), so it is not restored. Best effort."""
    try:
        from planetgen.queue import planrows
        control = planrows.control_connection(conn, run_id)
        try:
            control.execute("DELETE FROM engine_run_records WHERE run_id = ?", (run_id,))
            control.commit()
        finally:
            control.close()
    except Exception as exc:  # noqa: BLE001 -- the Redis record is gone; a copy left is dropped when it is old
        log.swallowed("queue/runstore", "forget", exc)
