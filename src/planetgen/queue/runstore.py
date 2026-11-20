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
    A run that was held, done, failed or canceled comes back as it was and is left alone.

    Returns:
        list[str]: The ids of the runs it ended.
    """
    from planetgen.queue import runcontrol
    ended = []
    for run_id in restore_lost(conn, control):
        if runrecord.run_state(conn, run_id) in runrecord.LIVE_STATES and runcontrol.clear_run(conn, run_id, who="the supervisor", reason=LOST_TEXT) is not None:
            ended.append(run_id)
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
