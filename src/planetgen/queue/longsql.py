# planetgen/queue/longsql.py

"""
One long SQL statement of a command, shown in the beat (OPS.195; reports/job-engine-items-stage-1-2-2026-11-13.md, item P5).

`run(conn, statement, label)` sends one statement the way `conn.execute` does and, while it runs, names it in the beat: the label and the time it was sent (`db/statementmark.py`, which
`unitbeat.evidence` writes every 10 s), and the unit's phase (`unitbeat.phase`, so the Queue page says what the unit is doing). A clear or a drop that takes minutes is then a unit that
says "in SQL since 6:40: truncate scatter_masses", not one that went quiet. The engine's `ctx.sql(statement, label)` (OPS.191) is this function.

**MySQL's stage events are behind a check, and off as installed.** Where the account can read `performance_schema.events_stages_current` the beat also carries `sql_done` and `sql_estimated`
(MySQL's own work done and work estimated for the statement, which it keeps for some statements, such as an `ALTER TABLE` or an index build). The account INSTALL.md creates cannot read it
(`GRANT ALL ON planetgen%.*` gives no `SELECT` on `performance_schema`), and the stage consumers are off by default, so `stage_events_readable` says no and nothing more is asked of the server.
"""

import pymysql

from planetgen.db import statementmark, store
from planetgen.queue import unitbeat
from planetgen.util import log

_readable = {}      # database key -> whether the stage events can be read

STAGE_SQL = ("SELECT e.WORK_COMPLETED, e.WORK_ESTIMATED FROM performance_schema.events_stages_current e "
             "JOIN performance_schema.threads t ON t.THREAD_ID = e.THREAD_ID WHERE t.PROCESSLIST_ID = ?")


def stage_events_readable(config):
    """Whether the account of `config` can read MySQL's current stage events; asked once for each database and kept."""
    key = (config.host, config.database, config.user)
    if key not in _readable:
        conn = store.get_connection(config, ensure_schema=False)
        try:
            conn.execute(STAGE_SQL, (0,)).fetchall()
            _readable[key] = True
        except pymysql.err.MySQLError as exc:
            log.swallowed("queue", "reading MySQL's stage events (not readable as installed; a statement then shows its label and start only)", exc)
            _readable[key] = False
        finally:
            conn.close()
    return _readable[key]


def _stage_reader(config, connection_id):
    """A function that gives `(done, estimated)` for the statement the connection `connection_id` runs, from a second short connection (the statement's own is busy); `None` where it has none."""
    def read():
        other = store.get_connection(config, ensure_schema=False)
        try:
            row = other.execute(STAGE_SQL, (connection_id,)).fetchone()
        finally:
            other.close()
        if not row:
            return None
        done, estimated = list(row.values()) if hasattr(row, "values") else row
        return (int(done), int(estimated)) if estimated else None
    return read


def run(conn, statement, label, params=(), config=None):
    """
    Sends `statement` on `conn` (as `conn.execute(statement, params)`) with `label` shown in the beat while it runs; returns the cursor.

    Args:
        conn: A `db.store` connection.
        statement (str): The statement.
        label (str): What it does, in the words the Queue page shows ("truncate scatter_masses").
        params (tuple): Its parameters.
        config (MySQLConfig, optional): The database's config; where given and MySQL's stage events can be read, the beat carries the statement's work done and estimated.
    """
    progress = None
    if config is not None and stage_events_readable(config):
        progress = _stage_reader(config, conn.execute("SELECT CONNECTION_ID() AS id").fetchone()["id"])
    unitbeat.phase(label)
    with statementmark.labelled(label, progress):
        return conn.execute(statement, params)
