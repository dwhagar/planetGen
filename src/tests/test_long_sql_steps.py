"""OPS.195: a long statement of a command is named in the beat (`planetgen/queue/longsql.py`, `db/statementmark.py`)."""

import threading
import time

from planetgen.db import statementmark, store
from planetgen.queue import longsql, unitbeat


def test_longsql_names_the_statement_while_it_runs_and_not_after(mysql_config):
    conn = store.get_connection(mysql_config)
    seen = []

    def watch():
        deadline = time.monotonic() + 10
        while time.monotonic() < deadline and not seen:
            mark = statementmark.running()
            if mark is not None:
                seen.append(mark)
            time.sleep(0.02)

    watcher = threading.Thread(target=watch)
    try:
        watcher.start()
        longsql.run(conn, "SELECT SLEEP(0.6)", "Sleeping a moment")
        watcher.join()
        assert seen and seen[0][0] == "Sleeping a moment", "the statement is named for the beat while it runs"
        assert statementmark.running() is None, "and not after"
        assert longsql.stage_events_readable(mysql_config) in (True, False)
    finally:
        conn.close()


def test_the_beat_carries_the_work_done_and_estimated_where_the_server_gives_it():
    """Where MySQL's stage events can be read a statement's `progress` gives `(done, estimated)`; a reading that fails or gives nothing adds nothing and never stops the beat."""
    with statementmark.labelled("Building an index", lambda: (30, 120)), statementmark.statement("ALTER TABLE x ADD INDEX i (a)"):
        found = unitbeat.evidence(100.0)
    assert found["sql"] == "Building an index" and (found["sql_done"], found["sql_estimated"]) == (30, 120)
    for progress in (lambda: None, lambda: 1 / 0):
        with statementmark.labelled("Building an index", progress), statementmark.statement("ALTER TABLE x ADD INDEX i (a)"):
            found = unitbeat.evidence(100.0)
        assert "sql_done" not in found and found["sql"] == "Building an index"
