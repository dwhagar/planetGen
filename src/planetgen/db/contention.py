"""
planetgen/db/contention.py

Running one transaction of a unit that may meet a deadlock (MySQL 1213) or a lock wait timeout (1205) (OPS.261).

A unit used to try such a transaction again itself: eight tries with random waits scaled by the process's trouble, then one more with an exclusive grant of the contended step (DB.113). The
job engine does that now for every unit (OPS.220): the failure is classified `passing` (`queue/failures.py`), the unit is run again after a doubling wait inside the run's budget, from the
piece it had finished (OPS.219), and the unit's own retry, its waits and the exclusive grant were deleted. What stays is what a transaction needs to be run again safely: it starts from a clean
transaction (`begin_read_committed`), and a failed one is rolled back and, for a lost connection, explained (`run_once`).

A stage uses it by giving `attempt(conn)`, one whole transaction that commits at its end: `run_once(conn, attempt)`.
"""
import pymysql

from planetgen.db import connection_loss, store


def begin_read_committed(conn):
    """
    Starts a transaction at READ COMMITTED (DB.113, DB.114): the Object ID step of the scatter tables (`scatter_uids.assign`) reads every row that
    has no ID yet, and at the default level that read takes locks that meet the other writers'. Call it first in an `attempt`.
    """
    conn.commit()           # a connection that just read has a transaction open, and MySQL refuses the next line inside one
    conn.execute("SET TRANSACTION ISOLATION LEVEL READ COMMITTED")


def run_once(conn, attempt):
    """
    Runs `attempt(conn)` (one whole transaction, committed at its end) and returns what it returns. After a database error the transaction is rolled back and the error goes on to the engine, which runs
    the unit again when it is a deadlock or a lock wait timeout (module note); a lost connection is explained first (OPS.171: why the server dropped it).

    Raises:
        pymysql.err.OperationalError: Whatever the attempt raised.
    """
    try:
        return attempt(conn)
    except pymysql.err.OperationalError as exc:
        conn.rollback()
        if exc.args and exc.args[0] in store.LOST_CONNECTION_ERRORS and getattr(conn, "_config", None) is not None:
            connection_loss.explain_lost(conn._config, exc)
        raise
