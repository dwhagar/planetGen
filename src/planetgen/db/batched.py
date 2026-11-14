# db/batched.py

"""
Long statements in batches that count (OPS.185).

A clear of a scatter is a few statements over tables of millions of rows (the `sector_stats` levels, the stars' member trees, the spheres the star
draw kept). Each one as a single statement is minutes in which a unit's count does not move, and the heartbeat of the unit (`queue/unitbeat.py`) can say
only that nothing happened. `repeat` runs such a statement a batch of rows at a time: each batch commits, tells the running unit it moved
(`unitbeat.moving`, which does nothing outside a unit) and feeds a `Ticker` that says the rows done so far, so the unit's line is a real count and a unit that
works never looks quiet.

    removed = batched.repeat(conn, "DELETE FROM sky_tier_draws", label="Emptying the sky tier draws")
"""

from planetgen.util import ticker

BATCH_ROWS = 100_000
"""int: How many rows one batch of `repeat` changes before it commits and says where it is (defined by the program, not a physical value: a few seconds of a statement on the tables of a
real galaxy, short enough that a reader of the table never waits and long enough that the commits are no cost)."""


def moving():
    """The running unit has made progress its count does not show: the beat's quiet clock starts again (nothing happens outside a unit)."""
    from planetgen.queue import unitbeat      # imported here, as `db/contention.py` does: the queue package reads the database layer

    unitbeat.moving()


def repeat(conn, sql, args=(), *, label, batch=BATCH_ROWS, pace=None):
    """
    Runs `sql` again and again until it changes no row, `batch` rows at a time: a `DELETE`, or an `UPDATE` whose change takes the row out of its own
    `WHERE` (`SET level = -1 WHERE level > 0`), written without a `LIMIT` (it is added). Each batch is committed, moves the running unit's beat and, at most every
    `PROGRESS_LOG_SECONDS`, says `label` with the rows changed so far.

    Args:
        conn (Connection): An open connection.
        sql (str): The statement, with `?` marks, no `LIMIT`.
        args (sequence): Its values.
        label (str): What the line says is going on, `Emptying the sky tier draws`.
        batch (int): Rows per batch.
        pace (Ticker, optional): The ticker to feed, for a caller whose several statements are one count; a new one by default.

    Returns:
        int: The rows changed in all.
    """
    pace = pace if pace is not None else ticker.Ticker(label)
    done = 0
    while True:
        changed = conn.execute(f"{sql} LIMIT {int(batch)}", tuple(args)).rowcount
        if not changed:
            return done
        done += changed
        conn.commit()
        moving()
        pace.tick(done)
