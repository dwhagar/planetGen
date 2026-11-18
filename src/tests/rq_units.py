# tests/rq_units.py

"""
A run's units as RQ jobs, for the tests of what still reads them there: the stuck-unit watch (`unitstuck`), the cancel of a run's jobs (`runcontrol`) and the Queue page's old job view. A run's units are rows of the job
ledger that the supervisor's launcher forks (OPS.228, OPS.271) and no product code queues them on RQ any more; these tests build the RQ job by hand until OPS.234 and OPS.235 move what they test to the launcher.
"""

from planetgen.queue import redisqueue, units

UNIT_FUNCTION = units.UNIT_FUNCTION


def enqueue_units(conn, run_id, stage, kind, unit_list, queue_name=None, expected_seconds=None, rows=None):
    """Queues the units of stage `stage` of the run as RQ jobs named by `UNIT_FUNCTION`, with the arguments the old `units.enqueue_units` gave them (`run_id, stage, index, kind, args, cost, row`)."""
    queue = redisqueue.queue(queue_name or redisqueue.BULK_QUEUE, conn)
    timeout = units.unit_timeout(expected_seconds)
    queue.enqueue_many([queue.prepare_data(UNIT_FUNCTION, args=(run_id, stage, index, kind, args, cost, None if rows is None else rows[index]), job_id=units.job_id(run_id, stage, index),
                                           timeout=timeout, result_ttl=0, failure_ttl=86400, description=f"{kind} {index}")
                        for index, (args, cost) in enumerate(unit_list)])


def unit_failed(job, connection, exc_type, exc_value, traceback):
    """What RQ's `on_failure` of a unit did, for the tests of what a unit's failure does to its run: `units.fail_unit` with the arguments of the job, as `forkedjob.run` gives them."""
    run_id, stage, unit = job.args[:3]
    kind = job.args[3] if len(job.args) > 3 else "?"
    row = job.args[6] if len(job.args) > 6 else None
    return units.fail_unit(connection, run_id, stage, unit, kind, getattr(job, "retries_left", None), exc_type, exc_value, traceback, row=row)


def ledger_units(conn, run_id, stage, kind, unit_list, expected_seconds=None, rows=None):
    """Cuts stage `stage` of the run into the ledger's unit rows and releases them, as `wrapper._open` does, for the tests that run the units with the inline launcher (`tests.inline_queue`)."""
    from planetgen.queue import planrows
    planrows.record_cut(conn, run_id, stage, kind, unit_list, sum(cost for _args, cost in unit_list), {})
    planrows.release_units(conn, run_id, stage)
