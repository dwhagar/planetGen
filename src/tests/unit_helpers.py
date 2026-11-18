# tests/unit_helpers.py

"""
The units of a run for the tests that build them by hand (OPS.235: RQ is gone; a run's units are rows of the ledger that the supervisor's launcher forks): `ledger_units` cuts a stage into unit rows and releases them for the inline
launcher (`tests.inline_queue`), and `unit_failed` is what a unit's failure does to its run (`units.fail_unit`).
"""

from planetgen.queue import units

UNIT_FUNCTION = units.UNIT_FUNCTION


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
