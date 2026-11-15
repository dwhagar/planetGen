# planetgen/queue/forkedjob.py

"""
A unit as a job process the template forked (OPS.226): the runner `template.run_job` calls in the child, in place of the RQ job `units.run_unit` is on main.

The job is `{"run_id", "stage", "unit", "kind", "args", "cost", "row", "redis_url"}`, what the RQ job carries. The process opens its own Redis connection (the template holds none), runs
the unit through the same code (`units.run_unit_on`: the attempt in the ledger, the fenced done mark, the end of the stage) and says on its pipe how it went. Between two pieces it asks the pipe:
a stop request from the supervisor ends the unit there, the pieces done so far kept, and a supervisor that has gone makes the job finish the piece it is in and exit
(`template.EXIT_SUPERVISOR_GONE`). A failure fails the run as RQ's `on_failure` does (`units.fail_unit`) and ends the job with `template.EXIT_FAILED`.

Built and tested, not yet used: the bulk lane moves to it in OPS.228.
"""

import sys

from planetgen.queue import redisqueue, template, units


def _halt(link):
    """The `halt` of the engine for a job on `link`: why the unit is to stop before its next piece, or `None`."""
    def halt():
        link.requests()
        if link.gone:
            return "the supervisor is gone"
        if link.stop_requested:
            return "the supervisor asked it to stop"
        return None
    return halt


def run(job, link):
    """Runs the unit `job` names and returns the exit code of the job process."""
    conn = redisqueue.connect(job.get("redis_url"))
    run_id, stage, unit, kind = job["run_id"], job["stage"], job["unit"], job["kind"]
    try:
        result = units.run_unit_on(conn, run_id, stage, unit, kind, job["args"], job.get("cost"), job.get("row"), halt=_halt(link))
    except template.SupervisorGone:
        raise
    except Exception:  # noqa: BLE001 -- a unit's failure is the run's failure, as RQ's `on_failure` made it
        exc_type, exc_value, tb = sys.exc_info()
        units.fail_unit(conn, run_id, stage, unit, kind, job.get("retries_left"), exc_type, exc_value, tb, row=job.get("row"))
        link.send({"event": "failed", "error": f"{exc_type.__name__}: {exc_value}"[:500]})
        return template.EXIT_FAILED
    if link.gone:
        return template.EXIT_SUPERVISOR_GONE
    link.send({"event": "unit", "run_id": run_id, "stage": stage, "unit": unit, "result": result})
    return 0
