"""OPS.226: a unit as a job process the template forked (`queue/forkedjob.py`): the same run ends the same way, with the units run in forked jobs and not in RQ's horses, and a unit stopped or
failed ends its job as the engine says. Real MySQL, Redis and processes."""

import pytest
import rq

from planetgen.queue import allowlist, jobhost, memorygate, redisqueue, runrecord, units, workers
from planetgen.cli import reset as reset_cli
from tests import inline_queue
from tests.test_one_unit_commands import _place_sector

pytestmark = [pytest.mark.fresh_database]

R = "tests.template_runners"


@pytest.fixture
def forking(monkeypatch):
    """The test's queue with the units of a run taken off it and run in jobs forked from a template; the other jobs (the openings) run in a worker as before."""
    monkeypatch.setenv(allowlist.EXTRA_ENV, f"{R}:*")
    host = jobhost.JobHost()
    host.start(ready_seconds=60)
    forked = []
    inline = workers.ensure_workers          # the test queue's stand-in for a worker (`tests/inline_queue.py`)

    def ensure(name, connection=None, url=None, mysql_host=None, count=None):
        connection = connection or redisqueue.connect(url)
        queue = redisqueue.queue(name, connection)
        while True:
            waiting = [job for job in queue.jobs if job.func_name == units.UNIT_FUNCTION]
            if not waiting:
                break
            for job in waiting:
                queue.remove(job)
                run_id, stage, unit, kind, args, cost, row = job.args
                handed = {"run_id": run_id, "stage": stage, "unit": unit, "kind": kind, "args": args, "cost": cost, "row": row, "redis_url": url,
                          "prefix": runrecord.PREFIX, "gate_prefix": memorygate.PREFIX, "bulk_queue": redisqueue.BULK_QUEUE, "maintenance_queue": redisqueue.MAINTENANCE_QUEUE}
                held = host.fork(f"{R}.unit_in_the_tests_queue", handed)
                exit_ = []
                while not exit_:
                    exit_ = host.poll(0.5)
                forked.append((unit, exit_[0].code, held.messages()))
                held.close()
            inline(name, connection, url, mysql_host, count)           # whatever else is waiting: the opening the unit's end queued
        return 0

    monkeypatch.setattr(workers, "ensure_workers", ensure)
    yield forked
    host.stop()


def test_a_run_whose_units_are_forked_jobs_ends_as_a_run_in_a_worker_does(mysql_config, redis_server, forking):
    _place_sector(mysql_config, "One", (5.0, 5.0, 5.0))
    assert reset_cli.reset_database(mysql_config, assume_yes=True)
    conn = redisqueue.connect(redis_server)
    ours = [run for run in (runrecord.read_run(conn, run_id) for run_id in runrecord.recent_runs(conn)) if run and run["database"] == mysql_config.database]
    assert [(run["title"], run["state"]) for run in ours] == [("Reset the galaxy", "done")]
    assert [(unit, code) for unit, code, _said in forking] == [(0, 0)], "the unit ran in a forked job, once, and ended well"
    said = forking[0][2]
    assert said and said[0]["event"] == "unit" and said[0]["unit"] == 0 and "seconds" in said[0]["result"]


class _Link:
    """A job's pipe as a test holds it: what the supervisor asked and whether it is gone."""

    def __init__(self, stop=False, gone=False):
        self.stop_requested, self.gone, self.said = stop, gone, []

    def requests(self):
        return []

    def send(self, message):
        self.said.append(message)


def test_the_engine_stops_a_unit_between_two_pieces_when_halt_says_why():
    from planetgen.queue import commands, engine
    done = []
    spec = commands.CommandSpec("pieces-test", "A test", pieces=(__name__, "_pieces"), do=(__name__, "_do"))   # named for the guards; the functions are patched in
    asked = iter([None, None, "the supervisor asked it to stop", None])
    import pytest as _pytest
    with _pytest.MonkeyPatch.context() as patch:
        patch.setitem(commands.REGISTRY, "pieces-test", spec)
        patch.setattr(commands, "pieces_function", lambda kind: lambda args, cursor: [0, 1, 2, 3])
        patch.setattr(commands, "do_function", lambda kind: lambda args, piece, ctx: done.append(piece))
        patch.setattr(commands, "gather_function", lambda kind: None)
        with pytest.raises(engine.UnitStopped, match="asked it to stop"):
            engine.run("pieces-test", {}, lambda done_, total: None, halt=lambda: next(asked))
    assert done == [0, 1], "the two pieces before the request were done, and the third was not begun"


def test_a_job_halts_for_a_stop_request_and_for_a_supervisor_that_is_gone():
    from planetgen.queue import forkedjob
    assert forkedjob._halt(_Link())() is None
    assert forkedjob._halt(_Link(stop=True))() == "the supervisor asked it to stop"
    assert forkedjob._halt(_Link(gone=True))() == "the supervisor is gone"


def test_a_unit_that_raised_fails_its_run_as_rqs_failure_callback_does_and_the_job_ends_failed(monkeypatch, redis_server):
    from planetgen.queue import forkedjob, template
    failed = []

    def boom(*args, **kwargs):
        raise ValueError("the unit broke")

    monkeypatch.setattr(units, "run_unit_on", boom)
    monkeypatch.setattr(units, "fail_unit", lambda conn, run_id, stage, unit, kind, retries, exc_type, exc_value, tb: failed.append((run_id, stage, unit, kind, exc_type.__name__, str(exc_value))) or True)
    link = _Link()
    code = forkedjob.run({"run_id": "r1", "stage": 2, "unit": 5, "kind": "scatter-mass", "args": {}, "redis_url": redis_server}, link)
    assert code == template.EXIT_FAILED and failed == [("r1", 2, 5, "scatter-mass", "ValueError", "the unit broke")]
    assert link.said == [{"event": "failed", "error": "ValueError: the unit broke"}]


def test_a_unit_stopped_because_the_supervisor_went_ends_its_job_with_the_gone_code(monkeypatch, redis_server):
    from planetgen.queue import forkedjob, template
    monkeypatch.setattr(units, "run_unit_on", lambda *args, **kwargs: {"stopped": True})
    assert forkedjob.run({"run_id": "r1", "stage": 0, "unit": 0, "kind": "x", "args": {}, "redis_url": redis_server}, _Link(gone=True)) == template.EXIT_SUPERVISOR_GONE
