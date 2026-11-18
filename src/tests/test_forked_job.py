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


class _RealForker:
    """The forker of `inline_queue.run_ready_rows` for this test: each ready row is forked as a job process from a real template (`jobhost.JobHost`) and waited for, so the launcher finds it ended at its next look."""

    def __init__(self, host, record, url):
        self.host, self.record, self.url, self.exits, self.made = host, record, url, [], 0

    def fork(self, kind, runner, job=None):
        handed = dict(job, redis_url=self.url, prefix=runrecord.PREFIX, gate_prefix=memorygate.PREFIX, bulk_queue=redisqueue.BULK_QUEUE, maintenance_queue=redisqueue.MAINTENANCE_QUEUE)
        held = self.host.fork(f"{R}.unit_in_the_tests_queue", handed)      # the key prefix and queue names of the test's own queue go with the job
        exit_ = []
        while not exit_:
            exit_ = self.host.poll(0.5)
        said = held.messages()
        self.record.append((job.get("unit"), exit_[0].code, said))        # an opening has no unit
        self.exits.append(("t", exit_[0]))
        return "t", _Said(held, said)

    def poll(self, timeout=0.0):
        found, self.exits = self.exits, []
        return found


class _Said:
    """A forked job as the launcher reads it after it ended: what it said is what was read when it ended."""

    def __init__(self, held, said):
        self.held, self.said, self.pid, self.id = held, said, held.pid, held.id

    def messages(self):
        said, self.said = self.said, []
        return said

    def close(self):
        self.held.close()


@pytest.fixture
def forking(monkeypatch, redis_server):
    """The units of a run handed out by the launcher and forked as job processes from a real template (OPS.226, OPS.228); each is recorded as `(unit, exit code, what it said)`."""
    from tests import inline_queue
    monkeypatch.setenv(allowlist.EXTRA_ENV, f"{R}:*")
    from planetgen.db import store
    control = store.get_control_connection()
    try:        # rows another test of the session left ready are not this test's units (every run is a ledger run, OPS.271)
        control.execute("UPDATE engine_jobs SET state = 'canceled' WHERE state IN ('ready', 'waiting', 'running')")
        control.commit()
    finally:
        control.close()
    host = jobhost.JobHost()
    host.start(ready_seconds=60)
    forked = []
    monkeypatch.setattr(inline_queue, "FORKER", lambda: _RealForker(host, forked, redis_server))
    yield forked
    host.stop()


def test_a_run_whose_units_are_forked_jobs_ends_as_a_run_in_a_worker_does(mysql_config, redis_server, forking):
    _place_sector(mysql_config, "One", (5.0, 5.0, 5.0))
    assert reset_cli.reset_database(mysql_config, assume_yes=True)
    conn = redisqueue.connect(redis_server)
    ours = [run for run in (runrecord.read_run(conn, run_id) for run_id in runrecord.recent_runs(conn)) if run and run["database"] == mysql_config.database]
    assert [(run["title"], run["state"]) for run in ours] == [("Reset the galaxy", "done")]
    units_run = [entry for entry in forking if entry[0] is not None]
    assert [(unit, code) for unit, code, _said in units_run] == [(0, 0)], "the unit ran in a forked job, once, and ended well"
    said = units_run[0][2]
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
    monkeypatch.setattr(units, "fail_unit", lambda conn, run_id, stage, unit, kind, retries, exc_type, exc_value, tb, row=None: failed.append((run_id, stage, unit, kind, exc_type.__name__, str(exc_value))) or True)
    link = _Link()
    code = forkedjob.run({"run_id": "r1", "stage": 2, "unit": 5, "kind": "scatter-mass", "args": {}, "redis_url": redis_server}, link)
    assert code == template.EXIT_FAILED and failed == [("r1", 2, 5, "scatter-mass", "ValueError", "the unit broke")]
    assert link.said == [{"event": "failed", "error": "ValueError: the unit broke"}]


def test_a_unit_stopped_because_the_supervisor_went_ends_its_job_with_the_gone_code(monkeypatch, redis_server):
    from planetgen.queue import forkedjob, template
    monkeypatch.setattr(units, "run_unit_on", lambda *args, **kwargs: {"stopped": True})
    assert forkedjob.run({"run_id": "r1", "stage": 0, "unit": 0, "kind": "x", "args": {}, "redis_url": redis_server}, _Link(gone=True)) == template.EXIT_SUPERVISOR_GONE


def test_the_opening_job_runs_the_opening_and_says_how_it_went(monkeypatch, redis_server):
    from planetgen.queue import forkedjob, template
    seen = []
    monkeypatch.setattr(units, "run_opening_on", lambda conn, run_id, index: seen.append((run_id, index)) or {"opened": True})
    link = _Link()
    assert forkedjob.run({"opening": True, "run_id": "r1", "stage": 3, "redis_url": redis_server}, link) == 0
    assert seen == [("r1", 3)] and link.said == [{"event": "unit", "run_id": "r1", "stage": 3, "unit": None, "result": {"opened": True}}]

    def boom(conn, run_id, index):
        raise ValueError("the cut broke")

    monkeypatch.setattr(units, "run_opening_on", boom)
    link = _Link()
    assert forkedjob.run({"opening": True, "run_id": "r1", "stage": 3, "redis_url": redis_server}, link) == template.EXIT_FAILED
    assert link.said == [{"event": "failed", "error": "ValueError: the cut broke"}]
