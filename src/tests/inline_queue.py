# tests/inline_queue.py

"""
The work queue run inside the test process (PERF.97, PERF.106): a context that gives the code in it a queue and a key
prefix of its own and has its `ensure_workers` run the waiting units right here, so no worker process is left running
when the test's database goes and a patch the test made reaches the units. The workers read the MySQL server's
credentials from the settings, so the environment carries the test server's. Nothing is patched, and nothing is
touched on Redis, without a Redis server.

The autouse fixture of `conftest.py` applies it to every test; a fixture of wider scope than a test (a module's galaxy,
built by `run_cli`) runs before that one, so `bughunt_support.run_cli` applies it itself.
"""

import contextlib
import os
import signal
import uuid

import pymysql

_active = []
"""list: One entry per context entered; a context inside another does nothing."""


def _server():
    from tests.conftest import _test_server_kwargs
    return _test_server_kwargs()


class _Link:
    """What `forkedjob.run` takes of its pipe, for a unit run in this process: nothing asks it to stop, the supervisor is never gone, what it says is kept."""

    gone = False
    stop_requested = False

    def __init__(self):
        self.sent = []

    def requests(self):
        return []

    def send(self, message):
        self.sent.append(message)


class _InlineJob:
    def __init__(self, job_id, link):
        self.id, self.pid, self.link = job_id, os.getpid(), link

    def messages(self):
        said, self.link.sent = list(self.link.sent), []
        return said

    def close(self):
        pass


class InlineForker:
    """The forker of a `launcher.Launcher` for the tests (OPS.228): `fork` runs the unit at once, in this process, through `forkedjob.run` as a job process would, so a ledger run's units are done where an RQ
    unit is: when the follower asks for workers."""

    def __init__(self):
        self.exits, self.made = [], 0

    def fork(self, kind, runner, job=None, template=None):
        from planetgen.queue import forkedjob, jobhost
        self.made += 1
        link = _Link()
        paused = signal.alarm(0)            # a test's own time limit (tests/test_fuzz_cli.py) is for the command, as it was with the RQ worker beside it: the unit's time is not counted
        try:
            code = forkedjob.run(job, link)
        finally:
            if paused:
                signal.alarm(paused)
        self.exits.append(("inline", jobhost.Exit(self.made, os.getpid(), code=code)))
        return "inline", _InlineJob(self.made, link)

    def poll(self, timeout=0.0):
        found, self.exits = self.exits, []
        return found


class InlineLightForker:
    """The forker of a `lightjobs.LightJobs` for the tests (OPS.230): `fork` runs the light job at once, in this process, through `lightjobs.run` as the job process would, so an API job or a slow page's
    answer is done where an RQ job was: when whoever submitted it, or asked after it, asks for the supervisor."""

    def __init__(self):
        self.exits, self.made = [], 0

    def fork(self, kind, runner, job=None, owner=None, template=None):
        from planetgen.queue import jobhost, lightjobs
        self.made += 1
        link = _Link()
        paused = signal.alarm(0)            # a test's own time limit is for the command: the job's time is not counted
        try:
            code = lightjobs.run({**(job or {}), "nice": 0}, link)       # this process is the test's: it is not made nicer
        finally:
            lightjobs._current.update(id=None, url=None)
            if paused:
                signal.alarm(paused)
        self.exits.append(("inline", jobhost.Exit(self.made, os.getpid(), code=code)))
        return "inline", _InlineJob(self.made, link)

    def poll(self, timeout=0.0, owner=None):
        found, self.exits = self.exits, []
        return found

    def stop(self):
        pass


def run_light_jobs(connection=None, url=None):
    """Hands out the pending light jobs and runs each here until none is waiting (the supervisor's `LightJobs` duty with an inline forker); returns how many it ran."""
    from planetgen.queue import lightjobs, memorygate, redisqueue
    connection = connection or redisqueue.connect(url)
    duty = lightjobs.LightJobs(connection, InlineLightForker(), read_host=lambda: memorygate.Host(10 ** 6, 10 ** 6, 64, 0.0), held_mb=lambda pid: 0.0)
    ran = 0
    while True:
        started = duty.turn()
        ran += len(started)
        if not started:
            break
    duty.ended()
    return ran


def run_ready_rows():
    """Hands out the ready rows of the job ledger and runs each unit here until none is ready (the launcher with an inline forker); returns how many units it ran. Nothing happens for a run with no ledger."""
    from planetgen.db import store
    from planetgen.queue import launcher, memorygate
    try:
        control = store.get_control_connection()
    except Exception:  # noqa: BLE001 -- a test with no control database has no ledger rows to hand out
        return 0
    forker = InlineForker()
    one = launcher.Launcher(control, forker, read_host=lambda: memorygate.Host(10 ** 6, 10 ** 6, 64, 0.0))
    ran = 0
    try:
        while True:
            started = one.turn()
            ran += len(started)
            if not started:
                break
        one.ended()
    except pymysql.err.ProgrammingError as exc:
        if exc.args[0] != 1146:             # the control database of this test has no ledger tables: no run of it has rows
            raise
    finally:
        control.close()
    return ran


@contextlib.contextmanager
def inline_queue(monkeypatch=None):
    """
    The units of a run started inside the context run in this process, on a queue of the context's own.

    Args:
        monkeypatch (pytest.MonkeyPatch, optional): The test's own `monkeypatch` fixture, which then undoes the patches
            too, in the order they were made: a test that patches the same names after this has the fixture put back
            what it found, and a patch undone earlier would leave that behind. By default the context has its own.
    """
    url = os.environ.get("PLANETGEN_REDIS_URL")
    if not url or _active:
        yield
        return
    import pytest
    import redis
    import rq

    from planetgen.queue import allowlist, lightjobs, memorygate, redisqueue, runrecord, workers, wrapper

    token = uuid.uuid4().hex[:12]
    prefix = f"planetgen-test:{token}"
    queue_name = f"planetgen-bulk-test-{token}"
    maintenance_name = f"planetgen-maintenance-test-{token}"
    used = []

    def run_waiting_units(name, connection=None, url=None, mysql_host=None, count=None):
        used.append(True)
        connection = connection or redisqueue.connect(url)
        rq.SimpleWorker([redisqueue.queue(name, connection)], connection=connection, serializer=redisqueue.JsonSerializer,
                        job_class=allowlist.GuardedJob).work(burst=True, logging_level="WARNING")     # SEC.45: the real worker's guard
        run_ready_rows()            # OPS.228, OPS.229: the units a ledger run released for the launcher, whatever its kind
        return 0

    real_stage_ensure = wrapper.stage_ensure

    def inline_stage_ensure(conn, url, kind, step_kind, workers_count=None, run_id=None):
        """The `ensure` of a stage, and, for a run whose units go to the launcher (no worker to ask for), the hand-out of its ready rows here: the follower's ask for workers is when the units run."""
        ensure = real_stage_ensure(conn, url, kind, step_kind, workers_count, run_id)

        def ask():
            started = ensure()
            run_ready_rows()                # OPS.229: a warm-up or a check has its units handed out by the launcher too
            used.append(True)
            return started
        return ask

    server = _server()
    _active.append(True)
    try:
        with contextlib.ExitStack() as stack:
            patch = monkeypatch or stack.enter_context(pytest.MonkeyPatch.context())
            patch.setattr(runrecord, "PREFIX", prefix)
            patch.setattr(memorygate, "PREFIX", f"{prefix}:gate")
            patch.setattr(redisqueue, "BULK_QUEUE", queue_name)
            patch.setattr(redisqueue, "MAINTENANCE_QUEUE", maintenance_name)
            patch.setattr(workers, "ensure_workers", run_waiting_units)
            patch.setattr(wrapper, "stage_ensure", inline_stage_ensure)
            patch.setattr(lightjobs, "ask_supervisor", lambda connection=None, url=None: run_light_jobs(connection, url))        # OPS.230: the supervisor's duty, in this process, at once
            for name, value in (("HOST", server["host"]), ("PORT", server["port"]), ("USER", server["user"]),
                                ("PASSWORD", server["password"])):
                patch.setenv(f"PLANETGEN_MYSQL_{name}", str(value))
            yield
    finally:
        _active.pop()
        if used:
            client = redis.Redis.from_url(url)
            for key in list(client.scan_iter(match=f"{prefix}:*")) + list(client.scan_iter(match=f"rq:*{token}*")):
                client.delete(key)
            client.delete(f"rq:queue:{queue_name}", f"rq:queue:{maintenance_name}")
