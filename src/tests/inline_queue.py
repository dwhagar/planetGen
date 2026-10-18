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
import uuid

_active = []
"""list: One entry per context entered; a context inside another does nothing."""


def _server():
    from tests.conftest import _test_server_kwargs
    return _test_server_kwargs()


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

    from planetgen.queue import redisqueue, runrecord, workers

    token = uuid.uuid4().hex[:12]
    prefix = f"planetgen-test:{token}"
    queue_name = f"planetgen-bulk-test-{token}"
    used = []

    def run_waiting_units(name, connection=None, url=None, mysql_host=None, count=None):
        used.append(True)
        connection = connection or redisqueue.connect(url)
        rq.SimpleWorker([redisqueue.queue(name, connection)], connection=connection,
                        serializer=redisqueue.JsonSerializer).work(burst=True, logging_level="WARNING")
        return 0

    server = _server()
    _active.append(True)
    try:
        with contextlib.ExitStack() as stack:
            patch = monkeypatch or stack.enter_context(pytest.MonkeyPatch.context())
            patch.setattr(runrecord, "PREFIX", prefix)
            patch.setattr(redisqueue, "BULK_QUEUE", queue_name)
            patch.setattr(workers, "ensure_workers", run_waiting_units)
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
            client.delete(f"rq:queue:{queue_name}")
