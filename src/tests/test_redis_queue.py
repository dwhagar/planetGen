# tests/test_redis_queue.py

"""
The work queue's Redis side (PERF.24, step 1): `planetgen.queue.
redisqueue` and the command line of the workers it starts.
"""

from planetgen.queue import redisqueue
from planetgen.util import settings


def test_redis_url_prefers_the_environment(monkeypatch):
    monkeypatch.setenv("PLANETGEN_REDIS_URL", "redis://example:1/2")
    assert redisqueue.redis_url() == "redis://example:1/2"
    monkeypatch.delenv("PLANETGEN_REDIS_URL")
    assert redisqueue.redis_url() == settings.get_settings().redis.url


def test_worker_argv_names_the_worker():
    argv = redisqueue.worker_argv(["q"], "redis://h:1/0", "q-1")
    assert argv[1:] == ["-m", "planetgen.cli.worker", "--burst", "--url", "redis://h:1/0", "--name", "q-1", "q"]
