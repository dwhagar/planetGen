# tests/test_private_redis.py

"""The Redis database a test worker runs its tests and the processes they start on (TEST.266, TEST.280) is its own and never the server's."""

import pytest

from tests.private_redis import url_for_worker


@pytest.mark.parametrize("worker, database", [("gw0", 1), ("gw1", 2), ("gw3", 4), ("master", 1), ("gw15", 1)])
def test_a_worker_takes_the_database_numbered_after_it(worker, database):
    assert url_for_worker("redis://127.0.0.1:6379/0", worker) == f"redis://127.0.0.1:6379/{database}"


def test_the_server_urls_own_database_is_never_chosen():
    assert url_for_worker("redis://127.0.0.1:6379/2", "gw1") != "redis://127.0.0.1:6379/2"
    assert url_for_worker("redis://127.0.0.1:6379", "gw0") == "redis://127.0.0.1:6379/1"


def test_two_workers_never_share_a_database():
    urls = {url_for_worker("redis://127.0.0.1:6379/0", f"gw{n}") for n in range(8)}
    assert len(urls) == 8


def test_the_tests_of_a_worker_run_on_the_database_numbered_after_it():
    """TEST.280: a process a test starts reads `PLANETGEN_REDIS_URL`, so it is the worker's own database, not the server's."""
    import os

    server = os.environ.get("PLANETGEN_TEST_REDIS_URL")
    if not server:
        pytest.skip("PLANETGEN_TEST_REDIS_URL is not set")
    assert os.environ["PLANETGEN_REDIS_URL"] == url_for_worker(server)
    assert os.environ["PLANETGEN_REDIS_URL"] != server
