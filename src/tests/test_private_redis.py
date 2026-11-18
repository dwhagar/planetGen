# tests/test_private_redis.py

"""The Redis database a test worker runs its tests and the processes they start on (TEST.266, TEST.280) is its own and never the server's."""

import pytest

from tests import private_redis
from tests.private_redis import url_for_worker


@pytest.mark.parametrize("worker, database", [("gw0", 1), ("gw1", 2), ("gw3", 4), ("master", 1), ("gw15", 1)])
def test_a_worker_takes_the_database_numbered_after_it(worker, database, monkeypatch):
    monkeypatch.delenv(private_redis.SLOT_ENVIRONMENT, raising=False)
    assert url_for_worker("redis://127.0.0.1:6379/0", worker) == f"redis://127.0.0.1:6379/{database}"


def test_the_server_urls_own_database_is_never_chosen(monkeypatch):
    monkeypatch.delenv(private_redis.SLOT_ENVIRONMENT, raising=False)
    assert url_for_worker("redis://127.0.0.1:6379/2", "gw1") != "redis://127.0.0.1:6379/2"
    assert url_for_worker("redis://127.0.0.1:6379", "gw0") == "redis://127.0.0.1:6379/1"


def test_two_workers_never_share_a_database(monkeypatch):
    monkeypatch.delenv(private_redis.SLOT_ENVIRONMENT, raising=False)
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


def test_a_slot_holds_three_databases_and_the_slots_never_overlap():
    """TEST.724: five slots of three databases make the fifteen a stock server has besides the machine's own."""
    assert private_redis.SLOTS == 5 and private_redis.SLOT_WIDTH == 3
    every = [database for slot in range(private_redis.SLOTS) for database in private_redis.slot_databases(slot)]
    assert every == list(range(1, 16))
    for slot in (-1, 5):
        with pytest.raises(ValueError):
            private_redis.slot_databases(slot)


def test_a_run_with_a_slot_takes_the_databases_of_that_slot_only(monkeypatch):
    server = "redis://127.0.0.1:6379/0"
    monkeypatch.setenv(private_redis.SLOT_ENVIRONMENT, "2")      # databases 7, 8, 9
    taken = {url_for_worker(server, f"gw{n}") for n in range(6)}
    assert taken == {f"redis://127.0.0.1:6379/{database}" for database in (7, 8, 9)}
    monkeypatch.setenv(private_redis.SLOT_ENVIRONMENT, "0")
    assert {url_for_worker(server, f"gw{n}") for n in range(3)} == {f"redis://127.0.0.1:6379/{database}" for database in (1, 2, 3)}


def test_two_runs_with_different_slots_never_meet(monkeypatch):
    server = "redis://127.0.0.1:6379/0"
    seen = []
    for slot in ("0", "1", "4"):
        monkeypatch.setenv(private_redis.SLOT_ENVIRONMENT, slot)
        seen.append({url_for_worker(server, f"gw{n}") for n in range(3)})
    assert not (seen[0] & seen[1]) and not (seen[0] & seen[2]) and not (seen[1] & seen[2])


def test_the_slot_comes_from_the_environment_and_a_bad_value_is_none():
    assert private_redis.slot_of({private_redis.SLOT_ENVIRONMENT: "3"}) == 3
    assert private_redis.slot_of({}) is None and private_redis.slot_of({private_redis.SLOT_ENVIRONMENT: "x"}) is None
