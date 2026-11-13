"""OPS.157: the stage worker caps and the adaptive starting limit read the job profile (workers.stage_worker_cap, adaptive.start_limit)."""

import pytest

from planetgen.queue import adaptive, jobprofiles, memorygate, profileclaim, redisqueue, runrecord, workers
from planetgen.queue.memorygate import Host

KIND = "link-neighbors"
TOTAL_MB = 16000.0


@pytest.fixture
def conn(redis_server, key_prefix, monkeypatch):
    monkeypatch.setattr(memorygate, "PREFIX", f"{key_prefix}:gate")
    monkeypatch.setattr(runrecord, "PREFIX", key_prefix)
    monkeypatch.setenv(profileclaim.ENV_VAR, "1")
    monkeypatch.setitem(workers.STAGE_WORKER_MEMORY_MB, "generate-sectors", 1000)         # eight workers in half of 16 GB
    monkeypatch.setattr(memorygate, "read_host", lambda: Host(TOTAL_MB, 8000.0, 8, 0.1))
    monkeypatch.setattr(workers, "_physical_memory_bytes", lambda: int(TOTAL_MB * 2 ** 20))
    from planetgen.util import generator_epoch
    monkeypatch.setattr(generator_epoch, "generator_epoch", lambda: 146)
    return redisqueue.connect(redis_server)


def history(peak, count=30):
    return lambda kind, **k: [{"peak_mb": peak, "outcome": "done", "host_memory_mb": TOTAL_MB, "host_cores": 8, "generator_epoch": 146}] * count


def test_with_no_history_the_cap_is_the_typed_one(conn, monkeypatch):
    monkeypatch.setattr(jobprofiles, "peaks", history(0.0, 0))
    assert workers.stage_worker_cap(KIND, conn, "r1") == workers.stage_worker_cap(KIND)
    assert workers.stage_worker_mb(KIND, conn, "r1") == workers.stage_worker_mb(KIND)


def test_a_kind_that_is_one_unit_by_nature_has_no_cap_whatever_the_profile_says(conn, monkeypatch):
    monkeypatch.setattr(jobprofiles, "peaks", history(100.0))
    assert workers.stage_worker_cap("not-a-staged-command", conn, "r1") is None


def test_history_that_shows_smaller_workers_lets_more_of_them_run_within_a_quarter(conn, monkeypatch):
    typed = workers.stage_worker_mb("generate-sectors")
    monkeypatch.setattr(jobprofiles, "peaks", history(typed / 4))
    cap = workers.stage_worker_cap("generate-sectors", conn, "r1")
    assert cap > workers.stage_worker_cap("generate-sectors")
    assert workers.stage_worker_mb("generate-sectors", conn, "r1") == pytest.approx(typed * 0.75)       # one run moves a figure by at most a quarter
    assert workers.stage_worker_cap("generate-sectors", conn, "r1") == cap                              # held for the run


def test_history_that_shows_bigger_workers_lowers_the_cap(conn, monkeypatch):
    typed = workers.stage_worker_mb("generate-sectors")
    monkeypatch.setattr(jobprofiles, "peaks", history(typed * 1.2))
    assert workers.stage_worker_cap("generate-sectors", conn, "r1") < workers.stage_worker_cap("generate-sectors")


def test_off_gives_the_typed_cap(conn, monkeypatch):
    monkeypatch.setenv(profileclaim.ENV_VAR, "0")
    monkeypatch.setattr(jobprofiles, "peaks", history(10.0))
    assert workers.stage_worker_cap("generate-sectors", conn, "r1") == workers.stage_worker_cap("generate-sectors")


def test_a_memory_wait_that_lowers_the_limit_is_remembered_for_the_kind(conn, monkeypatch):
    monkeypatch.setattr(adaptive, "lower", lambda *a, **k: 3)
    assert adaptive.on_memory_wait(conn, "w", f"unit:{KIND}", 31.0) == 3
    assert adaptive.settled(conn, f"unit:{KIND}") == {"limit": 3, "memory": TOTAL_MB, "cores": 8}


def test_the_next_run_of_the_kind_starts_at_the_settled_limit_and_logs_it(conn, caplog):
    import logging
    caplog.set_level(logging.INFO, logger="planetgen")
    adaptive._remember(conn, f"unit:{KIND}", 3)
    assert adaptive.start_limit(conn, KIND, 8, now=1000.0) == 3
    assert adaptive.limit_of(conn) == 3 and adaptive.recent(conn)[0]["from"] == 8
    assert f"start limit of unit:{KIND} was 8.0, now 3.0" in caplog.text
    (change,) = profileclaim.changes(conn)
    assert change["what"] == "start limit" and change["new"] == 3


def test_a_limit_already_set_is_left_alone_and_so_is_a_cap_below_it(conn):
    adaptive._remember(conn, f"unit:{KIND}", 6)
    assert adaptive.start_limit(conn, KIND, 6) is None                      # not under the cap
    adaptive._set(conn, 2, 1000.0, "lowered_at")
    adaptive._remember(conn, f"unit:{KIND}", 4)
    assert adaptive.start_limit(conn, KIND, 8) is None and adaptive.limit_of(conn) == 2


def test_a_limit_from_another_machine_shape_counts_half(conn):
    conn.hset(adaptive._key("settled"), f"unit:{KIND}", '{"limit": 2, "memory": 64000, "cores": 32}')
    assert adaptive.start_limit(conn, KIND, 8) == 5                          # (2 + 8) / 2


def test_no_memory_of_a_kind_starts_nothing_and_off_forgets_nothing(conn, monkeypatch):
    assert adaptive.start_limit(conn, KIND, 8) is None
    monkeypatch.setenv(profileclaim.ENV_VAR, "0")
    adaptive._remember(conn, f"unit:{KIND}", 3)
    assert adaptive.settled(conn, f"unit:{KIND}") is None


def test_a_limit_that_has_risen_to_the_cap_is_forgotten(conn):
    adaptive._remember(conn, f"unit:{KIND}", 3)
    adaptive._remember(conn, f"unit:{KIND}", None)
    assert adaptive.settled(conn, f"unit:{KIND}") is None
