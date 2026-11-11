"""PERF.228: the habitability re-score's units, run by four workers at once, no longer die with a MySQL deadlock (error 1213)."""
import threading

import pymysql
import pytest

from planetgen.db import store
from planetgen.generation import run_plan
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.queue import wrapper
from planetgen.util import draw
from tests.fuzz_support import deterministic_entropy

WORKERS = 4
SYSTEMS_PER_WORKER = 4
ROUNDS = 12


def _save(mysql_config, count):
    ids = []
    for seed in range(count):
        with deterministic_entropy(seed):
            draw.set_run_seed(seed)
            cfg = SystemConfig()
            system = StarSystem(cfg, galactic_center_dist_ly=26000.0)
        ids.append(store.save_system(system, cfg, config=mysql_config))
    conn = store.get_connection(mysql_config)
    with conn:      # the interest rows (MAP.305) only exist for a system with a galaxy position
        conn.execute("UPDATE star_systems SET galaxy_x_pc = id * 10, galaxy_y_pc = 1, galaxy_z_pc = 1")
    conn.close()
    return ids


def _unit(mysql_config, low, high):
    """`all` because GEN.1479 skips a system the running code already scored, and the contention these tests make needs every round to re-score."""
    return run_plan.rescore_habitability_unit({"database": wrapper.database_args(mysql_config), "low": low, "high": high, "all": True}, None)


def test_four_workers_re_scoring_neighbouring_ranges_all_finish(mysql_config):
    ids = _save(mysql_config, WORKERS * SYSTEMS_PER_WORKER)
    failures = []

    def work(first):
        for _ in range(ROUNDS):
            try:
                _unit(mysql_config, ids[first], ids[first + SYSTEMS_PER_WORKER - 1])
            except pymysql.err.OperationalError as exc:
                failures.append(exc.args[0])

    threads = [threading.Thread(target=work, args=(worker * SYSTEMS_PER_WORKER,)) for worker in range(WORKERS)]
    [thread.start() for thread in threads]
    [thread.join() for thread in threads]
    assert failures == []


def test_a_unit_that_meets_a_deadlock_starts_over_and_gives_the_same_answer(mysql_config, monkeypatch):
    ids = _save(mysql_config, 2)
    expected = _unit(mysql_config, ids[0], ids[-1])
    real = store.refresh_habitability
    hits = []

    def deadlocked_once(conn, system_id):
        if not hits:
            hits.append(system_id)
            real(conn, system_id)           # work already done in the transaction the deadlock takes back
            raise pymysql.err.OperationalError(1213, "Deadlock found when trying to get lock; try restarting transaction")
        return real(conn, system_id)

    monkeypatch.setattr(store, "refresh_habitability", deadlocked_once)
    monkeypatch.setattr(run_plan.time, "sleep", lambda seconds: None)
    assert _unit(mysql_config, ids[0], ids[-1]) == expected
    assert hits


def test_a_unit_that_keeps_meeting_a_deadlock_gives_up_after_the_attempt_limit(mysql_config, monkeypatch):
    ids = _save(mysql_config, 1)
    calls = []

    def always(conn, system_id):
        calls.append(system_id)
        raise pymysql.err.OperationalError(1213, "Deadlock found when trying to get lock; try restarting transaction")

    monkeypatch.setattr(store, "refresh_habitability", always)
    monkeypatch.setattr(run_plan.time, "sleep", lambda seconds: None)
    with pytest.raises(pymysql.err.OperationalError):
        _unit(mysql_config, ids[0], ids[0])
    assert len(calls) == store.SECTOR_SAVE_ATTEMPTS


def test_any_other_database_error_is_not_retried(mysql_config, monkeypatch):
    ids = _save(mysql_config, 1)
    calls = []

    def broken(conn, system_id):
        calls.append(system_id)
        raise pymysql.err.OperationalError(1146, "Table doesn't exist")

    monkeypatch.setattr(store, "refresh_habitability", broken)
    with pytest.raises(pymysql.err.OperationalError):
        _unit(mysql_config, ids[0], ids[0])
    assert len(calls) == 1
