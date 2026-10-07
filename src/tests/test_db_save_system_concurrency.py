# tests/test_db_save_system_concurrency.py

"""
`save_system` with several one-off systems saved at once (found by
TEST.37): systems with the same name saved concurrently are named
"Alpha X", "Beta X", ... like a sector's (the first used to keep its bare
name), and a deadlock or lock wait timeout is retried as `save_sector`
retries it.
"""

import threading

import pymysql
import pytest

from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation.system import StarSystem

WRITERS = 4


def _system(name):
    cfg = SystemConfig()
    cfg.STAR_TYPE = "M2V"
    cfg.BINARY_SYSTEM = False
    system = StarSystem(system_config=cfg)
    system.name = name
    return system, cfg


def _save_at_once(mysql_config, name):
    systems = [_system(name) for _ in range(WRITERS)]
    start = threading.Barrier(WRITERS)
    errors, ids = [], []

    def save(system, cfg):
        start.wait()
        try:
            ids.append(store.save_system(system, cfg, config=mysql_config))
        except Exception as exc:  # reported below, with the others
            errors.append(exc)

    threads = [threading.Thread(target=save, args=pair) for pair in systems]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    assert not errors, errors
    conn = store.get_connection(mysql_config)
    try:
        return sorted(row["name"] for row in conn.execute(
            f"SELECT name FROM star_systems WHERE id IN ({', '.join('?' * len(ids))})", tuple(ids)).fetchall())
    finally:
        conn.close()


@pytest.mark.parametrize("round_", range(3))
def test_same_named_systems_saved_at_once_are_all_lettered(mysql_config, round_):
    names = _save_at_once(mysql_config, "Halveth")
    assert names == ["Alpha Halveth", "Beta Halveth", "Delta Halveth", "Gamma Halveth"]


def test_same_named_systems_saved_at_once_after_a_matching_sector(mysql_config):
    store.save_sector(SpaceSector("Halveth", edge_ly=11.5), config=mysql_config)
    names = _save_at_once(mysql_config, "Halveth")
    assert len(set(names)) == WRITERS
    # One takes a diminutive ("Little Halveth"); the rest would need a
    # Greek letter too, three words, so they draw fresh names (GEN.46).
    assert sum(name.endswith(" Halveth") for name in names) == 1
    assert "Halveth" not in names
    assert all(len(name.split()) <= 2 for name in names)


def test_save_system_retries_a_deadlock_with_its_name_restored(mysql_config, monkeypatch):
    system, cfg = _system("Retrivel")
    real_insert = store.insert_star_system
    seen = []

    def deadlock_once(conn, star_system, system_config, *args, **kwargs):
        seen.append(star_system.name)
        if len(seen) == 1:
            star_system.name = "Changed Before The Deadlock"
            raise pymysql.err.OperationalError(1213, "Deadlock found when trying to get lock")
        return real_insert(conn, star_system, system_config, *args, **kwargs)

    monkeypatch.setattr(store, "insert_star_system", deadlock_once)
    monkeypatch.setattr(store.time, "sleep", lambda _seconds: None)
    system_id = store.save_system(system, cfg, config=mysql_config)
    assert seen == ["Retrivel", "Retrivel"]
    conn = store.get_connection(mysql_config)
    try:
        assert conn.execute("SELECT name FROM star_systems WHERE id = ?", (system_id,)).fetchone()["name"] == "Retrivel"
    finally:
        conn.close()


def test_save_system_gives_up_after_the_last_attempt(mysql_config, monkeypatch):
    system, cfg = _system("Stuckvel")
    calls = []

    def always_deadlocks(*_args, **_kwargs):
        calls.append(1)
        raise pymysql.err.OperationalError(1213, "Deadlock found when trying to get lock")

    monkeypatch.setattr(store, "insert_star_system", always_deadlocks)
    monkeypatch.setattr(store.time, "sleep", lambda _seconds: None)
    with pytest.raises(pymysql.err.OperationalError):
        store.save_system(system, cfg, config=mysql_config)
    assert len(calls) == store.SECTOR_SAVE_ATTEMPTS
