# tests/test_name_claims.py

"""
PERF.49: system names are claimed in a short transaction of their own, so a
writer saving a dense sector no longer holds the name registry's row locks
until its whole save commits and the other workers no longer queue for them.

Every database test takes the `mysql_config` fixture and is skipped, not
failed, when no MySQL test server is reachable.
"""

import threading

import pymysql
import pytest

from planetgen.db import store
from planetgen.galaxy.sector import SpaceSector

from tests.test_parallel_names import _names, _save_in_sector


@pytest.fixture
def mysql_config(mysql_config):
    store.get_connection(mysql_config).close()
    return mysql_config


def _registry(config, base):
    conn = store.get_connection(config)
    try:
        return conn.execute(
            "SELECT occurrence_count, first_star_system_id FROM system_name_registry WHERE base_name = ?",
            (base,)).fetchone()
    finally:
        conn.close()


def test_a_save_does_not_wait_for_another_that_claimed_the_same_name(mysql_config, monkeypatch):
    """A claims "Vega" and pauses before it commits; B claims it, saves and finishes meanwhile (B would have
    waited for A's commit when the claim was part of A's transaction). A then names itself "Alpha Vega"."""
    real_confirm = store.confirm_system_names
    at_confirm, release = threading.Event(), threading.Event()
    first = threading.current_thread()

    def confirm(conn, confirmations):
        if threading.current_thread() is not first:
            at_confirm.set()
            assert release.wait(timeout=60)
        return real_confirm(conn, confirmations)

    monkeypatch.setattr(store, "confirm_system_names", confirm)
    errors = []

    def slow_save():
        try:
            _save_in_sector(mysql_config, 0, "Vega")
        except BaseException as exc:  # noqa: BLE001 -- re-raised below
            errors.append(exc)

    writer = threading.Thread(target=slow_save)
    writer.start()
    try:
        assert at_confirm.wait(timeout=60)
        _save_in_sector(mysql_config, 1, "Vega")  # returns while the first save is still open
        assert _names(mysql_config, "star_systems") == ["Beta Vega"]
    finally:
        release.set()
        writer.join(timeout=60)
    if errors:
        raise errors[0]
    assert sorted(_names(mysql_config, "star_systems")) == ["Alpha Vega", "Beta Vega"]
    registry = _registry(mysql_config, "Vega")
    assert registry["occurrence_count"] == 2 and registry["first_star_system_id"] is not None
    conn = store.get_connection(mysql_config)
    try:
        holder = conn.execute("SELECT name FROM star_systems WHERE id = ?", (registry["first_star_system_id"],)).fetchone()
    finally:
        conn.close()
    assert holder["name"] == "Alpha Vega"


def test_a_failed_save_gives_its_claims_back(mysql_config, monkeypatch):
    sector = SpaceSector(name="Claimed")
    from tests.test_parallel_names import _system
    system = _system("Quillon")
    sector.add_system(system, position=(0.0, 0.0, 0.0), system_config=system.system_config)

    def broken(conn, confirmations):
        raise pymysql.err.ProgrammingError(1064, "boom")

    monkeypatch.setattr(store, "confirm_system_names", broken)
    with pytest.raises(pymysql.err.ProgrammingError):
        store.save_sector(sector, config=mysql_config)
    monkeypatch.undo()
    assert _registry(mysql_config, "Quillon") is None
    assert system.name == "Quillon"
    store.save_sector(sector, config=mysql_config)
    assert _names(mysql_config, "star_systems") == ["Quillon"]
    assert _registry(mysql_config, "Quillon")["occurrence_count"] == 1


def test_a_claim_someone_else_followed_stays_counted(mysql_config):
    """A count is only given back while nobody has claimed after it, so a later holder keeps its place."""
    conn = store.get_connection(mysql_config)
    claim = store.get_connection(mysql_config)
    try:
        conn.claim_connection = claim
        reserved = store.reserve_system_names(conn, ["Zorrel"])
        assert reserved == [("Zorrel", "Zorrel", None)]
        other = store.get_connection(mysql_config)
        try:
            other.claim_connection = store.get_connection(mysql_config)
            assert store.reserve_system_names(other, ["Zorrel"]) == [("Beta Zorrel", "Zorrel", None)]
        finally:
            other.claim_connection.close()
            other.close()
        store._release_name_claims(claim, conn.name_claims)
    finally:
        claim.close()
        conn.close()
    assert _registry(mysql_config, "Zorrel")["occurrence_count"] == 2
