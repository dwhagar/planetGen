# tests.test_warm_start

"""PERF.141: the warm-up starts by itself after a data change (`queue/warmstart.py`, `planetgen warm --on-change`)."""

import os
import types
import uuid

import pytest

from planetgen.admin import updateguard
from planetgen.cli import warm as warm_cli
from planetgen.db import dataversion
from planetgen.queue import redisclient, warmstart


@pytest.fixture
def database(redis_server, monkeypatch):
    """A database name of this test's own, with the automatic start on and the process it would start faked."""
    name = f"planetgen_warmstart_{uuid.uuid4().hex[:8]}"
    monkeypatch.setenv("PLANETGEN_WARM_ON_CHANGE", "1")
    monkeypatch.setattr(updateguard, "is_updating", lambda: False)
    started = []
    monkeypatch.setattr(warmstart, "spawn", lambda db: started.append(db) or os.getpid())
    client = redisclient.client(decode_responses=True)
    yield types.SimpleNamespace(name=name, started=started, client=client)
    client.delete(warmstart.wanted_key(name), warmstart.lock_key(name), dataversion.key(name))


def test_a_change_starts_one_warm_up_and_the_next_one_joins_it(database):
    assert warmstart.request(database.name, "sector-set") == "started"
    assert warmstart.request(database.name, "system-data") == "running"
    assert database.started == [database.name]
    assert warmstart.pending(database.name) == ["sector-set", "system-data"]
    assert database.client.get(warmstart.lock_key(database.name)) == str(os.getpid())


def test_a_holder_that_is_gone_is_taken_over(database):
    database.client.set(warmstart.lock_key(database.name), 2 ** 22 + 12345)
    assert warmstart.request(database.name, "orbits") == "started"
    assert database.started == [database.name]


def test_nothing_starts_when_it_is_off_or_an_update_runs(database, monkeypatch):
    monkeypatch.setenv("PLANETGEN_WARM_ON_CHANGE", "0")
    assert warmstart.request(database.name, "orbits") == "off"
    monkeypatch.setenv("PLANETGEN_WARM_ON_CHANGE", "1")
    monkeypatch.setattr(updateguard, "is_updating", lambda: True)
    assert warmstart.request(database.name, "orbits") == "update"
    assert warmstart.request(None, "orbits") == "off"
    assert database.started == [] and warmstart.pending(database.name) == []


def test_a_start_that_fails_lets_go_of_the_lock_and_says_so(database, monkeypatch, capsys):
    def refuse(_db):
        raise OSError("no such interpreter")

    monkeypatch.setattr(warmstart, "spawn", refuse)
    assert warmstart.request(database.name, "orbits") == "failed"
    assert "did not start" in capsys.readouterr().err
    assert database.client.get(warmstart.lock_key(database.name)) is None
    assert warmstart.pending(database.name) == ["orbits"]


def test_a_data_version_bump_asks_for_the_warm_up(database):
    assert dataversion.bump(database.name, dataversion.TAG_SECTOR_SET, dataversion.TAG_SECTOR_DATA) is True
    assert database.started == [database.name]
    assert warmstart.pending(database.name) == [dataversion.TAG_SECTOR_SET]
    dataversion.record(database.name, "orbits")
    assert database.started == [database.name], "the one already waiting takes it"
    assert dataversion.TAG_POSITIONS in warmstart.pending(database.name)


def test_the_loop_warms_once_for_a_burst_and_again_for_a_change_during_the_run(database):
    warmstart.request(database.name, "a")
    warmstart.request(database.name, "b")
    runs, sleeps = [], []

    def run(kinds):
        runs.append(kinds)
        if len(runs) == 1:
            warmstart.request(database.name, "c")

    assert warmstart.drain(database.name, run, settle=7, sleep=sleeps.append) == 2
    assert runs == [["a", "b"], ["c"]]
    assert sleeps and set(sleeps) == {7}
    assert database.client.get(warmstart.lock_key(database.name)) is None
    assert warmstart.pending(database.name) == []


def test_the_loop_stops_at_its_most_and_leaves_the_rest_for_the_next_request(database):
    warmstart.request(database.name, "a")

    def run(_kinds):
        warmstart.request(database.name, "again")

    assert warmstart.drain(database.name, run, settle=0, max_rounds=3, sleep=lambda _s: None) == 3
    assert warmstart.pending(database.name) == ["again"]
    assert database.client.get(warmstart.lock_key(database.name)) is None
    assert warmstart.request(database.name, "next") == "started"


def test_a_failed_run_lets_go_of_the_lock_and_raises(database):
    warmstart.request(database.name, "a")

    def run(_kinds):
        raise warm_cli.Failed("a family failed")

    with pytest.raises(warm_cli.Failed):
        warmstart.drain(database.name, run, settle=0, sleep=lambda _s: None)
    assert database.client.get(warmstart.lock_key(database.name)) is None


def test_the_process_runs_as_the_owner_of_the_cache_directory_when_root_starts_it(monkeypatch):
    seen = {}

    def popen(argv, **options):
        seen.update(argv=argv, options=options)
        return types.SimpleNamespace(pid=4242)

    monkeypatch.setattr(warmstart.subprocess, "Popen", popen)
    monkeypatch.setattr(warmstart, "_cache_owner", lambda: (33, 34))
    assert warmstart.spawn("galaxy") == 4242
    assert seen["argv"][-3:] == ["warm", "--on-change", "galaxy"]
    assert seen["options"]["user"] == 33 and seen["options"]["group"] == 34 and seen["options"]["cwd"] == "/"
    monkeypatch.setattr(warmstart, "_cache_owner", lambda: None)
    warmstart.spawn("galaxy")
    assert "user" not in seen["options"]


def _on_change_args(changed, site):
    parser, _parsers = __import__("planetgen.cli.generate", fromlist=["x"]).build_parser()
    args = parser.parse_args(["warm", "--on-change", changed])
    args.database = site
    return args


def test_the_warm_up_after_a_change_runs_the_whole_warm_up_for_the_site_database(database, monkeypatch):
    warmed = []
    monkeypatch.setattr(warm_cli, "_site_database", lambda: database.name)
    monkeypatch.setattr(warm_cli, "_warm_once", lambda args: warmed.append(args.database))
    monkeypatch.setattr(warmstart, "settle_seconds", lambda: 0.0)
    warmstart.request(database.name, "orbits")
    warm_cli.run_warm(_on_change_args(database.name, database.name))
    assert warmed == [database.name]
    assert database.client.get(warmstart.lock_key(database.name)) is None


def test_a_change_on_another_database_warms_nothing(database, monkeypatch):
    monkeypatch.setattr(warm_cli, "_site_database", lambda: "the_site")
    monkeypatch.setattr(warm_cli, "_warm_once", lambda args: pytest.fail("nothing to warm"))
    warmstart.request(database.name, "orbits")
    warm_cli.run_warm(_on_change_args(database.name, "the_site"))
    assert database.client.get(warmstart.lock_key(database.name)) is None


def test_a_failed_warm_up_after_a_change_exits_non_zero(database, monkeypatch, caplog):
    monkeypatch.setattr(warm_cli, "_site_database", lambda: database.name)
    monkeypatch.setattr(warmstart, "settle_seconds", lambda: 0.0)

    def fail(_args):
        raise warm_cli.Failed("systems did not answer")

    monkeypatch.setattr(warm_cli, "_warm_once", fail)
    warmstart.request(database.name, "orbits")
    with pytest.raises(SystemExit) as exit_info:
        warm_cli.run_warm(_on_change_args(database.name, database.name))
    assert exit_info.value.code == 1
    assert "the cache warm-up after a change stopped: systems did not answer" in caplog.text
