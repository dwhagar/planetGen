"""PERF.123 (slice A): `planetgen link`, `backfill` and `settle` each do their one stage for the galaxy sectors saved since a
time, as a run of one command through the wrapper."""

import datetime

import pytest

from planetgen.cli import generate
from planetgen.db import sector_paths, store
from planetgen.queue import commands, redisqueue, runrecord, wrapper
from tests import test_sector_paths_db as helper
from tests.bughunt_support import mysql_argv

pytestmark = pytest.mark.db

LONG_AGO = "2000-01-01 00:00:00"
FAR_AHEAD = "2999-01-01 00:00:00"


def _args(*argv):
    return generate.build_parser()[0].parse_args(list(argv))


def _run(monkeypatch, mysql_config, *argv):
    monkeypatch.setattr("sys.argv", ["generate", *argv, *mysql_argv(mysql_config)])
    generate.main()


@pytest.mark.parametrize("command, kind", [("link", "link-neighbours"), ("backfill", "scatter-backfill"),
                                           ("settle", "settle-sector-paths")])
def test_each_command_is_one_stage_and_one_registered_command(command, kind):
    args = _args(command, "--since", LONG_AGO)
    assert args.since == datetime.datetime(2000, 1, 1)
    assert [step.kind for step in wrapper.steps_for(command, args)] == [kind]
    assert kind in commands.REGISTRY


@pytest.mark.parametrize("command", ["link", "backfill", "settle"])
def test_the_time_is_required_and_must_be_a_time(command, capsys):
    with pytest.raises(SystemExit):
        _args(command)
    with pytest.raises(SystemExit):
        _args(command, "--since", "yesterday")
    assert "not a time" in capsys.readouterr().err


def test_settle_saves_the_paths_of_the_sectors_since_the_time_and_leaves_older_ones(mysql_config, monkeypatch, redis_server):
    a_id, _center = helper._save(mysql_config, systems=1)
    conn = store.get_connection(mysql_config)
    try:
        sector_paths.settle_sectors(conn, [a_id])
        conn.execute("UPDATE sector_paths SET duration_years = -1")
        conn.commit()
    finally:
        conn.close()
    _run(monkeypatch, mysql_config, "settle", "--since", FAR_AHEAD)
    conn = store.get_connection(mysql_config)
    try:
        marked = conn.execute("SELECT COUNT(*) AS n FROM sector_paths WHERE duration_years = -1").fetchone()["n"]
    finally:
        conn.close()
    assert marked > 0, "a time after every sector settles nothing"
    _run(monkeypatch, mysql_config, "settle", "--since", LONG_AGO)
    conn = store.get_connection(mysql_config)
    try:
        assert conn.execute("SELECT COUNT(*) AS n FROM sector_paths WHERE duration_years = -1").fetchone()["n"] == 0
    finally:
        conn.close()
    run = runrecord.read_run(redisqueue.connect(redis_server), runrecord.recent_runs(redisqueue.connect(redis_server), 5)[0])
    assert run["state"] == "done" and run["title"] == "Sector paths"


def test_link_says_how_many_sectors_it_linked(mysql_config, monkeypatch, redis_server, capsys):
    helper._save(mysql_config, systems=1)
    _run(monkeypatch, mysql_config, "link", "--since", LONG_AGO)
    assert "Linked 1 sectors" in capsys.readouterr().out


@pytest.mark.parametrize("command", ["link", "backfill", "settle"])
def test_a_time_after_every_sector_starts_no_run(mysql_config, monkeypatch, redis_server, capsys, command):
    helper._save(mysql_config, systems=1)
    before = runrecord.recent_runs(redisqueue.connect(redis_server), 50)
    _run(monkeypatch, mysql_config, command, "--since", FAR_AHEAD)
    assert "nothing to" in capsys.readouterr().out
    assert runrecord.recent_runs(redisqueue.connect(redis_server), 50) == before
