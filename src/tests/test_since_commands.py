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
def test_the_time_is_optional_and_must_be_a_time(command, capsys):
    assert _args(command).since is None        # the newest sectors or galaxy run's start (PERF.123)
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


def _path_rows(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        return [tuple(row.values()) for row in conn.execute(
            "SELECT p.sector_id, p.object_table, p.object_id, p.exited, p.duration_years, k.position, k.t_years, k.x_pc,"
            " k.y_pc, k.z_pc, k.vx_kms, k.vy_kms, k.vz_kms FROM sector_paths p JOIN sector_path_knots k ON k.path_id = p.id"
            " ORDER BY p.sector_id, p.object_table, p.object_id, k.position").fetchall()]
    finally:
        conn.close()


def test_the_settle_stage_is_cut_into_static_chunks_that_save_the_same_paths(mysql_config, monkeypatch, redis_server):
    from planetgen.generation import run_galaxy

    ids = []
    for wedge in (700, 701, 702, 703):
        monkeypatch.setattr(helper, "ADDRESS", (1500, 2, wedge))
        ids.append(helper._save(mysql_config, systems=2)[0])
    params = {"database": wrapper.database_args(mysql_config), "since": LONG_AGO, "workers": 4}
    whole = run_galaxy.expand_save_sector_paths(params)
    assert len(whole.units) == 1 and whole.units[0][0]["workers"] == 4, "a short stage is one unit"
    _run(monkeypatch, mysql_config, "settle", "--since", LONG_AGO)
    once = _path_rows(mysql_config)
    assert once

    monkeypatch.setattr(run_galaxy, "SETTLE_UNIT_SECTORS", 1)
    cut = run_galaxy.expand_save_sector_paths(params)
    assert len(cut.units) == cut.plan["sectors"] >= len(ids) and all(args["workers"] == 1 for args, _cost in cut.units)
    assert sorted(sector for args, _cost in cut.units for sector in args["sector_ids"]) == sorted(
        {sector for args, _cost in whole.units for sector in args["sector_ids"]})
    assert [args["sector_ids"] for args, _cost in run_galaxy.expand_save_sector_paths(params).units] == [
        args["sector_ids"] for args, _cost in cut.units], "the same chunks every time"
    conn = store.get_connection(mysql_config)
    try:
        conn.execute("DELETE FROM sector_paths")
        conn.commit()
    finally:
        conn.close()
    _run(monkeypatch, mysql_config, "settle", "--since", LONG_AGO)
    assert _path_rows(mysql_config) == once, "the paths are those of the one-unit stage"


def test_a_settle_that_names_no_sector_is_one_empty_unit(monkeypatch):
    from planetgen.generation import run_galaxy

    expansion = run_galaxy.expand_save_sector_paths({"database": {}, "sector_ids": [], "expand": False, "workers": 3})
    assert [args["sector_ids"] for args, _cost in expansion.units] == [[]] and expansion.plan == {"sectors": 0}
    monkeypatch.setattr(run_galaxy, "SETTLE_UNIT_SECTORS", 2)
    expansion = run_galaxy.expand_save_sector_paths(
        {"database": {}, "sector_ids": [9, 3, 5, 1, 7, 5], "expand": False, "workers": 3})
    assert [args["sector_ids"] for args, _cost in expansion.units] == [[1, 3], [5, 7], [9]]
    assert [cost for _args, cost in expansion.units] == [2.0, 2.0, 1.0]


def _link_rows(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        return [tuple(row.values()) for row in conn.execute(
            "SELECT sector_id, object_table, object_id, neighbor_rank, neighbor_system_id, distance_pc FROM nearest_systems"
            " ORDER BY sector_id, object_table, object_id, neighbor_rank").fetchall()]
    finally:
        conn.close()


def test_the_link_stage_is_cut_into_static_chunks_that_store_the_same_links(mysql_config, monkeypatch, redis_server):
    from planetgen.generation import run_galaxy

    for wedge in (700, 701, 702, 703):
        monkeypatch.setattr(helper, "ADDRESS", (1500, 2, wedge))
        helper._save(mysql_config, systems=2)
    saved = _link_rows(mysql_config)
    assert saved, "saving a sector links it"
    params = {"database": wrapper.database_args(mysql_config), "since": LONG_AGO, "workers": 4}
    assert len(run_galaxy.expand_link_neighbours(params).units) == 1, "a short stage is one unit"

    monkeypatch.setattr(run_galaxy, "LINK_UNIT_SECTORS", 1)
    cut = run_galaxy.expand_link_neighbours(params)
    assert len(cut.units) == cut.plan["sectors"] == 4
    assert [unit["chunk"] for unit, _cost in cut.units] == [unit["chunk"] for unit, _cost in
                                                           run_galaxy.expand_link_neighbours(params).units]
    conn = store.get_connection(mysql_config)
    try:
        conn.execute("DELETE FROM nearest_systems")
        conn.commit()
    finally:
        conn.close()
    _run(monkeypatch, mysql_config, "link", "--since", LONG_AGO)
    assert _link_rows(mysql_config) == saved, "the links are those of the one-unit stage"
