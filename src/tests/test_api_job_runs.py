"""PERF.119 (F3): the API's jobs that generate sectors (neighborhood, regenerate, settle) start runs through the wrapper: each
stage is a run on the run record."""

import itertools

import pytest

from planetgen.db import store
from planetgen.generation import run_galaxy
from planetgen.queue import api_jobs, commands, redisqueue, runrecord
from tests.test_galaxy_gen import _all_sectors, _seed_skeleton, _shells

pytestmark = [pytest.mark.db, pytest.mark.usefixtures("unit_density_plans")]


@pytest.fixture
def conn(redis_server, monkeypatch):
    """The Redis connection; run ids count up, because ids made in the same second sort by their random tail."""
    numbers = itertools.count()
    monkeypatch.setattr(runrecord, "new_run_id", lambda now=None: f"20261022-000000-{next(numbers):08x}")
    return redisqueue.connect(redis_server)


def _titles(conn):
    """The titles of the runs on the record, oldest first."""
    runs = [runrecord.read_run(conn, run_id) for run_id in runrecord.recent_runs(conn, 50)]
    return [run["title"] for run in reversed(runs)], runs


def _a_galaxy_sector(mysql_config):
    _seed_skeleton(mysql_config, shells=_shells(2, 1))
    first = run_galaxy.ensure_sector_generated(0, 0, 0, config=mysql_config, backfill=False, settle=False)
    assert first["created"]
    return first["sector_id"]


def test_regenerating_a_sector_is_a_run_of_the_sector_a_backfill_and_a_settle(mysql_config, conn):
    old_id = _a_galaxy_sector(mysql_config)
    result = api_jobs.regenerate_sector(old_id, mysql_config)
    assert result["sector_id"] is not None and result["sector_id"] != old_id
    assert [row["id"] for row in _all_sectors(mysql_config)] == [result["sector_id"]]
    titles, runs = _titles(conn)
    assert len(titles) == 3 and titles[0].startswith("Sector (") and titles[1:] == ["Bright-star backfill", "Sector paths"]
    assert {run["state"] for run in runs} == {"done"} and {run["kind"] for run in runs} == {runrecord.GENERATION}


def test_a_neighborhood_job_runs_its_sectors_through_the_wrapper(mysql_config, conn):
    center = _a_galaxy_sector(mysql_config)
    result = api_jobs.generate_neighborhood(center, 25.0, mysql_config)
    titles, runs = _titles(conn)
    # TEST.289: this failed once under a gate's four workers with no text; the answer, the sectors and the runs go in the message.
    saw = (result, len(_all_sectors(mysql_config)), titles, [(run["title"], run["state"]) for run in runs])
    assert result["generated"] >= 1 and len(_all_sectors(mysql_config)) == 1 + result["generated"], saw
    assert titles[0] == f"Sectors (around sector {center})", saw
    assert titles[-2:] == ["Bright-star backfill", "Sector paths"], saw
    assert {run["state"] for run in runs} == {"done"}, saw


def test_the_settle_job_is_a_run_of_the_named_sectors(mysql_config, conn):
    center = _a_galaxy_sector(mysql_config)
    assert isinstance(api_jobs.settle_sectors([center], mysql_config)["saved"], int)
    titles, runs = _titles(conn)
    assert titles == ["Sector paths"] and runs[0]["state"] == "done"


def test_the_named_sectors_are_the_ones_the_settle_stage_works_on():
    assert run_galaxy._sectors_since({"sector_ids": [4, 2], "since": "2000-01-01 00:00:00"}) == [4, 2]
    assert commands.REGISTRY["settle-sector-paths"].unit and commands.REGISTRY["scatter-backfill"].expand


def test_a_settle_that_can_be_neither_queued_nor_run_says_so_and_does_not_fail_the_caller(monkeypatch):
    def no_queue(*args, **options):
        raise api_jobs.NoQueue("an update is running")

    def no_redis(*args, **options):
        raise redisqueue.Unavailable("no Redis server answers")

    monkeypatch.setattr(api_jobs, "submit", no_queue)
    monkeypatch.setattr(api_jobs, "settle_sectors", no_redis)
    assert run_galaxy.queue_settle(None, [3, 1]) is None
