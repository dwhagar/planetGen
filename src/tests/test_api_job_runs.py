"""PERF.119 (F3): the API's jobs that generate sectors (neighborhood, regenerate, settle) start runs through the wrapper: each
stage is a run on the run record."""

import itertools

import pytest

from planetgen.generation import run_galaxy, run_sector
from planetgen.queue import api_jobs, commands, redisqueue, runrecord
from tests import first_fill
from tests.test_galaxy_gen import _all_sectors, _seed_skeleton, _shells
from tests.sector_support import make_sector

pytestmark = [pytest.mark.db, pytest.mark.usefixtures("unit_density_plans")]


@pytest.fixture(autouse=True)
def few_systems(monkeypatch):
    """Each sector these jobs fill holds one system, and the first fill does not draw down to the backfill's deepest floor
    (GEN.633): a sector of the toy galaxy holds hundreds of systems, so a test that makes a few sectors ran past its time,
    and what is tested is the runs a job starts, not what a sector holds (TEST.429)."""
    first_fill.without_the_deepest_floor(monkeypatch)
    monkeypatch.setattr(run_sector, "_sample_poisson_count", lambda mean, rng=None: 1)


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
    first = make_sector(0, 0, 0, config=mysql_config, backfill=False, settle=False)
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


def test_a_block_job_runs_the_sectors_of_the_block_through_the_wrapper(mysql_config, conn):
    """GEN.310: the sectors the outline holds in a block, as runs; a second job finds nothing left to make."""
    from planetgen.galaxy import drill
    _seed_skeleton(mysql_config, shells=[(0, 0)])           # the outline holds the core and nothing else: one sector to make
    block = drill.drill_chain_of(0, 0, 0)[-2]
    key = drill.format_drill_key(block)
    estimate = run_galaxy.generate_block(key, config=mysql_config, estimate_only=True)
    assert estimate["generated"] == 0 and estimate["candidates"] == 1 and estimate["already_existed"] == 0
    assert estimate["outside_galaxy"] >= 1, "the block holds sectors past the outline too"
    assert _all_sectors(mysql_config) == [], "an estimate writes nothing"
    result = api_jobs.generate_block(key, mysql_config)
    titles, runs = _titles(conn)
    saw = (result, len(_all_sectors(mysql_config)), titles, [(run["title"], run["state"]) for run in runs])
    assert result["generated"] == estimate["candidates"] == len(_all_sectors(mysql_config)), saw
    assert titles[0] == f"Sectors (block {key})" and titles[-2:] == ["Bright-star backfill", "Sector paths"], saw
    assert {run["state"] for run in runs} == {"done"}, saw
    again = run_galaxy.generate_block(key, config=mysql_config, estimate_only=True)
    assert again["already_existed"] == estimate["candidates"] and again["candidates"] == estimate["candidates"]


def test_a_block_that_holds_more_than_a_request_may_list_is_refused(mysql_config, monkeypatch):
    from planetgen.galaxy import drill
    _seed_skeleton(mysql_config, shells=_shells(2, 1))
    monkeypatch.setattr(run_galaxy.run_common, "WEB_ENUMERATION_LIMIT", 2)
    with pytest.raises(run_galaxy.EnumerationTooLarge):
        run_galaxy.generate_block(drill.format_drill_key(drill.drill_chain_of(2, 0, 0)[-2]), config=mysql_config,
                                  estimate_only=True)


def test_settling_the_named_sectors_is_a_run_of_them(mysql_config, conn):
    center = _a_galaxy_sector(mysql_config)
    assert isinstance(run_galaxy.settle_named_sectors([center], mysql_config), int)
    titles, runs = _titles(conn)
    assert titles == ["Sector paths"] and runs[0]["state"] == "done"


def test_the_named_sectors_are_the_ones_the_settle_stage_works_on():
    assert run_galaxy._sectors_since({"sector_ids": [4, 2], "since": "2000-01-01 00:00:00"}) == [4, 2]
    assert commands.REGISTRY["settle-sector-paths"].unit and commands.REGISTRY["scatter-backfill"].expand
