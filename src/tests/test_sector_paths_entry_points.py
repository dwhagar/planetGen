# tests/test_sector_paths_entry_points.py

"""
GEN.126: every place that changes the masses in a sector saves the sector
paths of that sector and its neighbours as its last step: the admin API's
queued jobs (delete a sector or a phenomenon, regenerate a phenomenon or a
sector, change a system's star) and `planetgen phenomenon --sector-id`.
Each test marks every saved path (duration -1), runs the entry point and
looks for the marks it should have replaced. Needs MySQL.
"""

import sys

import pytest

from planetgen.api import edits as api_edits
from planetgen.db import sector_paths, store
from planetgen.generation import run_galaxy
from planetgen.queue import api_jobs
from tests import test_sector_paths_db as helper
from tests.bughunt_support import mysql_argv, run_cli

pytestmark = pytest.mark.db

MARK = -1.0


def _two_sectors(mysql_config, monkeypatch):
    """Sector A (a black hole and a system) and its neighbour B (a system), paths saved and marked."""
    here = sys.modules[helper.__name__]
    a_id, _center = helper._save(mysql_config, systems=1)
    monkeypatch.setattr(here, "ADDRESS", (1500, 2, 701))
    b_id, _center = helper._save(mysql_config, with_black_hole=False, systems=1)
    conn = store.get_connection(mysql_config)
    try:
        sector_paths.settle_sectors(conn, [a_id])
        conn.execute("UPDATE sector_paths SET duration_years = ?", (MARK,))
        conn.commit()
    finally:
        conn.close()
    return a_id, b_id


def _marked(mysql_config, sector_id):
    """How many of the sector's paths still carry the mark."""
    conn = store.get_connection(mysql_config)
    try:
        return conn.execute("SELECT COUNT(*) AS n FROM sector_paths WHERE sector_id = ? AND duration_years = ?",
                            (sector_id, MARK)).fetchone()["n"]
    finally:
        conn.close()


def _count(mysql_config, sector_id):
    conn = store.get_connection(mysql_config)
    try:
        return conn.execute("SELECT COUNT(*) AS n FROM sector_paths WHERE sector_id = ?", (sector_id,)).fetchone()["n"]
    finally:
        conn.close()


def _black_hole_id(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        return conn.execute("SELECT MIN(id) AS id FROM black_holes").fetchone()["id"]
    finally:
        conn.close()


def test_marks_are_set_up(mysql_config, monkeypatch):
    a_id, b_id = _two_sectors(mysql_config, monkeypatch)
    assert _marked(mysql_config, a_id) and _marked(mysql_config, b_id)


def test_deleting_a_phenomenon_settles_its_sector_and_the_neighbours(mysql_config, monkeypatch):
    a_id, b_id = _two_sectors(mysql_config, monkeypatch)
    result = api_edits.delete_phenomenon_job(mysql_config, "black_hole", _black_hole_id(mysql_config))
    assert result["status"] == "ok"
    assert _marked(mysql_config, a_id) == 0 and _marked(mysql_config, b_id) == 0
    assert _count(mysql_config, a_id) and _count(mysql_config, b_id)


def test_regenerating_a_phenomenon_settles_its_sector_and_the_neighbours(mysql_config, monkeypatch):
    a_id, b_id = _two_sectors(mysql_config, monkeypatch)
    api_edits.regenerate_phenomenon_job(mysql_config, "black_hole", _black_hole_id(mysql_config))
    assert _marked(mysql_config, a_id) == 0 and _marked(mysql_config, b_id) == 0


def test_deleting_a_sector_settles_the_neighbours(mysql_config, monkeypatch):
    a_id, b_id = _two_sectors(mysql_config, monkeypatch)
    result = api_edits.delete_sector_job(mysql_config, a_id)
    assert result["systems"] == 1 and "Sector deleted" in result["summary"]
    assert _marked(mysql_config, b_id) == 0 and _count(mysql_config, b_id) == 1
    assert _count(mysql_config, a_id) == 0


def test_regenerating_a_sector_settles_the_neighbours(mysql_config, monkeypatch):
    a_id, b_id = _two_sectors(mysql_config, monkeypatch)
    monkeypatch.setattr(run_galaxy, "ensure_sector_generated",
                        lambda *address, config=None, settle=True: {"sector_id": None, "sector_name": None})
    result = api_jobs.regenerate_sector(a_id, mysql_config)
    assert result["deleted"]["systems"] == 1
    assert _marked(mysql_config, b_id) == 0


def test_changing_a_systems_star_settles_its_sector_and_the_neighbours(mysql_config, monkeypatch):
    a_id, b_id = _two_sectors(mysql_config, monkeypatch)
    conn = store.get_connection(mysql_config)
    try:
        system_id = conn.execute("SELECT id FROM star_systems WHERE sector_id = ?", (a_id,)).fetchone()["id"]
    finally:
        conn.close()
    api_edits.change_star_job(mysql_config, system_id, "K2V", True)
    assert _marked(mysql_config, a_id) == 0 and _marked(mysql_config, b_id) == 0


def test_the_phenomenon_command_settles_its_sector_unless_told_not_to(mysql_config, monkeypatch):
    a_id, b_id = _two_sectors(mysql_config, monkeypatch)
    run_cli("phenomenon", ["--type", "black-hole", "--sector-id", str(a_id), "--no-settle", "--quiet"]
            + mysql_argv(mysql_config))
    assert _marked(mysql_config, b_id) == 1
    run_cli("phenomenon", ["--type", "black-hole", "--sector-id", str(a_id), "--quiet"] + mysql_argv(mysql_config))
    assert _marked(mysql_config, a_id) == 0 and _marked(mysql_config, b_id) == 0


def test_a_sector_generated_on_the_spot_queues_a_settle_job(mysql_config, monkeypatch):
    a_id, b_id = _two_sectors(mysql_config, monkeypatch)
    submitted = []
    monkeypatch.setattr(api_jobs, "submit", lambda function, *args: submitted.append((function, args)) or "job1")
    assert run_galaxy.queue_settle(mysql_config, [b_id, a_id]) == "job1"
    assert submitted == [(api_jobs.settle_sectors, ([a_id, b_id], mysql_config))]
    assert _marked(mysql_config, a_id) and _marked(mysql_config, b_id)  # nothing ran yet


def test_the_settle_job_saves_the_paths_and_without_a_queue_it_runs_inline(mysql_config, monkeypatch):
    a_id, b_id = _two_sectors(mysql_config, monkeypatch)
    assert api_jobs.settle_sectors([a_id], mysql_config) == {"saved": 2}
    assert _marked(mysql_config, a_id) == 0 and _marked(mysql_config, b_id) == 0

    conn = store.get_connection(mysql_config)
    try:
        conn.execute("UPDATE sector_paths SET duration_years = ?", (MARK,))
        conn.commit()
    finally:
        conn.close()

    def no_queue(function, *args):
        raise api_jobs.NoQueue("no Redis")

    monkeypatch.setattr(api_jobs, "submit", no_queue)
    assert run_galaxy.queue_settle(mysql_config, [a_id]) is None
    assert _marked(mysql_config, a_id) == 0 and _marked(mysql_config, b_id) == 0
