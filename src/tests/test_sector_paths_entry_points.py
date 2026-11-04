# tests/test_sector_paths_entry_points.py

"""
GEN.126: every place that changes the masses in a sector saves the sector
paths of that sector and its neighbours as its last step: the admin API's
queued jobs (delete a sector or a phenomenon, regenerate a phenomenon or a
sector, change a system's star) and `planetgen settle --sector-id` (the step after `planetgen phenomenon --sector-id`).
Each test marks every saved path (duration -1), runs the entry point and
looks for the marks it should have replaced. Needs MySQL.
"""

import sys

import pytest

from planetgen.api import edits as api_edits
from planetgen.cli import generate as generate_cli
from planetgen.db import sector_paths, store
from planetgen.generation import run_galaxy
from planetgen.queue import api_jobs
from tests import sector_support
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
    monkeypatch.setattr(run_galaxy, "regenerate_sector_cell",
                        lambda address, config: {"sector_id": None, "sector_name": None})
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


def test_the_phenomenon_command_only_generates_and_the_settle_command_saves_the_paths(mysql_config, monkeypatch):
    """PERF.124: `phenomenon` saves the phenomenon and nothing else; `settle --sector-id` is the step that follows it."""
    a_id, b_id = _two_sectors(mysql_config, monkeypatch)
    run_cli("phenomenon", ["--type", "black-hole", "--sector-id", str(a_id), "--quiet"] + mysql_argv(mysql_config))
    assert _marked(mysql_config, a_id) and _marked(mysql_config, b_id), "the phenomenon command settled nothing"
    run_cli("settle", ["--sector-id", str(a_id), "--quiet"] + mysql_argv(mysql_config))
    assert _marked(mysql_config, a_id) == 0 and _marked(mysql_config, b_id) == 0


def test_the_phenomenon_command_has_no_settle_option_and_settle_names_sectors_or_a_time_not_both(capsys):
    parser, commands = generate_cli.build_parser()
    assert "--no-settle" not in commands["phenomenon"].format_help()
    with pytest.raises(SystemExit):
        parser.parse_args(["phenomenon", "--type", "black-hole", "--no-settle"])
    assert parser.parse_args(["settle", "--sector-id", "3", "4"]).sector_id == [3, 4]
    with pytest.raises(SystemExit):
        run_cli("settle", ["--sector-id", "3", "--since", "2026-10-25 01:00:00"])
    assert "cannot be combined" in capsys.readouterr().err


def test_settle_refuses_a_sector_that_does_not_exist(mysql_config, monkeypatch, capsys):
    a_id, _b_id = _two_sectors(mysql_config, monkeypatch)
    with pytest.raises(SystemExit):
        run_cli("settle", ["--sector-id", str(a_id), "999999", "--quiet"] + mysql_argv(mysql_config))
    assert _marked(mysql_config, a_id), "nothing was settled for a request that named an unknown sector"


def test_settling_named_sectors_saves_the_paths_of_them_and_their_neighbours(mysql_config, monkeypatch):
    a_id, b_id = _two_sectors(mysql_config, monkeypatch)
    assert run_galaxy.settle_named_sectors([a_id], mysql_config) == 2
    assert _marked(mysql_config, a_id) == 0 and _marked(mysql_config, b_id) == 0

    conn = store.get_connection(mysql_config)
    try:
        conn.execute("UPDATE sector_paths SET duration_years = ?", (MARK,))
        conn.commit()
    finally:
        conn.close()

    sector_support.settle_sectors(mysql_config, [a_id])
    assert _marked(mysql_config, a_id) == 0 and _marked(mysql_config, b_id) == 0
