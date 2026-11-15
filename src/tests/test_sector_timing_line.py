# tests/test_sector_timing_line.py

"""PERF.270: a stage of sectors says each sector's cost on a line of its own, so a long stage says which sectors took the time (the stage 5 of 21,919 s of PERF.267 had no such line)."""

import argparse
import os

from planetgen.generation import run_galaxy, run_sector

ADDRESS = (211, -53, 1022)


def _saved(**more):
    saved = {"address": list(ADDRESS), "sector_id": 7, "name": "Velnalo Dranhas", "systems": 37, "stars": 54, "planets": 452, "density": 3.5, "phenomena": 16, "summary": "",
             "memory_mb": 2181.4, "seconds": 853.0, "suffix": ""}
    saved.update(more)
    return saved


def _follower(monkeypatch, count):
    """A follower of `count` sectors whose plan has been announced and whose lines are caught."""
    lines = []
    monkeypatch.setattr(run_galaxy.log, "normal", lines.append)
    monkeypatch.setattr(run_galaxy.run_common, "_count_reported", lambda *_args: None)
    follower = run_galaxy.SectorsFollower(argparse.Namespace(), "neighborhood", "Neighborhood", 4.0, announce=lambda *_args: None)
    follower.on_plan({"detail": {"sectors": count}})
    return follower, lines


def test_the_line_says_the_sector_its_systems_planets_time_and_memory():
    line = run_galaxy.timing_line(_saved(), ADDRESS)
    assert "37 systems, 452 planets" in line
    assert "in 14.2 minutes," in line      # the time in the project's duration form
    assert "2,181 MB" in line


def test_a_stage_of_few_sectors_says_each_sector_and_names_it_in_the_plan_line(monkeypatch):
    follower, lines = _follower(monkeypatch, 109)
    assert follower.timed and "each sector says its time" in lines[0]
    follower.on_result("unit", {"pid": os.getpid(), "sectors": [_saved()]})
    assert len(lines) == 2 and "37 systems" in lines[1]


def test_a_stage_of_very_many_sectors_keeps_the_line_in_the_debug_log(monkeypatch):
    follower, lines = _follower(monkeypatch, run_galaxy.SECTOR_TIMING_LINES_MAX + 1)
    assert not follower.timed and "debug log" in lines[0]
    follower.on_result("unit", {"pid": os.getpid(), "sectors": [_saved()]})
    assert len(lines) == 1


def test_a_result_without_planets_or_memory_still_says_its_line():
    line = run_galaxy.timing_line({"systems": 3, "seconds": 5.0}, ADDRESS)
    assert "3 systems, 0 planets" in line and "0 MB" in line


def test_the_planets_of_a_sector_do_not_count_its_belts():
    belt, planet = argparse.Namespace(body_type="a"), argparse.Namespace(body_type="p")
    entry = argparse.Namespace(star_system=argparse.Namespace(planets=[planet, belt, planet]))
    assert run_sector.sector_planet_count(argparse.Namespace(entries=[entry, entry])) == 4
