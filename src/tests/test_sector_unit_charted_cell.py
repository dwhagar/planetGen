"""
DB.131: a sector unit that finds its cell charted does not fill it again. A unit the job engine runs again after its first try saved the sector (a connection dropped at the
commit, a worker killed just after it) would insert the sector's systems a second time, and the second `INSERT` stops on the object IDs the first gave them
(`star_systems.uq_star_systems_uid`: the sector, serial 0, body 0). `run_galaxy.generate_sectors_do` reports what the cell holds instead (`store.saved_sector_report`).
"""

import argparse
import math

import pytest

from planetgen import tuning
from planetgen.db import store
from planetgen.galaxy import spherical
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation import run_common, run_galaxy
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem

pytestmark = pytest.mark.db

ADDRESS = (41, 3, 5)


def _sector(count=2):
    sector = SpaceSector("Charted twice", edge_ly=20.0)
    for n in range(count):
        config = SystemConfig()
        config.STAR_TYPE = "G2V"
        config.PLANETS = False
        config.BINARY_SYSTEM = False
        sector.add_system(StarSystem(system_config=config), position=(1.0 + n, 1.0, 1.0), system_config=config)
    return sector


def _position(address=ADDRESS):
    center = spherical.sector_position_pc(*address, float(tuning.DEFAULT_SECTOR_EDGE_PC))
    return {"center_x_pc": float(center[0]), "center_y_pc": float(center[1]), "center_z_pc": float(center[2]),
            "galactic_radius_pc": math.dist(center, (0.0, 0.0, 0.0)),
            "shell_index": address[0], "band_index": address[1], "wedge_index": address[2]}


def _save(config, sector=None):
    return store.save_sector(sector or _sector(), config=config, galaxy_position=_position())


def test_a_second_fill_of_a_charted_cell_finds_the_first_fill_and_adds_nothing(mysql_config):
    """The failure of Boss's 8.0.4317 run, made on purpose: since OPS.268 the second save finds the cell charted (under the cell's lock) and returns the first save's sector, instead of stopping on the object IDs
    the first fill gave its systems (`uq_star_systems_uid`) after a whole sector's work."""
    first = _save(mysql_config)
    second = _save(mysql_config)
    assert second == first
    conn = store.get_connection(mysql_config)
    try:
        assert conn.execute("SELECT COUNT(*) AS n FROM star_systems WHERE sector_id = ?", (first,)).fetchone()["n"] == 2
    finally:
        conn.close()


def test_the_report_of_a_cell_says_what_it_holds_and_an_uncharted_cell_has_none(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        assert store.saved_sector_report(conn, ADDRESS) is None
    finally:
        conn.close()
    sector_id = _save(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        report = store.saved_sector_report(conn, ADDRESS)
        assert report["sector_id"] == sector_id and report["name"] and report["systems"] == 2
        assert report["stars"] == 2 and report["planets"] == 0 and report["phenomena"] == 0
        assert store.saved_sector_report(conn, (41, 3, 6)) is None
    finally:
        conn.close()


class _Context:
    def kept(self, key, make):
        return make()


class _Density:
    def __init__(self, config):
        pass

    def resolve(self, args, address, position_pc):
        return args


def _unit(config):
    return {"database": {"host": config.host, "port": config.port, "user": config.user, "database": config.database},
            "args": {"density": 1.0, "num_systems": None}, "edge_pc": float(tuning.DEFAULT_SECTOR_EDGE_PC)}


@pytest.fixture
def unit_args(mysql_config, monkeypatch):
    unit = _unit(mysql_config)
    monkeypatch.setattr(run_galaxy, "_BatchDensity", _Density)
    monkeypatch.setattr(run_galaxy, "_sector_args_of", lambda unit_args: argparse.Namespace(
        density=1.0, num_systems=None, mysql_host=mysql_config.host, mysql_port=mysql_config.port, mysql_user=mysql_config.user,
        mysql_password=mysql_config.password, mysql_database=mysql_config.database))
    return unit


def test_a_unit_that_finds_its_cell_charted_reports_it_and_does_not_fill_it_again(mysql_config, unit_args, monkeypatch):
    sector_id = _save(mysql_config)

    def fill(*args, **kwargs):
        raise AssertionError("a charted cell was filled again")

    monkeypatch.setattr(run_galaxy, "generate_and_save_sector_at", fill)
    result = run_galaxy.generate_sectors_do(unit_args, {"address": list(ADDRESS), "suffix": " (again)", "weight": 1}, _Context())
    assert result["address"] == list(ADDRESS) and result["sector_id"] == sector_id and result["systems"] == 2 and result["stars"] == 2
    assert result["stored"] is True and result["summary"] == [] and result["suffix"] == " (again)" and result["phenomena"] == 0
    assert {"name", "planets", "density", "memory_mb", "seconds"} <= set(result)


def test_a_unit_whose_cell_is_not_charted_fills_it(unit_args, monkeypatch):
    seen = []

    def fill(args, address, position_pc, edge_pc):
        seen.append(tuple(address))
        sector = _sector(1)
        sector.name = "Made"
        return 7, "Made", sector

    monkeypatch.setattr(run_galaxy, "generate_and_save_sector_at", fill)
    result = run_galaxy.generate_sectors_do(unit_args, {"address": list(ADDRESS), "suffix": "", "weight": 1}, _Context())
    assert seen == [ADDRESS] and result["sector_id"] == 7 and not result.get("stored")


def test_a_sector_found_charted_is_counted_for_the_run_but_adds_no_speed(monkeypatch):
    recorded = []
    monkeypatch.setattr(run_common, "_record_sector", lambda args, saved, seconds: recorded.append(saved["address"]))
    before = dict(run_common.RUN_COUNTS)
    sectors, systems = before.get("sectors", 0), before.get("systems", 0)
    try:
        run_common._count_reported(argparse.Namespace(), {"address": [1, 2, 3], "systems": 2, "phenomena": 1, "seconds": 0.0, "stored": True}, True)
        run_common._count_reported(argparse.Namespace(), {"address": [1, 2, 4], "systems": 3, "phenomena": 0, "seconds": 9.0}, True)
        assert recorded == [[1, 2, 4]]
        assert run_common.RUN_COUNTS["sectors"] == sectors + 2 and run_common.RUN_COUNTS["systems"] == systems + 5
    finally:
        run_common.RUN_COUNTS.clear()
        run_common.RUN_COUNTS.update(before)
