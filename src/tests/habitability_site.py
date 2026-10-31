# tests/habitability_site.py

"""
The shared fixture of the habitability system's end-to-end browser tests (TEST.300): a throwaway database with
a sector that holds a seeded system with Blue or Green worlds, a second sector whose system has none, the
Galaxy Map's stats and block totals recorded for both, and a local server over that database. Each lane
tests its own piece against fixture data; the tests that use this check the pieces together on stored scores.

`HabitableGalaxy` names what the tests look for: the ids and names, the number of Blue or Green bodies the stored
scores give each system, and the address of each sector.
"""

import uuid
from dataclasses import dataclass

import pymysql
import pytest

from planetgen.api.config import Config
from planetgen.api.limiter import API_LIMIT_OFF, PAGE_LIMITS_OFF
from planetgen.db import block_stats, interests, objectindex, store
from planetgen.db.store import MySQLConfig
from planetgen.galaxy import spherical
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import habitability
from planetgen.util import draw
from planetgen.web.app import create_app
from tests import dbpool
from tests.conftest import _test_server_kwargs
from tests.served import served

HABITABLE_ADDRESS = (1, 0, 0)
BARREN_ADDRESS = (1, 0, 1)
EDGE_PC = 4.0
"""float: The sector edge of the galaxy, parsecs (the default of the galaxy model)."""
EDGE_MPC = 4000
HABITABLE_CENTER_PC = spherical.sector_position_pc(*HABITABLE_ADDRESS, EDGE_PC)
BARREN_CENTER_PC = spherical.sector_position_pc(*BARREN_ADDRESS, EDGE_PC)
SYSTEM_OFFSET_PC = (0.5, 0.5, 0.5)
"""tuple: Where a sector's one system sits from its center, parsecs."""


@dataclass(frozen=True)
class HabitableGalaxy:
    """What the seeded database holds."""
    habitable_sector_id: int
    habitable_system_id: int
    habitable_system_name: str
    habitable_worlds: int      # the Blue and Green planets and moons: the sector's stored count
    marker_worlds: int         # the Yellow and better ones the map's habitable marker counts for the system (Boss, 7:31 PM EDT 2026-10-24)
    barren_sector_id: int
    barren_system_id: int
    barren_system_name: str
    habitable_bodies: tuple     # (kind, id, name) of every Blue or Green planet and moon of the habitable system
    base_url: str
    config: object = None   # the MySQLConfig of the database


def _config_for(site_config):
    class SiteConfig(Config):
        MYSQL_CONFIG = site_config
        WRITE_MYSQL_CONFIG = site_config
        CONTROL_MYSQL_CONFIG = site_config
        WEB_DATABASE = ""
        SESSION_COOKIE_SECURE = False
        SECRET_KEY = "habitability-site-secret"
        RATELIMIT_PAGES = PAGE_LIMITS_OFF
        RATELIMIT_DEFAULT = API_LIMIT_OFF   # fixtures and pages call the API from one address (TEST.385)

    return SiteConfig


def _system(seed, want_habitable):
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.PLANETS = True
    cfg.MAX_PLANETS = True
    cfg.BINARY_SYSTEM = False
    cfg.HABITABLE_WORLD = want_habitable
    with draw.bound(seed):
        return StarSystem(system_config=cfg), cfg


def _tiers(conn, system_id):
    """`[(kind, id, name, tier)]` of the planets and moons of one system with a stored overall score."""
    rows = []
    for kind, table in (("planet", "planets"), ("moon", "moons")):
        for row in conn.execute(f"SELECT id, name, phi4 FROM {table} WHERE star_system_id = ? AND phi4 IS NOT NULL", (system_id,)).fetchall():
            rows.append((kind, row["id"], row["name"], habitability.tier(row["phi4"])))
    return rows


def _add_sector(conn, name, address, center_pc, want_habitable):
    """A sector at `address` with one system that has Blue or Green worlds (`want_habitable`) or none, its stats recorded;
    returns (sector id, system id, system name, the Blue or Green bodies)."""
    conn.execute("INSERT INTO sectors (name, edge_mpc, shell_index, band_index, wedge_index) VALUES (?, ?, ?, ?, ?)",
                 (name, EDGE_MPC, *address))
    sector_id = conn.execute("SELECT id FROM sectors WHERE name = ?", (name,)).fetchone()["id"]
    for seed in range(400):
        system, cfg = _system(seed, want_habitable)
        system_id = store.insert_star_system(conn, system, cfg, sector_id=sector_id)
        good = [row for row in _tiers(conn, system_id) if row[3] in habitability.POTENTIALLY_HABITABLE_TIERS]
        if bool(good) == want_habitable and _tiers(conn, system_id):
            break
        conn.execute("DELETE FROM star_systems WHERE id = ?", (system_id,))
    else:
        pytest.fail(f"no system {'with' if want_habitable else 'without'} a Blue or Green world in 400 seeds")
    conn.execute("UPDATE star_systems SET galaxy_x_pc = ?, galaxy_y_pc = ?, galaxy_z_pc = ? WHERE id = ?",
                 (*(c + o for c, o in zip(center_pc, SYSTEM_OFFSET_PC)), system_id))
    store.assign_uids(conn, system_ids=[system_id])
    store.recount_habitable_worlds(conn, [system_id])
    store.record_sector_stats(conn, sector_id, address, center_pc)
    block_stats.refresh_chain(conn, address)
    name = conn.execute("SELECT name FROM star_systems WHERE id = ?", (system_id,)).fetchone()["name"]
    return sector_id, system_id, name, tuple((kind, body_id, body) for kind, body_id, body, _tier in good)


def _renamed(conn, added):
    """`added` (from `_add_sector`) with the names the database holds now."""
    sector_id, system_id, _name, bodies = added
    name = conn.execute("SELECT name FROM star_systems WHERE id = ?", (system_id,)).fetchone()["name"]
    names = tuple((kind, body_id, conn.execute(f"SELECT name FROM {'planets' if kind == 'planet' else 'moons'} WHERE id = ?",
                                               (body_id,)).fetchone()["name"]) for kind, body_id, _old in bodies)
    return sector_id, system_id, name, names


@pytest.fixture(scope="module")
def habitable_galaxy(_mysql_server_available):
    """The seeded database, served: a `HabitableGalaxy`."""
    kwargs = _test_server_kwargs()
    name = f"planetgen_habsite_{uuid.uuid4().hex[:12]}"
    admin = pymysql.connect(**kwargs)
    try:
        with admin.cursor() as cur:
            cur.execute(f"CREATE DATABASE `{name}`")
        admin.commit()
    finally:
        admin.close()
    config = MySQLConfig(database=name, **kwargs)
    try:
        conn = store.get_connection(config)
        try:
            habitable = _add_sector(conn, "Habit", HABITABLE_ADDRESS, HABITABLE_CENTER_PC, True)
            barren = _add_sector(conn, "Barren", BARREN_ADDRESS, BARREN_CENTER_PC, False)
            objectindex.build(conn)
            conn.commit()
            habitable, barren = (_renamed(conn, sector) for sector in (habitable, barren))     # two systems of one seed get an Alpha and a Beta
            marker_worlds = sum(conn.execute(f"SELECT COUNT(*) AS n FROM {table} WHERE star_system_id = ? AND phi4 >= ?",
                                             (habitable[1], interests.HABITABLE_MIN_SCORE)).fetchone()["n"] for table in ("planets", "moons"))
        finally:
            conn.close()
        app = create_app(_config_for(config))
        app.testing = True
        from planetgen import tuning
        kept = tuning.SLOW_PAGE_SECONDS
        tuning.SLOW_PAGE_SECONDS = float("inf")      # a loaded run must not turn the first build of a page into the loading state
        try:
            with served(app) as url:
                yield HabitableGalaxy(
                    habitable_sector_id=habitable[0], habitable_system_id=habitable[1], habitable_system_name=habitable[2],
                    habitable_worlds=len(habitable[3]), marker_worlds=marker_worlds, barren_sector_id=barren[0], barren_system_id=barren[1],
                    barren_system_name=barren[2], habitable_bodies=habitable[3], base_url=url, config=config)
        finally:
            tuning.SLOW_PAGE_SECONDS = kept
    finally:
        store.close_pool(config)
        dbpool.drop_database(kwargs, name)


__all__ = ["habitable_galaxy", "HabitableGalaxy", "HABITABLE_ADDRESS", "BARREN_ADDRESS", "HABITABLE_CENTER_PC", "BARREN_CENTER_PC", "EDGE_PC"]
