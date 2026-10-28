# tests/scored_site.py

"""
UX.214, UX.232: a small database served over HTTP for the tests that follow search results: one sector holding three
generated systems (fixed seeds, so the run is repeatable and fast), whose planets and moons are then given scores that
cover every PHI-6 colour, every domain colour and every life score, so each filter and each sort has bodies to pick and
order. `scored_galaxy` is the fixture; it holds the database config and the server's address.
"""

import uuid
from dataclasses import dataclass

import pymysql
import pytest

from planetgen.db import objectindex, store
from planetgen.db.store import MySQLConfig
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics.habitability_world import DOMAINS
from planetgen.util import draw
from planetgen.web.app import create_app

from tests import dbpool
from tests.conftest import _test_server_kwargs
from tests.habitability_site import EDGE_MPC, HABITABLE_ADDRESS, _config_for
from tests.served import served

SEEDS = (11, 12, 13)
"""tuple: The seed of each of the three systems."""

PHI4_CYCLE = (1.2, 0.9, 0.5, 0.1)
"""tuple: The overall PHI-6 score given to the bodies in turn: Blue, Green, Yellow, Red."""


@dataclass(frozen=True)
class ScoredGalaxy:
    config: object
    base_url: str


def _score(conn, table):
    """Gives every planet and moon of `table` scores that cycle through the colours and the life scores, so every filter has
    a match and the sorts have different values to order (the scores a generated body has are replaced)."""
    ids = [row["id"] for row in conn.execute(f"SELECT id FROM {table} ORDER BY id").fetchall()]
    for index, body_id in enumerate(ids):
        sets = {"phi4": PHI4_CYCLE[index % 4]}
        for step, domain in enumerate(DOMAINS):
            sets[f"tier_{domain}"] = (index + step) % 4
        sets.update(life_microbial=index % 10, life_flora=(index * 3) % 10, life_fauna=(index * 2) % 5,
                    life_sentient=1 if index % 7 == 0 else 0)
        marks = ", ".join(f"{column} = ?" for column in sets)
        conn.execute(f"UPDATE {table} SET {marks} WHERE id = ?", (*sets.values(), body_id))


@pytest.fixture(scope="module")
def scored_galaxy(_mysql_server_available):
    kwargs = _test_server_kwargs()
    name = f"planetgen_scored_{uuid.uuid4().hex[:12]}"
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
            conn.execute("INSERT INTO sectors (name, edge_mpc, shell_index, band_index, wedge_index) VALUES (?, ?, ?, ?, ?)",
                         ("Scored", EDGE_MPC, *HABITABLE_ADDRESS))
            sector_id = conn.execute("SELECT id FROM sectors WHERE name = 'Scored'").fetchone()["id"]
            system_ids = []
            for seed in SEEDS:
                cfg = SystemConfig()
                cfg.STAR_TYPE = "G2V"
                cfg.PLANETS = True
                cfg.MAX_PLANETS = True
                cfg.BINARY_SYSTEM = False
                with draw.bound(seed):
                    system = StarSystem(system_config=cfg)
                system_ids.append(store.insert_star_system(conn, system, cfg, sector_id=sector_id))
            store.assign_uids(conn, system_ids=system_ids)
            for table in ("planets", "moons"):
                _score(conn, table)
            objectindex.build(conn)
            conn.commit()
        finally:
            conn.close()
        app = create_app(_config_for(config))
        app.testing = True
        from planetgen import tuning
        kept = tuning.SLOW_PAGE_SECONDS
        tuning.SLOW_PAGE_SECONDS = float("inf")
        try:
            with served(app) as url:
                yield ScoredGalaxy(config=config, base_url=url)
        finally:
            tuning.SLOW_PAGE_SECONDS = kept
    finally:
        store.close_pool(config)
        dbpool.drop_database(kwargs, name)
