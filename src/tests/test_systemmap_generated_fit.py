# tests/test_systemmap_generated_fit.py

"""
MAP.88 over real generated systems: every element the System Map draws
(stars and their names, planets with their rings, moons, belts, orbits,
facilities and labels) sits inside its scene's viewBox, for single stars,
close pairs and wide pairs straight from generation, saved and read back
the way the system page reads them (`query.system_detail`).

Takes the `mysql_config` fixture (see `conftest.py`): skipped, not failed,
without a MySQL test server, like every other database-backed test.
"""

import os
import random
import sys

import pytest

from planetgen.db import query
from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "html", "lib"))

import systemmap as sm  # noqa: E402
from tests.test_systemmap import drawn_outside_view  # noqa: E402


def _config(kind):
    cfg = SystemConfig()
    if kind:
        cfg.BINARY_SYSTEM = True
        cfg.WIDE_BINARY = kind == "wide"
    return cfg


@pytest.mark.parametrize("kind", [None, "close", "wide"])
def test_generated_systems_fit_the_system_map(mysql_config, kind):
    rng_state = random.getstate()
    random.seed(8800 + (0 if kind is None else 1 if kind == "close" else 2))
    try:
        systems = [StarSystem(system_config=_config(kind)) for _ in range(8)]
    finally:
        random.setstate(rng_state)
    conn = store.get_connection(mysql_config)
    try:
        for system in systems:
            with conn:
                system_id = store.insert_star_system(conn, system, _config(kind))
            detail = query.system_detail(conn, system_id)
            html = sm.render_system_map_panel(detail, detail["stars"], detail["planets"], detail["belts"])
            assert drawn_outside_view(html) == [], (kind, system_id)
    finally:
        conn.close()
