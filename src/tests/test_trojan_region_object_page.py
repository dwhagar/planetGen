# tests/test_trojan_region_object_page.py

"""
NAV.80: a giant's trojan region (GEN.1150) takes no belt number, so the object resolver that read a number for every belt raised KeyError for it. It names the region for its
giant now, as every other list of belts does (`belt_names.names_for_belts`), and so does the chain of a facility built on one.
"""

import pytest

from planetgen.db import query, store
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.generation.trojan_region import TrojanRegion
from planetgen.physics import trojan_swarm
from planetgen.util import draw

pytestmark = pytest.mark.fresh_database


def _saved_system_with_regions(conn, config):
    cfg = SystemConfig()
    cfg.PLANETS = True
    cfg.MAX_PLANETS = True
    cfg.BINARY_SYSTEM = False
    cfg.FORCE_CLASS = "J"
    with draw.bound(3):
        made = StarSystem(system_config=cfg)
    giant = next(planet for planet in made.primary.orbiters if planet.body_type == "g")
    regions = [TrojanRegion(cfg, made.star, giant, region)
               for region in trojan_swarm.draw_regions(made.star.mass, giant.mass, 1.0e22, draw.Stream("test:trojan-page"))]
    assert len(regions) == 2
    made.primary.orbiters.extend(regions)
    with conn:
        system_id = store.insert_star_system(conn, made, cfg)
    return system_id, giant


def test_a_trojan_region_resolves_with_its_giants_name(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        system_id, giant = _saved_system_with_regions(conn, mysql_config)
        regions = conn.execute("SELECT id, lagrange_point FROM asteroid_belts WHERE star_system_id = ? AND lagrange_point IS NOT NULL ORDER BY lagrange_point", (system_id,)).fetchall()
        giant_row = conn.execute("SELECT name FROM planets WHERE id = (SELECT parent_id FROM asteroid_belts WHERE id = ?)", (regions[0]["id"],)).fetchone()
        for row in regions:
            found = query.resolve_object(conn, "belt", row["id"])
            assert found["name"] == f"{giant_row['name']} {row['lagrange_point']} trojan region"
    finally:
        conn.close()
