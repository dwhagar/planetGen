# tests/test_container_plan.py

"""
DB.85: `planetgen plan` stores the galaxy's containers (its globular clusters and its superbubbles) beside the
components, drawn from the same seed, so planning again writes the same rows.
"""

import pytest

from planetgen import tuning
from planetgen.db import containers as containersDb, store
from tests.test_galaxy_gen import _build_real_skeleton

pytestmark = pytest.mark.fresh_database


def _stored(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        return containersDb.listing(conn)
    finally:
        conn.close()


def test_a_plan_stores_a_container_for_each_globular_cluster(mysql_config):
    _build_real_skeleton(mysql_config)
    found = _stored(mysql_config)
    clusters = [c for c in found if c["kind"] == "globular_cluster"]
    assert len(clusters) == tuning.CLUSTER_COUNT_MILKY_WAY
    assert all(c["name"] and c["component_id"] and c["extent_pc"] > 0 for c in clusters)
    assert all(c["kind"] != "superbubble" or c["bubble"] for c in found)
    assert len({c["uid"] for c in found}) == len(found)


def test_planning_again_stores_the_same_containers(mysql_config):
    _build_real_skeleton(mysql_config)
    first = _stored(mysql_config)
    _build_real_skeleton(mysql_config)
    again = _stored(mysql_config)
    strip = lambda rows: [{k: v for k, v in row.items() if k not in ("id", "component_id")} for row in rows]
    assert strip(first) == strip(again)


def test_a_cluster_removed_from_the_components_takes_its_container_with_it(mysql_config):
    _build_real_skeleton(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        before = len(containersDb.listing(conn, "globular_cluster"))
        conn.execute("DELETE FROM galaxy_components WHERE id = (SELECT component_id FROM galaxy_containers"
                     " WHERE kind = 'globular_cluster' ORDER BY id LIMIT 1)")
        conn.commit()
        assert len(containersDb.listing(conn, "globular_cluster")) == before - 1
    finally:
        conn.close()
