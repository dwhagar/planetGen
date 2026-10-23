# tests/test_cluster_plan.py

"""
GEN.160: `planetgen plan` draws the galaxy's globular clusters from its seed before it walks the outline, stores them in
`galaxy_components`, and the shape read back carries them (`GalaxyShape.components`) for the sector gate.
"""

import pytest

from planetgen import tuning
from planetgen.db import store
from planetgen.galaxy import components, mass_chain
from planetgen.generation import clusters
from tests.test_galaxy_gen import _build_real_skeleton

pytestmark = pytest.mark.fresh_database


def _rows(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        return store.get_galaxy_components(conn, kind="cluster"), store.get_galaxy_shape(conn)
    finally:
        conn.close()


def planned_cluster_count(mysql_config, info):
    """How many globular clusters the stored plan's seed and halo draw (GEN.403: the count follows the halo mass)."""
    conn = store.get_connection(mysql_config)
    try:
        seed = store.get_galaxy_seed(conn)
    finally:
        conn.close()
    chain = mass_chain.draw_mass_chain(seed, info.shape.structure.stellar_mass_solar)
    return clusters.cluster_count(seed, chain.halo_mass_solar)


def _satellite_rows(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        return store.get_galaxy_components(conn, kind="satellite")
    finally:
        conn.close()


def test_a_plan_stores_the_clusters_and_the_shape_carries_them(mysql_config):
    _build_real_skeleton(mysql_config)
    found, info = _rows(mysql_config)
    globular = [c for c in found if not components.is_open_cluster(c)]
    assert len(globular) == planned_cluster_count(mysql_config, info)
    assert len(found) > len(globular)              # the young open clusters (GEN.244) are stored with them
    assert {c.source for c in found} == {"generated"}
    assert info.shape.components is not None and len(info.shape.components) == len(found) + len(_satellite_rows(mysql_config))
    assert info.shape.components.edge_pc == info.edge_pc
    assert info.shape.components.expected_at_density_1 == info.expected_system_count_at_density_1


def test_the_clusters_scale_with_the_galaxy_so_a_toy_keeps_them_inside_it(mysql_config):
    _build_real_skeleton(mysql_config)         # a toy galaxy of 40 pc scale length
    found, _info = _rows(mysql_config)
    scale = 40.0 / tuning.HALO_SCALE_LENGTH_PC
    farthest = max(sum(v * v for v in c.center_pc) ** 0.5 for c in found if not components.is_open_cluster(c))
    assert farthest <= tuning.CLUSTER_HALO_RADIUS_RANGE_PC[1] * scale * 1.001


def test_planning_again_keeps_the_seed_and_so_the_same_clusters(mysql_config):
    _build_real_skeleton(mysql_config)
    first, _ = _rows(mysql_config)
    _build_real_skeleton(mysql_config)
    again, info = _rows(mysql_config)
    assert [(c.name, c.center_pc, c.mass_solar, c.king_w0) for c in first] == [
        (c.name, c.center_pc, c.mass_solar, c.king_w0) for c in again]
    assert info.shape.components is not None


def test_a_galaxy_without_components_has_none_on_its_shape(mysql_config):
    _build_real_skeleton(mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        conn.execute("DELETE FROM galaxy_components")
        conn.commit()
        assert store.get_galaxy_shape(conn).shape.components is None
    finally:
        conn.close()
