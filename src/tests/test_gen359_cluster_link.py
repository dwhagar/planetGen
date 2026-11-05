# tests/test_gen359_cluster_link.py

"""
GEN.359: a star of a cluster points to the cluster. The component pass writes the cluster's `galaxy_components` row on every star it places for a
cluster (none for a satellite galaxy's, none when the component has no row yet), and a system built around such a pre-placed star passes the
cluster to all of its stars: through `mark_bright_star_filled` for a single uncharted system and through the sector save for a sector's fill.
"""

import dataclasses

from planetgen.db import store
from tests.bright_star_rows import B2V, star_row
from tests.test_db_persistence import _named_system
from tests.test_gen_component_pass import SMALL_SATELLITE_MASS, _cluster, _satellite, _stars


def test_the_component_pass_writes_the_cluster_on_each_star_of_a_cluster_and_none_on_a_satellites():
    column = store.BRIGHT_STAR_COLUMNS.index("cluster_id")
    cluster = dataclasses.replace(_cluster(), id=42)
    rows = _stars([cluster, _satellite(mass=SMALL_SATELLITE_MASS)])
    by_population = {}
    for row in rows:
        by_population.setdefault(row[6], set()).add(row[column])
    assert by_population["cluster"] == {42}
    assert all(ids == {None} for population, ids in by_population.items() if population != "cluster")
    assert {row[column] for row in _stars([_cluster()])} == {None}          # a component with no row yet links to nothing


def test_the_link_does_not_move_the_draw():
    column = store.BRIGHT_STAR_COLUMNS.index("cluster_id")
    plain = _stars([_cluster()])
    linked = _stars([dataclasses.replace(_cluster(), id=7)])
    assert [row[:column] for row in plain] == [row[:column] for row in linked]


def test_a_system_built_around_a_clusters_star_passes_the_cluster_to_its_stars(mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            conn.execute("INSERT INTO galaxy_shape (id, disk_scale_length_pc, disk_scale_height_pc,"
                         " bulge_scale_radius_pc, bulge_amplitude, arm_count, pitch_angle_rad, arm_amplitude,"
                         " spiral_reference_radius_pc, spiral_reference_angle_rad, k_norm, edge_pc,"
                         " expected_system_count_at_density_1, outer_shell_index)"
                         " VALUES (1, 1, 1, 1, 1, 2, 0.2, 0.3, 1, 0, 1, 4, 10, 5)")
            store.replace_galaxy_components(conn, [_cluster("Linked")])
            cluster_id = store.get_galaxy_components(conn)[0].id
            store.insert_bright_stars(conn, [
                star_row((3, 0, 0), (12000, 3000, 0), luminosity_w=3.0e29, cluster_id=cluster_id, **B2V),
                star_row((3, 0, 1), (12000, 3000, 0), luminosity_w=3.0e29, **B2V)])
        linked, field = sorted(store.bright_stars_for_sector(conn, 3, 0, 0) + store.bright_stars_for_sector(conn, 3, 0, 1),
                               key=lambda row: row["wedge_index"])
        assert linked["cluster_id"] == cluster_id and field["cluster_id"] is None
        for star, name in ((linked, "Clusterholm"), (field, "Fieldholm")):
            system, cfg = _named_system(name)
            with conn:
                system_id = store.insert_star_system(conn, system, cfg)
                store.mark_bright_star_filled(conn, star["id"], system_id)
            ids = {row["cluster_id"] for row in conn.execute("SELECT cluster_id FROM stars WHERE star_system_id = ?", (system_id,)).fetchall()}
            assert ids == {star["cluster_id"]}
    finally:
        conn.close()


def test_a_sectors_save_passes_the_cluster_of_a_pre_placed_star_to_its_stars(mysql_config):
    from planetgen.generation.config import SystemConfig
    from planetgen.generation.system import StarSystem
    from planetgen.galaxy.sector import SpaceSector

    conn = store.get_connection(mysql_config)
    try:
        with conn:
            store.replace_galaxy_components(conn, [_cluster("Saved")])
            cluster_id = store.get_galaxy_components(conn)[0].id
    finally:
        conn.close()
    sector = SpaceSector("Cluster Sector", edge_ly=20.0)
    for name, position, cluster in (("Memberholm", (1.0, 2.0, 3.0), cluster_id), ("Fieldfall", (-4.0, 0.0, 5.5), None)):
        cfg = SystemConfig()
        cfg.STAR_TYPE = "G2V"
        cfg.PLANETS = False
        cfg.BINARY_SYSTEM = False
        system = StarSystem(system_config=cfg)
        system.name = name
        entry = sector.add_system(system, position=position, system_config=cfg)
        entry.cluster_id = cluster
    store.save_sector(sector, config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        found = {row["name"]: row["cluster_id"] for row in conn.execute(
            "SELECT s.name AS name, st.cluster_id AS cluster_id FROM star_systems s JOIN stars st ON st.star_system_id = s.id").fetchall()}
    finally:
        conn.close()
    assert found == {"Memberholm": cluster_id, "Fieldfall": None}
