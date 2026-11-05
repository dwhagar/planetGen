"""scripts/census_keep_out_overlap.py (NAV.56): the overlap census of keep-out spheres."""

import importlib.util
import pathlib

import numpy as np
import pytest

from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.asteroid_field import AsteroidField
from tests.test_api import client  # noqa: F401

_PATH = pathlib.Path(__file__).resolve().parents[2] / "scripts" / "census_keep_out_overlap.py"
_SPEC = importlib.util.spec_from_file_location("census_keep_out_overlap", _PATH)
census_script = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(census_script)


def test_no_overlap_means_no_cluster_fallback():
    found = census_script.census([[0, 0, 0], [10, 0, 0]], [1.0, 2.0], ["system", "system"])
    assert found["pairs"] == 0 and found["clusters"] == {} and found["in_overlap"] == 0
    assert "never fires" in census_script.render(found)
    assert census_script.census(np.zeros((0, 3)), [], [])["spheres"] == 0


def test_pairs_nested_spheres_and_clusters_are_counted():
    centers = [[0, 0, 0], [1.5, 0, 0], [3.0, 0, 0],       # a chain of three: 0-1 and 1-2 overlap, 0-2 do not
               [100, 0, 0], [100.2, 0, 0],                  # two nested
               [200, 0, 0]]                                  # alone
    radii = [1.0, 1.0, 1.0, 1.0, 0.1, 1.0]
    kinds = ["system", "system", "black_hole", "system", "rogue_planet", "quasar"]
    found = census_script.census(centers, radii, kinds)
    assert found["pairs"] == 3 and found["nested"] == 1
    assert found["in_overlap"] == 5
    assert found["clusters"] == {2: 1, 3: 1} and found["largest_cluster"] == 3
    assert found["pairs_by_kinds"] == {"black_hole+system": 1, "rogue_planet+system": 1, "system+system": 1}
    assert found["by_kind"]["system"] == 3


def test_the_reader_finds_the_placed_objects_of_a_database(client, mysql_config):
    conn = store.get_connection(mysql_config)
    try:
        field = AsteroidField(SystemConfig())
        store.insert_asteroid_field(conn, field, placement={
            "center_x_pc": 5.0, "center_y_pc": 1.0, "center_z_pc": 2.0, "galactic_radius_pc": 5.5})
        conn.commit()
        centers, radii, kinds = census_script.load_spheres(conn)
    finally:
        conn.close()
    assert "asteroid_field" in kinds
    index = kinds.index("asteroid_field")
    assert centers[index].tolist() == [5.0, 1.0, 2.0]
    assert radii[index] == pytest.approx(field.radius_ly / 3.261563777)


def test_the_reader_gives_compact_objects_the_galactic_hill_radius_of_the_query(client, mysql_config):
    from planetgen.db import query
    from planetgen.generation.phenomena.compact_remnant import BlackHole, NeutronStar
    from planetgen.generation.phenomena.rogue import RoguePlanet

    placement = {"center_x_pc": 8000.0, "center_y_pc": 0.0, "center_z_pc": 0.0, "galactic_radius_pc": 8000.0}
    conn = store.get_connection(mysql_config)
    try:
        ids = {"black_hole": store.insert_black_hole(conn, BlackHole(SystemConfig()), placement=placement),
               "neutron_star": store.insert_neutron_star(conn, NeutronStar(SystemConfig()), placement=placement),
               "rogue_planet": store.insert_rogue_planet(conn, RoguePlanet(SystemConfig()), placement=placement)}
        center = {"center_x_pc": 0.0, "center_y_pc": 0.0, "center_z_pc": 0.0, "galactic_radius_pc": 0.0}
        store.insert_black_hole(conn, BlackHole(SystemConfig(), mass_class="supermassive"), placement=center)
        conn.commit()
        centers, radii, kinds = census_script.load_spheres(conn)
        expected = {kind: query.keep_out_radius(conn, kind, object_id).radius_km / census_script.KM_PER_PC
                    for kind, object_id in ids.items()}
    finally:
        conn.close()
    for kind, want in expected.items():
        assert radii[kinds.index(kind)] == pytest.approx(want, rel=0.02), kind
    assert "supermassive_black_hole" in kinds
