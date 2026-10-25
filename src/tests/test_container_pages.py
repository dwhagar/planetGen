# tests/test_container_pages.py

"""
UX.184: a globular cluster and a superbubble are objects a reference can name (`globular_cluster:<ID>`,
`superbubble:<ID>`): the reference resolves by the container's own object ID and by its kind, its property sheet is
served, and its page (`/object/<ref>`) shows where it lies, what it is and what sits inside it.
"""

import math

import pytest

from planetgen.api import ids
from planetgen.api.config import Config
from planetgen.db import containers as containersDb, store
from planetgen.galaxy import object_uid, objectref
from planetgen.galaxy.density import build_galaxy_shape
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation import clusters, containers
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.web import item_pages
from planetgen.web.app import create_app

SHAPE = build_galaxy_shape(
    disk_scale_length_pc=2600.0, disk_scale_height_pc=300.0, bulge_scale_radius_pc=1580.0, bulge_amplitude=3.11,
    arm_count=2, pitch_angle_rad=math.radians(15.0), arm_amplitude=0.4)
SEED = bytes(range(16))


def _printed(container):
    return object_uid.format_id(object_uid.from_bytes(container.uid))


@pytest.fixture
def world(mysql_config):
    """A client over a galaxy with stored containers and one system at the first cluster's centre, and the printed IDs
    of the first cluster and the first superbubble."""
    conn = store.get_connection(mysql_config)
    try:
        store.replace_galaxy_components(conn, clusters.generate_clusters(SEED, count=3), kind="cluster")
        conn.commit()
        made = containers.build(SEED, SHAPE, 50.0, 6000.0, store.get_galaxy_components(conn, kind="cluster"))
        containersDb.replace(conn, made)
        conn.commit()
    finally:
        conn.close()
    sector = SpaceSector("Cluster Sector", edge_ly=20.0)
    cfg = SystemConfig()
    cfg.BINARY_SYSTEM = False
    sector.add_system(StarSystem(system_config=cfg), position=(1.0, 1.0, 1.0), system_config=cfg)
    store.save_sector(sector, config=mysql_config)
    cluster = next(c for c in made if c.kind == containers.GLOBULAR_CLUSTER)
    bubble = next(c for c in made if c.kind == containers.SUPERBUBBLE)
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            x, y, z = cluster.center_pc
            conn.execute("UPDATE star_systems SET galaxy_x_pc = ?, galaxy_y_pc = ?, galaxy_z_pc = ?", (x, y, z))
    finally:
        conn.close()

    class TestConfig(Config):
        MYSQL_CONFIG = mysql_config
        WRITE_MYSQL_CONFIG = mysql_config
        CONTROL_MYSQL_CONFIG = mysql_config

    app = create_app(TestConfig)
    app.testing = True
    return app.test_client(), cluster, bubble


def test_a_cluster_and_a_superbubble_are_kinds_a_reference_can_name():
    for kind in ("globular_cluster", "superbubble"):
        assert kind in objectref.KINDS and kind in item_pages.page_kinds() and kind in ids.ACTIVE_KINDS
        assert objectref.parse_public(f"{kind}:FE81000A2B-0000005-007") == (kind, "FE81000A2B-0000005-007")


def test_a_cluster_resolves_by_its_object_id_and_never_as_a_superbubble(world):
    client, cluster, bubble = world
    found = client.get(f"/api/objects/globular_cluster:{_printed(cluster)}").get_json()
    assert found["kind"] == "globular_cluster" and found["name"] == cluster.name
    assert [parent["kind"] for parent in found["parents"]] == ["galaxy"]
    assert found["positions"]["galaxy_pc"] == list(cluster.center_pc)
    assert found["id"] == _printed(cluster)
    assert client.get(f"/api/objects/superbubble:{_printed(cluster)}").status_code == 404
    assert client.get(f"/api/objects/globular_cluster:{_printed(bubble)}").status_code == 404
    assert client.get(f"/api/objects/superbubble:{_printed(bubble)}").get_json()["kind"] == "superbubble"


def test_the_sheet_serves_the_stored_extent_and_a_container_holds_no_children(world):
    client, cluster, _bubble = world
    ref = f"globular_cluster:{_printed(cluster)}"
    sheet = client.get(f"/api/objects/{ref}/sheet").get_json()
    rows = {row["key"]: row for group in sheet["groups"] for row in group["rows"]}
    assert sheet["kind"] == "globular_cluster" and rows["extent_pc"]["stored"] is True
    children = client.get(f"/api/objects/{ref}/children").get_json()
    assert children["total"] == 0 and children["items"] == []


def test_the_container_api_carries_a_clusters_component_and_the_printed_ids_of_what_it_holds(world):
    client, cluster, bubble = world
    one = client.get(f"/api/containers/{_printed(cluster)}").get_json()
    assert one["component"]["age_gy"] > 0 and one["component"]["mass_solar"] > 0
    assert "component" not in client.get(f"/api/containers/{_printed(bubble)}").get_json()
    held = client.get(f"/api/containers/{_printed(cluster)}/contents").get_json()
    system = next(item for item in held["items"] if item["kind"] == "system")
    assert "-" in str(system["id"])


def test_the_page_of_a_cluster_says_what_it_is_and_links_what_sits_inside_it(world):
    client, cluster, _bubble = world
    page = client.get(f"/object/globular_cluster:{_printed(cluster)}")
    assert page.status_code == 200
    html = page.get_data(as_text=True)
    assert cluster.name in html and "Globular cluster" in html and "What is inside it" in html
    assert 'id="item-properties"' in html and 'id="location"' not in html    # UX.200: nothing above a cluster but the galaxy
    assert 'href="/system/' in html and "billion years old" in html


def test_the_page_of_a_superbubble_opens_and_a_wrong_kind_is_a_404(world):
    client, cluster, bubble = world
    page = client.get(f"/object/superbubble:{_printed(bubble)}")
    assert page.status_code == 200 and "Superbubble" in page.get_data(as_text=True)
    assert client.get(f"/object/superbubble:{_printed(cluster)}").status_code == 404
    assert client.get("/object/globular_cluster:FFFFFFFFFF-0000000-000").status_code == 404
