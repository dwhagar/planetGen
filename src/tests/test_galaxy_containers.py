# tests/test_galaxy_containers.py

"""
DB.85: the globular clusters and superbubbles of a galaxy are rows of `galaxy_containers`, named by an object ID of the
`SERIAL_REGION` kind. The draw (`generation/containers.py`) is deterministic; the rows answer which containers hold a
point and what a container holds (`db/containers.py`), and the API speaks printed object IDs.
"""

import dataclasses
import math

import pytest

from planetgen.api.config import Config
from planetgen.db import containers as containersDb, store
from planetgen.galaxy import object_uid, spherical, superbubble_field
from planetgen.galaxy.density import build_galaxy_shape
from planetgen.galaxy.sector import SpaceSector
from planetgen.generation import clusters, containers
from planetgen.generation.config import SystemConfig
from planetgen.generation.phenomena.nebula import Nebula
from planetgen.generation.system import StarSystem
from planetgen.web.app import create_app

SHAPE = build_galaxy_shape(
    disk_scale_length_pc=2600.0, disk_scale_height_pc=300.0, bulge_scale_radius_pc=1580.0, bulge_amplitude=3.11,
    arm_count=2, pitch_angle_rad=math.radians(15.0), arm_amplitude=0.4)
SEED = bytes(range(16))
EDGE_PC = 50.0
REACH_PC = 6000.0


def _clusters(count=4):
    return [dataclasses.replace(c, id=index + 1) for index, c in enumerate(clusters.generate_clusters(SEED, count=count))]


def _built(count=4):
    return containers.build(SEED, SHAPE, EDGE_PC, REACH_PC, _clusters(count))


def test_the_same_seed_gives_the_same_containers_in_the_same_order():
    found = _built()
    assert found == _built()
    assert [c.kind for c in found].count(containers.GLOBULAR_CLUSTER) == 4
    assert found[0].kind == containers.GLOBULAR_CLUSTER and found[-1].kind == containers.SUPERBUBBLE
    assert containers.build(bytes(range(1, 17)), SHAPE, EDGE_PC, REACH_PC, _clusters()) != found


def test_a_bubble_container_is_a_bubble_of_the_field_and_a_cluster_container_is_a_king_cluster():
    found = _built()
    bubbles = {(b.index[0], b.index[1], b.n): b for b in superbubble_field.all_bubbles(SEED, SHAPE, REACH_PC)}
    assert sum(1 for c in found if c.kind == containers.SUPERBUBBLE) == len(bubbles) > 0
    for c in found:
        if c.kind == containers.SUPERBUBBLE:
            bubble = bubbles[c.bubble]
            assert c.center_pc == tuple(bubble.center_pc) and c.extent_pc == bubble.radius_pc
            assert c.name is None and c.component_id is None
        else:
            assert c.bubble is None and c.name and c.component_id in (1, 2, 3, 4) and c.extent_pc > 0


def test_an_open_cluster_is_not_a_container():
    drawn = _clusters(2)
    plummer = dataclasses.replace(drawn[0], profile="plummer", king_w0=None, name="Open cluster 00001")
    found = containers.build(None, None, EDGE_PC, REACH_PC, [plummer, drawn[1]])
    assert [c.name for c in found] == [drawn[1].name]


def test_an_object_id_is_the_birth_sector_the_region_kind_and_the_rank_in_that_sector():
    found = _built()
    seen = set()
    for c in found:
        assert c.birth == spherical.sector_address_at(c.center_pc, EDGE_PC)
        unpacked = object_uid.unpack(object_uid.from_bytes(c.uid))
        assert (unpacked.shell, unpacked.band, unpacked.wedge) == c.birth
        assert unpacked.serial_kind == object_uid.SERIAL_REGION
        seen.add(c.uid)
    assert len(seen) == len(found)


def test_two_containers_born_in_one_sector_take_the_next_rank():
    drawn = _clusters(2)
    near = dataclasses.replace(drawn[1], center_pc=(drawn[0].center_pc[0] + 0.5, drawn[0].center_pc[1],
                                                    drawn[0].center_pc[2]))
    found = containers.build(None, None, EDGE_PC, REACH_PC, [drawn[0], near])
    assert found[0].birth == found[1].birth
    ranks = [object_uid.unpack(object_uid.from_bytes(c.uid)).serial for c in found]
    assert sorted(ranks) == [0, 1]


@pytest.fixture
def galaxy(mysql_config):
    """Components and containers stored; one system at the first cluster's centre, one far off, one nebula."""
    drawn = clusters.generate_clusters(SEED, count=3)
    conn = store.get_connection(mysql_config)
    try:
        store.replace_galaxy_components(conn, drawn, kind="cluster")
        conn.commit()
        stored = store.get_galaxy_components(conn, kind="cluster")
        made = containers.build(SEED, SHAPE, EDGE_PC, REACH_PC, stored)
        containersDb.replace(conn, made)
        conn.commit()
    finally:
        conn.close()
    sector = SpaceSector("Container Sector", edge_ly=20.0)
    cfg = SystemConfig()
    cfg.BINARY_SYSTEM = False
    sector.add_system(StarSystem(system_config=cfg), position=(1.0, 1.0, 1.0), system_config=cfg)
    sector.add_system(StarSystem(system_config=cfg), position=(-2.0, 0.5, 3.0), system_config=cfg)
    sector_row = store.save_sector(sector, config=mysql_config)
    first = next(c for c in made if c.kind == containers.GLOBULAR_CLUSTER)
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            inside, outside = [row["id"] for row in conn.execute("SELECT id FROM star_systems ORDER BY id").fetchall()]
            x, y, z = first.center_pc
            conn.execute("UPDATE star_systems SET galaxy_x_pc = ?, galaxy_y_pc = ?, galaxy_z_pc = ? WHERE id = ?",
                         (x, y, z, inside))
            conn.execute("UPDATE star_systems SET galaxy_x_pc = ?, galaxy_y_pc = ?, galaxy_z_pc = ? WHERE id = ?",
                         (x + 10.0 * first.extent_pc, y, z, outside))
            veil = store.insert_nebula(conn, Nebula(SystemConfig(), name="Cluster Veil"), sector_id=sector_row,
                                       placement={"center_x_pc": x + 0.1 * first.extent_pc, "center_y_pc": y,
                                                  "center_z_pc": z, "galactic_radius_pc": 3.7})
        row = conn.execute("SELECT id FROM galaxy_containers WHERE uid = ?", (first.uid,)).fetchone()["id"]
    finally:
        conn.close()
    return {"made": made, "first": first, "row": row, "inside": inside, "outside": outside, "veil": veil}


def test_the_stored_rows_are_the_built_containers(mysql_config, galaxy):
    conn = store.get_connection(mysql_config)
    try:
        listed = containersDb.listing(conn)
        assert len(listed) == len(galaxy["made"])
        assert len(containersDb.listing(conn, containers.GLOBULAR_CLUSTER)) == 3
        assert {c["uid"] for c in listed} == {object_uid.format_id(object_uid.from_bytes(c.uid)) for c in galaxy["made"]}
        first = containersDb.get(conn, container_id=galaxy["row"])
        assert first["kind"] == "globular_cluster" and first["name"] == galaxy["first"].name and "bubble" not in first
        assert containersDb.get(conn, uid=first["uid"]) == first
        assert containersDb.get(conn, container_id=10 ** 9) is None
        with pytest.raises(ValueError):
            containersDb.listing(conn, "galaxy")
        containersDb.replace(conn, galaxy["made"][:2])
        assert len(containersDb.listing(conn)) == 2
        conn.rollback()
    finally:
        conn.close()


def test_a_cluster_holds_the_system_at_its_centre_and_not_the_one_far_off(mysql_config, galaxy):
    conn = store.get_connection(mysql_config)
    try:
        found = containersDb.contents(conn, galaxy["row"])
        assert found["counts"]["system"] == 1 and found["counts"]["nebula"] == 1 and found["total"] == 2
        assert {(i["kind"], i["id"]) for i in found["items"]} == {("system", galaxy["inside"]),
                                                                  ("nebula", galaxy["veil"])}
        assert len(containersDb.contents(conn, galaxy["row"], limit=1)["items"]) == 1
        with pytest.raises(ValueError):
            containersDb.contents(conn, 10 ** 9)
    finally:
        conn.close()


def test_the_containers_holding_an_object_are_found_deepest_first(mysql_config, galaxy):
    conn = store.get_connection(mysql_config)
    try:
        held = containersDb.holding_object(conn, "system", galaxy["inside"])
        assert held and held[0]["id"] == galaxy["row"] and held[0]["depth"] == pytest.approx(0.0, abs=1e-9)
        assert [h["depth"] for h in held] == sorted(h["depth"] for h in held)
        assert all(h["id"] != galaxy["row"] for h in containersDb.holding_object(conn, "system", galaxy["outside"]))
        assert containersDb.holding_object(conn, "system", 10 ** 9) == []
        assert containersDb.holding_object(conn, "nebula", galaxy["veil"])[0]["id"] == galaxy["row"]
        only = containersDb.holding(conn, galaxy["first"].center_pc, ["superbubble"])
        assert all(h["kind"] == "superbubble" for h in only)
        with pytest.raises(ValueError):
            containersDb.holding(conn, (0, 0, 0), ["moon"])
        with pytest.raises(ValueError):
            containersDb.holding_object(conn, "moon", 1)
    finally:
        conn.close()


@pytest.fixture
def client(mysql_config, galaxy):
    class TestConfig(Config):
        MYSQL_CONFIG = mysql_config
        WRITE_MYSQL_CONFIG = mysql_config
        CONTROL_MYSQL_CONFIG = mysql_config

    app = create_app(TestConfig)
    app.testing = True
    return app.test_client()


def test_the_routes_speak_printed_object_ids(client, galaxy):
    listed = client.get("/api/containers").get_json()["containers"]
    assert len(listed) == len(galaxy["made"])
    printed = object_uid.format_id(object_uid.from_bytes(galaxy["first"].uid))
    one = client.get(f"/api/containers/{printed}").get_json()
    assert one["uid"] == printed and one["kind"] == "globular_cluster"
    assert client.get("/api/containers?kind=globular_cluster").get_json()["containers"]
    assert client.get("/api/containers?kind=moon").status_code == 400
    assert client.get("/api/containers/FFFFFFFFFF-0000000-000").status_code == 404
    assert client.get("/api/containers/not-an-id").status_code == 404
    contents = client.get(f"/api/containers/{printed}/contents?limit=1").get_json()
    assert contents["total"] == 2 and len(contents["items"]) == 1
    assert client.get("/api/containers/FFFFFFFFFF-0000000-000/contents").status_code == 404


def test_the_at_route_lists_the_containers_around_a_point(client, galaxy):
    x, y, z = galaxy["first"].center_pc
    found = client.get(f"/api/containers/at?x={x}&y={y}&z={z}").get_json()["containers"]
    assert found and found[0]["uid"] == object_uid.format_id(object_uid.from_bytes(galaxy["first"].uid))
    assert client.get(f"/api/containers/at?x={x}&y={y}").status_code == 400
    assert client.get(f"/api/containers/at?x=nan&y={y}&z={z}").status_code == 400
    assert client.get(f"/api/containers/at?x={x}&y={y}&z={z}&kind=moon").status_code == 400
