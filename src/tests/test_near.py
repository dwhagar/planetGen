# tests/test_near.py

"""
NAV.43: `planetgen.db.near.objects_within`, `GET /api/near` and `planetgen.cli.query near`
against a small generated galaxy, checked against distances worked out
straight from the stored rows.
"""

import math

import pytest

from planetgen.api.config import Config
from planetgen.db import near, query, store
from planetgen.galaxy.geometry import sector_position_pc
from planetgen.generation import run_galaxy
from planetgen.physics.units import ly_to_pc, mpc_to_pc
from planetgen.web.app import create_app
from planetgen import tuning
from planetgen.galaxy.density import build_galaxy_shape

EDGE_PC = ly_to_pc(tuning.DEFAULT_SECTOR_EDGE_LY)
SHAPE = build_galaxy_shape(
    disk_scale_length_pc=40.0, disk_scale_height_pc=12.0, bulge_scale_radius_pc=10.0,
    bulge_amplitude=2.0, arm_count=2, pitch_angle_rad=math.radians(15), arm_amplitude=0.4,
)
E_VALUE = 2000.0
EXTENTS = [(1, 6), (0, 8), (-1, 6)]

ADDRESSES = [(2, 0, 0), (2, 0, 1), (3, 0, 4)]


@pytest.fixture
def galaxy(mysql_config):
    store.save_galaxy_shape(SHAPE, edge_pc=EDGE_PC, outer_ring_index=8,
                            expected_system_count_at_density_1=E_VALUE, config=mysql_config)
    store.replace_galaxy_layers(EXTENTS, config=mysql_config)
    args = run_galaxy._default_generation_args(config=mysql_config)
    args.num_systems = 4
    for address in ADDRESSES:
        run_galaxy.generate_and_save_sector_at(args, address, sector_position_pc(*address, EDGE_PC), EDGE_PC)
    conn = store.get_connection(mysql_config)
    try:
        sector_id = conn.execute("SELECT id FROM sectors ORDER BY id LIMIT 1").fetchone()["id"]
        store.add_facility(conn, "Far Light", "station", "standalone", "space", sector_id,
                           offset_ly=(1.0, 0.0, 0.0))
        conn.commit()
    finally:
        conn.close()
    return mysql_config


def _conn(config):
    return query.open_readonly(config)


def _systems(conn):
    found = {}
    for row in conn.execute(
            "SELECT ss.id, ss.name, sec.center_x_pc cx, sec.center_y_pc cy, sec.center_z_pc cz,"
            " ss.position_x_mpc px, ss.position_y_mpc py, ss.position_z_mpc pz"
            " FROM star_systems ss JOIN sectors sec ON sec.id = ss.sector_id").fetchall():
        found[row["id"]] = (row["name"], (row["cx"] + mpc_to_pc(row["px"]), row["cy"] + mpc_to_pc(row["py"]),
                                         row["cz"] + mpc_to_pc(row["pz"])))
    return found


def _everything(conn, place, distance, **kwargs):
    rows = []
    offset = 0
    while True:
        page = near.objects_within(conn, place, distance, limit=near.MAX_LIMIT, offset=offset, **kwargs)
        rows += page["rows"]
        offset += near.MAX_LIMIT
        if offset >= page["total"]:
            return page, rows


def test_systems_within_a_distance_are_exactly_the_ones_inside_the_sphere(galaxy):
    conn = _conn(galaxy)
    try:
        systems = _systems(conn)
        center = sector_position_pc(*ADDRESSES[0], EDGE_PC)
        distance = 6.0
        page, rows = _everything(conn, near.place_from_point(center), distance, kinds=["system"])
        expected = {f"system:{i}" for i, (_n, point) in systems.items() if math.dist(center, point) <= distance}
        assert {row["ref"] for row in rows} == expected
        assert expected and len(expected) < len(systems)  # across a sector boundary, but not everything
        away = [row["distance_pc"] for row in rows]
        assert away == sorted(away)
        for row in rows:
            assert row["distance_pc"] == pytest.approx(math.dist(center, systems[row["id"]][1]))
        assert page["total"] == len(expected) == page["by_kind"]["system"]
    finally:
        conn.close()


def test_bodies_ride_with_their_system_and_pages_add_up(galaxy):
    conn = _conn(galaxy)
    try:
        center = sector_position_pc(*ADDRESSES[0], EDGE_PC)
        place = near.place_from_point(center)
        page, rows = _everything(conn, place, 8.0)
        for kind, table in (("star", "stars"), ("planet", "planets"), ("moon", "moons"),
                            ("belt", "asteroid_belts"), ("comet", "comets")):
            in_range = {row["id"] for row in rows if row["kind"] == kind}
            ids = {int(r["ref"].split(":")[1]) for r in rows if r["kind"] == "system"}
            expected = {r["id"] for r in conn.execute(
                f"SELECT id FROM {table} WHERE star_system_id IN ({','.join(map(str, ids or [0]))})").fetchall()}
            assert in_range == expected, kind
        assert page["total"] == len(rows)
        assert len({row["ref"] for row in rows}) == len(rows)

        # Paging by small pages gives the same rows in the same order.
        paged = []
        for offset in range(0, page["total"], 7):
            paged += near.objects_within(conn, place, 8.0, limit=7, offset=offset)["rows"]
        assert [row["ref"] for row in paged] == [row["ref"] for row in rows]
        assert [row["distance_pc"] for row in rows] == sorted(row["distance_pc"] for row in rows)
    finally:
        conn.close()


def test_from_a_system_leaves_it_out_and_finds_its_neighbours(galaxy):
    conn = _conn(galaxy)
    try:
        systems = _systems(conn)
        system_id = next(iter(systems))
        place = near.place_from_reference(conn, str(system_id))
        assert place["ref"] == f"system:{system_id}" and place["name"] == systems[system_id][0]
        _page, rows = _everything(conn, place, 10.0, kinds=["system"])
        refs = {row["ref"] for row in rows}
        assert f"system:{system_id}" not in refs
        expected = {f"system:{i}" for i, (_n, point) in systems.items()
                    if i != system_id and math.dist(place["point_pc"], point) <= 10.0}
        assert refs == expected
    finally:
        conn.close()


def test_a_stand_alone_facility_is_found_by_its_own_position(galaxy):
    conn = _conn(galaxy)
    try:
        row = conn.execute("SELECT id, center_x_pc x, center_y_pc y, center_z_pc z FROM facilities").fetchone()
        place = near.place_from_point((row["x"], row["y"], row["z"]))
        page = near.objects_within(conn, place, 0.5, kinds=["facility"])
        assert [r["ref"] for r in page["rows"]] == [f"facility:{row['id']}"]
        assert page["rows"][0]["distance_pc"] == pytest.approx(0.0)
        assert page["rows"][0]["parent"]["kind"] == "sector"
    finally:
        conn.close()


def test_the_search_reports_sectors_it_did_not_find_generated(galaxy):
    conn = _conn(galaxy)
    try:
        place = near.place_from_point(sector_position_pc(*ADDRESSES[0], EDGE_PC))
        page = near.objects_within(conn, place, 12.0, kinds=["system"])
        assert page["sectors_in_range"] > page["sectors_generated"] >= len(ADDRESSES) - 1
        before = conn.execute("SELECT COUNT(*) n FROM sectors").fetchone()["n"]
        near.objects_within(conn, place, 12.0)
        assert conn.execute("SELECT COUNT(*) n FROM sectors").fetchone()["n"] == before
    finally:
        conn.close()


@pytest.mark.parametrize("distance", [0, -1, "x", float("nan"), float("inf"), near.MAX_DISTANCE_PC + 1])
def test_a_bad_distance_is_refused(galaxy, distance):
    conn = _conn(galaxy)
    try:
        with pytest.raises(near.NearError):
            near.objects_within(conn, near.place_from_point((0, 0, 0)), distance)
    finally:
        conn.close()


def test_other_bad_input_is_refused(galaxy):
    conn = _conn(galaxy)
    try:
        place = near.place_from_point((0, 0, 0))
        with pytest.raises(near.NearError):
            near.objects_within(conn, place, 5, kinds=["spaceship"])
        with pytest.raises(near.NearError):
            near.objects_within(conn, place, 5, limit=0)
        with pytest.raises(near.NearError):
            near.place_from_point((1, 2))
        with pytest.raises(near.NearError):
            near.place_from_point((1, 2, float("nan")))
        with pytest.raises(near.NearError):
            near.place_from_reference(conn, "system:99999999")
        with pytest.raises(near.NearError):
            near.place_from_reference(conn, "not a ref")
    finally:
        conn.close()


def test_the_largest_distance_searches_without_trouble(galaxy):
    conn = _conn(galaxy)
    try:
        place = near.place_from_point(sector_position_pc(*ADDRESSES[0], EDGE_PC))
        page = near.objects_within(conn, place, near.MAX_DISTANCE_PC, limit=10)
        assert page["total"] >= len(page["rows"]) == 10
    finally:
        conn.close()


@pytest.fixture
def client(galaxy):
    class TestConfig(Config):
        MYSQL_CONFIG = galaxy
        WRITE_MYSQL_CONFIG = galaxy
        CONTROL_MYSQL_CONFIG = galaxy

    app = create_app(TestConfig)
    app.testing = True
    return app.test_client()


def test_api_near_answers_and_refuses(client, galaxy):
    conn = _conn(galaxy)
    try:
        system_id = next(iter(_systems(conn)))
    finally:
        conn.close()
    response = client.get(f"/api/near?from=system:{system_id}&distance=8&kinds=system,planet&limit=5")
    assert response.status_code == 200
    body = response.get_json()
    assert body["place"]["ref"] == f"system:{system_id}" and len(body["rows"]) <= 5
    assert body["total"] >= len(body["rows"])
    assert {row["kind"] for row in body["rows"]} <= {"system", "planet"}

    point = client.get("/api/near?point=0,0,0&distance=2")
    assert point.status_code == 200 and point.get_json()["place"]["ref"] is None

    for bad in ("/api/near?distance=5", "/api/near?point=0,0,0", "/api/near?from=system:1&point=0,0,0&distance=5",
                "/api/near?point=0,0,0&distance=51", "/api/near?point=0,0,0&distance=nan",
                "/api/near?point=0,0&distance=5", "/api/near?point=0,0,0&distance=5&kinds=ship",
                "/api/near?point=0,0,0&distance=5&limit=x", "/api/near?from=system:99999999&distance=5"):
        assert client.get(bad).status_code == 400, bad
