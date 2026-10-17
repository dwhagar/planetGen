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
from planetgen.galaxy.spherical import local_to_galaxy_pc, sector_position_pc
from planetgen.generation import run_galaxy
from planetgen.physics.units import ly_to_pc, mpc_to_pc
from planetgen.web.app import create_app
from planetgen import tuning
from planetgen.galaxy.density import build_galaxy_shape
from tests.publicids import pid, pids

EDGE_PC = ly_to_pc(tuning.DEFAULT_SECTOR_EDGE_LY)
SHAPE = build_galaxy_shape(
    disk_scale_length_pc=40.0, disk_scale_height_pc=12.0, bulge_scale_radius_pc=10.0,
    bulge_amplitude=2.0, arm_count=2, pitch_angle_rad=math.radians(15), arm_amplitude=0.4,
)
E_VALUE = 2000.0
EXTENTS = [(0, 0)] + [(shell, 1) for shell in range(1, 9)]

ADDRESSES = [(2, 0, 0), (2, 0, 1), (3, 0, 4)]


@pytest.fixture
def galaxy(mysql_config):
    store.save_galaxy_shape(SHAPE, edge_pc=EDGE_PC, outer_shell_index=8,
                            expected_system_count_at_density_1=E_VALUE, config=mysql_config, cell_star_ceiling=1e12)
    store.replace_galaxy_shells(EXTENTS, config=mysql_config)
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
        # The stored offset follows the sector's own axes (MAP.186).
        found[row["id"]] = (row["name"], local_to_galaxy_pc(
            (row["cx"], row["cy"], row["cz"]), (mpc_to_pc(row["px"]), mpc_to_pc(row["py"]), mpc_to_pc(row["pz"]))))
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


def _scatter_rows(config, center):
    """Two pre-placed bright stars (one built into a system) and an unbuilt scattered black hole near `center`, and
    a bright star and a black hole 40 pc away."""
    from planetgen.galaxy.spherical import sector_address_at
    conn = store.get_connection(config)
    try:
        for offset in (0.0, 1.0, 40.0):
            point = (center[0] + offset, center[1], center[2])
            store.insert_bright_stars(conn, [(*sector_address_at(point, EDGE_PC), *(int(round(v * 1000)) for v in point),
                                              "young", "B", "III", 1e31, 1e7, 15000.0, 1e28, 0.05, 0.1, 8.0, None, 1, None)])
        for offset in (2.0, 40.0):
            point = (center[0] + offset, center[1], center[2])
            store.insert_phenomenon_scatter(conn, [(7, 2, 3, "black-hole", "stellar",
                                                    *(int(round(v * 1000)) for v in point), None, None, None, 99,
                                                    None, None)])
        conn.commit()
    finally:
        conn.close()


def test_the_search_counts_what_is_not_generated_yet_and_never_reads_the_sectors_to_search(galaxy):
    center = sector_position_pc(*ADDRESSES[0], EDGE_PC)
    _scatter_rows(galaxy, center)
    conn = _conn(galaxy)
    try:
        place = near.place_from_point(center)
        statements = []
        real = conn.execute
        conn.execute = lambda sql, *args, **kwargs: statements.append(sql) or real(sql, *args, **kwargs)
        page = near.objects_within(conn, place, 12.0)
        conn.execute = real
        assert page["not_generated"] == {"bright_star": 2, "black-hole": 1} and page["not_generated_total"] == 3
        searched = [sql for sql in statements if "sectors" in sql.replace("\n", " ")]
        assert all("WHERE id IN" in sql and "VALUES" not in sql for sql in searched), "the sectors are named by id only"
        before = conn.execute("SELECT COUNT(*) n FROM sectors").fetchone()["n"]
        near.objects_within(conn, place, 12.0)
        assert conn.execute("SELECT COUNT(*) n FROM sectors").fetchone()["n"] == before
        wide = near.objects_within(conn, place, 50.0, kinds=["system"])
        assert wide["not_generated"] == {"bright_star": 3, "black-hole": 2}
    finally:
        conn.close()
    conn = store.get_connection(galaxy)
    try:
        system_id = conn.execute("SELECT MIN(id) id FROM star_systems").fetchone()["id"]
        conn.execute("UPDATE bright_stars SET star_system_id = ? ORDER BY id LIMIT 1", (system_id,))
        conn.commit()
    finally:
        conn.close()
    conn = _conn(galaxy)
    try:
        assert near.objects_within(conn, place, 12.0)["not_generated"] == {"bright_star": 1, "black-hole": 1}
    finally:
        conn.close()


def test_a_thing_outside_the_key_box_is_found_by_the_null_key_branch(galaxy):
    conn = store.get_connection(galaxy)
    try:
        conn.execute("UPDATE facilities SET center_x_pc = 300000.0, center_y_pc = 0.0, center_z_pc = 0.0, sector_id = NULL")
        conn.commit()
        assert conn.execute("SELECT cell_key FROM facilities").fetchone()["cell_key"] is None
    finally:
        conn.close()
    conn = _conn(galaxy)
    try:
        page = near.objects_within(conn, near.place_from_point((300003.0, 0.0, 0.0)), 5.0, kinds=["facility"])
        assert [row["kind"] for row in page["rows"]] == ["facility"] and page["rows"][0]["parent"] is None
        assert page["rows"][0]["distance_pc"] == pytest.approx(3.0)
        assert near.objects_within(conn, near.place_from_point((300010.0, 0.0, 0.0)), 5.0, kinds=["facility"])["total"] == 0
    finally:
        conn.close()


def test_a_system_with_no_sector_is_found_and_has_no_parent(galaxy):
    conn = store.get_connection(galaxy)
    try:
        row = conn.execute("SELECT id, galaxy_x_pc x, galaxy_y_pc y, galaxy_z_pc z FROM star_systems ORDER BY id LIMIT 1").fetchone()
        conn.execute("UPDATE star_systems SET sector_id = NULL WHERE id = ?", (row["id"],))
        conn.commit()
    finally:
        conn.close()
    conn = _conn(galaxy)
    try:
        page = near.objects_within(conn, near.place_from_point((row["x"], row["y"], row["z"])), 0.5, kinds=["system"])
        found = [r for r in page["rows"] if r["id"] == row["id"]]
        assert found and found[0]["parent"] is None
    finally:
        conn.close()


def test_the_statement_reads_the_covering_index_with_the_run_list_first(galaxy):
    """The plan a MySQL upgrade must keep: the run list drives, each table is read by its covering cell index (a range
    per run), and the NULL-key branch reads the same covering index -- never the table itself."""
    import json
    from planetgen.db import region
    conn = _conn(galaxy)
    try:
        for table, position in (("star_systems", ("p.galaxy_x_pc", "p.galaxy_y_pc", "p.galaxy_z_pc")),
                                ("facilities", ("p.center_x_pc", "p.center_y_pc", "p.center_z_pc")),
                                ("nebulae", ("p.center_x_pc", "p.center_y_pc", "p.center_z_pc"))):
            sql, runs = region.region_statement(table, "p.id", position, (8200.0, 0.0, 0.0), 20.0)
            assert 1 <= len(runs) <= 256
            plan = json.loads(list(conn.execute("EXPLAIN FORMAT=JSON " + sql).fetchone().values())[0])
            first, second = plan["query_block"]["union_result"]["query_specifications"]
            order = [entry["table"] for entry in first["query_block"]["nested_loop"]]
            assert [t["table_name"] for t in order] == ["rg", "p"], table
            probe = order[1]
            assert probe["possible_keys"] == [region.index_name(table)], table
            assert "range_checked_for_each_record" in probe or probe.get("key") == region.index_name(table), \
                f"{table}: a full scan of the table"
            null_branch = second["query_block"]["table"]
            assert null_branch["access_type"] in ("range", "ref", "index") and \
                null_branch["key"] == region.index_name(table) and null_branch.get("using_index"), table
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
    response = client.get(f"/api/near?from=system:{pid('system', system_id)}&distance=8&kinds=system,planet&limit=5")
    assert response.status_code == 200
    body = response.get_json()
    assert body["place"]["ref"] == f"system:{pid('system', system_id)}" and len(body["rows"]) <= 5
    assert body["total"] >= len(body["rows"])
    assert {row["kind"] for row in body["rows"]} <= {"system", "planet"}

    point = client.get("/api/near?point=0,0,0&distance=2")
    assert point.status_code == 200 and point.get_json()["place"]["ref"] is None

    for bad in ("/api/near?distance=5", "/api/near?point=0,0,0", "/api/near?from=system:1&point=0,0,0&distance=5",
                "/api/near?point=0,0,0&distance=51", "/api/near?point=0,0,0&distance=nan",
                "/api/near?point=0,0&distance=5", "/api/near?point=0,0,0&distance=5&kinds=ship",
                "/api/near?point=0,0,0&distance=5&limit=x", "/api/near?from=system:99999999&distance=5"):
        assert client.get(bad).status_code == 400, bad
