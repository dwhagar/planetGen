# tests/test_api_galaxy_profile.py

"""API.27: `GET /api/galaxy/profile` answers for any cell of the grid, charted or not, and changes nothing."""

import json
import math

import pytest

from planetgen import tuning
from planetgen.db import store as _db
from planetgen.galaxy import superbubble_field
from planetgen.galaxy.density import build_galaxy_shape
from planetgen.galaxy.profile import POPULATIONS, sector_profile
from planetgen.galaxy.spherical import MAX_ADDRESSABLE_SHELL, sector_address_at, sector_position_pc
from tests.publicids import pid

# Imported fixtures (see test_api_thin_routes.py).
from tests.test_api import _place_sector, client  # noqa: F401

SEED = bytes(range(16))


def _plan(mysql_config, **kwargs):
    shape = build_galaxy_shape(disk_scale_length_pc=2600.0, disk_scale_height_pc=300.0, bulge_scale_radius_pc=1580.0,
                               bulge_amplitude=3.11, arm_count=2, pitch_angle_rad=math.radians(15), arm_amplitude=0.4)
    _db.save_galaxy_shape(shape, edge_pc=float(tuning.DEFAULT_SECTOR_EDGE_PC), outer_shell_index=1500,
                          expected_system_count_at_density_1=600.0, config=mysql_config, galaxy_seed=SEED, **kwargs)
    return _db.get_galaxy_shape(_db.get_connection(mysql_config))


def _tables(mysql_config):
    conn = _db.get_connection(mysql_config)
    try:
        return {row["TABLE_NAME"]: conn.execute(f"SELECT COUNT(*) AS n FROM `{row['TABLE_NAME']}`").fetchone()["n"]
                for row in conn.execute("SELECT TABLE_NAME FROM information_schema.TABLES WHERE TABLE_SCHEMA = DATABASE()"
                                        " AND TABLE_TYPE = 'BASE TABLE'").fetchall()}
    finally:
        conn.close()


def test_the_profile_of_an_uncharted_cell_is_the_programs_own_profile(client, mysql_config):
    skeleton = _plan(mysql_config)
    address = (1500 // 4, 3, 17)  # nothing was generated here
    body = client.get("/api/galaxy/profile?shell=%d&band=%d&wedge=%d" % address).get_json()
    expected = sector_profile(address, skeleton)
    assert (body["shell_index"], body["band_index"], body["wedge_index"]) == address
    assert body["sector_id"] is None and body["in_galaxy"] in (True, False)
    assert body["region"] == expected.region
    assert body["density"] == pytest.approx(expected.density) and body["expected_systems"] == pytest.approx(
        expected.expected_systems)
    assert body["center_pc"] == pytest.approx(list(expected.position_pc))
    assert body["volume_pc3"] == pytest.approx(expected.volume_pc3) and body["edge_pc"] == skeleton.edge_pc
    assert set(body["populations"]) == set(POPULATIONS) == set(body["mix"])
    assert sum(body["mix"].values()) == pytest.approx(body["density"], rel=1e-9)
    young = body["populations"]["young"]
    assert set(young) == {"density", "share", "age_gy", "fe_h", "dispersion_kms", "rotation"}
    assert len(young["age_gy"]) == 2 and all(len(peak) == 3 for peak in young["fe_h"])
    assert set(body["gas"]) == {"cold_fraction", "ambient", "ambient_density_cm3", "mass_solar", "superbubble"}
    assert body["dark_matter"]["mass_solar"] == pytest.approx(expected.dark_matter.mass_solar)
    assert body["kind_factors"] == pytest.approx(expected.kind_factors)
    assert "unresolved" in body
    json.dumps(body, allow_nan=False)  # every number is finite


def test_a_point_names_the_cell_that_holds_it_and_a_charted_cell_names_its_sector(client, mysql_config):
    skeleton = _plan(mysql_config)
    edge_pc = skeleton.edge_pc
    address = (3, -1, 5)
    center = sector_position_pc(*address, edge_pc)
    sector_id = _place_sector(mysql_config, "Profiled", center, edge_ly=tuning.DEFAULT_SECTOR_EDGE_LY, address=address)
    by_address = client.get("/api/galaxy/profile?shell=3&band=-1&wedge=5").get_json()
    by_point = client.get(f"/api/galaxy/profile?x={center[0] + 0.3}&y={center[1]}&z={center[2] - 0.2}").get_json()
    assert by_address["sector_id"] == pid("sector", sector_id)
    assert by_point == by_address
    assert sector_address_at(center, edge_pc) == address


def test_a_place_inside_a_superbubble_reports_it(client, mysql_config):
    skeleton = _plan(mysql_config)
    bubble = next(b for b in superbubble_field.all_bubbles(SEED, skeleton.shape, 5000.0) if b.radius_pc > 200.0)
    x, y, z = bubble.center_pc
    body = client.get(f"/api/galaxy/profile?x={x}&y={y}&z={z}").get_json()
    found = body["gas"]["superbubble"]
    assert found is not None and found["radius_pc"] == pytest.approx(bubble.radius_pc)
    assert body["gas"]["cold_fraction"] == 0.0


def test_the_lookup_reads_and_writes_nothing(client, mysql_config):
    _plan(mysql_config)
    before = _tables(mysql_config)
    first = client.get("/api/galaxy/profile?shell=40&band=0&wedge=3").get_json()
    second = client.get("/api/galaxy/profile?shell=40&band=0&wedge=3").get_json()
    assert first == second and _tables(mysql_config) == before


def test_the_profile_of_the_outermost_cells_of_the_grid_is_answered(client, mysql_config):
    _plan(mysql_config)
    for address in ((0, 0, 0), (MAX_ADDRESSABLE_SHELL, 0, 0), (MAX_ADDRESSABLE_SHELL, -1, 7)):
        response = client.get("/api/galaxy/profile?shell=%d&band=%d&wedge=%d" % address)
        if response.status_code == 400:  # an address that is not a cell is refused, never a 500
            assert "error" in response.get_json()
        else:
            assert response.status_code == 200
            json.dumps(response.get_json(), allow_nan=False)


def test_bad_queries_are_clean_400s_and_an_unplanned_galaxy_is_a_404(client, mysql_config):
    assert client.get("/api/galaxy/profile?shell=3&band=0&wedge=1").status_code == 404
    _plan(mysql_config)
    for query in ("", "shell=1", "shell=1&band=0", "shell=a&band=0&wedge=0", "x=nan&y=0&z=0", "x=1&y=2",
                  "shell=-1&band=0&wedge=0", "shell=0&band=0&wedge=4", f"shell={MAX_ADDRESSABLE_SHELL + 1}&band=0&wedge=0",
                  "x=1e300&y=0&z=0", "x=inf&y=0&z=0"):
        response = client.get("/api/galaxy/profile?" + query)
        assert response.status_code == 400, query
        assert "error" in response.get_json(), query
