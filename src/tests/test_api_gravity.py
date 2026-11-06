# src/tests/test_api_gravity.py

"""
API.24: `GET /api/sectors/<id>/gravity`, the per-zone gravity grid as JSON (`api/gravity.py`, `galaxy/gravity_field.py`, `db/gravity_stamp.py`).

* `gravity_field`: the smooth galaxy's tidal tensor against the Poisson density, the three modes' numbers against closed forms, the saddle between two stars;
* `gravity_stamp`: the numba-free half of the grid job names the file the job writes;
* the route: its refusals, the `202` job and its reuse, the `200` answer of a cached grid, the total and the local field, and a stale grid built again.
"""

import math

import numpy as np
import pytest

from planetgen.api import gravity as gravity_api
from planetgen.api.config import Config
from planetgen.db import gravity_grid, gravity_stamp, store
from planetgen.galaxy import gravity_cache, gravity_field, gravity_kernel, potential
from planetgen.queue import api_jobs
from planetgen.web.app import create_app
from planetgen.web.lib import tilecache
from tests.publicids import pid
from tests.test_gravity_grid import EDGE_PC, POT, _save
from tests.toy_shape import toy_shape

G = potential.G_PC_KMS2_PER_MSUN
HERNQUIST = potential.Potential(spheres=(potential.Hernquist("bulge", 1.0e10, 700.0),))


def test_the_smooth_tidal_tensor_has_the_trace_of_minus_four_pi_g_rho():
    points = np.array([[1500.0, 400.0, -200.0], [3000.0, 0.0, 100.0], [800.0, 800.0, 800.0]])
    tensor = gravity_field.smooth_fields(points, HERNQUIST)
    trace = tensor[:, gravity_field.T_XX] + tensor[:, gravity_field.T_YY] + tensor[:, gravity_field.T_ZZ]
    np.testing.assert_allclose(trace, -4.0 * math.pi * G * potential.stellar_density(HERNQUIST, *points.T), rtol=2e-4)
    np.testing.assert_allclose(tensor[:, gravity_field.POTENTIAL], potential.potential(HERNQUIST, *points.T), rtol=1e-12)


def test_the_smooth_tensor_is_the_derivative_of_the_smooth_pull():
    point = np.array([[2500.0, 900.0, 300.0]])
    tensor = gravity_field.smooth_fields(point, POT)[0]
    h = 2.0
    for column, offset in enumerate(np.eye(3) * h):
        plus = np.array(potential.acceleration(POT, *(point[0] + offset)))
        minus = np.array(potential.acceleration(POT, *(point[0] - offset)))
        derivative = (plus - minus) / (2.0 * h)
        expected = [tensor[gravity_field.T_XX], tensor[gravity_field.T_XY], tensor[gravity_field.T_XZ]] if column == 0 else (
            [tensor[gravity_field.T_XY], tensor[gravity_field.T_YY], tensor[gravity_field.T_YZ]] if column == 1 else
            [tensor[gravity_field.T_XZ], tensor[gravity_field.T_YZ], tensor[gravity_field.T_ZZ]])
        np.testing.assert_allclose(derivative, expected, rtol=1e-3, atol=1e-9)


def _two_stars(softening=1.0e-3):
    c = (np.arange(16) + 0.5) / 16 * 4.0 - 2.0
    x, y, z = np.meshgrid(c, c, c, indexing="ij")
    points = np.stack([x.ravel(), y.ravel(), z.ravel()], axis=1)
    stars = np.array([[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    return points, gravity_kernel.fields(points, stars, np.array([1.0, 1.0]), softening)


def test_the_three_modes_are_the_pull_the_depth_and_the_tidal_stretch_in_si():
    points, local = _two_stars()
    pull = gravity_field.scalar(local, "pull")
    near = int(np.argmin(np.linalg.norm(points - [-1.0, 0.25, 0.0], axis=1)))
    distance = np.linalg.norm(points[near] - [-1.0, 0.0, 0.0])
    expected = G * 1.0 / distance ** 2 * potential.ACCELERATION_TO_SI
    assert pull[near] == pytest.approx(expected, rel=0.1)           # the other star pulls a little the other way
    depth = gravity_field.scalar(local, "depth")
    assert depth.min() == 0.0 and depth.max() > 0.0 and (depth >= 0.0).all()
    assert depth[near] > depth[int(np.argmin(np.linalg.norm(points, axis=1)))]
    tidal = gravity_field.scalar(local, "tidal")
    assert tidal[near] == pytest.approx(2.0 * G / distance ** 3 * gravity_field.TENSOR_TO_SI, rel=0.2)
    with pytest.raises(ValueError):
        gravity_field.scalar(local, "colour")


def test_the_total_field_is_the_local_one_plus_the_smooth_galaxy():
    points, local = _two_stars()
    points = points + [8000.0, 0.0, 0.0]
    total = gravity_field.field(local, points, HERNQUIST)
    np.testing.assert_allclose(total, local + gravity_field.smooth_fields(points, HERNQUIST))
    np.testing.assert_array_equal(gravity_field.field(local, points, HERNQUIST, total=False), local)
    # the smooth pull swamps the stars' in the pull, the tidal stretch keeps them
    assert gravity_field.scalar(total, "pull").std() / gravity_field.scalar(total, "pull").mean() < 0.02
    assert gravity_field.scalar(total, "tidal").max() > 100.0 * abs(gravity_field.scalar(gravity_field.smooth_fields(points, HERNQUIST), "tidal")).max()


def test_the_zone_between_two_equal_stars_is_the_saddle():
    _points, local = _two_stars()
    found = gravity_field.saddles(local, 16)
    assert len(found) >= 1
    i, j, k = found[0]
    assert i in (7, 8) and j in (7, 8) and k in (7, 8)
    assert gravity_field.saddles(np.zeros((16 ** 3, 10)), 16) == []                    # nothing anywhere: no zone has both signs


def test_the_stamp_module_needs_no_compiled_code():
    import subprocess
    import sys
    code = ("import sys\nimport planetgen.api.gravity, planetgen.db.gravity_stamp, planetgen.galaxy.gravity_field\n"
            "assert 'numba' not in sys.modules and 'planetgen.galaxy.gravity_tree' not in sys.modules and 'planetgen.db.gravity_grid' not in sys.modules\n")
    done = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert done.returncode == 0, done.stderr


# --- The route -----------------------------------------------------------------------------------------------------------------------------------

@pytest.fixture
def galaxy(mysql_config, tmp_path, monkeypatch):
    shape = toy_shape(disk_scale_length_pc=3000.0, disk_scale_height_pc=300.0, bulge_scale_radius_pc=400.0, bulge_amplitude=2.0, arm_count=2,
                      pitch_angle_rad=math.radians(15), arm_amplitude=0.4)
    store.save_galaxy_shape(shape, edge_pc=EDGE_PC, outer_shell_index=1600, expected_system_count_at_density_1=2000.0, config=mysql_config,
                            cell_star_ceiling=1e12)
    near = _save(mysql_config, (1500, 2, 700), 2, True, "Near")
    beside = _save(mysql_config, (1500, 2, 701), 1, False, "Beside")
    monkeypatch.setattr(tilecache, "cache_dir", lambda: str(tmp_path))
    return mysql_config, near, beside, tmp_path


@pytest.fixture
def client(galaxy):
    config = galaxy[0]

    class TestConfig(Config):
        MYSQL_CONFIG = config
        WRITE_MYSQL_CONFIG = config
        CONTROL_MYSQL_CONFIG = config

    app = create_app(TestConfig)
    app.testing = True
    return app.test_client()


def _printed(config, row_id):
    conn = store.get_connection(config)
    try:
        return pid("sector", row_id, conn)
    finally:
        conn.close()


def _build(config, sector_id, folder):
    conn = store.get_connection(config)
    try:
        info = store.get_galaxy_shape(conn, light=True)
        pot = potential.for_galaxy(info.shape, info.galaxy_seed)
        return gravity_grid.grid(conn, sector_id, str(folder), pot=pot)
    finally:
        conn.close()


@pytest.fixture
def queued(monkeypatch):
    jobs = []

    def submit(function, *args, expected_seconds=None):
        jobs.append((function, args))
        return f"job{len(jobs)}"
    monkeypatch.setattr(api_jobs, "submit", submit)
    gravity_route = __import__("planetgen.api.routes", fromlist=["_gravity_jobs"])
    gravity_route._gravity_jobs.clear()
    return jobs


def test_a_mode_or_field_that_is_not_one_is_a_400_and_an_unknown_sector_a_404(galaxy, client, queued):
    config, near, *_ = galaxy
    name = _printed(config, near)
    assert client.get(f"/api/sectors/{name}/gravity?mode=speed").status_code == 400
    assert client.get(f"/api/sectors/{name}/gravity?field=both").status_code == 400
    assert client.get("/api/sectors/999999/gravity").status_code == 404
    assert queued == []


def test_with_no_cache_folder_the_route_says_so(galaxy, client, queued, monkeypatch):
    config, near, *_ = galaxy
    monkeypatch.setattr(tilecache, "cache_dir", lambda: None)
    assert client.get(f"/api/sectors/{_printed(config, near)}/gravity").status_code == 503


def test_a_grid_not_worked_out_yet_is_queued_once_and_answered_202(galaxy, client, queued):
    config, near, _beside, folder = galaxy
    url = f"/api/sectors/{_printed(config, near)}/gravity?mode=tidal"
    first = client.get(url)
    assert first.status_code == 202 and first.get_json()["job_id"] == "job1" and first.get_json()["status_url"] == "/api/jobs/job1"
    again = client.get(url)
    assert again.status_code == 202 and again.get_json()["job_id"] == "job1"            # the job is still running: no second one
    assert len(queued) == 1 and queued[0][0] is api_jobs.build_gravity_grid
    assert queued[0][1][0] == near and queued[0][1][2].endswith("gravity")


def test_a_cached_grid_is_answered_in_the_mode_asked(galaxy, client, queued):
    config, near, _beside, folder = galaxy
    built = _build(config, near, folder / gravity_api.SUBDIRECTORY)
    assert not built.hit
    name = _printed(config, near)
    body = client.get(f"/api/sectors/{name}/gravity?mode=tidal").get_json()
    assert body["sector_id"] == name and body["mode"] == "tidal" and body["field"] == "total" and body["unit"] == "1/s^2"
    assert body["zones"] == 16 and len(body["values"]) == 4096 and body["stamp"] == built.stamp
    assert body["edge_pc"] == pytest.approx(EDGE_PC) and body["zone_edge_pc"] == pytest.approx(EDGE_PC / 16)
    assert body["order"] == "x-major" and len(body["origin_pc"]) == 3 and len(body["center_pc"]) == 3
    assert body["minimum"] == min(body["values"]) and body["maximum"] == max(body["values"]) and body["maximum"] > body["minimum"]
    assert all(set(spot) == {"zone", "position_pc"} for spot in body["saddles"])
    assert queued == []
    pull = client.get(f"/api/sectors/{name}/gravity").get_json()
    assert pull["mode"] == "pull" and pull["unit"] == "m/s^2" and pull["minimum"] > 0.0
    depth = client.get(f"/api/sectors/{name}/gravity?mode=depth").get_json()
    assert depth["unit"] == "m^2/s^2" and depth["minimum"] == 0.0 and depth["maximum"] > 0.0


def test_the_total_field_adds_the_smooth_galaxy_to_the_local_one(galaxy, client, queued):
    config, near, _beside, folder = galaxy
    _build(config, near, folder / gravity_api.SUBDIRECTORY)
    name = _printed(config, near)
    total = client.get(f"/api/sectors/{name}/gravity?mode=pull").get_json()
    local = client.get(f"/api/sectors/{name}/gravity?mode=pull&field=local").get_json()
    assert local["field"] == "local"
    # the sector is kiloparsecs out in the toy galaxy: the smooth pull is the bulk of the total, the stars' a small part of it
    assert np.mean(total["values"]) > 20.0 * np.mean(local["values"])


def test_a_grid_made_stale_by_a_new_star_is_queued_again(galaxy, client, queued):
    config, near, beside, folder = galaxy
    _build(config, near, folder / gravity_api.SUBDIRECTORY)
    name = _printed(config, near)
    assert client.get(f"/api/sectors/{name}/gravity").status_code == 200
    _save(config, (1500, 2, 702), 1, False, "Newcomer")
    stale = client.get(f"/api/sectors/{name}/gravity")
    assert stale.status_code == 202 and len(queued) == 1


def test_the_job_stores_the_file_the_route_reads(galaxy, client, queued):
    config, near, _beside, folder = galaxy
    result = api_jobs.build_gravity_grid(near, config, str(folder / gravity_api.SUBDIRECTORY))
    assert result == {"sector_id": near, "stamp": result["stamp"], "built": True}
    assert api_jobs.build_gravity_grid(near, config, str(folder / gravity_api.SUBDIRECTORY))["built"] is False
    assert client.get(f"/api/sectors/{_printed(config, near)}/gravity").status_code == 200


def test_the_directory_of_a_cache_that_is_off_is_none():
    assert gravity_api.directory(None) is None and gravity_api.directory("/x").endswith("gravity")
