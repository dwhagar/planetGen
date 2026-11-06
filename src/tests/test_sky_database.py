"""VIEW.23 with VIEW.22: the sky of an observer from the stars a real database stores (`galaxy.sky.sky_from_database`)."""

import numpy as np
import pytest

from planetgen.db import store
from planetgen.galaxy import sky
from tests.test_sky_sphere import add_bright, add_system


@pytest.fixture
def conn(mysql_config):
    connection = store.get_connection(mysql_config)
    try:
        yield connection
    finally:
        connection.close()


def flat(radius):
    return np.full(np.shape(radius), 220.0)


def test_a_charted_star_and_a_bright_star_appear_in_the_sky_and_the_observers_own_system_does_not(conn):
    observer_point = (8128.0, 0.0, 0.0)
    own_id, _ = add_system(conn, observer_point, [(1.0, 5772.0)], name="Own")
    near_id, near_stars = add_system(conn, (8138.0, 0.0, 0.0), [(1.0, 5772.0)], velocity=(0.0, 0.0, 0.0), name="Near")
    bright_id = add_bright(conn, (8128.0, 60.0, 0.0), 2000.0)
    conn.commit()
    observer = sky.Observer(observer_point, 1.0e9, (0.0, 0.0, 0.0), own_id)
    result = sky.sky_from_database(conn, observer, 1.0e9, 7.5, circular_speed_kms=flat, dust={"av_per_kpc_at_sun": 0.0}, max_distance_pc=200.0)
    assert f"0:{near_stars[0]}" in list(result.uid)
    assert f"1:{bright_id}" in list(result.uid)
    assert not any(uid.startswith("0:") and uid != f"0:{near_stars[0]}" for uid in result.uid)       # the own system is not in its own sky
    assert result.report["tiers"]["charted"]["visible"] == 1 and result.report["tiers"]["bright_star"]["visible"] == 1
    assert "built_to_pc" in result.report


def test_a_galaxy_with_a_shape_is_seen_through_the_dust_of_its_own_gas_and_one_without_through_the_fixed_law(conn, monkeypatch):
    from types import SimpleNamespace
    from planetgen.galaxy import shape_draw
    from planetgen.physics import dust_model
    shape = shape_draw.default_shape(bytes.fromhex("5e7a1c3b9d2f4608a1b2c3d4e5f60718"))
    monkeypatch.setattr(store, "get_galaxy_shape", lambda connection, light=False: SimpleNamespace(shape=shape))
    given = sky.default_dust(conn)
    assert set(given) == {"profile"} and given["profile"] == dust_model.profile_for_shape(shape)
    monkeypatch.setattr(store, "get_galaxy_shape", lambda connection, light=False: None)
    assert sky.default_dust(conn) == {}


def test_the_sky_of_a_database_takes_the_default_dust_unless_the_caller_gives_one(conn, monkeypatch):
    observer_point = (8128.0, 0.0, 0.0)
    own_id, _ = add_system(conn, observer_point, [(1.0, 5772.0)], name="Own")
    bright_id = add_bright(conn, (8128.0, 60.0, 0.0), 2000.0)
    conn.commit()
    observer = sky.Observer(observer_point, 1.0e9, (0.0, 0.0, 0.0), own_id)
    monkeypatch.setattr(sky, "default_dust", lambda connection: {"av_per_kpc_at_sun": 1000.0})
    dim = sky.sky_from_database(conn, observer, 1.0e9, 7.5, circular_speed_kms=flat, max_distance_pc=200.0)
    clear = sky.sky_from_database(conn, observer, 1.0e9, 7.5, circular_speed_kms=flat, dust={"av_per_kpc_at_sun": 0.0}, max_distance_pc=200.0)
    assert f"1:{bright_id}" in list(clear.uid) and f"1:{bright_id}" not in list(dim.uid)
