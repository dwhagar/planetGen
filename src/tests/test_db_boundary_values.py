# tests/test_db_boundary_values.py

"""
TEST.12: values at the edges of what the columns hold. DOUBLE extremes
round-trip exactly; NaN and infinity are refused before anything is
written; a system name at the length limit saves with every planet and
moon named after it, and one past it is refused cleanly by `_db`, the
API and the CLI instead of failing mid-save with MySQL error 1406;
4-byte UTF-8 names round-trip; and the system config's three-way
(True/False/None) flags come back as they went in.
"""

import itertools
import math
import subprocess
import sys

import pymysql
import pytest

from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from tests.fuzz_support import deterministic_entropy
from tests.test_api import admin_client, client, default_admin_client, first_admin_password  # noqa: F401  (fixtures)
from tests.publicids import pid, pids

pytestmark = pytest.mark.db

_seeds = itertools.count(1201)


def _system(name=None, **flags):
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.BINARY_SYSTEM = False
    cfg.NAME = name
    for flag, value in flags.items():
        setattr(cfg, flag, value)
    with deterministic_entropy(next(_seeds)):
        return StarSystem(system_config=cfg), cfg


def _load(config, system_id):
    conn = store.get_connection(config)
    try:
        return store.load_star_system(conn, system_id)
    finally:
        conn.close()


def _count(config, table):
    conn = store.get_connection(config)
    try:
        return conn.execute(f"SELECT COUNT(*) AS n FROM {table}").fetchone()["n"]
    finally:
        conn.close()


@pytest.mark.parametrize("value", [
    1.7976931348623157e308, -1.7976931348623157e308, 5e-324, 2.2250738585072014e-308, -0.0, 1e-300])
def test_double_extremes_round_trip(mysql_config, value):
    system, cfg = _system(PLANETS=False)
    system.star.mass = value
    system_id = store.save_system(system, cfg, config=mysql_config)
    loaded = _load(mysql_config, system_id).star.mass
    assert loaded == value
    if value != 0:
        # Zero's sign carries no meaning in a stored value, and MariaDB
        # 11 reads a stored -0.0 back as +0.0 (DB.12); every other value
        # keeps its sign.
        assert math.copysign(1, loaded) == math.copysign(1, value)


@pytest.mark.parametrize("value", [math.nan, math.inf, -math.inf])
def test_nan_and_infinity_are_refused_with_nothing_written(mysql_config, value):
    system, cfg = _system(PLANETS=False)
    name = system.name
    system.star.mass = value
    with pytest.raises(pymysql.err.ProgrammingError):
        store.save_system(system, cfg, config=mysql_config)
    assert system.name == name
    assert (_count(mysql_config, "star_systems"), _count(mysql_config, "stars")) == (0, 0)
    assert _count(mysql_config, "system_name_registry") == 0


def test_longest_system_name_saves_with_every_body_named_after_it(mysql_config):
    name = "N" * store.SYSTEM_NAME_MAX_LENGTH
    ids = []
    for _ in range(2):  # the second collides, so both gain a Greek-letter prefix
        system, cfg = _system(name, PLANETS=True, MAX_PLANETS=True, MOONS=True)
        assert system.planets, "the test needs planets named after the system"
        ids.append(store.save_system(system, cfg, config=mysql_config))
    loaded = [_load(mysql_config, system_id) for system_id in ids]
    assert sorted(system.name for system in loaded) == [f"Alpha {name}", f"Beta {name}"]
    assert all(planet.name.startswith(loaded[0].name) for planet in loaded[0].planets if hasattr(planet, "name"))


def test_too_long_system_name_is_refused_before_saving(mysql_config):
    system, cfg = _system("N" * (store.SYSTEM_NAME_MAX_LENGTH + 1), PLANETS=False)
    with pytest.raises(ValueError, match="at most"):
        store.save_system(system, cfg, config=mysql_config)
    assert _count(mysql_config, "star_systems") == 0


def test_api_refuses_a_too_long_system_name(admin_client):
    too_long = "N" * (store.SYSTEM_NAME_MAX_LENGTH + 1)
    assert admin_client.post("/api/systems", json={"planets": False, "name": too_long}).status_code == 400
    response = admin_client.post("/api/systems", json={"planets": True, "moons": True})
    system_id = response.get_json()["id"]
    assert admin_client.patch(f"/api/systems/{pid('system', system_id)}", json={"name": too_long}).status_code == 400
    longest = "N" * store.SYSTEM_NAME_MAX_LENGTH
    response = admin_client.patch(f"/api/systems/{pid('system', system_id)}", json={"name": longest})
    assert response.status_code == 200, response.get_json()
    assert admin_client.get(f"/api/systems/{pid('system', system_id)}").get_json()["name"] == longest
    star_id = admin_client.get(f"/api/systems/{pid('system', system_id)}").get_json()["stars"][0]["id"]
    assert admin_client.patch(f"/api/stars/{star_id}", json={"name": too_long}).status_code == 400


def test_cli_refuses_a_too_long_name(tmp_path):
    result = subprocess.run(
        [sys.executable, "-m", "planetgen.cli.generate", "system", "--name", "N" * (store.SYSTEM_NAME_MAX_LENGTH + 1)],
        cwd=str(tmp_path), capture_output=True, text=True, timeout=120)
    assert result.returncode == 2
    assert "--name" in result.stderr and "Traceback" not in result.stderr


@pytest.mark.parametrize("name", ["Vega 🌌", "𝔙𝔢𝔤𝔞", "Ṽéga 星"])
def test_four_byte_names_round_trip(mysql_config, admin_client, name):
    system, cfg = _system(name, PLANETS=True)
    system_id = store.save_system(system, cfg, config=mysql_config)
    loaded = _load(mysql_config, system_id)
    assert loaded.name == name
    assert all(planet.name.startswith(name) for planet in loaded.planets if hasattr(planet, "name"))

    renamed = f"{name} ⭐"
    response = admin_client.patch(f"/api/systems/{pid('system', system_id)}", json={"name": renamed})
    assert response.status_code == 200, response.get_json()
    assert admin_client.get(f"/api/systems/{pid('system', system_id)}").get_json()["name"] == renamed
    assert _load(mysql_config, system_id).name == renamed


@pytest.mark.parametrize("value", [True, False, None])
def test_system_config_tristates_round_trip(mysql_config, value):
    cfg = SystemConfig()
    for flag in ("HABITABLE_WORLD", "ASTEROID_BELT", "LARGE_STAR", "MOONS", "MAX_PLANETS", "PLANETS", "COMETS"):
        setattr(cfg, flag, value)
    conn = store.get_connection(mysql_config)
    try:
        config_id = store.insert_system_config(conn, cfg)
        conn.commit()
        loaded = store.load_system_config(conn, config_id)
    finally:
        conn.close()
    for flag in ("HABITABLE_WORLD", "ASTEROID_BELT", "LARGE_STAR", "MOONS", "MAX_PLANETS", "PLANETS", "COMETS"):
        assert getattr(loaded, flag) is value, flag
