# tests/test_migration_aerosols.py

"""
GEN.414, DB.68: the haze, dust and cloud depths of a planet or moon. A system saved and loaded keeps them; its draws do not move any other
draw of the run (the same system generated with the aerosols switched off differs in these five fields alone); migration 0138 gives a
galaxy made at v137 five NULL columns on `planets` and `moons`, which read as no aerosol.
"""

import pytest

from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import aerosols
from planetgen.util import draw
from tests.db_schema_support import load_old_schema

pytestmark = pytest.mark.fresh_database


def _config():
    cfg = SystemConfig()
    cfg.STAR_TYPE = "G2V"
    cfg.PLANETS = True
    cfg.MAX_PLANETS = True
    cfg.BINARY_SYSTEM = False
    return cfg


def _bodies(system):
    return [b for p in system.planets if p.body_type != 'a' for b in [p, *p.moons]]


def _system(seed=None):
    for number in range(60) if seed is None else [seed]:
        cfg = _config()
        with draw.bound(number):
            system = StarSystem(system_config=cfg)
        drawn = [b for b in _bodies(system) if b.tau_cloud is not None]
        if len(drawn) >= 2 and any(b.tau_cloud > 0 for b in drawn):
            return system, cfg, number
    pytest.fail("no system with two bodies that have air and a cloud in 60 seeds")


READ_FROM_AEROSOLS = ("uv_surface_index", "phi4_radiation", "phi4", "phi4_suitability", "tier_radiation", "eq_body", "eq_body_kind", "atmosphere_colour", "eq_dust")
"""tuple: What the body stores that reads the aerosols (GEN.445: the cloud, dust and haze cut the ground UV, and the Radiation score is the lower of the dose's score and the ground UV's; the Suitability score, the overall score `phi4`, the radiation tier and the body-protection equipment level and kind are worked from the Radiation score, so a ground UV that moves a dose score moves them too, with no draw; GEN.418: the sky colour is seen through them; the dust equipment level follows the dust depth), so it moves with them without a draw."""


def _fields(body):
    """What the body stores except the aerosol values and what is derived from them, moons included."""
    return {name: ([_stripped(moon) for moon in value] if name == "moons" else value)
            for name, value in body.to_dict().items() if name not in aerosols.BODY_FIELDS + READ_FROM_AEROSOLS}


def _stripped(record):
    return {name: value for name, value in record.items() if name not in aerosols.BODY_FIELDS + READ_FROM_AEROSOLS}


def test_a_body_with_air_has_every_value_and_one_without_has_none():
    system, _, _ = _system()
    for body in _bodies(system):
        values = [getattr(body, field) for field in aerosols.BODY_FIELDS]
        if (body.atmospheric_pressure or 0.0) > 0.0 and body.surface_temperature:
            assert all(value is not None for value in values), body.name
            assert 0.0 <= body.cloud_cover <= 1.0 and 0.0 <= body.aerosol_draw < 1.0
            assert min(body.tau_haze, body.tau_dust, body.tau_cloud) >= 0.0
        else:
            assert all(value is None for value in values), body.name


def test_the_same_seed_draws_the_same_aerosols_and_they_move_no_other_draw(monkeypatch):
    system, _, seed = _system()
    again, _, _ = _system(seed)
    assert [aerosols.of_body(b) for b in _bodies(system)] == [aerosols.of_body(b) for b in _bodies(again)]
    monkeypatch.setattr(aerosols, "generate_aerosols", lambda body: [setattr(body, f, None) for f in aerosols.BODY_FIELDS])
    with draw.bound(seed):
        without = StarSystem(system_config=_config())
    assert [_fields(b) for b in _bodies(without)] == [_fields(b) for b in _bodies(system)]


def test_a_saved_system_keeps_its_aerosols(mysql_config):
    system, cfg, _ = _system()
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        loaded = store.load_star_system(conn, system_id)
        assert [aerosols.of_body(b) for b in _bodies(loaded)] == [aerosols.of_body(b) for b in _bodies(system)]
        rows = conn.execute("SELECT tau_haze, tau_dust, tau_cloud, cloud_cover, aerosol_draw FROM planets WHERE star_system_id = ?", (system_id,)).fetchall()
        assert any(row["tau_cloud"] is not None for row in rows)
    finally:
        conn.rollback()
        conn.close()


def test_a_v137_galaxy_keeps_its_bodies_with_no_aerosol(mysql_config):
    load_old_schema(mysql_config, 137)
    assert store.schema_status(mysql_config)[0] == 137
    assert store.migrate_database(mysql_config) == store.SCHEMA_VERSION
    conn = store.get_connection(mysql_config, ensure_schema=False)
    try:
        for table in ("planets", "moons"):
            found = conn.execute("SELECT COLUMN_NAME, IS_NULLABLE, DATA_TYPE FROM information_schema.COLUMNS WHERE TABLE_SCHEMA = DATABASE()"
                                 " AND TABLE_NAME = ? AND COLUMN_NAME IN ('tau_haze', 'tau_dust', 'tau_cloud', 'cloud_cover', 'aerosol_draw')", (table,)).fetchall()
            assert sorted(row["COLUMN_NAME"] for row in found) == sorted(aerosols.BODY_FIELDS)
            assert all(row["IS_NULLABLE"] == "YES" and row["DATA_TYPE"] == "double" for row in found)
        assert store._stored_version(conn, "schema_migrations") == store.SCHEMA_VERSION
    finally:
        conn.close()
    assert aerosols.of_body(type("Old", (), {})()) == aerosols.NONE
