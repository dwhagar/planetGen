# tests/test_exosphere_stored.py

"""GEN.455 part B: the exosphere of an airless body (`exosphere.FIELD`) is kept on the body, serialized with it and stored as the rows of `planet_exosphere` and `moon_exosphere`, and a body with any air has none."""

import pytest

from planetgen.db import store
from planetgen.generation.config import SystemConfig
from planetgen.generation.planet import Planet
from planetgen.generation.system import StarSystem
from planetgen.physics import exosphere
from planetgen.util import draw

SEED = 1
"""int: The seed of a system with planets and a moon."""
RECORD = [
    {"species": "Na", "column_cm2": 2.7e10, "source": "photon-stimulated desorption", "scale_height_km": 268.0},
    {"species": "SiO", "column_cm2": 1.0e9, "source": "micrometeoroid impact vaporisation", "scale_height_km": 407.5},
    {"species": "CO2", "column_cm2": 1.0e15, "source": "thermal release", "scale_height_km": 14.0},
]
"""list: A record of three species, each source of the sources named."""


def _config():
    return SystemConfig()


def _bodies(system):
    for planet in system.primary.orbiters:
        if not hasattr(planet, "moons") or getattr(planet, "body_type", None) == "a":        # a belt or a comet is not a body with a ground
            continue
        yield planet
        yield from planet.moons


def _system():
    config = _config()
    draw.set_run_seed(SEED)
    return StarSystem(system_config=config), config


def test_the_field_is_serialized_and_only_an_airless_ground_body_has_one():
    assert exosphere.FIELD in Planet.SERIALIZABLE_FIELDS
    system, _config_used = _system()
    for body in _bodies(system):
        entries = getattr(body, exosphere.FIELD)
        if entries is None:
            continue
        assert exosphere.has_exosphere(body) and (body.atmospheric_pressure or 0.0) == 0.0         # a body with any air never has one
        assert all(set(entry) == {"species", "column_cm2", "source", "scale_height_km"} and entry["column_cm2"] >= exosphere.COLUMN_FLOOR_CM2 for entry in entries)
        assert [entry["column_cm2"] for entry in entries] == sorted((entry["column_cm2"] for entry in entries), reverse=True)


def test_the_record_is_saved_with_the_body_and_loaded_with_it(mysql_config):
    system, config = _system()
    bodies = list(_bodies(system))
    assert any(body.is_moon for body in bodies) and any(not body.is_moon for body in bodies)
    for body in bodies:
        setattr(body, exosphere.FIELD, [dict(entry) for entry in RECORD] if (body.atmospheric_pressure or 0.0) <= 0.0 or body is bodies[0] else None)
    system_id = store.save_system(system, config, config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        loaded = store.load_star_system(conn, system_id)
        rows = conn.execute("SELECT COUNT(*) AS n FROM planet_exosphere").fetchone()["n"] + conn.execute("SELECT COUNT(*) AS n FROM moon_exosphere").fetchone()["n"]
    finally:
        conn.close()
    stored = {(body.mass, round(body.distance, 9)): getattr(body, exosphere.FIELD) for body in _bodies(system)}
    assert rows == sum(len(entries or ()) for entries in stored.values()) and rows >= len(RECORD)
    for body in _bodies(loaded):
        want, got = stored[(body.mass, round(body.distance, 9))], getattr(body, exosphere.FIELD)
        assert (got is None) == (want is None)
        for g, w in zip(got or [], want or []):
            assert g["species"] == w["species"] and g["source"] == w["source"]
            assert g["column_cm2"] == pytest.approx(w["column_cm2"], rel=1e-9) and g["scale_height_km"] == pytest.approx(w["scale_height_km"], rel=1e-9)
