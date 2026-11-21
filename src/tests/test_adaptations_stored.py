# tests/test_adaptations_stored.py

"""GEN.1809: the adaptations the life of a world needed (`adaptation_record`, GEN.1808) are kept on the body (`adaptation_record.FIELD`) and stored as the rows of `planet_adaptations` and `moon_adaptations`."""

import pytest

from planetgen.db import store
from planetgen.generation import life
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import adaptation_record
from tests.fuzz_support import deterministic_entropy

SEEDS = (7533,)
"""tuple: The seed of a system (with `HABITABLE_WORLD`) that holds a world with a record of adaptations (Irachas I: hours of liquid water at the extreme grade and two habitat rows)."""
RECORD = [
    {"axis": "t_cold", "habitat": "brine", "demand": 258.25, "grade": 3, "rows": ["C8"]},
    {"axis": "water_hours", "habitat": "melt_film", "demand": 4707.04, "grade": None, "rows": []},
    {"axis": None, "habitat": "melt_film", "demand": None, "grade": None, "rows": ["C6"]},
    {"axis": "light_floor", "habitat": "water_under_ice", "demand": 0.3, "grade": 2, "rows": ["C1", "C2"]},
]
"""list: A record with every kind of row: a grade, past the table (no grade, no rows), a habitat row and two rows that give one need."""


def _config():
    config = SystemConfig()
    config.HABITABLE_WORLD = True
    return config


def _bodies(system):
    for planet in system.planets:
        if getattr(planet, "body_type", None) == "a":           # a belt is not a body with an air
            continue
        yield planet
        yield from getattr(planet, "moons", ())


def _system_with_life():
    for seed in SEEDS:
        config = _config()
        with deterministic_entropy(seed):
            system = StarSystem(system_config=config)
        if any(body.life_class is not None for body in _bodies(system)):
            return system, config
    pytest.fail("no generated system with a world that has life")


def test_a_body_with_no_record_has_none_and_the_attribute_is_serialized():
    from planetgen.generation.planet import Planet
    assert adaptation_record.FIELD in Planet.SERIALIZABLE_FIELDS
    system, _config_used = _system_with_life()
    for body in _bodies(system):
        entries = getattr(body, adaptation_record.FIELD)
        assert entries is None or all(set(entry) == {"axis", "habitat", "demand", "grade", "rows"} for entry in entries)


def test_an_edit_keeps_the_record_the_body_has():
    system, _config_used = _system_with_life()
    body = next(body for body in _bodies(system) if body.life_class is not None)
    setattr(body, adaptation_record.FIELD, [dict(entry) for entry in RECORD])
    life.apply_life_data(body, keep_stored=True)
    assert getattr(body, adaptation_record.FIELD) == RECORD


def test_the_record_is_saved_with_the_body_and_loaded_with_it(mysql_config):
    system, config = _system_with_life()
    recorded = [body for body in _bodies(system) if getattr(body, adaptation_record.FIELD, None)]
    assert recorded, "the seed no longer makes a world with a record: find another"
    extra = next(body for body in _bodies(system) if body is not recorded[0] and body.life_class is not None) if len([b for b in _bodies(system) if b.life_class is not None]) > 1 else None
    if extra is not None and not getattr(extra, adaptation_record.FIELD, None):
        setattr(extra, adaptation_record.FIELD, [dict(entry) for entry in RECORD])       # every kind of row, with none past the table, and the empty row list
    system_id = store.save_system(system, config, config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        loaded = store.load_star_system(conn, system_id)
        rows = conn.execute("SELECT COUNT(*) AS n FROM planet_adaptations").fetchone()["n"] + conn.execute("SELECT COUNT(*) AS n FROM moon_adaptations").fetchone()["n"]
    finally:
        conn.close()
    stored = {(body.mass, round(body.distance, 9)): getattr(body, adaptation_record.FIELD) for body in _bodies(system)}
    assert rows == sum(len(entries or ()) for entries in stored.values()) and rows >= len(getattr(recorded[0], adaptation_record.FIELD))
    for body in _bodies(loaded):
        want, got = stored[(body.mass, round(body.distance, 9))], getattr(body, adaptation_record.FIELD)
        assert (got is None) == (want is None)
        for g, w in zip(got or [], want or []):
            assert g["axis"] == w["axis"] and g["habitat"] == w["habitat"] and g["grade"] == w["grade"] and g["rows"] == w["rows"]
            assert (g["demand"] is None) == (w["demand"] is None) and (w["demand"] is None or g["demand"] == pytest.approx(w["demand"], rel=1e-9))
