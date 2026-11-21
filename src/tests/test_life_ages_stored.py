# tests/test_life_ages_stored.py

"""GEN.1746 slice 3: the life ages of a world are made at the end of the life pass, kept on the body (`life_ages.FIELD`) and stored as the rows of `planet_life_ages` and `moon_life_ages`."""

import pytest

from planetgen.db import store
from planetgen.generation import life
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.physics import air_record, life_ages
from tests.fuzz_support import deterministic_entropy

SEEDS = range(7500, 7540)


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


def _systems():
    for seed in SEEDS:
        config = _config()
        with deterministic_entropy(seed):
            yield StarSystem(system_config=config), config


def _first(condition):
    for system, _config_used in _systems():
        for body in _bodies(system):
            if condition(body):
                return body
    pytest.fail("no generated body of the kind the test needs")


def test_the_entries_and_the_ages_go_into_each_other_and_none_is_none():
    ages = (life_ages.LifeAge(1, 0.5, 1.2, "impact", 0.4), life_ages.LifeAge(2, 1.21, None, None, None))
    entries = life_ages.entries_of(ages)
    assert entries[0] == {"number": 1, "start_gyr": 0.5, "end_gyr": 1.2, "cause": "impact", "kill_share": 0.4}
    assert life_ages.from_entries(entries) == ages
    assert life_ages.entries_of(()) is None and life_ages.from_entries(None) == ()


def test_a_world_with_life_keeps_the_ages_its_life_pass_made_and_the_record_agrees():
    body = _first(lambda b: b.life_class is not None and getattr(b, life_ages.FIELD, None))
    ages = life_ages.from_entries(getattr(body, life_ages.FIELD))
    assert ages[0].number == 1 and ages[-1].end_gyr is None and ages[-1].kill_share is None
    begin = min(one.started(body) for one in body.build.air_living)
    assert ages == tuple(life_ages.ages_of(body, begin))
    for before, after in zip(ages, ages[1:]):
        assert before.end_gyr is not None and before.kill_share >= life_ages.AGE_KILL_THRESHOLD and before.cause in ("flood basalt", "impact")
        assert after.start_gyr >= before.end_gyr
    if air_record.eligible(body):
        assert air_record.unpack(body.air_history).ages == ages


def test_a_world_with_no_life_has_no_ages():
    body = _first(lambda b: b.life_class is None and b.body_type != "g")
    assert getattr(body, life_ages.FIELD) is None


def test_an_edit_keeps_the_ages_the_body_has():
    body = _first(lambda b: b.life_class is not None and getattr(b, life_ages.FIELD, None))
    before = getattr(body, life_ages.FIELD)
    life.apply_life_data(body, keep_stored=True)
    assert getattr(body, life_ages.FIELD) == before


def test_the_ages_are_saved_with_the_body_and_loaded_with_it(mysql_config):
    for made, config in _systems():
        if any(getattr(b, life_ages.FIELD, None) and len(getattr(b, life_ages.FIELD)) > 1 for b in _bodies(made)):
            break
    else:
        pytest.fail("no system with a world that has more than one life age in these seeds")
    system_id = store.save_system(made, config, config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        loaded = store.load_star_system(conn, system_id)
        rows = conn.execute("SELECT COUNT(*) AS n FROM planet_life_ages").fetchone()["n"] + conn.execute("SELECT COUNT(*) AS n FROM moon_life_ages").fetchone()["n"]
    finally:
        conn.close()
    stored = {(b.mass, round(b.distance, 9)): getattr(b, life_ages.FIELD) for b in _bodies(made)}
    assert rows == sum(len(entries or ()) for entries in stored.values()) and rows > 1
    for body in _bodies(loaded):
        want, got = stored[(body.mass, round(body.distance, 9))], getattr(body, life_ages.FIELD)
        assert (got is None) == (want is None)
        for g, w in zip(got or [], want or []):
            assert g["number"] == w["number"] and g["cause"] == w["cause"]
            assert g["start_gyr"] == pytest.approx(w["start_gyr"], rel=1e-9)
            assert (g["end_gyr"] is None) == (w["end_gyr"] is None) and (w["end_gyr"] is None or g["end_gyr"] == pytest.approx(w["end_gyr"], rel=1e-9))
