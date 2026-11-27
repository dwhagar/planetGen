# tests/test_restored_tilt_and_fallen_moons.py

"""
GEN.1874: a restored planet has the two things the generator set and the stored galaxy lost: whether its tilt wanders chaotically (`oddity_quantities["tilt_chaotic"]`, kept as
`body_history.tilt_chaotic`) and the moons its tide brought inside the Roche limit (`planet.fallen_moons`, the `fallen_moons` table), which its rings are made from.
"""

import pytest

from planetgen.db import store
from planetgen.physics import body_history, ring_origin

from tests.test_restored_birth_numbers import _make, _rocks_and_giants


class _Body:
    def __init__(self, **quantities):
        self.build = type("Build", (), {"orbit_birth": None, "birth_spin": None, "resonance_e": None})()
        self.oddity_quantities = dict(quantities)


def _row(**values):
    row = {name: None for name in body_history.COLUMN_NAMES}
    row["formula_version"] = body_history.FORMULA_VERSION
    row.update(values)
    return row


def test_a_planet_with_a_chaotic_tilt_keeps_one_number_and_one_without_keeps_none():
    assert body_history.numbers_of(_Body(tilt_chaotic=1))["tilt_chaotic"] == 1
    assert "tilt_chaotic" not in body_history.numbers_of(_Body())
    assert "tilt_chaotic" not in body_history.numbers_of(_Body(tilt_chaotic=0))


def test_a_restored_planet_has_the_quantity_the_tilt_oddities_read():
    body = _Body(forced_eccentricity=0.01)
    assert body_history.attach(body, _row(tilt_chaotic=1)) is True
    assert body.oddity_quantities == {"forced_eccentricity": 0.01, "tilt_chaotic": 1}
    other = _Body()
    assert body_history.attach(other, _row()) is False and other.oddity_quantities == {}


def test_the_moons_that_fell_are_written_in_order_and_read_back_whole():
    planet = type("Planet", (), {})()
    planet.fallen_moons = (ring_origin.FallenMoon(3.0e20, 1.5e8), ring_origin.FallenMoon(1.0e19, 4.0e9))
    rows = body_history.fallen_rows(7, planet, 11)
    assert [row["fall_index"] for row in rows] == [0, 1] and all(row["planet_id"] == 11 and row["star_system_id"] == 7 for row in rows)
    restored = type("Planet", (), {})()
    body_history.attach_fallen(restored, list(reversed(rows)))
    assert restored.fallen_moons == planet.fallen_moons


def test_a_planet_with_no_fallen_moons_keeps_what_it_has():
    planet = type("Planet", (), {})()
    assert body_history.fallen_rows(7, planet, 11) == []
    body_history.attach_fallen(planet, [])
    assert not hasattr(planet, "fallen_moons")


@pytest.mark.fresh_database
def test_a_saved_system_gives_back_the_fallen_moons_and_the_chaotic_tilt_of_each_planet(mysql_config):
    for seed in range(60):
        system, cfg = _make(seed)
        if any(getattr(planet, "fallen_moons", None) for planet in _rocks_and_giants(system)):
            break
    else:
        pytest.fail("no system with a fallen moon in 60 seeds")
    chosen = next(planet for planet in _rocks_and_giants(system) if getattr(planet, "fallen_moons", None))
    chosen.oddity_quantities = {**(getattr(chosen, "oddity_quantities", None) or {}), "tilt_chaotic": 1}
    conn = store.get_connection(mysql_config)
    try:
        with conn:
            system_id = store.insert_star_system(conn, system, cfg)
        loaded = store.load_star_system(conn, system_id)
        written = conn.execute("SELECT COUNT(*) AS n FROM fallen_moons WHERE star_system_id = ?", (system_id,)).fetchone()["n"]
    finally:
        conn.close()
    made = {planet.name: planet for planet in _rocks_and_giants(system)}
    assert written == sum(len(getattr(planet, "fallen_moons", None) or ()) for planet in made.values()) and written >= 1
    for again in _rocks_and_giants(loaded):
        original = made[again.name]
        assert tuple(getattr(again, "fallen_moons", None) or ()) == tuple(getattr(original, "fallen_moons", None) or ())
        assert bool((getattr(again, "oddity_quantities", None) or {}).get("tilt_chaotic")) == bool((getattr(original, "oddity_quantities", None) or {}).get("tilt_chaotic"))
