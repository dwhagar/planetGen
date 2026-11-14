# tests/test_flavor_scopes.py

"""GEN.1433, GEN.1438: the roots of the flavor tree for the bodies that are not planets of a star: a rogue planet and the moons. Each has a line for a plain body, a body with no star is offered none of the lines
that need a sun, and a rogue planet and its moons draw a line on their own streams."""

from types import SimpleNamespace

import pytest

from planetgen import tuning
from planetgen.generation import life
from planetgen.physics import flavor_tree
from planetgen.util import draw

YEAR_H = 8766.0


def _body(starless=False, is_moon=False, body_type="t", locked=False, **extra):
    star = SimpleNamespace(starless=True, age=3.0, mass=2.0e30, temperature=5000.0, radius=6.0e5) if starless else SimpleNamespace(age=3.0, mass=2.0e30, temperature=5700.0, radius=7.0e5)
    period = 0.01 if is_moon else 1.0
    spin = period * YEAR_H if locked else period * YEAR_H * 0.37
    return SimpleNamespace(body_type=body_type, star=star, is_moon=is_moon, surface_liquid=None, life_class=None, life_stage=None, mass=5.0e24, radius=6400.0, distance=1.0, period=period,
                           rotation_period_hours=spin, atmospheric_pressure=0.0, flavor_key=None, flavor_text_count=0, **extra)


def _ids(body):
    return {o.id for o in flavor_tree.walk("planet", body)}


def test_a_rogue_planet_is_offered_a_line_for_a_plain_body_and_none_that_need_a_sun():
    options = flavor_tree.walk("planet", _body(starless=True))
    assert len(options) >= 3
    assert not any(o.id.startswith(("planet:", "orbital:")) for o in options)     # the ground and orbit lists have beacons, foil sheets and orbits round a sun


def test_a_rogue_gas_giant_and_a_rogue_world_with_water_get_their_own_lines():
    plain = _ids(_body(starless=True))
    giant = _ids(_body(starless=True, body_type="g"))
    assert len(giant - plain) == 2 and len(giant) >= 3
    wet = _ids(_body(starless=True, water_mass_fraction=0.05, has_ground=1))
    assert wet >= plain


def test_a_moon_of_a_star_gets_a_line_for_a_plain_moon_and_the_one_for_its_spin():
    locked, turning = _ids(_body(is_moon=True, locked=True)), _ids(_body(is_moon=True, locked=False))
    assert locked - turning and turning - locked
    assert any(flavor_tree.text_of(i) == "a small world in orbit around a larger one." for i in locked & turning)
    assert not any(flavor_tree.text_of(i).startswith("a small world") for i in _ids(_body(is_moon=False)))


def test_a_world_that_is_not_a_moon_gets_no_moon_line():
    assert not any(i for i in _ids(_body()) if flavor_tree.text_of(i) in ("a host that hangs in one spot in the sky, never rising and never setting.", "a host that rises and sets as the moon turns."))


def test_a_rogue_body_draws_its_own_line_and_the_system_bookkeeping_is_not_touched(monkeypatch):
    monkeypatch.setattr(tuning, "FLAVOR_CHANCE_PLANET", 1.0)
    body = _body(starless=True)
    life.decide_flavor_text_alone(body)
    first = body.flavor_key
    assert first and first.startswith(("planet:", "orbital:", "life:", "habitable:")) and body.flavor_text_count == 1
    other = _body(starless=True)
    life.decide_flavor_text_alone(other)
    assert other.flavor_key == first                      # the same body, the same line, whatever else was drawn before it


def test_a_roll_that_fails_leaves_the_body_without_a_line(monkeypatch):
    monkeypatch.setattr(tuning, "FLAVOR_CHANCE_PLANET", 0.0)
    body = _body(starless=True)
    life.decide_flavor_text_alone(body)
    assert body.flavor_key is None and body.flavor_text_count == 0


def _belt(**extra):
    values = dict(mass_kg=3.0e21, lower_limit=2.2, upper_limit=3.6, largest_body_km=400.0, stirred=False, distance=2.7, flavor_key=None, flavor_text_count=0,
                  star=SimpleNamespace(age=3.0, mass=2.0e30))
    values.update(extra)
    return SimpleNamespace(**values)


def test_a_plain_belt_is_offered_the_three_lines_for_a_plain_belt_and_the_gated_ones_add_to_them():
    plain = {o.id for o in flavor_tree.walk("belt", _belt())}
    assert len(plain) >= 3 and all(i.startswith("belt:") for i in plain)
    assert {o.id for o in flavor_tree.walk("belt", _belt(stirred=True))} - plain
    assert {o.id for o in flavor_tree.walk("belt", _belt(largest_body_km=900.0))} - plain
    assert {o.id for o in flavor_tree.walk("belt", _belt(mass_kg=1.0e23))} - plain
    assert {o.id for o in flavor_tree.walk("belt", _belt(upper_limit=2.8))} - plain


def test_a_belt_draws_its_own_line(monkeypatch):
    monkeypatch.setattr(tuning, "FLAVOR_CHANCE_PLANET", 1.0)
    belt = _belt()
    life.decide_flavor_text_alone(belt, "belt")
    assert belt.flavor_key and belt.flavor_key.startswith("belt:") and flavor_tree.text_of(belt.flavor_key)
    other = _belt()
    life.decide_flavor_text_alone(other, "belt")
    assert other.flavor_key == belt.flavor_key
