# tests/test_class_outcome.py

"""GEN.1315: a class counts as an outcome only for a body whose stored match is high enough (`physics/class_outcome`), in a request for a class, the count of Class M worlds and the search."""

from types import SimpleNamespace

import pytest

from planetgen import tuning
from planetgen.db import search_values
from planetgen.generation.system import StarSystem
from planetgen.physics import class_outcome


def _body(letter, match):
    return SimpleNamespace(planet_class=letter, class_match=match, moons=[], body_type="t", phi4=None)


def test_a_match_counts_from_the_threshold_up():
    assert class_outcome.fits(1.0) and class_outcome.fits(tuning.CLASS_MATCH_THRESHOLD)
    assert not class_outcome.fits(tuning.CLASS_MATCH_THRESHOLD - 0.01) and not class_outcome.fits(0.0) and not class_outcome.fits(None)


def test_a_body_has_a_class_only_with_a_match_that_counts():
    assert class_outcome.has_class(_body("M", 1.0), "M")
    assert not class_outcome.has_class(_body("M", 0.2), "M")
    assert not class_outcome.has_class(_body("H", 1.0), "M")
    assert not class_outcome.has_class(SimpleNamespace(planet_class="M"), "M")        # a body with no stored match


def test_the_class_m_count_leaves_out_a_poor_fit():
    planets = [_body("M", 1.0), _body("M", 0.3), _body("H", 1.0)]
    system = StarSystem.__new__(StarSystem)
    _, m_count = system.count_habitable(planets)
    assert m_count == 1


def test_the_search_takes_a_range_on_the_class_match():
    key = search_values.KEYS["class_match"]
    assert set(key.entities) == {"planet", "moon"} and key.column == "class_match"
    ranges = search_values.ranges_from_args({"planet_min_class_match": "0.5"})
    clauses, params = [], []
    search_values.clauses("planet", "p", ranges["planet"], clauses, params)
    assert clauses == ["p.class_match >= ?"] and params == [pytest.approx(0.5)]
