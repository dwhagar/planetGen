# tests/test_day_rows.py

"""
GEN.1153: a star's sheet says `no night`, and a free brown dwarf's `no solid surface` (with no solar day), as derived rows; a rogue planet the planet chain built keeps the state it stores.
"""

import pytest

from planetgen.describe import clock_rows, day_rows, model
from planetgen.physics import day_night


class _Conn:
    def __init__(self, row):
        self.row = row

    def execute(self, sql, params=()):
        return self

    def fetchone(self):
        return self.row


def _text(row):
    return model.inline_text(row[3]) if hasattr(model, "inline_text") else row[3][0]


def test_a_stars_sheet_says_it_has_no_night():
    (group, key, label, value, raw), = day_rows.rows(_Conn(None), "star", 1)
    assert (group, key, label) == ("rotation", day_rows.STAR_KEY, "State of its day and night")
    assert value == [day_night.NO_NIGHT] and raw is None


def test_a_free_brown_dwarf_with_no_stored_state_has_no_solid_surface_and_no_sun_to_make_a_day():
    (group, key, label, value, raw), = day_rows.rows(_Conn({"day_night_state": None, "planet_type": "g"}), "rogue_planet", 3)
    assert key == day_rows.FREE_KEY and value == [day_night.NO_SURFACE]


def test_a_rogue_planet_that_stores_its_state_keeps_it_and_other_kinds_get_nothing():
    assert day_rows.rows(_Conn({"day_night_state": "no sun", "planet_type": "t"}), "rogue_planet", 3) == []
    for kind in ("planet", "moon", "belt", "comet", "nebula"):
        assert day_rows.rows(_Conn(None), kind, 1) == []


def test_the_keys_are_derived_keys_of_the_sheet():
    assert all(clock_rows.is_derived_key(key) and day_rows.is_derived_key(key) for key in day_rows.KEYS)
    assert not clock_rows.is_derived_key("day_night_state")
