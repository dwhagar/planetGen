# tests/test_deck_rows.py

"""GEN.1577: the sheet rows of a giant's stored cloud decks and trace-gas ceilings (`describe/deck_rows.py`)."""

import pytest

from planetgen.describe import deck_rows
from planetgen.physics import hot_cloud_deck as hcd
from planetgen.physics import hot_condensates as hc
from planetgen.physics import stellar_abundances as sa


class _Conn:
    def __init__(self, text):
        self.text = text

    def execute(self, sql, params):
        assert "cloud_decks" in sql and "planets" in sql
        return self

    def fetchone(self):
        return {"cloud_decks": self.text}


def _text():
    host = sa.composition(0.0, {})
    return hcd.pack(hcd.deck(host, lambda p: 1000.0 * p ** 0.15), hcd.trace_gas_ceilings(host))


def test_a_giant_gets_a_row_for_each_cloud_and_each_gas():
    found = deck_rows.rows(_Conn(_text()), "planet", 1)
    keys = [row[1] for row in found]
    assert keys[0].startswith(deck_rows.CLOUD_KEY_PREFIX) and len(keys) == len(set(keys))
    assert sum(key.startswith(deck_rows.CEILING_KEY_PREFIX) for key in keys) == len(hc.HOT_TRACE_GASES)
    assert all(row[0] == "atmosphere" and deck_rows.is_derived_key(row[1]) for row in found)


def test_a_cloud_deeper_than_the_grid_says_so():
    text = hcd.pack(hcd.deck(sa.composition(0.0, {}), lambda p: 100.0), {})
    row = deck_rows.rows(_Conn(text), "planet", 1)[0]
    assert row[3][0] == "deeper than " and row[4] == pytest.approx(hcd.GRID_BAR[-1], rel=1e-3)


def test_no_stored_deck_and_other_kinds_give_no_rows():
    assert deck_rows.rows(_Conn(None), "planet", 1) == []
    assert deck_rows.rows(_Conn(_text()), "moon", 1) == []
    assert not deck_rows.is_derived_key("surface_pressure")
