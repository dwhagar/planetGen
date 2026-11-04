# tests/test_belt_history_page.py

"""
GEN.1232: what the star's evolution did to a belt (`physics/belt_history.py`, run by `remnant_orbits`) is stored with the belt
(`asteroid_belts.history_text`) and shown on its page in plain words. A belt nothing happened to stores NULL and shows no history.
"""

import types

from planetgen.db import object_rows, store
from planetgen.describe import build, sheet
from planetgen.generation import remnant_orbits
from planetgen.generation.belt import AsteroidBelt
from planetgen.physics import belt_history

LINES = ("A giant engulfed its inner 0.4 au.", "The star lost mass and a share of the belt was thrown inward.")


def _belt(history):
    belt = types.SimpleNamespace(distance=2.0, lower_limit=1.5, upper_limit=3.0, density="typical", mass_kg=1.5e21, chondrite_group="CM",
                                 get_composition_summary=lambda: "rock", history=history)
    return belt


def _stored(belt):
    columns = object_rows.columns_of(store.BELT_FIELDS)
    values = object_rows.values_of(store.BELT_FIELDS, belt, {"star_system_id": 1, "star_id": None, "orbital_index": 0})
    return dict(zip(columns, values))


def test_the_history_is_one_paragraph_in_the_column_and_null_when_nothing_happened():
    assert _stored(_belt(LINES))["history_text"] == " ".join(LINES)
    assert _stored(_belt(()))["history_text"] is None


def test_the_column_reads_back_as_the_history_of_the_belt():
    row = {"history_text": " ".join(LINES)}
    field = next(field for field in store.BELT_FIELDS if field.columns == ("history_text",))
    assert field.load(row) == {"history": [" ".join(LINES)]}
    assert field.load({"history_text": None}) == {"history": []}


def test_a_belt_keeps_its_history_through_its_dict_and_an_untouched_belt_has_none():
    assert AsteroidBelt.history == ()
    belt = object.__new__(AsteroidBelt)
    assert "history" in AsteroidBelt.SERIALIZABLE_FIELDS and belt.history == ()


def test_the_events_of_the_giant_branch_make_the_lines_the_belt_stores():
    belt = types.SimpleNamespace(lower_limit=1.0, upper_limit=3.0, distance=2.0, mass_kg=1.0e21)
    remnant_orbits.engulf_belt(belt, 1.5)
    assert belt.history and belt_history.history_text(belt) == " ".join(belt.history)
    assert _stored(_belt(belt.history))["history_text"] == " ".join(belt.history)


def test_the_belt_page_shows_the_history_after_its_composition_and_nothing_for_a_quiet_belt():
    def section(history):
        belt = _belt(history)
        belt.composition = [("iron", "high")]
        belt.lower_limit, belt.upper_limit = 1.5, 3.0
        return build._belt_section(belt)

    shown = [block for block in section(LINES)["blocks"] if block.get("role") == "history"]
    assert len(shown) == 1 and shown[0]["type"] == "prose"
    assert not [block for block in section(())["blocks"] if block.get("role") == "history"]


def test_the_property_sheet_lists_the_history_row_of_a_belt():
    assert "history_text" in [spec.key for spec in sheet.SPECS["belt"]]
