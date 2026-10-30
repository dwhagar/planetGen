"""
The burn-down register (GEN.946, docs/design/research/burn-down-register.md): every stand-in row at tweak level L3 or L4 says what it
needs before it can go, and a line that names a `replace` item names one that is still open.
"""

import pathlib
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import build_stand_in_values as builder  # noqa: E402
import wiring_guard  # noqa: E402

PLAN = (ROOT / "docs" / "design" / "research" / "gap-filling-plan.md").read_text(encoding="utf-8")
REGISTER = (ROOT / "docs" / "design" / "research" / "burn-down-register.md").read_text(encoding="utf-8")
TODO = (ROOT / "docs" / "TODO.md").read_text(encoding="utf-8")
LINES = builder.burn_down_rows(REGISTER)


def hand_tuned_rows(tables):
    """`{section}:{fill_id}` of every registry row at tweak level L3 or L4."""
    found = set()
    for number, (header, rows) in tables.items():
        tweak_at = [name.lower() for name in header].index("tweak")
        found |= {f"{number}:{row[0]}" for row in rows if row[tweak_at] in ("L3", "L4")}
    return found


TUNED = hand_tuned_rows(builder.registry_tables(PLAN))


def test_every_hand_tuned_row_has_a_line():
    assert not TUNED - set(LINES), "rows at L3 or L4 without a line in burn-down-register.md: " + ", ".join(sorted(TUNED - set(LINES)))


def test_a_line_names_a_row_that_is_still_hand_tuned_or_says_it_was_replaced():
    stale = [key for key, (disposition, _, _) in LINES.items() if key not in TUNED and disposition != "replaced"]
    assert not stale, "lines for rows that are no longer at L3 or L4 (delete them, or mark them replaced): " + ", ".join(stale)


def test_a_replace_line_names_an_item_that_is_open():
    opened = wiring_guard.open_items(TODO)
    for key, (disposition, item, _) in LINES.items():
        if disposition == "replace":
            assert item in opened, f"{key}: item {item} is not open in docs/TODO.md"


@pytest.mark.parametrize("key", sorted(LINES))
def test_a_line_has_a_reason_and_only_replace_names_an_item(key):
    disposition, item, reason = LINES[key]
    assert reason, key
    assert (item == "-") is (disposition != "replace"), key


def test_the_register_shows_in_the_last_column_of_the_list():
    listed = builder.render(PLAN, "", REGISTER)
    assert builder.BURN_DOWN_COLUMN in listed
    assert listed.count("| no data |") + listed.count("| replace (") >= len(LINES) - listed.count("replaced")
