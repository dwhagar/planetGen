"""
GEN.1002: every registered row that rests on a Solar System body has a row in the stand-in registry (docs/design/research/gap-filling-plan.md
section 6.56) that names it as `table:row` and the item that replaces it, and the registered row carries the same item.

`test_calibration_fence.py` finds the rows (TEST.215) and `tests/data/solar_system_rows.json` lists them: the rows flagged as stand-ins and the
rows still pending the flag. The state of a pending row changes with the digest of its table, which rides a merge train (GEN.1002), so the
registry rows exist first and this test holds from then on.
"""

import json
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import build_stand_in_values as builder  # noqa: E402
from tests.audit import sweep_tables  # noqa: E402

SECTION = "6.56"
DATA = json.loads((pathlib.Path(__file__).resolve().parent / "data" / "solar_system_rows.json").read_text(encoding="utf-8"))
TABLES = builder.registry_tables((ROOT / "docs" / "design" / "research" / "gap-filling-plan.md").read_text(encoding="utf-8"))


def _registry():
    """`{(table, row): (fill_id, replaced_by)}` of section 6.56, from the `constant` cell that opens with `table:row`."""
    header, rows = TABLES[SECTION]
    found = {}
    for row in rows:
        cells = dict(zip(header, row))
        opening = re.match(r"([a-z0-9-]+):([a-z0-9_]+)\b", cells["constant"])
        assert opening, f"{cells['fill_id']}: the constant does not open with table:row"
        key = (opening.group(1), opening.group(2))
        assert key not in found, f"{key} has two rows"
        found[key] = (cells["fill_id"], cells["replaced_by"].strip())
    return found


def _listed(*sections):
    return {(table, row) for section in sections for table, rows in DATA[section].items() for row in rows}


def test_every_pending_row_has_a_registry_row_and_every_registry_row_is_listed():
    """The rows flagged earlier have rows in the sections of their own items; this section holds the rows of GEN.1002, pending or flagged by it."""
    registry = set(_registry())
    assert _listed("pending_decision") - registry == set(), f"pending rows that have no row in section {SECTION}"
    assert registry - _listed("pending_decision", "stand_ins") == set(), f"section {SECTION} rows that tests/data/solar_system_rows.json does not list"


def test_a_registry_row_names_a_replacement_item_and_the_registered_row_agrees_once_flagged():
    rows = {(table.name, row.name): row for table in sweep_tables.load_tables().tables() for row in table.rows}
    for key, (fill_id, replaced_by) in _registry().items():
        assert key in rows, f"{fill_id}: {key} is not a registered row"
        assert re.fullmatch(r"(GEN|MAP|TEST|PERF)\.\d+", replaced_by), f"{fill_id}: replaced_by {replaced_by!r} is not one item"
        if rows[key].state == "stand-in":
            shown = rows[key].replacement
            if re.fullmatch(r"(GEN|MAP|TEST|PERF)\.\d+", shown):      # TEST.487: a row shown on a page says what replaces it in words, not by item code
                assert shown == replaced_by, f"{fill_id}: the table row names {shown!r}, the registry {replaced_by!r}"
            else:
                assert shown and not re.search(r"\b(GEN|MAP|TEST|PERF)\.\d+", shown), f"{fill_id}: the shown replacement {shown!r} is not plain words"
