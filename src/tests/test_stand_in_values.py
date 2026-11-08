"""
Every stand-in carries a source tier and a tweak level, and the list of them is the registry's (GEN.725,
docs/design/research/stand-in-values.md, scripts/build_stand_in_values.py).
"""

import pathlib
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import build_stand_in_values as builder  # noqa: E402

TEXT = (ROOT / "docs" / "design" / "research" / "gap-filling-plan.md").read_text(encoding="utf-8")
WIRING = (ROOT / "docs" / "design" / "research" / "wiring-register.md").read_text(encoding="utf-8")
BURN_DOWN = (ROOT / "docs" / "design" / "research" / "burn-down-register.md").read_text(encoding="utf-8")
EVIDENCE = (ROOT / "docs" / "design" / "research" / "evidence-register.md").read_text(encoding="utf-8")
TABLES = builder.registry_tables(TEXT)


def test_every_registry_table_is_found_and_has_rows():
    assert set(TABLES) == {key for key, *_rest in builder.sections()}
    assert all(rows for _, rows in TABLES.values())


def test_every_registry_row_has_a_source_tier_and_a_tweak_level():
    assert builder.problems(TABLES) == []


@pytest.mark.parametrize("cell,ok", [("T1", True), ("T7 (core), T2 (peak)", True), ("T9", True), ("T0", False), ("T10", False),
                                     ("", False), ("verified", False)])
def test_a_tier_cell_must_open_with_a_tier_from_one_to_nine(cell, ok):
    tables = {key: (["Id", "Tier", "Tweak"], [["X", cell, "L0"]]) for key, *_rest in builder.sections()}
    assert (builder.problems(tables) == []) is ok


@pytest.mark.parametrize("cell,ok", [("L0", True), ("L4", True), ("L5", False), ("", False), ("L1, L2", False)])
def test_a_tweak_cell_is_one_level_from_zero_to_four(cell, ok):
    tables = {key: (["Id", "Tier", "Tweak"], [["X", "T1", cell]]) for key, *_rest in builder.sections()}
    assert (builder.problems(tables) == []) is ok


def test_a_row_without_the_columns_is_named():
    tables = {key: (["Id", "Kind"], [["X", "D"]]) for key, *_rest in builder.sections()}
    assert all("no Tier and Tweak columns" in message for message in builder.problems(tables))


def test_the_list_is_the_registry_as_built():
    path = ROOT / "docs" / "design" / "research" / "stand-in-values.md"
    assert path.read_text(encoding="utf-8") == builder.render(TEXT, WIRING, BURN_DOWN, EVIDENCE), \
        "run python scripts/build_stand_in_values.py to rebuild docs/design/research/stand-in-values.md"


def test_every_row_of_the_registry_is_in_the_list():
    listed = builder.render(TEXT, WIRING, BURN_DOWN, EVIDENCE)
    for number, (header, rows) in TABLES.items():
        for row in rows:
            assert f"| {row[0]} |" in listed, (number, row[0])


def test_the_rows_the_wiring_register_names_show_their_item_in_the_list():
    listed = builder.render(TEXT, WIRING, BURN_DOWN, EVIDENCE)
    unread = builder.unread_rows(WIRING)
    assert unread, "the wiring register names no registry row"
    for key, item in unread.items():
        section, fill_id = key.split(":")
        line = next(line for line in listed.splitlines() if line.startswith(f"| {fill_id} |"))
        assert [cell.strip() for cell in line.strip("|").split("|")][-2] == item, key
    assert builder.UNREAD_COLUMN in listed


def test_every_row_has_an_evidence_status_and_none_is_blank():
    """GEN.1428: the list names one of the four statuses for every row, the default being predicted only."""
    names = {name for name, _ in builder.STATUSES}
    listed = builder.render(TEXT, WIRING, BURN_DOWN, EVIDENCE)
    count = 0
    for key, number, title, columns in builder.sections():
        block = listed.split(f"## {title}\n", 1)[1].split("\n## ", 1)[0]
        lines = [line for line in block.splitlines() if line.startswith("| ") and not line.startswith("| ---") and not line.startswith("|---")]
        header, body = lines[0], lines[1:]
        assert builder.EVIDENCE_COLUMN in header, number
        for line in body:
            status = [cell.strip() for cell in line.strip().strip("|").split("|")][-3]
            assert status in names, (number, line[:60])
            count += 1
    assert count == sum(len(rows) for _, rows in TABLES.values())


def test_the_evidence_register_names_only_real_rows_and_sections():
    assert builder.evidence_problems(TABLES, builder.evidence_rows(EVIDENCE)) == []


def test_the_materials_the_report_names_are_ordinary_observed_rows():
    """GEN.1428: BC8 carbon, superionic ice, MgO2, Na2He, FeHe, Xe-Fe and Xe-Ni, ambient cg-N are observed, not predicted."""
    register = builder.evidence_rows(EVIDENCE)
    for fill_id in ("NM1", "NM2", "NM3", "NM4", "NM5", "NM6", "NM7"):
        assert register[f"6.237:{fill_id}"][0] == "observed in a lab", fill_id
    assert register["6.229:CV1"][0] == "observed in a planet"
    assert register["6.229:CV4"][0] == "predicted only"
