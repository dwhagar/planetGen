"""
The stand-in registry rows of the life class rules (GEN.893, docs/design/research/gap-filling-plan.md section 6.19): every class of the
rule file has a row, the row names every gate, window, solvent and supply of its class, and its tier and tweak are those of the
least certain number in the class.
"""

import json
import pathlib

import pytest

from planetgen.physics import life_classes

PLAN = pathlib.Path(__file__).resolve().parents[2] / "docs" / "design" / "research" / "gap-filling-plan.md"


def _rows():
    block = PLAN.read_text().split("<!-- life-class-standins:begin -->")[1].split("<!-- life-class-standins:end -->")[0]
    rows = [[cell.strip() for cell in line.strip().strip("|").split("|")] for line in block.strip().splitlines()]
    return [row for row in rows[2:] if len(row) == 10]


def _classes():
    with open(life_classes.DATA_FILE, encoding="utf-8") as handle:
        return json.load(handle)["classes"]


def _numbers(row):
    return [row["solvent"], row["energy"], *row["gates"], *row["windows"], *([row["supply"]] if row.get("supply") else [])]


def _by_class():
    return {row[4].split(":")[0].removeprefix("life_classes.json Class "): row for row in _rows() if row[4].startswith("life_classes.json Class ")}


def test_every_class_of_the_rule_file_has_one_row_and_no_row_names_a_class_that_is_gone():
    assert sorted(_by_class()) == sorted(c["id"] for c in _classes())


@pytest.mark.parametrize("rule", _classes(), ids=lambda c: c["id"])
def test_the_row_of_a_class_has_the_tier_and_tweak_of_its_least_certain_number(rule):
    row = _by_class()[rule["id"]]
    assert row[2] == "T" + str(max(int(n["tier"][1:]) for n in _numbers(rule)))
    assert row[3] == "L" + str(max(int(n["tweak"][1:]) for n in _numbers(rule)))
    assert row[1] == ("D" if row[3] == "L0" else "P")


@pytest.mark.parametrize("rule", _classes(), ids=lambda c: c["id"])
def test_the_row_of_a_class_names_every_gate_window_solvent_and_supply_of_it(rule):
    text = _by_class()[rule["id"]][4]
    assert rule["solvent"]["medium"] in text
    for number in [*rule["gates"], *rule["windows"], *([rule["supply"]] if rule.get("supply") else [])]:
        assert number["quantity"] in text, (rule["id"], number["quantity"])


def test_the_shared_thresholds_and_references_have_rows():
    text = " ".join(row[4] for row in _rows())
    for word in ("complexity thresholds", "hydrogen_inventory_floor", "reference_star_k"):
        assert word in text
