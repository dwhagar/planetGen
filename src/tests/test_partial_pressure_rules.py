"""
Oxygen and carbon dioxide are judged by partial pressure wherever they gate class or life (GEN.1320, GEN.1321; the sweep GEN.1331).

`physics/data/partial_pressure_rules.json` lists every rule that was converted, every rule that deliberately stays a share (and why) and the
constants that must stay pressures. The tests read every entry of it and every rule of the three rule files, so a share-based oxygen or carbon
dioxide threshold in a converted rule, or in a rule the list does not name, fails here.
"""

import importlib
import json
import os
import pathlib
import re

import pytest

from planetgen.physics import elements, planet_classes
from planetgen.web.lib import life_help

ROOT = pathlib.Path(__file__).resolve().parents[2]
REGISTRY = json.load(open(os.path.join(elements.DATA_DIRECTORY, "partial_pressure_rules.json"), encoding="utf-8"))
PRESSURE = tuple(REGISTRY["pressure_quantities"])
SHARE = tuple(REGISTRY["share_quantities"])
GAS_OF = {"o2_partial_pa": "o2", "co2_partial_pa": "co2", "share_o2": "o2", "share_co2": "co2"}
TABLE_KEY = {"planet_class_fit.json": "classes", "life_classes.json": "classes", "age_needs.json": "rows"}


def _document(filename):
    return json.load(open(os.path.join(elements.DATA_DIRECTORY, filename), encoding="utf-8"))


def _rules(filename):
    """`(row key, kind, quantity, condition)` for every condition of a rule file that reads a gas quantity of the registry."""
    table = _document(filename)[TABLE_KEY[filename]]
    found = []
    for row in table:
        key = row.get("letter") or row.get("id") or row.get("name")
        conditions = list(row.get("when", [])) + list(row.get("gates", [])) + list(row.get("windows", []))
        for condition in conditions:
            if condition.get("quantity") in PRESSURE + SHARE:
                found.append((key, row.get("kind"), condition["quantity"], condition))
    return found


def test_the_registry_names_a_row_for_every_converted_rule_and_it_reads_a_pressure():
    for entry in REGISTRY["converted"]:
        for quantity in entry["quantities"]:
            assert quantity in PRESSURE, entry
            assert any(key == entry["row"] and quantity == q and entry.get("kind", kind) == kind for key, kind, q, _c in _rules(entry["file"])), (
                f"{entry['file']} row {entry['row']} has no {quantity} condition")


def test_no_converted_rule_reads_a_share_of_the_gas_it_judges_by_pressure():
    kept = {(e["file"], e["row"], e["quantity"]) for e in REGISTRY["kept_as_share"]}
    for entry in REGISTRY["converted"]:
        gases = {GAS_OF[q] for q in entry["quantities"]}
        for key, kind, quantity, _condition in _rules(entry["file"]):
            if key == entry["row"] and entry.get("kind", kind) == kind and quantity in SHARE and GAS_OF[quantity] in gases:
                assert (entry["file"], key, quantity) in kept, f"{entry['file']} row {key} still judges {GAS_OF[quantity]} by a share ({quantity})"


def test_every_share_of_oxygen_or_carbon_dioxide_in_a_rule_file_is_a_listed_exception_with_a_reason():
    kept = {(e["file"], e["row"], e["quantity"]): e["reason"] for e in REGISTRY["kept_as_share"]}
    seen = set()
    for filename in TABLE_KEY:
        for key, _kind, quantity, _condition in _rules(filename):
            if quantity in SHARE:
                assert (filename, key, quantity) in kept, f"{filename} row {key} reads {quantity}, which the registry does not list"
                seen.add((filename, key, quantity))
    assert seen == set(kept), f"listed but not found: {sorted(set(kept) - seen)}"
    assert all(len(reason) > 40 for reason in kept.values())


def test_every_pressure_condition_of_a_rule_file_is_a_listed_converted_rule():
    converted = {(e["file"], e["row"], q) for e in REGISTRY["converted"] for q in e["quantities"]}
    for filename in TABLE_KEY:
        for key, _kind, quantity, _condition in _rules(filename):
            if quantity in PRESSURE:
                assert (filename, key, quantity) in converted, f"{filename} row {key} reads {quantity}, which the registry does not list"


def test_a_pressure_threshold_is_in_pascals_and_not_a_share():
    """A share is at most 1; no oxygen or carbon dioxide pressure edge is under a pascal."""
    for filename in TABLE_KEY:
        for key, _kind, quantity, condition in _rules(filename):
            for bound in ("min", "lt", "low", "high"):
                if bound in condition and quantity in PRESSURE:
                    assert condition[bound] >= 1.0, f"{filename} row {key}: {quantity} {bound} {condition[bound]} looks like a share, not Pa"
                if bound in condition and quantity in SHARE:
                    assert 0.0 <= condition[bound] <= 1.0, f"{filename} row {key}: {quantity} {bound} {condition[bound]} is not a share"


def test_the_pressure_constants_are_pressures():
    for entry in REGISTRY["code_constants"]:
        value = getattr(importlib.import_module(entry["module"]), entry["name"])
        assert entry["unit"] == "Pa" and value >= entry["floor"], f"{entry['module']}.{entry['name']} = {value}"


def test_the_fit_widths_of_the_pressures_are_in_logs_of_a_ratio_and_those_of_the_shares_are_shares():
    """GEN.1314: the fudge table (GEN.1313) gives a pressure a width as the natural log of a ratio and a share a width in share points."""
    widths = _document("fudge_table.json")["quantities"]
    for quantity in PRESSURE:
        assert widths[quantity]["unit"] == "ln", quantity
    for quantity in SHARE:
        if quantity in widths:
            assert widths[quantity]["unit"] == "plain" and all(edge["width"] <= 1.0 for edge in widths[quantity]["edges"]), quantity


def test_the_class_pages_give_the_converted_gases_as_pressures_and_the_sentences_hold_no_number():
    """GEN.1432: a class sentence holds no number; the page's Conditions line (read from planet_class_fit.json) gives a converted gas as a pressure, never a share."""
    from planetgen.physics import class_fit

    sentences = {row["letter"]: row["sentence"] for row in _document("planet_classes.json")["rows"] if row.get("kind", "t") == "t"}
    for entry in REGISTRY["converted"]:
        if entry["file"] != "planet_class_fit.json":
            continue
        text = sentences[entry["row"]]
        assert not re.search(r"\d", text), f"class {entry['row']}: {text}"
        fit = next(row for row in class_fit.CLASSES if row.letter == entry["row"] and row.kind == "t")
        words = "; ".join(class_fit.rule_words(fit))
        for quantity in entry["quantities"]:
            word = {"o2": "oxygen", "co2": "carbon dioxide"}[GAS_OF[quantity]]
            assert f"{word} pressure" in words, f"class {entry['row']}: {words}"
            assert f"{word} share" not in words, f"class {entry['row']}: {words}"


def test_the_class_pages_read_the_pressures_of_the_life_gates_as_pressures():
    for gate in ({"quantity": "o2_partial_pa", "lt": 3000.0}, {"quantity": "o2_partial_pa", "min": 100.0}, {"quantity": "co2_partial_pa", "min": 5.0}):
        words = life_help._gate_words(gate)
        assert "pressure" in words and "percent" not in words and "%" not in words, words


def test_the_class_quantities_hold_the_pressure_beside_the_share():
    from tests.test_planet_classes import EARTH

    quantities = planet_classes.quantities(dict(EARTH))
    for name in PRESSURE + SHARE:
        assert name in quantities
    assert quantities["o2_partial_pa"] == pytest.approx(quantities["share_o2"] * quantities["pressure_pa"], rel=1e-6)
