# tests/test_design_flags.py

"""
The design-flag registry (GEN.525): the life-linked gases m3, m4 and m11 are labelled design brackets that are off unless a named
flag allows them, and never tuned (decision D4 of docs/design/composition/tunable-gaps.md).
"""

import copy
import json

import pytest

from planetgen.generation import tuning_harness as th
from planetgen.generation.possibility_table import Subject
from planetgen.physics import design_flags as df
from planetgen.physics.atmosphere_recipes import RECIPES, Recipes
from planetgen.physics.elements import DataError


def test_the_recipes_honour_the_registry():
    assert df.check() == ()


def test_the_registry_lists_the_three_life_linked_gases():
    assert set(df.FLAGS) == {"m3", "m4", "m11"}
    assert df.flag("m3").enabled_by == ("life_chemical", "sep_reducing")
    assert df.flag("m4").enabled_by == ("life_chemical",)
    assert df.flag("m11").enabled_by == ("life_chemical", "lightning")


@pytest.mark.parametrize("row", sorted(df.FLAGS))
def test_a_flagged_row_is_off_without_a_flag_and_on_within_its_bracket_with_each_flag(row):
    item = df.flag(row)
    low, high, _ = RECIPES.entry(row).constants[item.constant]
    off = Subject("A", {q: False for q in item.enabled_by})
    assert RECIPES.raw_weight(RECIPES.entry(row), off, gated=False) == 0.0
    for q in item.enabled_by:
        on = Subject("A", {k: k == q for k in item.enabled_by})
        weight = RECIPES.raw_weight(RECIPES.entry(row), on, gated=False)
        assert low <= weight <= high, "the bracket midpoint lies inside the labelled bracket"


def test_a_flagged_row_is_neither_tunable_nor_a_stand_in_nor_in_the_target_table():
    for row in df.FLAGS:
        entry = RECIPES.entry(row)
        assert not entry.tunable and entry.state == "partly verified" and "design flag" in entry.limit
        assert all(t.row != f"atmosphere-m:{row}" for t in th.TARGETS)


def test_the_check_names_a_row_that_loses_its_gate_or_gains_a_target():
    ungated = copy.deepcopy(RECIPES.document)
    next(e for e in ungated["entries"] if e["id"] == "m4")["weight"] = "c_m4"
    assert any("m4" in p and "not zero" in p for p in df.check(recipes=Recipes(ungated, register=False)))
    tunable = copy.deepcopy(RECIPES.document)
    next(e for e in tunable["entries"] if e["id"] == "m3")["tunable"] = True
    assert any("m3" in p and "never tunable" in p for p in df.check(recipes=Recipes(tunable, register=False)))


def test_the_loader_refuses_a_flag_that_is_not_labelled_design(tmp_path):
    document = json.load(open(df.DATA_FILE, encoding="utf-8"))
    document["flags"][0]["label"] = "READ"
    path = tmp_path / "flags.json"
    path.write_text(json.dumps(document), encoding="utf-8")
    with pytest.raises(DataError, match="DESIGN"):
        df.load(str(path))
    document["flags"][0]["label"] = "DESIGN"
    document["flags"][0]["enabled_by"] = []
    path.write_text(json.dumps(document), encoding="utf-8")
    with pytest.raises(DataError, match="at least one flag"):
        df.load(str(path))
