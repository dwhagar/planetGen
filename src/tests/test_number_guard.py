# tests/test_number_guard.py

"""
GEN.772, Boss 9:03 PM EDT 2026-10-24: "we should never have numbers without units or context". The guard has three
layers, and this file holds the ones that are not about the property sheet's own rows (`test_property_sheet.py`):

* a `ratio` quantity of a Document cannot be built without a `note` (its unit, scale or noun), and a facts row whose
  value reads as a number alone fails `describe.model.check`;
* every Document of a generated system of each variant and of each phenomenon passes that check;
* a body's page, as the site writes it, holds no table cell that is a number alone (`bare_cells`).
"""

import re

import pytest

from planetgen import describe
from planetgen.db import column_units
from planetgen.describe import model, sheet
from planetgen.describe.quantity import inline_text
from planetgen.generation.phenomena.asteroid_field import AsteroidField
from planetgen.generation.phenomena.compact_remnant import BlackHole, NeutronStar
from planetgen.generation.phenomena.nebula import Nebula
from planetgen.generation.phenomena.quasar import Quasar
from planetgen.generation.phenomena.rogue import InterstellarComet, RoguePlanet
from planetgen.generation.phenomena.supernova_remnant import SupernovaRemnant
from planetgen.generation.config import SystemConfig
from planetgen.util import draw
from tests.test_describe import _systems

CELL = re.compile(r"<t[dh][^>]*>(.*?)</t[dh]>", re.S)
TAG = re.compile(r"<[^>]+>")


def bare_cells(page):
    """The table cells of an HTML page whose whole text is a measured number (it has a decimal point or a power of ten):
    a figure with no unit, scale or noun beside it. A whole number is a count, and its row's label is its noun."""
    found = []
    for cell in CELL.findall(page):
        text = TAG.sub("", cell).replace("&nbsp;", " ").strip()
        if model.is_bare_number(text) and ("." in text or "×" in text):
            found.append(text)
    return found


def test_a_ratio_cannot_be_built_without_saying_what_it_is():
    with pytest.raises(ValueError, match="bare number"):
        model.qty("ratio", 0.5)
    assert inline_text([model.qty("ratio", 0.0167, note="(0 is a circle, 1 is the escape orbit)")]) \
        == "0.017 (0 is a circle, 1 is the escape orbit)"


def test_a_facts_row_that_is_a_number_alone_fails_the_document_check():
    """The guard fails on a deliberately bare row: this is what a new fact with a bare number would hit."""
    bare = model.document("system", "", "A", [model.section("system", "", "A", [model.facts("x", [model.fact("v", "Value", "0.71")])])])
    with pytest.raises(ValueError, match="no unit, scale or noun"):
        model.check(bare)
    fine = model.document("system", "", "A", [model.section("system", "", "A", [
        model.facts("x", [model.fact("v", "Value", "0.71 of 1"), model.fact("c", "Planets", model.qty("count", 4))])])])
    model.check(fine)


def test_the_property_sheet_refuses_a_column_declared_without_a_unit(monkeypatch):
    spec = next(spec for spec in sheet.SPECS["planet"] if spec.key == "porosity")
    monkeypatch.setitem(column_units.COLUMNS, ("planets", "porosity"), "ratio")
    with pytest.raises(ValueError, match="bare number"):
        sheet._value(spec, "planets", 0.3)


def test_every_generated_system_has_no_bare_number_in_its_facts():
    for _overrides, system in _systems(count=2):
        model.check(describe.describe(system))


@pytest.mark.parametrize("cls", [Nebula, Quasar, RoguePlanet, InterstellarComet, AsteroidField, SupernovaRemnant, BlackHole,
                                 NeutronStar])
def test_every_phenomenon_has_no_bare_number_in_its_facts(cls):
    for seed in range(8):
        draw.set_run_seed(seed)
        model.check(describe.describe_phenomenon(cls(SystemConfig())))


def test_a_bare_cell_is_found_in_a_page():
    assert bare_cells("<table><tr><th>Water activity</th><td>0.71</td></tr><tr><td>400 ppm</td><td><b>1,200</b></td></tr></table>") \
        == ["0.71"]
