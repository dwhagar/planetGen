# tests/test_chem_formulas.py

"""
UX.149: the one chemical formula writer (`planetgen.util.chem`): counts as subscripts, charges and isotope mass numbers
as superscripts, hydrates and brackets, the plain-text twin, the registry-only marking of a text that names chemicals,
the Document's `chem` item in all four renderers, and no template that prints a bare formula.
"""

import os
import re

import pytest

from planetgen.describe import html as describe_html, markdown, model, text as describe_text
from planetgen.describe.quantity import inline_text
from planetgen.util import chem

TEMPLATES = os.path.join(os.path.dirname(chem.__file__), os.pardir, "web", "templates")


@pytest.mark.parametrize("formula, html, plain", [
    ("H2O", "H<sub>2</sub>O", "H₂O"),
    ("C2H6", "C<sub>2</sub>H<sub>6</sub>", "C₂H₆"),
    ("CO2", "CO<sub>2</sub>", "CO₂"),
    ("HDO", "HDO", "HDO"),
    ("CxHy", "C<sub>x</sub>H<sub>y</sub>", "CₓHy"),
    ("(Mg0.88,Fe0.12)SiO3", "(Mg<sub>0.88</sub>,Fe<sub>0.12</sub>)SiO<sub>3</sub>", "(Mg₀.₈₈,Fe₀.₁₂)SiO₃"),
    ("Fe(1-x)O", "Fe<sub>(1-x)</sub>O", "Fe₍₁₋ₓ₎O"),
])
def test_counts_are_subscripts(formula, html, plain):
    assert chem.html(formula) == html and chem.text(formula) == plain


@pytest.mark.parametrize("formula, html, plain", [
    ("(CH3)2S", "(CH<sub>3</sub>)<sub>2</sub>S", "(CH₃)₂S"),
    ("Ca5(PO4)3F", "Ca<sub>5</sub>(PO<sub>4</sub>)<sub>3</sub>F", "Ca₅(PO₄)₃F"),
    ("Al2Si2O5(OH)4", "Al<sub>2</sub>Si<sub>2</sub>O<sub>5</sub>(OH)<sub>4</sub>", "Al₂Si₂O₅(OH)₄"),
    ("[Fe(CN)6]^4-", "[Fe(CN)<sub>6</sub>]<sup>4−</sup>", "[Fe(CN)₆]⁴⁻"),
])
def test_brackets_take_their_own_count(formula, html, plain):
    assert chem.html(formula) == html and chem.text(formula) == plain


@pytest.mark.parametrize("formula, html, plain", [
    ("CaSO4.2H2O", "CaSO<sub>4</sub>·2H<sub>2</sub>O", "CaSO₄·2H₂O"),
    ("CaSO4·2H2O", "CaSO<sub>4</sub>·2H<sub>2</sub>O", "CaSO₄·2H₂O"),
    ("Ca(UO2)2(PO4)2.11H2O", "Ca(UO<sub>2</sub>)<sub>2</sub>(PO<sub>4</sub>)<sub>2</sub>·11H<sub>2</sub>O",
     "Ca(UO₂)₂(PO₄)₂·11H₂O"),
])
def test_a_hydrate_keeps_its_coefficient_on_the_line(formula, html, plain):
    assert chem.html(formula) == html and chem.text(formula) == plain


@pytest.mark.parametrize("formula, html, plain", [
    ("3He", "<sup>3</sup>He", "³He"),
    ("14C", "<sup>14</sup>C", "¹⁴C"),
    ("13CO", "<sup>13</sup>CO", "¹³CO"),
    ("14CO2", "<sup>14</sup>CO<sub>2</sub>", "¹⁴CO₂"),
])
def test_an_isotope_s_mass_number_is_raised_before_its_symbol(formula, html, plain):
    assert chem.html(formula) == html and chem.text(formula) == plain


@pytest.mark.parametrize("formula, html, plain", [
    ("Na+", "Na<sup>+</sup>", "Na⁺"),
    ("Cl-", "Cl<sup>−</sup>", "Cl⁻"),
    ("OH-", "OH<sup>−</sup>", "OH⁻"),
    ("NH4+", "NH<sub>4</sub><sup>+</sup>", "NH₄⁺"),
    ("e-", "e<sup>−</sup>", "e⁻"),
    ("Fe2+", "Fe<sup>2+</sup>", "Fe²⁺"),
    ("Fe3+", "Fe<sup>3+</sup>", "Fe³⁺"),
    ("O2+", "O<sup>2+</sup>", "O²⁺"),
    ("O2^+", "O<sub>2</sub><sup>+</sup>", "O₂⁺"),
    ("SO4^2-", "SO<sub>4</sub><sup>2−</sup>", "SO₄²⁻"),
])
def test_a_charge_is_raised_after_the_formula(formula, html, plain):
    assert chem.html(formula) == html and chem.text(formula) == plain


@pytest.mark.parametrize("bad", ["", " H2O", "H2 O", "h2o", "Xx2", "(H2", "H2)", "H2O!", None])
def test_what_is_no_formula_is_refused(bad):
    with pytest.raises(chem.FormulaError):
        chem.parse(bad)


def test_html_is_escaped_and_every_registry_formula_reads():
    assert "<" not in chem.html("H2O").replace("<sub>", "").replace("</sub>", "")
    known = chem.known()
    assert {"H2O", "CO2", "CH4", "H2SO4", "3He", "CaSO4.2H2O", "(CH3)2S"} <= known
    assert not {"CO", "He", "Fe", "N"} & known, "a formula with nothing to write out is not marked"
    for formula in known:
        assert chem.text(formula) and chem.html(formula)


def test_a_text_marks_only_registry_formulas_and_charged_elements():
    assert chem.split("H+, e-, O2+, N+ and S+") == [
        chem.Formula("H+"), ", ", chem.Formula("e-"), ", ", chem.Formula("O2+"), ", ", chem.Formula("N+"), " and ",
        chem.Formula("S+")]
    assert chem.mark_html("H2, He, CO, PAHs, silicates") == "H<sub>2</sub>, He, CO, PAHs, silicates"
    assert chem.mark_text("electron donors like Hydrogen Sulfide (H2S), Ferrous Iron (Fe2+).") == (
        "electron donors like Hydrogen Sulfide (H₂S), Ferrous Iron (Fe²⁺).")
    for words in ("a G2V star", "2MASS", "the H II region", "Herbig-Haro jets", "3D view", "12L"):
        assert chem.split(words) == [words], words
    assert chem.mark_html("<b> H2O") == "&lt;b&gt; H<sub>2</sub>O"


def test_the_document_s_chem_item_is_written_by_every_renderer():
    value = model.chem_text("H2, CO and Fe2+")
    assert value == [{"chem": "H2"}, ", CO and ", {"chem": "Fe2+"}]
    assert model.chem_text(None) == [] and model.chem_text("") == []
    with pytest.raises(chem.FormulaError):
        model.chem("not a formula")
    assert inline_text(value, "text") == "H₂, CO and Fe²⁺"
    assert inline_text(value, "markdown") == "H<sub>2</sub>, CO and Fe<sup>2+</sup>"
    assert describe_html.inline_html(value) == "H<sub>2</sub>, CO and Fe<sup>2+</sup>"
    document = model.document("phenomenon", "nebula:1", "A Cloud", [
        model.section("nebula", "nebula:1", "A Cloud", [model.facts("Facts", [model.fact("species", "Dominant Gas", value)])])])
    model.check(document)
    assert {"chem": "H2"} in list(model.walk_inline(document))
    assert "H₂, CO and Fe²⁺" in describe_text.render(document)
    assert "H<sub>2</sub>, CO and Fe<sup>2+</sup>" in markdown.render(document)
    assert "H<sub>2</sub>, CO and Fe<sup>2+</sup>" in describe_html.render(document)


def _template_text(source):
    """What a template writes itself: its text with every Jinja expression, statement, comment and tag taken out."""
    for pattern in (r"\{#.*?#\}", r"\{%.*?%\}", r"\{\{.*?\}\}", r"<script\b.*?</script>", r"<style\b.*?</style>",
                    r"<[^>]*>"):
        source = re.sub(pattern, " ", source, flags=re.S)
    return source


def test_no_template_prints_a_bare_formula():
    """A formula a template writes itself must go through `|formula` (or a caller's `chem`): none is typed bare."""
    found = []
    for folder, _dirs, names in os.walk(TEMPLATES):
        for name in names:
            if name.endswith(".html"):
                with open(os.path.join(folder, name), encoding="utf-8") as handle:
                    words = _template_text(handle.read())
                found += [(name, part) for part in chem.split(words) if isinstance(part, chem.Formula)]
    assert not found, found


# --- The pages that print formulas today ------------------------------------------------------------------------

def test_a_cloud_class_s_contents_and_summary_write_their_formulas_out():
    from planetgen.web.class_pages import _plain
    from planetgen.web.lib.classref import class_entry
    dark = class_entry("nebula", "M")
    assert "Dark nebula: H<sub>2</sub>, He, CO, PAHs, silicates." == dark["summary"]
    assert dict(dark["facts"])["Contents"] == "H<sub>2</sub>, He, CO, PAHs, silicates"
    assert _plain(dark["summary"]) == "Dark nebula: H2, He, CO, PAHs, silicates."
    emission = class_entry("nebula", "D")
    assert "H<sup>+</sup>, e<sup>−</sup>, O<sup>2+</sup>, N<sup>+</sup> and S<sup>+</sup>" in emission["summary"]


def test_a_nebula_page_s_dominant_species_and_a_species_limit_keep_their_formulas():
    from planetgen.web.population_pages import _limit_label, _limit_rows
    from planetgen.web.system_pages import FIELD_SPECS
    species = dict((column, formatter) for column, _label, formatter in FIELD_SPECS["nebula"])["dominant_species"]
    assert species("H2, CO, cold dust") == "H<sub>2</sub>, CO, cold dust" and species("") is None
    assert _limit_label("inspired_o2_min_kpa") == "Inspired O<sub>2</sub> min"
    assert _limit_label("comfort_high_c") == "Comfort high" and _limit_label("pressure_no_gear_kpa") == "Pressure no gear"
    (row,) = _limit_rows([{"name": "inspired_o2_min_kpa", "human": 8.0, "value": 8.0, "unit": "kPa", "changed": False,
                           "meaning": "least inspired O2 for consciousness"}])
    assert row["meaning"] == "least inspired O<sub>2</sub> for consciousness"
