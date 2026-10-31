# tests/test_life_help.py

"""
DOC.84: the help pages about life (`/help/life`): the kinds of life, the four scales (abundance, trend, security, suitability), the
stages and animal grades, and a page for each of the nine life classes. Every number on them is read from the code that uses it,
and every sentence keeps the plain-words style guide.
"""

import ast
import os
import re

import pytest

from planetgen import tuning
from planetgen.physics import life_air, life_classes, life_readings
from planetgen.util import plainwords
from planetgen.web.lib import life_help
from tests.test_class_pages import _html, app, client  # noqa: F401

LEVEL_CHECK = {key: name.lower() for key, name in life_help.LEVEL_NAMES.items()}

ALL_PAGES = [("kinds", None)] + [(slug, None) for slug, _t, _s in life_help.TOPICS if slug != "kinds"] + \
            [("classes", number) for number in range(1, 10)]


def _texts(page):
    """Every sentence-like text of a page: its intro, its paragraphs, the cells of its tables and the values of its facts."""
    found = list(page["intro"])
    for section in page["sections"]:
        found += section["paragraphs"]
        for row in (section["table"] or {"rows": []})["rows"]:
            found += [cell[0] if isinstance(cell, tuple) else cell for cell in row]
        for _label, value in section["facts"] or []:
            found += [value] if isinstance(value, str) else list(value)
    return found


@pytest.mark.parametrize("topic, number", ALL_PAGES)
def test_every_life_help_page_renders_with_the_other_pages_listed(client, topic, number):
    path = "/help/life" if topic == "kinds" else f"/help/life/classes/{number}" if number else f"/help/life/{topic}"
    html = _html(client, path)
    assert 'id="life-help"' in html
    for slug, title, _summary in life_help.TOPICS:
        assert f">{title}</a>" in html
    assert 'aria-current="page"' in html


def test_a_page_that_does_not_exist_is_a_404(client):
    for path in ("/help/life/nope", "/help/life/kinds", "/help/life/classes/10", "/help/life/classes/0"):
        assert client.get(path).status_code == 404, path
    assert life_help.page("nope") is None and life_help.page("classes", 12) is None


def test_the_help_menu_has_a_life_entry(client):
    html = _html(client, "/help/life")
    assert re.search(r'<sl-menu-item data-href="/help/life"[^>]*>Life</sl-menu-item>', html)
    assert re.search(r'<sl-menu-item data-href="/help/life"[^>]*aria-current="page"', _html(client, "/help/life/trend"))


@pytest.mark.parametrize("topic, number", ALL_PAGES)
def test_every_sentence_keeps_the_style_guide(topic, number):
    page = life_help.page(topic, number)
    broken = [(text, rule) for text in _texts(page) for sentence in plainwords.sentences(text) for rule, _detail in plainwords.check_sentence(sentence)]
    assert not broken, broken


def test_the_abundance_page_has_the_bands_of_the_code():
    table = life_help.page("abundance")["sections"][0]["table"]
    assert [row[0] for row in table["rows"]] == ["Absent", "Rare", "Occasional", "Frequent", "Abundant", "Dominant"]
    assert [row[0] for row in table["rows"]] == [word for _low, word in life_readings.ABUNDANCE_BANDS]
    assert [row[1] for row in table["rows"]] == ["0", "1 to 2", "3 to 4", "5 to 6", "7 to 8", "9"]


def test_the_trend_and_security_pages_state_the_numbers_the_code_uses():
    trend = " ".join(_texts(life_help.page("trend")))
    assert f"{life_readings.TREND_SPAN_GYR:g} billion years" in trend and f"{life_readings.TREND_BAND:.0%}" in trend
    assert f"{life_readings.SPREAD_GYR:g} billion years" in trend
    security = " ".join(_texts(life_help.page("security")))
    for value in (life_readings.DYING_GYR, life_readings.PRECARIOUS_GYR, life_readings.MARGINAL_GYR, life_readings.STEP_GYR):
        assert f"{value:g} billion years" in security
    assert f"{life_readings.MARGINAL_SCORE} or less" in security


def test_the_kinds_and_suitability_pages_state_the_tuning_they_use():
    kinds = " ".join(_texts(life_help.page("kinds")))
    assert f"{tuning.LIFE_SCORE_MAX} (abundant)" in kinds and f"never score above {tuning.LIFE_SCORE_SIMPLE_FLORA_MAX}" in kinds
    assert f"{tuning.HABITABILITY_LIFE_SCORE_MIN:.1%}" in kinds
    assert f"{int(tuning.LIFE_SCORE_SENTIENT_FULL_GYR * 1000):,} million years" in kinds
    suitability = " ".join(_texts(life_help.page("suitability")))
    assert f"{tuning.LIFE_SCORE_MAX} times that number" in suitability and f"{tuning.HABITABILITY_LIFE_SCORE_MIN:.1%}" in suitability


def test_the_stages_page_has_every_rung_grade_and_pace():
    sections = {section["id"]: section for section in life_help.page("stages")["sections"]}
    assert len(sections["stages"]["table"]["rows"]) == 7 and [row[0] for row in sections["grades"]["table"]["rows"]] == ["A1", "A2", "A3"]
    assert {row[0] for row in sections["pace"]["table"]["rows"]} == set(tuning.EVOLUTIONARY_PACE_LABELS.values())


def test_every_class_has_a_page_with_all_its_variants_and_the_rule_files_words():
    assert {row.number for row in life_classes.CLASSES} == set(range(1, 10))
    for number in range(1, 10):
        page = life_help.page("classes", number)
        headings = " ".join(section["heading"] for section in page["sections"])
        for row in life_help.family(number):
            assert f"Class {row.id}" in headings
            facts = dict(next(s for s in page["sections"] if s["id"] == f"class-{row.id}")["facts"])
            assert facts["Drinks"] == row.drinks and facts["Eats"] == row.eats and facts["Breathes"] == row.breathes
            assert facts["Gives off"] == row.gives_off and facts["How far it can go"] == row.reaches
            assert facts["Confidence it can exist"].startswith(row.tier)
            assert isinstance(facts["Needs"], list) and facts["Needs"]


def test_a_class_is_marked_speculative_when_its_air_chemistry_is_at_level_l3_or_l4():
    for number in range(1, 10):
        entry = life_air.CLASSES[str(number)]
        text = " ".join(next(s for s in life_help.page("classes", number)["sections"] if s["id"] == "sources")["paragraphs"])
        assert text.startswith("Speculative chemistry") == (entry["tweak"] in ("L3", "L4"))
        assert entry["tier"] in text and entry["tweak"] in text and LEVEL_CHECK[entry["tweak"]] in text.lower()
    rows = next(s for s in life_help.page("classes")["sections"] if s["id"] == "sources")["table"]["rows"]
    assert [row[3] for row in rows] == ["Yes" if life_air.CLASSES[str(n)]["tweak"] in ("L3", "L4") else "No" for n in range(1, 10)]


def test_the_needs_of_a_class_come_from_its_rule_row():
    needs = " ".join(life_help.needs(life_classes.BY_ID["1a"]))
    assert "255 to 395 K" in needs and "10 Pa" in needs and "5 Pa" in needs and "1% of Earth's light" in needs
    assert "518 kPa" in " ".join(life_help.needs(life_classes.BY_ID["6"]))
    assert "sulfuric acid" in " ".join(life_help.needs(life_classes.BY_ID["5b"]))
    assert any("hydrogen" in need for need in life_help.needs(life_classes.BY_ID["3a"]))


def test_the_tier_names_are_the_registrys():
    path = os.path.join(os.path.dirname(__file__), "..", "..", "scripts", "build_stand_in_values.py")
    with open(path, encoding="utf-8") as handle:
        tree = ast.parse(handle.read())
    tiers = next(ast.literal_eval(node.value) for node in tree.body if isinstance(node, ast.Assign) and node.targets[0].id == "TIERS")
    for key, name, _meaning in tiers:
        assert life_help.TIER_NAMES[key].startswith(name.split(" (")[0]), key
