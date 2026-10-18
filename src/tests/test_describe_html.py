# tests/test_describe_html.py

"""
UX.119: the HTML renderer of a Document (`planetgen.describe.html`): a snapshot for every block type, every quantity as a
`qty` span the unit preference can change, and every piece of text escaped.
"""

import re

import pytest

from planetgen import describe
from planetgen.describe import html, model
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.util import draw

HOSTILE = '<script>alert("x")</script> & "quoted"'


def _document(*blocks, children=(), title="Sol <b>"):
    return model.document("system", "system:1", title, [model.section("planet", "planet:7", title, blocks, children)])


def test_a_facts_block_is_a_table_with_a_quantity_span_in_each_measured_value():
    out = html.render(_document(model.facts("Planet Data", [
        model.fact("distance", "Distance", model.qty("length", 1.496e11)),
        model.fact("class", "Class", "Terrestrial"),
        None])))
    assert 'class="doc-facts"' in out and '<span class="qty" data-q="length" data-si="149600000000.0">' in out
    assert "<td>Terrestrial</td>" in out and "<h2>Sol &lt;b&gt;</h2>" in out


def test_a_power_of_ten_is_marked_up_and_a_quantity_without_a_ladder_is_plain_text():
    out = html.inline_html([model.qty("length", 3.0e-9), " ", model.qty("radius", 6.371e6), " ", model.qty("count", 12000)])
    assert "<sup>" in out and "×" in out
    assert out.count('class="qty"') == 1  # the length; a radius and a count are not on a ladder
    assert "12,000" in out and "km" in out
    assert html.inline_html([model.qty("length", None)]) == "–"


def test_a_ratings_block_names_each_tier_with_a_colour_class():
    block = model.ratings({"name": "PHI-6", "score": 0.82, "tier": "Green"},
                          [{"domain": "Pressure", "tier": "Blue", "score": 1.0, "text": ["Pressure is ideal."]},
                           {"domain": "Air", "tier": "Red", "score": 0.1, "text": ["Too much carbon dioxide."]}],
                          ["a mask with a scrubber"])
    out = html.render(_document(block))
    assert 'class="doc-ratings"' in out
    assert '<span class="tier tier-green">Green</span>' in out
    assert '<span class="tier tier-blue">Blue</span>' in out and '<span class="tier tier-red">Red</span>' in out
    assert "a visitor needs a mask with a scrubber." in out
    assert "<td>0.82</td>" not in out and "<td>1.00</td>" in out and "<td>0.10</td>" in out


def test_a_composition_block_lists_shares_or_a_summary():
    block = model.composition("Composition", [{"name": "iron", "share": 0.321}, {"name": "ice", "level": "trace"}],
                              summary="mostly rock")
    out = html.render(_document(block))
    assert '<strong>Composition</strong>: mostly rock.' in out
    assert '<th scope="row">Iron</th><td>32.1%</td>' in out and '<th scope="row">Ice</th><td>trace</td>' in out


def test_prose_carries_its_role_and_a_link_only_when_the_page_gives_one():
    block = model.prose("life", ["Orbits ", model.ref("planet:9", "Terra"), " and ", model.strong("glows"), "."])
    plain = html.render(_document(block))
    assert '<p class="doc-prose doc-life">Orbits Terra and <strong>glows</strong>.</p>' in plain
    linked = html.render(_document(block), link=lambda ref: f"/planet/{ref.split(':')[1]}")
    assert 'Orbits <a href="/planet/9">Terra</a> and' in linked
    assert "Terra" in html.render(_document(block), link=lambda ref: None)


def test_a_sections_blocks_alone_have_no_heading_and_no_children():
    child = model.section("moon", "moon:1", "Luna", [model.prose("summary", "A moon.")])
    top = _document(model.prose("summary", "A star.", ), children=[child])["sections"][0]
    out = html.blocks_html(top["blocks"])
    assert out == '<p class="doc-prose doc-summary">A star.</p>'
    assert html.blocks_html([]) == ""


def test_children_nest_one_heading_deeper():
    child = model.section("moon", "moon:1", "Luna", [model.prose("summary", "A moon.")])
    out = html.render(_document(children=[child]), first_heading=3)
    assert "<h3>Sol &lt;b&gt;</h3>" in out and "<h4>Luna</h4>" in out
    assert out.count("<section") == 2 == out.count("</section>")


def test_nothing_the_document_holds_comes_out_unescaped():
    hostile = _document(
        model.facts(HOSTILE, [model.fact(HOSTILE, HOSTILE, HOSTILE, model.ref(HOSTILE, HOSTILE), model.strong(HOSTILE))]),
        model.prose("summary", HOSTILE),
        model.composition(HOSTILE, [{"name": HOSTILE, "level": HOSTILE}], summary=HOSTILE),
        model.ratings({"name": HOSTILE, "score": 0.5, "tier": HOSTILE},
                      [{"domain": HOSTILE, "tier": HOSTILE, "score": 0.5, "text": [HOSTILE]}], [HOSTILE]),
        title=HOSTILE)
    out = html.render(hostile, link=lambda ref: ref)
    assert "<script" not in out and "alert(" not in out.replace("alert(&quot;x&quot;)", "")
    assert not re.search(r'<(?!/?(?:section|h\d|div|table|caption|thead|tbody|tr|th|td|p|span|strong|a|sup)\b)', out)
    assert "&lt;script&gt;" in out


@pytest.mark.parametrize("overrides", [
    dict(BINARY_SYSTEM=False, MOONS=True, COMETS=True, ASTEROID_BELT=True, HABITABLE_WORLD=True),
    dict(BINARY_SYSTEM=True, WIDE_BINARY=False, MOONS=True, COMETS=True),
])
def test_a_generated_document_renders_every_title_and_balanced_markup(overrides):
    draw.set_run_seed(119)
    config = SystemConfig()
    for key, value in overrides.items():
        setattr(config, key, value)
    document = describe.describe(StarSystem(system_config=config))
    out = html.render(document)
    for section in model.sections(document):
        assert html_escape(section["title"]) in out
    for tag in ("section", "table", "tbody", "div", "p"):
        assert len(re.findall(rf"<{tag}\b", out)) == out.count(f"</{tag}>"), tag
    assert 'class="qty"' in out


def html_escape(text):
    import html as stdlib
    return stdlib.escape(text)
