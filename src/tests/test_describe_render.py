# tests/test_describe_render.py

"""
`planetgen.describe` renderers (API.25): the Document written as text, Markdown and wikitext.

All three carry the same facts and sentences as the Document, a power of ten is marked the page
format's way, and the plain text has no markup at all.
"""

import re

import pytest

from planetgen import describe
from planetgen.db import render as system_render, store
from planetgen.describe import model, quantity
from planetgen.describe.render import share_text
from planetgen.generation.config import SystemConfig
from planetgen.generation.system import StarSystem
from planetgen.util import draw

_VARIANTS = [
    dict(BINARY_SYSTEM=False, MOONS=True, COMETS=True, ASTEROID_BELT=True, HABITABLE_WORLD=True),
    dict(BINARY_SYSTEM=True, WIDE_BINARY=False, MOONS=True, COMETS=True),
    dict(BINARY_SYSTEM=True, WIDE_BINARY=True, MOONS=True, COMETS=True),
]


def _documents():
    for index, overrides in enumerate(_VARIANTS):
        draw.set_run_seed(100 + index)
        config = SystemConfig()
        for key, value in overrides.items():
            setattr(config, key, value)
        yield describe.describe(StarSystem(system_config=config))


def _facts(document):
    for section in model.sections(document):
        for block in section["blocks"]:
            if block["type"] == "facts":
                yield section, block


def test_the_formats_are_the_three_renderers():
    assert set(describe.RENDERERS) == {"text", "markdown", "wikitext"}
    assert system_render.FORMATS == tuple(describe.RENDERERS)


def test_an_unknown_format_is_refused():
    with pytest.raises(ValueError):
        describe.render({"kind": "system", "ref": None, "title": "x", "sections": []}, "html")


@pytest.mark.parametrize("fmt", ["text", "markdown", "wikitext"])
def test_every_format_carries_every_title_and_every_fact_label_or_key(fmt):
    for document in _documents():
        out = describe.render(document, fmt)
        for section in model.sections(document):
            assert section["title"] in out
        for section, block in _facts(document):
            for row in block["rows"]:
                if fmt == "wikitext" and (section["kind"] in ("star", "planet", "moon") or block["caption"].endswith("binary system")
                                          or block["caption"] == "Binary system"):
                    assert f"|{row['key']}=" in out
                else:
                    assert row["label"] in out


def test_prose_is_the_same_sentences_in_every_format():
    for document in _documents():
        outputs = {fmt: describe.render(document, fmt) for fmt in describe.RENDERERS}
        for section in model.sections(document):
            for block in section["blocks"]:
                if block["type"] == "prose":
                    for sentence in block["text"]:
                        for fmt, out in outputs.items():
                            assert quantity.inline_text(sentence, fmt) in out, (fmt, sentence)


def test_plain_text_has_no_markup():
    for document in _documents():
        out = describe.render(document, "text")
        for mark in ("**", "'''", "<sup>", "{{", "[[", "|---", "{|"):
            assert mark not in out, mark
        assert not re.search(r"^#", out, re.M)


def test_a_power_of_ten_is_marked_the_page_formats_way():
    radius = {"qty": {"kind": "radius", "value": 1.71e8, "unit": "m"}}
    assert quantity.inline_text([radius], "text") == "1.71 × 10⁵ km"
    assert quantity.inline_text([radius], "markdown") == "1.71 × 10<sup>5</sup> km"
    assert quantity.inline_text([radius], "wikitext") == "{{Exp|1.71|5}} km"
    tiny = {"qty": {"kind": "radius", "value": 2e-3, "unit": "m"}}
    assert quantity.inline_text([tiny], "wikitext") == "{{Exp|2.00|-6}} km"


def test_emphasis_is_written_per_format():
    item = [model.strong("Blue")]
    assert quantity.inline_text(item, "text") == "Blue"
    assert quantity.inline_text(item, "markdown") == "**Blue**"
    assert quantity.inline_text(item, "wikitext") == "'''Blue'''"


def test_quantities_read_like_the_unit_ladders():
    assert quantity.inline_text([model.qty("pressure", 70400.0)], "text") == "70.4 kPa (0.695 atm, 10.2 psi)"
    assert quantity.inline_text([model.qty("count", 3)], "text") == "3"
    assert quantity.inline_text([model.qty("length", None)], "text") == "–"
    assert share_text({"share": 0.5}) == "50%"
    assert share_text({"share": 0.0001}) == "< 0.1%"
    assert share_text({"level": "trace"}) == "trace"


def test_headings_follow_the_tree_in_each_format():
    document = next(iter(_documents()))
    markdown = describe.render(document, "markdown")
    wikitext = describe.render(document, "wikitext")
    assert markdown.startswith(f"# {document['title']}\n")
    assert wikitext.startswith(f"= {document['title']} =\n")
    assert wikitext.rstrip().endswith("[[Category:Star Systems]]")
    moons = [s for s in model.sections(document) if s["kind"] == "moon"]
    assert moons and f"\n### {moons[0]['title']}\n" in markdown
    assert f"\n=== {moons[0]['title']} ===\n" in wikitext


def test_the_wiki_templates_are_the_old_ones():
    document = next(iter(_documents()))
    wikitext = describe.render(document, "wikitext")
    for name in ("Star Data", "Planet Data", "Class Data"):
        assert "{{" + name + "\n" in wikitext
    binary = describe.render(list(_documents())[1], "wikitext")
    assert "{{Binary System Data\n" in binary
    assert "{{Wide Binary System Data\n" in describe.render(list(_documents())[2], "wikitext")


def test_a_stored_system_renders_like_the_generated_one(mysql_config):
    draw.set_run_seed(7)
    config = SystemConfig()
    config.MOONS = True
    system = StarSystem(system_config=config)
    system_id = store.save_system(system, config, config=mysql_config)
    conn = store.get_connection(mysql_config)
    try:
        for fmt in describe.RENDERERS:
            stored = system_render.render_system_text(conn, system_id, fmt)
            assert stored.splitlines()[:3] == system_render.render_star_system(system, fmt).splitlines()[:3]
        assert system_render.describe_system(conn, system_id)["title"] == system.name
        with pytest.raises(ValueError):
            system_render.render_system_text(conn, system_id, "html")
    finally:
        conn.close()
