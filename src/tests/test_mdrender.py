"""
planetgen/web/lib/mdrender.py regression tests.

Covers exactly the narrow Markdown subset `StarSystem.__str__` actually
generates: ATX headers, GFM-style pipe tables, plain paragraphs, and the
one legitimate raw-HTML pattern (`<sup>...</sup>`) -- plus the escaping
behavior that keeps anything else (e.g. a `--name`-injected `<script>`)
from being rendered as live HTML.

Run with: pytest src/tests/test_mdrender.py
"""

from planetgen.web.lib.mdrender import markdown_to_html, markdown_to_html_with_headings  # noqa: E402


def test_empty_input_returns_empty_string():
    assert markdown_to_html("") == ""
    assert markdown_to_html(None) == ""


def test_headers_render_at_correct_level():
    html = markdown_to_html("# Sol\n\n## Data\n\n### Detail")
    assert '<h1 id="sol">Sol</h1>' in html
    assert '<h2 id="data">Data</h2>' in html
    assert '<h3 id="detail">Detail</h3>' in html


def test_headers_get_unique_ids_for_table_of_contents():
    html, headings = markdown_to_html_with_headings("# Kepler-42\n\n## Kepler-42 b\n\n### I\n\n### I")
    assert headings == [
        {"level": 1, "text": "Kepler-42", "id": "kepler-42"},
        {"level": 2, "text": "Kepler-42 b", "id": "kepler-42-b"},
        {"level": 3, "text": "I", "id": "i"},
        {"level": 3, "text": "I", "id": "i-2"},
    ]
    assert 'id="i"' in html
    assert 'id="i-2"' in html


def test_with_headings_returns_empty_list_for_empty_input():
    assert markdown_to_html_with_headings("") == ("", [])
    assert markdown_to_html_with_headings(None) == ("", [])


def test_pipe_table_renders_as_html_table():
    md = "| Property | Value |\n|---|---|\n| Type | G2V |\n| Mass | 1.0 |"
    html = markdown_to_html(md)
    assert "<table>" in html
    assert "<th>Property</th>" in html
    assert "<th>Value</th>" in html
    assert "<td>Type</td>" in html
    assert "<td>G2V</td>" in html


def test_paragraph_is_wrapped_and_escaped():
    html = markdown_to_html("A plain sentence about the system.")
    assert html == "<p>A plain sentence about the system.</p>"


def test_sup_tag_survives_while_everything_else_is_escaped():
    md = "Radius 4.38 x 10<sup>5</sup> km, injected <script>alert(1)</script>."
    html = markdown_to_html(md)
    assert "<sup>5</sup>" in html
    assert "<script>" not in html
    assert "&lt;script&gt;" in html


def test_ampersand_is_escaped_in_paragraphs_and_tables():
    html = markdown_to_html("The star's wind & heliosphere.")
    assert "star&#x27;s wind &amp; heliosphere" in html

    md = "| Name | Note |\n|---|---|\n| Ka'Iara | AT&T style name |"
    table_html = markdown_to_html(md)
    assert "<td>Ka&#x27;Iara</td>" in table_html
    assert "AT&amp;T style name" in table_html


def test_entities_and_markdown_syntax_in_names_stay_literal():
    html = markdown_to_html("&lt;b&gt; *star* [link](http://x) `code`\n\n- item\n\n> quote\n\n---")
    assert "&amp;lt;b&amp;gt; *star* [link](http://x) `code`" in html
    assert "<ul>" not in html and "<blockquote>" not in html and "<hr" not in html and "<code>" not in html
    assert "<p>- item</p>" in html and "<p>&gt; quote</p>" in html


def test_heading_list_has_the_text_as_written():
    html, headings = markdown_to_html_with_headings("# A & B")
    assert headings == [{"level": 1, "text": "A & B", "id": "a-b"}]
    assert '<h1 id="a-b">A &amp; B</h1>' in html


def test_table_alignment_uses_an_attribute_not_an_inline_style():
    html = markdown_to_html("| a | b |\n|:--|--:|\n| 1 | 2 |")
    assert 'align="left"' in html and 'align="right"' in html
    assert "style=" not in html


def test_table_is_wrapped_in_a_scrollable_div():
    html = markdown_to_html("| a |\n|---|\n| 1 |")
    assert html.startswith('<div class="table-scroll" tabindex="0"><table>')


def test_a_generated_system_page_renders_headings_tables_and_exponents():
    from planetgen.db import render
    from planetgen.generation.system import StarSystem
    from tests.test_systems import make_config

    text = render.render_star_system(StarSystem(system_config=make_config("G2V")), "markdown")
    html, headings = markdown_to_html_with_headings(text)
    assert headings and headings[0]["level"] == 1
    assert html.count("<table>") == html.count("</table>") >= 1
    assert "<sup>" in html
    assert "<script" not in html


def test_multiple_blocks_are_joined_in_order():
    md = "# Title\n\nFirst paragraph.\n\nSecond paragraph."
    html = markdown_to_html(md)
    assert html.index('<h1 id="title">Title</h1>') < html.index("First paragraph")
    assert html.index("First paragraph") < html.index("Second paragraph")
