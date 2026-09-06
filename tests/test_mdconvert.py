"""
html/lib/mdconvert.py regression tests.

Covers exactly the narrow Markdown subset `StarSystem.__str__` actually
generates: ATX headers, GFM-style pipe tables, plain paragraphs, and the
one legitimate raw-HTML pattern (`<sup>...</sup>`) -- plus the escaping
behavior that keeps anything else (e.g. a `--name`-injected `<script>`)
from being rendered as live HTML.

`html/lib` isn't part of the installed `stellarObjects` package (CGI-only
plumbing), so it's added to `sys.path` here the same way the CGI scripts
themselves do.

Run with: pytest tests/test_mdconvert.py
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "html", "lib"))

from mdconvert import markdown_to_html  # noqa: E402


def test_empty_input_returns_empty_string():
    assert markdown_to_html("") == ""
    assert markdown_to_html(None) == ""


def test_headers_render_at_correct_level():
    html = markdown_to_html("# Sol\n\n## Data\n\n### Detail")
    assert "<h1>Sol</h1>" in html
    assert "<h2>Data</h2>" in html
    assert "<h3>Detail</h3>" in html


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


def test_apostrophe_and_ampersand_are_escaped_in_paragraphs_and_tables():
    html = markdown_to_html("The star's wind & heliosphere.")
    assert "&#x27;" in html
    assert "&amp;" in html

    md = "| Name | Note |\n|---|---|\n| Ka'Iara | AT&T style name |"
    table_html = markdown_to_html(md)
    assert "&#x27;" in table_html
    assert "&amp;" in table_html


def test_multiple_blocks_are_joined_in_order():
    md = "# Title\n\nFirst paragraph.\n\nSecond paragraph."
    html = markdown_to_html(md)
    assert html.index("<h1>Title</h1>") < html.index("First paragraph")
    assert html.index("First paragraph") < html.index("Second paragraph")
