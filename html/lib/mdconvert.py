# html/lib/mdconvert.py

"""
A small, purpose-built Markdown-to-HTML converter for exactly the
Markdown `StarSystem.__str__` generates (see
`stellarObjects/utils.py:properties_to_string` and `systemData.py`) --
not a general-purpose Markdown implementation.

That output uses a narrow, fixed subset: ATX headers (`#`/`##`/`###`,
never trailing-hash-closed in Markdown mode), GFM-style pipe tables
(`| Property | Value |` + a `|---|---|` separator row), plain-text
paragraphs, and exactly one form of inline raw HTML --
`<sup>exponent</sup>`, from `to_scientific_notation`. No lists, links,
emphasis, blockquotes, or code blocks appear anywhere in generated
output (confirmed by grepping the generator for markdown syntax markers).
Reaching for a full CommonMark parser for that surface would be the
wrong tool -- and `html/` is deliberately dependency-free (standard
library only), so pulling one in isn't an option anyway.

Security note: `markdown_content` is trusted-ish generated data, but not
fully -- a `--name`/`--star-type` override on the CLI can inject
arbitrary characters into it. Every block is HTML-escaped in full before
any markup is emitted, and the one legitimate raw-HTML pattern
(`<sup>...</sup>`) is re-enabled afterward via a narrow, specific regex
over the *escaped* text -- never by leaving arbitrary input unescaped.
"""

import html
import re

_HEADER_RE = re.compile(r'^(#{1,6})\s+(.*)$')
_TABLE_ROW_RE = re.compile(r'^\|(.+)\|\s*$')
_TABLE_SEPARATOR_RE = re.compile(r'^\|?\s*:?-{2,}:?\s*(\|\s*:?-{2,}:?\s*)*\|?$')
_ESCAPED_SUP_RE = re.compile(r'&lt;sup&gt;(.*?)&lt;/sup&gt;')


def _restore_safe_sup_tags(escaped_text):
    """
    Re-enables the one raw-HTML pattern generated output ever contains --
    `<sup>...</sup>` around a plain exponent -- after the surrounding text
    has already been fully HTML-escaped. Operates only on the escaped
    text, so anything else that happened to contain literal `<`/`>` stays
    safely escaped.
    """
    return _ESCAPED_SUP_RE.sub(r'<sup>\1</sup>', escaped_text)


def _render_table(lines):
    rows = [line.strip() for line in lines if line.strip()]
    header_cells = [c.strip() for c in rows[0].strip('|').split('|')]
    body_rows = [
        [c.strip() for c in row.strip('|').split('|')]
        for row in rows[2:]  # rows[1] is the '|---|---|' separator
    ]

    out = ['<table>', '<thead><tr>']
    for cell in header_cells:
        out.append(f'<th>{_restore_safe_sup_tags(html.escape(cell))}</th>')
    out.append('</tr></thead><tbody>')
    for row in body_rows:
        out.append('<tr>')
        for cell in row:
            out.append(f'<td>{_restore_safe_sup_tags(html.escape(cell))}</td>')
        out.append('</tr>')
    out.append('</tbody></table>')
    return ''.join(out)


def markdown_to_html(text):
    """
    Converts one `markdown_content` string (a full rendered system page)
    into an HTML fragment suitable for embedding inside a container
    element -- no `<html>`/`<body>` wrapper.

    Args:
        text (str): The Markdown source, e.g. `star_systems.markdown_content`.

    Returns:
        str: HTML markup. Empty string for `None`/empty input.
    """
    if not text:
        return ""

    blocks = re.split(r'\n\s*\n', text.strip())
    html_parts = []

    for block in blocks:
        lines = block.split('\n')

        header_match = _HEADER_RE.match(lines[0]) if len(lines) == 1 else None
        if header_match:
            level = len(header_match.group(1))
            content = _restore_safe_sup_tags(html.escape(header_match.group(2).strip()))
            html_parts.append(f'<h{level}>{content}</h{level}>')
            continue

        if (
            len(lines) >= 2
            and _TABLE_ROW_RE.match(lines[0])
            and _TABLE_SEPARATOR_RE.match(lines[1].strip())
        ):
            html_parts.append(_render_table(lines))
            continue

        paragraph = _restore_safe_sup_tags(html.escape(block)).replace('\n', '<br>')
        html_parts.append(f'<p>{paragraph}</p>')

    return '\n'.join(html_parts)
