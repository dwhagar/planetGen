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
_SLUG_INVALID_RE = re.compile(r'[^a-z0-9]+')


def _slugify(raw_text, used_ids):
    """
    Turns one heading's raw text into a unique `id` attribute value, for
    the anchors the system page's floating table-of-contents links to --
    lowercased, non-alphanumeric runs collapsed to a single `-`, and
    disambiguated with a `-2`/`-3`/... suffix against every id already
    handed out for this same document (two planets can share a name-ish
    heading, e.g. two moons both just called "I").

    Args:
        raw_text (str): The heading's text, before HTML-escaping.
        used_ids (set): Every id already returned for this conversion --
                        mutated in place to record the new one.

    Returns:
        str: A unique, non-empty id.
    """
    slug = _SLUG_INVALID_RE.sub('-', raw_text.lower()).strip('-') or 'section'
    candidate = slug
    suffix = 2
    while candidate in used_ids:
        candidate = f"{slug}-{suffix}"
        suffix += 1
    used_ids.add(candidate)
    return candidate


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


def _convert(text):
    """
    Shared implementation behind `markdown_to_html` and
    `markdown_to_html_with_headings` -- does the actual block-by-block
    conversion once, so the two never drift into assigning different ids
    to the same headings.

    Returns:
        tuple: `(html, headings)` -- `html` as described on
               `markdown_to_html`; `headings` a list of
               `{"level", "text", "id"}` dicts, one per header block, in
               document order, with `id` matching the `id="..."` embedded
               in that same heading's `<hN>` tag in `html`.
    """
    if not text:
        return "", []

    blocks = re.split(r'\n\s*\n', text.strip())
    html_parts = []
    headings = []
    used_ids = set()

    for block in blocks:
        lines = block.split('\n')

        header_match = _HEADER_RE.match(lines[0]) if len(lines) == 1 else None
        if header_match:
            level = len(header_match.group(1))
            raw_text = header_match.group(2).strip()
            content = _restore_safe_sup_tags(html.escape(raw_text))
            heading_id = _slugify(raw_text, used_ids)
            html_parts.append(f'<h{level} id="{heading_id}">{content}</h{level}>')
            headings.append({"level": level, "text": raw_text, "id": heading_id})
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

    return '\n'.join(html_parts), headings


def markdown_to_html(text):
    """
    Converts one `markdown_content` string (a full rendered system page)
    into an HTML fragment suitable for embedding inside a container
    element -- no `<html>`/`<body>` wrapper. Every heading gets an `id`
    (see `_slugify`) so it can be linked to directly; callers that need
    the heading list itself (e.g. to build a table of contents) want
    `markdown_to_html_with_headings` instead.

    Args:
        text (str): The Markdown source, e.g. `star_systems.markdown_content`.

    Returns:
        str: HTML markup. Empty string for `None`/empty input.
    """
    html_text, _headings = _convert(text)
    return html_text


def markdown_to_html_with_headings(text):
    """
    Same conversion as `markdown_to_html`, but also returns the document's
    headings for building a table-of-contents sidebar alongside the
    rendered content.

    Args:
        text (str): The Markdown source, e.g. `star_systems.markdown_content`.

    Returns:
        tuple: `(html, headings)` -- `html` as `markdown_to_html` returns;
               `headings` a list of `{"level", "text", "id"}` dicts in
               document order (empty list for `None`/empty input).
    """
    return _convert(text)
