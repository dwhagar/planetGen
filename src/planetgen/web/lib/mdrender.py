# planetgen/web/lib/mdrender.py

"""
Renders the Markdown `StarSystem.__str__` generates into HTML with the
`markdown` package (UX.39), cut down to the narrow subset that output uses:
ATX headers, GFM-style pipe tables, plain paragraphs (a single newline is a
`<br>`), and exactly one inline tag, `<sup>exponent</sup>`
(`to_scientific_notation`).

Everything else stays plain text. Names and other fields reach this text
from the CLI and the database, so no Markdown syntax other than the above
is interpreted (no emphasis, links, lists, code, quotes, raw HTML or
entities), and the serializer escapes every `<`, `>` and `&`. Headers get
unique `id`s for the system page's table of contents, and tables are wrapped
in a scrollable `div`.
"""

import re
import threading
import xml.etree.ElementTree as etree

import markdown
from markdown.blockprocessors import BlockProcessor
from markdown.extensions import Extension
from markdown.inlinepatterns import InlineProcessor
from markdown.treeprocessors import Treeprocessor

_HEADER_RE = re.compile(r'^(#{1,6})\s+(.*)$')
_SLUG_INVALID_RE = re.compile(r'[^a-z0-9]+')

_REMOVED_PREPROCESSORS = ('html_block',)
_REMOVED_BLOCKS = ('indent', 'code', 'hashheader', 'setextheader', 'hr', 'olist', 'ulist', 'quote', 'reference')
_REMOVED_INLINE = ('backtick', 'escape', 'reference', 'link', 'image_link', 'image_reference', 'short_reference',
                   'short_image_ref', 'autolink', 'automail', 'linebreak', 'html', 'entity', 'not_strong',
                   'em_strong', 'em_strong2')


def _escape_text(text):
    """Escapes `&`, `"` and `'` in the source. The serializer escapes `<` and
    `>` itself but leaves anything shaped like an entity alone, so this makes
    a literal `&lt;` in a name show as `&lt;`, and keeps quotes safe even if
    the HTML lands inside an attribute."""
    return text.replace('&', '&amp;').replace('"', '&quot;').replace("'", '&#x27;')


def _unescape_text(text):
    """Undoes `_escape_text`."""
    return text.replace('&#x27;', "'").replace('&quot;', '"').replace('&amp;', '&')


def _slugify(raw_text, used_ids):
    """
    Turns one heading's raw text into a unique `id` attribute value, for
    the anchors the system page's floating table-of-contents links to --
    lowercased, non-alphanumeric runs collapsed to a single `-`, and
    disambiguated with a `-2`/`-3`/... suffix against every id already
    handed out for this same document (two planets can share a name-ish
    heading, e.g. two moons both just called "I").

    Args:
        raw_text (str): The heading's text.
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


class _HeaderProcessor(BlockProcessor):
    """`#`..`######` followed by whitespace, as a block of its own."""

    def __init__(self, parser, extension):
        super().__init__(parser)
        self._extension = extension

    def test(self, parent, block):
        return '\n' not in block and _HEADER_RE.match(block) is not None

    def run(self, parent, blocks):
        match = _HEADER_RE.match(blocks.pop(0))
        shown = match.group(2).strip()
        # `shown` went through `_escape_text`; the heading list wants the
        # text as written.
        raw_text = _unescape_text(shown)
        level = len(match.group(1))
        heading_id = _slugify(raw_text, self._extension.used_ids)
        element = etree.SubElement(parent, f'h{level}')
        element.set('id', heading_id)
        element.text = shown
        self._extension.headings.append({"level": level, "text": raw_text, "id": heading_id})


class _SupProcessor(InlineProcessor):
    """`<sup>...</sup>` around plain text."""

    def handleMatch(self, m, data):
        element = etree.Element('sup')
        element.text = m.group(1)
        return element, m.start(0), m.end(0)


class _TableScrollProcessor(Treeprocessor):
    """Wraps each top-level table in a keyboard-scrollable `div`."""

    def run(self, root):
        for position, child in enumerate(list(root)):
            if child.tag == 'table':
                wrapper = etree.Element('div')
                wrapper.set('class', 'table-scroll')
                wrapper.set('tabindex', '0')
                wrapper.append(child)
                root[position] = wrapper


class _GeneratedMarkdown(Extension):
    """Cuts `markdown` down to the subset this module's docstring lists."""

    def __init__(self):
        super().__init__()
        self.headings = []
        self.used_ids = set()

    def extendMarkdown(self, md):
        md.registerExtension(self)
        for name in _REMOVED_PREPROCESSORS:
            md.preprocessors.deregister(name)
        for name in _REMOVED_BLOCKS:
            md.parser.blockprocessors.deregister(name)
        for name in _REMOVED_INLINE:
            md.inlinePatterns.deregister(name)
        md.parser.blockprocessors.register(_HeaderProcessor(md.parser, self), 'generated_header', 70)
        md.inlinePatterns.register(_SupProcessor(r'<sup>(.*?)</sup>', md), 'generated_sup', 190)
        md.treeprocessors.register(_TableScrollProcessor(md), 'generated_table_scroll', 5)

    def reset(self):
        self.headings = []
        self.used_ids = set()


_local = threading.local()


def _converter():
    """This thread's converter (building one per call is the slow part)."""
    converter = getattr(_local, 'converter', None)
    if converter is None:
        extension = _GeneratedMarkdown()
        # `use_align_attribute`: the default `style="text-align:..."` would
        # be blocked by the page's Content-Security-Policy.
        converter = markdown.Markdown(
            extensions=[extension, 'tables', 'nl2br'],
            extension_configs={'tables': {'use_align_attribute': True}},
            output_format='html',
        )
        converter.generated = extension
        _local.converter = converter
    return converter


def _convert(text):
    """
    Shared implementation behind `markdown_to_html` and
    `markdown_to_html_with_headings`, so the two never assign different ids
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
    converter = _converter()
    converter.reset()
    html_text = converter.convert(_escape_text(text.strip()))
    return html_text, list(converter.generated.headings)


def markdown_to_html(text):
    """
    Converts one Markdown string (a rendered system page, or one piece of it)
    into an HTML fragment suitable for embedding inside a container
    element -- no `<html>`/`<body>` wrapper. Every heading gets an `id`
    (see `_slugify`) so it can be linked to directly; callers that need
    the heading list itself (e.g. to build a table of contents) want
    `markdown_to_html_with_headings` instead.

    Args:
        text (str): The Markdown source, e.g. from `planetgen.db.render`.

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
        text (str): The Markdown source, e.g. from `planetgen.db.render`.

    Returns:
        tuple: `(html, headings)` -- `html` as `markdown_to_html` returns;
               `headings` a list of `{"level", "text", "id"}` dicts in
               document order (empty list for `None`/empty input).
    """
    return _convert(text)
