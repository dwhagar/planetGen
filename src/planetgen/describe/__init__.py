# planetgen/describe/__init__.py

"""
One description of a star system, for every output that shows one
(docs/design/architecture/system-description.md, UX.104).

`describe(system)` turns a `StarSystem` (generated, or rebuilt from the
database) into a Document: a tree of sections whose blocks are fact tables,
rating tables, compositions and short prose, with every number a typed
quantity and no Markdown, wikitext or HTML in it. The web page, the API, the
command line and the wiki exports render that one Document.
"""

from planetgen.describe import html, markdown, text, wikitext
from planetgen.describe.build import describe
from planetgen.describe.phenomena import describe_phenomenon

DOCUMENT_FORMAT = "document"
"""str: The name the command line gives the Document itself, written as JSON, beside the formats of `RENDERERS`."""

RENDERERS = {"text": text.render, "markdown": markdown.render, "wikitext": wikitext.render}
"""dict: The formats a Document is written in on request, each to its renderer. HTML is not one: the web page calls
`html.render` itself, with the links of its own routes."""


def render(document, fmt):
    """The Document as a string in `fmt` (a key of `RENDERERS`); raises `ValueError` for any other."""
    if fmt not in RENDERERS:
        raise ValueError(f"unknown format {fmt!r} (expected one of {', '.join(RENDERERS)})")
    return RENDERERS[fmt](document)


__all__ = ["DOCUMENT_FORMAT", "RENDERERS", "describe", "describe_phenomenon", "html", "render"]
