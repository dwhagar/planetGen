# planetgen/describe/__init__.py

"""
One description of a star system, for every output that shows one
(docs/design/system-description.md, UX.104).

`describe(system)` turns a `StarSystem` (generated, or rebuilt from the
database) into a Document: a tree of sections whose blocks are fact tables,
rating tables, compositions and short prose, with every number a typed
quantity and no Markdown, wikitext or HTML in it. The web page, the API, the
command line and the wiki exports render that one Document.
"""

from planetgen.describe.build import describe

__all__ = ["describe"]
