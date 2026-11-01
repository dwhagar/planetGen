# planetgen/web/lib/citations.py

"""
Citations on the References pages (DOC.97): a citation written as a source ID becomes a link to its entry on the reference list.

Everywhere else in the project a citation is only a source ID in square brackets (`[R0123]`, `[R0123, R0456]`, `[R0123: Table 2]`; DOC.98, the numbered
style of NASA documents). The References pages, and only they, show the source by name and make it a link: `link(text, data)` turns each bracket of a
sentence into links, one for each source, to `/references/papers#paper-r0123` (the entry, with its DOI, arXiv, ADS and the document information).
The names come from `docs/references/citations.json`, which `scripts/build_reference_pages.py` builds from the reference list
(`docs/reference/references-data.json`), so a name cannot differ from the list.

A bracket that names no source, a source that the list does not hold, and an author and year written out in the text are errors, never a dead
link or a guess: the build and the tests fail. Standard library only.
"""

import html
import re

from .reference_urls import PAPER_ID, reference_path

_ID = PAPER_ID.pattern
_BRACKET = re.compile(r"\[(" + _ID + r"(?:,\s*" + _ID + r")*)(?::\s*([^\]]+))?\]")
_ONE = re.compile(_ID)
_NAME = r"[A-Z][A-Za-zÀ-ſ'’-]{2,}"
_AUTHOR_YEAR = re.compile(r"\b" + _NAME + r"(?:,? (?:and|&) " + _NAME + r"| et al\.?)?,? \(?(?:19|20)\d\d[a-z]?\b")


class CitationError(ValueError):
    """A citation on a References page that cannot become a link."""


def label(entry):
    """The words a link shows for a source, as the reference list prints them: family name, year and the letter that tells two apart."""
    return entry["label"]


def link(text, data, strict=True):
    """
    `text` as HTML with each citation bracket replaced by one link for every source in it, the rest of the text escaped. `data` is the content of
    `citations.json` (`{"sources": {"R0123": {"label", "title", "kind"}}, "aliases": {...}}`). A place inside the source (`[R0123: Table 2]`) follows its
    link after a comma. With `strict` (the build) an author and year written out in the text raises `CitationError`.
    """
    if strict:
        stripped = _BRACKET.sub("", text)
        bare = _AUTHOR_YEAR.search(stripped)
        if bare:
            raise CitationError(f"{bare.group(0)!r} is written out; write the source ID, for example [R0123]")
    out, pos = [], 0
    sources, aliases = data["sources"], data.get("aliases", {})
    for m in _BRACKET.finditer(text):
        out.append(html.escape(text[pos:m.start()]))
        links = []
        for sid in _ONE.findall(m.group(1)):
            sid = aliases.get(sid, sid)
            entry = sources.get(sid)
            if entry is None:
                raise CitationError(f"{sid} is not a source of the reference list")
            page = "data" if entry.get("kind") == "data" else "papers"
            links.append('<a class="citation" href="%s" title="%s">%s</a>' % (reference_path(page, sid), html.escape(entry.get("title") or label(entry), quote=True), html.escape(label(entry))))
        place = (", " + html.escape(m.group(2).strip())) if m.group(2) else ""
        out.append("; ".join(links) + place)
        pos = m.end()
    out.append(html.escape(text[pos:]))
    return "".join(out)


def ids_in(text):
    """Every source ID that a text cites, in order, once each."""
    seen = []
    for m in _BRACKET.finditer(text):
        for sid in _ONE.findall(m.group(1)):
            if sid not in seen:
                seen.append(sid)
    return seen
