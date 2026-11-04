# planetgen/web/lib/references_site.py

"""
The text of the References pages (DOC.88, References R2): the one section of the site that speaks out of universe.

`/references` (Help > References) says where each number came from, what is a stand-in, what was fitted to a result and how sure each
part is. Its pages are static: the numbers are read from `docs/references/*.json`, which `scripts/build_reference_pages.py` builds from
the registries (the reference list and the stand-in registry), so a count on a page cannot differ from the registry. `page()` builds the
start page; the pages of the later items (papers, data, stand-ins, method, speculative chemistry) add their own builders here.

The pages are held to the plain-words guide like every page (short sentences, American spelling, the Oxford comma); only the
out-of-universe words are allowed here (`plainwords.exempt`). Standard library only.
"""

import functools
import html
import json
import os
import re
import urllib.parse

from markupsafe import Markup

from . import citations, references_usage
from .reference_urls import anchor, paper_anchor, reference_path

DOCS = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..", "docs", "references"))
"""str: The folder the build writes (`docs/references/`)."""

BANNER = "Out of universe"
"""str: The words of the banner every References page opens with; the banner's id is `out-of-universe-banner`."""

BANNER_TEXT = ("You have left the story. Everywhere else the site speaks as a working navigation system that simply knows things. "
               "This section says how the numbers were made.")
"""str: The sentence after the banner words."""


@functools.lru_cache(maxsize=8)
def _read(path, stamp):
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def load(name, folder=None):
    """
    The data of one built file (`overview.json`), or `None` when the build has not been run. The file is read again when it changes.
    """
    path = os.path.join(folder or DOCS, name)
    try:
        stamp = os.stat(path).st_mtime_ns
    except OSError:
        return None
    return _read(path, stamp)


def cite(text, folder=None, strict=True):
    """
    A sentence of a References page as HTML, each citation ID in it (`[R0123]`) replaced by a link to the source's entry (DOC.97). Only these pages
    name a source; everywhere else a citation stays an ID. A source the list does not hold, or a name written out, raises `citations.CitationError`.
    `strict=False` skips the check for a name written out: the registry rows (DOC.90) carry decisions ("Boss, 2026-10-21") that are not citations,
    and the build has already refused an unknown ID in them.
    """
    data = load("citations.json", folder)
    if data is None:
        raise citations.CitationError("docs/references/citations.json is missing: run python scripts/build_reference_pages.py")
    return citations.link(text, data, strict)


def count(number):
    """A whole number with thousands separators (`2,444`)."""
    return f"{number:,}"


def _percent_row(label, found, meaning):
    return [label, count(found["rows"]), f"{found['percent']}%", meaning]


def overview(data):
    """
    The start page, `/references`, as `{"title", "intro", "sections"}`. Each section is `{"id", "heading", "paragraphs", "table"}`; a
    table is `{"head", "rows"}` and may be `None`. `data` is `overview.json`, or `None` before the build has run.
    """
    page = {
        "title": "References",
        "intro": [
            "This section says where the numbers on the other pages come from. It lists the papers and data sets behind them, "
            "the values that stand in for a measurement, and how sure each part is.",
            "Each stand-in has a source tier, from verified to chosen for the output. It also has a tweak level, "
            "from used as its source gives it to fitted with no source.",
        ],
        "sections": [],
        "pages": [{"text": "Papers", "href": reference_path("papers")}, {"text": "Data sets and archives", "href": reference_path("data")}],
    }
    if data is None:
        page["sections"].append({
            "id": "not-built", "heading": "The pages are not built", "table": None,
            "paragraphs": ["The numbers on this page are built at setup by scripts/build_reference_pages.py. Run it, or run update.sh, and reload."]})
        return page
    papers, sources, stand_ins, verdict = data["papers"], data["data"], data["stand_ins"], data["verdict"]
    page["sections"].append({
        "id": "standing", "heading": "Where the numbers stand", "table": None,
        "paragraphs": [
            f"The registry holds {count(stand_ins['rows'])} stand-ins: values that are not a verified measurement of the thing they stand for. "
            f"{verdict['read_from_a_source']['percent']}% are read from a source we trust, "
            f"and {verdict['unverified']['percent']}% come from a source we doubt or from none.",
            f"{verdict['educated_guess_from_physics']['percent']}% are an educated guess from physics, "
            f"and {verdict['chosen_for_the_output']['percent']}% are chosen so the output comes out right.",
            f"{verdict['fitted_to_a_result']['percent']}% were bent toward a result. "
            f"They were tuned to a published statistic, tuned against a sweep, or fitted with no source."]})
    page["sections"].append({
        "id": "verdict", "heading": "How sure the stand-ins are", "paragraphs": [
            "Every stand-in is counted once. A row that is a family of values counts at the tier of its first part."],
        "table": {"head": ["Kind of value", "Stand-ins", "Share", "What it means"], "rows": [
            _percent_row("Read from a source", verdict["read_from_a_source"], "tiers T1 and T2: verified, or from a source we trust"),
            _percent_row("Unverified", verdict["unverified"], "tiers T3 and T4: a source we doubt, or no source and a reason for the range"),
            _percent_row("Educated guess from physics", verdict["educated_guess_from_physics"], "tiers T5 to T8: physics fixes the form"),
            _percent_row("Chosen for the output", verdict["chosen_for_the_output"], "tier T9: chosen so the result is right"),
            _percent_row("Fitted to a result", verdict["fitted_to_a_result"], "tweak levels L2 to L4, of any tier"),
        ]}})
    page["sections"].append({
        "id": "sources", "heading": "What the numbers are read from", "paragraphs": [
            f"The reference list holds {count(papers['entries'])} papers, books, and data sets. {count(papers['with_doi'])} have a DOI page "
            f"and {count(papers['with_arxiv'])} have an arXiv copy. {count(papers['without_doi'])} have no DOI.",
            f"{count(sources['downloaded_at_setup'])} data sets are downloaded when the site is set up, and {count(sources['archives'])} "
            f"archives and tools are named."],
        "table": None})
    page["sections"].append({
        "id": "tiers", "heading": "The source tiers", "paragraphs": [
            "A tier says where a value came from, from the best source to a value chosen for the output."],
        "table": {"head": ["Tier", "Name", "Stand-ins", "Meaning"], "rows": [
            [code, name, count(stand_ins["by_tier"][code]), meaning] for code, name, meaning in data["tiers"]]}})
    page["sections"].append({
        "id": "tweaks", "heading": "The tweak levels", "paragraphs": [
            "A tweak level says how far a value was bent from its source, from not at all to fitted with no source."],
        "table": {"head": ["Level", "Name", "Stand-ins", "Meaning"], "rows": [
            [code, name, count(stand_ins["by_tweak"][code]), meaning] for code, name, meaning in data["tweaks"]]}})
    return page


# ------------------------------------------------------------------------------------------------------ papers and data (DOC.89)

PAPER_VIEWS = (("az", "A to Z"), ("area", "By area"), ("field", "By field"), ("year", "By year"), ("type", "By type"), ("totals", "Totals"))
"""tuple: The views of the Papers page and their names; the first is the default and lists every source once, with its anchor."""

UNUSED = "none"
UNUSED_TITLE = "Not cited in code, data, or notes"
"""The group of the sources that no code file, data file, registry row, or design note cites (reports and older notes may)."""

SHOWN = 8
"""int: How many places of one kind an entry lists before it says how many more there are."""

MOST_CITED = 15
SEARCH_LIMIT = 300
ADS = "https://ui.adsabs.harvard.edu/abs/"


def paper_links(row):
    """
    The links of one source, `[{"label", "url", "kind"}]`, the one to follow first at the head: its DOI page, else its arXiv page, else the link
    the record carries, else a web search for its title and first author. A source always has at least one (a test fails when it does not). The
    ADS page of a DOI or an arXiv id follows.
    """
    links = []
    if row.get("doi"):
        links.append({"label": "DOI", "url": "https://doi.org/" + urllib.parse.quote(row["doi"], safe="/()"), "kind": "doi"})
    if row.get("arxiv"):
        links.append({"label": "arXiv", "url": "https://arxiv.org/abs/" + row["arxiv"], "kind": "arxiv"})
    if row.get("link"):
        label = "ADS search" if "adsabs" in row["link"] else "Link"
        links.append({"label": label, "url": row["link"], "kind": "ads" if "adsabs" in row["link"] else "link"})
    if row.get("doi"):
        links.append({"label": "ADS", "url": ADS + urllib.parse.quote(row["doi"], safe="/()") + "/abstract", "kind": "ads"})
    elif row.get("arxiv"):
        links.append({"label": "ADS", "url": ADS + "arXiv:" + row["arxiv"] + "/abstract", "kind": "ads"})
    if not links:
        words = " ".join(part for part in (row.get("title", ""), row.get("who", "").split(",")[0]) if part)
        links.append({"label": "Search", "url": "https://www.google.com/search?q=" + urllib.parse.quote(words), "kind": "search"})
    return links


def _place(text, strip):
    return text[len(strip):] if text.startswith(strip) else text


def _used(row_id, usage):
    entry = usage["uses"].get(row_id)
    return entry or {"values": [], "notes": [], "stand_ins": [], "areas": []}


def _listed(names):
    names = list(names)
    return {"names": names[:SHOWN], "more": max(0, len(names) - SHOWN)}


@functools.lru_cache(maxsize=2)
def _keys_of(stamp):
    return frozenset(row["key"] for _part, row in rows_of(load("stand-ins.json")))


def _stand_in_keys():
    """The registry keys that have a page, read again when `stand-ins.json` changes; empty before the build has run."""
    try:
        stamp = os.stat(os.path.join(DOCS, "stand-ins.json")).st_mtime_ns
    except OSError:
        return frozenset()
    return _keys_of(stamp)


def _stand_in_places(keys):
    """
    The stand-ins that cite a source, as `{"names": [{"text", "url"}], "more": number}`: each name is the row ID and links to the page of the row (DOC.90)
    when the stand-in pages are built and hold that row (the registry key is `5.1:GR1`).
    """
    known = _stand_in_keys()
    keys = list(keys)
    names = [{"text": key.split(":", 1)[-1], "url": reference_path("stand-in", key) if key in known else ""} for key in keys[:SHOWN]]
    return {"names": names, "more": max(0, len(keys) - SHOWN)}


def paper_entry(row, usage, anchor_id=True, compact=False):
    """
    One source as the page shows it: its ID, name, authors, title, where it was published, its links, and the places that cite it. A compact entry
    (the A to Z list, which holds every source so that each citation link finds its anchor) keeps the first link and the areas, and leaves out the
    lists of places; the search for its ID gives the whole entry.
    """
    used = _used(row["id"], usage)
    links = paper_links(row)
    return {
        "compact": compact,
        "search": _paper_href("az", q=row["id"]),
        "anchor": paper_anchor(row["id"]) if anchor_id else "",
        "id": row["id"],
        "label": row["label"],
        "who": row.get("who", ""),
        "year": row.get("year") or "",
        "title": row.get("title", ""),
        "venue": row.get("venue", ""),
        "type": row["type"],
        "field": row["field"],
        "links": links[:1] if compact else links,
        "areas": [references_usage.AREA_TITLES[area] for area in used["areas"]],
        "code": _listed(_place(path, "src/planetgen/") for path in used["values"]),
        "stand_ins": _stand_in_places(used["stand_ins"]),
        "notes": _listed(_place(path, "docs/design/") for path in used["notes"]),
        "cited_by": references_usage.cited_by(used) if used["areas"] else 0,
    }


def compact_html(entry):
    """
    A compact entry as the HTML of one list item. The A to Z list holds every source, so it is written here and not by the template, which would
    take seconds to fill five thousand items.
    """
    escape = lambda text: html.escape(str(text), quote=True)
    first = entry["links"][0]
    used = (". Used in " + escape(", ".join(entry["areas"]))) if entry["areas"] else ""
    anchor_id = f' id="{entry["anchor"]}"' if entry["anchor"] else ""
    venue = (". " + escape(entry["venue"])) if entry["venue"] else ""
    return (f'<li class="reference-entry reference-compact"{anchor_id}><strong>{escape(entry["label"])}</strong> '
            f'<a class="reference-id" href="{escape(entry["search"])}">{entry["id"]}</a>. {escape(entry["title"])}{venue}. '
            f'<a href="{escape(first["url"])}" rel="noopener noreferrer">{escape(first["label"])}</a>{used}.</li>')


def _slug(name):
    return anchor(name)


def paper_groups(rows, usage, axis):
    """
    `{group slug: (name, [rows])}` for one axis of the Papers page (`area`, `field`, `year`, `type`). A source cited in two areas is in both; a
    source that nothing cites is in the group `none` of the area axis.
    """
    groups = {}
    for row in rows:
        if axis == "area":
            areas = _used(row["id"], usage)["areas"] or [UNUSED]
            names = [(area, references_usage.AREA_TITLES.get(area, UNUSED_TITLE)) for area in areas]
        elif axis == "field":
            names = [(_slug(row["field"]), row["field"])]
        elif axis == "year":
            decade = f"{row['year'] // 10 * 10}s" if row.get("year") else "Not stated"
            names = [(_slug(decade), decade)]
        else:
            names = [(_slug(row["type"]), row["type"])]
        for slug, name in names:
            groups.setdefault(slug, (name, []))[1].append(row)
    return groups


def _ordered(groups, axis):
    if axis == "area":
        order = [area for area in list(references_usage.AREA_TITLES) if area in groups] + ([UNUSED] if UNUSED in groups else [])
    elif axis == "year":
        order = sorted(groups, key=lambda slug: (slug == "not-stated", slug))
    else:
        order = sorted(groups, key=lambda slug: (-len(groups[slug][1]), groups[slug][0]))
    return order


def _link(text, href):
    return {"text": text, "href": href}


def _paper_href(view, group=None, q=None):
    query = {}
    if view != "az":
        query["view"] = view
    if group:
        query["group"] = group
    if q:
        query["q"] = q
    return reference_path("papers") + (("?" + urllib.parse.urlencode(query)) if query else "")


def _group_table(rows, usage, axis, name):
    groups = paper_groups(rows, usage, axis)
    body = []
    for slug in _ordered(groups, axis):
        title, members = groups[slug]
        body.append({"id": "", "cells": [_link(title, _paper_href(axis, slug)), count(len(members)),
                                         count(sum(1 for row in members if row.get("doi"))), count(sum(1 for row in members if row.get("arxiv")))]})
    return {"heading": "", "head": [name, "Sources", "With a DOI", "With an arXiv copy"], "rows": body}


def _number_rows(pairs):
    return [{"id": "", "cells": [name, count(number)]} for name, number in pairs]


def _source_cell(item):
    return _link(f"{item['label']}", reference_path("papers", item["id"]))


def totals_sections(data, usage):
    """The Totals view: the reference list in numbers, the sources we cite most, and the counts by area, field, decade, and type."""
    totals, rows = data["totals"], data["rows"]
    cited = [row for row in rows if _used(row["id"], usage)["areas"]]
    ranked = sorted(cited, key=lambda row: (-references_usage.cited_by(_used(row["id"], usage)), row["id"]))[:MOST_CITED]
    numbers = [
        ["Papers, books, and data sets", count(totals["papers"])],
        ["With a DOI page", count(totals["with_doi"])],
        ["With an arXiv copy", count(totals["with_arxiv"])],
        ["Books", count(totals["books"])],
        ["Data sets with a DOI", count(totals["datasets"])],
        ["Cited in code, data, registry rows, or design notes", count(len(cited))],
        ["Oldest", _link(totals["oldest"]["label"], reference_path("papers", totals["oldest"]["id"]))],
        ["Newest", _link(totals["newest"]["label"], reference_path("papers", totals["newest"]["id"]))],
    ]
    sections = [
        {"id": "numbers", "heading": "The reference list in numbers", "paragraphs": [
            "Every count on this page is made from the reference list and from the files that cite it, when the page is built or asked for."],
         "tables": [{"heading": "", "head": ["Count", "Number"], "rows": [{"id": "", "cells": row} for row in numbers]}]},
        {"id": "most-cited", "heading": "Most cited by us", "paragraphs": [
            "These sources are cited in the most places: code files, data files, stand-in rows, and design notes."],
         "tables": [{"heading": "", "head": ["Source", "Year", "Cited in", "Areas"], "rows": [
             {"id": "", "cells": [_source_cell(row), str(row.get("year") or ""), count(references_usage.cited_by(_used(row["id"], usage))),
                                  ", ".join(references_usage.AREA_TITLES[a] for a in _used(row["id"], usage)["areas"])]} for row in ranked]}]},
    ]
    for axis, name, heading in (("area", "Area", "Where we use them"), ("field", "Field", "The field of the journal"),
                                ("year", "Decade", "The decade of the paper"), ("type", "Type", "The type of source")):
        sections.append({"id": f"by-{axis}", "heading": heading, "paragraphs": [], "tables": [_group_table(rows, usage, axis, name)]})
    return sections


def _matches(row, words):
    text = " ".join(str(row.get(key, "")) for key in ("id", "label", "who", "title", "doi", "arxiv", "venue")).lower()
    return all(word in text for word in words)


def _letter(row):
    first = (row.get("label") or "#")[:1].upper()
    return first if first.isalpha() and first.isascii() else "#"


def papers_page(data, usage, view="az", group="", q=""):
    """
    The Papers page, `/references/papers`, as a dict the template shows: `title`, `intro`, `views` (the tabs), `query`, `sections` (tables), `blocks`
    (lists of entries), `note`. `view` is one of `PAPER_VIEWS`; `group` names one group of an axis; `q` is a search over the title, authors, DOI,
    arXiv id, and ID. A source has an anchor (`paper-r0123`) in the first place the page lists it, so a citation link always lands on it.
    """
    view = view if view in dict(PAPER_VIEWS) else "az"
    page = {
        "title": "Papers",
        "intro": [
            "This is every paper, book, and data set that the numbers, laws, and notes of the project cite. Each source has an ID, such as R0123, "
            "that the other pages and the code use.",
            "A source links to its DOI page, else its arXiv page, else a search. An area is read from the name of the file that cites a source, "
            "and a source cited in two areas is in both. The field is read from the journal.",
        ],
        "views": [], "query": q, "sections": [], "blocks": [], "note": "", "view": view, "count": 0,
    }
    for key, name in PAPER_VIEWS:
        page["views"].append({"key": key, "label": name, "href": _paper_href(key), "current": key == view and not q})
    if data is None:
        page["sections"].append({"id": "not-built", "heading": "The pages are not built", "tables": [],
                                 "paragraphs": ["The papers are built at setup by scripts/build_reference_pages.py. Run it, or run update.sh, and reload."]})
        return page
    rows, seen = data["rows"], set()

    def block(heading, members, block_id="", compact=False):
        entries = []
        for row in members:
            entries.append(paper_entry(row, usage, anchor_id=row["id"] not in seen, compact=compact))
            seen.add(row["id"])
        if compact:
            return {"id": block_id, "heading": heading, "entries": [], "html": Markup("".join(compact_html(entry) for entry in entries))}
        return {"id": block_id, "heading": heading, "entries": entries, "html": ""}

    if q.strip():
        words = q.lower().split()
        found = [row for row in rows if _matches(row, words)]
        page["count"] = len(found)
        page["note"] = (f"{count(len(found))} sources match." if found else "No source matches.") + (
            f" The first {SEARCH_LIMIT} are shown." if len(found) > SEARCH_LIMIT else "")
        page["blocks"].append(block("", found[:SEARCH_LIMIT]))
        return page
    if view == "totals":
        page["sections"] = totals_sections(data, usage)
        return page
    if view == "az":
        page["count"] = len(rows)
        page["letters"] = []
        by_letter = {}
        for row in rows:
            by_letter.setdefault(_letter(row), []).append(row)
        for letter in sorted(by_letter, key=lambda l: (l == "#", l)):
            page["letters"].append({"letter": letter, "id": f"letter-{letter.lower() if letter != '#' else 'other'}"})
            page["blocks"].append(block(letter if letter != "#" else "Other", by_letter[letter], page["letters"][-1]["id"], compact=True))
        page["note"] = f"{count(len(rows))} sources."
        return page
    axis = view
    groups = paper_groups(rows, usage, axis)
    if group in groups:
        name, members = groups[group]
        page["count"] = len(members)
        page["note"] = f"{count(len(members))} sources in {name}."
        page["blocks"].append(block(name, members, f"group-{group}"))
        return page
    names = {"area": "Area", "field": "Field", "year": "Decade", "type": "Type"}
    page["sections"].append({"id": f"by-{axis}", "heading": dict(PAPER_VIEWS)[view], "tables": [_group_table(rows, usage, axis, names[axis])],
                             "paragraphs": ["Choose a group to list its sources."]})
    return page


LICENCE_FILE = os.path.join(references_usage.ROOT, "DATA-LICENSE.md")
"""str: The file that records where each data file came from and what is known of its terms (read by the Data page, never copied)."""

_MARKUP = re.compile(r"`|\*\*|\\(?=\|)")
_MD_LINK = re.compile(r"\[([^\]]+)\]\([^)]*\)")


def _plain_cell(cell):
    return _MD_LINK.sub(r"\1", _MARKUP.sub("", cell)).strip()


def _american(text):
    return re.sub(r"\bLicence\b", "License", re.sub(r"\blicence\b", "license", text))


def licence_tables(text):
    """
    The tables of `DATA-LICENSE.md` as `[{"heading", "head", "rows"}]`: the heading above each table (backticks removed), its column names, and
    its rows with the same number of cells as the columns. Text between the tables is not read.
    """
    tables, heading, lines, index = [], "", text.splitlines(), 0
    while index < len(lines):
        line = lines[index]
        if line.startswith("#"):
            heading = _american(_plain_cell(line.lstrip("#")))
        elif line.startswith("|") and index + 1 < len(lines) and re.match(r"^\|\s*:?-{3,}", lines[index + 1]):
            head = [_american(_plain_cell(cell)) for cell in line.strip().strip("|").split("|")]
            rows, index = [], index + 2
            while index < len(lines) and lines[index].startswith("|"):
                cells = [_plain_cell(cell) for cell in lines[index].strip().strip("|").split("|")]
                rows.append((cells + [""] * len(head))[:len(head)])
                index += 1
            tables.append({"heading": heading, "head": head, "rows": rows})
            continue
        index += 1
    return tables


def licence_text(path=None):
    """The text of `DATA-LICENSE.md`, or `""` when the file is not there."""
    try:
        with open(path or LICENCE_FILE, encoding="utf-8") as handle:
            return handle.read()
    except OSError:
        return ""


def _setup_terms(tables):
    """`{dataset key: terms recorded}` from the table of `DATA-LICENSE.md` that lists the data sets fetched at setup."""
    terms = {}
    for table in tables:
        if table["head"][:1] == ["Dataset"] and "setup" in table["heading"].lower():
            for name, recorded, *_rest in table["rows"]:
                terms[re.match(r"[A-Za-z0-9_]+", name).group(0) if re.match(r"[A-Za-z0-9_]+", name) else name] = recorded
    return terms


def _row(source_id, cells, seen):
    """A row of a table of the Data page; two data sets may share one source ID, and only the first carries its anchor."""
    first = source_id not in seen
    seen.add(source_id)
    return {"id": paper_anchor(source_id) if first else "", "cells": cells}


def data_page(data, licences=""):
    """
    The Data page, `/references/data`, as a dict the template shows: the data sets fetched at setup, the archives and catalogs, the standards and
    web sources, the data sets with a DOI, and the terms of use that `DATA-LICENSE.md` records. Each source has a row with the anchor
    `paper-r0123`, so a citation link lands on it.
    """
    page = {
        "title": "Data sets and archives",
        "intro": [
            "This page lists the data behind the numbers. It covers the data sets the site downloads when it is set up, the archives and catalogs "
            "that were read, and the standards and web pages that were cited. It also gives what is known of the terms of use.",
            "The raw data sets are never stored in the project and are never passed on. Only the values read from them are kept.",
        ],
        "views": [], "query": "", "sections": [], "blocks": [], "note": "",
    }
    if data is None:
        page["sections"].append({"id": "not-built", "heading": "The pages are not built", "tables": [],
                                 "paragraphs": ["The data is built at setup by scripts/build_reference_pages.py. Run it, or run update.sh, and reload."]})
        return page
    tables = licence_tables(licences)
    terms, seen = _setup_terms(tables), set()
    page["sections"] = [
        {"id": "setup", "heading": "Data sets downloaded at setup", "paragraphs": [
            f"{count(len(data['setup']))} data sets are downloaded when the site is set up. A step keeps the values it needs and removes the rest."],
         "tables": [{"heading": "", "head": ["ID", "Data set", "Source", "Terms recorded"], "rows": [
             _row(row["id"], [row["id"], _link(row["name"], row["link"]), row["source"], terms.get(row["key"], "")], seen) for row in data["setup"]]}]},
        {"id": "archives", "heading": "Archives, catalogs, and tools", "paragraphs": [
            f"{count(len(data['archives']))} archives, catalogs, and tools are named in the notes and the code. The last column counts the places that name each one."],
         "tables": [{"heading": "", "head": ["ID", "Name", "Provider", "Places"], "rows": [
             _row(row["id"], [row["id"], _link(row["name"], row["link"]), row["provider"], count(row["mentions"])], seen) for row in data["archives"]]}]},
        {"id": "standards", "heading": "Standards and web pages", "paragraphs": [
            f"{count(len(data['web']))} standards, regulations, reports, and web pages are cited. Each was opened on the date shown."],
         "tables": [{"heading": "", "head": ["ID", "Name", "Provider", "Kind", "Identifier", "Year", "Opened"], "rows": [
             _row(row["id"], [row["id"], _link(row["name"], row["link"]), row["provider"], row["kind"], row["identifier"], row["year"], row["checked"]], seen)
             for row in data["web"]]}]},
        {"id": "doi-data", "heading": "Data sets that have a DOI", "tables": [], "paragraphs": [
            "A data set that has its own DOI is listed with the papers, because it is cited the same way. The Papers page lists them under the type Data set."],
         "links": [_link("Data sets on the Papers page", _paper_href("type", "data-set"))]},
    ]
    if tables:
        page["sections"].append({
            "id": "terms", "heading": "Terms of use", "paragraphs": [
                "The program's own code and text are released under CC0. The tables below come from DATA-LICENSE.md. They record where each data file "
                "came from and what is known of its terms. Where no statement was found, the entry says so."],
            "tables": [{"heading": table["heading"], "head": table["head"], "rows": [{"id": "", "cells": row} for row in table["rows"]]} for table in tables]})
    page["sections"].append({"id": "thanks", "heading": "With thanks", "tables": [], "paragraphs": [
        "This project would not be possible without the researchers, observers, instrument teams, and data archives whose work it builds on. "
        "Their papers are on the Papers page."]})
    return page

# --- The stand-in pages (DOC.90, References R4) -------------------------------------------------------------------------------------

FITTED_TWEAKS = ("L2", "L3", "L4")
"""tuple: The tweak levels that were bent toward a result; the Fitted to a result page lists these rows."""

SUMMARY_WORDS = 28
"""int: How many words of a row's description the index shows."""

_BRACKET = re.compile(r"\s*\[R\d{4,}(?:,\s*R\d{4,})*(?::[^\]]*)?\]")

_UNREAD = {
    "by design": "Nothing in the app reads this value, by design: it is a record of a choice and not an input.",
}


def summary(text, words=SUMMARY_WORDS):
    """The first `words` words of a description with its citations left out, and an ellipsis where it was cut."""
    plain = _BRACKET.sub("", text).split()
    return " ".join(plain[:words]) + (" ..." if len(plain) > words else "")


def _lookup(table, code):
    """The `(name, meaning)` of a tier or tweak code in the legend rows `[[code, name, meaning]]`."""
    return next(((name, meaning) for found, name, meaning in table if found == code), ("", ""))


def rows_of(data):
    """Every row of `stand-ins.json` as `(section, row)`, in the order the registry lists them."""
    return [(section, row) for section in data["sections"] for row in section["rows"]]


def find(data, slug):
    """The `(section, row)` whose page address is `slug`, or `None`."""
    for section, row in rows_of(data):
        if row["slug"] == slug:
            return section, row
    return None


def _filtered(data, tier=None, tweak=None, section=None, text=None):
    """The sections of `data` cut down to the rows that match every filter given; a section with no row left is dropped."""
    words = (text or "").lower().split()
    kept = []
    for part in data["sections"]:
        if section and section not in (part["slug"], part["number"], part["key"]):
            continue
        rows = []
        for row in part["rows"]:
            if tier and row["tier_code"] != tier:
                continue
            if tweak and row["tweak"] != tweak:
                continue
            if words:
                hay = " ".join([row["id"], row["what"], part["title"]] + [value for _label, value in row["fields"]]).lower()
                if not all(word in hay for word in words):
                    continue
            rows.append(row)
        if rows:
            kept.append(dict(part, rows=rows))
    return kept


def standins(data, overview, tier=None, tweak=None, section=None, text=None):
    """
    The registry page, `/references/stand-ins`, as `{"title", "intro", "filters", "found", "total", "sections", "tiers", "tweaks"}`. `data` is
    `stand-ins.json` and `overview` is `overview.json`; the filters narrow the rows to one tier, one tweak level, one section, or the words
    in a text. Every row is `{"slug", "id", "what", "tier", "tweak", "kind"}` and links to its own page.
    """
    kept = _filtered(data, tier, tweak, section, text)
    sections = [{"slug": part["slug"], "number": part["number"], "title": part["title"], "rows": [
        {"slug": row["slug"], "url": reference_path("stand-in", row["key"]), "id": row["id"], "what": summary(row["what"]), "tier": row["tier"], "tweak": row["tweak"],
         "kind": data["kinds"].get(row["kind"], row["kind"])} for row in part["rows"]]} for part in kept]
    found = sum(len(part["rows"]) for part in sections)
    total = sum(len(part["rows"]) for part in data["sections"])
    return {
        "title": "Stand-ins",
        "intro": [
            f"This is the registry of every stand-in: {count(total)} values that are not a verified measurement of the thing they stand for. "
            "Each has its own page that says what it is, how it was made, and what would replace it.",
            "A tier says where the number came from. A tweak level says how far it was bent. The legends are on the start page and below.",
        ],
        "filters": {"tier": tier or "", "tweak": tweak or "", "section": section or "", "text": text or ""},
        "found": found, "total": total, "sections": sections,
        "tiers": [list(row) for row in overview["tiers"]], "tweaks": [list(row) for row in overview["tweaks"]],
        "tier_options": [{"value": "", "label": "Any tier"}] + [{"value": code, "label": f"{code}: {name}"} for code, name, _meaning in overview["tiers"]],
        "tweak_options": [{"value": "", "label": "Any tweak level"}] + [{"value": code, "label": f"{code}: {name}"} for code, name, _meaning in overview["tweaks"]],
        "section_options": [{"value": "", "label": "Every section"}] + [{"value": part["slug"], "label": f"{part['number']} {part['title']}"}
                                                                         for part in data["sections"]],
    }


def _near(row):
    """The next or the last stand-in of a section, for the links at the foot of a page."""
    return {"id": row["id"], "url": reference_path("stand-in", row["key"])}


def standin(data, overview, slug):
    """
    The page of one stand-in, `/references/stand-ins/<slug>`, or `None` when the registry has no such row: `{"title", "section", "row", "head",
    "paragraphs", "fields", "before", "after"}`. Each field is `[label, html]`; the html links its citations.
    """
    found = find(data, slug)
    if found is None:
        return None
    section, row = found
    rows = section["rows"]
    at = rows.index(row)
    tier_name, tier_meaning = _lookup(overview["tiers"], row["tier_code"])
    tweak_name, tweak_meaning = _lookup(overview["tweaks"], row["tweak"])
    kind = data["kinds"].get(row["kind"], row["kind"])
    paragraphs = [
        f"This is stand-in {row['id']}, row {row['key']} of the registry.",
        f"Its source tier is {row['tier']}: {tier_name.lower()}. {tier_meaning[:1].upper()}{tier_meaning[1:]}." if tier_name else "",
        f"Its tweak level is {row['tweak']}: {tweak_name.lower()}. {tweak_meaning[:1].upper()}{tweak_meaning[1:]}.",
    ]
    if kind:
        paragraphs.append(f"It is a {kind} fill.")
    if row["unread"] == "by design":
        paragraphs.append(_UNREAD["by design"])
    elif row["unread"]:
        paragraphs.append(f"Nothing in the app reads this value yet. The item {row['unread']} will wire it.")
    if row["burn_down"]:
        disposition, item, reason = row["burn_down"]
        lead = {"no data": "Nothing we hold or can reach measures it, so it stays as it is.",
                "replace": f"The item {item} replaces it with data.",
                "replaced": "It has been replaced by data; this line stays for one release so the change can be found."}[disposition]
        paragraphs.append(f"{lead} {reason[:1].upper()}{reason[1:]}.")
    return {
        "title": f"Stand-in {row['id']}", "section": {"slug": section["slug"], "number": section["number"], "title": section["title"]},
        "row": row, "head": cite(row["what"], strict=False), "paragraphs": [text for text in paragraphs if text],
        "fields": [[label, cite(text, strict=False)] for label, text in row["fields"]],
        "before": _near(rows[at - 1]) if at else None, "after": _near(rows[at + 1]) if at + 1 < len(rows) else None,
        "all": reference_path("stand-in"),
    }


def fitted(data, targets, overview):
    """
    The Fitted to a result page, `/references/fitted`: every stand-in at a tweak level of L2 to L4, grouped by level, and the table of tuning targets
    the L2 rows were tuned to. `{"title", "intro", "levels", "targets", "total"}`.
    """
    levels = []
    for code in FITTED_TWEAKS:
        name, meaning = _lookup(overview["tweaks"], code)
        rows = [{"slug": row["slug"], "url": reference_path("stand-in", row["key"]), "id": row["id"], "section": f"{part['number']} {part['title']}", "what": summary(row["what"]), "tier": row["tier"]}
                for part, row in rows_of(data) if row["tweak"] == code]
        levels.append({"code": code, "name": name, "meaning": meaning, "rows": rows})
    return {
        "title": "Fitted to a result",
        "intro": [
            "A stand-in is fitted when it was bent toward a result. Level L2 was tuned until a published count or distribution came out. "
            "Level L3 was chosen by running the generator over many draws and looking at the outcome. Level L4 is whatever made the output right, with no independent source.",
            f"{count(sum(len(level['rows']) for level in levels))} stand-ins are at these levels. The tuning targets they were tuned to are listed last.",
        ],
        "levels": levels,
        "targets": [dict(target, statistic=cite(target["statistic"] or "", strict=False), source=cite(target["source"] or "", strict=False),
                         falsifier=cite(target["falsifier"] or "", strict=False))
                    for target in (targets["targets"] if targets else [])],
        "total": sum(len(level["rows"]) for level in levels),
    }


def _tier_name(data, code):
    return next(name for tier, name, _meaning in data["tiers"] if tier == code)


def _level_name(data, code):
    return next(name for level, name, _meaning in data["tweaks"] if level == code)


def speculative(data, folder=None):
    """
    `/references/speculative` (DOC.92, References R6) as `{"title", "intro", "sections"}`, in the shape of `overview`: how sure the air
    chemistry of each life class is, which classes are speculative, and how sure it is that each class can exist. `data` is `speculative.json`,
    or `None` before the build has run. A basis that cites a source (`[R2336]`) is shown with the citation as a link (`cite`).
    """
    page = {
        "title": "Speculative chemistry",
        "intro": [
            "This page says how sure the air chemistry of each life class is. The planet pages show the result: which gases a class "
            "adds to the air and takes out of it. Here is how well each part is known.",
            "A class is speculative when its air numbers are hand-tuned against many worlds (level L3) or fitted with no independent source (level L4).",
        ],
        "sections": [],
    }
    if data is None:
        page["sections"].append({
            "id": "not-built", "heading": "The pages are not built", "table": None,
            "paragraphs": ["The numbers on this page are built at setup by scripts/build_reference_pages.py. Run it, or run update.sh, and reload."]})
        return page
    classes = data["classes"]
    page["sections"].append({
        "id": "classes", "heading": "All the life classes", "paragraphs": [
            "Each class has two tiers. The first says how sure we are that the class can exist at all. The second, with its level, says how sure we are "
            "of how the class changes the air."],
        "table": {"head": ["Class", "Can exist", "Air chemistry tier", "Adjustment level", "Speculative"], "rows": [
            [f"Class {row['id']}", f"{row['exists_tier']}: {_tier_name(data, row['exists_tier'])}",
             f"{row['air_tier']}: {_tier_name(data, row['air_tier'])}", f"{row['air_tweak']}: {_level_name(data, row['air_tweak'])}",
             "Yes" if row["speculative"] else "No"] for row in classes]}})
    page["sections"].append({
        "id": "several", "heading": "Several classes at once", "paragraphs": [
            "A world can hold more than one class of life. Its air is then the sum of what each class adds and takes out. "
            "The combined air is only as sure as its least sure class."],
        "table": None})
    seen = set()
    for row in classes:
        paragraphs = [
            f"Air chemistry: tier {row['air_tier']} ({_tier_name(data, row['air_tier']).lower()}), level {row['air_tweak']} "
            f"({_level_name(data, row['air_tweak']).lower()}). " + ("This class is speculative." if row["speculative"] else "This class is not speculative."),
            f"Can exist: tier {row['exists_tier']} ({_tier_name(data, row['exists_tier']).lower()}). How far it can go: {row['reaches'].lower()}.",
        ]
        if row["number"] not in seen:
            seen.add(row["number"])
            paragraphs.append(f"The air numbers of the whole Class {row['number']} family ({row['air_name'].lower()}) rest on: {row['air_basis']}.")
        page["sections"].append({"id": f"class-{row['id'].lower()}", "heading": f"Class {row['id']}: {row['name']}", "paragraphs": paragraphs, "table": None})
    for section in page["sections"]:
        section["paragraphs"] = [cite(text, folder) if "[R" in text else text for text in section["paragraphs"]]
    return page
