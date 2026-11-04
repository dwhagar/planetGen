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
import json
import os
import re

from . import citations
from .reference_urls import reference_path

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
