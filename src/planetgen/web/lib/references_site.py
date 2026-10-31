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
