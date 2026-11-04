#!/usr/bin/env python3
"""
Builds the data of the References pages (DOC.88, References R2) from the registries they describe.

    python scripts/build_reference_pages.py            # write docs/references/*.json
    python scripts/build_reference_pages.py --check    # build in memory and fail if the committed files differ

The References section (Help > References, `/references`, `docs/design/interface/references.md`) is the one place that speaks out of
universe: where each number came from, what is a stand-in, what was fitted to a result and how sure each part is. Its pages are static and
need no database. Every number on them is read from the registries here, never typed by hand, so it cannot drift; `--check` is the
test's way to prove it (`src/tests/test_references_pages.py`), and `update.sh` runs the build at setup like the other generated pages.

Inputs:

* `docs/reference/references-data.json`: the papers, books and datasets (`scripts/build_references.py`), the datasets fetched at setup and
  the archives;
* `docs/design/research/gap-filling-plan.md`: the stand-in registry (GEN.725), read with `scripts/build_stand_in_values.py`'s own parser.

Output, in `docs/references/`:

* `overview.json`: the counts the start page shows (papers, links, datasets, archives, stand-ins by tier and by tweak), the legends of the
  tiers and tweak levels, and the headline verdict, computed from the counts.

* `stand-ins.json` (DOC.90): every row of the stand-in registry, section by section, with the legend text of each field, the item that will
  wire it (`wiring-register.md`) and what a hand-tuned row needs (`burn-down-register.md`). One page of `/references/stand-ins/<slug>`
  is built from each row.
* `fitted.json` (DOC.90): the tuning targets (`src/planetgen/physics/data/tuning_targets.json`) that the rows tuned to a published
  statistic point at, for `/references/fitted`.

Later items add their files here (the papers, data, method and speculative pages: DOC.89, DOC.91 and DOC.92); each build function is added
to `BUILDERS`, and `--check` covers it at once.
"""

import argparse
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))

import build_references as refs  # noqa: E402
import build_stand_in_values as stand_ins  # noqa: E402

DATA = os.path.join(ROOT, "docs", "reference", "references-data.json")
PLAN = os.path.join(ROOT, "docs", "design", "research", "gap-filling-plan.md")
WIRING = os.path.join(ROOT, "docs", "design", "research", "wiring-register.md")
BURN_DOWN = os.path.join(ROOT, "docs", "design", "research", "burn-down-register.md")
TARGETS = os.path.join(ROOT, "src", "planetgen", "physics", "data", "tuning_targets.json")
OUT = os.path.join(ROOT, "docs", "references")

READ_FROM_A_SOURCE = ("T1", "T2")
"""tuple: The tiers whose value is read from a source (verified, or from a source we trust)."""

DERIVED_FROM_PHYSICS = ("T5", "T6", "T7", "T8")
"""tuple: The tiers whose value is an educated guess from physics."""

CHOSEN_FOR_THE_OUTPUT = ("T9",)
"""tuple: The tiers whose value is chosen so the output comes out right."""

FITTED_TWEAKS = ("L2", "L3", "L4")
"""tuple: The tweak levels that were bent toward a result: a published statistic, a sweep, or the outcome alone."""


def paper_counts(data):
    """
    How many papers and books the reference list holds and how each can be looked up: `{"entries", "with_doi", "without_doi", "with_arxiv",
    "datasets", "books"}`. A dataset or a book is counted among the entries too.
    """
    entries = data["entries"].values()
    kinds = [entry.get("type") for entry in entries]
    return {
        "entries": len(data["entries"]),
        "with_doi": sum(1 for entry in data["entries"].values() if entry.get("doi")),
        "without_doi": sum(1 for entry in data["entries"].values() if not entry.get("doi")),
        "with_arxiv": sum(1 for entry in data["entries"].values() if entry.get("arxiv")),
        "datasets": kinds.count("dataset"),
        "books": sum(kinds.count(kind) for kind in ("book", "monograph", "edited-book")),
    }


def stand_in_counts(plan_text):
    """
    The stand-in registry in numbers: `{"rows", "by_tier", "by_tweak"}` (a code to a count, every code present), counted from the registry
    tables of the gap-filling plan by the parser of `build_stand_in_values.py`.
    """
    tables = stand_ins.registry_tables(plan_text)
    found = stand_ins.problems(tables)
    if found:
        raise ValueError("the stand-in registry has rows without a tier or a tweak: " + "; ".join(found[:3]))
    total, tier_count, tweak_count, _cross = stand_ins._counts(tables)
    return {"rows": total, "by_tier": tier_count, "by_tweak": tweak_count}


def share(part, whole):
    """`part` as a whole percent of `whole` (an integer, rounded half up; 0 for an empty whole)."""
    return (200 * part + whole) // (2 * whole) if whole else 0


def verdict(counts):
    """
    The headline of the start page, computed from the stand-in counts and never rated by hand: how many values are read from a source we
    trust, how many are an educated guess from physics, how many are chosen for the output, and how many were bent toward a result.
    """
    total, tiers, tweaks = counts["rows"], counts["by_tier"], counts["by_tweak"]
    read = sum(tiers[code] for code in READ_FROM_A_SOURCE)
    physics = sum(tiers[code] for code in DERIVED_FROM_PHYSICS)
    chosen = sum(tiers[code] for code in CHOSEN_FOR_THE_OUTPUT)
    fitted = sum(tweaks[code] for code in FITTED_TWEAKS)
    return {
        "read_from_a_source": {"rows": read, "percent": share(read, total)},
        "educated_guess_from_physics": {"rows": physics, "percent": share(physics, total)},
        "chosen_for_the_output": {"rows": chosen, "percent": share(chosen, total)},
        "unverified": {"rows": total - read - physics - chosen, "percent": share(total - read - physics - chosen, total)},
        "fitted_to_a_result": {"rows": fitted, "percent": share(fitted, total)},
    }


def build_overview():
    """The data of the start page, `/references`."""
    with open(DATA, encoding="utf-8") as handle:
        data = json.load(handle)
    with open(PLAN, encoding="utf-8") as handle:
        counts = stand_in_counts(handle.read())
    return {
        "papers": paper_counts(data),
        "data": {"downloaded_at_setup": len(data["setup_datasets"]), "archives": len(data["archives"])},
        "stand_ins": counts,
        "verdict": verdict(counts),
        "tiers": [list(row) for row in stand_ins.TIERS],
        "tweaks": [list(row) for row in stand_ins.TWEAKS],
    }


def source_label(e, suffix):
    """The words a citation link shows for a paper or a book: family name or names, a year and a letter where two share one."""
    names = [family.title() if family.isupper() and len(family) > 3 else family for family, _ in e["authors"]]
    if not names:
        who = (e.get("title") or "").strip()
    elif len(names) == 1:
        who = names[0]
    elif len(names) == 2:
        who = names[0] + " and " + names[1]
    else:
        who = names[0] + " et al."
    return ("%s %s%s" % (who, e.get("year") or "", suffix)).strip()


def build_citations():
    """
    The citation links of the References pages (DOC.97): `{"sources": {id: {"label", "title", "kind"}}, "aliases": {merged id: id}}`.

    A text on a References page cites `[R0123]`; `planetgen.web.lib.citations.link` turns it into the label and a link to the entry. The label
    is built the way the reference list prints the author and year, with the same `a`, `b` letter, so the link and the entry cannot differ.
    """
    with open(DATA, encoding="utf-8") as handle:
        cache = json.load(handle)
    keyed, suffix = refs.sorted_entries({k: e for k, e in cache["entries"].items() if e.get("authors") or e.get("title")})
    sources = {}
    for key, e in keyed:
        sources[refs.id_of(cache, key)] = {
            "kind": "book" if e.get("type") in ("monograph", "book", "edited-book") else "paper",
            "label": source_label(e, suffix.get(key, "")),
            "title": e.get("title") or "",
        }
    for a in cache.get("archives", []):
        sources[cache["ids"]["archive:" + a["url"].lower()]] = {"kind": "data", "label": a["name"], "title": a["name"]}
    for w in cache.get("web_sources", []):
        sources[cache["ids"]["web:" + w["url"].lower()]] = {"kind": "web", "label": w["name"], "title": w["name"]}
    for d in cache.get("setup_datasets", []):
        sources[cache["ids"]["setup:" + d["url"].lower()]] = {"kind": "data", "label": d["title"], "title": d["title"]}
    return {"aliases": dict(sorted(cache.get("id_aliases", {}).items())), "sources": sources}


KINDS = {"D": "derived from a source", "B": "bracket between sources", "T": "theory", "P": "placeholder"}
"""dict: The kind letters of a registry row (`gap-filling-plan.md`, section 5): a value derived from a source, a bracket between sources, a
theory fill, or a placeholder with no source."""

FIELD_LABELS = {
    "Fill (ranked)": "How it is filled", "Draw": "How it is drawn", "Wrong if": "It is wrong if", "Result lands as": "The result lands as",
    "constant": "What it is", "knobs": "The range it moves in", "validity": "Where it holds", "falsifier": "What would show it wrong",
    "replaced_by": "What replaces it", "evidence": "Evidence", "stand-in": "What it is", "sources": "Sources", "Range": "The range",
    "Label": "How it is labeled", "Source": "Source",
}
"""dict: The words a field of a registry table goes by on a stand-in's page."""

HEADLINE = ("constant", "Fill (ranked)", "stand-in", "Lever")
"""tuple: The columns that say what a row is, the first one a table has; the index shows its first words."""

SKIP_FIELDS = ("Id", "fill_id", "Lever", "Kind", "kind", "Tier", "tier", "Tweak", "tweak")
"""tuple: The columns a stand-in's page shows in its head and not in its body."""


def clean(text):
    """A registry cell as the page shows it: no markdown ticks, one line."""
    return " ".join(text.replace("`", "").split())


def build_stand_ins():
    """
    The stand-in pages (DOC.90): `{"sections": [{"key", "number", "title", "slug", "rows": [row]}], "kinds"}`; a row is `{"key", "slug", "id", "kind", "tier",
    "tier_code", "tweak", "what", "fields": [[label, text]], "unread", "burn_down"}`. The key is the registry's own (`6.7:AM1`, the section and the id)
    and the slug is the page's address (`6-7-am1`). Every text is passed through the citation links of the References pages, so an unknown ID
    or an author and year written out stops the build and names the row.
    """
    from planetgen.web.lib import citations, reference_urls
    cites = build_citations()
    with open(PLAN, encoding="utf-8") as handle:
        plan = handle.read()
    tables = stand_ins.registry_tables(plan)
    found = stand_ins.problems(tables)
    if found:
        raise ValueError("the stand-in registry is not complete: " + "; ".join(found[:3]))
    with open(WIRING, encoding="utf-8") as handle:
        unread = stand_ins.unread_rows(handle.read())
    with open(BURN_DOWN, encoding="utf-8") as handle:
        burn = stand_ins.burn_down_rows(handle.read())
    sections, slugs = [], {}
    for key, number, title, _columns in stand_ins.sections():
        header, rows = tables[key]
        lowered = [name.lower() for name in header]
        out = []
        for row in rows:
            cells = dict(zip(header, row))
            row_key = f"{number}:{row[0].strip('`')}"
            slug = reference_urls.anchor(row_key)
            if slug in slugs:
                raise ValueError(f"{row_key} and {slugs[slug]} are the same page address {slug}")
            slugs[slug] = row_key
            tier_cell = clean(cells[header[lowered.index("tier")]])
            head = next((clean(cells[name]) for name in HEADLINE if name in cells and cells[name]), clean(row[0]))
            if "Fill (ranked)" in cells and head == clean(cells["Fill (ranked)"]):
                head = clean(stand_ins._first_option(cells["Fill (ranked)"]))
            by_label = "Kind" not in cells and "kind" not in cells
            kind = clean(cells.get("Kind") or cells.get("kind") or cells.get("Label") or "")
            kind = {name.split()[0]: letter for letter, name in KINDS.items()}.get(kind, kind)
            fields = [[FIELD_LABELS.get(name, name), clean(cells[name])] for name in header
                      if name not in SKIP_FIELDS and not (by_label and name == "Label") and cells[name].strip() not in ("", "-")]
            if head == clean(row[0]) and fields:
                head = f"{head}: {fields[0][1]}"
            for text in [head] + [text for _label, text in fields]:
                try:
                    citations.link(text, cites, strict=False)
                except citations.CitationError as error:
                    raise ValueError(f"stand-in {row_key}: {error}") from error
            disposition = burn.get(row_key)
            out.append({
                "key": row_key, "slug": slug, "id": clean(row[0]),
                "kind": kind,
                "tier": tier_cell, "tier_code": stand_ins.TIER_CODE.match(tier_cell).group(0),
                "tweak": clean(cells[header[lowered.index("tweak")]]), "what": head, "fields": fields,
                "unread": unread.get(row_key), "burn_down": list(disposition) if disposition else None,
            })
        sections.append({"key": key, "number": number, "title": re.sub(r"^\d+\.\d+a? ", "", title), "slug": reference_urls.anchor(title), "rows": out})
    return {"kinds": KINDS, "sections": sections}


def build_fitted():
    """The tuning targets for `/references/fitted` (DOC.90), as `{"targets": [target]}`; each target keeps what the table states about it."""
    from planetgen.web.lib import citations
    cites = build_citations()
    with open(TARGETS, encoding="utf-8") as handle:
        table = json.load(handle)
    targets = []
    for target in table["targets"]:
        entry = {name: target.get(name) for name in ("id", "row", "statistic", "kind", "value", "sigma", "n", "label", "source", "falsifier", "role")}
        for name in ("statistic", "source", "falsifier"):
            try:
                citations.link(entry[name] or "", cites, strict=False)
            except citations.CitationError as error:
                raise ValueError(f"tuning target {entry['id']}: {error}") from error
        targets.append(entry)
    return {"targets": targets}


BUILDERS = {"overview.json": build_overview, "citations.json": build_citations, "stand-ins.json": build_stand_ins, "fitted.json": build_fitted}
"""dict: `{file under docs/references/: function that returns its data}`."""


def render(data):
    """A data file's text: sorted keys, two-space indent, one trailing newline, so a build is the same bytes every time."""
    return json.dumps(data, indent=2, sort_keys=True, ensure_ascii=False) + "\n"


def build():
    """`{file name: text}` for every file of `BUILDERS`."""
    return {name: render(builder()) for name, builder in BUILDERS.items()}


def stale(built, folder=None):
    """The names of the files in `built` that are missing from `folder` (`docs/references/`) or differ from it, and the files there that nothing builds."""
    folder = folder or OUT
    names = [name for name, text in built.items() if not _same(os.path.join(folder, name), text)]
    if os.path.isdir(folder):
        names += [f"{name} (nothing builds it)" for name in sorted(os.listdir(folder)) if name not in built]
    return names


def _same(path, text):
    if not os.path.isfile(path):
        return False
    with open(path, encoding="utf-8") as handle:
        return handle.read() == text


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--check", action="store_true", help="build in memory and fail if the committed files differ")
    args = parser.parse_args(argv)
    built = build()
    if args.check:
        differ = stale(built)
        if differ:
            print("docs/references is out of date: " + ", ".join(differ) + "\nrun: python scripts/build_reference_pages.py", file=sys.stderr)
            return 1
        print(f"docs/references is up to date ({len(built)} file{'s' if len(built) != 1 else ''}).")
        return 0
    os.makedirs(OUT, exist_ok=True)
    for name, text in built.items():
        with open(os.path.join(OUT, name), "w", encoding="utf-8") as handle:
            handle.write(text)
    print(f"wrote {len(built)} file{'s' if len(built) != 1 else ''} to docs/references/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
