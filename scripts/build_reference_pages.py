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

Later items add their files here (the papers, data, stand-in, method and speculative pages: DOC.89 to DOC.92); each build function is added
to `BUILDERS`, and `--check` covers it at once.
"""

import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

import build_references as refs  # noqa: E402
import build_stand_in_values as stand_ins  # noqa: E402

DATA = os.path.join(ROOT, "docs", "reference", "references-data.json")
PLAN = os.path.join(ROOT, "docs", "design", "research", "gap-filling-plan.md")
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
    for d in cache.get("setup_datasets", []):
        sources[cache["ids"]["setup:" + d["url"].lower()]] = {"kind": "data", "label": d["title"], "title": d["title"]}
    return {"aliases": dict(sorted(cache.get("id_aliases", {}).items())), "sources": sources}


BUILDERS = {"overview.json": build_overview, "citations.json": build_citations}
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
