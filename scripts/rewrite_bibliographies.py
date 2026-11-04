#!/usr/bin/env python3
"""
Turns the APA-form bibliography lists of the imported research notes into citation IDs (DOC.98 follow-up).

    python scripts/rewrite_bibliographies.py report           # what each entry would become
    python scripts/rewrite_bibliographies.py apply            # write the files and add the sources found to the list

The imported notes (`docs/design/research/imported/`, `docs/design/interface/3d-rendering-camera-research.md`) end in lists such as
`* Atri, D., & Melott, A. L. (2011). Title. _Journal_, 35(1), 31-37.` which `scripts/rewrite_citations.py` cannot read as a mention. An entry is
tied to the reference list (`docs/reference/references-data.json`) when the first author, the year and most of the list's title words are in
the entry; otherwise Crossref is searched with the entry, and a record whose first author, year (within one) and title words agree is added to
the list (`src="crossref-search"`), as `build_references.py add` does. A tied entry becomes `* [R0123]`. An entry that nothing verifies is left
as it is and listed, because an unverifiable paper is treated as not real (the text is removed by hand, or the source added with `build_references.py add`).
After `apply`, run `python scripts/build_references.py render` and `python scripts/build_reference_pages.py`.
"""

import argparse
import glob
import json
import os
import re
import sys
import time
import unicodedata

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

import build_references as br  # noqa: E402

FILES = ("docs/design/research/imported/*.md", "docs/design/interface/3d-rendering-camera-research.md")
APA = re.compile(r"^(?P<lead>\s*(?:[*\-]|\d+\.)?\s*)(?P<first>[A-Z][A-Za-zÀ-ſ'’\- ]*?),\s+(?:[A-Z]\.[ -]?)+.*?\((?P<year>(?:19|20)\d\d)[a-z]?\)\.\s*(?P<rest>.*)$")
IDS = re.compile(r"\[(R\d{4,})\]")
COVER = 0.7
"""float: The share of a record's title words that must be in the entry for the two to be one paper."""
STOP = {"the", "and", "for", "with", "from", "that", "this", "are", "was", "into", "over", "under", "its", "their", "of", "in", "on", "an", "a", "to"}


def title_words(text):
    """The words (three letters or more, no stop words) of a title, lower case and without accents."""
    plain = unicodedata.normalize("NFKD", text.replace("\u00df", "ss")).encode("ascii", "ignore").decode().lower()
    return {w for w in re.findall(r"[a-z0-9]{3,}", plain) if w not in STOP}


def covers(entry_text, title):
    """True when most of the record's title words are in the entry, or most of the entry's own title (its words up to the first journal or
    publisher mark) are in the record's."""
    words = title_words(title)
    if len(words) >= 3 and len(words & title_words(entry_text)) / len(words) >= COVER:
        return True
    own = title_words(re.split(r"[_*]|\. [A-Z]", entry_text.lstrip("_* "), maxsplit=1)[0])
    return len(own) >= 3 and len(own & words) / len(own) >= COVER


def surname(first):
    return br.fold(first.split()[-1]) if first.split() else ""


def local_match(cache, first, year, rest):
    """The key of the listed source this entry names, or `None`."""
    found = []
    for key, e in cache["entries"].items():
        if e.get("authors") and e.get("year") and surname(e["authors"][0][0]) == surname(first) and abs(int(e["year"]) - int(year)) <= 1:
            if covers(rest, e.get("title") or ""):
                found.append(key)
    return found[0] if len(found) == 1 else None


def crossref_match(first, year, rest):
    """A Crossref record for this entry (first author, year within one, title words), or `None`."""
    query = "%s %s %s" % (first, year, re.sub(r"[*_]", "", rest)[:240])
    url = "https://api.crossref.org/works?rows=5&mailto=%s&query.bibliographic=%s" % (br.MAILTO, br.urllib.parse.quote(query))
    raw = br.http_get(url)
    if not raw:
        return None
    try:
        items = json.loads(raw)["message"]["items"]
    except (ValueError, KeyError):
        return None
    for item in items:
        title = (item.get("title") or [""])[0]
        authors = item.get("author") or []
        issued = ((item.get("issued") or {}).get("date-parts") or [[None]])[0][0]
        if not (title and authors and issued and abs(int(issued) - int(year)) <= 1):
            continue
        if surname(authors[0].get("family", "")) == surname(first) and covers(rest, title):
            return item["DOI"].lower()
    return None


def entries(path):
    with open(os.path.join(ROOT, path), encoding="utf-8") as handle:
        lines = handle.read().split("\n")
    for i, line in enumerate(lines):
        m = APA.match(line)
        if m:
            yield i, lines, m


def cmd(args):
    cache = br.load_cache()
    br.assign_ids(cache)
    tied = unresolved = added = 0
    left = []
    for pattern in FILES:
        for path in sorted(glob.glob(os.path.join(ROOT, pattern))):
            rel = os.path.relpath(path, ROOT)
            lines, changed, drop = None, False, set()
            for i, lines, m in entries(rel):
                lead, first, year, rest = m.group("lead"), m.group("first"), m.group("year"), m.group("rest")
                have = sorted(set(IDS.findall(lines[i])))
                key = None
                if not have:
                    key = local_match(cache, first, year, rest)
                    if key is None:
                        doi = crossref_match(first, year, rest)
                        time.sleep(0.2)
                        if doi:
                            if doi not in cache["entries"]:
                                r = br.fetch_crossref(doi)
                                if r:
                                    r["src"] = "crossref-search"
                                    cache["entries"][doi] = r
                                    added += 1
                            if doi in cache["entries"]:
                                br.assign_ids(cache)
                                key = doi
                    have = [br.id_of(cache, key)] if key else []
                if have:
                    tied += 1
                    new = lead + ", ".join("[%s]" % h for h in have)
                    if new != lines[i]:
                        lines[i] = new
                        changed = True
                else:
                    unresolved += 1
                    left.append("%s:%d: %s" % (rel, i + 1, lines[i].strip()))
                    if args.drop_unverified:
                        drop.add(i)
                        changed = True
            if drop:
                kept = []
                for i, line in enumerate(lines):
                    if i in drop or (i - 1 in drop and not line.strip() and kept and not kept[-1].strip()):
                        continue
                    kept.append(line)
                lines = kept
            if changed and args.cmd == "apply":
                with open(os.path.join(ROOT, rel), "w", encoding="utf-8") as handle:
                    handle.write("\n".join(lines))
    if args.cmd == "apply":
        br.assign_ids(cache)
        br.save_cache(cache)
    print("entries tied to a source: %d (new sources found at Crossref: %d); left unverified: %d" % (tied, added, unresolved))
    for row in left:
        print("  ", row[:150])
    if args.log and args.cmd == "apply" and args.drop_unverified:
        with open(args.log, "w", encoding="utf-8") as handle:
            handle.write("\n".join(left) + "\n")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=("report", "apply"))
    ap.add_argument("--drop-unverified", action="store_true", help="remove the entries nothing verifies (an unverifiable paper is treated as not real)")
    ap.add_argument("--log", help="write the removed entries to this file")
    cmd(ap.parse_args())
