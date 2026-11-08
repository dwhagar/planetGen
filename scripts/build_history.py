#!/usr/bin/env python3
"""
Builds the project history, `docs/history/`: what was chosen, what was dropped and how the work moved, read from the
repository itself (git log, CHANGELOG.md and the Markdown documents).

Boss (2026-11-08, 6:13 PM EST): "read through every single commit, todo, choice, issue ... alive or dead, and add past
choices and paths and trends to our design documents."

This is a SNAPSHOT, not a global document: it reads git history, which a shallow clone does not have, so no test or gate
compares it with its sources. Run it (in a full clone) when the record should be brought up to date, and commit the result.

What it writes
==============
* `README.md`: the index, the method and the as-of commit.
* `timeline.md`: commits, releases, items filed, items retired and pull requests per week, by category of work.
* `releases-YYYY-MM.md`: every release headline of CHANGELOG.md, day by day (the headline is "the single best feature
  that release brought").
* `items-<CATEGORY>.md`: the life of every item ID that appears in a commit title: the first and last date, whether it
  is open in `docs/TODO.md` now, the pull requests, and the commit titles in order.
* `dropped-and-closed.md`: every commit that retired, dropped, closed, rejected, reverted or replaced something, with the
  items named, so a path that was tried and left is findable.
* `decisions-YYYY-MM.md`: every sentence in the documents that names Boss and a date (the choices as they were made),
  with the file and line, oldest first.
* `by-note.md`: for every design note, the items it names and the span of dates those items were worked.

Usage
=====
    python scripts/build_history.py
"""

import collections
import datetime
import glob
import os
import re
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "docs", "history")
ID = re.compile(r"\b(UX|MAP|NAV|GEN|PERF|DB|API|ADM|SEC|TEST|USR|OPS|DOC|VIEW|POP)\.(\d+)\b")
PR = re.compile(r"\(#(\d+)\)\s*$")
DROP = re.compile(r"\b(retire[ds]?|drop(?:ped|s)?|close[ds]?|reject(?:ed|s)?|revert(?:ed|s)?|replace[ds]?|supersede[ds]?|abandon(?:ed)?|withdraw[ns]?|"
                  r"delete[ds]?|remove[ds]?|not built|as a finding|won'?t|no longer)\b", re.I)
CATEGORIES = {"UX": "Web pages", "MAP": "Maps", "NAV": "Navigation", "GEN": "Generation and physics", "PERF": "Speed",
              "DB": "Database", "API": "API", "ADM": "Admin", "SEC": "Security", "TEST": "Tests", "USR": "Accounts",
              "OPS": "Operations", "DOC": "Documentation", "VIEW": "Views", "POP": "Population"}
DECISION = re.compile(r"(Boss[^.\n]{0,100}20\d\d-\d\d-\d\d|20\d\d-\d\d-\d\d[^.\n]{0,100}Boss)")
DATE = re.compile(r"20\d\d-\d\d-\d\d")


def citations_by_id(name, text):
    """`text` with every source written by name (author and year, DOI, arXiv number) replaced by its `[R0123]` ID (DOC.99:
    a citation is written only in the registry); a name the list does not hold is dropped by `rewrite_citations`."""
    sys.path.insert(0, os.path.join(ROOT, "scripts"))
    import rewrite_citations as rc
    reg = rc.Registry(all_candidates=True)
    text = rc.rewrite_file(name, text, reg, rc.Counter(), strings=True)
    for _key, mention in reg.unresolved:  # a source the list does not hold stays unnamed
        text = text.replace(mention, "a source")
    return text


def git(*args):
    return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True, check=True).stdout


def commits():
    """`[(sha, date, subject, body)]`, oldest first."""
    raw = git("log", "--reverse", "--format=%H%x1f%ad%x1f%s%x1f%b%x1e", "--date=short")
    found = []
    for record in raw.split("\x1e"):
        parts = record.strip("\n").split("\x1f")
        if len(parts) == 4:
            found.append(tuple(parts))
    return found


def write(name, lines):
    os.makedirs(OUT, exist_ok=True)
    text = citations_by_id("docs/history/" + name, "\n".join(lines) + "\n")
    with open(os.path.join(OUT, name), "w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)


def esc(text):
    return text.replace("|", "\\|").replace("<", "&lt;")


def week_of(date):
    d = datetime.date.fromisoformat(date)
    start = d - datetime.timedelta(days=d.weekday())
    return start.isoformat()


def open_items():
    """Item IDs written as an open item in docs/TODO.md."""
    with open(os.path.join(ROOT, "docs", "TODO.md"), encoding="utf-8") as handle:
        return set(m.group(1) for m in re.finditer(r"^\s*- \[[ x]\] \*\*([A-Z]+\.\d+) ", handle.read(), re.M))


def build():
    log = commits()
    live = open_items()
    head = git("rev-parse", "--short", "HEAD").strip()
    today = log[-1][1]
    items = collections.defaultdict(list)  # id -> [(date, subject, pr)]
    for _sha, date, subject, _body in log:
        pr = PR.search(subject)
        for m in sorted(set(f"{a}.{b}" for a, b in ID.findall(subject))):
            items[m].append((date, subject, pr.group(1) if pr else ""))

    # --- timeline
    weekly = collections.OrderedDict()
    for _sha, date, subject, _body in log:
        w = weekly.setdefault(week_of(date), collections.Counter())
        w["commits"] += 1
        if subject.startswith("Release "):
            w["releases"] += 1
        if PR.search(subject):
            w["prs"] += 1
        if re.match(r"(File|Filed)\b", subject):
            w["filed"] += 1
        if re.match(r"Retire\b", subject):
            w["retired"] += 1
        if re.match(r"(Train \d+|.*\btrain \d+)", subject, re.I):
            w["train"] += 1
    lines = ["# Timeline", "",
             f"Commits per week, from the git history as of `{head}` ({today}). `Pull requests` counts squash-merged titles ending in `(#N)`; "
             "`Releases` counts `Release x.y.z` commits (one per changelog entry); `Filed` and `Retired` count commits whose title starts with "
             "`File` or `Retire` (an item opened or closed in docs/TODO.md); `Train` counts commits that name a merge train.", "",
             "| Week starting | Commits | Pull requests | Releases | Filed | Retired | Train |", "|---|---|---|---|---|---|---|"]
    for week, c in weekly.items():
        lines.append(f"| {week} | {c['commits']} | {c['prs']} | {c['releases']} | {c['filed']} | {c['retired']} | {c['train']} |")
    lines += ["", "## Work by category", "",
              "Item IDs named in commit titles, by category: how many different items, how many commits, and the first and last date.", "",
              "| Category | What | Items | Commits | First | Last |", "|---|---|---|---|---|---|"]
    by_cat = collections.defaultdict(list)
    for item, events in items.items():
        by_cat[item.split(".")[0]].append((item, events))
    for cat in CATEGORIES:
        rows = by_cat.get(cat, [])
        if rows:
            dates = [d for _i, ev in rows for d, _s, _p in ev]
            lines.append(f"| {cat} | {CATEGORIES[cat]} | {len(rows)} | {sum(len(ev) for _i, ev in rows)} | {min(dates)} | {max(dates)} |")
    write("timeline.md", lines)

    # --- releases
    with open(os.path.join(ROOT, "CHANGELOG.md"), encoding="utf-8") as handle:
        changelog = handle.read()
    entries = re.findall(r"^## \[([^\]]+)\] - (\d{4}-\d\d-\d\d)\n+(?:\*\*Headline:\*\* (.*)\n)?", changelog, re.M)
    by_month = collections.defaultdict(list)
    for version, date, headline in entries:
        by_month[date[:7]].append((date, version, headline))
    for month, rows in sorted(by_month.items()):
        out = [f"# Release headlines, {month}", "",
               f"{len(rows)} releases. Each headline is the single best feature that release brought (CHANGELOG.md, newest first). "
               "Releases without a headline predate the rule (2026-10-10).", ""]
        day = None
        for date, version, headline in sorted(rows, key=lambda r: (r[0], r[1]), reverse=True):
            if date != day:
                day = date
                out += ["", f"## {date}", ""]
            out.append(f"- **{version}**: {esc(headline) if headline else '(no headline)'}")
        write(f"releases-{month}.md", out)

    # --- items
    for cat in CATEGORIES:
        rows = sorted(by_cat.get(cat, []), key=lambda r: int(r[0].split(".")[1]))
        if not rows:
            continue
        out = [f"# Items: {cat}, {CATEGORIES[cat]}", "",
               f"The life of every {cat} item that a commit title names, from the git history as of `{head}`: first and last date, whether it is "
               "open in [docs/TODO.md](../TODO.md) now, and its commit titles in order (the pull request numbers are the GitHub pull requests). "
               "An item not open in TODO.md has shipped, been closed as a finding or been dropped; "
               "[dropped-and-closed.md](dropped-and-closed.md) lists the closings.", ""]
        for item, events in rows:
            state = "open" if item in live else "not open"
            out += [f"## {item}", "", f"{events[0][0]} to {events[-1][0]}; {len(events)} commits; {state}.", ""]
            shown = events if len(events) <= 12 else events[:6] + [None] + events[-5:]
            for event in shown:
                out.append("- ..." if event is None else f"- {event[0]}: {esc(event[1])}")
            out.append("")
        write(f"items-{cat}.md", out)

    # --- dropped and closed
    out = ["# Dropped, retired, closed and replaced", "",
           f"Every commit title (as of `{head}`) that has one of these words in it: retire, drop, close, reject, revert, replace, supersede, abandon, "
           "withdraw, delete, remove, \"as a finding\", \"not built\", \"no longer\". A `Retire GEN.n` commit usually means the item shipped; "
           "read the title for \"finding\", \"dropped\" or \"closed\". Newest first.", ""]
    hits = [(d, s) for _sha, d, s, _b in reversed(log) if DROP.search(s) and not s.startswith("Release ")]
    month = None
    for date, subject in hits:
        if date[:7] != month:
            month = date[:7]
            out += ["", f"## {month}", ""]
        out.append(f"- {date}: {esc(subject)}")
    write("dropped-and-closed.md", out)

    # --- decisions
    found = collections.defaultdict(list)
    for path in sorted(glob.glob(os.path.join(ROOT, "docs", "**", "*.md"), recursive=True)):
        rel = os.path.relpath(path, ROOT).replace(os.sep, "/")
        if rel.startswith(("docs/history/", "docs/generated/", "docs/site/")) or rel in ("docs/TODO.md", "docs/plan/todo-number-map.md"):
            continue
        with open(path, encoding="utf-8") as handle:
            for number, line in enumerate(handle, 1):
                for sentence in re.split(r"(?<=[.!?])\s+(?=[A-Z(])", re.sub(r"\[([^\]]*)\]\([^)]*\)", r"\1", line.strip())):
                    if DECISION.search(sentence) and len(sentence) > 30:
                        dates = DATE.findall(sentence)
                        if dates:
                            found[dates[0][:7]].append((dates[0], rel, number, sentence[:420]))
    total = sum(len(v) for v in found.values())
    for month, rows in sorted(found.items()):
        out = [f"# Decisions, {month}", "",
               f"{len(rows)} sentences in the documents (outside TODO.md and the number map) that name Boss and a date in {month}; the choices as "
               "they were made, oldest first, with the file and line they are written on.", ""]
        for date, rel, number, sentence in sorted(rows):
            out.append(f"- {date} ([{rel}:{number}](../{rel[len('docs/'):]})): {esc(sentence)}")
        write(f"decisions-{month}.md", out)

    # --- by note
    out = ["# Design notes and the work behind them", "",
           f"For each design note, the items it names, with the dates those items were worked (git history as of `{head}`) and where they stand. "
           "Use it to jump from a note to the history of its items in [items-<category>.md](README.md).", ""]
    for path in sorted(glob.glob(os.path.join(ROOT, "docs", "design", "**", "*.md"), recursive=True)):
        rel = os.path.relpath(path, ROOT).replace(os.sep, "/")
        with open(path, encoding="utf-8") as handle:
            ids = sorted(set(f"{a}.{b}" for a, b in ID.findall(handle.read())), key=lambda i: (i.split(".")[0], int(i.split(".")[1])))
        known = [i for i in ids if i in items]
        if not known:
            continue
        dates = [d for i in known for d, _s, _p in items[i]]
        opened = [i for i in known if i in live]
        out += [f"## [{rel[len('docs/design/'):]}](../{rel[len('docs/'):]})", "",
                f"{len(known)} items with commits, {min(dates)} to {max(dates)}; {len(opened)} still open"
                + (f" ({', '.join(opened[:12])}{', ...' if len(opened) > 12 else ''})" if opened else "") + ".", ""]
    write("by-note.md", out)

    # --- readme
    pages = ["timeline.md", "dropped-and-closed.md", "by-note.md"]
    out = ["# Project history", "",
           f"What was chosen, what was dropped and how the work moved, read from the repository as of commit `{head}` ({today}): "
           f"{len(log)} commits, {len(entries)} releases, {len(items)} items named in commit titles, {total} dated decisions in the documents. "
           "Built by `scripts/build_history.py`; this folder is a snapshot (git history is not in a shallow clone, so no test compares it "
           "with its sources). Rebuild it in a full clone when it should catch up, and commit the result.", "",
           "## Where to look", "",
           "- [timeline.md](timeline.md): commits, releases, items filed and retired per week; work by category.",
           "- [dropped-and-closed.md](dropped-and-closed.md): paths that were retired, dropped, closed as findings, reverted or replaced.",
           "- [by-note.md](by-note.md): each design note with the dates and state of the items behind it.",
           "- Release headlines by month (the best feature of each release): " + ", ".join(f"[{m}](releases-{m}.md)" for m in sorted(by_month)) + ".",
           "- Decisions by month (Boss's choices as the documents record them): " + ", ".join(f"[{m}](decisions-{m}.md)" for m in sorted(found)) + ".",
           "- Item histories: " + ", ".join(f"[{c}](items-{c}.md)" for c in CATEGORIES if by_cat.get(c)) + ".", "",
           "## How the design evolved, in four eras", "",
           "Read from the timeline, the release headlines and the decision log; each era's decisions are in the month pages above.", "",
           "- **Before 2026-09**: a first generator (changelog 0.1 to 7.x, June 2024), then a pause.",
           "- **2026-09 to 2026-10-09 (Phase 0 and 1)**: the package layout, MySQL, Alembic, per-unit seeds, the map engine and Redis queues; "
           "release numbers 7.x. Phase 0 closed on 2026-10-09 and the version became 8.0.",
           "- **2026-10-10 to 2026-10-31**: the version map (8.1 to 10.0, 2026-10-16), the lane model (cores, bugs, aux), the TODO thread as sole owner of items "
           "and global documents, the stand-in and source rules, condition-based planet classes and the nine life classes.",
           "- **2026-11 onward (8.1 build-out)**: merge trains (train 1 on 2026-10-29 to train 19 on 2026-11-08), the gate, lanes renumbered "
           "on 2026-11-07, the everything-affects-everything principle (2026-11-08) and real materials.", ""]
    write("README.md", out)
    return len(log), len(items), total


def main(argv=None):
    commits_count, item_count, decision_count = build()
    print(f"docs/history: {commits_count} commits, {item_count} items, {decision_count} decisions")
    return 0


if __name__ == "__main__":
    sys.exit(main())
