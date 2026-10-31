#!/usr/bin/env python3
"""
Builds the five browsable TODO pages from the TODO sources.

Why this exists
===============
Boss (2026-10-07 16:22Z) asked for the plan and a searchable page that
explains every TODO item in detail, rebuilt whenever `docs/TODO.md`
changes. Boss (2026-10-16, 1:07 PM EDT) turned the tiers into a version map:
one version per release from 8.1, each with one or more phases (a theme
each), and a blocker sits in the same thread, an earlier phase or an earlier
version. The pages are generated from the files the TODO tasks thread
already keeps, so one command refreshes them:

- `docs/generated/version-plan.html`: every open item in build order, version by
  version, phase by phase and thread by thread (from
  `docs/plan/releases/version-*.md`), with a progress bar per version; each thread lists
  its finished items checked off and crossed out above its open ones, and the
  finished threads that left the plan are listed per version (Boss, 2026-10-21,
  9:15 PM EDT; the finished items come from `docs/plan/progress-ledger.md`).
- `docs/generated/todo-reference.html`: every ID ever issued, open or closed,
  with its full text, status, version phase, prerequisites and what it
  unblocks, plus search and filters.
- `docs/generated/todo-tree.html`: every open item as a tree, each item under the
  prerequisite that finishes last, with version and bug marks, search and
  filters (Boss 2026-10-09 07:03Z: "a TODO as a tree structure").
- `docs/generated/execution-plan.html` (DOC.68): where each build lane is (`docs/plan/lanes.md`), the headline
  panel, the open items by version and 8.1 phase, the goal issues, and the conflicts and needs-Boss sections
  (`docs/plan/execution-notes.md`) with the older decisions (`scripts/todo_docs/execution-history.snippet.html`).
  Every part comes from a file in the repo, none from the clock, so a test can tell when it is stale.
- `docs/generated/release-roadmap.html`: the releases from 8.1 to 9.3, each with
  its headline, its phases, the open items it holds, what it needs, its share
  finished and the suggestions not yet filed (Boss, 2026-10-16).

Sources
=======
- `docs/TODO.md`: open items, their full text, nesting and prerequisites.
- `docs/plan/releases/version-*.md`: each version's goal, phases (theme) and threads,
  with "Needs" and a note per open item.
- `docs/plan/notes.md`: "Decisions for Boss".
- `docs/plan/todo-number-map.md`: the status of every ID ever issued
  (its "Reverse lookup" table).
- `docs/plan/issue-groups.md`: the GitHub issues (goal groups of items and the
  issues Boss opened) and the items each one covers; every page links them.
- `docs/plan/progress-ledger.md`: the version each finished item belonged to,
  for the progress bar on every version (Boss, 2026-10-16).
- `docs/plan/release-roadmap.md`: the release roadmap's order, headlines and
  suggestions; open counts come from TODO.md on every build.

Usage
=====
    python scripts/build_todo_docs.py          # rewrite all four pages
    python scripts/build_todo_docs.py --check  # fail if a page is stale

Each page carries a "Last updated" banner at the top (Boss, 3:50 PM EDT
2026-10-31). The stamp is the Eastern time of the last rebuild that changed
the page, kept in the page itself: a rebuild of unchanged sources keeps the
stamp and gives byte-identical pages (so `--check` stays exact), and a
rebuild that changes the page writes the new time. Templates live in
`scripts/todo_docs/`.
"""

from __future__ import annotations

import argparse
import datetime
import html
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TEMPLATES = ROOT / "scripts" / "todo_docs"
OUTPUTS = {
    "version-plan.template.html": ROOT / "docs" / "generated" / "version-plan.html",
    "todo-reference.template.html": ROOT / "docs" / "generated" / "todo-reference.html",
    "todo-tree.template.html": ROOT / "docs" / "generated" / "todo-tree.html",
    "release-roadmap.template.html": ROOT / "docs" / "generated" / "release-roadmap.html",
}
EXECUTION_TEMPLATE = TEMPLATES / "execution-plan.template.html"
EXECUTION_OUTPUT = ROOT / "docs" / "generated" / "execution-plan.html"
ISSUES_URL = "https://github.com/dwhagar/planetGen/issues/"
ID_RE = re.compile(r"^[A-Z]+\.\d+$")
HEADLINE_RE = re.compile(r"Headline:\s*(\d+\.\d+)\s*\(([^)]*)\)")
ITEM_RE = re.compile(r"^(\s*)- \[[ x]\] \*\*([A-Z]+\.\d+) (.*?)\*\*\s*$")
VERSION0 = "8.0"
VERSION0_TITLE = "Every bug and the groundwork (Phase 0, finished 2026-10-09)"
VERSION_RE = r"(?:\d+\.\d+|later)"
CATEGORIES = {
    "GEN": "Generation", "MAP": "Maps", "NAV": "Navigation", "UX": "User experience",
    "ADM": "Admin", "OPS": "Operations", "API": "API", "TEST": "Tests",
    "PERF": "Performance", "USR": "Accounts", "DB": "Database", "SEC": "Security",
    "VIEW": "View from a planet", "POP": "Population", "DOC": "Documentation",
}


def id_key(item_id: str) -> tuple:
    cat, num = item_id.split(".")
    return (cat, int(num))


def split_ids(cell: str) -> list[str]:
    return [x.strip() for x in cell.split(",") if ID_RE.match(x.strip())]


def parse_todo(text: str) -> dict:
    """Open items in TODO.md: title, body, parent, prerequisites, section."""
    items: dict[str, dict] = {}
    stack: list[tuple[int, str]] = []
    section = ""
    current = None
    body: list[str] = []

    def flush():
        if current:
            joined = "\n".join(body).strip()
            items[current]["text"] = re.sub(r"\n{2,}", "\n\n", joined)

    for line in text.splitlines():
        if line.startswith("## "):
            flush()
            current, body, stack = None, [], []
            section = line[3:].strip()
            continue
        m = ITEM_RE.match(line)
        if m:
            flush()
            indent, item_id, title = len(m.group(1)), m.group(2), m.group(3).strip()
            while stack and stack[-1][0] >= indent:
                stack.pop()
            parent = stack[-1][1] if stack else ""
            stack.append((indent, item_id))
            items[item_id] = {"title": title, "parent": parent, "section": section, "text": ""}
            current, body = item_id, []
            continue
        if current is None:
            continue
        if not line.strip():
            body.append("")
        elif re.match(r"^(#|\S)", line):
            flush()
            current, body = None, []
        else:
            body.append(line.strip())
    flush()
    for it in items.values():
        prereq = re.findall(r"Prerequisites?:\s*([^\n]*?)\.(?:\s|$)", it["text"].replace("\n", " "))
        it["prereqs"] = [x for p in prereq for x in split_ids(p)]
        it["design"] = re.findall(r"Design:\s*\[([^\]]+)\]\(([^)]+)\)", it["text"])
        found = HEADLINE_RE.search(it["text"])
        it["headline"] = (found.group(1), re.sub(r"\s+", " ", found.group(2)).strip()) if found else None
    return items


def version_key(path: Path) -> str:
    """docs/plan/releases/version-8-1-name.md -> "8.1"; version-later-name.md -> "later"."""
    m = re.match(r"version-(\d+)-(\d+)-|version-(later)-", path.name)
    return m.group(3) or f"{m.group(1)}.{m.group(2)}"


def version_sort(key: str) -> tuple:
    return (1, 0, 0) if key == "later" else (0, *[int(x) for x in key.split(".")])


def phase_label(key: str) -> str:
    """"8.1-2" -> "8.1 phase 2"."""
    version, _, n = key.partition("-")
    return f"{version} phase {n}" if n else version


def parse_phases(plan_dir: Path) -> tuple[list[dict], dict, list[dict]]:
    """The version files: versions, their phases (themes) and threads, and each item's row."""
    versions, phases, rows = [], [], {}
    files = sorted(plan_dir.glob("version-*.md"), key=lambda p: version_sort(version_key(p)))
    for path in files:
        text = path.read_text(encoding="utf-8")
        version = version_key(path)
        title = re.sub(r"^Version \S+:\s*", "", text.splitlines()[0].lstrip("# ").strip())
        goal = ""
        if "## Goal" in text:
            goal = text.split("## Goal", 1)[1].split("\n## ", 1)[0].strip()
        entry = {"key": version, "file": path.name, "title": title, "goal": goal, "phases": []}
        versions.append(entry)
        for ph in re.finditer(r"^## Phase (\d+): (.+?)\n(.*?)(?=^## |\Z)", text, re.S | re.M):
            n, theme, body = int(ph.group(1)), ph.group(2).strip(), ph.group(3)
            key = f"{version}-{n}"
            groups = []
            for g in re.finditer(r"^### (.+?)\n(.*?)(?=^### |\Z)", body, re.S | re.M):
                name, chunk = g.group(1).strip(), g.group(2)
                ids = []
                for line in chunk.splitlines():
                    if not line.startswith("| "):
                        continue
                    cells = [c.strip() for c in line.strip().strip("|").split("|")]
                    if not ID_RE.match(cells[0]):
                        continue
                    ids.append(cells[0])
                    rows[cells[0]] = {
                        "phase": key, "version": version, "thread": name,
                        "needs": split_ids(cells[2]) if len(cells) > 2 else [],
                        "note": cells[3] if len(cells) > 3 else "",
                    }
                intro = re.split(r"\n\|", chunk)[0].strip()
                groups.append({"name": name, "ids": ids, "intro": intro})
            intro = re.split(r"^### ", body, maxsplit=1, flags=re.M)[0].strip()
            phases.append({"key": key, "version": version, "n": n, "theme": theme, "intro": intro,
                           "label": phase_label(key), "title": title, "file": path.name, "groups": groups})
            entry["phases"].append(key)
    return phases, rows, versions


def parse_decisions(notes: str) -> list[str]:
    if "### Decisions for Boss" not in notes:
        return []
    chunk = notes.split("### Decisions for Boss", 1)[1]
    chunk = re.split(r"\n#{2,3} ", chunk, maxsplit=1)[0]
    return [line[2:].strip() for line in chunk.splitlines() if line.startswith("- ")]


def parse_number_map(text: str) -> dict:
    """Every ID ever issued: its title and recorded status."""
    out = {}
    # The "Reverse lookup" table holds most IDs; a few (TEST.1 onward) sit in
    # a later table with the same four columns, so read every such row and
    # let the first one for an ID win.
    for line in text.splitlines():
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) >= 4 and ID_RE.match(cells[0]) and cells[0] not in out:
            out[cells[0]] = {"title": cells[1], "old": cells[2], "status": cells[3]}
    return out


def status_bucket(status: str) -> str:
    s = status.lower()
    if s.startswith("open"):
        return "open"
    if s.startswith("done"):
        return "done"
    if s.startswith(("dropped", "merged", "replaced", "folded")):
        return "closed"
    return "other"


ISSUE_URL = "https://github.com/dwhagar/planetGen/issues/"
ISSUE_ROW_RE = re.compile(r"^\|\s*\[#(\d+)\]\([^)]*\)\s*\|\s*(open|closed)\s*\|([^|]*)\|([^|]*)\|\s*$")


def parse_issue_groups(text: str) -> tuple[dict[str, list[int]], dict[int, dict]]:
    """The issue table: item -> issue numbers, and issue number -> title and state."""
    by_item: dict[str, list[int]] = {}
    issues: dict[int, dict] = {}
    for line in text.splitlines():
        m = ISSUE_ROW_RE.match(line)
        if not m:
            continue
        n = int(m.group(1))
        issues[n] = {"title": m.group(3).strip(), "state": m.group(2)}
        for item_id in re.findall(r"[A-Z]+\.\d+", m.group(4)):
            if n not in by_item.setdefault(item_id, []):
                by_item[item_id].append(n)
    return by_item, issues


LEDGER_ROW_RE = re.compile(r"^\|\s*([A-Z]+\.\d+)\s*\|\s*(" + VERSION_RE + r")\s*\|([^|]*)\|\s*$")
LEDGER_BASELINE_RE = re.compile(r"^Baseline: PR #(\d+)\.", re.M)


def parse_ledger(text: str) -> tuple[dict[str, str], int]:
    """Finished item -> the version it belonged to, and the baseline PR number."""
    rows = {m.group(1): m.group(2) for m in map(LEDGER_ROW_RE.match, text.splitlines()) if m}
    base = LEDGER_BASELINE_RE.search(text)
    return rows, int(base.group(1)) if base else 0


def parse_ledger_rows(text: str) -> list[tuple[str, str, str]]:
    """Every ledger row as (item, version, thread), in file order."""
    return [(m.group(1), m.group(2), m.group(3).strip()) for m in map(LEDGER_ROW_RE.match, text.splitlines()) if m]


def build_finished(items: dict, rows: list[tuple[str, str, str]], phases: list[dict]) -> tuple[list[dict], dict]:
    """The finished subtasks of every thread (Boss, 2026-10-21, 9:15 PM EDT): each thread of a version file gets the
    finished items whose ledger row names it, to be shown checked off and crossed out above its open items; finished
    items whose thread no longer exists in the version file are listed per version under their old thread name.
    Dropped items (superseded, folded or not reproduced) are listed too, marked dropped, so a thread reads whole."""
    live = {(p["version"], g["name"]) for p in phases for g in p["groups"]}
    by_thread: dict[tuple[str, str], list[dict]] = {}
    for item_id, version, thread in rows:
        it = items.get(item_id)
        if not it or it["status"] not in ("done", "dropped"):
            continue
        by_thread.setdefault((version, thread), []).append(
            {"id": item_id, "title": it["title"], "dropped": it["status"] == "dropped"})
    plan_phases = [{**p, "groups": [{**g, "done": by_thread.get((p["version"], g["name"]), [])} for g in p["groups"]]}
                   for p in phases]
    retired: dict[str, list[dict]] = {}
    for (version, thread), done in by_thread.items():
        if (version, thread) not in live:
            retired.setdefault(version, []).append({"thread": thread, "items": done})
    for v in retired.values():
        v.sort(key=lambda r: r["thread"].lower())
    return plan_phases, retired


def build_progress(items: dict, ledger: dict[str, str], baseline: int, versions: list[dict], phases: list[dict]) -> dict:
    """Finished, open and dropped items per version, with the percentage finished and each phase's open count."""
    missing = []
    for item_id, it in items.items():
        pr = re.search(r"PR #(\d+)", it["statusText"])
        if it["status"] == "done" and pr and int(pr.group(1)) > baseline and item_id not in ledger:
            missing.append(item_id)
    if missing:
        raise SystemExit("docs/plan/progress-ledger.md has no row for " + ", ".join(sorted(missing, key=id_key)) +
                         ": add one (ID | version | thread) for each, copied from the version file row before it was deleted")
    keys = [VERSION0] + [v["key"] for v in versions]
    titles = {VERSION0: VERSION0_TITLE, **{v["key"]: v["title"] for v in versions}}
    out = {k: {"key": k, "title": titles.get(k, ""), "done": 0, "open": 0, "dropped": 0} for k in keys}
    before = 0
    for item_id, it in items.items():
        if it["status"] == "open":
            if it["version"] in out:
                out[it["version"]]["open"] += 1
        elif item_id in ledger:
            out[ledger[item_id]]["done" if it["status"] == "done" else "dropped"] += 1
        elif it["status"] == "done":
            before += 1
    for t in out.values():
        t["total"] = t["done"] + t["open"]
        t["pct"] = round(100 * t["done"] / t["total"], 1) if t["total"] else 0.0
        t["phases"] = [{"key": p["key"], "label": p["label"], "theme": p["theme"],
                        "open": sum(1 for it in items.values() if it["status"] == "open" and it["phase"] == p["key"])}
                       for p in phases if p["version"] == t["key"]]
    done, total = sum(t["done"] for t in out.values()), sum(t["total"] for t in out.values())
    return {"versions": [out[k] for k in keys], "before": before,
            "done": done, "total": total, "pct": round(100 * done / total, 1) if total else 0.0}


def build_blockers(items: dict, phases: list[dict], top: int = 5) -> dict:
    """Boss (2026-10-16): the top blockers overall and for each phase: the open items that hold up the most others.

    An item's weight is the number of open items that wait on it directly or through a chain; ties break on the
    direct count, then the ID. An item that waits on, or is waited on by, a blocker already listed, or holds up almost the same items (80% overlap), is part of that
    blocker's chain and is not listed again (Boss, 2026-10-16: an umbrella and its pieces count once); the row
    names those chain members in `chain`."""
    open_ids = {i for i, it in items.items() if it["status"] == "open"}
    waiting: dict[str, set] = {}

    def below(i: str, stack: tuple = ()) -> set:
        if i in waiting:
            return waiting[i]
        out: set = set()
        for u in items[i]["unblocks"]:
            if u in open_ids and u not in stack:
                out.add(u)
                out |= below(u, stack + (i,))
        if not stack:
            waiting[i] = out
        return out

    rows = {i: {"id": i, "blocks": len(below(i)),
                "direct": sum(1 for u in items[i]["unblocks"] if u in open_ids)} for i in open_ids}

    def related(a: str, b: str) -> bool:
        """One chain: one waits on the other, or the two hold up almost the same set of items."""
        if b in below(a) or a in below(b):
            return True
        shared = len(below(a) & below(b))
        return shared >= 0.8 * min(len(below(a)), len(below(b))) > 0

    def best(ids) -> list[dict]:
        ranked = sorted((rows[i] for i in ids if rows[i]["blocks"]), key=lambda r: (-r["blocks"], -r["direct"], id_key(r["id"])))
        picked: list[dict] = []
        for r in ranked:
            if len(picked) < top and not any(related(r["id"], p["id"]) for p in picked):
                picked.append(dict(r, chain=[]))
        for r in ranked:
            if all(r["id"] != p["id"] for p in picked):
                for p in picked:
                    if related(r["id"], p["id"]):
                        p["chain"].append(r["id"])
                        break
        for p in picked:
            p["chain"] = p["chain"][:6]
        return picked

    return {"overall": best(open_ids),
            "phases": [{"key": p["key"], "label": p["label"], "theme": p["theme"],
                        "rows": best(i for i in open_ids if items[i]["phase"] == p["key"])} for p in phases]}


def check_order(items: dict) -> list[dict]:
    """Boss (2026-10-16): a blocker sits in the same thread (earlier) or in an earlier phase, never later."""
    order = {}
    for it in items.values():
        if it["status"] == "open" and it["phase"]:
            order[it["id"]] = (version_sort(it["version"]), int(it["phase"].partition("-")[2] or 0))
    bad = []
    for item_id, it in items.items():
        if item_id not in order:
            continue
        for dep in it["deps"]:
            if dep not in order or items[dep]["thread"] == it["thread"] and items[dep]["phase"] == it["phase"]:
                continue
            if order[dep] >= order[item_id]:
                bad.append({"id": item_id, "phase": it["phase"], "needs": dep, "needsPhase": items[dep]["phase"]})
    return bad


ROADMAP_RELEASE_RE = re.compile(r"^Release (" + VERSION_RE + r"): (.*)$")
ROADMAP_FIELD_RE = re.compile(r"^(Status|Headline|Today|Change|Needs|Lanes|Alternative)(?: \(([^)]+)\))?:\s*(.*)$")


def parse_roadmap(text: str) -> dict:
    """docs/plan/release-roadmap.md: each version's headline and needs."""
    out = {"title": "", "intro": [], "how": [], "releases": [], "questions": []}
    cur, section, last = None, "intro", None
    for line in text.splitlines():
        if line.startswith("# "):
            out["title"] = line[2:].strip()
            continue
        if line.startswith("## "):
            title = line[3:].strip()
            found = ROADMAP_RELEASE_RE.match(title)
            if found:
                cur = {"version": found.group(1), "name": found.group(2), "fields": {},
                       "alternatives": []}
                out["releases"].append(cur)
                section = "release"
            else:
                cur, section = None, title
            last = None
            continue
        if not line.strip():
            continue
        if section == "intro":
            out["intro"].append(line.strip())
        elif section == "How to read this":
            out["how"].append(line[2:].strip() if line.startswith("- ") else line.strip())
        elif section == "Questions for Boss":
            out["questions"].append(line[2:].strip() if line.startswith("- ") else line.strip())
        elif section == "release" and cur is not None:
            field = ROADMAP_FIELD_RE.match(line)
            if line.startswith("  ") and last is not None:
                last["text"] += " " + line.strip()
            elif field:
                key, label, value = field.groups()
                if key == "Alternative":
                    last = {"label": label or "", "text": value}
                    cur["alternatives"].append(last)
                else:
                    last = {"text": value}
                    cur["fields"][key] = last
    return out


def build_roadmap(root: Path, items: dict, progress: dict, phases: list[dict], versions: list[dict]) -> dict:
    """The roadmap: each version's narrative from the markdown plus its phases and threads from the version files."""
    path = root / "docs" / "plan" / "release-roadmap.md"
    road = parse_roadmap(path.read_text(encoding="utf-8") if path.exists() else "")
    prog = {t["key"]: t for t in progress["versions"]}
    seen = set()
    for rel in road["releases"]:
        if rel["version"] in seen:
            raise SystemExit(f"docs/plan/release-roadmap.md: release {rel['version']} appears twice")
        seen.add(rel["version"])
        rel["fields"] = {k: v["text"].strip() for k, v in rel["fields"].items()}
        rel["alternatives"] = [{"label": a["label"], "text": a["text"].strip()} for a in rel["alternatives"]]
        ver = next((v for v in versions if v["key"] == rel["version"]), None)
        rel["title"] = ver["title"] if ver else rel["name"]
        rel["goal"] = ver["goal"] if ver else ""
        rel["phases"] = []
        open_ids = []
        for ph in (p for p in phases if p["version"] == rel["version"]):
            groups = []
            for g in ph["groups"]:
                rows = [{"id": i, "title": items[i]["title"], "bug": items[i]["bug"],
                         "deps": items[i]["deps"]} for i in g["ids"] if i in items and items[i]["status"] == "open"]
                open_ids += [r["id"] for r in rows]
                groups.append({"name": g["name"], "items": rows})
            rel["phases"].append({"key": ph["key"], "label": ph["label"], "theme": ph["theme"], "intro": ph["intro"],
                                  "groups": groups})
        t = prog.get(rel["version"], {"done": 0, "open": 0, "pct": 0.0, "dropped": 0})
        rel.update({
            "open": len(open_ids), "ready": sum(1 for i in open_ids if not items[i]["deps"]),
            "bugs": sum(1 for i in open_ids if items[i]["bug"]),
            "done": t["done"], "share": t["pct"] if t["done"] + len(open_ids) else None,
        })
    road["violations"] = check_order(items)
    road["unplaced"] = sorted((i for i, it in items.items() if it["status"] == "open" and not it["phase"]), key=id_key)
    road["unlisted"] = [v["key"] for v in versions if v["key"] not in seen]
    return road


def build_data(root: Path) -> dict:
    todo = parse_todo((root / "docs" / "TODO.md").read_text(encoding="utf-8"))
    phases, rows, versions = parse_phases(root / "docs" / "plan" / "releases")
    decisions = parse_decisions((root / "docs" / "plan" / "notes.md").read_text(encoding="utf-8"))
    numbers = parse_number_map((root / "docs" / "plan" / "todo-number-map.md").read_text(encoding="utf-8"))
    groups_file = root / "docs" / "plan" / "issue-groups.md"
    item_issues, issues = parse_issue_groups(groups_file.read_text(encoding="utf-8") if groups_file.exists() else "")

    items: dict[str, dict] = {}
    for item_id in sorted(set(todo) | set(numbers), key=id_key):
        t, n, r = todo.get(item_id), numbers.get(item_id, {}), rows.get(item_id, {})
        title = t["title"] if t else n.get("title", "")
        status = "open" if t else status_bucket(n.get("status", ""))
        if status == "open" and not t:
            status = "other"  # the map says open but TODO.md has no entry
        deps = []
        for x in (r.get("needs", []) + (t["prereqs"] if t else [])):
            if x != item_id and x not in deps and x in todo:
                deps.append(x)
        items[item_id] = {
            "id": item_id, "cat": item_id.split(".")[0], "title": title,
            "bug": title.rstrip().endswith("(bug)"),
            "status": status, "statusText": n.get("status", "open" if t else ""),
            "old": n.get("old", ""),
            "phase": r.get("phase", ""), "version": r.get("version", ""), "thread": r.get("thread", ""), "note": r.get("note", ""),
            "text": t["text"] if t else "", "parent": t["parent"] if t else "",
            "section": t["section"] if t else "",
            "design": [{"label": a, "href": b} for a, b in (t["design"] if t else [])],
            "deps": deps, "unblocks": [], "children": [],
            "issues": item_issues.get(item_id, []),
            "headline": t["headline"][0] if t and t["headline"] else "",
        }
    for it in items.values():
        for d in it["deps"]:
            items[d]["unblocks"].append(it["id"])
        if it["parent"] in items:
            items[it["parent"]]["children"].append(it["id"])
    for it in items.values():
        it["unblocks"].sort(key=id_key)

    open_ids = [i for i, it in items.items() if it["status"] == "open"]
    plan_items = {i: {
        "id": i, "title": items[i]["title"], "phase": items[i]["phase"], "version": items[i]["version"], "was": "",
        "bug": items[i]["bug"], "thread": items[i]["thread"], "deps": items[i]["deps"],
        "unblocks": items[i]["unblocks"], "note": items[i]["note"], "text": items[i]["text"],
        "headline": items[i]["headline"], "issues": items[i]["issues"],
    } for i in open_ids}
    headlines = build_headlines(items, todo)
    ledger_file = root / "docs" / "plan" / "progress-ledger.md"
    ledger, baseline = parse_ledger(ledger_file.read_text(encoding="utf-8") if ledger_file.exists() else "")
    progress = build_progress(items, ledger, baseline, versions, phases)
    plan_phases, retired = build_finished(items, parse_ledger_rows(ledger_file.read_text(encoding="utf-8") if ledger_file.exists() else ""), phases)
    roadmap = build_roadmap(root, items, progress, phases, versions)
    summary = (f"{len(open_ids)} open items, {sum(items[i]['bug'] for i in open_ids)} bugs; "
               f"{len(items)} IDs issued in all")
    return {
        "plan": {"items": plan_items, "phases": plan_phases, "retired": retired, "versions": versions, "decisions": decisions, "generated": summary,
                 "headlines": headlines, "issues": {str(k): v for k, v in issues.items()},
                 "progress": progress, "blockers": build_blockers(items, phases)},
        "roadmap": roadmap,
        "reference": {"items": items, "phases": [{"key": p["key"], "label": p["label"], "theme": p["theme"], "version": p["version"]} for p in phases],
                      "categories": CATEGORIES, "summary": summary, "headlines": headlines,
                      "issues": {str(k): v for k, v in issues.items()}},
    }


def build_headlines(items: dict, todo: dict) -> list[dict]:
    """Each headline feature with every open item it needs, in build order."""
    out = []
    for item_id, t in todo.items():
        if not t["headline"]:
            continue
        version, name = t["headline"]
        closure: dict[str, int] = {}

        def walk(i: str, depth_seen: tuple = ()) -> int:
            if i in closure:
                return closure[i]
            if i in depth_seen:
                return 0
            below = [walk(d, depth_seen + (i,)) for d in items[i]["deps"] if d in items]
            closure[i] = 1 + max(below) if below else 0
            return closure[i]

        walk(item_id)
        rows = []
        for i in sorted(closure, key=lambda x: (closure[x], id_key(x))):
            it = items[i]
            rows.append({"id": i, "title": it["title"], "phase": it["phase"], "bug": it["bug"],
                         "deps": [d for d in it["deps"] if d in closure],
                         "ready": not [d for d in it["deps"] if d in closure]})
        out.append({"version": version, "name": name, "id": item_id, "title": t["title"],
                    "remaining": len(rows), "ready": sum(r["ready"] for r in rows),
                    "bugs": sum(r["bug"] for r in rows), "items": rows})
    return sorted(out, key=lambda h: tuple(int(x) for x in h["version"].split(".")))


def _inline(text: str) -> str:
    """`text` as HTML: escaped, with `**bold**` and `code` spans."""
    text = html.escape(text, quote=False)
    text = re.sub(r"`([^`]+)`", r"<code>\1</code>", text)
    return re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", text)


def parse_lanes(text: str) -> dict:
    """`docs/plan/lanes.md`: the `Updated:` stamp, the table rows (`Lane | Doing now | Gate | Next`) and the `Merged:`
    line."""
    stamp = re.search(r"^Updated: (.+)$", text, re.M)
    merged = re.search(r"^Merged: (.+)$", text, re.M)
    rows = []
    for line in text.splitlines():
        cells = [cell.strip() for cell in line.strip().strip("|").split("|")] if line.startswith("|") else []
        if len(cells) == 4 and cells[0] != "Lane" and not set(cells[0]) <= {"-", ":", " "}:
            rows.append(cells)
    return {"updated": stamp.group(1).strip() if stamp else "", "rows": rows,
            "merged": merged.group(1).strip() if merged else ""}


def render_notes(text: str) -> str:
    """The sections of `docs/plan/execution-notes.md` below its first paragraph block as HTML: `##` and `###` headings,
    `- ` bullets, paragraphs (one wholly in `*...*` is muted). Everything before the first `##` is the file's own note."""
    text = text[text.index("\n## "):] if "\n## " in text else ""
    out, bullets = [], []

    def close():
        if bullets:
            out.append('<ul class="tight">' + "".join(f"<li>{_inline(b)}</li>" for b in bullets) + "</ul>")
            bullets.clear()

    for line in text.splitlines():
        line = line.rstrip()
        if line.startswith("- "):
            bullets.append(line[2:])
            continue
        close()
        if line.startswith("### "):
            out.append(f"<h3>{_inline(line[4:])}</h3>")
        elif line.startswith("## "):
            out.append(f"<h2>{_inline(line[3:])}</h2>")
        elif len(line) > 2 and line.startswith("*") and line.endswith("*") and not line.startswith("**"):
            out.append(f'<p class="muted">{_inline(line[1:-1])}</p>')
        elif line:
            out.append(f"<p>{_inline(line)}</p>")
    close()
    return "\n".join(out)


def build_execution_plan(root: Path, data: dict) -> str:
    """`docs/generated/execution-plan.html` (DOC.68): the lanes board from `docs/plan/lanes.md`, the headline panel, the open
    items by version and by 8.1 phase and the goal issues from the TODO data, the conflicts and needs-Boss sections from
    `docs/plan/execution-notes.md`, and the older decisions from `scripts/todo_docs/execution-history.snippet.html`.
    Nothing is read from the clock, so the page is the same on every run until a source changes."""
    esc = lambda value: html.escape(str(value), quote=True)
    reference, plan = data["reference"], data["plan"]
    items, issues = reference["items"], reference["issues"]
    opened = {i for i, item in items.items() if item["status"] == "open"}

    def idspan(item_id):
        item = items.get(item_id)
        if not item:
            return f'<span class="id">{esc(item_id)}</span>'
        title = esc(item["title"])
        if item["issues"]:
            number = item["issues"][0]
            return f'<a class="id" href="{ISSUES_URL}{number}" title="{title} (#{number})">{esc(item_id)}</a>'
        return f'<span class="id" title="{title}">{esc(item_id)}</span>'

    lanes = parse_lanes((root / "docs" / "plan" / "lanes.md").read_text(encoding="utf-8"))
    lane_html = "".join(f"<tr><td><b>{esc(row[0])}</b></td><td>{esc(row[1])}</td><td>{esc(row[2])}</td><td>{esc(row[3])}</td></tr>"
                        for row in lanes["rows"])
    headlines = ""
    for headline in plan["headlines"]:
        entries = "".join(
            f'<li>{idspan(row["id"])} {esc(row["title"])} <span class="muted">'
            f'{esc(row["phase"].replace("-", " phase ") if row["phase"] else "unplaced")} · '
            f'{"ready" if row["ready"] else "needs " + ", ".join(row["deps"])}</span></li>'
            for row in headline["items"])
        headlines += (f'<details><summary><b>{esc(headline["version"])}</b> {esc(headline["name"])} <span class="muted">'
                      f'({headline["remaining"]} open items needed, {headline["ready"]} ready to start)</span></summary>'
                      f'<ol class="tight">{entries}</ol></details>\n')
    per_version, bugs = {}, {}
    for item_id in opened:
        item = items[item_id]
        key = item["version"] or "unplaced"
        per_version[key] = per_version.get(key, 0) + 1
        bugs[key] = bugs.get(key, 0) + (1 if item["bug"] else 0)
    order = sorted(per_version, key=lambda v: (v == "later", v == "unplaced",
                                                [int(x) for x in v.split(".")] if v[0].isdigit() else [99]))
    version_rows = "".join(f'<tr><td>{esc(v)}</td><td class="k">{per_version[v]}</td><td class="k">{bugs[v]}</td></tr>'
                           for v in order)
    per_phase = {}
    for item_id in opened:
        if items[item_id]["version"] == "8.1":
            per_phase[items[item_id]["phase"]] = per_phase.get(items[item_id]["phase"], 0) + 1
    labels = {p["key"]: f'{p["label"]}: {p["theme"]}' for p in reference["phases"]}
    phase_rows = "".join(f'<tr><td>{esc(labels.get(k, k or "unplaced"))}</td><td class="k">{per_phase[k]}</td></tr>'
                         for k in sorted(per_phase, key=lambda k: (k == "", k or "")))
    goal_rows = ""
    for number in sorted(issues, key=int, reverse=True):
        issue = issues[number]
        everyone = [i for i, item in items.items() if int(number) in item["issues"]]
        members = [i for i in everyone if items[i]["status"] == "open"]
        if issue["state"] != "open" or not members:
            continue
        # Boss (2026-10-21, 9:15 PM EDT): the goal's subtask list shows the finished ones checked off and crossed out.
        finished = [i for i in everyone if items[i]["status"] == "done"]
        dropped = [i for i in everyone if items[i]["status"] == "dropped"]
        pct = round(100 * len(finished) / (len(finished) + len(members)), 1)
        chips = " ".join(f'<span class="id" title="{esc(items[i]["title"])}">{esc(i)}</span>' for i in members)
        chips += " " + " ".join(f'<span class="id dn" title="finished: {esc(items[i]["title"])}">&#10003; {esc(i)}</span>' for i in finished)
        chips += " " + " ".join(f'<span class="id dn dropped" title="dropped: {esc(items[i]["title"])}">&#10007; {esc(i)}</span>' for i in dropped)
        meter = (f'<span class="meter" role="progressbar" aria-label="#{number} progress" aria-valuemin="0" aria-valuemax="100" '
                 f'aria-valuenow="{pct}"><i style="width:{pct}%"></i></span>')
        goal_rows += (f'<tr><td><a href="{ISSUES_URL}{number}">#{number}</a></td><td>{esc(issue["title"])}</td>'
                      f'<td class="k">{len(members)}<div class="gp">{meter}<span>{len(finished)} of {len(finished) + len(members)} finished</span></div></td>'
                      f'<td><div class="chips">{chips}</div></td></tr>')
    goals = goal_rows.count("<tr>")
    bug_count = sum(1 for i in opened if items[i]["bug"])
    notes = render_notes((root / "docs" / "plan" / "execution-notes.md").read_text(encoding="utf-8"))
    history = (TEMPLATES / "execution-history.snippet.html").read_text(encoding="utf-8").rstrip()
    body = f'''<header>
  <div class="eyebrow">planetGen · version map 8.1 to 10.0, Core lanes 1 to 5 and Bug Hunt · updated {esc(lanes["updated"])}</div>
  <div class="meta"><a href="../index.html">Docs index</a> · <a href="version-plan.html">Version plan</a></div>
  <h1>Execution Plan</h1>
  <p class="lede">{len(opened)} open items, {bug_count} of them bugs. Version 8.1 (the infinite zoom map and the new generate scatters, halo and satellite galaxies included) is the headline in progress; its stamp comes when every item of its six phases is done. Five Core lanes and a bug lane build it in parallel, each through its own gate; the TODO lane keeps this page, the plan documents, the number map and the GitHub issues true, and every item ID below links to the issue that groups it with the other items of its goal.</p>
</header>

<h2>Where each lane is</h2>
<p class="muted">Updated {esc(lanes["updated"])} from the lanes' own status lines (<code>docs/plan/lanes.md</code>); times are the lanes' clocks as reported. A gate is the <code>scripts/test.py changed</code> run that has to be green before a lane's PR merges (Boss: merges need no human review).</p>
<div class="tw"><table><thead><tr><th>Lane</th><th>Doing</th><th>Gate</th><th>Next</th></tr></thead><tbody>{lane_html}</tbody></table></div>
<p class="note">Merged: {esc(lanes["merged"])}</p>

<h2>Headline features and what each needs</h2>
<p class="muted">Boss, 2026-10-10: every release has one headline feature. The open items each one needs, in build order (generated from docs/TODO.md; the TODO tree, version plan and TODO reference carry the same panel).</p>
{headlines}
<h2>Open items by version</h2>
<div class="split">
<div><div class="tw"><table><thead><tr><th>Version</th><th>Open items</th><th>Bugs</th></tr></thead><tbody>{version_rows}</tbody></table></div></div>
<div><div class="tw"><table><thead><tr><th>8.1 phase</th><th>Open items</th></tr></thead><tbody>{phase_rows}</tbody></table></div></div>
</div>

<h2>Goals and their GitHub issues</h2>
<p class="muted">{goals} open goal issues. Each groups the TODO items of one goal and is what we are working toward; a goal issue carries the status label "in progress" only while code for one of its items is being written, plus a dated comment with where it stands. Findings that belong to a goal (a failed test, a bug noticed in passing) are comments on its issue. Closed goals are on GitHub.</p>
<div class="tw"><table><thead><tr><th>Issue</th><th>Goal</th><th>Open items</th><th>Items</th></tr></thead><tbody>{goal_rows}</tbody></table></div>

{notes}

<h2>Earlier decisions</h2>
{history}
<p class="muted">Links: <a href="version-plan.html">version plan</a> · <a href="todo-tree.html">TODO tree</a> · <a href="todo-reference.html">every ID explained</a> · <a href="release-roadmap.html">release roadmap</a> · <a href="https://github.com/dwhagar/planetGen/issues?q=is%3Aissue+is%3Aopen+label%3Agoal">goal issues on GitHub</a></p>'''
    return EXECUTION_TEMPLATE.read_text(encoding="utf-8").replace("__BODY__", body)


def render(template: str, data: dict) -> str:
    payload = json.dumps(data, ensure_ascii=False, sort_keys=True).replace("</", "<\\/")
    snippet = (TEMPLATES / "headlines.snippet.html").read_text(encoding="utf-8")
    return template.replace("__HEADLINES__", snippet).replace("__DATA__", payload)


BANNER_RE = re.compile(r'\n<div class="updated-banner"[^>]*>Last updated: (.*?)</div>')


def banner(stamp: str) -> str:
    style = ("font:600 13px/1.4 system-ui,-apple-system,'Segoe UI',sans-serif;padding:6px 14px;text-align:center;"
             "background:#F6E7D6;color:#152030;border-bottom:1px solid #D8DEE7")
    return f'\n<div class="updated-banner" style="{style}">Last updated: {stamp}</div>'


def stamp_of(text: str) -> str:
    """The stamp a page carries, `""` if it has none."""
    found = BANNER_RE.search(text)
    return found.group(1) if found else ""


def with_stamp(page: str, stamp: str) -> str:
    """`page` with its banner set to `stamp`: just after `<body>`, or after the last `</style>` in a fragment."""
    page = BANNER_RE.sub("", page)
    mark = "<body>" if "<body>" in page else "</style>"
    at = page.index(mark) if mark == "<body>" else page.rindex(mark)
    at += len(mark)
    return page[:at] + banner(stamp) + page[at:]


def now_stamp() -> str:
    """The Eastern time now, for example `2026-10-31 4:05 PM EDT`."""
    from zoneinfo import ZoneInfo
    moment = datetime.datetime.now(ZoneInfo("America/New_York"))
    return f"{moment:%Y-%m-%d} {moment.hour % 12 or 12}:{moment:%M %p} {moment:%Z}"


def build(root: Path = ROOT) -> dict[Path, str]:
    data = build_data(root)
    out = {}
    for name, target in OUTPUTS.items():
        template = (TEMPLATES / name).read_text(encoding="utf-8")
        key = ("reference" if name.startswith("todo-reference") else
               "roadmap" if name.startswith("release-roadmap") else "plan")
        out[target] = render(template, data[key])
    out[EXECUTION_OUTPUT] = build_execution_plan(root, data)
    # each page keeps the stamp it already carries, so an unchanged page compares equal; `main` writes a new one
    return {target: with_stamp(page, stamp_of(target.read_text(encoding="utf-8")) if target.exists() else "")
            for target, page in out.items()}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--check", action="store_true", help="fail if a page is out of date")
    args = parser.parse_args(argv)
    pages = build()
    stale = [p for p, html in pages.items() if not p.exists() or p.read_text(encoding="utf-8") != html]
    if args.check:
        for p in stale:
            print(f"{p.relative_to(ROOT)} is out of date: run python scripts/build_todo_docs.py")
        return 1 if stale else 0
    now = now_stamp()
    for p, page_html in pages.items():
        if p in stale:
            page_html = with_stamp(page_html, now)
        p.write_text(page_html, encoding="utf-8")
        print(f"wrote {p.relative_to(ROOT)}" + (f" (updated {now})" if p in stale else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
