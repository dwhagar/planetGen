#!/usr/bin/env python3
"""
Builds the four browsable TODO pages from the TODO sources.

Why this exists
===============
Boss (2026-10-07 16:22Z) asked for the tier plan and a searchable page that
explains every TODO item in detail, rebuilt whenever `docs/TODO.md`
changes. Both pages are generated from the files the TODO tasks thread
already keeps, so one command refreshes them:

- `docs/plan/tier-plan.html`: every open item in build order, phase by
  phase and thread by thread (from `docs/plan/phase-*.md`).
- `docs/plan/todo-reference.html`: every ID ever issued, open or closed,
  with its full text, status, phase, prerequisites and what it unblocks,
  plus search and filters.
- `docs/plan/todo-tree.html`: every open item as a tree, each item under the
  prerequisite that finishes last, with phase and bug marks, search and
  filters (Boss 2026-10-09 07:03Z: "a TODO as a tree structure").
- `docs/plan/release-roadmap.html`: the proposed releases from v8.0 to v10.0,
  each with its headline, the open items it holds, what it needs, its share
  built and the suggestions not yet filed (Boss, 2026-10-16).

Sources
=======
- `docs/TODO.md`: open items, their full text, nesting and prerequisites.
- `docs/plan/phase-*.md`: phase, thread, "Needs" and notes per open item.
- `docs/plan/notes.md`: "Decisions for Boss".
- `docs/design/todo-number-map.md`: the status of every ID ever issued
  (its "Reverse lookup" table).
- `docs/plan/issue-groups.md`: the GitHub issues (goal groups of items and the
  issues Boss opened) and the items each one covers; every page links them.
- `docs/plan/progress-ledger.md`: the phase each finished item belonged to, for
  the progress bar on every tier of the tier plan (Boss, 2026-10-16).
- `docs/plan/release-roadmap.md`: the release roadmap's order, headlines and
  item IDs; open counts come from TODO.md on every build.

Usage
=====
    python scripts/build_todo_docs.py          # rewrite all four pages
    python scripts/build_todo_docs.py --check  # fail if a page is stale

The output carries no dates, so rebuilding unchanged sources gives
byte-identical pages. Templates live in `scripts/todo_docs/`.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TEMPLATES = ROOT / "scripts" / "todo_docs"
OUTPUTS = {
    "tier-plan.template.html": ROOT / "docs" / "plan" / "tier-plan.html",
    "todo-reference.template.html": ROOT / "docs" / "plan" / "todo-reference.html",
    "todo-tree.template.html": ROOT / "docs" / "plan" / "todo-tree.html",
    "release-roadmap.template.html": ROOT / "docs" / "plan" / "release-roadmap.html",
}
ID_RE = re.compile(r"^[A-Z]+\.\d+$")
HEADLINE_RE = re.compile(r"Headline:\s*(\d+\.\d+)\s*\(([^)]*)\)")
ITEM_RE = re.compile(r"^(\s*)- \[[ x]\] \*\*([A-Z]+\.\d+) (.*?)\*\*\s*$")
PHASE_ORDER = ["1", "2", "3", "3+"]
TIER_KEYS = ["0", "1", "2", "3", "3+"]
PHASE0_TITLE = "Every bug and the groundwork (finished 2026-10-09)"
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
        it["headline"] = (found.group(1), found.group(2).strip()) if found else None
    return items


def phase_key(path: Path) -> str:
    m = re.search(r"phase-(1|2|3plus|3)-", path.name)
    return {"3plus": "3+"}.get(m.group(1), m.group(1))


def parse_phases(plan_dir: Path) -> tuple[list[dict], dict]:
    phases, rows = [], {}
    files = sorted(plan_dir.glob("phase-*.md"), key=lambda p: PHASE_ORDER.index(phase_key(p)))
    for path in files:
        text = path.read_text(encoding="utf-8")
        key = phase_key(path)
        title = text.splitlines()[0].lstrip("# ").strip()
        goal = ""
        if "## Goal" in text:
            goal = text.split("## Goal", 1)[1].split("\n## ", 1)[0].strip()
        groups = []
        for g in re.finditer(r"^### (.+?)\n(.*?)(?=^### |^## |\Z)", text, re.S | re.M):
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
                    "phase": key, "thread": name,
                    "needs": split_ids(cells[2]) if len(cells) > 2 else [],
                    "note": cells[3] if len(cells) > 3 else "",
                }
            intro = re.split(r"\n\|", chunk)[0].strip()
            groups.append({"name": name, "ids": ids, "intro": intro})
        questions = ""
        if "## Open questions for Boss" in text:
            questions = text.split("## Open questions for Boss", 1)[1].strip()
        phases.append({"key": key, "file": path.name, "title": title, "goal": goal,
                       "groups": groups, "questions": questions})
    return phases, rows


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


LEDGER_ROW_RE = re.compile(r"^\|\s*([A-Z]+\.\d+)\s*\|\s*(0|1|2|3\+?)\s*\|([^|]*)\|\s*$")
LEDGER_BASELINE_RE = re.compile(r"^Baseline: PR #(\d+)\.", re.M)


def parse_ledger(text: str) -> tuple[dict[str, str], int]:
    """Finished item -> the phase it belonged to, and the baseline PR number."""
    rows = {m.group(1): m.group(2) for m in map(LEDGER_ROW_RE.match, text.splitlines()) if m}
    base = LEDGER_BASELINE_RE.search(text)
    return rows, int(base.group(1)) if base else 0


def build_progress(items: dict, ledger: dict[str, str], baseline: int, phases: list[dict]) -> dict:
    """Finished, open and dropped items per tier, with the percentage finished."""
    missing = []
    for item_id, it in items.items():
        pr = re.search(r"PR #(\d+)", it["statusText"])
        if it["status"] == "done" and pr and int(pr.group(1)) > baseline and item_id not in ledger:
            missing.append(item_id)
    if missing:
        raise SystemExit("docs/plan/progress-ledger.md has no row for " + ", ".join(sorted(missing, key=id_key)) +
                         ": add one (ID | phase | thread) for each, copied from the phase file row before it was deleted")
    titles = {"0": PHASE0_TITLE}
    titles.update({p["key"]: re.sub(r"^Phase [0-9+]+:\s*", "", p["title"]) for p in phases})
    tiers = {k: {"key": k, "title": titles.get(k, ""), "done": 0, "open": 0, "dropped": 0} for k in TIER_KEYS}
    before = 0
    for item_id, it in items.items():
        if it["status"] == "open":
            if it["phase"] in tiers:
                tiers[it["phase"]]["open"] += 1
        elif item_id in ledger:
            tiers[ledger[item_id]]["done" if it["status"] == "done" else "dropped"] += 1
        elif it["status"] == "done":
            before += 1
    for t in tiers.values():
        t["total"] = t["done"] + t["open"]
        t["pct"] = round(100 * t["done"] / t["total"], 1) if t["total"] else 0.0
    done, total = sum(t["done"] for t in tiers.values()), sum(t["total"] for t in tiers.values())
    return {"tiers": [tiers[k] for k in TIER_KEYS], "before": before,
            "done": done, "total": total, "pct": round(100 * done / total, 1) if total else 0.0}


ROADMAP_RELEASE_RE = re.compile(r"^Release (\d+\.\d+): (.*)$")
ROADMAP_FIELD_RE = re.compile(r"^(Status|Headline|Today|Change|Needs|Lanes|Built before|Items|Alternative)(?: \(([^)]+)\))?:\s*(.*)$")
ROADMAP_SUGGEST_RE = re.compile(r"^- (S\d+\.\d+-[a-z]): (.*)$")


def parse_roadmap(text: str) -> dict:
    """docs/plan/release-roadmap.md: the proposed releases, their item IDs and suggestions."""
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
                cur = {"version": found.group(1), "name": found.group(2), "fields": {}, "groups": [],
                       "alternatives": [], "suggested": []}
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
            sug = ROADMAP_SUGGEST_RE.match(line)
            field = ROADMAP_FIELD_RE.match(line)
            if sug:
                cur["suggested"].append({"code": sug.group(1), "text": sug.group(2)})
                last = None
            elif line.startswith("  ") and last is not None:
                last["text"] += " " + line.strip()
            elif field:
                key, label, value = field.groups()
                if key == "Items":
                    last = {"label": label or "", "text": value}
                    cur["groups"].append(last)
                elif key == "Alternative":
                    last = {"label": label or "", "text": value}
                    cur["alternatives"].append(last)
                else:
                    last = {"text": value}
                    cur["fields"][key] = last
            elif line.strip() == "Suggested:":
                last = None
    return out


def build_roadmap(root: Path, items: dict, progress: dict) -> dict:
    """The roadmap with every listed item resolved against TODO.md, counts and order checks."""
    path = root / "docs" / "plan" / "release-roadmap.md"
    road = parse_roadmap(path.read_text(encoding="utf-8")) if path.exists() else parse_roadmap("")
    rank = {r["version"]: n for n, r in enumerate(road["releases"])}
    placed: dict[str, str] = {}
    for rel in road["releases"]:
        fields = {k: v["text"].strip() for k, v in rel["fields"].items()}
        built = re.match(r"\d+", fields.get("Built before", ""))
        rel["fields"] = fields
        rel["built"] = int(built.group(0)) if built else None
        for grp in rel["groups"]:
            ids = split_ids(grp.pop("text"))
            for item_id in ids:
                if item_id not in items:
                    raise SystemExit(f"docs/plan/release-roadmap.md: {item_id} in release {rel['version']} is not an item")
                if item_id in placed:
                    raise SystemExit(f"docs/plan/release-roadmap.md: {item_id} is in both {placed[item_id]} and {rel['version']}")
                placed[item_id] = rel["version"]
            grp["ids"] = ids
        rel["alternatives"] = [{"label": a["label"], "text": a["text"].strip()} for a in rel["alternatives"]]
    tiers = {t["key"]: t for t in progress["tiers"]}
    violations, unplaced = [], []
    for rel in road["releases"]:
        open_ids = []
        for grp in rel["groups"]:
            rows = []
            for item_id in grp["ids"]:
                it = items[item_id]
                is_open = it["status"] == "open"
                if is_open:
                    open_ids.append(item_id)
                    for dep in it["deps"]:
                        if dep in placed and rank[placed[dep]] > rank[rel["version"]]:
                            violations.append({"id": item_id, "release": rel["version"], "needs": dep, "needsRelease": placed[dep]})
                rows.append({"id": item_id, "title": it["title"], "phase": it["phase"], "bug": it["bug"],
                             "open": is_open, "deps": it["deps"] if is_open else []})
            grp["items"] = rows
            del grp["ids"]
        done = sum(1 for g in rel["groups"] for r in g["items"] if not r["open"])
        mix: dict[str, int] = {}
        for g in rel["groups"]:
            for r in g["items"]:
                if r["open"]:
                    mix[r["phase"] or "?"] = mix.get(r["phase"] or "?", 0) + 1
        finished = (rel["built"] or 0) + done
        total = finished + len(open_ids)
        rel.update({
            "open": len(open_ids), "ready": sum(1 for g in rel["groups"] for r in g["items"] if r["open"] and not r["deps"]),
            "bugs": sum(1 for g in rel["groups"] for r in g["items"] if r["open"] and r["bug"]),
            "doneListed": done, "mix": mix,
            "share": round(100 * finished / total, 1) if rel["built"] is not None and total else None,
        })
    for item_id, it in items.items():
        if it["status"] == "open" and item_id not in placed:
            unplaced.append(item_id)
    road["violations"], road["unplaced"] = violations, sorted(unplaced, key=id_key)
    road["tierShare"] = {k: t["pct"] for k, t in tiers.items()}
    return road



def build_data(root: Path) -> dict:
    todo = parse_todo((root / "docs" / "TODO.md").read_text(encoding="utf-8"))
    phases, rows = parse_phases(root / "docs" / "plan")
    decisions = parse_decisions((root / "docs" / "plan" / "notes.md").read_text(encoding="utf-8"))
    numbers = parse_number_map((root / "docs" / "design" / "todo-number-map.md").read_text(encoding="utf-8"))
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
            "phase": r.get("phase", ""), "thread": r.get("thread", ""), "note": r.get("note", ""),
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
        "id": i, "title": items[i]["title"], "phase": items[i]["phase"], "was": "",
        "bug": items[i]["bug"], "thread": items[i]["thread"], "deps": items[i]["deps"],
        "unblocks": items[i]["unblocks"], "note": items[i]["note"], "text": items[i]["text"],
        "headline": items[i]["headline"], "issues": items[i]["issues"],
    } for i in open_ids}
    headlines = build_headlines(items, todo)
    ledger_file = root / "docs" / "plan" / "progress-ledger.md"
    ledger, baseline = parse_ledger(ledger_file.read_text(encoding="utf-8") if ledger_file.exists() else "")
    progress = build_progress(items, ledger, baseline, phases)
    roadmap = build_roadmap(root, items, progress)
    summary = (f"{len(open_ids)} open items, {sum(items[i]['bug'] for i in open_ids)} bugs; "
               f"{len(items)} IDs issued in all")
    return {
        "plan": {"items": plan_items, "phases": phases, "decisions": decisions, "generated": summary,
                 "headlines": headlines, "issues": {str(k): v for k, v in issues.items()},
                 "progress": progress},
        "roadmap": roadmap,
        "reference": {"items": items, "phases": [{"key": p["key"], "title": p["title"]} for p in phases],
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


def render(template: str, data: dict) -> str:
    payload = json.dumps(data, ensure_ascii=False, sort_keys=True).replace("</", "<\\/")
    snippet = (TEMPLATES / "headlines.snippet.html").read_text(encoding="utf-8")
    return template.replace("__HEADLINES__", snippet).replace("__DATA__", payload)


def build(root: Path = ROOT) -> dict[Path, str]:
    data = build_data(root)
    out = {}
    for name, target in OUTPUTS.items():
        template = (TEMPLATES / name).read_text(encoding="utf-8")
        key = ("reference" if name.startswith("todo-reference") else
               "roadmap" if name.startswith("release-roadmap") else "plan")
        out[target] = render(template, data[key])
    return out


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
    for p, html in pages.items():
        p.write_text(html, encoding="utf-8")
        print(f"wrote {p.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
