#!/usr/bin/env python3
"""
Builds the two browsable TODO pages from the TODO sources.

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

Sources
=======
- `docs/TODO.md`: open items, their full text, nesting and prerequisites.
- `docs/plan/phase-*.md`: phase, thread, "Needs" and notes per open item.
- `docs/plan/notes.md`: "Decisions for Boss".
- `docs/design/todo-number-map.md`: the status of every ID ever issued
  (its "Reverse lookup" table).

Usage
=====
    python scripts/build_todo_docs.py          # rewrite both pages
    python scripts/build_todo_docs.py --check  # fail if either page is stale

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
}
ID_RE = re.compile(r"^[A-Z]+\.\d+$")
ITEM_RE = re.compile(r"^(\s*)- \[[ x]\] \*\*([A-Z]+\.\d+) (.*?)\*\*\s*$")
PHASE_ORDER = ["0", "1", "2", "3", "3+"]
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
    return items


def phase_key(path: Path) -> str:
    m = re.search(r"phase-(0|1|2|3plus|3)-", path.name)
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


def build_data(root: Path) -> dict:
    todo = parse_todo((root / "docs" / "TODO.md").read_text(encoding="utf-8"))
    phases, rows = parse_phases(root / "docs" / "plan")
    decisions = parse_decisions((root / "docs" / "plan" / "notes.md").read_text(encoding="utf-8"))
    numbers = parse_number_map((root / "docs" / "design" / "todo-number-map.md").read_text(encoding="utf-8"))

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
    } for i in open_ids}
    summary = (f"{len(open_ids)} open items, {sum(items[i]['bug'] for i in open_ids)} bugs; "
               f"{len(items)} IDs issued in all")
    return {
        "plan": {"items": plan_items, "phases": phases, "decisions": decisions, "generated": summary},
        "reference": {"items": items, "phases": [{"key": p["key"], "title": p["title"]} for p in phases],
                      "categories": CATEGORIES, "summary": summary},
    }


def render(template: str, data: dict) -> str:
    payload = json.dumps(data, ensure_ascii=False, sort_keys=True).replace("</", "<\\/")
    return template.replace("__DATA__", payload)


def build(root: Path = ROOT) -> dict[Path, str]:
    data = build_data(root)
    out = {}
    for name, target in OUTPUTS.items():
        template = (TEMPLATES / name).read_text(encoding="utf-8")
        key = "plan" if name.startswith("tier-plan") else "reference"
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
