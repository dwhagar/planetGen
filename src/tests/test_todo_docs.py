"""The TODO pages in docs/plan/ match their sources.

`scripts/build_todo_docs.py` builds docs/plan/version-plan.html,
docs/plan/todo-reference.html, docs/plan/todo-tree.html and
docs/plan/release-roadmap.html from docs/TODO.md, the version files, the
plan notes and the number map; docs/plan/execution-plan.html (DOC.68) comes from the same data plus
docs/plan/lanes.md, docs/plan/execution-notes.md and the history snippet. Boss (2026-10-07) asked for them to be rebuilt
with every TODO.md change, so a change that forgets the rebuild fails here
with the command to run. Boss (2026-10-16) turned the tiers into a version
map: versions with phases, a blocker never in a later phase.
"""

import importlib.util
import os

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
_spec = importlib.util.spec_from_file_location("build_todo_docs", os.path.join(ROOT, "scripts", "build_todo_docs.py"))
build_todo_docs = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(build_todo_docs)


def test_pages_are_up_to_date():
    stale = [str(p.relative_to(build_todo_docs.ROOT)) for p, html in build_todo_docs.build().items()
             if not p.exists() or p.read_text(encoding="utf-8") != html]
    assert not stale, f"{stale} out of date: run python scripts/build_todo_docs.py"


def test_every_open_item_is_in_the_reference_with_its_text():
    data = build_todo_docs.build_data(build_todo_docs.ROOT)
    todo = build_todo_docs.parse_todo(open(os.path.join(ROOT, "docs", "TODO.md"), encoding="utf-8").read())
    items = data["reference"]["items"]
    for item_id, it in todo.items():
        assert items[item_id]["status"] == "open"
        assert items[item_id]["text"] == it["text"]
    assert set(data["plan"]["items"]) == set(todo)


def test_needs_and_unblocks_mirror_each_other():
    items = build_todo_docs.build_data(build_todo_docs.ROOT)["reference"]["items"]
    for item_id, it in items.items():
        for dep in it["deps"]:
            assert item_id in items[dep]["unblocks"]


def test_progress_counts_every_open_item_once_and_matches_the_ledger():
    data = build_todo_docs.build_data(build_todo_docs.ROOT)["plan"]
    progress = data["progress"]
    assert sum(t["open"] for t in progress["versions"]) == len(data["items"])
    for t in progress["versions"]:
        assert t["total"] == t["done"] + t["open"]
        assert t["pct"] == (round(100 * t["done"] / t["total"], 1) if t["total"] else 0.0)
        assert sum(p["open"] for p in t["phases"]) == t["open"]
    assert progress["versions"][0]["key"] == "8.0" and progress["versions"][0]["open"] == 0


def test_every_open_item_sits_in_one_phase_of_one_version():
    data = build_todo_docs.build_data(build_todo_docs.ROOT)["plan"]
    phases = {p["key"]: p for p in data["phases"]}
    placed = [i for p in data["phases"] for g in p["groups"] for i in g["ids"]]
    assert len(placed) == len(set(placed)), "an item is listed twice in the version files"
    for item_id, it in data["items"].items():
        assert it["phase"] in phases, f"{item_id} is in no phase"
        assert it["version"] == phases[it["phase"]]["version"]


def test_a_finished_item_without_a_ledger_row_fails_the_build():
    items = {"GEN.1": {"status": "done", "statusText": "done, PR #9999", "phase": "", "version": ""}}
    versions = [{"key": "8.1", "title": "A", "phases": []}]
    try:
        build_todo_docs.build_progress(items, {}, 1000, versions, [])
    except SystemExit as err:
        assert "GEN.1" in str(err) and "progress-ledger.md" in str(err)
    else:
        raise AssertionError("a retired item with no ledger row must fail the build")
    done = build_todo_docs.build_progress(items, {"GEN.1": "8.1"}, 1000, versions, [])
    assert [t["done"] for t in done["versions"] if t["key"] == "8.1"] == [1]


def test_the_release_roadmap_is_in_build_order_with_a_headline_per_release():
    data = build_todo_docs.build_data(build_todo_docs.ROOT)
    road = data["roadmap"]
    assert road["violations"] == [], f"items placed before something they need: {road['violations']}"
    assert road["unplaced"] == [] and road["unlisted"] == []
    versions = [r["version"] for r in road["releases"]]
    assert versions[0] == "8.1" and versions[-1] == "later" and len(set(versions)) == len(versions)
    assert "10.0" in versions, "Boss (2026-10-16, 1:43 PM EDT): 10.0 is roles and permissions"
    for rel in road["releases"]:
        assert rel["fields"].get("Headline"), f"release {rel['version']} has no headline"
        prefix = "SL-" if rel["version"] == "later" else f"S{rel['version']}-"
        for sug in rel["suggested"]:
            assert sug["code"].startswith(prefix), sug


def test_every_headline_tag_names_a_version_of_the_map():
    data = build_todo_docs.build_data(build_todo_docs.ROOT)
    keys = {v["key"] for v in data["plan"]["versions"]}
    for h in data["plan"]["headlines"]:
        assert h["version"] in keys | {"8.0"}, h


def _roadmap_error(tmp_path, text, items):
    plan = tmp_path / "docs" / "plan"
    plan.mkdir(parents=True)
    (plan / "release-roadmap.md").write_text(text, encoding="utf-8")
    try:
        build_todo_docs.build_roadmap(tmp_path, items, {"versions": []}, [], [])
    except SystemExit as err:
        return str(err)
    return ""


def test_the_release_roadmap_refuses_a_version_twice(tmp_path):
    twice = "## Release 8.1: A\nHeadline: x\n\n## Release 8.1: B\nHeadline: y\n"
    assert "8.1" in _roadmap_error(tmp_path, twice, {})


def test_a_blocker_in_a_later_phase_is_reported():
    def items(*rows):
        return {i: {"id": i, "status": "open", "phase": ph, "version": ph.partition("-")[0],
                    "thread": thread, "deps": deps} for i, ph, deps, thread in rows}
    ok = items(("GEN.1", "8.1-1", [], "t"), ("GEN.2", "8.1-2", ["GEN.1"], "t"), ("GEN.3", "8.2-1", ["GEN.2"], "t"))
    assert build_todo_docs.check_order(ok) == []
    same_phase = items(("GEN.1", "8.1-1", [], "a"), ("GEN.2", "8.1-1", ["GEN.1"], "b"))
    assert [b["id"] for b in build_todo_docs.check_order(same_phase)] == ["GEN.2"]
    later = items(("GEN.1", "8.1-2", [], "t"), ("GEN.2", "8.1-1", ["GEN.1"], "t"))
    assert [b["needs"] for b in build_todo_docs.check_order(later)] == ["GEN.1"]


def test_the_tree_lists_the_top_five_blockers_overall_and_per_phase():
    plan = build_todo_docs.build_data(build_todo_docs.ROOT)["plan"]
    blockers = plan["blockers"]
    assert 0 < len(blockers["overall"]) <= 5
    weights = [r["blocks"] for r in blockers["overall"]]
    assert weights == sorted(weights, reverse=True)
    assert [p["key"] for p in blockers["phases"]] == [p["key"] for p in plan["phases"]]
    for p in blockers["phases"]:
        assert len(p["rows"]) <= 5
        for r in p["rows"]:
            assert plan["items"][r["id"]]["phase"] == p["key"] and r["blocks"] >= 1


def test_blocker_weight_counts_a_chain_once_per_item():
    def item(*unblocks):
        return {"status": "open", "unblocks": list(unblocks), "phase": "8.1-1"}
    items = {"A.1": item("A.2", "A.3"), "A.2": item("A.4"), "A.3": item("A.4"), "A.4": item()}
    rows = {r["id"]: r for r in build_todo_docs.build_blockers(items, [])["overall"]}
    assert rows["A.1"]["blocks"] == 3 and rows["A.1"]["direct"] == 2 and "A.4" not in rows


def test_blocker_ranking_lists_a_chain_once():
    def item(*unblocks):
        return {"status": "open", "unblocks": list(unblocks), "phase": "8.1-1"}
    items = {"A.1": item("A.3"), "A.2": item("A.3"), "A.3": item("A.4", "A.5"), "A.4": item(), "A.5": item(),
             "B.1": item("B.2"), "B.2": item()}
    rows = build_todo_docs.build_blockers(items, [])["overall"]
    ids = [r["id"] for r in rows]
    assert ids == ["A.1", "B.1"] and rows[0]["chain"] == ["A.2", "A.3"]


# --- The execution plan page (DOC.68) --------------------------------------------------------------------------

def test_the_execution_plan_is_generated_and_keeps_what_the_docs_index_needs():
    page = build_todo_docs.build()[build_todo_docs.EXECUTION_OUTPUT]
    assert '<meta name="description" content="' in page
    assert '<a href="../index.html">Docs index</a>' in page
    assert page == build_todo_docs.build()[build_todo_docs.EXECUTION_OUTPUT], "the page must not depend on the clock"
    assert "__BODY__" not in page


def test_the_lanes_file_gives_the_stamp_the_rows_and_the_merged_line():
    lanes = build_todo_docs.parse_lanes("""# Lanes board

Updated: 2026-10-18 1:02 AM EDT

| Lane | Doing now | Gate | Next |
|---|---|---|---|
| Core 1 | a | none | b |
| Bug Hunt | c | gate 2 | d |

Merged: GEN.1 (#1), GEN.2 (#2).
""")
    assert lanes["updated"] == "2026-10-18 1:02 AM EDT"
    assert lanes["rows"] == [["Core 1", "a", "none", "b"], ["Bug Hunt", "c", "gate 2", "d"]]
    assert lanes["merged"] == "GEN.1 (#1), GEN.2 (#2)."
    real = build_todo_docs.parse_lanes((build_todo_docs.ROOT / "docs" / "plan" / "lanes.md").read_text(encoding="utf-8"))
    assert real["updated"] and real["merged"] and len(real["rows"]) >= 5


def test_the_notes_become_headings_bullets_and_paragraphs_with_escaping():
    page = build_todo_docs.render_notes("""# Notes

An opening note that is not shown.

## Needs Boss

*Muted lead.*

- **GEN.1:** a <script> and `code`
- two

### Sub

A paragraph.
""")
    assert "opening note" not in page
    assert "<h2>Needs Boss</h2>" in page and "<h3>Sub</h3>" in page
    assert '<p class="muted">Muted lead.</p>' in page
    assert "<li><b>GEN.1:</b> a &lt;script&gt; and <code>code</code></li><li>two</li>" in page
    assert "<p>A paragraph.</p>" in page
