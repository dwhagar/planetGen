"""The two TODO pages in docs/plan/ match their sources.

`scripts/build_todo_docs.py` builds docs/plan/tier-plan.html and
docs/plan/todo-reference.html from docs/TODO.md, the phase plans, the plan
notes and the number map. Boss (2026-10-07) asked for both to be rebuilt
with every TODO.md change, so a change that forgets the rebuild fails here
with the command to run.
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
    assert sum(t["open"] for t in progress["tiers"]) == len(data["items"])
    for t in progress["tiers"]:
        assert t["total"] == t["done"] + t["open"]
        assert t["pct"] == (round(100 * t["done"] / t["total"], 1) if t["total"] else 0.0)
    assert progress["tiers"][0]["key"] == "0" and progress["tiers"][0]["open"] == 0


def test_a_finished_item_without_a_ledger_row_fails_the_build():
    items = {"GEN.1": {"status": "done", "statusText": "done, PR #9999", "phase": ""}}
    try:
        build_todo_docs.build_progress(items, {}, 1000, [])
    except SystemExit as err:
        assert "GEN.1" in str(err) and "progress-ledger.md" in str(err)
    else:
        raise AssertionError("a retired item with no ledger row must fail the build")
    done = build_todo_docs.build_progress(items, {"GEN.1": "2"}, 1000, [])
    assert done["tiers"][2]["done"] == 1
