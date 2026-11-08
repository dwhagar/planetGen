# The plan: where each kind of fact lives

This page is the map of everything under `docs/plan/` and its neighbours (DOC.115). It names the **one source** for each kind of fact, so a fact is written once and every other page points at it. When two pages disagree, the source named here wins. Nothing was merged away to make this page: every plan file keeps its item numbers, its links and its history; the rule below is only about which one to trust and edit.

## One source for each fact

| Fact | Its one source | Written by | Pages built from it |
|---|---|---|---|
| What work is open: the item, what is wrong, where to start, what done means, the default taken | [../TODO.md](../TODO.md) (one block per ID) | TODO thread (lane 01) only | [todo-tree](../generated/todo-tree.html), [todo-reference](../generated/todo-reference.html), [version-plan](../generated/version-plan.html) |
| The next free ID of each category and the old-number to ID map | [todo-number-map.md](todo-number-map.md) | TODO thread | [todo-reference](../generated/todo-reference.html) |
| Which version an open item belongs to, and in what order | the file for that version in [releases/](releases/) | TODO thread | [version-plan](../generated/version-plan.html) |
| What a version is for: headline, phases, needs, share finished | [release-roadmap.md](release-roadmap.md) | TODO thread | [release-roadmap](../generated/release-roadmap.html) |
| The version of an item that has shipped | [progress-ledger.md](progress-ledger.md) | TODO thread | progress bars on the roadmap |
| Which lane (thread) holds which item, and each lane's state | [lanes.md](lanes.md); the lane numbers and old names are in [lane-designations.md](lane-designations.md) | TODO thread, from the lanes' own status lines | [execution-plan](../generated/execution-plan.html) |
| What two lanes might both break, and what needs Boss | [execution-notes.md](execution-notes.md) | TODO thread | [execution-plan](../generated/execution-plan.html) |
| Which GitHub issue (goal) groups which items | [issue-groups.md](issue-groups.md) | TODO thread | every page that links items to goals |
| The plan for space beyond the galaxy | [intergalactic-space-plan.md](intergalactic-space-plan.md) | TODO thread | none (read as written) |
| Old phase notes and dependency notes, kept as written | [notes.md](notes.md) (history; TODO.md wins where they disagree) | no one (frozen) | none |
| How a part of the program is designed | the note in [../design/](../design/README.md), named by its row in the design index | the lane that builds it | none |
| What was chosen, dropped or replaced, and when | [../history/](../history/README.md) (snapshot) and the design notes' own status lines | `scripts/build_history.py` | none |
| What a release changed | [../../CHANGELOG.md](../../CHANGELOG.md) (from `changes/*.md`) | the release stamp | none |
| What the code does today | the code, then [../reference/](../reference/), then [../generated/database-columns.md](../generated/database-columns.md) and [../generated/api-routes.md](../generated/api-routes.md) | the lanes | none |

## Rules that hold for every plan page

- **Lanes never edit the TODO files.** A lane tells TODO; TODO files, changes and retires items, edits the bookkeeping files in this folder and rebuilds the generated pages (`python scripts/build_global_docs.py`).
- **IDs are permanent.** An ID is never renumbered or reused; a shipped item is deleted from TODO.md (its ID stays taken, and the progress ledger keeps its version). Never invent an ID: ask TODO.
- **Every open item sits in a version file** (or the roadmap's "later" list) so it has a place; the audit of 2026-11-08 found all 620 open items placed.
- **Everything in a system affects everything else in the system** (Boss, 2026-11-08): every body acts on every other body of its system by 3D position, averaged over the orbits; a host-only or pair-only effect needs a reason registered in GEN.725. The rule and its reasons are in [../design/generation/system-interactions.md](../design/generation/system-interactions.md) and at the top of [../design/README.md](../design/README.md).
- **All jobs run under one memory and load gate** (PERF.254, OPS.124, 2026-11-08): the work queue, the Generate page, command-line runs and the warm-up all start only when the gate says the machine has room. The design is in [../design/operations/redis-and-the-work-queue.md](../design/operations/redis-and-the-work-queue.md).
- **Every merge goes through the gate and, for generation changes, a merge train** (lane 02, lane 03); a merge tests only what changed (`scripts/test.py changed --resume`) and the full run is kept for the end-of-version bug hunt.
- **Data and designs are archived, never deleted** without Boss's say: replaced designs go to [../archive/](../archive/) or stay in place with a `history` status.
- **Times to Boss are EST clock times.** Plans written before 2026-11-01 use EDT; both are kept as written.

## Reading order for someone new

1. [release-roadmap.md](release-roadmap.md): where the project is going, version by version.
2. [lanes.md](lanes.md) and [execution-notes.md](execution-notes.md): who is building what now.
3. The version file of the version you care about in [releases/](releases/).
4. [../TODO.md](../TODO.md): the item itself.
5. The design note the item points at, then [../history/by-note.md](../history/by-note.md) for how it got there.
