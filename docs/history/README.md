# Project history

What was chosen, what was dropped and how the work moved, read from the repository as of commit `d02661568` (2026-11-08): 12144 commits, 3399 releases, 3349 items named in commit titles, 2034 dated decisions in the documents. Built by `scripts/build_history.py`; this folder is a snapshot (git history is not in a shallow clone, so no test compares it with its sources). Rebuild it in a full clone when it should catch up, and commit the result.

## Where to look

- [timeline.md](timeline.md): commits, releases, items filed and retired per week; work by category.
- [dropped-and-closed.md](dropped-and-closed.md): paths that were retired, dropped, closed as findings, reverted or replaced.
- [by-note.md](by-note.md): each design note with the dates and state of the items behind it.
- Release headlines by month (the best feature of each release): [2024-06](releases-2024-06.md), [2024-07](releases-2024-07.md), [2026-06](releases-2026-06.md), [2026-07](releases-2026-07.md), [2026-09](releases-2026-09.md), [2026-10](releases-2026-10.md), [2026-11](releases-2026-11.md).
- Decisions by month (Boss's choices as the documents record them): [2026-09](decisions-2026-09.md), [2026-10](decisions-2026-10.md), [2026-11](decisions-2026-11.md).
- Item histories: [UX](items-UX.md), [MAP](items-MAP.md), [NAV](items-NAV.md), [GEN](items-GEN.md), [PERF](items-PERF.md), [DB](items-DB.md), [API](items-API.md), [ADM](items-ADM.md), [SEC](items-SEC.md), [TEST](items-TEST.md), [USR](items-USR.md), [OPS](items-OPS.md), [DOC](items-DOC.md), [VIEW](items-VIEW.md), [POP](items-POP.md).

## How the design evolved, in four eras

Read from the timeline, the release headlines and the decision log; each era's decisions are in the month pages above.

- **Before 2026-09**: a first generator (changelog 0.1 to 7.x, June 2024), then a pause.
- **2026-09 to 2026-10-09 (Phase 0 and 1)**: the package layout, MySQL, Alembic, per-unit seeds, the map engine and Redis queues; release numbers 7.x. Phase 0 closed on 2026-10-09 and the version became 8.0.
- **2026-10-10 to 2026-10-31**: the version map (8.1 to 10.0, 2026-10-16), the lane model (cores, bugs, aux), the TODO thread as sole owner of items and global documents, the stand-in and source rules, condition-based planet classes and the nine life classes.
- **2026-11 onward (8.1 build-out)**: merge trains (train 1 on 2026-10-29 to train 19 on 2026-11-08), the gate, lanes renumbered on 2026-11-07, the everything-affects-everything principle (2026-11-08) and real materials.

