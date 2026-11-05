# Version 9.2: Politics and territories

Rebuilt on 2026-10-16 as a version map (Boss, 1:07 PM EDT: one version per headline, phases inside a version, blockers one phase earlier).
[docs/TODO.md](../../TODO.md) is the master file: it holds each item's full text, and its "Plan: versions" section indexes every version.
This file gives the version's goal, its phases and each phase's build threads with every item's prerequisites in order. A phase holds independent threads: an item's blockers sit in the same thread (earlier in the list) or in an earlier phase. Where the two disagree, TODO.md wins; when an item ships, it leaves TODO.md and its row here is deleted in the same PR (the progress ledger keeps its version). Research notes and judgment calls are in [notes.md](../notes.md); the headlines and the ideas not yet filed are in [release-roadmap.md](../release-roadmap.md).

Version: 9.2

Status: paused until 8.1 is complete (Boss, 7:44 PM EDT 2026-10-22: the five cores work 8.1 only; nothing in this version is started or owned until then, except the gates already running for UX.140, UX.128 and NAV.17).

## Goal

Territories that grow and shift over time on the map, with borders drawn and polities that trade. Boss (2026-10-16): politics and territories come before ships and stations. No items exist yet beyond the population design note.

Release gate (standing, Boss, 2026-10-21, 5:10 PM EDT): this version is not stamped until SEC.33 (the full security code audit) and TEST.261 (the full bug hunt audit) each have a dated report for it and no finding they mark blocking is open. Both are standing items: they stay open and are re-run for every major or minor version change (required whenever the major or minor version changes (8.1, 8.2 and on to 9.0 and 10.0), never for build numbers).

## Phase 1: Dynamic settlements

The fixed settlements of the 8.x chain become dynamic (Boss, 9:27 PM EDT 2026-10-24). The item below holds the place; the breakdown follows Boss's answers.

### Dynamic settlements (Boss approved the 9.2 design at 6:00 PM EST on 2026-11-05)

| ID | Item | Needs | Note |
|---|---|---|---|
| POP.35 | Dynamic world populations: populations step forward with the nightly update from the stored growth state | POP.30 | Boss, 9:38 PM EDT 2026-10-24: trade routes and world populations, done by 8.9, dynamic in 9.1. |
| POP.48 | Claim strength and frontier ownership in the political step: strength is population times a tech factor times cohesion, a system goes to the polity with the largest strength over distance, and it changes hands only when a challenger beats the holder by the hysteresis margin | GEN.1293, POP.46, DB.110, POP.35, DB.104 |  |
| POP.49 | Contested systems and the derived war state: a system whose two best claims are within the contested ratio is marked contested, and two polities with a contested frontier are at war, which closes trade between them and changes nothing else | POP.48 |  |
| POP.36 | Dynamic trade routes: routes grow, shift and close as populations, ports and territories change, using the same hop pathfinding | POP.33, POP.35, MAP.309, POP.49 | Boss, 9:38 PM EDT 2026-10-24: trade routes and world populations, done by 8.9, dynamic in 9.1. |
| POP.50 | Cohesion and splits: a polity's cohesion rises on contested frontier systems and falls with the share of territory far from its capital, and below the split threshold the polity splits along its weakest line into polities of the same species | POP.49 |  |
| POP.51 | Planet level in the political step: nations drift in dominance by a slow power score, cities change nation by the same weighted draw as at generation re-run with current dominance (the dominant nation keeps its capital), city sizes grow by a proportional law with a spread capped by carrying capacity, and ports follow their city's owner | POP.48, POP.35 |  |
| ADM.138 | The political step skips locked polities, nations and cities and shows the lock in the editors; releasing a lock lets the step move it again | DB.110, POP.48, ADM.131 |  |
| MAP.390 | Past borders: the Territories layer gets a time control that shows the borders at any stored epoch, from the history rows the step writes | MAP.389, POP.48 |  |
| UX.271 | The polity history page: a timeline of its territory, cohesion, splits and trade totals, and the journal feed on the API | POP.50, POP.51 |  |
| TEST.624 | Politics tests: determinism over 100 steps at different batchings, no flicker, a split, a locked case, and polity count and area in a sane range on a census galaxy | POP.48, POP.49, POP.50, POP.51, ADM.138 |  |
| DOC.105 | Docs for politics: the help page for borders, contested systems, splits and the map time control, and the design note | MAP.390, UX.271 |  |
| POP.24 | Settlements become dynamic: nations and cities grow, borders and dominance shift and cities change hands over time | ADM.131, TEST.381, DB.104, POP.26, POP.27, POP.48, POP.49, POP.50, POP.51, ADM.138, MAP.390, UX.271, TEST.624, DOC.105 | Boss, 9:27 PM EDT 2026-10-24: settlements tree; fixed until 9.2. |
