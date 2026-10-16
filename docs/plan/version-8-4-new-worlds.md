# Version 8.4: New worlds

Rebuilt on 2026-10-16 as a version map (Boss, 1:07 PM EDT: one version per headline, phases inside a version, blockers one phase earlier).
[docs/TODO.md](../TODO.md) is the master file: it holds each item's full text, and its "Plan: versions" section indexes every version.
This file gives the version's goal, its phases and each phase's build threads with every item's prerequisites in order. A phase holds independent threads: an item's blockers sit in the same thread (earlier in the list) or in an earlier phase. Where the two disagree, TODO.md wins; when an item ships, it leaves TODO.md and its row here is deleted in the same PR (the progress ledger keeps its version). Research notes and judgment calls are in [notes.md](notes.md); the headlines and the ideas not yet filed are in [release-roadmap.md](release-roadmap.md).

Version: 8.4

## Goal

The planet class refactor around the habitability index (classes R to Z), nebula conditions on planets, asteroid fields as object systems and the habitability symbols on the pages. Proposed by the TODO thread; not yet approved by Boss.

## Phase 1: New worlds

Planet classes, nebula conditions, asteroid fields and habitability on the pages.

### Planet classes

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.33 | One class per PR, each with its tests |  | Moved to phase 2 under GEN.90, after the habitability score. One class per PR (R and S first). Built on fixed physics so new classes aren't tuned to wrong masses, moons or zones. |
| GEN.27 | Class P (glaciated world) only in the habitable zone, and fitting there | GEN.33 | Under GEN.90. Same reconcile/zone code as phase 0's physics fixes. |
| GEN.28 | Seven new planet classes in the letter gaps (R, S, U, W, X, Y, Z) | GEN.33 | Under GEN.90; class Z is Boss's Earth-size world that never had life (2026-10-03). Class S landed with GEN.38 (PR #415); the other six classes here. Closes with GEN.33's PRs; PLANET_CLASSES in program_constants.py. |
| GEN.91 | Classes like S and V in the hot and cold zones | GEN.33 |  |
| GEN.92 | Life and its highest stage follow the habitability score | GEN.28 |  |
| GEN.29 | Sweep every planet class for sense once the new ones are in (bug) | GEN.28, GEN.27, GEN.91, GEN.92 | Under GEN.90: the refactor is the sweep. Bug, but by definition a sweep after the new classes; it can't go earlier. Includes rocky rogues of 10-16 Earth masses (up to 17,600 km) that get S as nearest fit (PR #415). |
| GEN.90 | Refactor the planet classes around the habitability index | GEN.33, GEN.28, GEN.27, GEN.29, GEN.91, GEN.92 | Parent; takes GEN.33, GEN.28, GEN.27 and GEN.29 as its subitems. |
| GEN.145 | Class S atmosphere rule: S keeps air unless the shoreline ratio is over 30 | GEN.91 | Research follow-up to class S (built). |
| GEN.146 | Teff-dependent habitable zone from the Kopparapu table | GEN.91 | Research: in GEN.91's dependency chain. |

### Nebula conditions on planets

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.95 | Nebula conditions applied when planets and surfaces are generated |  |  |
| GEN.93 | Nebula conditions in planet generation | GEN.95 | Parent. |

### Habitability on the pages

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.143 | Color sectors by their number of habitable locations |  | GitHub issue [#717](https://github.com/dwhagar/planetGen/issues/717). |
| UX.32 | Planet rows show the class only, without the type and moon labels |  | Same rows as UX.105 (absorbed UX.30); one thread. |
| UX.100 | Planet and moon lists show habitability as symbols, not notes |  |  |
| UX.101 | PHI-6 symbols in the planet and moon lists, coloured by rating | UX.100 |  |

### Asteroid fields

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.112 | Plan asteroid fields and belts as object systems for rendering |  | Plan only. |
