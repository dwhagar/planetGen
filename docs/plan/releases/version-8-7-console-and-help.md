# Version 8.7: The pilot's console and built-in help

Rebuilt on 2026-10-16 as a version map (Boss, 1:07 PM EDT: one version per headline, phases inside a version, blockers one phase earlier).
[docs/TODO.md](../../TODO.md) is the master file: it holds each item's full text, and its "Plan: versions" section indexes every version.
This file gives the version's goal, its phases and each phase's build threads with every item's prerequisites in order. A phase holds independent threads: an item's blockers sit in the same thread (earlier in the list) or in an earlier phase. Where the two disagree, TODO.md wins; when an item ships, it leaves TODO.md and its row here is deleted in the same PR (the progress ledger keeps its version). Research notes and judgment calls are in [notes.md](../notes.md); the headlines and the ideas not yet filed are in [release-roadmap.md](../release-roadmap.md).

Version: 8.7

Status: paused until 8.1 is complete (Boss, 7:44 PM EDT 2026-10-22: the five cores work 8.1 only; nothing in this version is started or owned until then, except the gates already running for UX.140, UX.128 and NAV.17).

## Goal

A built-in help section and the pilot's starmap look. Proposed by the TODO thread; not yet approved by Boss.

## Phase 1: Help and console look

The help pages and the console design.

Release gate (standing, Boss, 2026-10-21, 5:10 PM EDT): this version is not stamped until SEC.33 (the full security code audit) and TEST.261 (the full bug hunt audit) each have a dated report for it and no finding they mark blocking is open. Both are standing items: they stay open and are re-run for every major or minor version change (required whenever the major or minor version changes (8.1, 8.2 and on to 9.0 and 10.0), never for build numbers).

### Help section

| ID | Item | Needs | Note |
|---|---|---|---|
| DOC.7 | The Galaxy Map help page: layers, zoom, fly-through, Color by, select modes, bookmarks and the locate box |  | Boss 23:53Z; Phase 2. |
| DOC.8 | The sector help pages: the sector list, a sector page, the sector map and the sector scene |  | Boss 23:53Z; Phase 2. |
| DOC.9 | The star system help pages: the system list, a system page, the system map and the planets, moons and belts shown |  | Boss 23:53Z; Phase 2. |
| DOC.10 | The search and navigation help pages: search, nearby, the nav page and routes |  | Boss 23:53Z; Phase 2. |
| DOC.11 | The reference browser help pages: species, polities, object classes and phenomena |  | Boss 23:53Z; Phase 2. |
| DOC.12 | The account help pages: signing in, two-factor, the account page, API keys and bookmarks |  | Boss 23:53Z; Phase 2. |
| DOC.13 | The Generate page help pages: layer specs, spans, radial fills, directives, one-off systems and jobs |  | Boss 23:53Z; Phase 2. |
| DOC.14 | The admin help pages: the queue, the stats page, settings, lockouts and the naming key |  | Boss 23:53Z; Phase 2. |
| DOC.15 | A glossary and units help page: coordinates, scales, sector paths, object IDs, time and the in-universe wording |  | Boss 23:53Z; Phase 2. |
| DOC.16 | The API help page for visitors: what the API is, how to get a key, and where the reference lives |  | Boss 23:53Z; Phase 2. |
| DOC.58 | Help coverage test: every user-facing route has a help page and every help page has a route | DOC.7, DOC.8, DOC.9, DOC.10, DOC.11, DOC.12, DOC.13, DOC.14, DOC.15, DOC.16, DOC.54 | Split from DOC.6. |

### Console look

| ID | Item | Needs | Note |
|---|---|---|---|
| UX.82 | Theme checks after PR #800: SVG currentColor, two Shoelace contrast failures, alpha in --bg-subtle |  | Research follow-up to UX.76 (built). |
| UX.130 | Console design brief: colours, type, panel shapes and iconography, approved by Boss |  | Split from UX.43. |
| UX.131 | Console tokens, type and base template in both themes | UX.130, UX.82 | Split from UX.43. |
| UX.132 | Console look for the shared components: buttons, tables, panels, status chips, icons and the maps' chrome | UX.131 | Split from UX.43. |
| UX.43 | A visual design built like a pilot's starmap and navigation console | UX.130, UX.131, UX.132 |  |

### Resource help

| ID | Item | Needs | Note |
|---|---|---|---|
| DOC.81 | Help page: resource locations, what the scores mean and how to find them on the map | UX.223, NAV.64 | Boss, 10:29 PM EDT 2026-10-24: resource locations (valuable belts and fields as points of interest), foundations 8.2 to 9.0. |
