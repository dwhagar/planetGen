# Version 8.7: The pilot's console and built-in help

Rebuilt on 2026-10-16 as a version map (Boss, 1:07 PM EDT: one version per headline, phases inside a version, blockers one phase earlier).
[docs/TODO.md](../TODO.md) is the master file: it holds each item's full text, and its "Plan: versions" section indexes every version.
This file gives the version's goal, its phases and each phase's build threads with every item's prerequisites in order. A phase holds independent threads: an item's blockers sit in the same thread (earlier in the list) or in an earlier phase. Where the two disagree, TODO.md wins; when an item ships, it leaves TODO.md and its row here is deleted in the same PR (the progress ledger keeps its version). Research notes and judgment calls are in [notes.md](notes.md); the headlines and the ideas not yet filed are in [release-roadmap.md](release-roadmap.md).

Version: 8.7

## Goal

A built-in help section and the pilot's starmap look. Proposed by the TODO thread; not yet approved by Boss.

## Phase 1: Help and console look

The help pages and the console design.

### Help section

| ID | Item | Needs | Note |
|---|---|---|---|
| DOC.6 | A static help section in the web interface: page template, index, per-page help links and a coverage test |  | Boss 23:53Z; the framework the per-feature help pages use. |
| DOC.7 | The Galaxy Map help page: layers, zoom, fly-through, Color by, select modes, bookmarks and the locate box | DOC.6 | Boss 23:53Z; Phase 2. |
| DOC.8 | The sector help pages: the sector list, a sector page, the sector map and the sector scene | DOC.6 | Boss 23:53Z; Phase 2. |
| DOC.9 | The star system help pages: the system list, a system page, the system map and the planets, moons and belts shown | DOC.6 | Boss 23:53Z; Phase 2. |
| DOC.10 | The search and navigation help pages: search, nearby, the nav page and routes | DOC.6 | Boss 23:53Z; Phase 2. |
| DOC.11 | The reference browser help pages: species, polities, object classes and phenomena | DOC.6 | Boss 23:53Z; Phase 2. |
| DOC.12 | The account help pages: signing in, two-factor, the account page, API keys and bookmarks | DOC.6 | Boss 23:53Z; Phase 2. |
| DOC.13 | The Generate page help pages: layer specs, spans, radial fills, directives, one-off systems and jobs | DOC.6 | Boss 23:53Z; Phase 2. |
| DOC.14 | The admin help pages: the queue, the stats page, settings, lockouts and the naming key | DOC.6 | Boss 23:53Z; Phase 2. |
| DOC.15 | A glossary and units help page: coordinates, scales, sector paths, object IDs, time and the in-universe wording | DOC.6 | Boss 23:53Z; Phase 2. |
| DOC.16 | The API help page for visitors: what the API is, how to get a key, and where the reference lives | DOC.6 | Boss 23:53Z; Phase 2. |

### Console look

| ID | Item | Needs | Note |
|---|---|---|---|
| UX.43 | A visual design built like a pilot's starmap and navigation console |  |  |
| UX.82 | Theme checks after PR #800: SVG currentColor, two Shoelace contrast failures, alpha in --bg-subtle |  | Research follow-up to UX.76 (built). |
