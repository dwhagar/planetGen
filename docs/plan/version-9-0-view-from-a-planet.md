# Version 9.0: The view from a planet

Rebuilt on 2026-10-16 as a version map (Boss, 1:07 PM EDT: one version per headline, phases inside a version, blockers one phase earlier).
[docs/TODO.md](../TODO.md) is the master file: it holds each item's full text, and its "Plan: versions" section indexes every version.
This file gives the version's goal, its phases and each phase's build threads with every item's prerequisites in order. A phase holds independent threads: an item's blockers sit in the same thread (earlier in the list) or in an earlier phase. Where the two disagree, TODO.md wins; when an item ships, it leaves TODO.md and its row here is deleted in the same PR (the progress ledger keeps its version). Research notes and judgment calls are in [notes.md](notes.md); the headlines and the ideas not yet filed are in [release-roadmap.md](release-roadmap.md).

Version: 9.0

## Goal

Stand on a world and look up: the sky from any planet or moon. Proposed headline of 9.0 (recommended by the TODO thread; Boss has not chosen). A major version needs Boss's word.

## Phase 1: Sky data

The sky data and the research items that settle it.

### Sky data

| ID | Item | Needs | Note |
|---|---|---|---|
| VIEW.4 | Constellation names in the name generator |  | Constellation names from the codec under the naming key, not sliced word lists. Floats: no blockers, can run any time. |
| VIEW.6 | Settle the handedness of the generated galaxy before mapping real sky coordinates |  | Research: needed before real sky coordinates. |
| VIEW.9 | Replace the BC_V table with a published relation (Flower 1996 or Torres 2010) |  | Research. |
| VIEW.10 | Declare `pillow` in `setup.py` if it is used, and add a HYG calibration test |  | Research. |

## Phase 2: Sky rendering

The view itself, the star map and the PNG.

### Sky rendering

| ID | Item | Needs | Note |
|---|---|---|---|
| VIEW.1 | View from a planet |  | Research session with Boss. |
| VIEW.2 | A starmap seen from a planet. RESEARCH WITH BOSS FIRST | VIEW.1 |  |
| VIEW.3 | Render the view as a PNG, with constellations | VIEW.2, VIEW.4 |  |
| VIEW.7 | Sky tiers and the backfill cost measurement | VIEW.1 | Research: the sky is about 8% to 10% complete today. |
| VIEW.8 | A dust model for the sky: the galactic extinction law and `nebulae.extinction_av` | VIEW.1 | Research. |
