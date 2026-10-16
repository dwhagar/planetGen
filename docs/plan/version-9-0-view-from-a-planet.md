# Version 9.0: The view from a planet

Rebuilt on 2026-10-16 as a version map (Boss, 1:07 PM EDT: one version per headline, phases inside a version, blockers one phase earlier).
[docs/TODO.md](../TODO.md) is the master file: it holds each item's full text, and its "Plan: versions" section indexes every version.
This file gives the version's goal, its phases and each phase's build threads with every item's prerequisites in order. A phase holds independent threads: an item's blockers sit in the same thread (earlier in the list) or in an earlier phase. Where the two disagree, TODO.md wins; when an item ships, it leaves TODO.md and its row here is deleted in the same PR (the progress ledger keeps its version). Research notes and judgment calls are in [notes.md](notes.md); the headlines and the ideas not yet filed are in [release-roadmap.md](release-roadmap.md).

Version: 9.0

## Goal

Stand on a world and look up: the sky from any planet or moon. Proposed headline of 9.0 (recommended by the TODO thread; Boss has not chosen). A major version needs Boss's word. The version also carries the backups and repair chain Boss targeted for 9.0 (Boss, 2026-10-16, 2:00 PM EDT; ADM.91): phases 3 and 4.

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

## Phase 3: Backup parts and restore

Target 9.0 (Boss, 2:00 PM EDT): the control database and parity parts of the backup archive, and a restore for each.

### Control database

| ID | Item | Needs | Note |
|---|---|---|---|
| ADM.88 | Control database backup part: a dump of the planetgen_control schema (admin logins, API keys, run records) in the backup archive | ADM.87 | Boss 1:59 PM EDT; part of the ADM.87 archive. |
| ADM.89 | Control database restore: load a backed-up planetgen_control into a scratch schema, check it, then replace the live schema, from the command and from the settings page | ADM.88, ADM.85, OPS.41 | Boss 2:00 PM EDT. |

### Parity data

| ID | Item | Needs | Note |
|---|---|---|---|
| DB.40 | Parity data backup part: the parity file of the DB repair copied into the backup archive | ADM.87, DB.9 | Boss 1:59 PM EDT; Core 1 owns DB.9. |
| DB.41 | Parity data restore: put a backed-up parity file back and check it against the database | DB.40, DB.8 | Boss 2:00 PM EDT; Core 1. |

## Phase 4: Repair from backups, the repair page and proof

The repair of the galaxy from an archive, the one page for check, repair and restore, and the end-to-end test, guide and umbrella.

### Repair and proof

| ID | Item | Needs | Note |
|---|---|---|---|
| DB.42 | Repair the galaxy from an archive: `repair-db --from-archive` uses the parity part of a backup when the live parity file is missing or damaged | DB.41, DB.9, DB.17 | Boss 2:00 PM EDT. |
| ADM.90 | Database health and repair page: run the check, see damaged and unprotected sectors, repair from the live parity or from an archive, restore parts, as a queue job | DB.42, ADM.89, DB.41 | Boss 2:00 PM EDT. |
| TEST.180 | Backup and repair end to end: take a full archive, lose the control schema and the live parity file, restore and repair, and prove the galaxy and the logins come back | ADM.89, DB.42, ADM.90 | Boss 2:00 PM EDT. |
| DOC.56 | Backups, restore and repair: the deployment guide and the help page | TEST.180 | Boss 2:00 PM EDT. |
| ADM.91 | Backups and repair, target 9.0: one bz2 archive of the settings, the control database and the parity data, and the ways to use it (restore the control database, restore parity, repair the galaxy, the repair page) | ADM.88, ADM.89, DB.40, DB.41, DB.42, ADM.90, TEST.180, DOC.56 | Umbrella of the backups and repair chain; done when its items are. |
