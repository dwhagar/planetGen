# Version 8.6: Edit everything

Rebuilt on 2026-10-16 as a version map (Boss, 1:07 PM EDT: one version per headline, phases inside a version, blockers one phase earlier).
[docs/TODO.md](../TODO.md) is the master file: it holds each item's full text, and its "Plan: versions" section indexes every version.
This file gives the version's goal, its phases and each phase's build threads with every item's prerequisites in order. A phase holds independent threads: an item's blockers sit in the same thread (earlier in the list) or in an earlier phase. Where the two disagree, TODO.md wins; when an item ships, it leaves TODO.md and its row here is deleted in the same PR (the progress ledger keeps its version). Research notes and judgment calls are in [notes.md](notes.md); the headlines and the ideas not yet filed are in [release-roadmap.md](release-roadmap.md).

Version: 8.6

## Goal

The admin has complete control: move stars, change orbits, change densities, regenerate any star to any specification, and the world adjusts after the admin approves the generation plan (ADM.57). In every regen screen the admin picks the aspects to keep and the aspects to regenerate. The item breakdown is filed in the next change; this file holds the items that already exist.

## Phase 1: The editors that exist

The admin shell and the admin editors already on the list.

### Admin shell and audit

| ID | Item | Needs | Note |
|---|---|---|---|
| ADM.55 | Admin usability audit: walk every admin page and task, rank what is hard to use, and propose one layout | UX.115 | Queued admin chain (Boss, 12:52 AM EDT); first of the chain |
| ADM.56 | Streamline and fix the admin pages from the audit: one admin shell, one navigation and one pattern for tables, forms, confirmations and job logs | ADM.55, ADM.50 | Queued admin chain; the shell the rest of the chain uses |
| TEST.169 | Admin page browser tests: every admin page and task works end to end after the streamline | ADM.56 | Quality; queued admin chain |

### Existing admin editors

| ID | Item | Needs | Note |
|---|---|---|---|
| ADM.32 | Add a star system to a sector: at the emptiest spot, at given coordinates, or at random outside every Hill sphere |  | Merges two asks (one system added by the computer; placement chosen by admin). |
| ADM.36 | Change an object's trajectory vector | GEN.109 |  |
| ADM.52 | A star system editor: plan a system as a 2D diagram or from the command line |  |  |
| MAP.119 | Expected star density editable by admins on the Galaxy Map |  |  |
| MAP.120 | Bright-star backfill from the Galaxy Map's block, slab and wedge menus |  |  |
| ADM.35 | Full control from every screen: edit anything, regenerate with every input, backfill or erase what is in view | MAP.120 | Parent; MAP item MAP.120 is the map part. |
| POP.13 | A species editor: create and edit species and their environment limits |  |  |

## Phase 2: The headline

The umbrella that closes when the edit-everything items are done.

### Edit everything: the headline

| ID | Item | Needs | Note |
|---|---|---|---|
| ADM.57 | Headline of 8.6: edit everything, the admin can change anything and the world adjusts after the admin approves a generation plan | ADM.32, ADM.35, ADM.36, ADM.52, MAP.119, MAP.120, POP.13 |  |
