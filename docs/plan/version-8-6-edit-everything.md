# Version 8.6: Edit everything

Rebuilt on 2026-10-16 as a version map (Boss, 1:07 PM EDT: one version per headline, phases inside a version, blockers one phase earlier).
[docs/TODO.md](../TODO.md) is the master file: it holds each item's full text, and its "Plan: versions" section indexes every version.
This file gives the version's goal, its phases and each phase's build threads with every item's prerequisites in order. A phase holds independent threads: an item's blockers sit in the same thread (earlier in the list) or in an earlier phase. Where the two disagree, TODO.md wins; when an item ships, it leaves TODO.md and its row here is deleted in the same PR (the progress ledger keeps its version). Research notes and judgment calls are in [notes.md](notes.md); the headlines and the ideas not yet filed are in [release-roadmap.md](release-roadmap.md).

Version: 8.6

## Goal

The admin has complete control: move stars, change orbits, change densities, regenerate any star to any specification, and the world adjusts after the admin approves the generation plan (ADM.57). In every regen screen the admin picks the aspects to keep and the aspects to regenerate. Every change is a plan the admin sees and approves before anything is written; the design is [edit-everything.md](../design/edit-everything.md).


## Phase 1: Plan model, aspects and regions

The foundations every editor stands on: the change plan and its store, the unit finalizer, apply and undo, the aspect catalogue, the keep-or-regenerate resolution and the recipe schema, regions, and the admin editors that already exist. The threads are independent of each other, so they can be built side by side. The unit finalizer and the recipe schema are shared with version 8.9 (generate anywhere), so they are built once here.


### Admin shell and audit

| ID | Item | Needs | Note |
|---|---|---|---|
| ADM.55 | Admin usability audit: walk every admin page and task, rank what is hard to use, and propose one layout | UX.115 | Queued admin chain (Boss, 12:52 AM EDT); first of the chain |
| ADM.92 | Admin shell: one layout, one navigation and the shared table, form, confirmation and job-log components | ADM.55 | Split from ADM.56. |

### Existing admin editors

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.120 | Bright-star backfill from the Galaxy Map's block, slab and wedge menus |  |  |

### Change plans

| ID | Item | Needs | Note |
|---|---|---|---|
| DB.36 | Tables for change plans, their operations and the edit journal (migration only) |  | Core 1 owns migrations; schema is not a constraint |
| ADM.58 | Change plan engine: operations, the effect diff and the plan states, as a library |  | Builds on planetgen/admin/edits.py |
| ADM.59 | Change plan store: save, load and list plans, and mark them stale | DB.36, ADM.58 |  |
| ADM.82 | Unit finalizer: validate a system or sector, write it in one transaction, correct only a clashing name and derived columns, or reject with reasons |  | Shared with 8.9 (API.8) |
| ADM.60 | Apply, undo and audit a plan: one job, one activity-log entry per operation | ADM.59, ADM.82 | Uses the update guard (OPS.41, OPS.53, OPS.54) |
| ADM.61 | Edits survive repair and regeneration: replay the journal after a sector is rebuilt from its seed or from parity | ADM.60, DB.17, DB.9, GEN.57 | Boss dropped GEN.59 (daily merge into a JSON file); this keeps the journal in the database |

### Aspects

| ID | Item | Needs | Note |
|---|---|---|---|
| ADM.62 | Aspect catalogue for systems, stars, planets, moons and belts: every stored field and every generation input, with how they depend on each other |  |  |
| ADM.63 | Aspect catalogue for phenomena: black holes, neutron stars, nebulae, remnants, rogue planets, comets, asteroid fields and quasars | ADM.62 |  |
| ADM.64 | Keep-or-regenerate resolution: from the admin's marks to the exact set of aspects redrawn | ADM.62 |  |
| API.28 | Recipe schema and runner: a recipe fixes, ranges or randomises any aspect, and a library turns it into objects | ADM.62, ADM.64 | Shared with 8.9 (API.18, API.19) |
| DB.37 | Lock marks table (migration only) |  | Core 1 owns migrations |

### Regions

| ID | Item | Needs | Note |
|---|---|---|---|
| DB.38 | Region density overrides table (migration only) |  | Core 1 owns migrations |
| GEN.302 | A region object: a sphere, box, shell, band, wedge or set of sectors, and the systems and phenomena inside it |  | Reads by the cell key built in 8.1 |
| MAP.212 | A region selection tool on the Galaxy Map and the Sector Map, for admin actions | GEN.302, MAP.217 |  |
| DB.60 | Migration: record the density multiplier used on each filled sector (migration only) | DB.38 | Split from MAP.119. Core 1 owns migrations. |

## Phase 2: Interfaces: approval, locks, regen screens, density overrides

The screens and marks the editors use: the approval screen, the plan history, locks, the keep-or-regenerate matrix and specification form, and density overrides for any region.


### Approval and history

| ID | Item | Needs | Note |
|---|---|---|---|
| ADM.66 | Plan approval screen: every operation with its before and after, approve, reject one, or cancel | ADM.60 |  |
| ADM.52 | A star system editor: plan a system as a 2D diagram or from the command line | ADM.92, ADM.58, ADM.66, ADM.82, API.28 |  |
| ADM.67 | Plan diagram: the system before and after as the 2D diagram, and the sector with moved stars marked | ADM.52, ADM.66, MAP.217 |  |
| ADM.68 | Plan history: every plan with who asked, approved, applied and undid it, and an undo button | ADM.66, ADM.60 |  |

### Locks and regen screens

| ID | Item | Needs | Note |
|---|---|---|---|
| ADM.65 | Locks: keep marks that stay on an object or one of its aspects, shown on the pages and honoured by every plan | DB.37, ADM.62, ADM.58, MAP.217 |  |
| ADM.69 | Keep-or-regenerate matrix: the one screen component every regen screen uses | ADM.64, ADM.65 |  |
| ADM.70 | Specification form: "regenerate to this specification" built from the aspect catalogue and the system generator options | ADM.62, ADM.69, API.28 | Reuses the system file of the system command |

### Density overrides

| ID | Item | Needs | Note |
|---|---|---|---|
| ADM.94 | Density override store and the fill read: set, undo and combine overrides for any region, and have later fills use them | DB.38, DB.60, GEN.302 | Split from MAP.119. |
| MAP.218 | Galaxy Map density control: a readout for everyone, an editor for admins over any region, with preview, confirm and undo | ADM.94, DB.38, GEN.302, MAP.212 | Split from MAP.119. |
| MAP.119 | Expected star density editable by admins on the Galaxy Map | DB.60, ADM.94, MAP.218 |  |

### Admin shell

| ID | Item | Needs | Note |
|---|---|---|---|
| ADM.93 | Move every existing admin page onto the shell and fix the defects the audit ranks above its line | ADM.92 | Split from ADM.56. |
| ADM.56 | Streamline and fix the admin pages from the audit: one admin shell, one navigation and one pattern for tables, forms, confirmations and job logs | ADM.92, ADM.93 | Queued admin chain; the shell the rest of the chain uses |
| TEST.169 | Admin page browser tests: every admin page and task works end to end after the streamline | ADM.56, TEST.168 | Quality; queued admin chain |

### Help framework

| ID | Item | Needs | Note |
|---|---|---|---|
| DOC.57 | Help section framework: page template, index, '?' links, help search and the link check |  | Split from DOC.6. |

### Admin pages in the shell

| ID | Item | Needs | Note |
|---|---|---|---|
| POP.13 | A species editor: create and edit species and their environment limits | ADM.92 |  |

## Phase 3: Editors that make plans

Each editor produces a plan: move a system, edit an orbit, adjust the world to a density change, and the rearranger that fixes a system around a change.


### Move and orbit

| ID | Item | Needs | Note |
|---|---|---|---|
| ADM.32 | Add a star system to a sector: at the emptiest spot, at given coordinates, or at random outside every Hill sphere | ADM.92, ADM.58, ADM.65, ADM.66, ADM.82 | Merges two asks (one system added by the computer; placement chosen by admin). |
| ADM.71 | Move a star system or phenomenon within its sector or to another sector, as a plan | ADM.32, ADM.58, ADM.65 |  |
| ADM.36 | Change an object's trajectory vector | GEN.109, ADM.92, ADM.58, ADM.65, ADM.66 |  |
| ADM.72 | Edit a body's orbit: distance, eccentricity, inclination and angles, with its moons following, as a plan | ADM.36, ADM.52, ADM.58, ADM.65, GEN.289 |  |

### Density plans

| ID | Item | Needs | Note |
|---|---|---|---|
| ADM.74 | Adjust the world to a density change: a plan that adds or removes systems so a region matches its new density | ADM.58, ADM.65, MAP.218 |  |

### The rearranger

| ID | Item | Needs | Note |
|---|---|---|---|
| ADM.75 | Rearrange a system around a change to a planet, moon or belt: the plan that fewest changes make valid again | ADM.58, ADM.64, ADM.65 | Boss's own example: change an A to an M planet |
| ADM.76 | Rearrange a system around a change to a star: class, mass, luminosity or companions | ADM.75 |  |

## Phase 4: Regeneration and editing any field

Regeneration to a specification for every kind of object and region, with the keep set, and the generic field editor.


### Regenerate objects

| ID | Item | Needs | Note |
|---|---|---|---|
| ADM.77 | Regenerate a planet, moon or belt to a specification, with a keep set, rearranging the system if it must | ADM.70, ADM.69, ADM.75, ADM.66 |  |
| ADM.78 | Regenerate a star or a whole system to a specification, with a keep set | ADM.70, ADM.69, ADM.76, ADM.66 |  |
| ADM.79 | Regenerate a phenomenon to a specification, with a keep set | ADM.63, ADM.70, ADM.69, ADM.66 |  |

### Regenerate regions

| ID | Item | Needs | Note |
|---|---|---|---|
| ADM.80 | Regenerate a sector or a region keeping chosen aspects: locked objects, positions, names | ADM.69, MAP.212, ADM.74, ADM.66, ADM.65 |  |

### Edit any field

| ID | Item | Needs | Note |
|---|---|---|---|
| ADM.81 | Edit any stored field of any object as a plan: one generic editor driven by the aspect catalogue | ADM.63, ADM.66, ADM.69 |  |

## Phase 5: Close-out: proof, help and the headline

End-to-end and property tests, the admin help, and the items that close ADM.35 and ADM.57.


### Proof and help

| ID | Item | Needs | Note |
|---|---|---|---|
| TEST.175 | Edit everything end to end: the headline scenarios with approval and undo | ADM.77, ADM.78, ADM.71, ADM.72, ADM.74, ADM.80, TEST.168 |  |
| TEST.176 | Plan engine properties: any sequence of edits leaves every system valid, and undo restores every row | ADM.76, ADM.60, ADM.77 |  |
| TEST.177 | Edited sectors survive repair, regeneration and the golden seed | ADM.61, ADM.80 |  |
| DOC.54 | Admin help: how to edit, regenerate, lock, approve and undo | TEST.175, ADM.68, DOC.57 |  |
| ADM.35 | Full control from every screen: edit anything, regenerate with every input, backfill or erase what is in view | MAP.120, ADM.77, ADM.78, ADM.79, ADM.80, ADM.81 | Closes when the regen and edit screens of phase 4 are done; MAP.120 is the map part. |
| ADM.57 | Headline of 8.6: edit everything, the admin can change anything and the world adjusts after the admin approves a generation plan | ADM.32, ADM.35, ADM.36, ADM.52, MAP.119, MAP.120, POP.13, ADM.68, ADM.71, ADM.72, ADM.74, ADM.77, ADM.78, ADM.79, ADM.80, ADM.81, TEST.175, TEST.176, TEST.177, DOC.54 | The headline; done when its items are |
