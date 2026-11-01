# Version 8.6: Edit everything

Rebuilt on 2026-10-16 as a version map (Boss, 1:07 PM EDT: one version per headline, phases inside a version, blockers one phase earlier).
[docs/TODO.md](../../TODO.md) is the master file: it holds each item's full text, and its "Plan: versions" section indexes every version.
This file gives the version's goal, its phases and each phase's build threads with every item's prerequisites in order. A phase holds independent threads: an item's blockers sit in the same thread (earlier in the list) or in an earlier phase. Where the two disagree, TODO.md wins; when an item ships, it leaves TODO.md and its row here is deleted in the same PR (the progress ledger keeps its version). Research notes and judgment calls are in [notes.md](../notes.md); the headlines and the ideas not yet filed are in [release-roadmap.md](../release-roadmap.md).

Version: 8.6

Status: paused until 8.1 is complete (Boss, 7:44 PM EDT 2026-10-22: the five cores work 8.1 only; nothing in this version is started or owned until then, except the gates already running for UX.140, UX.128 and NAV.17).

## Goal

The admin has complete control: move stars, change orbits, change densities, regenerate any star to any specification, and the world adjusts after the admin approves the generation plan (ADM.57). In every regen screen the admin picks the aspects to keep and the aspects to regenerate. Every change is a plan the admin sees and approves before anything is written; the design is [edit-everything.md](../../design/interface/edit-everything.md).


## Phase 1: Plan model, aspects and regions

The foundations every editor stands on: the change plan and its store, the unit finalizer, apply and undo, the aspect catalogue, the keep-or-regenerate resolution and the recipe schema, regions, and the admin editors that already exist. The threads are independent of each other, so they can be built side by side. The unit finalizer and the recipe schema are shared with version 8.9 (generate anywhere), so they are built once here.


### Admin shell and audit

| ID | Item | Needs | Note |
|---|---|---|---|

### Existing admin editors

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.120 | Bright-star backfill from the Galaxy Map's block, slab and wedge menus |  |  |

### Change plans

| ID | Item | Needs | Note |
|---|---|---|---|
| ADM.59 | Change plan store: save, load and list plans, and mark them stale |  |  |
| ADM.82 | Unit finalizer: validate a system or sector, write it in one transaction, correct only a clashing name and derived columns, or reject with reasons |  | Shared with 8.9 (API.8) |
| ADM.60 | Apply, undo and audit a plan: one job, one activity-log entry per operation | ADM.59, ADM.82 | Uses the update guard (OPS.41, OPS.53, OPS.54) |
| ADM.61 | Edits survive repair and regeneration: replay the journal after a sector is rebuilt from its seed or from parity | ADM.60, DB.17 | Boss dropped GEN.59 (daily merge into a JSON file); this keeps the journal in the database |

### Aspects

| ID | Item | Needs | Note |
|---|---|---|---|
| ADM.62 | Aspect catalogue for systems, stars, planets, moons and belts: every stored field and every generation input, with how they depend on each other |  |  |
| ADM.63 | Aspect catalogue for phenomena: black holes, neutron stars, nebulae, remnants, rogue planets, comets, asteroid fields and quasars | ADM.62 |  |
| ADM.64 | Keep-or-regenerate resolution: from the admin's marks to the exact set of aspects redrawn | ADM.62 |  |
| API.28 | Recipe schema and runner: a recipe fixes, ranges or randomises any aspect, and a library turns it into objects | ADM.62, ADM.64 | Shared with 8.9 (API.18, API.19) |

### Regions

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.302 | A region object: a sphere, box, shell, band, wedge or set of sectors, and the systems and phenomena inside it |  | Reads by the cell key built in 8.1 |
| MAP.212 | A region selection tool on the Galaxy Map and the Sector Map, for admin actions | GEN.302 |  |

## Phase 2: Interfaces: approval, locks, regen screens, density overrides

The screens and marks the editors use: the approval screen, the plan history, locks, the keep-or-regenerate matrix and specification form, and density overrides for any region.


### Approval and history

| ID | Item | Needs | Note |
|---|---|---|---|
| ADM.66 | Plan approval screen: every operation with its before and after, approve, reject one, or cancel | ADM.60 |  |
| ADM.52 | A star system editor: plan a system as a 2D diagram or from the command line | ADM.66, ADM.82, API.28 |  |
| ADM.67 | Plan diagram: the system before and after as the 2D diagram, and the sector with moved stars marked | ADM.52, ADM.66 |  |
| ADM.68 | Plan history: every plan with who asked, approved, applied and undid it, and an undo button | ADM.66, ADM.60 |  |

### Locks and regen screens

| ID | Item | Needs | Note |
|---|---|---|---|
| ADM.65 | Locks: keep marks that stay on an object or one of its aspects, shown on the pages and honoured by every plan | ADM.62 |  |
| ADM.69 | Keep-or-regenerate matrix: the one screen component every regen screen uses | ADM.64, ADM.65 |  |
| ADM.70 | Specification form: "regenerate to this specification" built from the aspect catalogue and the system generator options | ADM.62, ADM.69, API.28 | Reuses the system file of the system command |

### Density overrides

| ID | Item | Needs | Note |
|---|---|---|---|
| ADM.94 | Density override store and the fill read: set, undo and combine overrides for any region, and have later fills use them | GEN.302 | Split from MAP.119. |
| MAP.218 | Galaxy Map density control: a readout for everyone, an editor for admins over any region, with preview, confirm and undo | ADM.94, GEN.302, MAP.212 | Split from MAP.119. |
| MAP.119 | Expected star density editable by admins on the Galaxy Map | ADM.94, MAP.218 |  |

### Admin shell

| ID | Item | Needs | Note |
|---|---|---|---|
| TEST.169 | Admin page browser tests: every admin page and task works end to end after the streamline | TEST.168 | Quality; queued admin chain |

### Help framework

| ID | Item | Needs | Note |
|---|---|---|---|

### Admin pages in the shell

| ID | Item | Needs | Note |
|---|---|---|---|
| POP.13 | A species editor: create and edit species and their environment limits |  |  |

## Phase 3: Editors that make plans

Each editor produces a plan: move a system, edit an orbit, adjust the world to a density change, and the rearranger that fixes a system around a change.


### Move and orbit

| ID | Item | Needs | Note |
|---|---|---|---|
| ADM.32 | Add a star system to a sector: at the emptiest spot, at given coordinates, or at random outside every Hill sphere | ADM.65, ADM.66, ADM.82 | Merges two asks (one system added by the computer; placement chosen by admin). |
| ADM.71 | Move a star system or phenomenon within its sector or to another sector, as a plan | ADM.32, ADM.65 |  |
| ADM.36 | Change an object's trajectory vector | GEN.109, ADM.65, ADM.66 |  |
| ADM.72 | Edit a body's orbit: distance, eccentricity, inclination and angles, with its moons following, as a plan | ADM.36, ADM.52, ADM.65 |  |

### Density plans

| ID | Item | Needs | Note |
|---|---|---|---|
| ADM.74 | Adjust the world to a density change: a plan that adds or removes systems so a region matches its new density | ADM.65, MAP.218 |  |

### The rearranger

| ID | Item | Needs | Note |
|---|---|---|---|
| ADM.75 | Rearrange a system around a change to a planet, moon or belt: the plan that fewest changes make valid again | ADM.64, ADM.65 | Boss's own example: change an A to an M planet |
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

### Civilization editors

| ID | Item | Needs | Note |
|---|---|---|---|
| ADM.130 | Aspects for civilizations: species, nations, cities, ports and stations in the aspect catalogue, with keep marks and locks | POP.22, ADM.62, ADM.63, DB.104, POP.26 | Boss, 9:27 PM EDT 2026-10-24: settlements tree; fixed until 9.2. |
| ADM.131 | Editable civilizations: edit a species, rename or split a nation, move a city, transfer a city between nations, set the dominant nation and place or remove a port, as approved plans | ADM.130, ADM.66, ADM.69 | Boss, 9:27 PM EDT 2026-10-24: settlements tree; fixed until 9.2. |

### Resource editors

| ID | Item | Needs | Note |
|---|---|---|---|
| ADM.133 | Admin edits of resource locations: change a belt, field or body's resource grades and scores as an aspect of the edit-everything plan, with a lock so regeneration keeps them | DB.106, GEN.799 | Boss, 10:29 PM EDT 2026-10-24: resource locations (valuable belts and fields as points of interest), foundations 8.2 to 9.0. |


## Phase 5: Close-out: proof, help and the headline

End-to-end and property tests, the admin help, and the items that close ADM.35 and ADM.57.

Release gate (standing, Boss, 2026-10-21, 5:10 PM EDT): this version is not stamped until SEC.33 (the full security code audit) and TEST.261 (the full bug hunt audit) each have a dated report for it and no finding they mark blocking is open. Both are standing items: they stay open and are re-run for every major or minor version change (required whenever the major or minor version changes (8.1, 8.2 and on to 9.0 and 10.0), never for build numbers).


### Proof and help

| ID | Item | Needs | Note |
|---|---|---|---|
| TEST.175 | Edit everything end to end: the headline scenarios with approval and undo | ADM.77, ADM.78, ADM.71, ADM.72, ADM.74, ADM.80, TEST.168 |  |
| TEST.176 | Plan engine properties: any sequence of edits leaves every system valid, and undo restores every row | ADM.76, ADM.60, ADM.77 |  |
| TEST.177 | Edited sectors survive repair, regeneration and the golden seed | ADM.61, ADM.80 |  |
| TEST.443 | Reclassification after an edit: an edited planet or moon changes class when its conditions cross a rule, and keeps it when they do not | DB.90, ADM.80 | Remainder of TEST.301 |
| DOC.54 | Admin help: how to edit, regenerate, lock, approve and undo | TEST.175, ADM.68 |  |
| ADM.35 | Full control from every screen: edit anything, regenerate with every input, backfill or erase what is in view | MAP.120, ADM.77, ADM.78, ADM.79, ADM.80, ADM.81 | Closes when the regen and edit screens of phase 4 are done; MAP.120 is the map part. |
| ADM.57 | Headline of 8.6: edit everything, the admin can change anything and the world adjusts after the admin approves a generation plan | ADM.32, ADM.35, ADM.36, ADM.52, MAP.119, MAP.120, POP.13, ADM.68, ADM.71, ADM.72, ADM.74, ADM.77, ADM.78, ADM.79, ADM.80, ADM.81, TEST.175, TEST.176, TEST.177, DOC.54 | The headline; done when its items are |
