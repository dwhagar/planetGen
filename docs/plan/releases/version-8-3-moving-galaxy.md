# Version 8.3: The moving galaxy

Rebuilt on 2026-10-16 as a version map (Boss, 1:07 PM EDT: one version per headline, phases inside a version, blockers one phase earlier).
[docs/TODO.md](../../TODO.md) is the master file: it holds each item's full text, and its "Plan: versions" section indexes every version.
This file gives the version's goal, its phases and each phase's build threads with every item's prerequisites in order. A phase holds independent threads: an item's blockers sit in the same thread (earlier in the list) or in an earlier phase. Where the two disagree, TODO.md wins; when an item ships, it leaves TODO.md and its row here is deleted in the same PR (the progress ledger keeps its version). Research notes and judgment calls are in [notes.md](../notes.md); the headlines and the ideas not yet filed are in [release-roadmap.md](../release-roadmap.md).

Version: 8.3

Status: paused until 8.1 is complete (Boss, 7:44 PM EDT 2026-10-22: the five cores work 8.1 only; nothing in this version is started or owned until then, except the gates already running for UX.140, UX.128 and NAV.17).

## Goal

The orbital update (GEN.105): every object moves on its orbit each day with n-body influence, collisions and ejections, built on one point-in-space object for every kind of object, with the gravity field evaluator, the Quantity class, long steps, events and the nightly run. The orbital motion update (four steps, forks F1 to F5) is approved to build (Boss, card taps of 2026-10-18, 5:30 to 5:33 PM EDT, as the coordinator reported them); the rest of the version is as proposed by the TODO thread.

## Phase 1: Positions and the nightly run

Every object as one SpatialPosition3D, the nightly schedule and the movement layers on the maps.

### Positions as one object

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.375 | Field kernel subtracts a local origin before its loop in `physics/sector_path.py` |  | GEN.372 audit. |
| MAP.307 | Stations and facilities become points of interest: interest reasons for the map mark and its menu filter |  |  |

### Nightly maintenance

| ID | Item | Needs | Note |
|---|---|---|---|
| OPS.16 | A daily maintenance script for Linux and macOS, with its lock helper |  | scripts/maintenance.sh and .ps1: positional update, delta merge, rotation; lock; optional OPS.15 check. |
| OPS.17 | Install and update set up the daily maintenance schedule | OPS.16 | Same scripts as OPS.7/OPS.8/OPS.13 (install/update, deploy-common), after them. |

### Movement on the maps

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.193 | Show every object's axis of rotation, direction of rotation and axial tilt |  |  |
| MAP.235 | Superbubbles as a registered Galaxy Map layer: a toggle, a hover name and radius, and a key row |  | Split from the Galaxy Map key item. |

### Galactic gravity

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.285 | Nuclear mass components in the galactic potential for the orbital update |  | Core 3; step 5 of the galactic core note; moved from 8.1 (no potential module until GEN.338); Boss 8:41 PM EDT: stays in 8.3. |

## Phase 2: The update itself

N-body influence, collisions, encounters and the Quantity class.

### Orbital update core

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.139 | Orbit-update thresholds: per-object epoch, path-length rule and what the 0.01 mpc applies to (GEN.106 built) |  | Research follow-up to GEN.106 (built). |
| GEN.339 | Generated and stored galactic velocities follow the potential: the old curve is replaced, dispersions are seeded, stored velocities are reseeded |  | Split from GEN.115. |
| GEN.340 | One star stepper: kick-drift-kick in the potential, a step per object, on arrays | GEN.139 | Split from GEN.115. |
| GEN.115 | The galaxy's own gravity: a smooth disk, bulge and halo potential | GEN.339, GEN.340 | Model from Boss's "Computational Astrodynamics.md" (2026-10-07): bulge, disk and halo potential. |
| GEN.212 | Galactic move in the orbital update as NumPy batches (3 to 4 times on the whole update) | GEN.340 |  |
| GEN.343 | The perturbed pass for a system: host, siblings, moons and outside perturbers integrated through an encounter |  | Split from GEN.109. |
| GEN.344 | Hill-radius warnings by the gravitational criterion, to the debug and activity logs |  | Split from GEN.109. |
| GEN.109 | N-body influence from every object inside the largest nearby Hill sphere plus the galactic gradient, with a Hill-radius warning | GEN.343, GEN.344 |  |
| GEN.197 | Orbital event record: the `record_event` writer and the object-ID rules for ejection, capture, merger and split |  |  |
| GEN.346 | Applying a collision to the galaxy: new asteroid field, merged giant, merger-remnant star, retired bodies, events | GEN.197 | Split from GEN.110. |
| GEN.110 | Rogue planet collisions: asteroid fields, merged giants and new stars | GEN.346 |  |
| GEN.290 | An encounter search that cannot miss a pass at large steps (100 years per second) | GEN.343 |  |
| GEN.291 | A bound planet becomes a rogue planet when an encounter unbinds it: the unbind and capture test and the event log | GEN.290, GEN.343, GEN.346, GEN.197 |  |
| GEN.105 | Orbital updates | GEN.115, GEN.109, GEN.110, GEN.290, GEN.291 | Parent of the orbital update work. |

### The Quantity class

| ID | Item | Needs | Note |
|---|---|---|---|

## Phase 3: Events, long steps and the evaluator

The event records, the long-step rules and the gravity field evaluator the courses reuse.

Release gate (standing, Boss, 2026-10-21, 5:10 PM EDT): this version is not stamped until SEC.33 (the full security code audit) and TEST.261 (the full bug hunt audit) each have a dated report for it and no finding they mark blocking is open. Both are standing items: they stay open and are re-run for every major or minor version change (required whenever the major or minor version changes (8.1, 8.2 and on to 9.0 and 10.0), never for build numbers).

### Events and ejections

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.111 | Email the admin when two objects are inside each other's Hill radius | USR.3, GEN.344 | Needs USR.3's SMTP settings. |
| GEN.347 | The last-collisions admin view and the test that runs the whole collision path | GEN.346 | Split from GEN.143. |
| GEN.143 | The orbital-event admin report, the last-collisions view and a test that runs the whole collision path | GEN.347 | Research: makes the collision code testable. |

### Long steps

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.292 | Run the orbital update for a stated span in one go from the command line (elapsed years) | GEN.291 |  |
| GEN.293 | Long-step rules for hypervelocity stars that leave the sector grid and for parabolic comets | GEN.290 |  |
| GEN.295 | Nebulae as point groups whose shape follows the orbital update (design and research) | GEN.109, GEN.292 | Boss 2026-10-16 1:01 AM EDT. Design only, no lane. |

### Gravity field evaluator

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.198 | The gravity field evaluator: pull, potential and tidal gradient at any point from the point masses and the galaxy potential |  | Gravity map foundation (Boss 2026-10-10). |

### What moved in the last update

| ID | Item | Needs | Note |
|---|---|---|---|
| DB.57 | The update records what it moved: a replace-each-run table of previous and new positions | GEN.212, GEN.340 | Split from MAP.192. Core 1 owns migrations. |
| MAP.216 | Admin page of the objects that moved in the last update, with the movement vector and a map ghost of the old position | DB.57 | Split from MAP.192. |
| MAP.192 | See the objects whose positions changed in the last orbital update | DB.57, MAP.216 |  |

## Phase 4: World population

The population simulation from habitability scores (Boss, 9:38 PM EDT 2026-10-24).

Release gate (standing, Boss, 2026-10-21, 5:10 PM EDT): this version is not stamped until SEC.33 (the full security code audit) and TEST.261 (the full bug hunt audit) each have a dated report for it and no finding they mark blocking is open. Both are standing items: they stay open and are re-run for every major or minor version change (required whenever the major or minor version changes (8.1, 8.2 and on to 9.0 and 10.0), never for build numbers).

### World population

| ID | Item | Needs | Note |
|---|---|---|---|
| POP.30 | World population simulation from habitability scores: population, carrying capacity and growth state as its own command and galaxy stage | DB.105, POP.29, GEN.774 | Boss, 9:38 PM EDT 2026-10-24: trade routes and world populations, done by 8.9, dynamic in 9.1. |
| POP.31 | Population on the pages and API: world population, capacity and share of capacity on the planet, species and polity pages | POP.30 | Boss, 9:38 PM EDT 2026-10-24: trade routes and world populations, done by 8.9, dynamic in 9.1. |

### Resource locations

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.799 | Resource scan: compute the grade, mass and score of every resource for belts, fields, giants and icy bodies from their stored composition, as its own command and galaxy stage | DB.106, GEN.798 | Boss, 10:29 PM EDT 2026-10-24: resource locations (valuable belts and fields as points of interest), foundations 8.2 to 9.0. |
| MAP.312 | Resource locations become points of interest: interest reasons for valuable asteroid belts and fields, icy bodies and giants, each with its own row in the map menu filter | GEN.799 | Boss, 10:29 PM EDT 2026-10-24: resource locations (valuable belts and fields as points of interest), foundations 8.2 to 9.0. |
