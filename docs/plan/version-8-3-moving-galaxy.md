# Version 8.3: The moving galaxy

Rebuilt on 2026-10-16 as a version map (Boss, 1:07 PM EDT: one version per headline, phases inside a version, blockers one phase earlier).
[docs/TODO.md](../TODO.md) is the master file: it holds each item's full text, and its "Plan: versions" section indexes every version.
This file gives the version's goal, its phases and each phase's build threads with every item's prerequisites in order. A phase holds independent threads: an item's blockers sit in the same thread (earlier in the list) or in an earlier phase. Where the two disagree, TODO.md wins; when an item ships, it leaves TODO.md and its row here is deleted in the same PR (the progress ledger keeps its version). Research notes and judgment calls are in [notes.md](notes.md); the headlines and the ideas not yet filed are in [release-roadmap.md](release-roadmap.md).

Version: 8.3

## Goal

The orbital update (GEN.105): every object moves on its orbit each day with n-body influence, collisions and ejections, built on one point-in-space object for every kind of object, with the gravity field evaluator, the Quantity class, long steps, events and the nightly run. Proposed by the TODO thread; not yet approved by Boss.

## Phase 1: Positions and the nightly run

Every object as one SpatialPosition3D, the nightly schedule and the movement layers on the maps.

### Positions as one object

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.115 | The galaxy's own gravity: a smooth disk, bulge and halo potential |  | Model from Boss's "Computational Astrodynamics.md" (2026-10-07): bulge, disk and halo potential. |
| GEN.142 | Peculiar velocity for rogue planets and asteroid fields |  | Research: prerequisite of GEN.110. |
| GEN.289 | Planets and moons kept as the universal point-in-space object, with point mass and movement vector stored, so perturbations have somewhere to land | GEN.115 |  |
| GEN.294 | Every object in the galaxy is held as one SpatialPosition3D with its position, point mass and movement vector stored as columns |  | Boss 2026-10-15 11:12 PM EDT. |

### Nightly maintenance

| ID | Item | Needs | Note |
|---|---|---|---|
| OPS.16 | A daily maintenance script for Linux, macOS and Windows |  | scripts/maintenance.sh and .ps1: positional update, delta merge, rotation; lock; optional OPS.15 check. |
| OPS.17 | Install and update set up the daily maintenance schedule | OPS.16 | Same scripts as OPS.7/OPS.8/OPS.13 (install/update, deploy-common), after them. |
| OPS.30 | A lock helper for the maintenance run |  | Research: used by OPS.16 (ADM.20's admin merge was dropped). |

### Movement on the maps

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.212 | Galactic move in the orbital update as NumPy batches (3 to 4 times on the whole update) |  |  |
| MAP.192 | See the objects whose positions changed in the last orbital update |  |  |
| MAP.193 | Show every object's axis of rotation, direction of rotation and axial tilt |  |  |

## Phase 2: The update itself

N-body influence, collisions, encounters and the Quantity class.

### Orbital update core

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.109 | N-body influence from every object inside the largest nearby Hill sphere plus the galactic gradient, with a Hill-radius warning | GEN.115, GEN.289, GEN.294 |  |
| GEN.110 | Rogue planet collisions: asteroid fields, merged giants and new stars | GEN.109, GEN.142 |  |
| GEN.139 | Orbit-update thresholds: per-object epoch, path-length rule and what the 0.01 mpc applies to (GEN.106 built) |  | Research follow-up to GEN.106 (built). |
| GEN.140 | Orbital math guards the edge-case table adds (GEN.108 built) |  | Research follow-up to GEN.108 (built). |
| GEN.141 | Faster Kepler solver (Mikkola or Markley) with brentq as fallback |  | Research: optional speed-up. |
| GEN.144 | One distance for the Sun from the galactic centre across the constants, the density model and the design docs |  | Research: three Sun distances today. |
| GEN.290 | An encounter search that cannot miss a pass at large steps (100 years per second) | GEN.109, GEN.289 |  |
| GEN.291 | A bound planet becomes a rogue planet when an encounter unbinds it: the unbind and capture test and the event log | GEN.109, GEN.289, GEN.290, GEN.110 |  |
| GEN.105 | Orbital updates | GEN.115, GEN.109, GEN.110, GEN.289, GEN.290, GEN.291 | Parent of the orbital update work. |

### The Quantity class

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.299 | The Quantity class (the unit class): one value, every unit version always available, all conversions and the shared ladders | GEN.294, GEN.135 | Boss 2026-10-16, 11:35 AM EDT; like SpatialPosition3D. |
| GEN.300 | Every quantity in the code moves onto the Quantity class, and hand conversions are removed | GEN.299, UX.107 | The call-site sweep. |

## Phase 3: Events, long steps and the evaluator

The event records, the long-step rules and the gravity field evaluator the courses reuse.

### Events and ejections

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.111 | Email the admin when two objects are inside each other's Hill radius | GEN.109, USR.3 | Needs USR.3's SMTP settings. |
| GEN.143 | A collision_events table, an admin report and a test that runs the whole collision path | GEN.110 | Research: makes the collision code testable. |
| GEN.197 | Object IDs on ejection, merger and split events | GEN.143 |  |

### Long steps

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.292 | Run the orbital update for a stated span in one go from the command line (elapsed years) | GEN.291 |  |
| GEN.293 | Long-step rules for hypervelocity stars that leave the sector grid and for parabolic comets | GEN.290 |  |
| GEN.295 | Nebulae as point groups whose shape follows the orbital update (design and research) | GEN.294, GEN.109, GEN.292 | Boss 2026-10-16 1:01 AM EDT. Design only, no lane. |

### Gravity field evaluator

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.198 | The gravity field evaluator: pull, potential and tidal gradient at any point from the point masses and the galaxy potential | GEN.115, GEN.109 | Gravity map foundation (Boss 2026-10-10). |
| TEST.131 | Tests for the gravity field evaluator: two bodies, Lagrange points, softening and the rotation curve | GEN.198 |  |
