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
| GEN.335 | One definition of the stored point columns, and a test that lists every object table |  | Split from GEN.294. |
| DB.55 | One additive migration for the 8.3 schema: point columns for every kind, heaviest object per sector, asteroid-field provenance |  | Split from GEN.294. Core 1 owns migrations. |
| GEN.336 | Compact and free objects hold a velocity vector and a mass: stars, black holes, neutron stars, rogue planets, interstellar comets, asteroid fields, nebulae, remnants | GEN.335, DB.55 | Split from GEN.294. |
| GEN.337 | Epoch, vector and mass for the rest: belts, quasars, comets, facilities, scatter points and members | GEN.335, DB.55, GEN.336 | Split from GEN.294. |
| GEN.294 | Every object in the galaxy is held as one SpatialPosition3D with its position, point mass and movement vector stored as columns | GEN.335, DB.55, GEN.336, GEN.337 | Boss 2026-10-15 11:12 PM EDT. |
| GEN.140 | Orbital math guards the edge-case table adds (GEN.108 built) |  | Research follow-up to GEN.108 (built). |
| GEN.289 | Planets and moons kept as the universal point-in-space object, with point mass and movement vector stored, so perturbations have somewhere to land | GEN.335, DB.55 |  |
| GEN.142 | Peculiar velocity for rogue planets and asteroid fields | GEN.336 | Research: prerequisite of GEN.110. |

### Nightly maintenance

| ID | Item | Needs | Note |
|---|---|---|---|
| OPS.30 | A lock helper for the maintenance run |  | Research: used by OPS.16 (ADM.20's admin merge was dropped). |
| OPS.16 | A daily maintenance script for Linux, macOS and Windows | OPS.30 | scripts/maintenance.sh and .ps1: positional update, delta merge, rotation; lock; optional OPS.15 check. |
| OPS.17 | Install and update set up the daily maintenance schedule | OPS.16 | Same scripts as OPS.7/OPS.8/OPS.13 (install/update, deploy-common), after them. |

### Movement on the maps

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.217 | Map layer interface: one registry the Galaxy stage view, the Sector Map and the System Map share for a toggled layer (menu entry, modes, URL key, legend, picker layer, fade rule) |  | Foundation found in the review. |
| MAP.193 | Show every object's axis of rotation, direction of rotation and axial tilt | MAP.217 |  |

### Galactic gravity

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.144 | One distance for the Sun from the galactic centre across the constants, the density model and the design docs |  | Research: three Sun distances today. |
| GEN.338 | The galactic potential as one pure module: Hernquist bulge, Miyamoto-Nagai disk and NFW halo, with settings and the rotation-curve tests | GEN.144 | Split from GEN.115. |
| GEN.370 | Hill spheres, keep-out spheres and the default tidal distance use the Milky Way's mass and Sol's distance for every galaxy: take the galaxy's own enclosed mass and the cell's own distance (bug) |  | Found by the TODO thread. |

## Phase 2: The update itself

N-body influence, collisions, encounters and the Quantity class.

### Orbital update core

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.139 | Orbit-update thresholds: per-object epoch, path-length rule and what the 0.01 mpc applies to (GEN.106 built) |  | Research follow-up to GEN.106 (built). |
| GEN.141 | Faster Kepler solver (Mikkola or Markley) with brentq as fallback |  | Research: optional speed-up. |
| GEN.339 | Generated and stored galactic velocities follow the potential: the old curve is replaced, dispersions are seeded, stored velocities are reseeded | GEN.338, GEN.294 | Split from GEN.115. |
| GEN.340 | One star stepper: kick-drift-kick in the potential, a step per object, on arrays | GEN.338, GEN.294, GEN.139 | Split from GEN.115. |
| GEN.115 | The galaxy's own gravity: a smooth disk, bulge and halo potential | GEN.338, GEN.339, GEN.340 | Model from Boss's "Computational Astrodynamics.md" (2026-10-07): bulge, disk and halo potential. |
| GEN.212 | Galactic move in the orbital update as NumPy batches (3 to 4 times on the whole update) | GEN.340 |  |
| GEN.334 | Pair closest-approach library: sweep, pair windows, window cap and chunking |  | Foundation found in the review. |
| GEN.341 | Influence-set search: every object inside the largest nearby Hill sphere, found by sector bounds, with cover lists and a brute-force check | GEN.294, GEN.289 | Split from GEN.109. |
| GEN.342 | The acceleration kernel: the pull of the influence set plus the galactic gradient, as one function, and the velocity update from it | GEN.341, GEN.338, GEN.340 | Split from GEN.109. |
| GEN.343 | The perturbed pass for a system: host, siblings, moons and outside perturbers integrated through an encounter | GEN.342, GEN.289, GEN.140 | Split from GEN.109. |
| GEN.344 | Hill-radius warnings by the gravitational criterion, to the debug and activity logs | GEN.341, GEN.342 | Split from GEN.109. |
| GEN.109 | N-body influence from every object inside the largest nearby Hill sphere plus the galactic gradient, with a Hill-radius warning | GEN.341, GEN.342, GEN.343, GEN.344 |  |
| DB.56 | The orbital_events table (migration only): collision outcomes, mergers, splits, ejections, captures, unbinds and warnings |  | Split from GEN.143. Core 1 owns migrations. |
| GEN.197 | Orbital event record: the `record_event` writer and the object-ID rules for ejection, capture, merger and split | DB.56 |  |
| GEN.345 | Collision outcomes as a pure function: regime table, giant merger, fusion check, conservation |  | Split from GEN.110. |
| GEN.346 | Applying a collision to the galaxy: new asteroid field, merged giant, merger-remnant star, retired bodies, events | GEN.345, GEN.334, GEN.341, GEN.142, DB.56, GEN.197, DB.55 | Split from GEN.110. |
| GEN.110 | Rogue planet collisions: asteroid fields, merged giants and new stars | GEN.345, GEN.346 |  |
| GEN.290 | An encounter search that cannot miss a pass at large steps (100 years per second) | GEN.289, GEN.341, GEN.343, GEN.334 |  |
| GEN.291 | A bound planet becomes a rogue planet when an encounter unbinds it: the unbind and capture test and the event log | GEN.289, GEN.290, GEN.343, GEN.346, GEN.197 |  |
| GEN.105 | Orbital updates | GEN.115, GEN.109, GEN.110, GEN.289, GEN.290, GEN.291 | Parent of the orbital update work. |

### The Quantity class

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.299 | The Quantity class (the unit class): one value, every unit version always available, all conversions and the shared ladders | GEN.294, GEN.135 | Boss 2026-10-16, 11:35 AM EDT; like SpatialPosition3D. |
| GEN.348 | Describer, CLI query lines and class page ranges take Quantity | GEN.299 | Split from GEN.300. |
| UX.128 | Map side panels and other browser call sites take Quantity | GEN.299, UX.107 | Split from GEN.300. |
| GEN.349 | API output and stored-column readers take Quantity; the grep test fails on a hand-written conversion factor; the old helpers are removed | GEN.348, UX.128 | Split from GEN.300. |
| GEN.300 | Every quantity in the code moves onto the Quantity class, and hand conversions are removed | GEN.348, UX.128, GEN.349 | The call-site sweep. |

## Phase 3: Events, long steps and the evaluator

The event records, the long-step rules and the gravity field evaluator the courses reuse.

### Events and ejections

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.111 | Email the admin when two objects are inside each other's Hill radius | USR.3, GEN.344 | Needs USR.3's SMTP settings. |
| GEN.347 | The last-collisions admin view and the test that runs the whole collision path | DB.56, GEN.346 | Split from GEN.143. |
| GEN.143 | A collision_events table, an admin report and a test that runs the whole collision path | DB.56, GEN.347 | Research: makes the collision code testable. |

### Long steps

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.292 | Run the orbital update for a stated span in one go from the command line (elapsed years) | GEN.291 |  |
| GEN.293 | Long-step rules for hypervelocity stars that leave the sector grid and for parabolic comets | GEN.290 |  |
| GEN.295 | Nebulae as point groups whose shape follows the orbital update (design and research) | GEN.294, GEN.109, GEN.292 | Boss 2026-10-16 1:01 AM EDT. Design only, no lane. |

### Gravity field evaluator

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.198 | The gravity field evaluator: pull, potential and tidal gradient at any point from the point masses and the galaxy potential | GEN.338, GEN.342 | Gravity map foundation (Boss 2026-10-10). |
| TEST.131 | Tests for the gravity field evaluator: two bodies, Lagrange points, softening and the rotation curve | GEN.198 |  |

### What moved in the last update

| ID | Item | Needs | Note |
|---|---|---|---|
| DB.57 | The update records what it moved: a replace-each-run table of previous and new positions | GEN.294, GEN.212, GEN.340 | Split from MAP.192. Core 1 owns migrations. |
| MAP.216 | Admin page of the objects that moved in the last update, with the movement vector and a map ghost of the old position | DB.57, MAP.217 | Split from MAP.192. |
| MAP.192 | See the objects whose positions changed in the last orbital update | DB.57, MAP.216 |  |
