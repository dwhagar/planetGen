# Edit everything: the admin changes anything and the world adjusts after a plan is approved

**Status (2026-10-16):** plan; not built. This is the design of the 8.6 headline (ADM.57). It builds on what exists: the admin edits in `planetgen/admin/edits.py` (ADM.6, ADM.7, ADM.8: regenerate or remove a body, change a class, change a single star), their writer in `planetgen/db/edits.py`, `validation.stabilize_star_system` (moves bodies until a system validates, never removes them) and the activity log.

Boss (2026-10-16, 1:07 PM EDT): "the admin has complete control they can edit and change anything how they want. They can move stars in a sector, change orbits in a system, change densities for a sector or region, regenerate any star they want to any specification they want. The world will adjust to what the admin wants. That means if I want to change an A to an M planet the system will figure out how to rearrange the star system to make that happen but it will get approval for the new generation plan first. In regen screens the admin should also be able to select aspects they want to keep and what they want regenerated."

## 1. The idea in one paragraph

Every admin change, small or large, becomes a **plan**: an ordered list of operations on stored objects. The planner computes the whole effect (every row that changes, including the bodies that have to move so the system stays valid), shows it, and nothing is written until the admin **approves**. An approved plan is applied as one visible job, journaled, and can be undone. The admin steers regeneration with **aspects**: for each aspect of an object (a star's class, a planet's orbit, a system's bodies) the admin marks keep, regenerate, or set to a value, and a **lock** is a keep mark that stays.

## 2. The parts

| Part | What it is | Items |
|---|---|---|
| Change plan | Operations, the effect diff, the states draft, proposed, approved, applied, undone, stale, rejected, and the stale check (an object changed since the plan was proposed) | ADM.58, DB.36, ADM.59 |
| Apply, undo, audit | One queue job, one transaction per sector, a journal row and an activity-log entry per operation, undo as the inverse plan from the stored before values; refuses while the update guard is on | ADM.60 |
| Aspect catalogue | Data: every stored column belongs to one aspect, each aspect names the generation input that draws it, what it derives from and what depends on it | ADM.62, ADM.63 |
| Keep or regenerate | Resolves the admin's marks into the redraw set, the fixed inputs and the conflicts (a kept radius with a regenerated class is a conflict to show) | ADM.64 |
| Locks | A standing keep mark on an object or one aspect; shown on the pages and maps; every plan refuses to change a locked aspect unless the plan lifts the lock | DB.37, ADM.65 |
| Unit finalizer | Validates a system or sector, writes it in one transaction, auto-fixes only a clashing name and derived columns, or rejects with reasons; plan apply and upload verification both end here (shared with 8.9) | ADM.82 |
| Recipe schema | For each aspect a fixed value, a range or "draw it"; the regeneration specification and generate-by-recipe are the same schema (shared with 8.9) | API.28 |
| Journal replay | The generated sector stays the base (GEN.57); the journal of applied operations replays over it after a rebuild from the seed (DB.17) or from parity (DB.9), so the rebuilt sector equals the edited one | ADM.61, TEST.177 |
| Regions | A region (sphere, box, shell, band, wedge, set of sectors), a selection tool, and density overrides for any region | GEN.302, MAP.212, DB.38, ADM.73 |
| Rearranger | The planner that makes a system valid again around a change to a planet, moon, belt or star, with the fewest changes | ADM.75, ADM.76 |
| Screens | Approval, before-and-after diagram, plan history with undo, the keep-or-regenerate matrix, the specification form, and one regen screen per kind | ADM.66 to ADM.70, ADM.77 to ADM.81 |
| Editors | Move a system or phenomenon, edit an orbit, adjust the world to a density change, regenerate a sector or region | ADM.71, ADM.72, ADM.74, ADM.80 |

## 3. The headline cases

**Change an A planet into an M planet.** The admin opens the planet, chooses class M in the regen screen and leaves the rest as marked. `rearrange` (ADM.75) takes the change and the keep marks: it asks which bodies conflict with an M planet at that distance (spacing, Hill spheres, the habitable zone, moons that cannot be kept), picks the smallest set of moves, resizes and regenerations that restore validity, and returns a plan. The approval screen shows the system before and after (ADM.67). Approve, and the job applies it.

**Move a star in a sector.** A plan operation by coordinates (the three coordinate forms of ADM.32), by drag on the Sector Map, or to the emptiest spot. It checks every Hill sphere and cover list, updates sector membership, and moves nothing else unless the admin approves that too.

**Change an orbit.** The elements of a body (distance, eccentricity, inclination and angles) with the moons following; the velocity-vector edit of ADM.36 is the same operation in another form. The next orbital update starts from the stored result (8.3).

**Change a density.** The admin selects a region, sets a multiplier (0.1 to 10) and previews it. If the region already holds systems, the world adjusts through a plan that adds systems from the fill, or removes the fewest, never a locked, named or populated one without approval (ADM.74).

**Regenerate to a specification.** A specification form, generated from the aspect catalogue, plus the keep-or-regenerate matrix: keep the name, the position and one chosen planet, redraw the rest. The same component serves the object, system, sector and region screens.

## 4. Decisions and their reasons

- **One mechanism for every edit.** Moves, orbits, densities and regenerations all produce plans, so approval, history, undo and audit are written once. The cost is that a one-field edit also goes through a plan; the screen approves a one-operation plan in one click.
- **The planner never applies.** A plan is data until approved; applying re-checks it is not stale. Two admins cannot silently overwrite each other.
- **Rearrange by fewest changes, deterministically.** The planner draws from the object's own seed stream, so the same request on the same stored state gives the same plan. Alternatives (two or three candidate plans to choose from) are a later idea (S8.6-b).
- **Edits live in the database as a journal over the generated base.** GEN.59 (admin changes as a net difference) and GEN.61 (a daily merge into a JSON file) were dropped on 2026-10-09; this keeps the journal in the database and uses it only so repair and regeneration give the edited rows back. A journal row that no longer applies after a generator epoch change is listed, never dropped.
- **Aspects are data.** A new column added to a body appears in the matrix, the specification form and the generic field editor (ADM.81) without a new screen; a test fails if a stored column belongs to no aspect.
- **Locks are the long-lived keep marks.** A regeneration of a sector keeps locked objects and aspects without the admin re-marking them every time.
- **The update guard applies.** Applying a plan refuses while maintenance is on, an update is running or the orbital update runs, as OPS.41, OPS.53 and OPS.54 do for other admin jobs.
- **Build once what 8.9 (generate anywhere) also needs** (Boss, 1:33 PM EDT). Two pieces were in both: the check-and-write of a unit (API.8 and plan apply) and the recipe (API.18 and the specification form). They are built here, as ADM.82 and API.28, and API.8 and API.18 in 8.9 keep only what is upload- or route-specific. A plan operation uses the upload envelope (`StarSystem.to_dict()`, API.13), the region object (GEN.302) is also what galaxy-scale recipes (API.19) apply to, and apply will refuse a sector an upload run has claimed (API.10).
- **Migrations are separate items.** DB.36, DB.37 and DB.38 are migration-only items for Core 1; the code that uses the tables follows in later phases.

## 5. Phases and threads (version 8.6)

| Phase | Theme | Threads |
|---|---|---|
| 1 | Plan model, aspects and regions | Admin shell and audit; Existing admin editors; Change plans (with the unit finalizer); Aspects (with the recipe schema); Regions |
| 2 | Interfaces: approval, locks, regen screens, density overrides | Approval and history; Locks and regen screens; Density overrides |
| 3 | Editors that make plans | Move and orbit; Density plans; The rearranger |
| 4 | Regeneration and editing any field | Regenerate objects; Regenerate regions; Edit any field |
| 5 | Close-out: proof, help and the headline | Proof and help, then ADM.35 and ADM.57 |

A blocker of an item sits in an earlier phase or earlier in its own thread; the threads of one phase are independent, so they can be built side by side.

## 6. Open questions (each with a default; work goes ahead on the default)

- **Undo window.** Undo is allowed while no newer plan touched the same objects. Default: yes; a later plan on the same objects blocks the undo and says which.
- **Who may approve.** Any admin may approve a plan another admin proposed. Default: yes (there is one admin class today); two-person approval waits for the roles of 8.2.
- **Removal in a density plan.** The planner may remove generated, unnamed, unpopulated systems without asking per system, inside the approved plan. Default: yes, listed on the approval screen with counts per sector.
- **Specification form limits.** The form accepts any value the generator can draw plus "force" values the generator would not draw (as ADM.6 does for classes). Default: forced values allowed, marked on the plan as forced.
