# Orbital update check against Boss's four steps (2026-10-18)

Read-only check, no code changes. Written after Boss's message of 5:24 PM EDT on 2026-10-18 (this lane's clock). Sources: docs on main at 440db88 (`docs/design/orbital-updates.md`, `orbital-solvers-and-integrators.md`, `Orbital Update Full Algorithm.md`, `galactic-potential.md`, `gravity-map.md`, `mass-object-class.md`, `docs/plan/version-8-3-moving-galaxy.md`, `docs/TODO.md`), the drafts in this folder (`unified-system-barycentre.md`, `unified-star-like-objects.md`, `mass-object-class-savings.md`) and the code in `src/planetgen/db/store.py`, `physics/position.py`, `physics/field_kernel.py`, `physics/sector_path.py`. Nothing was run; every code statement is from reading the source.

Boss's message, quoted:

> And let's also double-check that the orbital motion updates:
> 1. Checks for whatever records require an update.
> 2. New records without a next update field will always get an update, where the next update field will be set.
> 3. Calculate any required position changes first.
> 4. Change the motion vector of each object to account for gravitational effects, given the point mass, calculated gradient of the sector and area surrounding the sector, using the algorithm of: the largest hill radius present in the galaxy becomes the search radius for gravitational objects. So that we can calculate the vector field and the gravitational gradient at the exact location of the object to update its motion vector
>
> Based on speed the system will then set an entry in the database for that object for when it next needs a position or motion update. According to the floor, it is already approved

## Summary table

| Item | Verdict | Where |
|---|---|---|
| Step 1, find records that need an update | Covered, and built | orbital-updates.md section 3 and section 4 ("A ceiling on a night's work"); `store.py` selects `next_update_due IS NULL OR next_update_due <= end` ordered by `next_update_due` (lines 5373, 5456, 5517, 9581), indexed column, GEN.106 |
| Step 2, new records always get an update and a next-update value | Partly covers, and the built code differs | orbital-updates.md section 3 ("first works out a due time for rows that have none"); `store._schedule_updates` (9211) and `_due_interval_s` (9237) |
| Step 3, positions first | Partly covers; the ordering Boss states is not written anywhere | Full Algorithm "Master Simulation Workflow"; GEN.340 (kick-drift-kick); built code moves position and velocity in one analytic step |
| Step 4, vector change from point masses plus gradient, search radius = largest Hill radius in the galaxy | Partly covers, contradicts on two points, not built | GEN.109 children GEN.341, GEN.342, GEN.343; GEN.198; Boss 2026-10-09 rule; card S13 |
| Final entry, next update set from speed | Covered for position, silent on motion-vector updates | orbital-updates.md section 3; `position.update_interval_s` (physics/position.py:93); GEN.139 |
| "According to the floor, it is already approved" | Could not tell | see the last section |

## Step 1: check which records need an update. Covered, built.

Each moving row has `epoch_unix` and an indexed `next_update_due` (GEN.106, schema v66, listed in orbital-updates.md section 3 and the "Status" line). A run selects only rows with `next_update_due IS NULL OR next_update_due <= end of the span`, oldest due first. A row ceiling (PERF.103, `tuning.ORBIT_UPDATE_ROW_CEILING`, 10,000,000) can defer the rest; deferred rows stay due and go first next run. Nothing is skipped silently.

Two notes that touch the new standing rule (Boss, 5:20 PM EDT, GEN.565 and GEN.584):

- The check is not one function. It is repeated in `advance_galactic_positions`, `advance_orbital_phases` and `advance_comet_orbits`, and for about a dozen tables, each with its own epoch and due columns. GEN.584 (SY17) puts epoch and next-due in the core's kinematic record, so after GEN.565b and GEN.584 there is one check for every kind of motion. That agrees with Boss's rule. GEN.565a (#2274) only carries the galactic-motion fields, so the one-check outcome is not there until the later slices.
- The 8.3 plan page still says the whole orbital update is "Proposed by the TODO thread; not yet approved by Boss" (`docs/plan/version-8-3-moving-galaxy.md`, Goal).

## Step 2: records without a next-update field. Partly covers; the built code does something different.

Design text: a run "first works out a due time for rows that have none (new, or an orbit edited since)" (orbital-updates.md section 3). Built: `_schedule_updates` and `_due_interval_s` give such a row `epoch + interval` (epoch taken as the last run's time when the row has none), then the row is moved only if that due time is at or before the end of the span (`store.py:5384`, 5460, 5521, 5632). So today a new record is scheduled, not necessarily updated. Boss's wording is that it "will always get an update". That is a difference, not a contradiction of the design's intent.

One more point I could not test: a NULL epoch means "at the last run" (position.py docs; `clock.start_unix`). A record created after the last run has its position at its creation time, so treating its position as holding at the last run's time would advance it over time before it existed. Whether generation stamps an epoch on new rows was not checked here.

## Step 3: position changes first. Partly covers.

- The Full Algorithm (Boss's imported document) does half-kick, drift, then second half-kick per object, stages every new position and commits atomically, with the sector point-mass tables read at the start (so the pull comes from the old positions of the other objects).
- The design replaced that with one adaptive integration per object over its whole interval (orbital-updates.md section 4 and 10.4) and GEN.340 builds a kick-drift-kick stepper on arrays. Neither is "all positions first, then all vectors".
- The built code (`advance_galactic_positions`, `advance_orbital_phases`) moves position and velocity together in one analytic step (rotation about the galactic axis for stars, `kepler.universal_step` for planets, moons and comets) and applies no neighbour pull at all.
- A not-due object is not moved, so its stored position is up to one threshold old (0.01 mpc, 0.01 AU or 100,000 km). That is inside the design's own definition of "has not moved noticeably", so using stored positions of not-due objects as sources is consistent with Boss's thresholds.

## Step 4: vector change from gravity. Partly covers, contradicts on two points, not built.

Covered:

- "Point mass of the objects found plus the galactic gradient, summed at the object's own position": GEN.342 (acceleration kernel, "pairwise point masses ... plus the galactic gradient of GEN.115's module, leaving out the object itself", then the velocity update) and the potential in GEN.338 and `galactic-potential.md`. Vector fields are "never stored" (orbital-updates.md section 4), which matches "at the exact location of the object". GEN.198 is the same evaluator for any point (pull, potential, tidal gradient) and is meant to wrap GEN.342's kernel.
- Per-sector point-mass tables with the largest mass per sector exist (`sectors.max_mass_kg` and `max_mass_id`, DB.55, schema v102, `schema.sql:1469`), and the sector search that matches brute force is in `orbital-solvers-and-integrators.md` section 6.4 (GEN.341).
- Hill-radius warnings by the gravitational criterion: GEN.344.

Contradicts or differs:

1. **Search radius.** Boss, 2026-10-09 (relayed, TODO GEN.109): "the influence radius is the maximum Hill sphere of the largest nearby object". Solver doc section 6.3: take the largest mass within a first scan of about 6 pc, grow the scan to a fixed point, and add cover lists for objects whose Hill sphere is wider than 1.5 sector edges. Boss's message of 5:24 PM EDT says the largest Hill radius present in the galaxy is the search radius. That is a single global number, which is simpler (no fixed point, no cover lists) but much larger. Scaling the table in solver doc section 6.1 (1 solar mass: 1.37 pc Jacobi radius at 8.128 kpc, growing as mass^(1/3)), a 4.3e6 solar-mass Sgr A* treated like any object would give about 220 pc (my scaling, not a run), which is of order 700,000 sectors of 4 pc per object. The same section's cost table shows the explosion already at 1e4 solar masses (2,359 sectors). The central black hole is also already in the smooth potential (GEN.285, realistic-galactic-core.md), so counting it again as a point mass would double count. GEN.370 (open, "Hill spheres ... use the Milky Way's mass and Sol's distance for every galaxy") must be fixed first or the maximum is computed with the wrong enclosed mass.
2. **Which objects get the vector change.** Boss says "each object". Card S13 (tapped 12:30 AM EDT on 2026-10-18, recommended option, `unified-system-barycentre.md` section 10.5 and 13) says trajectories inside a system are an analytic two-body orbit per frame, and the outer system's pull is a validity condition (Hill fraction, triple and Kozai screens), "not a force in this model"; the orbital update "never integrates". The galactic pass in orbital-updates.md section 4 also applies influence only to root objects; planets and moons get Kepler steps and the perturbed pass (GEN.343) only when a perturber is inside the system (GEN.290). But Boss's own design of 2026-10-15, 11:06 PM EDT (quoted in GEN.290) is that "a planet or moon is moved by the summed gravitational influence of point masses ... plus passing objects", so the all-objects reading has precedent and S13 is the newer, narrower card. The two are not reconciled in the docs.

Silent:

- The word "gradient". In Boss's earlier texts (GEN.115, 2026-10-07) it means the galaxy's smooth gravitational field, which is what changes velocity. In `gravity-map.md` it means the gradient of the pull (the tidal tensor), which does not change a velocity at first order and is drawn as a heat map. I read step 4 as the first meaning (acceleration from the point masses plus the smooth galactic pull). The tidal tensor stays a display and screening quantity.
- The ordering between objects: whether the pull on one object uses the other objects' old or new positions (see fork 3).

## Final entry: next update from speed. Covered for position, a gap for vector updates.

Built and documented: after an object is moved, `position.update_interval_s(speed, scale)` (threshold divided by speed, capped at 1e9 years; `MAX_UPDATE_INTERVAL_S`, position.py:88) gives its next-due time (orbital-updates.md section 3). Gaps:

- The time is set from the speed. Step 4 changes the speed, so the next-due time has to be worked out after step 4 from the new vector. The docs do not say so, and the built code has no step 4.
- A speed-based interval says when a position crosses its threshold, not when the vector needs refreshing. An object at rest, or moving slowly, with a massive object approaching is not revisited until its own long interval ends. Solver doc section 5, recommended scheme item 2, makes next-due the earliest of sector exit, the time for the influence set to change (`R_infl / (10 v_rel)`), one period for a perturbed bound body, and the literal threshold time; GEN.139 (Boss, 2026-10-10 2:44 PM EDT, defaults approved) took the analytic-orbit part of that. Boss's message says "based on speed" only.
- Boss says "when it next needs a position or motion update" with one entry. The design has one `next_update_due` for both; if he means two separate times (position due, vector due) that is a different schema.
- Applied literally to the galaxy-frame speed, 0.01 mpc makes every star due every 16 days at 220 km/s (solver doc section 5). GEN.139 (c) leaves open whether the threshold uses galaxy-frame or peculiar speed.

## "According to the floor, it is already approved"

I could not tell what this means, and I did not treat it as approval of anything. Possible referents I can list without choosing: the word "floor" appears in this project as the 10 m representability floor of the mass object class (card E2), the 1 m/s focusing floor in the collision formula, and the per-body `min_update_interval_years` floor in the code. A voice-typing slip for another word (for example "flow", "file" or "plan") is also possible. What I could verify as approved by Boss: the 2026-10-09 "largest nearby object" influence rule, the 2026-10-10 defaults for GEN.109 and GEN.139 ("stars follow their analytic orbit and are rewritten at sector exit or when the influence set changes"), the real-time step, the 2026-10-15 design of planets and moons moved by summed pull, and cards S10 to S18. What I did not find recorded as approved: the global largest Hill radius as search radius, the "positions first, then vectors" order, step 2's "always get an update", and the speed-only next-due rule. The 8.3 plan page says the version is "not yet approved by Boss".

## Forks for Boss (recommendation first)

1. **Search radius for step 4.**
   - (a) Recommended: use Boss's rule as the broad-phase bound, one stored number, the largest Hill radius of any object in the galaxy; leave the central black hole out of that maximum because it is in the smooth potential (GEN.285); sum the nearest 64 exactly and the rest as one point mass per sector (the existing `max_exact_members` cap); fix GEN.370 first so the radius uses the galaxy's own enclosed mass.
   - (b) Literal: every object inside the global radius contributes exactly, central black hole included. Hundreds of pc around every object; not workable at the galaxy sizes in the plan.
   - (c) Keep the 2026-10-09 "largest nearby object" rule with cover lists.
2. **Who gets the vector change.**
   - (a) Recommended: every object goes through the same core function (Boss's 5:20 PM EDT rule). A bound planet or moon works in its host's frame with the host's pull taken out (GEN.342 already says so); when nothing but the host matters it takes the analytic Kepler step; S13's tide screens stay as validity checks.
   - (b) Stay with S13 as tapped: only root objects (stars, remnants, rogues) get the field; bound members stay analytic and use the perturbed pass only when a perturber is inside the system.
   - (c) Integrate every system numerically every update. Closest to Boss's literal wording, most expensive.
3. **Order of positions and vectors.**
   - (a) Recommended: two phases per run. Phase 1 advances every due object's position. Phase 2 computes each due object's pull at its new position from the sources' new positions (not-due sources at their stored positions), then changes the vector. Deterministic and independent of the order rows are visited; GEN.340's kick-drift-kick still serves as the integrator inside each phase.
   - (b) Per object, position then vector, using whatever the other objects' stored positions are at that moment. Order-dependent.
   - (c) Full Algorithm style: every pull read from the old positions, everything committed together.
4. **When an object is next due.**
   - (a) Recommended: Boss's speed rule applied to the vector after step 4, capped at the time for the influence set to change (so a slow object is revisited when a massive object approaches).
   - (b) Speed only, as Boss wrote it. Simplest; a nearly stationary object near a passing mass can miss the pass until its own interval ends.
   - (c) The full earliest-of rule in solver doc section 5.
5. **New records without a next-update value.**
   - (a) Recommended: always visited on the next run: next-due is set, step 4 runs once, and the position only advances over the time since the record's own epoch (generation stamps the epoch), never over time before it existed.
   - (b) As built: scheduled, moved only if the computed due time falls inside the span.
   - (c) Always advanced from the last run's time, as the NULL-epoch rule reads literally.

## Decided (Boss, 2026-10-18)

Card taps between 5:30 and 5:33 PM EDT, as the coordinator reported them to the TODO thread (the cards are in the project chat; their tap results are not readable by the TODO thread). The last card, approval of the orbital motion update to build, was answered "Approved as the cards say".

1. Search radius: option (a), the largest Hill radius of any object in the galaxy, central black hole left out of the maximum, 64 exact then one point mass per sector, GEN.370 first.
2. Who gets the vector change: option (a), every object through the one core function (amends card S13).
3. Order: per object steps (option (b)): each object moves, then gets its vector, in turn. The visit order is tier, then mass, and each chunk reads the state at its start; that chunk-start rule is the only rule (section "Decided: visit order"), so the visit order must be fixed and tested. Decided (Boss, contradiction sweep card, 7:30 AM EDT 2026-10-21): the chunk-start rule stands and the "mix of old and new positions" wording is removed. The recommended two-phase option was not chosen.
4. Next due: option (a), speed from the new vector, capped at the time for the influence set to change.
5. New records: option (a), always visited, position advanced only from the record's own epoch, then next-update set.

The "floor" remark stays unexplained; the approval card answered the question it raised.

## Decided: visit order (Boss, 2026-10-18)

Boss's words, 5:38 PM EDT: move every due object by mass first (a black hole before a planet), then hosts first; the chain is primary (the main object), secondary (planets), tertiary (moons), quinary (small moon satellites, if they exist); galactic-frame objects (a rogue planet, an interstellar comet, a single primary star) move before anything inside a system. His card tap at 5:41 PM EDT, as the coordinator reported it: "Tier first, then mass". The order is a stable sort of the due objects on (tier, mass heaviest first, ties by id). Tier is the number of hops from the galactic frame, derived from the parent link and not stored: 0 top-level objects (each system's primary included), 1 planets and companion stars, 2 moons, 3 small moon satellites; a rogue planet's moon is tier 1. Rows run in ordered chunks (tier, then mass band), each seeing the state at its start. A planet heavier than a companion star stays in its tier and is ordered by mass there; a capture in the middle of a run takes effect in the next run's order. Open measurement: sort cost on millions of due rows.
