# One galaxy-wide object store, region work units, and sectors as a navigation table (GEN.1118)

Status: final enough to file from (5:44 PM EDT 2026-10-30). Items filed as GEN.1123 to GEN.1129, TEST.501, with GEN.1111, GEN.1120 and GEN.1121 reworded. Boss's answers on tiles, bands, rogue ID after a move and the top-up ceiling are open, so the defaults of section 12 are taken. Section 11's item numbers are the working titles: 1 and 2 are GEN.1127, 3 is GEN.1111, 3b is GEN.1120 with GEN.1124 to GEN.1126, 3c is GEN.1121, 4 is GEN.1128, 5 is TEST.501. Item 6 (cross-face separation) waits on GEN.1006 and item 7 (partial sectors) is only built if Boss asks.

Design study, written about 5:35 PM EDT 2026-10-30 as the lane reported. Research only, nothing built. Tags: **[S]** read in the repository, **[C]** computed here by arithmetic from numbers in the repo, **[I]** inferred, **[U]** not measured or not read.

Boss's ideas, in order: (1) 5:15 PM EDT, store everything in one database so an object's place does not matter, sectors and regions are for navigation and are calculated when asked; build a whole region from one starting draw regardless of sectors, work units as bands within the sphere; a sector table with a charted flag and an uncharted boundary in spherical coordinates from the sector centre. (2) 5:16 PM EDT, a rogue body is stored in the sector where it ends up, that sector listed uncharted, and this applies to every object. (3) 5:18 PM EDT, generate the sectors rogues land in as a way to fill out the universe: charted sphere first, then phases that chart the sectors the earlier phase's rogues landed in, admin choosing how many phases.

## 1. Short answers

1. **The store is already one store.** Every positioned table carries galaxy-frame coordinates and a generated `cell_key` (Morton key, 0.25 pc cubes) with a covering index, and "what is within r of here" never reads `sectors` [S: `db/near.py`, `db/region.py`, spherical-region-queries.md section 6]. `star_systems.sector_id` is nullable with `ON DELETE SET NULL`. Two things still tie an object to its birth sector: `rogue_planets.sector_id` is `ON DELETE CASCADE` (deleting a sector deletes a rogue that may have drifted in), and the object ID's serial is the object's generation rank inside its birth sector [S: object-id-options.md section 3, `store._claim_object_ids`]. Both are fixable and the first is a one-line schema change.
2. **Uncharted space already holds known-unknowns.** `bright_stars` (every star above the luminosity floor, placed galaxy-wide by shell, with `star_system_id` NULL until a sector builds it), `phenomenon_scatter`, `galaxy_block_lights` and `sector_stats.level_before_fill_sol` are exactly "we know what might be there, not what is there" for stars and phenomena [S]. Boss's idea 2 extends that pattern to rogues. That is the new part, and it is cheap.
3. **The region is a request, the tile is the unit of randomness.** Determinism (one seed, same objects whatever order regions are charted) requires every draw to be keyed by a fixed global tile and candidate index, never by the region's start point. A region started from a random point is then only a choice of which fixed tiles to run, and "bands within the sphere" is a fine way to cut that list into work units (section 4). If the start point seeded the draw, two overlapping regions would disagree in the overlap.
4. **Recommended tile: keep the sector (whole sectors) for 8.1; do not make partial sectors.** Partial sectors need the sector fill to stop placing systems sequentially (it rejects against what is already placed in the sector, `SpaceSector._random_position` [S]) and the object serial to stop being "rank in the sector". That is a rewrite of the fill and the ID rule, and it buys only an exactly-spherical region edge, which rounding the region up to whole sectors gives to within 4 pc. Section 6 says what to do if Boss wants partial sectors anyway.
5. **Charted flag: derived, not drawn.** A `chart_ledger` of charted tiles (run, epoch, time) is the record; a sector's charted state and its uncharted remainder are computed from it. A stored spherical boundary per sector is not needed when tiles are whole sectors (a sector is wholly charted or wholly not), and is awkward when they are not (section 6).
6. **Boss's decision, 5:23 PM EDT: rogues come only from generated systems, and an uncharted sector stores every scatter and rogue object exactly as a charted one would.** This replaces virtual donors (GEN.1111 returns to real donors) and removes the rogue-landing phases entirely. Section 8 is rewritten for it. The one consequence to read first: the rogue population a region has is only what the charted stars have ejected (section 8.3).
7. **Generating an uncharted sector** fills around its existing residents and never moves or changes one (section 7.2). That is what `add_preplaced_system` already does for bright stars, so it is an extension, not a new mechanism.

## 2. What the code does today (the facts the answers rest on)

- **Two generators, one pipeline** [S]. A galaxy-wide, object-first scatter (`generation/bright_stars.py`, `object_first.py`): per shell chunk one Poisson count from a certified density bound, each candidate gets band, wedge and point, kept with probability true density over bound. Then a per-sector fill (`run_sector.py`, `galaxy/sector.py`): a Poisson count of dimmer systems for the sector, positions sampled uniformly in the sector's volume and rejected against systems already placed there (`_random_position`, Hill sphere separation); the sector fill first builds a system around each pre-placed bright star it contains.
- **Determinism rules already written** [S: sampling-backfill-and-resume.md 2.2]: the block grid must be fixed globally, nothing sequential may cross block boundaries, star parameters come from a stream keyed by (block key, candidate index), and "a region only chooses which blocks to read". That document tested bit-identical lists when processed in any order, split among four workers, or restricted to a sub-region.
- **Cross-sector coupling is small** [S]. Hill-sphere separation is enforced only inside a sector (`_massive_neighbors` reads the sector's own entries and phenomena). The one place a sector reads its neighbours is emission nebulae lit by pre-placed hot stars in neighbouring sectors (`neighbour_hosts`). Names and object-ID collision numbers are ranked by sector address whichever sector saves first (`_rank_chain`, GEN.307), so they are already order independent.
- **Extended objects** (nebulae, remnants, asteroid fields) are stored by their centre's sector and found with the widest placed radius [S]; they are the existing answer to "an object crosses a boundary".
- **Work units are not the sector** [S]: `generate_sectors_unit`, `expand_generate_sectors` and the queue already batch sectors into units, and `galaxy/drill.py` defines nested blocks of 243, 81, 27, 9, 3 and 1 sectors (`galaxy_block_stats`). The execution side of "work units" exists.
- **Charting state exists in pieces** [S]: `sectors` rows exist only for generated sectors, `sector_stats.filled_at`, `generation_runs`, `scatter_epochs`, `region_stats`, `galaxy_block_stats` and `galaxy_block_lights`.

## 3. Determinism: unit of work versus unit of randomness

| | Unit of randomness (the tile) | Unit of work (the region run) |
|---|---|---|
| Fixed by | the galaxy's grid, never by a request | the request (start point, radius, phases) |
| Seed | galaxy seed, version key, tile key; per object the candidate index | none of its own |
| Overlaps | cannot overlap | may overlap freely; a tile already in `chart_ledger` at this epoch is skipped |
| Order | irrelevant | irrelevant |

Rule for the build: **nothing in a draw may read a region parameter, the run id, the arrival order, or whether a neighbour tile exists.** The start point of Boss's "random draw" only picks the first tile; if "starting location for our random draw" meant seeding the stream from it, that has to become "choose the first tile from it".

Things that today depend on arrival and must be keyed before regions can run in any order [S]: the serial in the object ID (rank within the sector, assigned by the one worker that owns the sector, which stays correct while the tile is a whole sector); field-drawn serials for nebulae (already keyed); `id_counters` for run-time births (unaffected, they only grow).

## 4. Region work units: bands within the sphere

Two readings of "bands within the sphere", both fine:

- **A. Bands of the galaxy's own grid.** The region is the set of whole sectors a sphere around the start touches; a work unit is a run of wedges inside one (shell, band) arc, sized so a unit takes a few minutes. This is what `drill.py` blocks and the existing scatter chunks already are. Cost of a sphere in sectors [C]: radius 50 pc touches about 8,200 sectors, 100 pc about 65,000, 200 pc about 523,000 (volume over 64 pc³, before rounding up at the surface).
- **B. Concentric radial bands around the start point.** Easier to explain ("the inner 20 pc first"), gives a natural progress and resume order, but the tiles inside are still the fixed grid sectors, so it is only an ordering.

**Recommendation: A for the units, B as the ordering** (nearest bands first so a user can look while the outer bands run). Width of a unit: one block of 27 sectors (about 1,700 pc³) is a reasonable first guess [I]; size it from measured seconds per sector (PERF items own that number, not re-measured here) so a unit is minutes, not hours, because a unit is also the retry and resume grain (the resumable-run design in sampling-backfill-and-resume.md section 7).

What a region run gains over sector-by-sector, given tiles stay sectors: one density sweep for the whole region (the `_BatchDensity` pass in `run_galaxy.py` already shares it), one pass of the cross-face checks (section 5), one `chart_ledger` write per unit, and one place for the "how far do you want to spread" setting.

## 5. Consistency between adjacent regions

| Thing that crosses a boundary | Today | With region runs |
|---|---|---|
| Systems and their members (AU scale) | belong to the cell of the barycentre; members never leave it in practice [I] | same; ownership by barycentre |
| Hill-sphere separation between systems (about 1 pc scale [I, 4 pc sector, Hill radii of order a pc]) | enforced inside a sector only [S]; across faces not enforced | enforce across faces with a **priority rule**: a candidate is dropped if a higher-priority candidate (lower key) is closer than the required separation, where candidates come from the keyed draw of the neighbour tile, not from what is stored. A pure function, so order independent, at the cost of drawing the one-sector halo of candidates (about 3 x the draws at the surface of a region, none inside it). Optional; the current code does without it |
| Emission nebulae lit by neighbours | `neighbour_hosts` read from `bright_stars` of neighbours | unchanged; the scatter is global |
| Extended objects (nebulae, remnants, asteroid fields) | stored by centre sector, found by widest radius | unchanged; ownership by centre, so a region that charts a sector gets objects whose centre is in it |
| Stars already pre-placed | global scatter | unchanged |
| Rogue planets | drawn per sector | real donors write the rogue into the sector where it lands (section 8) |
| Names and collision numbers | ranked by address | unchanged |

Nothing found in the repo requires reading a neighbour's *accepted* objects, which is what would make order matter. That is the property to protect with a test: charting two adjacent regions in either order gives the same fingerprint (the per-sector digests of GEN.58 already exist for this; a new test charts A then B and B then A and compares the region digest).

## 6. The sector table, the charted flag and the uncharted boundary

**Tiles are whole sectors (recommended).** Then:

- `chart_ledger(tile address, epoch, run_id, charted_at)` is the truth; `sectors` keeps its role (name, centre, edge, stats) and `sector_stats.filled_at` is already this for filled sectors.
- A sector is **charted** if its tile is in the ledger at the current epoch, otherwise **uncharted**. No partial state, so no boundary to store.
- A block of sectors is charted when all its sectors are (`galaxy_block_stats` already counts them). "Fully within a previously charted area" is therefore a count comparison, which is how the Galaxy Map decides today [S: `region_stats`, `galaxy_block_stats`].
- The uncharted sector still lists what is known there: pre-placed bright stars, scattered phenomena, rogues that landed (section 7), all by position query.

**If Boss wants partial sectors** (a region edge cutting a sector), then:

1. The fill must be re-keyed: candidates drawn per fixed sub-tile (for example the 1 pc cube of the 0.25 pc cell grid, or a Morton block), kept or dropped by the priority rule of section 5, serial = (sub-tile key, candidate index) instead of rank. That is the block-first design already specified; it is not what the sector fill does [S].
2. The ID serial rule changes (the 28-bit serial is "rank in sector"; a partial sector needs a rank that does not depend on which part was charted first, which is the keyed index above). Object IDs of everything already generated change at the new epoch (testing mode allows it, no shim).
3. The uncharted part is the sector minus the union of charted sub-tiles. **Spherical coordinates from the centre** work only if the uncharted part is star-shaped about the centre. A sphere clipped to a cube is convex, but the uncharted remainder of a cube around a charted sphere that contains the centre is the cube minus a ball, which is star-shaped about the centre only when the ball is centred there; in general it is not. A stored polar function r(theta, phi) per sector would be wrong for those. The exact record is the sub-tile list, and a polar sketch can be derived from it for display. **Recommendation if partial sectors are built: store sub-tile coverage (a bitmap of 64 sub-tiles of 1 pc, 8 bytes per sector), derive any boundary from it.**

## 7. Every object in the sector where it ends up

Yes, as ownership by position, already how the data is laid out:

- **What belongs to a sector** is the answer to a position query on the sector's bounds, computed when asked. `sector_id` stays only as a cached birth sector for fast joins and for the ID prefix (the ID names the birth place and never changes, object-id-options.md principle 1); it never decides visibility or deletion.
- **Change:** `rogue_planets.sector_id` loses `ON DELETE CASCADE` (and every table like it is audited); deleting or regenerating a sector removes only objects whose *birth* tile is that sector and whose record says generated-by-it, never residents who drifted in.
- **A rogue that ends up in an uncharted sector is a full row** (Boss, 5:23 PM EDT): same table, same columns, same body chain as a rogue in a charted sector (GEN.1112 train 19), stored with its position. The sector gets a `sectors` row on demand with the charted flag off. There are no stubs and no second table.

### 7.2 Generating an uncharted sector that already has residents

- The residents (rogues, pre-placed bright stars, scattered phenomena, extended objects whose centre lies there) are **immutable**: position and attributes are kept.
- The fill treats them as already placed: they count as obstacles in the Hill-sphere placement (`_massive_neighbors` gains them) and are not drawn again. Rogues do not count toward the sector's system count (they are not systems), so the Poisson draw of systems is unchanged by them.
- The sector's own systems are drawn as today; their rogues leave for other sectors (section 8).
- Object-ID collision numbers are ranked by sector address, so a later-saved lower-address sector can bump a stored object's number up (`_rank_chain`). A resident's name can then change once. If "same attributes" must include the name, the rank rule has to stop moving stored objects; flagged in the questions.

### 7.3 Names and IDs may change; every designation is unique forever (Boss, 5:28 PM EDT)

Decision: an object may change its name and its ID, but it stores the history of its previous designations, and its current ID is always completely unique. This reverses the "assigned once, never recomputed" principle in object-id-options.md (principle 1) and makes the collision bump of `_rank_chain` legal (a stored object moving up a number is an ordinary rename).

Design:
- Every object keeps its **internal row id** (stable, never shown). Name and `uid` are attributes.
- A `designations` table, one row per designation an object has ever had: `designation` (the name or the uid in text form, UNIQUE across the whole galaxy), `table`, `object_row_id`, `kind` (name or uid), `valid_from_epoch`, `retired_epoch` (NULL while current), `reason` (collision bump, regeneration, moved, admin rename). A retired designation **stays reserved**, so it is never given to another object: uniqueness holds over current and past designations together.
- Lookup by an old designation (URL, bookmark, API, search) resolves through this table to the current object and returns a redirect, so nothing breaks when a name changes.
- The object page lists the previous designations.
- Fingerprints (GEN.58) digest the current designation only; the history is a separate digest, so a rename does not read as a changed object.
- An object moved by the nightly orbital update keeps its birth-position ID unless Boss wants rogue IDs to follow the current position; if so, the move writes a history row (cost: one row per ID change, moves that change the ID are rare, 1.2 in 10,000 for the cell key [S: region doc 6.6]; the 76-bit position ID changes far less often than the cell).
- Migration: backfill one current row per existing object; no old data to translate (testing mode).

- **Every object in the galaxy, same rule:** stars and phenomena already follow it (scatter tables with `star_system_id`/`built_at` NULL until built). Comets, asteroid fields, interstellar objects follow the same pattern when their items land; no new rule needed.

## 8. Rogues from real donors, stored where they land (GEN.1111, replaces the virtual-donor design)

Rule (Boss, 5:23 PM EDT): a rogue exists only because a generated system ejected it. When a system is built, its ejection events are drawn from its own seed, each rogue is drifted by the closed form (rogue-planet-drift.md section 3), and the rogue is saved at the place it ends up, in that sector, charted or not.

### 8.1 How a donor writes into other sectors

1. The donor's events are a pure function of the donor's seed: ejection time, kick, mass, moons kept, destination position. So the **set of rogues a given system produces is deterministic** and the same whatever else is charted.
2. The destination sector is found by `address_at` (closed form). Its row is created if missing (charted off). The rogue is saved with position, ID from that position, and a link to the donor record (donor system uid, ejection time, speed), so the object page can say where it came from.
3. A donor 7 kpc from its landing is far outside any charted region, so the write goes to isolated rows. Per charted sector about 37 rogues leave at solar density [C: 5.8 per star x 6.4 systems], each into a different sector [drift doc], so a 65,000 sector sphere writes about 2.4 million rogue rows and as many sector rows. Storage, not compute, is the cost: each rogue is a full body with moons (section 8.4).
4. **Idempotence:** regenerating a donor sector must not duplicate its rogues. Each rogue row carries the donor's uid and event index; saving a donor first deletes the rogues it made (by that key) then writes them again. Deleting a sector deletes the rogues its systems ejected, not the ones that landed in it (the cascade change in section 7).
5. A donor's destination can fall in a charted sector. That is the same write; the sector is simply already filled, and the rogue is added to its residents. The rogue count of a charted sector is therefore not fixed at its charting time (it grows as donors are charted). Fingerprints (GEN.58) must either exclude rogues whose donor is outside the sector or digest them separately; otherwise a sector's digest changes when an unrelated region is charted.

### 8.2 What this does to determinism

Each rogue is deterministic. **Which rogues exist in a place depends on which donor sectors are charted**: that is Boss's rule and it is physically honest (we only know the rogues of stars we generated). The test is: chart A then B, or B then A, and the *union* of objects is identical; and a sector's own systems and its own ejected rogues are identical in both. A sector's digest is split into "made here" and "arrived here".

### 8.3 The consequence to know about

Nearly all rogues in a charted sphere would physically have come from stars outside it (the median rogue is 7 kpc from its donor, and only 1.5% land within 100 pc [S: drift doc]). With rogues only from generated systems, a charted sphere of radius 100 pc holds only the rogues its own stars ejected and that did not leave, about 1.5% of its charted stars' yield [C], plus ones from other charted regions. So **the rogue population of a region is a tiny fraction of the physical one until most of the galaxy is charted**, and the rogue count in a sector says more about what is charted than about the sector. This matches the rule as given; if the shown rogue density should look physical, a smooth "expected rogues" number can be shown beside the real count (no objects). Listed as a question.

### 8.4 Cost and the other drift items

- Per donor: the events (cheap) plus one full rogue body per event through the same chain as a planet and moons (GEN.557). The body chain, not the drift, is the cost; a sphere of 8,200 sectors (50 pc) writes about 300,000 rogue bodies [C].
- GEN.1114 (closed-form drift) is needed as designed. GEN.1115 (smoothed density) is no longer needed by the sector draw; it stays useful only for the low-zoom map's expected-count display. GEN.1116 (regenerate a donor on demand) is replaced by a real link to the donor system. GEN.1113 (kick research) is unchanged.
- Kick direction decided (Boss, 5:29 PM EDT): isotropic, and the rogue's velocity is built from the origin planet's actual motion: the donor star's galactic velocity, plus the planet's orbital velocity vector at the moment of ejection (from the stored orbit, not a random draw), plus the scattering kick. The ejection vector therefore follows the planet's own movement; the kick is added to it, not applied perpendicular to galactic motion. Needs the orbit state at the ejection time from the architecture's ejection event (GEN.728 path); the "virtual donors" question is answered by this decision.

### 8.5 Three-stage rogue spread and the Rogues level (Boss, 5:36 PM EDT; replaces the single spread tweak)

Boss's rule, in order, for every sector that is generated:

**Stage 0, residents come off the budget.** Whatever the sector already holds when it is generated (arrived rogues, pre-placed bright stars, scattered phenomena, belts' parents, extended objects whose centre lies there) is subtracted from the sector's calculated mass before anything is drawn, in every phase of sector generation (the mass-backfill rings, the system draw, phenomena, remnants, belts). Today the expected figures live in `sector_stats` (`expected_systems`, `bright_mass_sol`, `total_mass_kg`, `max_mass_kg`) [S]; a resident-mass subtraction is the one new input. The residents themselves are untouched (section 7.2).

**Stage 1, real physics, rogues end up where they end up.** The sector's own systems eject rogues with the donor's real motion (section 8.4 below) and the unscaled closed-form drift (f = 1). Each rogue is saved in the sector where it lands, charted or not (section 8.1). Nothing is tweaked here, so on its own this gives a sparse, physical rogue population (section 8.3).

**Stage 2, the rest of the sector's own build.** Systems, black holes, neutron stars, pulsars, asteroid belts and the rogues that arrived in stage 1 of other sectors' builds are all in place.

**Stage 3, the top-up, runs once per region run, after every sector of the region has finished stages 0 to 2 (Boss, 5:37 PM EDT).** It is the last pass of the run, so the sector's mass shortfall is measured against everything the whole region naturally produced, including all rogues that landed from other sectors of the same run. Compare the sector's expected mass with what it actually holds (residents plus everything drawn). The shortfall is made up by a second round of ejections: **the sector's own systems donate another phase of rogues**, placed inside the sector. Here, and only here, the tweak is used so the rogues land in the sector: the displacement is scaled by the stand-in f so they stay inside (f about 0.014 for a 100 pc scale is no longer the question; inside one 4 pc sector the placement is simply uniform in the sector, with f an implementation detail). Once placed, each rogue is **given a realistic velocity and ejection history as if it came from another sector**: a velocity drawn from the physical distribution (donor-frame isotropic kick added to the galactic motion, as in 8.4) and an apparent origin (distance and direction back along its drifted path, ejection age) that the object page shows as "ejected about N Gyr ago from about D pc away". No donor row is claimed (the apparent origin is not a real system), so there is no fake link: the record says `origin: apparent`, with the sector's own system recorded as the donor of the mass only.

Rules around it:
- **The dropdown shows only the words** None, Few, Normal, Many; the per-star numbers are internal tuning values, never displayed (Boss, 5:39 PM EDT: the levels' numbers are accepted).
- **Cap by the Rogues level (a ceiling, so the top-up never adds too many).** The selector (None, Few, Normal, Many) sets the rogue amount the top-up may reach, so the shortfall never makes more than the level allows: None = 0 rogues, Few 2.5, Normal 4.3, Many 10 terrestrial rogues per star [S: observed 0.5 to 1.4 per pc³, 2 to 10 per star, 5.8 default; Normal is on the low side as Boss asked, Few below the observed range]. The shortfall is the mass the sector is missing, the level is the ceiling. If Boss means the shortfall may exceed the level, say so in the answers.
- **Pass structure.** The region run has a barrier: stages 0 to 2 for all units, then the top-up units. Each top-up unit is a block of sectors (same blocks as section 4), reads the final stored state, and writes its rogues; a sector's `chart_ledger` row gets `topped_up_at`. A later run does not top up a sector again (a sector keeps one top-up); rogues that arrive in it later from newly charted donors simply add to it. If a run is stopped before the barrier, resume restarts at stage 3 only after stages 0 to 2 are complete for the whole set.
- **Order dependence is stated, not hidden.** A sector's top-up depends on what had arrived before the end of its region run. Once generated its content is fixed (stored, residents never move), so a regenerated sector from the same stored state gives the same result, but charting two regions in the other order gives a different top-up. The order test in section 5 compares the union of real-donor objects only; the top-up is checked as "same state in, same objects out".
- **The tweak stays one registered stand-in** (GEN.725, L3, burn-down line), now confined to stage 3 and to the placement; the physics of stage 1 is untouched. Boss wants to play with it: it is a single tuning value `ROGUE_TOPUP_PLACEMENT` plus the Rogues level table, both in `tuning.py`, read from the version key.
- **Where the selector lives:** a row in `galaxy/dials.py` (flag, Generate field, help and the CLI table generated from it) and a sector option that overrides it for that sector, recorded on the sector like `density_multiplier_used` and `override_id` (DB.60) [S].
- **Remainder goes to asteroid fields (Boss, 5:41 PM EDT).** If the cap leaves mass unmatched, the rest is made up as interstellar asteroid fields in the sector, in the same last pass, using the field machinery (interstellar-asteroid-fields.md: ledger, mass cut rule, `<parent ID>.<rank hex>` designation, field radius from age). The fields get the same treatment as the rogues: placed inside the sector, then a realistic origin and velocity as if from elsewhere (apparent origin, marked). Fields retire past 1 ly as in that design. These top-up fields are a distinct, lighter kind (Boss, 5:42 PM EDT): smaller bodies, lower density, spread over a larger volume, and not given a full body history like a field with a real parent. Origin is a stand-in 'unbound remnant' (no collision event claimed), registered in GEN.725 (accepted). Sizes, density and spread are tuning values for the field items to set.
- **Rogue mass mix:** the existing mass bins (96% terrestrial, 4% gas giant) pick each rogue's mass; the rogue chain is the full planet chain (train 19).

What changes in the earlier sections: 8.3's "tiny fraction" holds for stage 1 only; with stage 3 a sector's rogue count reaches its level whatever is charted nearby, so the sector info box can show the count and the real/apparent split.

## 9. Rogue-landing phases

Cut by Boss, 5:23 PM EDT. Nothing is built for them. To fill out the universe, grow the charted sphere by shells (`--grow-shells K`, cost about r squared: roughly 3,100 sectors for the first 4 pc shell around a 50 pc sphere [C]); charting a landed rogue's sector remains possible by hand, the same as any sector.

## 10. Migration (testing mode: no shims, new epoch, wiped database allowed)

1. **Epoch with no generation change** (schema only): drop the cascade, add `chart_ledger`, backfill the ledger from `sector_stats.filled_at`. Alembic migration; update.sh applies it.
2. **Rogue rows written position-owned** with the GEN.1111 train (rogue train 19): `rogue_planets` gets the donor link and event index columns.
3. **If partial sectors are chosen** (section 6): a generation epoch changes the fill to keyed sub-tiles and the serial rule; every generated sector regenerates (wipe and regenerate, allowed). Sized as its own train after 8.1's physics build because it moves every ID.
4. Nothing keeps the old shape: no dual-write, no compat view; the old `sector_id` cascade is removed in step 1.

## 11. Proposed items (working titles; TODO assigns numbers and owners)

| # | Item | Needs | Notes |
|---|---|---|---|
| 1 | `chart_ledger` table, ledger written by sector fill units, charted state derived from it, backfilled from `sector_stats.filled_at` | none | Core 5 (sectors and fill) |
| 2 | Remove `ON DELETE CASCADE` on position-owned objects (rogues first, audit every table with `sector_id`), regeneration deletes by birth tile only | none | Core 5 with the DB lane |
| 3 | Rogue rows written by the donor into the sector where they land (full body, sector row created uncharted, donor link and event index, delete-and-rewrite on regeneration), fingerprint split into made here and arrived here | GEN.1114, 1 and 2 | Planets 1 (GEN.1111 train) |
| 3b | Rogues dial (None, Few, Normal, Many): `dials.py` row, Generate form and CLI, sector override recorded like DB.60, level table, stage 3 top-up with the placement stand-in registered in GEN.725 (L3), resident mass subtracted from the sector budget in every generation phase | 3 | Planets 1 and the dials owner |
| 3c | `designations` table, redirect on old designations, history on the object page, collision bump writes history | none | Core 5 with the DB lane |
| 4 | Region run command: sphere to whole sectors to units of a block, nearest bands first, `--grow-shells K` | 1 | Core 5 and Core 1 |
| 5 | Order-independence test: chart A then B, B then A, compare region digests; and chart a region twice, no change | 1 | Core 1 (reproducible galaxies) |
| 6 | Optional cross-face separation by priority rule (section 5) | 4, a design decision | Core 5, decide after the overlap-free work (GEN.1006) lands |
| 7 | Only if Boss wants partial sectors: keyed sub-tile fill, serial rule change, sub-tile coverage bitmap | 1 to 5 | own train, moves every ID |

## 12. Questions for Boss

1. **Tiles:** keep whole sectors as the unit of randomness and round a region up to whole sectors (recommended, small), or build partial sectors with sub-tile coverage (exact sphere edge, moves every object ID)?
2. **Rogue levels:** accepted by Boss, 5:39 PM EDT; shown to users as words only.
3. **Top-up ceiling:** may the stage 3 shortfall be filled only up to the Rogues level (my reading), or by any amount needed to reach the sector's expected mass?
4. **Bands:** the galaxy's own (shell, band) arcs cut into units (reading A, recommended), or concentric radial bands around the start (reading B, used only as the run order)?
5. **Rogue IDs after a move:** keep the ID from the birth position (recommended) or let it follow the current position, with history?

## 13. Not verified

- Seconds per sector and per unit (PERF items own this; the unit size is a guess).
- Whether every table with `sector_id` can drop its cascade without breaking a test or an admin tool (the audit is item 2).
- Hill separation across faces is currently unenforced; I did not measure how many systems would be dropped or moved by the priority rule.
- The rogue yield (5.8 per star) and the 7 kpc median come from stand-ins and a simulation in the drift study; the phase growth factor is only as good as they are.
- Partial-sector serial design is outlined, not designed.

## 14. Built

**GEN.1128, first slice.** `galaxy/region_run.py` is the plan of a region run, pure (no stream, no database): `plan(center_pc, radius_pc, edge_pc, block_sectors, grow_shells, skip)` takes the sectors a sphere touches (`Region`, reach = radius + `grow_shells` edges + the 0.75 edge allowance, so a sphere that touches a sector holds it) and cuts each (shell, band) arc into units of at most 27 consecutive wedges, nearest first; `skip` takes charted sectors out and splits the arcs they sat in, which is how a stopped run resumes and how an overlapping region costs only its new sectors; `next_units` is the stage barrier of 8.5 (every generate unit first, then the top-up units when a run has a stage 3); `cost` reports sectors, units, shells and the volume estimate (8,181 sectors for 50 pc at a 4 pc edge, the plan holds 10,395 with the allowance, in 530 units; one 4 pc shell around 50 pc adds 2,528). Tests: `tests/test_region_run.py`. The command that runs the units, and the chart ledger write (GEN.1127), follow.

**GEN.1128, second slice.** `planetgen sectors --around X Y Z --radius-pc R [--grow-shells K] [--unit-sectors N]` is the command: a mode of the sector stage (`run_galaxy._region_plan`, `_expand_region`), one unit of work per region unit, nearest first, units re-cut at each stage open less the sectors saved (a rerun resumes). `tests/test_region_mode.py`. The chart ledger write waits for GEN.1127; the stage 3 call rides GEN.1120.
