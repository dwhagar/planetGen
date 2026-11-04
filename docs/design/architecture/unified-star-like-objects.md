# One chain for everything that can be orbited: stars, star-like objects, planets, moons and wanderers

**Status:** plan (decided). Written after 11:20 PM EDT on 2026-10-17 from a read of `main` at bc97f75 (11:14 PM EDT); approved by Boss on cards D1 to D8, taps 11:28 to 11:31 PM EDT on 2026-10-17 (as the coordinator relayed), each the recommended option below. Nothing is built; the items are GEN.550 and its children (GEN.552 to GEN.560, DB.73, UX.138, TEST.239, PERF.136).

**Decided (Boss, cards, 11:28 to 11:31 PM EDT, 2026-10-17):** D1 a shared parent link for every kind, the kind tables stay (option 2); D2 derive the neutron star and black hole counts from the mass function and verified explosion maps, his densities become test targets, and the black hole share needs a verified source first (option 2); D3 rogue moons from both origins, ejection survivors and the isolated planet's disk, by origin row (option 3); D4 a rogue planet's composition from a notional formation host through the bound-planet code (option 1); D5 unbound comets as rows above a size cut with the rest a computed count (option 1); D6 a new part of 8.4 for the chain, GEN.315 stays in 8.1 (option 1); D7 Core 4 for the chain, Core 3 for the pair builder, Core 1 for storage (option 1); D8 rogue planets with moons on by default after a PERF measurement shows the fill stays inside its budget (option 1). Order rule from the coordinator (11:31 PM EDT): GEN.428 follows or joins the wanderer item (GEN.557). The group cap stays 12 (**decided: Boss tapped 'Keep 12' at 12:44 AM EDT on 2026-10-18, as the coordinator relayed; his 'group of up to 7' is read as a typical group, not a limit**) and moons stay the lowest level unless Boss says otherwise (asked at 11:30 PM EDT). **Superseded (Boss, card S6 of docs/design/architecture/unified-system-barycentre.md, 12:10 to 12:13 AM EDT, 2026-10-18, as the coordinator relayed): moons may have satellites through the same hosting gate, so moons are no longer the lowest level; the cap of 12 is unchanged.**

**Boss, 11:18 PM EDT (project thread):** "I want to unify generation of stars and star-like objects (former stars, including stellar black holes and neutron stars as well as brown dwarfs, etc). Then use the same generation chain we use for sector fill that we do with stars. We then use sector fill for unbound comets, rogue planets (which should have the chance to get moons)."

**Boss, 11:20 PM EDT:** "basically I want to unify the process of anything that can have something else orbiting it. Planets have moons at the lowest level, highest level a 7-star-like-object group with each having star systems. Then star system creation is recursive."

Evidence tags: [S] read in the repository; [C] computed here; [R] recalled and not verified (a lane verifies each before use, under the standing source rules). Standing rules kept: nothing fixed or copied from the Solar System or the Milky Way, the DB schema is not a constraint, no backward-compatibility shims, every draw through the one wrapper `util/draw.py` (GEN.56).

**Standing rule (Boss, 5:20 PM EDT on 2026-10-18):** every object with mass follows the same trajectory rules, galactic, planetary and lunar alike, built into the mass object functions; see [mass-object-class.md](mass-object-class.md).

## 1. Summary

1. **The ask in one sentence.** One recursive routine, `populate(host)`, that works for any body or group of bodies that can be orbited: a group of star-like objects, a star-like object, a planet, a rogue planet, a brown dwarf. It draws the host's own properties from one chain, then what orbits it, then calls itself for each thing that orbits. Moons were the lowest level, so the recursion had three working depths: group, star-like object or planet, moon; moons may have satellites through the same hosting gate (Boss, card S6 of unified-system-barycentre.md, 12:10 to 12:13 AM EDT, 2026-10-18), so a moon is a host too.
2. **It is mostly three half-built pieces joined, not new science.**
   - The scatter side already draws a hierarchy of up to twelve star-like members with stability checks (`multistar.draw_member_tree`) [S], but no fill code reads it (the fill builds only one or two stars with the older `StarSystem._companion`) [S]. GEN.264 and GEN.315 are the open items that join them.
   - The orbit tree tables (`orbit_nodes`, `stars.parent_node_id`, `planets.host_node_id`) exist since migration 0104 and nothing writes them [S].
   - Planets and moons already carry `parent_type` and `parent_id` (migration 0102 and 0109) [S], and the 96-bit object ID already counts bodies of one system on a single counter and says the kind is not in the ID, so a planet ejected from its star "keeps its ID" [S].
3. **What does not match the vision** (details in section 2):
   - Neutron stars and black holes never come out of the star chain. They are separate phenomenon classes with their own mass, age and rate, and they have three mass draws and two radius formulas copied into different files [S].
   - A free-floating brown dwarf is a `RoguePlanet` row with a canned composition sentence, an age uniform in 0.1 to 10 Gy that ignores the population, and no companions [S].
   - A rogue planet has `has_moons`, one boolean at 20%, used only to add tidal heating and one sentence. No moon object exists, and the `moons` table cannot hold one (`planet_id` and `star_system_id` are NOT NULL) [S].
   - An interstellar comet is a descriptive row: a diameter, a speed and three names drawn from a fixed list of ten [S].
4. **The recommended shape** (section 3): a registry of host classes, one gated table per class saying what can orbit it and where (the same gated-table pattern as GEN.472), one body draw (mass, age, [Fe/H], state, composition, motion) and one recursive orbiter routine. Remnants become evolved stars with a progenitor, so their counts follow from the same mass function and star formation history as every other star. Wanderers are bodies with no host, drawn by the same routine, so moons for a rogue planet are the ordinary moon process with a different outer limit.
5. **Where it sits.** Almost all of it changes stored values, so it belongs in 8.4 beside the composition chain, with one prerequisite already in 8.1 (GEN.315, Core 3). It should be built before GEN.428 reaches rogue planets, otherwise composition gets written twice (section 5). Ten choices for Boss are in section 6, each with a recommendation for a card.

## 2. How each kind is made and stored in the code

Counts per local 64 pc^3 sector use about 8 to 10 stars [S: `interstellar-object-rates.md`].

| Kind | Made by | Entry to sector fill | Stored in | What it lacks against the vision |
|---|---|---|---|---|
| Main-sequence, giant, white dwarf star | `Star` via `StarSystem`; `stellar_evolution.sample_living_star` draws mass (the primary mass function, Chabrier IMF less companions, GEN.928), age, then `evolve_star` picks the phase [S] | Stage A (pre-placed bright and massive stars, from `bright_stars` rows) and stage D (the cell's other systems) [S] | `stars`, `star_systems`, `bright_stars` [S] | At most two stars from the older `_companion` roll. `initial_mass_sol` is not stored on a built star [S]. |
| Neutron star, stellar and intermediate black hole | `BlackHole`, `NeutronStar` (`compact_remnant.py`, subclasses of `Star`) [S] | Stage B (scatter rows at or above the cut, default 14 Msun for stars), stage G (below the cut, `below_cut_draws`), stage J (ejected ones) [S] | `phenomenon_scatter` while only a point; `black_holes`, `neutron_stars` when built [S] | Mass uniform in 1.1 to 2.2 and 5 to 20 Msun with no progenitor; age uniform in 1e-5 to 10 Gy; rate a fixed density (NS 5.6e-4, BH 7e-5 per pc^3) [S]. |
| Remnant as a member of a multi-star tree | `multistar._remnant` and `_remnant_leaf`: progenitor from the IMF within a range, white dwarf mass from the initial-final relation, NS and BH mass uniform [S] | Scatter time only (`scatter_companions`); the fill does not read it [S] | `scatter_members`, `scatter_nodes` [S] | The progenitor is drawn here and not on the standalone path, so the same remnant has two draw rules. |
| White dwarf | Only as an evolved `Star` (`evolve_star`, line 495) and in trees [S] | Stage D | `stars` | The initial-final relation is applied in four places (`stellar_evolution`, `multistar`, `scatter_companions`, `place_terms`); a hot central star uses a separate Gaussian and a cool one a uniform draw [S]. |
| Brown dwarf, in a tree | `multistar._star_leaf` with `brown_dwarf.state` (GEN.252) [S] | Scatter time | `scatter_members.kind = 'brown-dwarf'` | Not read by the fill. |
| Brown dwarf, free-floating | `RoguePlanet(mass_bin="brown-dwarf")`; mass from `brown_dwarf.sample_mass_solar` [S] | Stage E, a Poisson count of 0.214 per star [S] | `rogue_planets` with `planet_class` NULL | Canned composition sentence; age uniform 0.1 to 10 Gy [S]; the rogue surface module re-types the Burrows luminosity law with a hard-coded 0.05 [S]; no companion (15% of brown dwarfs are pairs, GEN.331 is open). |
| Rogue planet | `RoguePlanet` [S] | Stage E, 5.8 per star over four mass bins, drawn from one mass function ([R3191], 2020; [R2178]) [S] | `rogue_planets` | Terrestrial bulk density uniform in a range, composition a canned sentence per bin, no link to a former host, velocity from the solar circle's speed wherever it sits, `has_moons` a boolean [S]. |
| Interstellar comet | `InterstellarComet` [S] | Stage E at a design rate of 0.007 per pc^3 (0.05 per star) against a measured 7e14 to 1.8e15 per pc^3 ([R2159], Do et al. 2018; corrected from 1e12 on 2026-10-23) shown as a computed count [S] | `interstellar_comets`, `interstellar_comet_composition` | No mass, no size distribution, activity a 50% coin, three components from a list of ten. |
| Planet and moon | `Planet`; moons only for a non-moon planet (`generate_moons`), stable zone out to 0.4895 of the Hill radius (`MOON_PROGRADE_STABLE_HILL_FRACTION`) [S] | Inside `StarSystem` | `planets`, `moons` | A moon never has a moon (nesting is fixed at two levels in the class) [S]. |
| Multi-star hierarchy | `draw_member_tree`, up to 12 members, Mardling-Aarseth, Eggleton-Kiseleva and Kozai screens, `group` nodes for young non-hierarchical systems [S] | Scatter time only | `scatter_nodes`, `scatter_members`; `orbit_nodes` exists and is empty [S] | Not consumed by the fill. |

Other facts the design depends on:

- The fill is a list of stages, each on its own random stream seeded from the sector's stream and the stage name (`fill_stages.py`, GEN.251) [S]. An object built from a scatter row draws under its own seed (`draw.bound(row.seed)`), and the system pass uses `Stream(f"{seed}:system-pass")` [S]. So a stage or a kind can be added without moving any other draw. This is the property the recursion needs.
- The system pass already is a small version of "anything that can hold something": `SYSTEM_PASS_ODDS` says a millisecond pulsar holds planets at 0.7% and a rogue planet a moon at 20% [S].
- The registry (`phenomena/registry.py`, GEN.250) has a `system_model` column with five fixed values (`none`, `star-system`, `central-star`, `rogue-moon`, `remnant-planets`) [S]. Those five values are what the host model generalises.
- The count of a non-star kind is a Poisson draw with mean `rate per star x realised star count x kind factor` (`run_sector.py:222`) [S]. The scatter and below-cut paths instead use the expected star density, on purpose, so their totals match the plan [S].
- A rogue planet or brown dwarf is placed at a uniform point of the cell. It may land inside a star's Hill sphere and is not counted among the massive neighbours [S].
- Galactic orbit fields: stars and remnants pass their galactocentric distance, rogue planets and interstellar comets call the function with none, so they get the solar circle's speed and period anywhere. GEN.142 (peculiar velocity for rogue planets) is open for the same gap [S].
- The default phenomenon mass cut in code is 14 solar masses (`tuning.PHENOMENON_MIN_MASS_SOLAR`); several sections of `phenomenon-scatter-mass-cut.md` still print 20 and 8 and are marked superseded there [S].

## 3. The unified chain

### 3.1 Vocabulary

- **Body:** anything with a mass that can be orbited or can orbit: a star-like object (protostar, pre-main-sequence star, main-sequence star, giant, white dwarf, neutron star, black hole, brown dwarf), a planet, a moon, a small body (comet, asteroid, belt object).
- **Node:** the barycentre of two children (the existing `close`, `wide`, `compact` and `group` nodes, with up to 12 members as decided earlier; Boss's "7" is read as "the largest group the generator makes", and the cap of 12 stands: Boss tapped 'Keep 12' at 12:44 AM EDT on 2026-10-18, as the coordinator relayed).
- **Host:** a body or a node that has orbiters. A host is either top-level (no parent: a wanderer or the head of a system) or it is itself an orbiter of another host.
- **Level:** group, member, planet, moon. A moon was the lowest level by Boss's earlier rule; superseded (Boss, card S6 of unified-system-barycentre.md, 12:10 to 12:13 AM EDT, 2026-10-18): moons may have satellites through the same hosting gate, so a moon can be a host.

### 3.2 One body draw

Every body, whatever its class, goes through the same five steps, each on a stream derived from the body's seed (`Stream(f"{seed}:birth")` and so on), so any body can be rebuilt alone from its seed and its host's state:

1. **Origin.** A gated table (the GEN.472 format, so rows have classes, weights and a source state) says how a body of this class came about: "formed in a system", "formed in isolation", "ejected from a system", "remnant of a star". The weights of the rows come from the cited rates.
2. **Birth.** Birth mass (the IMF, continued below 0.08 Msun by the brown-dwarf mass function; the planet mass function for planets), birth time from the population's star formation history (GEN.469, merged), [Fe/H] from the population and place (GEN.158, merged).
3. **State at the present.** Evolution to the age: pre-main-sequence, main sequence, giant branch and remnant for star-like objects (`evolve_star`, the neutron-star radius law, `brown_dwarf.state`); cooling and surface state for brown dwarfs and rogue planets (`rogue_surface`, kept, its duplicated Burrows constants replaced by the `brown_dwarf` module).
4. **Composition.** Drawn first, as decided on 1:21 AM EDT (composition model section 4a): stellar and remnant composition from GEN.439 and GEN.440, bulk composition of planets, moons and small bodies from GEN.428, GEN.442 and GEN.466. A wanderer needs the same inputs a bound body has (a formation distance, a snow line, a star's abundances), so it is given a notional formation host (section 3.5).
5. **Motion.** Position in the cell and velocity from the population's dispersion, plus the kick or ejection speed for remnants and ejected bodies (section 3.4).

What this removes: four copies of the white dwarf initial-final relation, three of the neutron star and black hole mass draw, two Schwarzschild formulas, three age rules, the separate canned composition sentences, and the rogue-only galactic orbit shortcut [S].

### 3.3 Host classes and what can orbit them

One registry row per host class, extending `phenomena/registry.py`. The row's `system_model` column becomes a pointer to an orbiter table, one gated table per class, in the GEN.472 format:

| Host class | Companion table (nodes) | Satellite table (zone rule) |
|---|---|---|
| Group node | members, 2 to 12, the existing count tables and shape draw | none (a group is not a body) |
| Star-like object | the pair rule of `multistar` (extended in 3.4) | planets, belts, comets: from the inner limit (engulfment, Roche) to the outer limit from the node tree (Holman-Wiegert for a binary, the Hill sphere in the galaxy for a single star) [S: both exist] |
| Remnant | as a star-like object, with the compact filters (Peters time, direct collapse, kick survival) [S] | millisecond pulsar planets 0.7%, white dwarf survivors outside 1.5 AU [S] |
| Brown dwarf | its 15% pairs and wide companions (GEN.331) | disk-limited zone; planets and moons, a few percent [R] |
| Planet (bound or free) | planet pairs (GEN.332) | moons: Roche limit to 0.4895 of the Hill radius, prograde [S]; for a free planet the limit is set by its history (3.5) |
| Moon | none | satellites through the same hosting gate (Boss, card S6, 12:10 to 12:13 AM EDT, 2026-10-18; it was "lowest level") |

The routine itself is short:

```
populate(host, depth):
    if host has companion table and depth allows: draw nodes; each member is a new host
    for each satellite row that passes the gate: draw orbiters inside the zone rule; each orbiter is a new host
    call populate(orbiter, depth + 1) for every orbiter with a table
```

Every orbiter is built from `Stream(f"{host_seed}:orbit:{slot}")`, so building a system from stored members (the GEN.264 and GEN.315 requirement: "the stars of a built system equal the stored members") is the same call with a different starting point, and a fill in any order gives the same system.

A moon of a free planet and a moon of a planet in a system are the same code. This is what makes "rogue planets with the chance of moons" a table row and not a feature.

### 3.4 Pairs of any mass

`multistar._Builder` already draws periods, eccentricities and orientations from a law per pair kind and checks stability (60 tries per pair) [S]. Two changes make it host-neutral:

- A pair kind for planet with planet (double planets, GEN.332's "rogue planet pairs") and brown dwarf with brown dwarf or planetary-mass companion (GEN.331), with period and mass-ratio laws taken from sources the lanes verify.
- A planet with a moon is a hierarchical pair whose mass ratio is extreme; the moon table of 3.3 is the same orbit draw in a smaller zone, so the stability test (Hill fraction) can be shared with the node builder.

Compact remnants stay in the tree builder with their filters. Nothing about the tree builder's output format changes, so the scatter-time companions (Boss, 5:51 PM EDT on 2026-10-15: companions of stars and star-like objects are made at scatter time, planets, moons and belts at sector fill) stand.

### 3.5 Remnants as evolved stars

A neutron star or black hole becomes a star-like object whose state is "remnant", with a progenitor mass and a birth time like every other star. The consequences:

- **Counts become derived.** A quick Kroupa integral using the generator's own breaks and slopes (`IMF_BREAKS_SOL = (0.08, 0.5, 150)`, slopes 1.3 and 2.3) gives the share of stars born above 0.08 Msun that have progenitors in 8 to 20 Msun: 0.45%, in 8 to 30 Msun: 0.53%, and in 20 to 100 Msun: 0.17% [C]. The fixed densities in the generator are 0.4% of stars for neutron stars and 0.05% for black holes [S]. The neutron star number is already the right size from the mass function alone; the black hole number needs a failed-supernova and direct-collapse share and a kick-escape share (about 40% of neutron stars and 2% of black holes leave the galaxy [S: `compact-remnant-regions.md`]). Boss's research rates then become test targets, as already done for the brown dwarfs in GEN.252, and not inputs. The progenitor ranges used by `multistar` (8 to 30, 20 to 100) are [S] but their sources are not named in the file; a lane must verify the explosion-model sources ([R1466, R3327, R4251] [R]) before they are inputs.
- **Space comes from a kick step.** `remnant_distribution.py` keeps its role (thicker layer, black hole core excess, the Lorimer pulsar profile) as the placement of the kick, so a remnant is born in the young population and moved [S].
- **The scatter and the mass cut stay.** A remnant above the cut is still a scatter row (kind, mass, seed); one below is still drawn at fill from the cell's stream (GEN.168 pattern). Only the draw function behind both changes: they call the body draw with the origin row "remnant of a star". Intermediate and supermassive black holes keep their own origin rows (cluster runaway growth, nuclear), since they do not come from the stellar mass function.
- **Planetary nebula central stars** are drawn once, from the star chain, and the Gaussian and uniform side paths for hot and cool white dwarfs go away.

### 3.6 Wanderers through sector fill

Stage E becomes "wanderers": the bodies whose origin row says "ejected from a system" or "formed in isolation", counted per star as in the code (Poisson, mean `rate per star x realised star count x factor`, anchored on the cited rates, which stay: 5.8 rogue planets, 0.214 brown dwarfs). What changes is what a wanderer is.

**Rogue planets.**

- **Origin.** Two rows in the origin table: ejected from a system (most terrestrial and sub-Neptune masses, because disk scattering ejects small bodies while giants stay bound [S: `interstellar-object-rates.md`]) and formed in isolation (planetary-mass objects continuing the star-like mass function). Weights by mass bin; where only a range is known the Dirichlet draw applies (decided 1:43 AM EDT).
- **Formation host.** An ejected planet is given a notional host to feed the composition draw: a host mass from the IMF at its population, a formation distance from the planet formation model, a snow line, and the host's abundances from the population's [Fe/H]. Nothing of the host is stored beyond the draw's inputs on the planet's own seed; when the orbital update unbinds a planet from a stored system (GEN.291), the real host plays this part, and the row simply changes its parent (the ID already allows it).
- **Age and motion.** The age is the population's, from the same star formation history as stars (not uniform 0.1 to 10 Gy); the velocity is the population's dispersion plus an ejection speed. This closes GEN.142 and the galactic orbit shortcut.
- **Moons.** The moon routine of 3.3 with the planet as host. For a bound planet the outer limit is a fraction of the Hill radius at its orbit; for an ejected one, only moons that were inside a smaller fraction of the Hill radius at ejection survive ([R1032, R1998] [R], to verify), so the draw takes the notional formation distance and gives a smaller zone, hence fewer and tighter moons. For a planet formed in isolation the zone comes from its circumplanetary disk (size limit from the planet's own rotation and tidal radius [R]). The fixed 20% in `SYSTEM_PASS_ODDS["rogue-planet"]` goes; the chance of at least one moon becomes an output, with 20% as a test target only until a source replaces it.
- **Pairs.** Planet pairs come from the pair builder (3.4), with their own rate (GEN.332).
- **Cost.** About 58 rogues per local sector (the code's own docstring figure) [S]; at 20% that is about 12 moon systems per sector. Rogue planets outnumber stars about 5 to 1, so this will be the main cost of the fill; no measurement was made here, so the cost is unknown and has to be measured with the PERF harness before it is switched on.

**Brown dwarfs.** Free-floating ones are star-like objects with no host and take the star chain's age, population and [Fe/H], plus the companion and planet tables of 3.3. They leave `rogue_planets`.

**Unbound comets.**

- **Class.** A small body (`comet`), origin "ejected from a system" (most of them) or "formed in isolation" (none expected, a table row with weight 0 until a source says otherwise).
- **Counts.** The measured density is about 7e14 to 1.8e15 per pc^3 ([R2159], Do et al. 2018, read by the bulk composition lane on 2026-10-23 in reports/small-body-size-laws-2026-10-23.md; this replaces the earlier 1e12, range 1e11 to 1e14), so most cannot be rows. Use the same device as the remnant mass cut: a lowest size, rows above it drawn at fill from the cell's stream, and the count below it stays a computed figure on the sector page (`interstellar_debris_count` exists [S]). The size-frequency law and the cut are a research item (candidates: the inferences from 'Oumuamua and Borisov in [R2159], which `tuning.py` already cites [S], and Do, [R1124] [R]). The design rate of 0.007 per pc^3 is a placeholder labelled as such until then.
- **Properties.** Mass and size from the size law, composition from the comet table (GEN.495, GEN.442; the within-system shape is a stand-in and the two interstellar objects bracket the offset, per decision D6 of the tunable gaps note [S]), activity derived from composition, size and distance to the nearest star (replacing the 50% coin), velocity from the population plus ejection speed. A comet is a body, so a star's captured comet or a comet bound to a wanderer is the same class with a parent.

### 3.7 Storage

"The DB schema is not a constraint", so three choices exist (section 6, card D1). In every case:

- The object ID is unchanged: body number 0 for the head of a system or a wanderer, 1 and up for the rest on one counter [S]. A planet leaving its star keeps its ID and changes its parent.
- `parent_type` and `parent_id` (on `planets` and `moons` since 0102 and 0109 [S]) widen to every orbiter and add `node` (using `orbit_nodes`, which GEN.315 begins to fill) and the free case (no parent).
- Free planets, brown dwarfs and comets need `moons` that point at them, so `moons.star_system_id` and `planet_id` stop being NOT NULL or the rogue gets a system id of its own (the ID scheme already treats it as a top-level object).
- The build that turns a scatter row into a system (GEN.315) deletes `scatter_members` and `scatter_nodes` in its own transaction; the same rule applies to every wanderer's stored seed.

## 4. Sector fill, as it would run

The stage list of `fill_stages.py` keeps its letters; the content changes.

| Stage | In the code | Unified |
|---|---|---|
| A pre-placed | bright and massive stars become systems | the same, built with `populate` from the stored members and nodes |
| B scattered | scatter rows become bare objects | rows become bodies by the body draw (remnants by progenitor) |
| C system pass | planets for ms pulsars only, from `SYSTEM_PASS_ODDS` | removed: `populate` of B and G objects does it, from each class's table |
| D systems | the cell's other systems, 1 to 2 stars | the same, with up to 12 members from the tree builder |
| E per-star kinds | rogue planets, brown dwarfs, interstellar comets, each a bare row | wanderers: bodies with no host, through `populate` |
| F nebulae | unchanged | unchanged (their required stars use the same star chain) |
| G below cut | NS and BH from independent draws | the same body draw as B |
| J ejecta | unchanged | unchanged (a kick strips planets, as stated in the catalogue) |
| H flags | unchanged | unchanged |

Each stage keeps its own stream, so adding wanderers' moons moves no draw of any other stage.

## 5. What it touches in TODO.md and docs/plan

Counted from `docs/TODO.md` on main; owner lines are those the items carry [S].

| Item | Relation | Note |
|---|---|---|
| GEN.315 (and GEN.264, Core 3) | **prerequisite, unchanged** | Fill builds from stored members; `StarSystem(member_tree=...)`. Ask Core 3 to keep its interface host-neutral (a function that takes a tree and returns bodies), at no extra cost. |
| GEN.331, GEN.332, GEN.327, GEN.273, GEN.269 to 271, GEN.129 (Core 3 train) | **keep, become rows** | Substellar companions, circumbinary planets, rogue planet pairs, group nodes: all are rows of the tables in 3.3 and builders of 3.4. GEN.332's "rogue planet pairs" is the pair builder with planets. |
| GEN.252, GEN.541 (merged), GEN.549 (Core 4, open) | **keep** | The brown dwarf state function is the brown-dwarf branch of step 3 of the body draw. GEN.549 changes only its cooling track. |
| GEN.439, GEN.440 (Blocker 1 then Core 4) | **consumed** | Stellar and remnant composition (including brown dwarfs) is step 4 for star-like objects. No conflict. |
| GEN.428, GEN.429, GEN.442, GEN.458, GEN.466, GEN.457, GEN.495, DB.69, DB.70 (Core 4, held for 8.4) | **conflict if built first** | They write composition for `Planet` objects. A `RoguePlanet` is not one, so rogue composition would be written a second time. Do the wanderer unification (item U8) before or together with GEN.428, or the rogue draw is built twice. |
| GEN.142 (rogue planet peculiar velocity) | **superseded by U8** | The wanderer's motion is a step of the body draw. Fold it in. |
| GEN.110, GEN.291, GEN.346 (collisions, unbinding) | **simpler** | A planet leaving its star changes its parent and keeps its row, ID and moons. GEN.291 said it becomes a "rogue-planet row"; with one body model there is no second table. |
| GEN.130, GEN.317, GEN.318, GEN.253, GEN.263, GEN.377 | **keep, one input changes** | Compact remnant systems, pulsar planets, white dwarf rows, scatter pairs and limits stay. U5 changes how a remnant's mass and age are drawn, not what systems it gets. |
| GEN.250, GEN.251 (built) | **extended** | The registry's `system_model` column is replaced by the orbiter table pointer; stage C is absorbed. |
| GEN.447 to 450, GEN.463, GEN.467, GEN.469 (merged or in 8.1) | **inputs** | Binary law, derived radius and lifetime, star formation history: the body draw reads them. |
| GEN.517, GEN.520, GEN.522 (tuned draw builders B9 brown dwarfs, B12 comet sampler, B14 moons) | **related** | Their tuned rows are inputs to the brown-dwarf, comet and moon tables; B12's two interstellar objects bracket is the offset in 3.6. |
| GEN.301 (headline 8.1) | **not a prerequisite; do not add one** | Only GEN.315 is in 8.1. The rest is 8.4. |
| MAP.82 to 84, UX.24, UX.25 | **touch** | Rogue planet display and the "Mark Rogue Planets" setting follow the new rows; map symbols for moons of a rogue are new. |
| ADM.62, ADM.63 | **touch** | Aspect catalogues: one catalogue for bodies and hosts is smaller than two. |
| TEST.159, TEST.160, TEST.162, TEST.211 | **extend** | Representation audit counts per kind and per level; multiplicity audit for planet and moon levels; star composition sweep includes the brown dwarf and remnant branches. |

Lanes (from memory of the board, 11:12 PM EDT as the lanes reported): Core 4 owns the composition and stellar chain, so it owns the body draw and the orbiter routine (U4 to U9); Core 3 owns the scatters and the multi-star train, so it owns U6 and the interface of GEN.315; Core 1 owns the identity and repair chain, so it owns the storage items (U11). Core 4's 8.4 items are held until its 8.1 chain is built (Boss, 10:17 AM EDT on 2026-10-17), which fits the order below.

Epoch and migrations: every item that changes a drawn value bumps the generator epoch (OPS.56 rule; epoch 21 is the last stamped [S]). The full wiped rescatter that already collects the epoch-changing items would be joined by the 8.4 ones, or this chain gets its own, since all of it follows the 8.1 stamp. Migration head is 0109 [S]; new revisions are numbered contiguously at merge. No compatibility code: rows in `rogue_planets` move to the new rows by regeneration, not conversion.

## 6. Choices for Boss (one card each)

**D1. How far the storage is unified.**
1. Code only: one builder and one in-memory tree, tables unchanged. Cheapest; every pair of tables stays and the rogue moon needs a second moon table.
2. Code plus one shared parent link for every kind (widen `parent_type`/`parent_id`, fill `orbit_nodes`), kind tables keep their own physical columns. **Recommended:** it gives the "anything can orbit anything" queries and the ID rule without moving 100-column tables.
3. One physical `bodies` table for everything. Cleanest queries, largest migration, and the highest risk for the pages, maps and the orbital update.

**D2. Remnant counts.**
1. Keep the fixed densities (NS 5.6e-4, BH 7e-5 per pc^3) as inputs.
2. Derive them from the mass function, the star formation history and verified explosion maps, with the densities as test targets. **Recommended:** the neutron star count already matches [C]; the black hole needs a verified failed-supernova share; it matches the rule that fixed content gives way to derivation.
3. Derive neutron stars and keep black holes fixed until their source is verified.

**D3. Where rogue moons come from.**
1. Survivors of ejection only (moons inside a small fraction of the original Hill radius).
2. Formed with an isolated planet only (circumplanetary disk).
3. Both, as rows of the origin table by mass bin. **Recommended:** ejected rogues and isolated planetary-mass objects are different populations in the cited rates; both mechanisms need sources checked first (section 3.6).

**D4. Notional formation host for a rogue planet's composition.**
1. Draw a notional host (mass from the IMF, a formation distance, a snow line) and use the bound-planet composition code. **Recommended:** one composition code for all planets.
2. Use the population's [Fe/H] and one global formation-distance distribution, no host.
3. Keep composition descriptive for rogues until 9.0.

**D5. Unbound comet rows.**
1. Rows above a size cut with a research-backed size law; the rest stays a computed count. **Recommended:** the same pattern as the remnant mass cut that Boss approved on 2026-10-09 (default 14 Msun in the code).
2. Keep the 0.007 per pc^3 design rate.
3. Raise the design rate and show every comet (not possible at 7e14 to 1.8e15 per pc^3, listed for completeness).

**D6. Version.**
1. A new phase of 8.4, after the composition fractions and before the surfaces; GEN.315 stays in 8.1. **Recommended.**
2. Do the structural parts (U2, U3, U6) in 8.1 phase 2 and the drawn-value parts in 8.4. Spreads the work but adds to the 8.1 count (686 open at 3:05 PM EDT on 2026-10-17).
3. After 9.0.

**D7. Lanes.**
1. Core 4 for the chain, Core 3 for the pair builder, Core 1 for storage. **Recommended:** each already owns the neighbouring code.
2. One new lane for the whole chain. More coordination, slower start.
3. Core 3 for everything.

**D8. Whether the rogue cost is switched on by default.**
1. On, after a PERF measurement shows the fill stays inside the budget. **Recommended.**
2. Off by default with a plan-configuration switch, like the halo and satellite switches (GEN.537).
3. On without a measurement.

Two points that I read as decided and did not make cards for: the group cap stays 12 (decided earlier and confirmed by Boss's 'Keep 12' card, 12:44 AM EDT on 2026-10-18; he wrote 7 as a typical group), and moons were the lowest level (Boss's words), so moons of moons were not generated; superseded (Boss, card S6 of unified-system-barycentre.md, 12:10 to 12:13 AM EDT, 2026-10-18): moons may have satellites through the same hosting gate, and the sub-moon stability question is open (unified-system-barycentre.md section 4.5).

## 7. Proposed items, in order

The TODO thread numbers them and assigns lanes. "8.4" means the composition and new-worlds version; none of these except U1 is needed by 8.1.

| # | Item | Needs | Owner (proposed) | Version |
|---|---|---|---|---|
| U1 | Approve this note and the cards of section 6 | none | Boss | before any build |
| U2 | Host registry: replace `system_model` by orbiter-table pointers; the orbiter tables as gated tables (GEN.472 format) with their sources and states; no behaviour change | U1, GEN.472 (built) | Core 4 | 8.1 phase 3 (back-scatter, Boss's card 11:38 to 11:40 PM EDT on 2026-10-17; filed as GEN.552) |
| U3 | GEN.315 as filed, with a host-neutral interface | none | Core 3 | 8.1 (unchanged) |
| U4 | One star-like body draw: origin, birth, state, [Fe/H]; the initial-final relation, remnant mass, radius and age rules in one place; the free-floating brown dwarf moves in | U1, GEN.252 (merged), GEN.439, GEN.440 | Core 4 | 8.4 phase 1 |
| U5 | Remnant counts derived from the mass function and star formation history, kick step on `remnant_distribution`, Boss's densities as tests; verified progenitor maps first | U4, GEN.463, GEN.467, GEN.469 (merged), D2 | Core 4 | 8.4 phase 1 |
| U6 | The pair builder for any masses: planet pairs, brown dwarf pairs, planet with moon share the stability test | U3, GEN.331, GEN.332 | Core 3 | 8.4 phase 1 |
| U7 | `populate(host)`: planets, belts, comets and moons from the class tables, replacing `StarSystem._generate_planets` host logic and `generate_moons`; stage C absorbed | U2, U3, U6 | Core 4 | 8.4 phase 2 |
| U8 | Wanderers: origin table, notional formation host, population age and motion, velocity (closes GEN.142), composition by the bound-planet code; rogue planet becomes a `Planet` with no parent | U4, GEN.428, GEN.458, D3, D4 | Core 4 | 8.4 phase 2 (with or before GEN.428) |
| U9 | Moons of wanderers (ejection survival and isolated-disk rows), the 20% boolean and `SYSTEM_PASS_ODDS["rogue-planet"]` removed | U7, U8, GEN.466 | Core 4 | 8.4 phase 2 |
| U10 | Unbound comets: size law and cut, composition from the comet table, derived activity and velocity | U8, GEN.495, GEN.442, D5 | Core 4 | 8.4 phase 2 |
| U11 | Storage: shared parent link, free parents, `orbit_nodes` filled, rogue moons storable, `rogue_planets` and `interstellar_comets` rows regenerated | D1 | Core 1 | 8.4 phase 1, with U7 |
| U12 | Pages, class pages, map symbols and descriptions for the new rows and moons of wanderers; GEN.446 text | U9, U10, MAP.82 to 84 | Interface | 8.4 phase 2 |
| U13 | Tests: determinism by host (a system rebuilt from seeds and stored members equals the stored one), representation counts per kind and level, rogue and comet rates against the cited ones, derived remnant shares against Boss's densities | U5, U9, U10 | Core 4 with TEST | 8.4 |
| U14 | PERF measurement of the fill with wanderers and moons, and the decision of D8 | U9, U10 | Core 5 | 8.4 |
| U15 | Research list for the lanes (each source verified before use): ejection survival of moons, circumplanetary disks of isolated planetary-mass objects, failed-supernova and direct-collapse shares, size-frequency law of interstellar objects | none | read-only research lane | 8.1 phase 1 (no behaviour change) |

Order: U1, then U15 and U3 in parallel with U2; then U4, U6, U11; then U5; then U7 and U8 together; then U9, U10; then U12, U13, U14.

## 8. What this does not do

- It does not change what a user sees in 8.1, and adds nothing to the 8.1 headline chain except GEN.315's interface note.
- It does not generate planets or moons at scatter time: companions of stars and star-like objects stay at scatter time, planets, moons and belts at fill (Boss, 5:51 PM EDT on 2026-10-15).
- It did not make moons of moons; superseded (Boss, card S6 of unified-system-barycentre.md, 12:10 to 12:13 AM EDT, 2026-10-18): moons may have satellites through the same hosting gate.
- It does not fix any number: every rate and law above is either a cited input (already in the repository), a value to derive, or a labelled placeholder until a verified source replaces it.
- It makes no claim about cost; the PERF item U14 exists because rogue planets with moons may be the heaviest new part of the fill and nothing here measured it.

## 9. Sources

In the repository, read for this note: `generation/object_first.py` (the scatter's object-first pass), `generation/run_sector.py`, `fill_stages.py`, `system_pass.py`, `multistar.py`, `scatter_companions.py`, `phenomena/registry.py`, `compact_remnant.py`, `rogue.py`, `builders.py`, `physics/brown_dwarf.py`, `physics/stellar_evolution.py`, `galaxy/object_uid.py`, `db/store.py`; docs `phenomena-catalogue-and-fill.md`, `interstellar-object-rates.md`, `compact-remnant-regions.md`, `phenomenon-scatter-mass-cut.md`, `brown-dwarf-state.md`, `composition-model.md`, `multistar-and-compact-systems.md`, `database-schema.md`. Recalled and still to verify: [R1032] (ejected planets with moons), [R1998] (exomoon dynamics in scattering), [R1466, R3327, R4251] (explosion and remnant maps), Do, [R1124] (interstellar object density). Already cited in the repository: [R3191] and 2020, [R0606, R1990, R2159, R2178, R2374]. The Kroupa share figures in 3.5 are a quick computation with the generator's own IMF breaks and slopes, not a sourced result, and ignore binary evolution and kicks.
