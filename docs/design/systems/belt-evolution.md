# Evolving asteroid belts, and a ledger of every body in a system

Status: plan (Boss said yes 4:27 PM EDT 2026-10-30; nothing built). Written 4:45 PM EDT 2026-10-30 as the lane reported. No generator code changed. Asked for by Boss (4:19 PM EDT 2026-10-30): evolve belts "so we know exactly what is in them, are they planets", and start an idea to "catalogue every single object in a star system individually with its own point of mass data". Companion: `reports/belt-root-cause-2026-10-30.md` (outside the repo) (what exists). Related: GEN.1014 (clearance ledger covers belts), GEN.1006, `systems/weird-worlds.md`, `docs/design/composition/small-bodies-composition.md` sections 13 and 14, `small-body-size-laws.md`.

Evidence tags: [S] standard result, number as published as far as recalled; [R] recalled, **no paper opened this session** (web research was not run), so the building lane must read the source before the number is used; [REPO] read in the repo; [C] computed from repo values. Anything still [R] when built is registered in GEN.725 as a stand-in (T2 or T3, L0 to L1).

## 1. The idea in one paragraph

A belt is not a ring with a word; it is the **surviving tail of a planetesimal population** that the system's planets, star and age have acted on. Generate it as: (1) what formed, from the disc solids over its annulus; (2) a size distribution with a known top end; (3) the history that removes mass (collisional grinding, giant-planet clearing and resonances, early instability, drift, sweeping, star changes); (4) a **ledger**: the largest bodies kept as individual objects with their own mass, the rest as a population with a size law and a total mass. "Is it a planet?" then falls out of the IAU test applied to each ledger body.

## 2. Physics to build on

### 2.1 What forms and where
* Solids over the annulus: `planet_mass.solids_surface_density_gcm2` integrated (`belt_mass.primordial_mass_earth`) [REPO].
* Belts sit where growth stalled: inside the snow line at the giant's inner edge (a gap left by an oligarch chain whose next mass would exceed the isolation mass), at the snow line (pressure bump, Ida and Guillot 2016 [R]), and beyond the last giant (the outer, Kuiper-type disc, Gladman and Volk 2021 review [R]). The architecture (GEN.728) already holds the chain; a belt is a **gap or an unfinished tail** in it, not a 10 percent roll.
* Occurrence of cold debris discs by star type: Sibthorpe 2018 [REPO, GEN.787]. Keep it as the occurrence for the outer disc; the warm belt is the architecture's.

### 2.2 Collisional grinding (built, not wired)
Wyatt et al. 2007 [S, READ in repo]: `M(t) = M0 / (1 + t / tc)`, late mass independent of M0, `M_max ~ rho r^3.5 (dr/r) Dc / (M*^0.5 t)`. `belt_mass_earth` implements it. The size distribution in equilibrium is Dohnanyi `q = 11/6` cumulative for the strength-dominated tail [S], steeper (q about 4) for the gravity-dominated bodies above ~100 km (Ivezic 2001 [REPO table]). Grinding needs `Dc`, the largest body, which the ledger supplies.

### 2.3 Giants clear gaps
* Mean-motion resonances: `a_res = a_giant (q/p)^(2/3)` for p:q. For the Solar System the strong ones are 4:1, 3:1 (2.50 AU), 5:2 (2.82), 7:3 (2.96), 2:1 (3.27) [S]. Small bodies at those semi-major axes have chaotic eccentricity and are lost in about 1 to 10 Myr [R, Wisdom 1982; Moons and Morbidelli 1995].
* Chaotic zone around a planet: half-width about `1.3 a_p (m_p / M*)^(2/7)` [S, Wisdom 1980] (the repo's mutual Hill spacing is the same family). A belt inside the chaotic zone of a giant is removed, a belt beyond it is stirred (more eccentricity, so less mass by Wyatt's `e` term).
* The giants' history matters more than their present orbits: Grand Tack and Nice-type instabilities depleted the Solar System's belt by about 99 percent [R, Walsh 2011; Raymond and Izidoro 2017; Morbidelli 2010]. `planet_architecture` already moves and scatters giants (hot Jupiters clear the inner chain); the same event must **sweep or excite** a belt it crosses. Hot Jupiter or migration through a belt: belt gone or scattered into a resonant fringe.

### 2.4 Drift: Yarkovsky and YORP
* Yarkovsky: thermal force changes `a` at `da/dt` of order 1e-4 AU/Myr for a 1 km body at 2.5 AU, **proportional to 1/D** [S, Bottke 2006 review]. Consequence: small bodies drift into resonances and are removed; bodies above ~40 km (belt ages of Gyr) do not move. This is the mechanism that keeps resonance gaps empty and turns a size distribution into a **size-dependent** radial distribution: the inner and outer edges get a size cut. Scales with `L*^-1`-type factors and `rho` [R]; a luminous star (A, giant) drives it hard.
* YORP: spin-up/down of small bodies over `~1e7 to 1e8 yr` for 1 to 10 km [R], leading to rubble-pile disruption and binaries (small-body-size-laws section 2 holds binary fractions). Needed only for the **spin and binary fraction of members**, not for belt mass.
* Combined: a **resonance-feeding loss rate** per size bin, `dN/dt = - N(<D) / t_Yark(D)`, applied to bodies below a cut; the cut grows with age as `D_cut ~ t^-1`... (to be derived from the da/dt law and the fractional resonance cross-section by the building lane).

### 2.5 Star and system events
* Main-sequence luminosity rise moves the snow line outward by `L^0.5` (the repo has `snow_line_au`); the composition of a belt is set at **formation**, so belts keep their formation snow-line position (store it), do not recompute from the present star.
* Giant branch: mass loss widens orbits by `M_form / M_now` (built for white dwarfs, `widen_for_mass_loss`). Add: an AGB/RGB star's luminosity drives Yarkovsky/YORP hard and ablates small bodies and melts ice [R]; a belt inside the engulfment or sublimation radius is lost partly, not whole (fix the inner-edge-only rule: keep the surviving outer part). Belts beyond a few AU survive and are **destabilised** by the sudden `M` drop: the resonance map moves, a share is scattered inward. The scattered share is what pollutes a white dwarf: 25 to 50 percent of white dwarfs show metals from accreted rock [S, Zuckerman 2003, 2010; Koester 2014], with accretion rates of 1e8 g/s typical [R]. The repo already has the `white-dwarf-pollution` row (multistar.py); link it to **whether a belt or planet survives** in the ledger.
* Binaries: a companion truncates a belt at about 1/3 of the periastron (Holman and Wiegert 1999 [REPO]); Kozai excitation of a belt by a distant companion; a circumbinary belt starts beyond the stability limit. Close-pair stars that merge or go through common envelope clear the belt (GEN.318 already removes planets there).
* Young star: debris-disc phase is stirred by self-stirring of Pluto-size bodies (Kenyon and Bromley 2004, 2008 [R]): the **number of dwarf planets** a young belt holds follows from the size law's top end.

## 3. What a belt is made of: the ledger

Three tiers per belt (and per outer disc). Sizes are in diameters; the cut is a **mass cut**, so it is the same rule for every system.

| Tier | Contents | Stored as | ID |
|---|---|---|---|
| L1 individual | every body with mass above `M_L1` (stand-in 1e19 kg, about 200 km rock; the whole Solar System's belt has ~30 above this, KBOs a few hundred [R]) | one row each: mass, radius, density, composition (from the belt makeup), orbit (a, e, i), spin, status | object ID, kind = small body, child of the belt |
| L2 reproducible | bodies from `M_L2` (stand-in 1e15 kg, ~10 km) to `M_L1`: drawn **on demand** from the belt seed and the size law, identical every time, never stored | none | serial under the belt, valid while the belt's version holds |
| L3 population | everything below: total mass, size-law slope and cut, radial profile, binary fraction | belt record | none (a belt is its own object) |

Rule: mass is conserved. `M_belt = sum(L1) + integral(L2) + L3`, and the belt's stored mass (Wyatt) is the budget the draw spends. The **point-mass data** Boss asked for is the L1 row plus the belt itself as one extended mass (position barycentre, total mass) for an N-body or barycentre calculation (`unified-system-barycentre.md`).

### 3.1 The L1 draw
1. Belt mass `M_belt(t)` from `belt_mass_earth` with `Dc` fixed by step 2.
2. Number above mass `m`: `N(>m) = N0 m^(-b)` from the size law (main belt cumulative q 3.0 from 5 to 40 km, steeper above (stand-in BR20), Ivezic 2001 [REPO]; Kuiper hot/cold Fraser 2014 [REPO]; built as `physics/belt_ledger.py`, GEN.1061). The top body is the one with `N(>m) = 1` (order statistic), jittered; the Solar System's Ceres (9.4e20 kg) is a test target, not an input.
3. Resolve `Dc` and the mass: iterate once (Dc changes `M_belt`, which changes the top end), as `belt_mass.py` notes.
4. Sample L1 bodies from the order statistics; assign orbits inside the belt, **excluding resonant gaps** and size-dependent Yarkovsky depletion; assign eccentricity and inclination from the stirring (`e = (upper - lower) / (upper + lower)` already used).
5. Compositions: the belt's taxon shares at the body's own `a` (the kernel makes a gradient, not one belt value), differentiated bodies (V, A, M) only if the system's 26Al ratio says so (GEN.473 step 4). Density from `compaction` (existing taper 1e19 to 1e21 kg).

### 3.2 Is it a planet? (status per body)
IAU 2006: a planet orbits the star, is in hydrostatic equilibrium, and has cleared its orbit. Apply as three computed flags, stored on every L1 body and every planet:
* **Round (hydrostatic).** Mass above the equilibrium mass for its composition: ice bodies roughly above 1e20 to 1e21 kg, rock above ~1e21 to 1e22 kg [R, Bierson and Nimmo 2019 for ice; Thomas 2010 for rock]. The repo's compaction taper (1e19 to 1e21 kg, GEN.442) is the same scale: reuse it. The existing `DWARF_PLANET_MAX_MASS_KG = 5e22` (a stand-in) is the upper bound of "dwarf".
* **Cleared.** Stern and Levison 2002 `Lambda = k m^2 / a^(3/2)` (>1 clears) [S]; Soter 2006 `mu = m / sum(other masses in its zone)` (Earth 1.7e6, Ceres 0.33, Pluto 0.08) [S]; Margot 2015 `Pi = k m / (M*^(5/2) a^(9/8)) >= 1`, `k = 807` [S]. Use Soter's mu (needs the ledger: it exists once every body has a mass) and Margot's Pi as the cross-check. This is what the ledger buys: the question "did it clear its orbit" has a number.
* **Orbit class.** star-orbiting (planet or minor planet), moon (host planet), or binary partner.
Status word: `planet`, `dwarf planet` (round, not cleared), `minor body`, `moon`. A "planet" that fails Soter's mu is a dwarf planet; a Mars-size body in a chain is cleared. **Decided by Boss (section 11): store both** the IAU word and the numbers.

## 4. Evolving a belt (the algorithm)

Order, per belt, once when the system is built (and again by the epoch step when a star changes):
1. **Birth**: annulus from the architecture's gap or tail; `M0` from disc solids; formation snow line stored.
2. **Giant acts**: list the giants' final and historic orbits; remove the chaotic zone of each giant at its **formation** and migration track; excite the rest (`e` up, mass down); mark resonance gaps for the present giants.
3. **Early depletion**: if the system had a giant instability (architecture event), multiply `M0` by a drawn survival (stand-in lognormal about 0.01 to 0.3, widest range until read: Walsh 2011 [R], Raymond and Izidoro 2017 [R]); record the fact.
4. **Grind**: `belt_mass_earth(M0, ..., age, Dc)`; top end from the order statistic.
5. **Drift**: apply Yarkovsky size-cut to the radial profile and L2/L3 slope; record the cut diameter.
6. **Star events**: engulf/sublimate part (keep survivors), mass-loss widening, scattering share (to pollution rate), companion truncation. Each event is a function on the belt and appends a line to its history (a text history the page can show: "the outer third was scattered when the star lost 46 percent of its mass").
7. **Ledger**: fix L1 bodies, compute status flags and mu.
Because grinding is `1/t`, a stored `age` is enough to reproduce any time: bring-to-epoch (GEN.1006) recomputes mass with the new age without redrawing members whose mass is above the grinding scale (the top bodies do not change; the tail does).

## 5. The ledger of the whole system (the idea Boss asked for)

* **Everything has a mass.** Star(s), planets, moons, rings (mass given by existing `draw_ring`), L1 bodies of every belt, comets above `M_L1` (a comet nucleus mass exists via `draw_comet`), the belt L3 tails as extended masses. A system view is then `sum(m_i) = M_total`, an N-body-ready table, and a **barycentre** that includes the belts (the Sun's barycentre wobble from the belts is a test).
* **One list, one rule.** A flat system table (`system_bodies`: id, parent id, kind, mass, radius, density, a, e, i, status, tier) with `kind` in star, planet, moon, belt-member, comet, ring, extended. Existing tables (`planets`, `moons`, `asteroid_belts`, `comets`) stay the typed store; the ledger is a view or a generated index (decision for the DB lane: view vs table; recommended: view first, no new storage).
* **Mass cut is the one rule** for what is its own object; the IDs of L1 bodies are stable (derived from belt seed and rank). L2 bodies are addressable by `(belt id, rank)`, rebuilt, never stored; a **generate anywhere** (8.9) path asks the belt for rank k.
* **Where it plugs in**: Mass-Object class (`mass_object.py`) already family-bases `SmallBody`; L1 members are `SmallBody`s with a parent link. GEN.1014 (clearance ledger) already needs "which radii are occupied": the belt's gap map (resonance gaps, chaotic zones) is exactly a clearance entry; reuse it so planets are not placed into a belt's gap and belts into a planet's.
* **Map and pages**: system map draws the belt as an annulus plus L1 dots at true orbits; a body page for each L1 object (status, mass, composition, history).

## 6. Weird-world hooks (for the weird-worlds table, kept separate)
Belt oddities this model makes possible without painting: a **heavy belt that never accreted** (a failed planet: a 5 Earth-mass belt at a gap, rare); **rubble belt of a destroyed planet** (collisional belt with one family composition, differentiated shards; Theia-like impactor remnants); **white dwarf with an intact outer belt** and ongoing pollution; **double belt with a giant between** (Epsilon Eridani-type, Su 2005 [R]); **Trojan swarms** (Jupiter Trojans outnumber the main belt above 10 km [R]); **dwarf-planet-dominated belt** (a belt whose mass is one Ceres, as ours is 1/3 Ceres); **circumbinary belt** beyond the stability limit; **resonance-locked belt** (Plutinos, Neptune 3:2); **belt around a pulsar** (the repo's pulsar planets already built from a disrupted disc). Feed the table with the trigger inputs listed above; no extra generator work.

## 7. Stand-ins this design needs (register in GEN.725)
Stand-in values to register when built: `M_L1` and `M_L2` cuts; the early-depletion survival lognormal; the Yarkovsky size-cut coefficient (until the da/dt law is read and calibrated, tier T2, tweak L1); the instability depletion fraction; the dwarf-planet roundness masses by composition; the Soter mu thresholds if a word cut is used. All carry tier and tweak level; none reads the Solar System as an input; the Solar System belt, Kuiper belt and the Ceres/Pluto/Eris masses are **test targets**.

## 8. What to read before building (research list)
Read in GEN.1060 (2026-10-30): see `belt-sources.md` for what was opened, what was not reachable and the numbers (registry rows BR1 to BR18, section 6.95). The [R] tags below and above are covered by those rows.

Wyatt 2007 (READ), Ivezic 2001 (READ), Fraser 2014 (READ), Singer 2019 (READ), DEEP V (READ). To read: Bottke 2006 (Yarkovsky/YORP review), Wisdom 1980/1982 and Moons-Morbidelli 1995 (resonance and chaotic zone), Walsh 2011 and Raymond-Izidoro 2017 (belt depletion), Stern and Levison 2002, Soter 2006, Margot 2015 (planet definition numbers), Bierson and Nimmo 2019 (ice hydrostatic mass, READ for density), Thomas 2010 (rock equilibrium), Kenyon and Bromley 2008 (self-stirring), Debes and Sigurdsson 2002, Zuckerman 2010, Koester 2014 (white dwarf pollution), Ida and Guillot 2016 (snow line pressure bump), Gladman and Volk 2021 (Kuiper belt review), Holman and Wiegert 1999 (READ in repo). Boss cannot access some papers: use second-hand if a source is closed.

## 9. Identity of belt members (Boss, 4:21 PM EDT 2026-10-30)

Asteroids get **no name**, only a designation: their object ID printed the way unnamed stellar objects print it (`FE81000A2B-0000005-000`: sector, serial, body, hex, `galaxy/object_uid.format_id`). The designation is a pure function of the ID, so one ID always gives the same text. Names need not be unique outside a system, which fits: the body number is only unique inside its system. Rules that follow [REPO: object-ids.md, object_uid.py]:
* No name draw for belt members. No salad name, no codec name (`naming_key` renames only IDs of the kinds it lists; belt members are left out of it).
* **Body number must be derived, not issued in insertion order.** Today a system's stars, planets, moons, belts and comets take body numbers 1, 2, 3 in insertion order (`store._UidIssuer`). Members of a regenerated or evolved belt must keep their number, so a member's body number is a function of `(belt body number, rank by mass, seed)`, never of save order.
* **Capacity problem.** The body field is **12 bits (4095 per system)**. A Kuiper-type belt with a hundred thousand L1-sized bodies is far over it, and even the tier cut at 1e19 kg may give thousands round a rich star. Options: (a) keep L1 small enough (a count cap per system, with the cut raised so no system exceeds a few hundred); (b) a member sub-number inside the belt's own body number (the belt's ID plus a rank field printed as an extra group, `...-000.0042`); (c) widen the body field in a new layout (costly: the layout is fixed at 96 bits). Recommended: **(b)**, because it keeps the system ID layout and gives every member a stable designation `<belt ID>.<rank hex>`; the evolved belt keeps ranks stable by sorting on a mass drawn once at birth. Decided (section 11).
* A body that a player or admin edits (8.6) is still an ID-designated member; only planets and moons take names.
* L2 bodies (drawn on demand) take the same designation with rank above the L1 count; they are addressable but not stored.

## 10. TODO items

Lock means it touches a lock (epoch, migration, golden, API version or shared generator path); the rest is separate-branch work that waits for no one. Owners and trains are the best guess when filed (4:35 PM EDT 2026-10-30) and can move; the TODO file and docs/plan/lanes.md hold the current ones.

| Item | Draft no. | What | Owner (guess) | When |
|---|---|---|---|---|
| GEN.1060 | 1 | Research read of the sources of section 8 | Sources 2 | none (starts now) |
| GEN.1061 | 2 | belt_ledger.py | Sources 2 | after 1 |
| GEN.1062 | 3 | resonance_gaps.py | Core 3 | after 1 |
| GEN.1063 | 4 | body_status.py | Aux 2 | after 1 |
| GEN.1064 | 5 | belt_history.py | Star 2 (or Sources 2) | after 1 |
| GEN.1065 | 6 | Member designation function | Core 1 | none |
| GEN.1066 | 7 | System body table or view and barycentre design | Planets 1 | after GEN.1014 wiring |
| GEN.1067 | 8 | Per-member composition from the taxon kernel | Sources 1 | after 2 and GEN.473 wiring |
| GEN.1068 | 9 | Belt birth from the architecture gap or tail (lock) | Planets 1 | train 15 |
| GEN.1069 | 10 | Store belt mass, Dc, history, snow line (lock, migration) | Sources 2 | train 15 |
| GEN.1070 | 11 | Store members, map dots, body pages, API bump (lock) | Core 2 with Sources 2 | train 15 |
| GEN.1071 | 12 | Outer Kuiper-type disc as a belt (lock) | Sources 2 | train 16 |
| GEN.1072 | 13 | White dwarf and giant-branch events (lock) | Star 2 | train 16 |
| GEN.1073 | 14 | Bring-to-epoch for belts (lock) | Bug 1 | train 16 after the GEN.1006 chain |
| GEN.1074 | 15 | Captured asteroid moons: research, then a build plan | research thread | now |
| GEN.1093 | belt_density.py: objects per square AU above a stated diameter, and mass per square AU | Aux 7 | after 2 |
| UX.240 | Belt page data builder | Core 2 | after 2 and GEN.1093 |
| UX.241 | Asteroid belt information page, route and API (lock) | Core 2 | train 16 |

Order: 1, then 2 to 6 in parallel, then 9 to 11 as train 15; 12 to 14 as train 16. Item 15 (captured asteroid moons, Boss: "Jupiter has so many moons because a lot of them are former asteroids from the belt, if gravitational bodies are near a belt, they should pick up asteroidal moons") is a research item beside the belt chain: a belt feeds the moons of nearby giants, so a captured moon leaves the belt's mass in the ledger.

## 11. Boss decisions (4:24 PM EDT 2026-10-30)
1. **Member numbering: decided, sub-number inside the belt's ID** (option b of section 9). Members print as `<belt ID>.<rank hex>`; the rank is stable because it sorts on a mass drawn once at birth.
2. **Planet words: decided, store both** the IAU word (planet, dwarf planet, minor body) and the numbers (roundness mass, Soter mu, Margot Pi, Stern-Levison Lambda), so either can be used.
3. **Mass cut for individual objects: decided, the rule (option E), 4:27 PM EDT 2026-10-30.** Store every body above 1e19 kg; if that is fewer than the bodies needed to cover 90 percent of the belt mass, add the next heaviest up to 300 per belt. The values are stand-ins registered in GEN.725. The options considered are kept below.

### 11.1 Options for the individual-object cut (all counts [R]/[C], order of magnitude, for a Sun-like star)
Diameter is for rock at 2,500 kg/m3 (m = rho pi D^3 / 6). The Solar System rows are the only worked check; they are test targets, not inputs. "Main-belt-like" is a belt of about 0.5 thousandths of an Earth mass; "Kuiper-like" is a disc of about 0.02 Earth masses [R].

| Option | Cut | About | Objects, main-belt-like | Objects, Kuiper-like | Belt mass catalogued | Storage and cost |
|---|---|---|---|---|---|---|
| A | 1e20 kg | 430 km | about 3 to 5 | tens | about 55 to 65 percent | a few rows per belt, negligible |
| B | 1e19 kg | 200 km | about 10 to 30 | hundreds to a few thousand | about 70 to 80 percent | under 1 MB per rich system; fast |
| C | 1e17 kg | 40 km | about 1,000 to 2,000 | tens of thousands | about 85 to 90 percent | tens of MB for a rich system; draw time seconds; the system page needs paging |
| D | 1e15 kg | 9 km | about 5,000 to 10,000 | hundreds of thousands | about 95 percent | too heavy to store per system; only as on-demand L2 |
| E | rule, not a number: catalogue from the top until 90 percent of the belt mass is covered, cap 300 per belt | varies | 10 to 300 | 300 (capped) | 90 percent where the cap allows | bounded: at most 300 rows per belt; a thin belt gets few, a heavy one never explodes |

Each row below the cut is still drawn on demand (tier L2) and is addressable, so nothing is lost by choosing a higher cut: the cut decides what is stored and listed, not what exists.

**Chosen: E with a floor of 1e19 kg (B) as the bound for the unit of "individual".** In words: store every body above 1e19 kg, and if that is fewer than the bodies needed to reach 90 percent of the belt's mass, keep adding the next heaviest up to 300. The Solar System belt falls in 10 to 30 objects (Ceres, Vesta, Pallas, Hygiea and neighbours) which a reader can name; a rich outer disc is bounded at 300; mass is always accounted for.

## 12. The asteroid belt information page (Boss, 4:37 PM EDT 2026-10-30)

A page per belt: (1) the individually tracked bodies by designation (no names), each linked to its member page; (2) the contents and sizes of the rest of the belt mass as diameter bins (count, mass share, typical composition) that sum with the listed bodies to the stored belt mass; (3) the density factor, objects per square astronomical unit above a stated diameter (1, 10 and 100 km), plus mass per square AU. The unit is registered in the unit registry so no number shows without it. Items: GEN.1093, UX.240, UX.241. Scatter placement: the scatter sweep runs after the giants are placed and before the belts are finalised (belt-origins.md).
