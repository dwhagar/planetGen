# Systems that cannot overlap by construction: one source of truth for every star and body (GEN.1006)

Boss (2026-10-30, 10:51 AM EDT, in the missing-third-star thread; asked again 10:53 AM EDT): for a multi-star system the star objects need to be aware of each other's rebuilding, so that in theory any number of stars or objects in one system never overlap by design, and any fitting and adjusting comes after. Also a test that measures how often items are redrawn to make them fit (TEST.476).

Informs: GEN.1006 (this design), GEN.1004 (the quick fix, replaced by this), TEST.476, TEST.162, TEST.430, GEN.725.

Status: design, revised 2026-10-30 after Boss asked (10:59 AM EDT) for the design to be checked against real data and papers, not memory. Nothing here is built. Evidence tags: **[S]** seen in this repository's code; **[C]** computed in this session by running the generator (scripts named in section 14); **[P]** read in the paper itself in this session (arXiv text); **[P2]** read in a second paper that quotes the source (the source is not open to this session); **[R]** recalled and not checked (these remain, listed in section 14, and stay below Verified in GEN.725); **[I]** inferred, with the step shown. A figure keeps the tag of its weakest input. Times to Boss are EDT only.

## Summary

- **The bug is two sources of truth, and I have now found the main difference between them [C].** `multistar.draw_member_tree` evolves each companion to the system's age and draws the orbits from those radii and masses (`_star_leaf`, `_Builder.clears`). The system build then makes each `Star` again from `(seed, initial_mass, age)`. The tree evolves the companions at [Fe/H] = 0, because `system._drawn_plan` makes the `Placed` without the primary's [Fe/H] (`Placed.tree_fe_h` falls back to 0.0), while the build gives each companion the primary's real [Fe/H]. Passing the primary's [Fe/H] into the tree (a one-line change, tried in a scratch edit and reverted) cut the first-draw failures of forced triples from 24 to 7 in 2,500, of quadruples from 17 to 8, of five-star systems from 50 to 22 and of six-star systems from 95 to 38. The rest come from a second, smaller difference: with the same [Fe/H] the built radius still differs from the tree's by up to about 4% (median under 1%; 400 stars compared), enough to break a pair drawn at the Roche floor. GEN.1004 added a second check after the build (`multistar.built_tree_clears`) and redraws the whole tree when it fails [S]. That finds the overlap; it does not remove its cause, and it cannot tell the tree about a rebuild before the orbits are placed.
- **Recommendation: evolve every member once, freeze the result as an immutable record, and make everything else read the record.** The tree draws orbits from the records; the build makes each `Star` from its record without evolving it again; planets, moons and belts take their room from a clearance ledger computed from the same records. Overlap then cannot arise from disagreement, because there is nothing left to disagree.
- **Draw orbits by truncated sampling, not by draw-and-redraw.** Each pair's period is drawn from its law cut to the range its members allow (Roche lobe, pericentre, hierarchical stability), by inverting the cut distribution. A pair with no allowed range is a tree that cannot exist, and the shape step handles that once, not by retry. Every retry loop in the stellar part then has a counted rate (TEST.476) that should read zero for overlap.
- **Add a history rule the code does not have.** A star that is a white dwarf today was a giant, maybe an AGB star, in the past. Its partner must have been outside that giant's Roche lobe at the separation of the time, or the pair went through a common envelope and is a post-common-envelope pair. The code handles this for white dwarfs only (`classify_pair`, `mass_loss.engulfment_radius_au`); the rule makes it general, for any pair and any depth of tree [I].
- **Prove it with an independent verifier, not with the generator's own screens.** A pure function reads the stored rows of a finished system (stars, nodes, planets, moons, belts) and checks every pair of bodies for overlap, Roche fill, Hill room and hierarchical stability, written without importing the generator's placement code. The test runs it over thousands of seeds with forced triples, quadruples and deeper trees and requires zero violations.
- **The quick fix goes away when this lands.** `built_tree_clears` and the redraw loop in `system.py` (`_companions_fit`, `_built_stars_clear`) are deleted, with no shim (Boss's standing rule).

## 1. What the code does now (the facts the design starts from)

All **[S]** unless marked.

1. **The tree draws leaves already evolved.** `_star_leaf(initial_mass, age_yr, rng, ...)` takes a birth mass to its state at the system's age: a brown dwarf from `brown_dwarf.state`, a young star from `pre_main_sequence.state`, otherwise `evolution.evolve_star(mass, age, rng, fe_h)`. If that returns None (collapsed) and `lighten` is on, the mass is cut by 20% and tried again. The leaf keeps `kind`, `mass`, `initial` and `radius_au`. The `rng` is passed into `evolve_star`, so the evolved state can depend on a draw.
2. **The orbits are drawn top down against those leaves** (`_Builder`). A pair of two single members gets a Roche floor on its axis (`axis_floor`, `roche.separation_for_filling` at `MAX_FILL_FACTOR` 1.2); any other pair gets an extent rule (`EXTENT_AXIS_FACTOR` 2.5 times the summed radii, pericentre at least `PERICENTRE_EXTENT_FACTOR` 2.0 times the extent). A pair inside a larger one is bounded by `stability.max_inner_axis` and then checked by `stability.is_triple_stable` (Mardling and Aarseth, Eggleton and Kiseleva criteria) and the Kozai-Lidov screen. Up to 60 tries a pair, `TREE_RESTARTS` restarts a tree, and then the lightest loose member is dropped.
3. **The build does not take the leaf's state.** `stored_trees.PlannedStar` carries `seed`, `initial_mass_solar` and `mass_solar`; the build constructs `Star` from the birth mass and the age. `MASS_TOLERANCE` (1e-6) guards only the mass: a pair that the evolution changed (transfer, merger) "does not build from a birth mass and an age" and falls back. The radius is not guarded.
4. **Members the build cannot make are cut out of the tree after it is drawn** (`without_members`: brown dwarfs, neutron stars, black holes in a system build). The cut keeps every other pair's drawn orbit, and re-computes the period of the nodes above from the mass that is left. The cut changes a node's mass and therefore the stability margin of every pair above it, and nothing re-screens them [S, I].
5. **The only protection after the build is GEN.1004's** `built_tree_clears`, which re-runs the pericentre and Roche screen on the built radii and masses, and `_companions_fit` redraws the whole companion set when it fails. The missing-third-star thread measured about 1 in 150 triples and 1 in 50 quadruples before the fix (seeds 241 and 297 reproduce a triple; one quadruple had a star 0.033 AU across with a partner at 0.0088 AU). I re-ran the draw on 2,500 seeds each (1000 to 3499) with the fix in, counting first-draw failures of `built_tree_clears`: 24 triples (0.96%), 17 quadruples (0.68%), 50 five-star (2.0%) and 95 six-star (3.8%) systems [C]. The numbers differ from the thread's because the fix redraws until a draw clears and I counted each first failure; the pattern that matters is that the rate rises with star count. The 0.033 AU giant of the quadruple appears again in my failures (a 1.125 solar-mass star of 0.0332 AU radius in a pair with a 0.6 solar-mass white dwarf) [C].
6. **Planets already have their own, separate, stability pieces:** `_orbit_floor_au` (Holman and Wiegert 1999 circumbinary limit for a close pair), a per-star S-type limit (`stellar_object.py`, `wide_binary.py`), `mutual_hill_radius_m` spacing between planets (`system.py _mutual_min_distance_au`), `min_orbit_distance` of `5 * hill_radius` (`planets.py`). They are computed after the stars are built, from the built stars, and do not share any structure with the tree.

So there are three separate "where is there room" computations (tree orbits, built-star check, planet floors) fed by two evolutions of the same star. The design merges them.

## 2. Why overlap happens, and the failure classes a design has to close

| Class | Cause | Today |
|---|---|---|
| A. Evolution mismatch | The tree and the build evolve the same member separately and get different radius or kind | GEN.1004 redraw |
| B. Cut-out mismatch | `without_members` removes a member and the masses above change | not screened |
| C. History mismatch | A pair that today is detached would have been inside the giant's lobe in the past | white dwarfs only |
| D. Depth mismatch | The check for a pair of two stars differs from the check for a pair containing a sub-tree (Roche versus extent rule) | two rules |
| E. Bodies after stars | Planets, moons and belts are placed against the built stars, with no tie to the tree's orbits other than a floor | separate floors |
| F. Mass changes after placing | Adiabatic mass loss widens orbits as `a (M_star + M_p)` constant (Jeans 1924; Hadjidemetriou 1963; Veras et al. 2011) [R]; the code already uses this for white-dwarf planets (`mass_loss.expansion_factor`) | planets only |

A design is complete when each class has a by-construction answer and a test.

## 3. Options compared

| Option | Idea | For | Against |
|---|---|---|---|
| 1. Keep the redraw (GEN.1004) | Build, check, redraw | Done, small, uses the real Roche physics | Finds a symptom, not the cause; the measured first-draw failure rate climbs with star count, 1.0% (3 stars), 0.7% (4), 2.0% (5), 3.8% (6) of 2,500 each [C]; biases which hierarchies survive; cannot reach class B, C, E; a retry loop is invisible until counted |
| 2. Build first, orbits second, one record | Evolve each member once into a frozen record; tree and build both read it | Removes class A at the source; no new physics; the build becomes simpler (it stops evolving); makes a stored tree rebuild exactly | Needs `Star` to accept a ready state (a constructor path, not a new class); merger and stripped members (no birth-mass recipe) need their own record fields |
| 3. Constraint pass over the whole tree | After drawing, solve for the periods that satisfy all constraints at once (a projection onto the feasible set) | Handles any depth and any mix; fits afterwards, as Boss asked | Needs a solver and a rule for which constraint yields; changes the drawn distributions in a way that is hard to state; tests are harder |
| 4. N-body check | Integrate every system for some Myr and keep survivors | The most direct stability evidence | Far too slow for millions of systems; needs a new dependency; only good as an offline audit |
| 5. Option 2 plus truncated sampling plus a ledger plus a verifier (recommended) | Frozen records; draw each period from its law cut to the allowed range; planets draw from allowed intervals of a ledger; a separate verifier proves it | Closes A to F; one measured retry rate; same code path for 2 or 12 stars; planets, moons and belts share the same structure | Largest change; output changes, so it rides trains in small parts |

Option 3 is what Boss's "fit and adjust after" would be, but it is the second step: once the records and ranges are shared, nothing needs fitting. Option 4 is kept as an optional offline audit only (section 9).

## 4. The design

### 4.1 The member record (one source of truth)

One function evolves a birth mass to a state at an age and returns a frozen `MemberRecord`. It is called once per member, by the roster step, and its result travels with the member everywhere (tree, stored rows, build).

Fields (all in the unit registry, GEN.568): `kind`, `state` (Yerkes or spectral type), `mass_solar` (present), `initial_mass_solar`, `radius_au` (present), `radius_max_au` (the largest radius over the past: the red-giant and AGB tips for a star that is now smaller, from the same evolution call), `mass_max_solar` (the largest mass it ever had, for the history rule), `teff_k`, `luminosity_sol`, `age_yr` (effective age for a merger product), `fe_h`, and the `seed` that produced it.

Rules:

- The roster draws everything random once (including whatever `evolve_star` takes the `rng` for). The record is then deterministic data.
- `Star` gets one new constructor path, `Star.from_record(record, ...)`. It does not evolve. A test requires that for every member kind the built star's radius, mass and Teff equal the record's within a stated tolerance (this replaces `MASS_TOLERANCE` as the guard, and it covers the radius).
- A stored tree stores the record fields it needs (most already live in `scatter_members`: mass, initial mass, radius, state, flags). Adding `radius_max_au` and `mass_max_solar` is the only schema change; the migration is a testing-mode epoch.
- Members the build cannot make as a `Star` (brown dwarf, neutron star, black hole) keep a record too, and are no longer cut out silently: see 4.4.

This closes class A. It is also the cheapest part and the one that makes the rest testable.

### 4.2 Orbits by truncated sampling

The top-down orbit draw already knows, for each pair, the closest it may be (Roche floor or extent floor) and the widest (`max_inner_axis` from the stability criterion). Today it draws a period from the law and retries when the draw falls outside. The design draws from the law **conditioned on the allowed range**:

1. Convert the allowed axis range `[a_min, a_max]` to a period range with the pair's total mass.
2. Take the law's cumulative distribution (log-normal, or log-uniform for the compact laws) at both ends, draw a uniform between them, and invert. One draw, no retry. This is the exact conditional distribution of the law restricted to that range.
3. Eccentricity next: the thermal law is cut so the pericentre clears the floor (`e <= 1 - a_floor/a`), drawn the same way (thermal `f(e) = 2e` has a closed-form inverse on a cut range), and capped by the existing caps.
4. If `a_min >= a_max` for a pair, no orbit exists for this shape and these members. That is decided in the shape step (try another shape, or drop the lightest as now), once per tree, and counted as its own named event. It is the one place a tree can fail; it is no longer a loop around the draw.

Is the bias acceptable? The observed period laws (Raghavan 2010; Moe and Di Stefano 2017 [R]) are laws of systems that exist today, so they already exclude pairs that merged or were unstable; cutting a draw to the allowed range is the same exclusion, done at the pair instead of the tree. The old whole-tree redraw conditioned on the joint outcome; the truncation conditions each pair on its own limit given the pairs above it. The difference is small and TEST.476 (section 8) measures how much probability is moved, which is the honest answer; if it matters, the shape weights absorb it.

The Kozai-Lidov and gravitational-wave screens stay as screens (they depend on the joint of inclination, eccentricity and mass and have no clean closed-form range). They should be counted by TEST.476 and expected to be rare after the truncation because the stability range already removed most close calls.

### 4.3 The history rule (class C, F)

For a pair of members `i` and `j` at separation `a` now, each with present mass `m` and largest mass `m_max`:

- Adiabatic mass loss gave a separation of the past `a_then = a * (m_i + m_j) / (m_i_max + m_j_max)` (a larger orbit grew, so the past orbit was smaller). Veras et al. (2011) treat the Jeans and Hadjidemetriou adiabatic case and show where it fails (fast loss and strong eccentricity, arXiv 1107.1239 [P]); for a binary the adiabatic form needs the mass loss to be slow against the orbital period, which holds for a giant branch loss but not for a supernova [I]. The same rule is already in `mass_loss.expansion_factor` for planets.
- The giant of the past had radius `radius_max_au`. If `radius_max_au` exceeds its Roche lobe at `a_then` (Eggleton 1983 [R]; `roche.roche_lobe_radius_au` already implements this) the pair went through mass transfer or a common envelope in the past. It is only allowed as a **post-common-envelope pair**, whose present axis comes from the law for that case (the existing `post-common-envelope` interaction and the 'stripped' laws), not from the general period law.
- Otherwise it is an ordinary detached pair and the draw above applies.

For a subtree the same test is applied to the subtree's outer orbit against the largest extent of the whole subtree in the past (its apocentre then, plus its members' `radius_max_au`).

This replaces the white-dwarf-only engulfment test with a rule that covers a stripped helium core, a neutron star's progenitor, and the third star in a triple. It also gives a principled value to a rule the code now holds as a lookup (`WD_ENGULFMENT_BASE_AU`). Mustill and Villaver (2012) show the engulfment limit is set by the star's maximum radius plus tidal decay (larger for Jovian planets, small for terrestrial ones), with the expansion from mass loss pushing some inner planets out to safety [P]; so the limit should come from `radius_max_au` times a tidal factor, with the lookup as the stand-in until then. Same principle for stars: the existing `GIANT_ENGULFMENT_RADIUS_FACTOR` of 2.0 times the radius is the same idea [S].

**Built (GEN.1012, first part).** `physics/pair_history.py` holds the rule as a pure function of numbers: `Member(mass, mass_max, radius_au, radius_max_au, sudden)`, `pair_history(axis, e, a, b)` (the orbit of the past, each member's filling of its lobe at that orbit's closest approach, whether a radius reached the closest approach) and `must_be_post_common_envelope`, with `widest_forced_axis_au` giving the present axis below which a pair is forced. A pair of members that never shrank is never forced: its present filling is `classify_pair`'s. A neutron star or a black hole reports `sudden` (the adiabatic form needs a slow loss). 

**Built (GEN.1012, wiring; stacked on GEN.1009's member record).** `multistar.classify_pair` reads the rule for every pair of two single members (`_history_member` builds a `Member` from the leaf's present and largest mass and radius; a neutron star or a black hole is `sudden`): a pair that fills no lobe now, but that `must_be_post_common_envelope`, is `post-common-envelope` (it was `detached`) and carries the interacted flag; a pair already `accreting`, `contact` or `semi-detached` keeps that label. The white-dwarf-only test (`mass_loss.engulfment_radius_au(dwarf.initial)` against the axis) is gone from the classifier. The axis is not redrawn: the pair keeps the axis the general period law or its own law gave, and only its label changes, as the white-dwarf test did. A 200-seed sweep of four primary masses labelled 64 more pairs of 2,674 as post-common-envelope (white dwarf with a white dwarf, a giant or a main-sequence star), and no pair changed any other way. For planets, `mass_loss.engulfment_radius_au` is now the lookup through the three Jovian points of Mustill and Villaver 2012 (`tuning.WD_ENGULFMENT_POINTS_AU`: 2.7, 3.8 and 10 AU at 1, 2.5 and 5 solar masses, straight lines between them), replacing `WD_ENGULFMENT_BASE_AU` and `WD_ENGULFMENT_AU_PER_PROGENITOR_MSUN`; registry GM1 is that lookup and GM2 is gone. Not built: the engulfment limit from the star's largest radius times a tidal factor, and a post-common-envelope axis drawn by the general law of that case for a pair the period law placed.

### 4.4 Consistent masses when members are cut out (class B)

Two choices, and the recommended one is the first:

1. **Do not cut them out of the orbit draw.** The record of a neutron star, black hole or brown dwarf stays in the tree during the orbit draw, with its real mass and radius (a neutron star's radius is tiny; its mass is what counts). The orbits and stability margins are then right for the full tree. The system build may still not make the member a `Star`; it builds the placed object and the others it can, and the member's orbit and mass remain in the stored tree and in the gravity (the rows are not dropped). This matches Boss's scatter rule that compact companions are represented because they impact gravity (the 2026-10-15 rule, multi-star-systems.md).
2. Keep the cut and re-run the pair screens above the cut. Works, but needs the screens run twice and leaves the cut-out member's gravity out of the tree.

**As built (GEN.1013, P7).** Option 1. `multistar.without_members` is removed; `stored_trees.drawn` and `plan` keep every member, and a member of `UNBUILT_KINDS` (brown dwarf, neutron star, black hole, stripped helium) is a `DarkMember`: `StarSystem.dark_members` holds its node, slot, kind, mass and radius, so every node keeps the mass of everything below it and the orbit draw's clearance is the tree the system is built in. `_add_companions` builds only the lit members. The pair is the node of two built stars (`StoredPlan.pair_node_no`); a tree with no such node (a star beside a sub-pair of a star and a brown dwarf, say) takes the two heaviest built stars and the lowest node above both, the pair classes holding the orbit of that node, turned half way round when the heavier hangs under its slot 1. A tree with fewer than two built stars is a single star with its dark members about it (`star_slots` of one place, `binary_type` `None`). The planet windows (`planet_zones`) and the star offsets count a dark member's mass; it gives no light (a brown dwarf's and a stripped helium star's own luminosity are not counted yet) and has no planets. Storage: table `dark_members` (migration 0194, schema v194), read and written by `db/orbit_trees.py`; `physics/clearance.stored_rows` reads it as `stars` rows with `dark` true, so the verifier screens the whole tree.

This needs one check in the code (whether a system build can carry a non-`Star` member in its node tree) before the lane commits. It does not block the other parts.

### 4.5 One rule at every depth (class D)

Replace the two rules (Roche for two singles, extent for anything else) with one: every node has a **keep-out sphere** about its barycentre.

- A single member's keep-out radius is `max(radius, Roche-lobe-based radius toward its partner)`; a subtree's is its outer orbit's apocentre plus the keep-out of its larger side (`a (1 + e)` plus the member's), which is what `built_tree_clears` already sums for nodes with children [S].
- A parent node's pericentre must exceed `PERICENTRE_EXTENT_FACTOR` times the sum of its two children's keep-outs. For two singles this reduces to the Roche rule.
- That one function is the floor for the truncated draw, the check in the verifier, and the base of the planet ledger. It can run over a tree of any size, so 12 members or 20 work the same way.

### 4.6 Planets, moons and belts: the clearance ledger (class E)

Build a **ledger** per system, from the records and the finished nodes: for every node and every single star, the intervals of semi-major axis (around that star or around that node's barycentre) where a body can orbit stably, and the intervals that are not allowed.

- **Around a single star in a binary (S-type):** inside the Holman and Wiegert (1999) critical semi-major axis `a_c` of that star toward its companion: `a_c / a_b = 0.464 - 0.380 mu - 0.631 e + 0.586 mu e + 0.150 e^2 - 0.198 mu e^2` [P]. I checked the repo's S-type and P-type (`1.60 + 5.10 e - 2.22 e^2 + 4.12 mu - 4.27 e mu - 5.09 mu^2 + 4.61 e^2 mu^2`) coefficients against the paper's Equations 1 and 3 and they match exactly, with the paper's stated errors of 0.006 to 0.36 on them [P]. The fit was made for **coplanar test particles, 0.1 <= mu <= 0.9 and e up to 0.7 to 0.8**, integrated for 1e4 binary periods [P]; so a companion lighter than 0.1 of the pair's mass (a 0.09 solar-mass star beside a solar-mass one has mu of 0.08) and a pair more eccentric than 0.8 are outside the fit, and the ledger must mark those intervals as extrapolated and keep a margin [I]. In a hierarchy, the star's limit is the minimum over its companions at every level above it, using the companion subtree's mass as the perturber.
- **Around a close pair (P-type):** outside the P-type critical axis (`1.60 + 5.10 e - 2.22 e^2 + 4.12 mu - 4.27 e mu - 5.09 mu^2 + 4.61 e^2 mu^2` times `a_b`) [R] and inside the S-type limit of that pair as a body toward the next companion up. That makes the allowed zone around a subtree the interval `[P-critical of its inner orbit, S-critical of its outer orbit]`, and the rule is the same recursion at each level.
- **Inside the star:** outside the star's keep-out sphere, plus the tidal and engulfment room for a star that will be a giant (the history rule: a planet inside the giant tip radius did not survive; one outside it has expanded by `expansion_factor`). The repo's white-dwarf engulfment radius (1.5 AU for a one solar-mass progenitor plus 0.8 AU per further solar mass, `tuning.WD_ENGULFMENT_*`) is **smaller than Mustill and Villaver (2012)**, who find the outermost engulfed planet at about 2.7 AU for a one solar-mass star, clearing Jovian planets out to about 10 AU and terrestrials out to 15 AU or more for a five solar-mass star, against 4.7 AU in the repo [P, S]. The constants are the lower side of the published range and need a registry row and a second look (section 7).
- **Planets against each other:** the mutual-Hill separation the code uses (`MUTUAL_HILL_RADII_SEPARATION`, Chambers, Wetherill and Boss 1996 [R]), inside the allowed interval.
- **Moons:** inside a fraction of the planet's Hill radius, the prograde limit `0.4895 R_H (1 - 1.0305 e_planet - 0.2738 e_moon)` and the retrograde `0.9309 R_H (1 - 1.0764 e_planet - 0.9812 e_moon)` (Domingos, Winter and Yokoyama 2006, read in the paper [P]) and outside the Roche limit of the planet. The repo already holds the 0.4895 (`tuning.MOON_PROGRADE_STABLE_HILL_FRACTION`, used in `orbit_elements.py` and `hosting_gate.py`) without the eccentricity terms and without the retrograde value, and `moon_orbits.HILL_ESCAPE_FRACTION` is 0.4 [S]. The ledger should use the eccentricity-dependent forms.
- **Belts:** an annulus that must lie inside one allowed interval and between two planets' Hill zones.

The planet placer then samples from the allowed intervals (truncated sampling again: log-uniform or the existing law restricted to the union of intervals, inverse over the union). A body whose law has nothing left in its zone is not placed, and that is a counted event, "no room". The existing floors (`_orbit_floor_au`, `min_orbit_distance`) become reads from the ledger, and the per-star `a_crit` code is removed with no shim.

The ledger is derived data: it can be rebuilt from the stored records and nodes, so the verifier can recompute it independently.

### 4.7 The verifier (the proof)

`physics/clearance.py` (name for the lane to settle): a pure function `audit(system_rows) -> list of violations`, reading only what a finished system stores. Checks, for every pair of bodies:

1. Stars: closest approach versus summed radii and versus the Roche lobe fill (cap `MAX_FILL_FACTOR`) except for pairs the record marks as interacting.
2. Subtrees: pericentre of each node versus the keep-out spheres below it.
3. Stability: the Mardling and Aarseth and Eggleton and Kiseleva criteria on each triple of the hierarchy (pairwise recursion for deeper trees), with the stored masses.
4. History: section 4.3 on each pair.
5. Planets, moons, belts: inside allowed intervals; no two within the mutual-Hill separation; no body within a star's keep-out or the Roche limit of a host.
6. No two bodies share a position at any epoch the test samples (a coarse time-sweep of the stored orbits).

It must **not** import the placement code; it re-derives each limit from the registry constants and the formulas, so that a bug in the generator's version of a limit cannot hide in the proof. The test generates N systems over seeds and forced star counts 1 to 12 and requires zero violations; every violation is a failing seed in the report.

## 5. Hierarchical stability for more than three bodies

What the repo already does [S, P]: `physics/stability.py` offers Mardling and Aarseth 2001 (MA01), Eggleton and Kiseleva 1995 (EK95) and Equation 4 of Vynatheya, Hamers and Mardling 2022 (arXiv 2207.03151), and `auto` uses Equation 4 inside its fitted range. I compared the code with the paper: EK95's `Y_crit` and Equation 4 (with the effective inner eccentricity of Equation 3) are the paper's forms term for term [P].

What that leaves open for deeper trees: the code applies the triple test to each node and its parent, a "nested triples" approach. Vynatheya, Mardling and Hamers 2023 (arXiv 2301.09930, read in full for the numbers below) integrated 5e5 quadruples of each type with a direct N-body code and found that nested triples are not enough [P]:

- For 2+2 quadruples the nested-triple classifier is right 88% of the time; for 3+1 only 66%, which they say is only slightly better than a random classifier for the false-stable cases.
- MA01 and the nested-triple approach **overestimate the number of stable systems**: the precision for "stable" is 0.85 (2+2) and 0.59 (3+1). In a 3+1 system the intermediate orbit's eccentricity is excited by the outermost star and that destabilises the inner triple, which a triple test cannot see.
- Their sets are coplanar, and they give simple multi-layer-perceptron models for 2+2 and 3+1 quadruples (94% and 93% accurate) as Python scripts on GitHub. I did not download them (a machine-learning dependency and a licence check would be needed first).

So the "keep the nested triple test with a margin" recommendation of the first draft is **not enough for 3+1 trees**, and deep trees (five or more stars) have no published criterion at all [I]. Revised recommendation: (a) keep the nested-triple test as the cheap first screen with a per-level safety margin registered as a stand-in; (b) treat a 3+1 node as a different case from 2+2: the intermediate orbit's pericentre must clear the inner triple's required ratio with the outer star's perturbation added (a stand-in factor until the 3+1 model is read), and the verifier flags every 3+1 node as "lower confidence"; (c) a research part reads the paper's model file and decides whether to port it as a table or a small formula (no runtime dependency on a machine-learning library). The verifier should use a different criterion from the generator where one exists (EK95 in the verifier, Equation 4 in the generator), so that a systematic miscoding shows as disagreement.

## 6. Time: what "never overlap" has to mean

The snapshot is the present epoch. A tree must also not overlap during the time it is valid: a pair whose pericentre is just outside the keep-out now but whose Kozai-Lidov cycle drives the eccentricity up is not safe. The Kozai screen (`kozai_verdict`, `kozai_max_eccentricity`) already enters the draw; the verifier uses the maximum eccentricity of the cycle for the pericentre test (`e_max` from the Kozai formula), not the present value, for any node flagged `kl_active`. Beyond that, long-term secular and tidal evolution is out of scope and belongs to the orbital update (version 8.3).

## 7. Stand-in values this design adds or moves (GEN.725)

Registered with a tier and a tweak level so research can replace them (Boss's relaxed source rule, 8:16 PM EDT 10-29):

| Value | Now | Proposed tier | Tweak |
|---|---|---|---|
| `MAX_FILL_FACTOR` 1.2 | design default | Marginal educated guess from physics | L1 |
| `WD_ENGULFMENT_POINTS_AU` (2.7, 3.8, 10 AU at 1, 2.5, 5 solar masses; GEN.1012, replaces the 1.5 + 0.8 lookup) | lookup through the paper's Jovian points [P] | Confident Unverified until re-fitted to their Figure 7 | L1 |
| Holman and Wiegert applied for mu < 0.1 or e > 0.8 | out of the fit range [P] | educated guess from physics (extrapolated) | L1 |
| 3+1 quadruple stability margin (new) | nested triples give 66% accuracy, over-predict stable [P] | educated guess from physics | L3 until the 2023 model is read |
| `EXTENT_AXIS_FACTOR` 2.5, `PERICENTRE_EXTENT_FACTOR` 2.0 | design default | Marginal educated guess from physics | L3 |
| Per-level stability margin (new) | none; EK95 is reliable to about 20% and its authors add 10 to 20% [P] | Marginal educated guess from physics, 10 to 20% | L1 |
| Engulfment tidal factor (replaces `WD_ENGULFMENT_*` as a lookup) | lookup | Confident Unverified, once checked against Mustill and Villaver 2012 | L1 |
| Moon Hill fractions: prograde 0.4895 in code without eccentricity terms; retrograde 0.9309 not in code | read in Domingos et al. 2006 [P] | **Verified** (T1); the code lacks the eccentricity terms and the retrograde limit | L0 |
| `MUTUAL_HILL_RADII_SEPARATION` 10 | above the Delta of about 8 to 9 where five-planet systems last 1e9 to 1e10 orbits [P, Obertas 2017] | Verified as consistent with the literature (T1 for the threshold) | L0 |
| S/P-type coefficients (Holman and Wiegert 1999) | in code, equal to the paper's Equations 1 and 3 [P] | **Verified** (T1) inside mu 0.1 to 0.9, e up to 0.8 | L0 |
| Eggleton Roche lobe, EK95, Equation 4 of Vynatheya 2022 | in code; EK95 and Equation 4 equal the 2022 paper's forms [P]; the Roche formula is the standard Eggleton form but his 1983 paper was not opened | Verified for Roche (Eggleton 1983 read), MA01 (read) and Equation 4; EK95 form Verified as printed in the 2022 paper, its 20% accuracy read in the original | L0 |
| Raghavan 2010 period law mean 5.03, sigma 2.28 (log10 days) | in code; equal to the paper [P] | Verified | L0 |

None of the papers cited here is on Boss's list of papers he cannot reach, but Holman and Wiegert, Mardling and Aarseth, Eggleton and Kiseleva, Raghavan, Moe and Di Stefano, Tokovinin, Domingos et al., Chambers et al., Veras et al., Mustill and Villaver, Villaver and Livio, Vynatheya and Hamers were all **recalled, not opened**. A research lane with access should confirm the formulas and numbers before the registry rows are marked Verified; where a number cannot be confirmed, use the second-hand rule.

## 8. TEST.476 (the redraw-rate test) and how it relates

Boss asked for a test that measures how often items are redrawn to make them fit. Here is how it fits this design.

**Counters.** A counter beside every retry loop, through the one random wrapper (`draw.py`), so counting does not change the stream. The retry loops found in the stellar part: pair period tries (`PAIR_RETRIES` 60), tree restarts (`TREE_RESTARTS`), drop-lightest, `lighten` mass cuts in `_star_leaf`, the GEN.1004 whole-set redraw (`_companions_fit`), the Kozai and gravitational-wave rejects, the hydrogen-limit redraws of the roster. In the planet part: mutual-Hill pushes, floor pushes, belt redraws. Each loop gets a named kind.

**Metrics.** For each kind: events per generated system, by star count (1, 2, 3, 4, 5+) and by placed-star class, over a seeded sample that includes forced triples and quadruples (`STARS=3`, `STARS=4`) so the rare cases appear. Report the rate and its count with a binomial interval; the census output of TEST.430 carries the table.

**Two kinds of gate.**
- **Ceilings** for loops that legitimately remain (Kozai reject, shape fail, no-room for a planet): pinned to a registered value with a tier; fails when a rate doubles (as the TODO text says).
- **Invariants** for the classes this design closes: the overlap redraw count is **exactly zero**, not "small". After the build, the only legal redraw reason for stars is "no allowed range for this shape" (4.2 step 4).

**Order.** TEST.476 goes first, before any output change, on today's code: it gives the baseline (how often each loop fires today, including the GEN.1004 redraw), which is the number the design's Done clause ("falls to a small stated rate") is measured against. It is lock-free (a counter, no output change).

**Built (TEST.476).** `draw.redrawn(kind)` and `draw.counting()` in `planetgen/util/draw.py`, 35 kinds in `planetgen/generation/redraw_kinds.py`, the harness `tests/audit/redraw_rates.py` and the tests in `tests/test_redraw_rates.py`; the report is `redraw-rate-latest.md` beside the TEST.430 reports. The baseline before this design's build is `tests/data/redraw_rate_expectations.json` (8,700 seeded systems measured on main after merge trains 1 to 5: 600 ordinary, 300 one-star, 600 two-star, 1,500 each of three to six stars, 600 each with a habitable world and belt requested and with habitable worlds barred), registered as stand-ins of tier T9 and tweak level L3. Ceilings are twice the baseline plus one event of the reference's size, held against the run's mean less two standard errors of the per-system counts (the events of a system come in clusters). The overlap invariant reads five known seeded systems, one or more of each overlap kind, because a random sample of a dozen seldom meets an overlap; it is a strict expected failure and GEN.1015 removes the mark, the sensitivity test beside it, and the three overlap kinds' call sites. Baseline overlap rates per generated system (3, 4, 5, 6 stars): `pair.overlap` 0.115, 0.057, 0.101, 0.838; `tree.restart_overlap` 0.011, 0.069, 0.240, 0.451; `companions.overlap` 0.045, 0.008, 0.021, 0.041.

## 9. Optional audit with an N-body code

Not part of the build. A script can integrate a sample of generated hierarchies for a few thousand outer periods with an open N-body code (REBOUND, IAS15 / WHFast; Rein and Liu 2012, Rein and Spiegel 2015 [R]) and report the fraction that disrupt. A new dependency and a long run; useful once, as evidence for the stability margin of section 5, and kept out of the repo's runtime and the test gate.

## 10. Build plan in small parts (for TODO)

Each part is one lane-sized item. "Lock" means it changes generated output (rides a train, epoch, golden); "free" means it does not.

| Part | What | Lock | Depends on |
|---|---|---|---|
| P0 | Give the tree the primary's [Fe/H] (`system._drawn_plan`: pass `fe_h` to `multistar.Placed`); also pass the same age and mass the build uses, remove the `lighten` and `draw_age_yr` differences or record them | lock | none |
| P1 | TEST.476 counters and baseline report (existing) | free | none |
| P2 | The verifier `audit()` and its seed-sweep test, run against today's code; known failures listed as expected-fail with their seeds (the GEN.1004 cases should already pass) | free | none |
| P3 | `MemberRecord` and `Star.from_record`; tree and build read the record; schema adds `radius_max_au`, `mass_max_solar`; the radius-equality test replaces `MASS_TOLERANCE` as the guard | lock | none |
| P4 | One keep-out rule at every depth (4.5) replacing the Roche-or-extent split | lock | P3 |
| P5 | Truncated sampling of period and eccentricity (4.2), deleting the 60-try loop and the restart loop; shape step reports "no allowed range" | lock | P3, P4 |
| P6 | History rule (4.3), replaces the white-dwarf-only engulfment test in `classify_pair` | lock | P3 |
| P7 | Compact and brown-dwarf members kept in the tree for the orbit draw (4.4), pending the code check | lock | P3 |
| P8 | The clearance ledger and planet, moon, belt sampling from it (4.6), replacing `_orbit_floor_au`, the per-star `a_crit` and `min_orbit_distance` reads | lock | P3, P5 |
| P9 | Delete `built_tree_clears`, `_companions_fit`'s clearance branch and `_built_stars_clear` (GEN.1004's fix), no shim; TEST.476 invariant set to zero | lock | P5 (and P8 for the planet part) |
| P10 | Stability margin per level; 2+2 quadruple fit if the verifier shows it is needed | lock | P2, P5 |

P1 and P2 can start now on any free lane and give the evidence base; P3 is the pivot. P4 to P7 are small and can run in parallel after P3. Parts P3 to P9 each need an epoch, golden refresh and census re-tag on a merge train (Boss's TEST GATE rule: each merge tests only what changed; `scripts/test.py changed --resume`), and each should carry `--affects` for the star and planet suites. Testing mode applies: new epochs and migrations are fine, no backward-compat shim.

## 11. Open points, answered

1. **Does `evolve_star` use its `rng` in a way that makes two calls with the same inputs differ?** Mostly no. With [Fe/H] given, an ordinary star's state is a function of mass, age and [Fe/H] (it reads a MIST track); the `rng` is used to draw the birth [Fe/H] when none is given (`birth_fe_h`), and inside the horizontal-branch/giant state and the Wolf-Rayet state (not examined further) [S]. The mismatch is in the **inputs**: the tree's [Fe/H] is 0 where the build's is the primary's [C] (section Summary); a smaller part is a radius scatter of up to about 4% between `_star_leaf` and the built `Star` [C]. The `lighten` mass cut (20% steps when a star has collapsed) and the `draw_age_yr` fallback also change the tree's mass or age from the build's; I did not measure those two.
2. **Can a system build carry a non-`Star` member (neutron star, black hole, brown dwarf) in its node tree?** Still not answered; `stored_trees` documents that a tree with such a member builds the placed object alone and marks `complete=False`, and a drawn tree has them cut out first [S]. The build lane must check what the node tree and `member_stars` accept before choosing 4.4 option 1 or 2.
3. **Does the moon code use the Domingos fractions?** Partly: it has the prograde 0.4895 without the eccentricity terms, and no retrograde fraction (section 4.6) [S].
4. **Planet ledger scope.** Unchanged: derived, not stored.
5. **New:** whether the 1-line [Fe/H] fix should go in at once as its own small item. It changes generated output (any system whose primary is not at [Fe/H] 0 and has a companion) so it rides a train with an epoch, but it removes about 60 to 70% of the redraws and is the cheapest win. It is the new first output-changing part, P0, in section 10.

## 12. Sources and what was checked

Checked means read in this session. arXiv text was fetched and read for the first group; the second group could not be opened (no arXiv copy or paywalled) and stays below Verified.

| Source | Used for | Status |
|---|---|---|
| Holman and Wiegert 1999, AJ 117, 621 (astro-ph/9809315) | S-type and P-type critical axes; fit range | **Read**; repo coefficients equal Equations 1 and 3; range mu 0.1 to 0.9, e to 0.8, coplanar, 1e4 periods |
| Vynatheya, Hamers and Mardling 2022 (arXiv 2207.03151) | Triple criteria: EK95, MA01, Equation 4 | **Read**; repo's EK95 and Equation 4 equal the paper's |
| Vynatheya, Mardling and Hamers 2023 (arXiv 2301.09930) | Quadruple stability: nested triples overestimate stable systems, 88% (2+2) and 66% (3+1) | **Read**. This replaces my earlier wrong citation "Hamers 2021; Vynatheya and Hamers 2022" |
| Raghavan et al. 2010, ApJS 190, 1 (arXiv 1007.0414) | Period law log P mean 5.03, sigma 2.28 | **Read** (the values) |
| Moe and Di Stefano 2017, ApJS 230, 15 (arXiv 1606.05347) | Massive-star period law peaks near log P = 3.5 | **Read** (abstract statement) |
| Tokovinin 2018, ApJS 235, 6 (arXiv 1712.04750) "The Updated Multiple Star Catalog" | Observed hierarchies | Exists, title and authors confirmed; contents not used for numbers |
| Veras, Wyatt, Mustill, Bonsor and Eldridge 2011, MNRAS 417, 2104 (publisher PDF uploaded by Boss, and arXiv 1107.1239) | Adiabatic mass-loss orbit growth, where it fails | **Read** (sections 1 and 2.4): in the adiabatic regime `a = a0 (1 - alpha t / mu0)^-1`, so a planet doubles its separation if its 1 solar-mass parent loses half its mass; holds when the mass-loss time-scale is long against the orbital period. Equal to the repo's `mass_loss.expansion_factor` for a test mass. **Verified** |
| Mustill and Villaver 2012, ApJ 761, 121 (arXiv 1210.0328) | Engulfment limits after the giant phase | **Read** in part: about 2.7 AU for a 1 solar-mass star; Jovian planets cleared to about 10 AU and terrestrial to 15 AU or more for 5 solar masses; repo's constants are smaller |
| Villaver and Livio 2009, ApJ 705, L81 (arXiv 0910.2396) | Orbital evolution of giant planets around giants | Exists, confirmed by title and authors; not used for numbers |
| Domingos, Winter and Yokoyama 2006, MNRAS 373, 1227 (publisher PDF uploaded by Boss) | Moon Hill-fraction fits | **Read**: prograde `a_E = 0.4895 (1 - 1.0305 e_P - 0.2738 e_sat)`, retrograde `a_E = 0.9309 (1 - 1.0764 e_P - 0.9812 e_sat)`, in units of the Hill radius, from the restricted elliptic three-body problem; planet orbits tested at 0.1 AU with e_P 0 to 0.9. Equal to the values quoted in 0912.3484. **Verified** |
| Obertas, Van Laerhoven and Tamayo 2017, Icarus (arXiv 1703.08426; uploaded by Boss in place of Chambers 1996) | Planet spacing in mutual Hill radii | **Read** (introduction, section 2, 5.2): the time to a first close encounter grows as `log(t_c/t_0) = b*Delta + c` for about 3.5 < Delta < 8 (Chambers et al. 1996, confirmed by later studies); five equal planets become long-lived beyond Delta of about 8 to 9 (stability to 1e9 and 1e10 orbits for over 5,000 systems with Delta > 8), and a two-planet Hill-stable limit exists (their Equation 5). The repo's `MUTUAL_HILL_RADII_SEPARATION` of 10 is above that threshold. Chambers 1996 itself is covered second-hand through this paper |
| Rein and Spiegel 2015 (arXiv 1409.4779) | IAS15 integrator | Exists, confirmed; optional audit only |
| Rein and Liu 2012 (REBOUND) | N-body package | Not located by search; the package is real, citation not checked |
| Mardling and Aarseth 2001, MNRAS 321, 398 (ADS scan) | Triple criterion, Equation 90 | **Read** (second pass): `R_p/a_in = C [(1+q_out)(1+e_out)/(1-e_out)^(1/2)]^(2/5)`, C = 2.8 empirical; "holds for q_out <= 5" (beyond it exchange, not escape, is the failure); inclined orbits can be up to about 30% closer. Repo's `mardling_aarseth_pericentre_ratio` equals it, including the (1 - 0.3 i/pi) factor. The code does not enforce the q_out <= 5 limit; a lane should check where `ma01` is used outside it |
| Eggleton 1983, ApJ 268, 368 (ADS scan, 2 pages) | Roche-lobe radius | **Read**: Equation 2 `r_L = 0.49 q^(2/3) / (0.6 q^(2/3) + ln(1 + q^(1/3)))`, better than 1% for all q. Repo's `roche_lobe_fraction` equals it. **Verified** |
| Eggleton and Kiseleva 1995, ApJ 455, 640 (ADS scan, image only) | Triple criterion `Y_crit` | **Read the conclusions page**: the criterion is "reliable to about 20%", weakest at extreme mass ratio (star plus planet), their own code adds 10% to 20% for caution. The formula's page was not transcribed; the repo's form equals the one printed in the 2022 paper [P]. The 20% accuracy is the source for the per-level stability margin of section 5 |
| Chambers, Wetherill and Boss 1996 (paywalled); Jeans 1924; Hadjidemetriou 1963; Tokovinin 2008 | Mutual-Hill spacing, mass-loss theory | **Not opened**; Chambers is covered through Obertas 2017, Jeans and Hadjidemetriou through Veras 2011; Tokovinin 2008 is not used for a number |
| In this repository | Code and docs named in the text | **Read** |

## 13. What to change in TODO (for the TODO thread)

- **GEN.1006:** keep the title; replace the "Do" text with this note's parts P0 and P3 to P10; Owner: Star 2 (or first free lane); Prerequisites: P1 and P2 (new items) for the baseline and the verifier, TEST.476 first; Related: GEN.1004, GEN.725, TEST.162, TEST.430, TEST.476. The Done clause stands, with "falls to a small stated rate" made "overlap redraw count is zero and every remaining redraw kind has a registered ceiling".
- **New items (suggested split, TODO to number):** P2 the verifier and seed-sweep test (lock-free, Bug 2 or any free lane); P3 `MemberRecord`; P4 keep-out rule; P5 truncated sampling; P6 history rule; P7 compact members in tree; P8 clearance ledger for planets, moons and belts; P9 removal of the GEN.1004 fix; P10 stability margin and quadruple fit.
- **TEST.476:** add the two gate kinds (ceilings and the zero invariant) and the baseline-first order.
- **GEN.725:** the stand-in rows of section 7, including the engulfment constants that sit below Mustill and Villaver 2012 and the 3+1 quadruple margin.
- **Research part (new, free):** read the 2023 quadruple models, decide how to port a 3+1 criterion without a machine-learning dependency; verify Domingos 2006, Eggleton 1983, MA01, EK95 and Chambers 1996 from the papers (a lane with access).
- **GEN.1004:** mark "replaced by GEN.1006 P9" when P9 is scheduled.

## 14. How this was checked (scripts, in the scratchpad of the research thread)

- Rates: a script forced `STARS=n` for seeds 1000 to 3499 with `PLANETS=False`, wrapped `multistar.built_tree_clears` to count calls and failures, and counted systems drawn more than once. Run for n = 3, 4, 5 and 6, once on current main and once with `fe_h` passed to `multistar.Placed` in `system._drawn_plan` (scratch edit, reverted; not committed). Environment: datasets MIST isochrones, white-dwarf models, Baraffe, spectra and Mamajek fetched with `python -m planetgen.cli.datasets fetch`, and the NLTK word list.
- Radius comparison: 400 forced binaries; for the second star the tree's `_star_leaf` radius (same birth mass, the primary's age and [Fe/H]) against the built `Star` radius; no star differed by more than 5%, the largest ratios 0.96 to 1.045.
- Papers: arXiv PDFs converted with pdftotext and searched for the formulas and quoted numbers.
- Not done: the quadruple of the thread's seed was not re-found by seed; the 0.033 AU giant shows in my failures (n = 4 and n = 6 lists). `lighten` and `draw_age_yr` effects not measured. No part of this was run on the test database.
