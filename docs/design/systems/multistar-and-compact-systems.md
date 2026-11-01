# Multi-star systems and compact-object primaries

How a system of up to twelve stars (the group cap `MAX_STARS`, GEN.267; Boss's 2026-10-09, 2:43 AM EDT cap of seven was superseded by the cap of 12, Boss, 2026-10-15, 6:08 PM EDT) is drawn, stored, named, shown on the maps and
updated, how often each size of system occurs, which hierarchies are dynamically
stable, and what a white dwarf, neutron star or black hole primary changes for
the bodies around it. It is the design note GEN.128 asks for and the input to
GEN.129 and GEN.130. How many neutron stars and black holes exist and where
they sit is not repeated here: those densities and the scatter are in
[interstellar-object-rates.md](interstellar-object-rates.md) (GEN.1, GEN.100),
the placement rules in [anomalies.md](anomalies.md) ("Placement by population"),
and the galaxy's stellar density in
[galaxy-disk-density.md](../galaxy/galaxy-disk-density.md). This note covers only the
company a compact object keeps and the planets that survive it.

Informs: GEN.128, GEN.129, GEN.130, GEN.62, GEN.71, GEN.84, GEN.85, GEN.86, GEN.87, GEN.88, GEN.89, GEN.100, GEN.103, GEN.109

Status: research, 2026-10-09; decisions marked "Boss" are his, everything else is a recommendation

Evidence tags: [S] seen in a search result (URL in Sources), [C] computed or
simulated in the research sandbox, [R] recalled and unconfirmed. The research
environment could only read search-result text, not the papers themselves; every
[R] is on the Evidence notes list. The [C] values in the probability table, the
period shares, the compact-object sizes and the Holman and Wiegert comparison
were rerun while writing this document.

## Decisions already taken

- Scope (Boss, GitHub issues #777 and #778, 2026-10-09 06:43Z): "scientifically
  accurate star systems with up to 7 stars, this is going to be complex but
  that is the highest number of stars we've seen in orbit around each other",
  and "exotic star systems that have black holes, neutron stars, or similar as
  the central star for systems, binary systems."
- Storage (Boss, GEN.128): a tree of pairs, each pair's orbit around its
  barycentre.
- Wide-pair planet names (Boss, 2026-10-07 17:11Z, GEN.71): a wide binary's
  planets are "<word 1> I", "<word 1> II" and the second star's "<word 2> I",
  never "A I" ([object-ids.md](../architecture/object-ids.md), "The naming key"). GEN.62 (done,
  PR #393) names the stars: `Voranthis A` and `Voranthis B` for a close pair,
  `Voranthis Kelmoor` and `Voranthis Pikkita` for a wide one.
- Influence rule and time step (Boss, 2026-10-09 and 2026-10-07 17:11Z): an
  object's influence radius is the Hill sphere of the largest nearby object, and
  each run advances by real elapsed time
  ([orbital-updates.md](../orbits/orbital-updates.md) sections 1 and 10.4). Section 2
  applies them to hierarchies.
- GEN.100 (done) scatters neutron stars and black holes galaxy-wide; GEN.130
  builds on that scatter instead of drawing its own.

## Summary

- **GEN.129: go, as a staged build.** A system is a tree of pairs (`orbit_nodes`,
  N-1 rows for N stars); one recursive generator makes every N from 3 to 7. The
  cost is plumbing: 34 source files and 42 test files mention the binary fields
  (`binary_*`, `is_binary`, `secondary_star`, `wide_binary`, by grep).
- **Frequency (3.3):** per system, single 69.1%, binary 25.6%, triple 4.0%,
  quadruple 1.07%, quintuple 0.24%, sextuple 0.029%, septuple 0.0030% (about 1 in
  33,000), mean 1.38 stars. Seven is the observational ceiling, with no seven-star
  system certain.
- **Stability needs a Kozai-Lidov screen (5.2).** Without it 16.7% of tested
  levels failed within 10^4 inner orbits; with it 598 of 600 held, and 30 of 30
  whole systems held [C].
- **The binary generator missed 43% of real binaries (4):** nothing was drawn
  between 0.26 and 50 AU. Built in GEN.448 (4).
- **GEN.130: go in slices, no-go on four items (8).**

## 1. What the repo does today (checked 2026-10-09)

- **Close (P-type) pair.** `BinaryStarProxy` (`generation/binary.py`) merges two
  stars into one effective star: separation `uniform(0.05, 0.25)` AU plus both
  radii, eccentricity 0, circular position, random orientation, mass ratio
  from the law of [R3108] over 0.1 to 1 (`generation/mass_ratio.py`, GEN.951, with `tuning.BINARY_MASS_RATIO_RANGE` the ends), same age as the primary
  (GEN.53). The orbit lives in `binary_*` columns on `star_systems`; the proxy has
  no `stars` row ([database-schema.md](../../reference/database-schema.md), "Binary
  configuration").
- **Wide (S-type) pair.** `WideBinaryPair` (`generation/wide_binary.py`):
  separation from Moe and Di Stefano's companions per decade of period, between
  0.25 AU and the place's survival limit (GEN.448; it was log-uniform from 50 to 10,000
  AU), thermal eccentricity capped at 0.8, live position circular ("a deliberate,
  documented simplification"). Each star keeps its planets and its Holman and
  Wiegert `wide_binary_a_crit_km`.
- **Choice.** `StarSystem._companion` rolls
  `BINARY_SYSTEM_PROBABILITY_BY_SPECTRAL_CLASS` (O 0.90, B 0.65, A 0.55, F 0.47,
  G 0.44, K 0.40, M 0.26) for "has a companion", then shares it between the close
  range (0.05 to 0.25 AU, a P-type pair), the inner range (0.25 to 10 AU) and the outer range (10 AU and out)
  by the companions per decade of period of Moe and Di Stefano's Table 13
  (`generation/companion_classes.py`, `tuning.COMPANIONS_PER_DECADE_BY_PRIMARY_MASS`), for the
  primary's mass. The close and inner ranges are the ones Moe, Kratter and Badenes measured (periods under 10^4 days), so
  since GEN.103 the odds of the two together are multiplied by `F([Fe/H]) / F(0)` of the primary's own [Fe/H]
  (every mass); the outer range is not, and a wide pair's axis is drawn up to the place's survival limit
  ([placement-by-population.md](../generation/placement-by-population.md) 3.1 and 3.3).
- **Holman and Wiegert** is implemented in `physics/orbits.py`
  (`holman_wiegert_critical_semimajor_axis`, `holman_wiegert_circumbinary_a_crit_au`);
  `StarSystem._orbit_floor_au` applies the P-type limit.
- **Stored.** `stars.role` is CHECK IN ('primary','secondary','single');
  `star_systems.binary_configuration` is 'close' or 'wide'; `planets.star_id` is
  the host for S-type planets and NULL for a close pair's planets.
- **Maps and update.** `web/maps/systemscene.py` gives each body an `around` of
  `barycenter` or `star:<id>`; `physics/body_positions.positions_at` balances a
  close pair on the barycentre by `secondary_mass_fraction`;
  `db/store.py advance_orbital_phases` advances phases with one set-based UPDATE
  per table ([orbital-updates.md](../orbits/orbital-updates.md), "Positions at any time").
- **Compact objects.** `generation/phenomena/compact_remnant.py` has `BlackHole`
  and `NeutronStar` as `Star` subclasses. They anchor a `StarSystem` only through
  the `compact_remnant` argument (`phenomenon --anchor-system`), and that path
  skips binary generation (`compact_remnant is None and
  self._companion()` in `StarSystem.__init__`). **The module
  docstring cites `docs/design/exotic-phenomena.md`, which does not exist**; a
  repo-wide grep finds the path only there (the words "exotic-phenomena" elsewhere
  are prose in `docs/reference/cli.md` and a test CLI). The design home is
  [anomalies.md](anomalies.md) plus this note; the docstring should point there.
- **Collapsed stars stay out of normal systems on purpose.**
  `physics/stellar_evolution.evolve_star` returns None for a star at or above 8
  Msun past its giant phase; `Star.generate` redraws such a primary and lightens a
  companion by 0.8 until it is alive. White dwarfs are ordinary class VII `Star`
  objects (a lighter white dwarf primary is swapped below its companion), so white
  dwarf primaries and companions occur already; neutron star and black hole ones
  cannot.

## 2. How the source documents update a multi-body system

"Orbital Update Full Algorithm.md" and "Orbital Position and Vector Update
Algorithms.md" are condensed in [orbital-updates.md](../orbits/orbital-updates.md); neither
has a notion of a bound group of stars. What they imply, and what changes:

| Source statement | Section | Status for hierarchies |
|---|---|---|
| Stage 1 integrates every entity with `parent_system_id IS NULL` (stars, remnants, rogues) on its own, under the point masses around it | Full Algorithm, Master Simulation Workflow | A star pair would be two independent roots: a partner at 0.1 AU pulls about 4e12 times harder than a star 1 pc away, and a period of days is far inside any galactic step. **A hierarchy is one root at its barycentre with the total mass**; its stars are leaves of the tree. |
| Stage 3 (`ExecuteStarSystemUpdates`): Phase A shifts every child by its star's displacement; Phase B advances each child's Kepler orbit about the star | Full Algorithm, Star System Engine | Stays, and generalizes. Phase A moves the root node by the barycentre's displacement; Phase B advances each node's relative orbit analytically (mean anomaly by `n dt mod 2 pi`, [orbital-updates.md](../orbits/orbital-updates.md) section 10.4); star offsets follow from the tree walk in 6.1. `parent_system_id` becomes "node or star". |
| `CheckPlanetaryHillCrossings`: occupants of a body's Hill sphere are captured, deflected or ejected | Full Algorithm, Encounters | Members of one hierarchy sit inside each other's Hill spheres by construction and are exempt from mutual Hill warnings and captures ([orbital-solvers-and-integrators.md](../orbits/orbital-solvers-and-integrators.md) section 8, "exclude hierarchical members"; [collisions-and-mergers.md](../orbits/collisions-and-mergers.md) section 3.3, "the same system"). Hill occupants are tracked per root. |
| Point-mass table with one row per body (id, mass, mu, position, velocity, tier) | Update Algorithms, Analytical Point-Mass Table Computation | One row per system at its barycentre with the total mass (GEN.109); tier from the total mass. Internal structure is read from `orbit_nodes`. |
| Top 10 influencers; `R_search = sqrt(G M_max / a_min)`; the `a_min` filter | Update Algorithms, Position and Vector Update; Full Algorithm | Superseded (Boss, 2026-10-09). A flyby inside the outermost orbit's Hill zone is a per-star event flagged "disrupts hierarchy", not integrated. |
| Velocity Verlet at fixed steps; 10 and 11.5 ly sector keys | both | Superseded ([orbital-updates.md](../orbits/orbital-updates.md) sections 4, 10.4, 11). Irrelevant inside a hierarchy, which is analytic. |

The stored pair columns already carry what a node needs
(`binary_mutual_orbital_period_years`, `_inclination_deg`, `_ascending_node_deg`,
`_phase_deg`, `_min_update_interval_years`, `_position_*_km`,
`binary_secondary_mass_fraction`). The 0.01 AU "companion stars" threshold of
[orbital-updates.md](../orbits/orbital-updates.md) section 3 is met by computing node
positions from the phase on read, as for planets: a pair at 0.1 AU moves 0.01 AU
in hours, so a next-due column would fire constantly.

## 3. Multiplicity statistics

### 3.1 What the literature says

| Quantity | Value | Tag and source |
|---|---|---|
| FGK dwarfs within 25 pc, single / double / triple / higher | 56% / 33% / 8% / 3% (454 stars) | [S] [R3646] |
| F and G dwarfs within 67 pc, n = 1 to 5 | 54 : 33 : 8 : 4 : 1; 2+2 quadruples 4% | [S] [R4384] |
| Bright systems (about V < 6), 1 to 7 components | 2718, 1437, 285, 86, 20, 11, 2 | [S] [R1223] |
| Updated Multiple Star Catalog | 17 systems with 6 components, 4 with 7, none of the 7s certain (65 UMa best) | [S] [R4385] |
| Companions per primary (q > 0.1, log P < 8) | 0.50 solar-type, 2.1 O-type; single most probable below about 2 Msun, triple above about 12, quadruple above about 23 | [S] [R3108] |
| M dwarfs | multiplicity about 26%, higher order about 3% of systems | [S] [R4724] summary; 26% is in `tuning.py` |
| O stars | 91 +- 3% multiple; over 70% exchange mass with a companion | [S] [R3892], SMaSH+ |
| O/B companion peak | log P (days) near 3.5 (about 10 AU) | [S] [R3108] |
| Solar-type period law | log-normal in log10 P (days), mean 5.03, sigma 2.28 | [R] [R3646]; every period share below depends on it |
| Close-binary fraction (a < 10 AU) against metallicity | about 10% at [Fe/H] +0.5, 24% at -0.5, 40% at -1, 53% at -3 | [R] [R3109] |

Theory allows several hundred components and observation confirms at most seven
[S, [R1529]]. The single fraction counts systems; by star count
about half of Sun-like stars sit in multiples (singles outnumber primaries of
multiples 1.28 to 1, arxiv.org/pdf/1502.04018 [S]). The generator draws a primary
from a per-star census and then adds companions, which slightly over-produces
low-mass stars: a refinement, not a blocker.

### 3.2 The extreme multiples

| System | Stars | Structure | Tag |
|---|---|---|---|
| Alpha Centauri | 3 | A and B (79.9 yr), Proxima far out | [S] |
| HD 98800 | 4 | 2+2, each a spectroscopic binary (262 d and 315 d) | [S] |
| Castor | 6 | three binaries, the outer one eclipsing | [S] |
| Mizar / Alcor | 6 | Mizar a quadruple, Alcor a binary; whether Alcor is bound has been debated | [S], [R] on the debate |
| TIC 168789840 | 6 | three eclipsing binaries, outer period about 2 kyr | [S] |
| Nu Scorpii | 7 (probably) | groups AB and CD 41 arcsec apart; AB a triple, CD a binary within a wider pair | [S] "most likely a septuple" |
| AR Cassiopeiae, 65 UMa | 7 (candidates) | structure not confirmed; 65 UMa has four hierarchy levels | [S] |

Boss's statement holds as an observational ceiling. A fair line for a page:
"seven components is the most ever proposed, with none beyond doubt".

### 3.3 Probability table the generator draws

Rows are the primary's spectral letter. These are design values, not fits: the
first three columns blend the sources above (M from Winters and Duchene and
Kraus, FGK from Raghavan and Tokovinin, A to O from Moe and Di Stefano's trend
toward more triples and quadruples with mass), and the N >= 4 mass is split by a
per-step ratio set so the FG row reproduces Tokovinin's 4% and 1%.

| Primary | N=1 | N=2 | N=3 | N=4 | N=5 | N=6 | N=7 | mean N |
|---|---|---|---|---|---|---|---|---|
| M | 0.73000 | 0.23500 | 0.02800 | 0.00560 | 0.00123 | 0.00015 | 0.00001 | 1.31 |
| K | 0.60000 | 0.31000 | 0.06500 | 0.02002 | 0.00440 | 0.00053 | 0.00005 | 1.52 |
| G | 0.54500 | 0.33500 | 0.08000 | 0.03202 | 0.00705 | 0.00085 | 0.00008 | 1.62 |
| F | 0.50000 | 0.34000 | 0.11000 | 0.04003 | 0.00881 | 0.00106 | 0.00011 | 1.72 |
| A | 0.40000 | 0.40000 | 0.15000 | 0.03741 | 0.01048 | 0.00189 | 0.00023 | 1.87 |
| B (2.1 to 5 Msun) | 0.25000 | 0.40000 | 0.22000 | 0.08962 | 0.03137 | 0.00784 | 0.00118 | 2.28 |
| B-mid (5 to 9 Msun) | 0.24000 | 0.36000 | 0.25143 | 0.10241 | 0.03585 | 0.00896 | 0.00135 | 2.37 |
| B-early (9 to 16 Msun) | 0.16000 | 0.32000 | 0.32686 | 0.13314 | 0.04660 | 0.01165 | 0.00175 | 2.63 |
| O | 0.10000 | 0.25000 | 0.30000 | 0.21354 | 0.09610 | 0.03363 | 0.00673 | 3.08 |
| Weighted by `SPECTRAL_PROBABILITIES_NORMAL` | 0.69070 | 0.25604 | 0.03987 | 0.01069 | 0.00238 | 0.00029 | 0.00003 | 1.38 |

The B-mid and B-early rows (GEN.449) take their first two columns from Moe and Di Stefano's Table 13 (single 0.24 and 0.16, binary 0.36 and 0.32,
triple or more 0.40 and 0.52; `tuning.MULTIPLICITY_BY_PRIMARY_MASS_BIN`) and split the rest across 3 to 7 stars in the B row's proportions. The chance of a companion
(1 - single) is 0.76 and 0.84 there against the letter's 0.65. The B row keeps its values for the late-B stars: the bin of 2 to 5 solar masses holds A and late B together at 0.59 +-0.11, and the letter's 0.65 sits inside that.
The weighted row below still reads the B row for every B star: B is 0.12% of the normal mix, so the two new rows move its figures by under 0.0005.

The N >= 4 total per class is 0.007 (M), 0.025 (K), 0.040 (G), 0.050 (F, A), 0.130
(B), 0.350 (O), split across N=4..7 by per-step ratios 0.22, 0.12, 0.10 (M to F),
0.28, 0.18, 0.12 (A), 0.35, 0.25, 0.15 (B), 0.45, 0.35, 0.20 (O). The bright-star
counts have gentler tails, but that sample favours rich A and B systems and
Tokovinin's volume-limited FG sample shows nothing beyond N=5, so the tails are
deliberately steep. For more exotic systems without new physics, scale the N >= 5
columns with a prevalence knob (`prevalence.scaled_chance`); to force N, use a
directive.

Against the code: M, K, G and F match `BINARY_SYSTEM_PROBABILITY_BY_SPECTRAL_CLASS`
within 0.03, but B (0.65) and A (0.55) were low against 0.75 and 0.60 (Moe and Di
Stefano: B about 0.75 to 0.8, A about 0.6 to 0.65). Table 13's bins settle it (GEN.449): A and late B
0.59 +-0.11 (the letters' 0.55 and 0.65 stand), mid B 0.76 and early B 0.84 (now the B-mid and B-early classes, by birth mass). The mean of 1.38 stars per
system replaces the "about 1.3" in
[interstellar-object-rates.md](interstellar-object-rates.md); those rates are per
star and `run_sector.sector_star_count` already sums `len(entry.star_system.stars)`,
so phenomenon counts scale correctly when N grows.

Mass ratios: keep `uniform(0.1, 1)` for companions of the primary (Moe and Di
Stefano find it close to flat). Outer companions (P above about 1e5 days) draw
their mass independently from the IMF, since those are near-random pairings [S].
Eccentricity: circular below about 7 days, otherwise thermal `f(e) = 2e` capped at
`1 - (P / 2 d)^(-2/3)` [R].

## 4. Binary period and mass-ratio redraw

With the log-normal law and a 1.5 Msun pair, the shares of binaries by separation
were [C]:

| Separation | Share of binaries | Generated before GEN.448 |
|---|---|---|
| below 0.26 AU | 6.6% (4.2% inside 0.05 to 0.26 AU) | `binary.py`, 0.05 to 0.25 AU |
| 0.26 to 50 AU | 43.3% | nothing |
| 50 to 10,000 AU | 43.6% | `wide_binary.py` |
| above 10,000 AU | 6.5% | nothing (unbound by passing stars) |

The 50 and 10,000 AU edges in this table are the old fixed ones: wide pairs start at 0.25 AU and the fixed 10,000 AU cap is replaced by the derived limit (GEN.448, Boss, 2026-10-16, 6:12 to 6:15 PM EDT).

The 0.26 to 50 AU range is where Holman and Wiegert planet zones of 0.1 to 15 AU
and the habitable zones interact most.

Built (GEN.448, no schema change): one distribution of companions over period for the
whole range replaces the close/wide coin flip (`WIDE_BINARY_DEFAULT_CHANCE` is gone) and
the nothing between 0.25 and 50 AU. It is Moe and Di Stefano's own (2017, Table 13: companions per decade of log P
at 1, 3, 5 and 7 for primaries of 1.0, 3.16, 6.71, 12 and 25 solar masses, straight lines between, flat past the ends)
rather than the log-normal law, because the lines give the massive stars their larger close share
(a log-normal fit left them about 2%). For a 1.55 solar mass pair around a G dwarf the shares of the companions
are 5.7% inside 0.25 AU, 24% from 0.25 to 10 AU, 17% from 10 to 50 AU, 31% from 50 to 1,000 AU and 22% beyond
1,000 AU up to the 11,400 AU the Sun's neighbourhood keeps [C]. A pair's axis is drawn from the same lines cut at 0.25 AU and the
place's limit, the part under 10 AU carrying the metallicity factor's weight, so the trend shapes the
separations as well as the odds. The limit is the place's alone (the survival limit of GEN.103 and GEN.447; a place that keeps no
pair wider than 0.25 AU loses the wide chance, and a forced companion is a close pair). The S-type stability limit of the pair
(`a_crit_au`, 5.6) is unchanged: a pair a few AU apart leaves each star only a fraction of an AU for planets and comets
(a forced comet is then an unmet requirement, like a forced habitable world), which is the physics, not a rule.
Eccentricity stays as in 3.3, and the `wide_binary` directive keeps its meaning (forced wide or close).
The absolute level inside 10 AU is about 0.6 of the pooled line's (a G dwarf 0.13 at [Fe/H] 0 against the line's 0.20 or 0.24 at -0.2),
(Table 13's 0.15 companions inside log P 3.7 are a count per primary and the generator's chance is the fraction of stars with one,
which may account for part of it; the cause is not settled [R]). TEST.196 checks the level against the published ranges.

## 5. Dynamical stability

### 5.1 Criteria

A pair node with children masses m1, m2 (m_in = m1 + m2) orbits a sibling of mass
m3 in its parent's orbit (a_out, e_out). q_out = m3 / m_in, q_in = m1 / m2, and i
is the mutual inclination of the two orbital planes (0 coplanar prograde, pi
coplanar retrograde).

- **[R2868]** [S]: `a_out (1 - e_out) / a_in > 2.8 *
  [(1 + q_out)(1 + e_out) / sqrt(1 - e_out)]^(2/5) * (1 - 0.3 i / pi)`. Valid for
  m3 / m_in below about 5.
- **[R1221]** [S, as quoted in arxiv.org/pdf/1408.5431]:
  `a_out (1 - e_out) / (a_in (1 + e_in)) > Y`, `Y = 1 + 3.7 / q^(1/3) - 2.2 /
  (1 + q^(1/3)) + 1.4 / q_in^(1/3) * (q^(1/3) - 1) / (q^(1/3) + 1)` with q = q_out.
- **Observed period ratios** [R]: P_out / P_in above about 5, mostly above 10; the
  generator gives a minimum of 4.8, a 5th percentile of 11.3 and a median of 233
  over 600 levels [C]. Mardling and Aarseth is stricter: a period ratio of 5 is an
  axis ratio of 2.9 for equal masses, below the 3.9 it asks for coplanar orbits.
- **Hill radius** `r_H = a (1 - e) (m / 3M)^(1/3)`: for an equal-mass circular pair
  0.55 a, and Holman and Wiegert's S-type limit 0.274 a is 0.50 r_H [C].
- **Several levels:** test every node against its own parent (the parent's orbit
  and the sibling's total mass), deeper subsystems as point masses. Per-level
  triples and whole-system integrations agree (5.4).
- **[R1990]** [S]: S-type `a_c / a_b = 0.464 - 0.380 mu - 0.631 e +
  0.586 mu e + 0.150 e^2 - 0.198 mu e^2` (mu = perturbing companion's mass
  fraction); P-type `a_c / a_b = 1.60 + 5.10 e - 2.22 e^2 + 4.12 mu - 4.27 e mu -
  5.09 mu^2 + 4.61 e^2 mu^2` (mu = the lighter star's fraction). Quoted accuracy:
  S-type about 4% typical and 11% worst for 0.1 <= mu <= 0.9, 0 <= e <= 0.8; P-type
  about 3% and 6%, with resonance "islands" of instability. 
- **Mardling and Aarseth, verified; updated criterion read from the paper.** The formula above was checked by Boss's engine on 2026-10-09 (18:45Z): P-type coplanar circular gives 2.80; S-type with `m1 = m3 = 1`, `e_out = 0.5` gives `R_p/a_in = 4.99` and `a_out/a_in = 9.98` (use as unit tests). [R4574] refine it (Equations 3 and 4 of the paper: an effective inner eccentricity `e~_in` and a new inclination bracket; 93% accuracy against 90%). Their equations were read from the paper Boss supplied at 18:51Z and are written out, with unit tests and the decided use (Boss, 18:56Z: Equation 4 for star-only triples inside the fitted mass-ratio range, this formula for planets), in [exotic-environments-planets-and-compact-binaries.md](exotic-environments-planets-and-compact-binaries.md) 2.3.

### 5.2 The criteria alone are not enough

Mardling-Aarseth and Eggleton-Kiseleva only guard against the outer body
disrupting the inner pair. A hierarchy that passes both can still destroy itself by
Kozai-Lidov cycles: at a mutual inclination between about 40 and 140 degrees the
inner eccentricity climbs toward 0.99 on a timescale about P_out^2 / P_in and the
pair merges. Direct integration of equal-mass triples (1, 1, 1 Msun; e_in = e_out
= 0.3; REBOUND IAS15 and WHFast agree) [C]:

| Mutual inclination | Result |
|---|---|
| coplanar | pericentre ratio 2 and 3 disrupted, 4 and above stable (the criterion asks 3.93: on the boundary) |
| 60 degrees | stable, e_in reaching 0.75 to 0.77 at ratio 4 to 10 (analytic maximum `sqrt(1 - 5/3 cos^2 i)` = 0.76) |
| 90 degrees | merged (e_in = 0.99) at every ratio from 4 to 10, stable only at 20, though the criterion asks only 3.34 |

Real systems avoid this because Kozai-Lidov plus tides makes tight binaries and
many real triples are near-coplanar at small period ratio. A generator must screen
or draw inclinations with a coplanar preference. The screen tested (quadrupole
maximum eccentricity, Kozai time against 10 Gyr, an octupole widening for eccentric
outer orbits, a pericentre floor of 4 times the children's extent):

```
kl_unsafe(a_in, e_in, m_in, parent, mut_incl, extent_sum, age_yr = 1e10):
    x = (5/3) (1 - e_in^2) cos^2(mut_incl)
    if x >= 1: return False                      # no Kozai window
    e_max = sqrt(1 - x)
    eps = (a_in / parent.a) * parent.e / (1 - parent.e^2)    # octupole strength
    if eps > 0.01 and 25 deg < min(mut, 180 deg - mut) < 90 deg:
        e_max = max(e_max, 0.98)
    t_kl = (2 / (3 pi)) (P_out^2 / P_in) (parent.mass / m3) (1 - parent.e^2)^1.5
    return t_kl < age_yr and a_in (1 - e_max) < 4 * extent_sum
```

"Extent" is the largest distance of any star of the subtree from its barycentre (a
star's radius, or `a (1 + e)` plus the children's extents). A first version that
used stellar radii where a subsystem's real extent applied let 7 of 180 levels
fail (mergers at inclination 73 to 127 degrees); using final extents and
re-checking after the draw took it to 0 of 180.

### 5.3 Generator algorithm

Pure Python, about 300 lines, no dependencies; 53 ms for N=7 (average 32
restarts), under 1 ms for N <= 4 [C].

1. Draw N from the table in 3.3 for the primary's class.
2. Draw the other N-1 masses as `m1 * uniform(0.1, 1)` (the primary is heaviest).
3. Draw the tree shape by random splitting: N=3 is 2+1; N=4 is 3+1 (60%) or 2+2
   (40%); N >= 5 splits off k stars with k=2 weighted 1.5 against 1.0 for the
   others (a design choice; the only published shape frequency is Tokovinin's 4%
   for 2+2 [S]). Shuffle the leaves over the shape.
4. Root pair: separation from the truncated log-normal (section 4) up to the derived limit (the birth distribution runs to 20,000 AU).
5. Each child pair, top-down: draw an isotropic orientation and the mutual
   inclination with the parent; set `a_hi = a_parent (1 - e_parent) / (safety *
   MA01_ratio)`; draw `a` from the log-normal truncated to `[2.5 * sum of the
   children's extents, a_hi]`; draw e (circular below 7 days); require the
   Eggleton-Kiseleva test and `kl_unsafe` false; retry up to 60 times, else restart
   the hierarchy.
6. Post-check every level with the final extents (pericentre at least twice the
   children's extents; Kozai screen).

Defaults: `safety = 1.0`, Kozai screen on, inner eccentricity capped at 0.8. All
draws come from the seeded stream in a fixed order (GEN.55/GEN.56).

### 5.4 N-body test summary

REBOUND 5.2.2 (WHFast, IAS15 as checker) in the research sandbox. Failure means a
body unbound or ejected, a semimajor axis changed by more than a factor of 2, or
the inner pericentre falling below the children's physical extent. No tides,
general relativity or mass loss.

| Test | Samples | Duration | Result |
|---|---|---|---|
| Screen off (MA01 + EK95 only), N=3 to 7, one triple per level | 180 levels | 10^4 inner orbits | 30 failed (16.7%): 29 Kozai-Lidov mergers, 1 unbound; at P_out/P_in from below 10 to about 1,000 |
| Screen on, final version | 180 levels, then 600 | 10^4 inner orbits | 0 of 180; 598 of 600 stable (failure 0.33%, 95% upper bound about 1.2%); the 2 failures: an unbound pair at inclination 95 degrees and ratio 9, a merger with e_in = 0.93 |
| Tightest margins, screen on | 30 levels at 1.00 to 1.45 times the criterion | 10^5 inner orbits (mean 5,000 outer orbits, minimum 46) | 30 of 30 stable |
| Whole system, all N stars together | 30 hierarchies (6 each for N=3 to 7), total period ratio below 3,300 | 10^4 innermost orbits (3 to 8 root orbits) | 30 of 30 stable (one WHFast "unbound" at dt = P/40 was a step artefact; stable at P/200 and with IAS15) |
| Placed at k times the criterion, coplanar, random masses and e | 40 per k | 10^4 inner orbits | k = 0.6: 10 of 40 fail (25%); k = 0.8, 1.0, 1.25, 1.5: 0 of 160 fail |
| Safety factor 0.7 (30% inside the criterion) | 180 levels | 10^4 inner orbits | 3 of 180 disrupted (1.7%) |

The criterion is about 25 to 40% conservative in semimajor axis for coplanar
orbits, the right side to err on. Whole-system runs only covered hierarchies whose
total period ratio fit the compute budget; per-level tests cover the rest. Gyr
stability of high-ratio levels rests on the criterion plus the Kozai screen and
was not integrated.

### 5.5 Holman and Wiegert, implemented and tested

The repo functions and an independent implementation from the paper's
coefficients agree at all test points (equal-mass circular: S-type 0.274, P-type
2.388). Test-particle integration (1,000 binary orbits, two starting phases,
coplanar prograde, grid steps 0.03 a_b for S and 0.1 a_b for P) [C]:

| Case | Formula | Integrated boundary |
|---|---|---|
| S-type mu 0.2, e 0.3 | 0.244 | 0.244 |
| S-type mu 0.5, e 0.0 | 0.274 | all stable to 0.394, the grid edge (formula conservative) |
| S-type mu 0.5, e 0.4 | 0.147 | 0.160 |
| P-type mu 0.5, e 0.0 | 2.388 | 2.188 |
| P-type mu 0.3, e 0.3 | 3.361 | 3.461 (a stable point inside at 3.261: an island) |
| P-type mu 0.5, e 0.4 | 3.403 | 3.503 |

Within one grid step in every case and never optimistic by more than the paper's
quoted P-type error (3%). The test is coarser than the original (10^4 orbits, more
phases): a sanity check, not a re-derivation.

### 5.6 Planet zones in a hierarchy

- **S-type around star X:** upper limit `HW_S(mu = m_partner / m_pair, e_pair) *
  a_pair`, the partner being X's sibling subtree as a point mass; lower limit the
  engulfment/Roche floor (`_engulfment_radius_au`).
- **P-type around pair node P:** lower limit `HW_P(mu = lighter fraction, e_P) *
  a_P`; upper limit the S-type limit of P as a body against its own sibling in its
  parent's orbit. The window is non-empty only when a_parent / a_P is above about 8,
  so many tight hierarchies have no circumbinary zone at some levels.
- Cap everything at a disk size (100 AU; `_disk_outer_edge_au` exists).

Over 3,000 draws per N with a 0.3 AU minimum window and a crude habitable-zone
test [C]:

| N | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|
| Stars with an S-type window above 0.3 AU | 77% | 40% | 28% | 20% | 15% | 13% |
| Pair nodes with a usable P-type window (earlier estimate, see GEN.543 below) | 34% | 35% | 28% | 23% | 20% | 18% |
| Pair nodes with a P-type window above 0.3 AU (measured, GEN.543) | 45% | 48% | 43% | 37% | 35% | 30% |

About 75 to 80% of systems at every N have a habitable-zone-sized stable window
somewhere: higher N removes zones around single stars but rarely all of them.

**As built (GEN.314).** `physics/planet_zones.py` is the one place these windows are worked out: `planet_zones(nodes, members, floors_au, disk_cap_au)` over the nodes and members of a tree (`OrbitNode` and `TreeMember`, or the generator's `Node` and `Member`) gives an S-type window for every member and a P-type window for every pair node, and `flux_at(zone, distance_au, luminosities, nodes, members)` sums the light on a planet of a window's host (6.5). The S-type ceiling is the Holman and Wiegert limit or the star's Hill radius in the pair, whichever is tighter, so two stars give exactly what `Star.a_crit_au` has always held, and a pair node's lower limit is the circumbinary fit or the highest engulfment floor of a member below it, which is `StarSystem._orbit_floor_au`. The disk cap is `tuning.PLANET_ZONE_DISK_CAP_AU`. On the generator's own trees (300 draws for each N, seeded) the share of stars with an S-type window wider than 0.3 AU is 73%, 39%, 27%, 21%, 15% and 15% for N = 2 to 7, the table above within 7 points at every N. The same count for pair nodes (width above 0.3 AU alone) gives 45%, 48%, 43%, 37%, 35% and 30%, above the table's 34% to 18%: the table's "usable" P-type window was not defined past the width, so the test held only the S-type row until GEN.543. No production code calls it yet (only its tests); GEN.315 and GEN.332 are its callers.

**As built (GEN.315, slice 1).** `generation/stored_trees.py` turns a row's stored tree (`db.scatter_systems.read_stored_trees`, every column of the member and node rows) into a `StoredPlan`, and `StarSystem(member_tree=plan)` builds from it. The pair classes hold exactly two ordinary stars, so a tree builds when it is one star (no companion is drawn, whatever `BINARY_SYSTEM` and `HABITABLE_WORLD` say) or two stars (kinds main-sequence, pre-main-sequence, protostar, giant, white-dwarf) joined by one `close` or `wide` node. The second star is built from its stored seed, its birth mass and the placed star's age and [Fe/H], and must come out at its stored mass (to 1e-6), else the placed star is built alone. `BinaryStarProxy` and `WideBinaryPair` take the stored node as `orbit=` (semi-major axis, eccentricity for a wide pair, plane, phase, period) and draw nothing for it; the stored node, carried to the database's orbit epoch (`stored_trees.node_at`: the same orbit at that instant), is `StarSystem.orbit_nodes`, so `orbit_nodes` equals the stored node. The old binary columns of a close pair still track a circle of the node's semi-major axis until the readers move to the nodes. Planets, moons, belts and comets are drawn as before around the host the pair classes give (the circumbinary floor and each star's stability ceiling); the zone windows of `physics/planet_zones.py` are still GEN.332's. A tree the classes cannot hold (3 or more members, a compact or brown-dwarf member, a `compact` or `group` node) builds the placed star alone, draws no companion, and keeps its stored rows (`StarSystem.tree_built` false) for the N-star build of GEN.316; a built tree's `scatter_members` and `scatter_nodes` are deleted in the sector's save (`_insert_sector_rows`) or the uncharted-system job, in the same transaction as the link to the system. Stage B (phenomenon rows) and the `stars.role` `member` migration are not part of this slice. Generator epoch 52.

**As built (GEN.630).** The N-star host and its storage, with the fill still building single stars and pairs (no seed changes, no epoch). `StarSystem.stars` is every star; `member_stars` is the ones after the pair, and `star_slots` says for each star in `stars` order the node and the slot (0 or 1) it hangs in, `None` for a single star and for the pair on the one root node. Migration 0133 (schema v133) lets `stars.role` be `member` and adds `stars.slot`, which a stored pair gets from its roles; `db/orbit_trees.star_places` gives each star `(parent_node_id, offset, slot)`, a system with members placing every star by the whole tree at the nodes' epoch (`physics.system_tree.star_offsets_au`), and `insert_star_system` writes a `member` row for each further star. `load_star_system` reads them back into `member_stars` and `star_slots`, and the JSON form carries both (serialization version 6). A member has no planets, comets or belts of its own yet: an S-type window around it is GEN.332's. The readers that took a system to have one or two stars: the sector list orders the stars primary, secondary, then members; the system scene's tree reads each star's stored node and slot instead of guessing the free one, and a member sits at its stored offset (the pair stars keep following the binary columns the orbit update moves; a member's offset is the one stored with the system until the update moves the whole tree); `galaxy/system_position` places every star of the tree; the system page's map draws the pair (the stars table and the scene list all of them); a facility around a close pair counts the pair only; the Galaxy Map tiles treat a member like a secondary (never a bright-stars duplicate); and the body slots that an ID follows number the members by their place among them. Nothing generates a third star until GEN.631 (stages D and G draw member trees).

**As built (GEN.631).** A system the fill makes from nothing draws its stars as a tree: `StarSystem._drawn_plan` calls `multistar.draw_member_tree` for the primary at the odds of the `star-count` table (`hosts.STAR_COUNT_PROBABILITIES`), keeps the members a build cannot make from a birth mass and an age (a brown dwarf, a neutron star, a black hole, a stripped helium star) in the tree as dark members (GEN.1013: `StarSystem.dark_members`, stored in `dark_members`; `multistar.without_members` is gone), and `StarSystem._add_companions` builds the others from their own seeds (`generation/stored_trees.py`: `drawn` for a tree just drawn, `plan` for a stored one, both a `StoredPlan` of every star and its nodes numbered from 0). The pair is the two stars of the node of the placed star when its children are both stars, else of the heaviest such node, the heavier its primary; every other star is one of `member_stars` (GEN.630), hung in the tree by `star_slots`. The pair is a P-type `close` pair (circular) when its axis is under `tuning.WIDE_BINARY_SEPARATION_MIN_AU` and an S-type `wide` pair otherwise (eccentricity held under `WIDE_BINARY_ECCENTRICITY_MAX`); the other nodes keep the kind the tree drew them with. The widest axis is the place's (`place_terms.wide_axis_limit_au`, passed as `max_axis_au`). The nodes hold their phases at `stored_trees.GENERATED_EPOCH_UNIX` until the store saves the system and puts the database's orbit epoch there; a stored row's tree builds the same way, three or more stars included. The directives: `BINARY_SYSTEM` false is one star, true two or more, `WIDE_BINARY` and a `wide_binary` prevalence ask for the kind of the pair, a `binary_system` prevalence moves the share with a companion, and a protostar's or Wolf-Rayet star's own chance of a companion (GEN.242, GEN.243) replaces the table's; a tree that is not what the directives and the place need is drawn again (`tuning.COMPANION_DRAW_ATTEMPTS`). `StarSystem._companion`, `BINARY_SYSTEM_PROBABILITY_BY_SPECTRAL_CLASS` and the [Fe/H] factor on the close and inner ranges (GEN.103) are gone from the fill; the tree does not take the factor yet. (GEN.899, which moved the draw: a companion's mass was drawn from 0.1 to 1 times the primary's, so more than half of a low-mass M dwarf's companions were brown dwarfs, cut out here, and an M2V primary kept a companion star in 16% of systems against the table's 27%. Now every other star of a tree has at least the hydrogen-burning limit, 0.08 solar masses (`multistar.companion_ratio_range`), so the table's count is a count of stars, and a brown dwarf is one more member in 3% of the drawn trees, whatever the star count (`tuning.BROWN_DWARF_COMPANION_CHANCE`; a directive's count adds none). The census draw gives 26.6% for an M2V primary.) The members have no planets of their own yet (GEN.332).

**As built (GEN.332, slice 1).** The stars beyond the pair have planets. `StarSystem.set_stability_limits` works out the tree's windows (`planet_zones.planet_zones`, the member number being the star's place in `stars`) and gives each of `member_stars` the outer edge of its S-type window as `a_crit_au`, which `validation.orbit_ceiling_au` reads like a wide pair's; a close pair beside other stars gets the outer edge of its P-type window the same way, so its circumbinary planets stop short of the stars beside it (a pair that is the whole tree keeps the Hill sphere against the galaxy alone). `StarSystem._generate_member_planets` then draws each member's planets and belts exactly as a wide pair's secondary's are drawn (`_generate_planets` with no guarantees, spaced, trimmed to the ceiling, aged by the star's history), into `member_planets`, a list per member in the members' order, with `all_planets` the whole system's and `member_hosts()` the `(star, planets)` pairs. All the stars of the tree share one age, so a member's planets set a minimum for it and its state a ceiling (`adjust_pair_age_for_planets` takes the members). A member's planets are named after it (`Voranthis C I`), saved under its own `stars` row (`planets.star_id`), read back by `load_star_system` and `from_dict` (serialization version 7), listed under it by `describe`, lit and placed like a wide secondary's (`insolation.refresh_system` takes the swing from the node and the stars of the subtree beside it, `planet_zones.sibling_of`), counted in the system's totals and given a reflex offset of the star. A member has no comets yet. The system page's map still draws the pair; a member's planets are in the table and the 3D scene, which goes round the member (`systemscene._star_around`). Not yet built: circumbinary planets at the inner nodes of a tree (`planets.host_node_id`), the pair's own planets taken from the tree rather than the pair classes (GEN.640 needs that), circumbinary disks and planet pairs.

**As built (GEN.709).** The places stored for the two stars of a pair are the places the tree puts them at its epoch. A close pair is circular, so its node folds the argument of periapsis into the mean anomaly (the periapsis is 0 and the anomaly is the angle from the ascending node, which `BinaryStarProxy` tracks as its phase); a wide pair built from a node sits where `system_tree.node_relative_au` puts the secondary (the real ellipse, not the circle its phase is advanced on). The system scene gives a close pair inside a larger tree the stars' stored offsets in the whole tree, as a member gets, instead of the pair's own frame. `tests/test_gen_pair_positions.py` pins both.

**As built (GEN.710).** `StarSystem.placed_star` is the star the build was made around (member 0 of the tree). It is `stars[0]` unless a companion outweighs it (GEN.631), when the heavier takes the primary's place. `run_sector.add_star_hosted_nebulae` rolls and names the nebula from `placed_star`, the star `nebula_contents.row_nebula` rolls for a pre-placed row, so the nebula built in the sector equals the row's prediction whatever the stored members weigh. A system loaded from the database has no `placed_star`; the nebula is a sector-build step only.

**"Usable" for a P-type window (GEN.543).** The table's row was an earlier estimate with no definition past the width and a "crude habitable-zone test"; it cannot be reproduced from the text, so two definitions were tried on the same seeded draw (300 trees per N, `draw.Stream(0..299)`, main-sequence members, the host's luminosity the sum of the members below the node) against the row's 34, 35, 28, 23, 20 and 18%:

| N | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|
| (a) width above 0.3 AU and outer edge at least twice the inner edge | 38% | 41% | 36% | 31% | 29% | 25% |
| distance from the row (points) | +4 | +6 | +8 | +8 | +9 | +7 |
| (b) width above 0.3 AU and the window overlaps the pair's habitable zone (`calculate_habitable_zone` of the summed luminosity) | 7% | 15% | 15% | 13% | 13% | 11% |
| distance from the row (points) | -27 | -20 | -13 | -10 | -7 | -7 |

Neither lands within 7 points at every N ((a) misses at N = 4, 5 and 6, (b) at N = 2 to 5), so by the rule written for the item the measured numbers stay: a P-type window is "usable" when it is wider than 0.3 AU, the same test as the S-type row, and the table's row is the earlier estimate. A planet needs a window, not a window that happens to hold the habitable zone: how much light a planet gets there is `flux_at`'s job (6.5). `test_planet_zones.py` pins the measured share.


## 6. Data model, naming, maps and updates

### 6.1 Tree of pairs

A flat star list with parent pointers cannot hold a pair's orbit, which belongs to
the pair. Proposed schema (Alembic migration, MySQL 8.4 and MariaDB 11.4):

- **`orbit_nodes`** (N-1 rows per system): `id`, `star_system_id` (FK, ON DELETE
  CASCADE), `parent_node_id` (NULL for the root), `slot` ('a' or 'b'),
  `designation`, `depth`, `path`, `mass_kg`, `separation_km`, `eccentricity`,
  `inclination_deg`, `ascending_node_deg`, `arg_periapsis_deg`, `mean_anomaly_deg`
  (the phase), `period_years`, `min_update_interval_years`,
  `secondary_mass_fraction`, relative position `x/y/z_km` of child b from child a,
  `critical_s_km`, `critical_p_km`, and `kind` ('close' or 'wide', for naming).
- **`stars`**: add `parent_node_id`, `slot`, `designation`; extend `role` with
  'member'. Keep 'primary'/'secondary'/'single' for N <= 2 so queries and tests
  keep working during the migration.
- **`planets`**: add nullable `host_node_id` beside `star_id` (P-type planets orbit
  a node's barycentre, S-type a star).
- **Migration of N=2 data:** one root node per existing pair, copied from the
  `binary_*` columns; keep `star_systems.binary_*`, `is_binary` and
  `binary_configuration` as a root summary until nothing reads them. Add the table
  to `_db.ID_BLOCK_TABLES` and the cascade lists ([object-ids.md](../architecture/object-ids.md)).
- A walk to the root is at most 6 steps, so `depth` and `path` avoid recursive
  queries.

A star's offset from the system barycentre is the sum, up the tree, of
`-/+ (other child's mass / node mass) * relative_vector`. The repo does this for
one level (`binary_primary_position_*`); a loop generalizes it.

### 6.2 Eccentric pair orbits

Pair orbits are circular today on purpose. Eccentricities up to 0.9 at inner levels
need real orbits: the stability criteria depend on pericentre, and a circular
picture of an e = 0.8 orbit shows stars where they never are. Store
`mean_anomaly_deg` as the phase: it still advances linearly, so the set-based
`UPDATE ... SET phase = phase + 360 * dt / period` in `advance_orbital_phases` stays
valid and only the position step changes, using `physics/kepler.py`
(`solve_eccentric_anomaly`, `true_anomaly_and_distance_elliptical`, already used for
comets): N-1 Kepler solves per system. `physics/body_positions.positions_at` and
`static/orbitpositions.js` change together (`tests/test_js_unit.py` checks the two
agree).

As built (TEST.477, the overlap verifier): `physics/clearance.py` `audit(rows)` reads only what a finished system stores (its stars,
the nodes of its tree, planets, moons and belts, `stored_rows` reads them with plain SELECTs) and lists every place two bodies
overlap or a limit the generator promises is broken (a **hard** violation: Roche fill at pericentre, overlap of a node's two sides,
stability by Mardling and [R4574] together, node masses, a planet inside its
star, the Holman and Wiegert S-type and P-type axes, planets crossing or under 10 mutual Hill radii apart, moons outside 0.4895 of
the Hill radius or touching their planet, a belt over a planet's orbit) and the physical limits it does not promise today (an
**advisory**: the rigid-body Roche limit of a planet or moon). It imports no placement or screening code and re-derives every form;
`tests/test_clearance.py` holds each form against the generator's own copy and sweeps stored systems of one to five stars
(`PLANETGEN_CLEARANCE_SCALE=12` gave 864: one hard failure, `triple-208`, a planet pair 9.93 mutual Hill radii apart where the
generator promises 10; it is the one entry of `KNOWN`, expected to fail). Not yet checked, because the rows do not hold what it needs:
the history rule, the Kozai-Lidov screen and the 3+1 quadruple margin (GEN.1006 4.3 and section 5).

As built (MAP.255): `orbitpositions.js` walks the scene's `tree` the way `system_tree.member_offsets_au` does (the slot-1 child of a node
on the node's Kepler orbit, both children about the node's centre by the other's mass, a group node's members at its centre, the
wide pair's first star shifted to the origin) and `relativeAt` answers each star relative to the star its scene orbit names, or to
the frame's origin. A scene with stars and no tree is refused.

As built (MAP.213, slice A, the 3D system view): `relativeAt` also answers a point for each node's centre of mass (`node:<n>`, `virtual`),
stars hang from their node's centre and nodes from their parent's, `nodePaths` gives each child's real ellipse about its node's centre
(the slot-0 child opposite the slot-1 child, each at the other's share), and the compressed layout (`systemscale.js`) scales each node on
its own: the node's separation is as many units as its two sides need to clear at periapsis (a star with bodies going round it needs
their whole span, a sub-pair the reach of its own orbit), the orbit itself drawn straight-scaled so its ellipse keeps its shape. True
scale draws the real distances. Slices B (the System Map diagram) and C (the system page lists the tree) follow.

As built (MAP.213, slice C, the system page): `db.query.system_detail` carries the stored tree as `orbits` (the `orbit_nodes` rows in node order) and each star's `parent_node_no` and `slot`, and the page draws "Orbits of the stars" from `web/system_orbits_view.orbit_rows`: one row per node naming what it joins (slot 0 first, a star's name or `orbit N` for a pair hung below), with the same kind, separation, period and eccentricity columns an unbuilt scatter row's Orbits table has, plus the node's mass and its Roche `interaction`. A single star, or a system built before it stored a tree, has no panel.

### 6.3 Naming (extends GEN.62)

Outside conventions are [R], not searched: the Washington Double Star catalogue
uses A, B, C and "AB" for combined light; the Washington Multiplicity Catalog and
Tokovinin's Multiple Star Catalog add lowercase letters for subsystem members (Aa,
Ab; Ba, Bb), write a pair of subsystems as "AB,C", and use digits at the next level
(Aa1, Aa2). Rule, keeping GEN.62's output for N <= 2:

- A split whose separation is below the wide threshold (50 AU; a pair this tight
  shares one planetary system) uses letters. A split at or above it is wide: **each
  branch gets its own word** (first word shared, second drawn from the star list,
  the lighter branch from the small/child list, as GEN.62 does).
- Inside a branch a star carries a path code: level 1 `A`/`B`, level 2 `a`/`b`,
  level 3 `1`/`2`, deeper levels `.1`/`.2` appended; a lone star in a branch has an
  empty code. N=2 is unchanged.
- Example, seven stars with a wide root split, branch X = ((x1, x2), x3), branch Y
  = ((y1, y2), (y3, y4)): `Voranthis Kelmoor Aa`, `Ab`, `B`; `Voranthis Pikkita Aa`,
  `Ab`, `Ba`, `Bb`. A close-only hierarchy with one word: `Voranthis A`,
  `Voranthis Ba`, `Voranthis Bb`.
- **Planets** follow Boss's rule (decisions above): the stem is a word, never a bare
  letter. The branch holding the primary uses the system's first word, any other
  branch its own word. S-type planets around a lone star in a branch: `Voranthis I`
  and `Pikkita I`. S-type planets around a coded star add the code (`Pikkita Ba I`),
  because one word would be ambiguous inside a branch; P-type planets around a node
  take the node's name (the branch word at the branch root, else `<branch word>
  <code>`). The code-bearing forms are a recommendation beyond Boss's two-star rule.
  The structure is displayed as a bracketed string such as `(Aa,Ab),B`.
- Branch words are drawn after the structure from the seeded stream in a fixed
  order (GEN.55/GEN.56), as `_draw_star_words` in `system.py` does for two words.
  GEN.71 touches only the system name, so there is no conflict.

### 6.4 Maps and the orbital update

- **Scene:** `around` becomes `node:<id>` as well as `star:<id>` and `barycenter`;
  the scene lists nodes with their orbits and every star's parent node.
  `_planet_center` in `systemscene.py` returns `star:<id>` or `node:<id>` by host.
- **Scale:** a seven-star system spans 0.01 AU to 10^4 AU, six orders of magnitude,
  which one linear map cannot show. The system map needs a per-node focus (click a
  pair to zoom into its two-body orbit at its own scale) and a log-radial overview;
  the bracketed structure string is the quick index.
- **Update:** each pair is an analytic Kepler orbit propagated by the clock, as
  `positions_at` does for a close pair. No mutual perturbation is needed for a
  stable hierarchy over game time spans (secular precession and Kozai-Lidov are
  ignored, a documented simplification).
- **GEN.109:** sector point-mass tables hold one entry per system at its barycentre
  with the total mass (section 2). `BinaryStarProxy._calculate_system_perimeter_static`
  should use the total mass.
- **Plumbing:** migrate N=2 through the new tables first, so the old code paths are
  proven equal before any N >= 3 data exists.

### 6.5 Habitability inputs (GEN.84 to GEN.89)

Insolation from several stars is a sum, `S = sum(L_i / d_i^2)`, with the
time-average over an eccentric orbit `<1/d^2> = 1 / (a^2 sqrt(1 - e^2))`. A
companion at 50 AU adds about 4e-4 of its luminosity to a planet at 1 AU, so for
S-type planets only a companion inside a few times the planet's orbit matters;
P-type planets see the combined luminosity with a flux that varies with the inner
pair's period. GEN.84's score needs a `flux_variation` input; GEN.86 (activity and
flares) is per star and the flux is summed.

**As built (GEN.148).** `physics/insolation.py` sets `flux_variation` on every planet
and moon (a gas giant keeps none, its moons take its value): the peak-to-peak swing of
the planet's total light as a fraction of its mean, 0 for a lone star. For a planet of
one star of a wide pair it is the companion's `L_c (1/periapsis^2 - 1/apoapsis^2)` over
the mean total `L_h/d^2 + L_c/(a^2 sqrt(1 - e^2))`; for a circumbinary planet it is
`4 delta / d` with `delta = |L1 r1 - L2 r2| / (L1 + L2)` the luminosity centroid's offset
from the barycentre (`r1 = a m2/M`, `r2 = a m1/M`). The Temperature domain scores the
worse of the cold and hot extremes `T (1 -+ f/2)^(1/4)` of the cycle, and the wet-bulb
test reads the hot one. The equilibrium temperature itself still comes from the host star's
light alone: the companion's mean share is at most a few percent of a planet's light
(under 0.1% for a companion beyond 100 AU) and a planetGen simplification kept for now.
Over 200 generated systems the median swing of a wide pair's planet is about 5e-7 (the
largest 0.23, a close companion) and of a circumbinary planet 0.05 (the largest 0.46).

### 6.6 Placement by population (GEN.103)

Multiplicity depends on the stellar population: the close-binary fraction rises
sharply as metallicity falls ([R], so more close pairs in the halo and thick disk),
massive multiples sit in young regions (arms, OB associations, clusters), and wide
pairs survive less in dense regions (bulge, clusters, halo), so the wide-pair limit (the derived limit that replaced the fixed 10,000 AU cap, GEN.448)
shrinks there. This is rate shaping on the population mix of
[galaxy-disk-density.md](../galaxy/galaxy-disk-density.md) (its section 4) and does not touch
GEN.100's densities.

**Research supplied (Boss, 2026-10-16, 5:55 PM EDT, for GEN.103).** Two constants Core 4 could not look up
arrived as research, and they are written here as rules the type draw derives, not as values it stores. The
sources Boss named are [R0223, R0229, R0281, R1235, R1880, R2166, R2474, R3109, R4204]. **The reference list is unverified; Core 4 checks each source and each
number before any of it is coded.**

1. **Close-binary fraction against metallicity.** For solar-type primaries (about 0.6 to 1.2 Msun) the fraction
   inside about 10 AU (P below about 10^3.5 days) falls as [Fe/H] rises: about 50% to 60% at [Fe/H] below -1
   and about 10% to 15% above +0.3 (the table in section 2 gives the finer steps; Core 4 reads the [R3109]
   points as 53%, 40%, 24% and 10% at [Fe/H] -3.0, -1.0, -0.2 and +0.5, and El-Badry and Rix as 2019, so the
   rounder figures here are summaries to be replaced by the points). The dependence vanishes
   for separations above about 200 AU. The cause is disk fragmentation: low metallicity means low dust
   opacity, so the outer disk cools fast (t_cool below about 3 over Omega, with Toomre Q near 1 or below),
   fragments into companions that migrate inward and stall inside about 10 AU, while metal-rich disks
   suppress fragmentation. **Rule:** the close-companion share of a draw is a function of the cell's [Fe/H]
   (from the population mix) and of nothing else, and the wide-companion share is not changed by [Fe/H].
2. **Encounter-disruption limit for wide pairs.** A wide pair loses its outer end to flybys and clouds
   (Heggie's law: soft binaries widen and dissociate) and, as the absolute outer envelope, to the Galactic
   tide. **Rule:** the survival limit is computed from the cell's own environment, never read from a table.
   (a) The diffusive limit: the largest separation whose binding energy survives the summed impulses of
   perturbers over the age of the population, from the local stellar number density, the mean perturber
   mass, the velocity dispersion, the Coulomb logarithm and the age ([R0229, R2166]; Core 4 found the formula as Boss typed it dimensionally inconsistent, so it is re-derived from those papers before it is coded). (b) The tidal limit:
   the Jacobi radius from the galactic potential at the cell, (G M / (4 Omega^2 - kappa^2))^(1/3) = (G M /
   (2 A (A - B)))^(1/3) for a pair of total mass M, which the generator evaluates from the same potential
   as the orbital update. A pair is kept only below the smaller of the two. **Check values, not inputs:** the
   solar neighbourhood should give a diffusive limit near 2 to 5 x 10^4 AU and a Jacobi radius near 1.5 to
   2 pc for 2 Msun, with Gaia wide-binary catalogues truncating near 0.1 to 0.25 pc. Where GEN.103 lands, this
   derived limit takes the place of the fixed 10,000 AU cap, and it shrinks in the bulge, clusters and halo
   because the density and dispersion there are higher.

**Boss's decisions (2026-10-16, 6:12 to 6:15 PM EDT, by cards in the project chat).** (1) The close-binary trend is
the pooled one, read at the star's own [Fe/H]. (2) The metallicity effect applies to every star mass, and the
generator is checked that it behaves as expected: TEST.196 compares the generated close-companion share, by primary-mass
bin and [Fe/H] bin, to the published ranges and flags any bin outside them. (3) Wide pairs have no metallicity
dependence. (4) The encounter limit is a sharp cut at a_1/2, the smaller of the two limits above. (5) The fixed
10,000 AU wide-pair cap is replaced by the derived limit. The last card, whether the metallicity effect fades above about 2 solar masses, was answered "all masses": there is no fade, and the check is the safeguard. A last card (2026-10-16, 7:19 PM EDT) chose the form: the factor F/F(0) is applied to the odds p/(1-p) of a close companion, not to the chance p, for every mass. It cannot pass 1, and a star already likely to have a companion rises less.

## 7. Compact-object primaries

### 7.1 How they arise

A compact remnant is the end state of one star of a system (white dwarf below about
8 Msun, neutron star or black hole above). When the progenitor was in a tree, the
rest stays bound only if the supernova's mass loss and natal kick do not unbind it.
A Monte Carlo of an instantaneous supernova on a circular orbit (Blaauw/Brandt and
Podsiadlowski style; kick and masses are assumptions) gives the bound fraction [C]:

| Case (pre-SN masses, companion) | 0.1 AU | 1 AU | 10 AU | 100 AU | 1,000 AU |
|---|---|---|---|---|---|
| NS, 12 to 1.4 Msun, companion 3 Msun, kick sigma 265 km/s (Hobbs-like [R]) | 10.8% | 0.8% | 0 | 0 | 0 |
| Same, low kick (sigma 30 km/s: electron capture or ultra-stripped) | 0.1% | 11.7% | 8.8% | 0.6% | 0 |
| BH, 30 to 10 Msun, companion 5 Msun, kick sigma 37 km/s (momentum conserving, 265 x 1.4/10) | 11.6% | 28.5% | 19.2% | 2.0% | 0.1% |
| Same, no mass-loss kick | 0 | 0 | 0 | 0 | 0 |

Wide orbits around fresh remnants almost never survive, and a collapse that loses
more than half the mass unbinds everything. Observed surviving binaries need an
earlier mass-transfer phase this generator does not model, so the table explains
why the companion fractions below are small; it is not a population synthesis. A
white dwarf has no kick, so its binary survives with the orbit widened by roughly
the mass-loss ratio.

### 7.2 What real systems show

| Population | Facts | Tag |
|---|---|---|
| Gaia BH1 | 9.62 +- 0.18 Msun hole, 0.93 Msun G dwarf, 480 pc; P = 185.6 d, e about 0.45, a = 1.40 AU (Kepler's third law) | [S] masses; [R] P and e |
| Gaia BH2, BH3 | BH2 about 9 Msun, P about 1,300 d, a about 5 AU; BH3 32.7 Msun with a 0.76 Msun giant, P about 4,000 d, e about 0.72 | [S] periods; [R] BH2 mass, orbit sizes |
| Dormant-BH puzzle | orbits too tight for non-interacting stars, too wide for common-envelope products | [S] |
| Pulsars in binaries | young: observed 1.2 to 1.6%, corrected below about 5% (conservative 8%); millisecond: about 79%; neutron stars born in binaries that end isolated: about 86% | [S] arxiv.org/pdf/2011.08075 and summaries |
| White dwarfs | about 22% of field white dwarfs in wide binaries; 25 to 50% of cool ones show metals (pollution) | [S] |
| Pulsar planets | PSR B1257+12: three planets at about 0.19, 0.36, 0.46 AU, about 0.02, 4.3, 3.9 Earth masses; PSR B1620-26 b in M4 (probably captured); PSR J1719-1438 b (about 1 Jupiter mass, stripped white dwarf companion); about 1% of millisecond pulsars have planets, 0 of 151 young pulsars surveyed | [S] existence and rates; [R] B1257 orbits and masses |
| WD 1856+534 b | transiting, P = 1.4 d, Jupiter-sized, at most about 6 Jupiter masses, about 186 K | [S]; host in a wide triple [R] |
| Compact triple | PSR J0337+1715: millisecond pulsar with 1.6 d and 327 d white dwarfs (ratio about 200) | [R] [R3664] |

### 7.3 Companion fractions, conditional on a compact primary

GEN.100 supplies how many neutron stars and black holes exist
([interstellar-object-rates.md](interstellar-object-rates.md)); the fractions below
say what company they keep. They are rounded design defaults from the table above,
for Boss's review. X-ray binaries get no separate rate, only the conditional
fraction applied to GEN.100's scattered objects.

| Primary | No stellar companion | Bound stellar companion | Typical companion |
|---|---|---|---|
| White dwarf | about 70% | about 30% (22% wide plus a few % close or double white dwarf) | main-sequence star at 1 to 10^4 AU; double white dwarf; post-common-envelope pair (hours to days) |
| Neutron star, young | 97 to 99% | 1 to 3% | Be or main-sequence star; wide orbits unlikely (kick) |
| Neutron star, recycled millisecond pulsar | 20 to 25% (field), lower in clusters | about 75 to 80% | low-mass white dwarf (days to years), "black widow" or "redback" of a few hundredths of a solar mass, rare double neutron star |
| Neutron star or black hole, active X-ray binary | n/a | the share of the above still accreting | low-mass donor at a period of hours to a day |
| Stellar black hole | 70 to 90% [R] | 10 to 30% [R] | Gaia-type dormant pair, a = 1 to 20 AU, Sun-like or giant star |
| Black hole or neutron star inside a triple | a few % of those with companions [R] | hierarchical, outer period at least 5 times larger | the companion's own inner pair (PSR J0337+1715) |

The black hole fraction rests on the three Gaia systems and is low confidence.

### 7.4 Radiation, tidal and size limits [C] unless marked

- **Black hole (10 Msun):** gravitational radius 14.8 km; innermost stable circular
  orbit 89 km non-rotating, about 18 km maximally spinning prograde [R, Kerr factor].
  Earth-density tidal-disruption radius about 9.5e5 km (0.0064 AU); for 33 Msun,
  1.4e6 km (0.0095 AU). Eddington luminosity 1.26e31 W per solar mass [R].
- **Accretion:** `L = 0.1 * Mdot * c^2`. A low-mass X-ray binary at 1e-9 Msun/yr
  gives 5.7e29 W, 2.0e6 W/m^2 at 1 AU (about 1,500 times the solar constant, mostly
  X-rays). A dormant pair (Gaia BH1) has none.
- **Neutron star / pulsar:** spin-down power `L = 4 pi^2 I Pdot / P^3` with I = 1e38
  kg m^2 [R]. PSR B1257+12 (P = 6.22 ms, Pdot = 1.14e-19): 1.9e27 W, 5.1e4 W/m^2 at
  0.36 AU if absorbed isotropically (38 times the solar constant). Most of it is
  relativistic wind and field, not heat, but lethal and atmosphere-stripping. A
  typical millisecond pulsar is similar (1.5e27 W); the Crab is 4.5e31 W; an old 1 s
  pulsar 4e24 W. Earth-density tidal radius for 1.4 Msun: 4.9e5 km (0.0033 AU).
- **White dwarf (0.6 Msun, R = 8,800 km, 4.2e8 kg/m^3):** Roche limit for a 5.5
  g/cm^3 rocky body 0.0031 AU rigid, 0.0061 AU fluid (3 g/cm^3: 0.0038 and 0.0074).
  Habitable zone with flux limits 1.1 and 0.35 of solar: 0.036 to 0.064 AU at
  Teff 10,000 K; 0.013 to 0.023 at 6,000 K; 0.0058 to 0.010 at 4,000 K; 0.0032 to
  0.0058 at 3,000 K. The zone moves inward as the star cools, reaches the Roche
  limit near 3,000 to 4,000 K, and has periods of 0.3 to 5 days, so planets there are
  tidally locked ([R0038, R1898]: continuously habitable within
  about 0.02 AU, locking in 10 to 1,000 years [S]).
- **Gaia-type stable zones** (Holman and Wiegert on the measured orbits, mass
  fraction clamped to the fit's range): BH1's Sun-like star keeps planets only
  inside 0.097 AU (its habitable zone is about 0.8 to 1.5 AU), circumbinary planets
  only beyond 5.1 AU. BH2's star: inside 0.31 AU, beyond 19 AU. BH3's star: 0.52 AU
  and 69 AU (approximate, the P-type fraction clamped to 0.1). The generator
  reproduces the literature: these companions cannot host habitable worlds.

### 7.5 What a planet around a black hole can be

Observationally nothing: no planet has been found around a stellar black hole.
Theory allows (a) planets of the companion star (S-type, inside the limit), (b)
circumbinary planets beyond the P-type limit, (c) captured rogues (negligible rate),
and (d) planets formed from a post-supernova fallback or debris disk. Case (d) is the
pulsar-planet analogy: PSR B1257+12 supports it for neutron stars and nothing
supports it for black holes. A rogue inside a hole's tidal radius is torn apart. The
generator puts (a) and (b) around black holes only through the companion and leaves
(d) as a hand-forced option.

### 7.6 Rules for the generator

| Primary | Put around it | Skip |
|---|---|---|
| **White dwarf** | normal planet draw with a floor at the fluid Roche limit; surviving outer planets expanded by the progenitor-to-remnant mass ratio (1 Msun to 0.54 Msun grows outer orbits by 1.85; the law `a (M_star + M_p) = constant` is verified for adiabatic mass loss, see [exotic-environments-planets-and-compact-binaries.md](exotic-environments-planets-and-compact-binaries.md) 1.2); engulfment floor from the giant phase (`_engulfment_radius_au`); optional debris belt; an extremely rare close giant (WD 1856+534 b type); pollution flag for a quarter to a half | habitable worlds inside the Roche limit; life on a planet younger than the zone's lifetime |
| **Neutron star, young** | nothing | all planets; the star's own fields only |
| **Millisecond pulsar** | 1% of them get 1 to 3 small planets at 0.1 to 0.5 AU (PSR B1257+12 type) or an ablated companion (black widow, J1719-1438 type); optionally a captured outer planet in a cluster | any habitable score above the lowest tier (radiation); atmospheres on short-period planets |
| **Stellar black hole, dormant** | the companion star with its S-type zone from Holman and Wiegert (often only inside 0.1 to 0.5 AU) and a P-type zone far out | planets in the hole's debris disk; a habitable zone around the hole |
| **Active X-ray binary** | a close donor, an accretion disk with `L = eta Mdot c^2`, eclipse/outburst flags | planets inside the disk; habitable zones (X-ray flux of 1e3 solar and more) |
| **Supermassive black hole** (`mass_class = "supermassive"`) | the galactic centre's S-star cluster as a phenomenon description only | any star system |

For the tree, a remnant is a leaf whose class is NS, BH or WD; mass and radius follow
the evolution rules, and the bound-companion probability comes from 7.1 and 7.3. The
hierarchy test (section 5) applies after the supernova, because the remaining orbits
must satisfy the criteria with the remnant's mass. A cheap rule that keeps draws
valid: draw the pre-collapse tree, evolve the leaves to the system's age, keep the
tree if a mass-loss and kick check (7.1 style) passes, else drop the leaf's
companions (the star becomes a free compact object, which GEN.100 already scatters).
The pulsing fraction and its age rule are in [anomalies.md](anomalies.md); that rule
keeps a floor for recycled pulsars in binaries, which is where the millisecond-pulsar
row here applies.

## 8. Go / no-go and build order

Boss's research of 2026-10-09 18:23Z refines the GEN.130 answer below (millisecond-pulsar planet rate and types, wide tertiaries only for direct-collapse holes, a gravitational-wave age check for compact binaries, Kozai screen with relativistic quenching): [exotic-environments-planets-and-compact-binaries.md](exotic-environments-planets-and-compact-binaries.md) sections 2 and 3.

**GEN.129: go.** The value is clear, the physics is tractable and verified, the risk
is plumbing. Do not start by adding code for N=7; start with the data model at N=2.

**GEN.130: go in three slices, no-go on four items.**

- Go: (a) allow collapsed leaves in the tree and a kick/survival check (a small
  change in `star.py`, where collapsed primaries are redrawn and companions
  lightened); (b) NS/BH/WD presets with the rules in 7.6, reusing `BlackHole` and
  `NeutronStar`; (c) radiation fields (spin-down power, X-ray luminosity) feeding
  GEN.86 and GEN.87.
- No-go: young-pulsar planets; debris-disk planets around black holes as a default;
  habitable scoring above the lowest tier around any remnant except a cool white
  dwarf; supermassive-hole stellar systems. Out of scope: kilonovae, mergers, and
  anything needing a binary-evolution code.

**Build order**

1. Redraw the existing binaries from one log-normal period law (fills the 0.26 to 50
   AU gap), eccentric pair orbits with `physics/kepler.py`, mass-ratio refinements.
   No schema change; tests the P-type/S-type switch by limits.
2. `orbit_nodes` and the N=2 migration; `stars.parent_node_id`;
   `planets.host_node_id`; readers and writers (`store.insert_star_system`,
   `_db.load_star_system`, `serialization`), keeping the old columns filled.
3. Naming (6.3), scene, maps (per-node zoom), `positions_at` and its JavaScript twin,
   `advance_orbital_phases` for nodes.
4. N=3 generator (the algorithm of 5.3). Tests: criterion unit tests; a property
   test (every generated tree passes MA01, EK95 and the Kozai screen, no pair inside
   another's extent); a small leapfrog or `DOP853` regression for a few seeded
   triples. REBOUND is GPL-3.0-only and the project is CC0, so the heavy N-body
   campaign stays an unshipped dev script
   ([orbital-solvers-and-integrators.md](../orbits/orbital-solvers-and-integrators.md), licence
   decision).
5. N=4 to 7: the same recursion, with tuning of shape weights and tail
   probabilities; map and search UI for deep trees.
6. Planet zones from 5.6, flux sums, GEN.84 to GEN.89 inputs.
7. GEN.130 slices (a) to (c).

## Evidence notes

The research environment could only read search-result text, not the papers. Claims
still [R], to verify when paper access is allowed:

- Log-normal period law for solar-type stars (mean 5.03, sigma 2.28 in log10 days);
  every period share in section 4 depends on it.
- Close-binary fraction against metallicity [R3109].
- Per-mass-bin companion, triple and quadruple fractions (Moe and Di Stefano table,): the table in 3.3 is a design blend, and its A, B and O rows and every
  N >= 4 split are the least certain.
- Tokovinin's "ratio above 5" rule; the WDS/WMC/MSC notation (6.3).
- Gaia BH1 period and eccentricity, BH2 mass and orbit, BH3 beyond the periods;
  PSR B1257+12 masses and orbits; PSR J0337+1715.
- Moment of inertia 1e38 kg m^2, Kerr ISCO factor, Eddington 1.26e31 W per Msun,
  neutron star kick sigma 265 km/s, white dwarf orbit expansion on mass loss.
- Black hole companion fraction (10 to 30%) and the compact-primary-in-a-triple
  fractions in 7.3: low confidence.
- Whether AR Cassiopeiae is septuple; whether Mizar and Alcor are bound.

The research scripts (generator and criteria `hier.py`, N-body drivers, Holman and
Wiegert and boundary scans, zone and table scripts) lived in a session scratchpad
outside the repository; GEN.129 step 4 recreates the generator from 5.1 to 5.3.

## Sources

Seen in search results (summaries only; no paper text was read):

- [R3646]: [R3646]; 2015 reanalysis
  [R4693]
- [R1175]: [R1175]
- [R3108]: [R3108, R3334]:
  [R3333]
- [R4384]: [R4384]; 2018 catalogue
  [R4385]; 2021 review [R4386]
- [R1223]: [R1223]; Gebrehiwot et al.
  2016: [R1529]
- Nu Scorpii and Mizar/Alcor (Wikipedia); HD 98800
  [R3603]; TIC 168789840
  [R4989]
- Stability: [R1537] (Eggleton and Kiseleva form); Mardling
  and Aarseth forms quoted in [R1855, R2519, R3864]; Holman and
  Wiegert coefficients quoted in [R0243, R3631]; Quarles et al. refit
  [R3633]
- Gaia black holes: [R1238] (BH1),
  [R2091] (BH3)
- Pulsars: [R0126, R2078, R2849, R2902] (natal kicks);
  X-ray binaries: [R0893, R2797]
- White dwarfs: [R0038, R4493],
  https://astrobites.org/2020/09/22/the-first-planet-found-orbiting-a-white-dwarf/,
  https://doaj.org/article/6ee22a7610b842f98a9ddef9f43c063e (WD 1856+534 b),
  [R3306] (pollution)
- Package page: https://pypi.org/pypi/rebound/json (rebound 5.2.2)
- Repo files read: `generation/binary.py`, `wide_binary.py`, `system.py`, `star.py`,
  `phenomena/compact_remnant.py`, `physics/orbits.py`, `physics/kepler.py`,
  `physics/stellar_evolution.py`, `names/bodies.py`, `tuning.py`,
  `web/maps/systemscene.py`, `docs/reference/database-schema.md`, and the design documents
  linked above, with "Orbital Update Full Algorithm.md" and "Orbital Position and
  Vector Update Algorithms.md"; TODO items GEN.128 to GEN.130, GEN.109, GEN.84 to
  GEN.89, GEN.103; GitHub issues #777 and #778.
