# Captured asteroid moons (irregular satellites): physics and a generator design

Status: plan (Boss said yes 4:27 PM EDT 2026-10-30 and confirmed both recommendations of section 8 at 4:34 PM EDT: capture windows only, and the belt storage rule with a 1e16 kg floor; nothing built). Written 4:40 PM EDT 2026-10-30 as the lane reported. No generator code changed. Asked for by Boss (4:27 PM EDT 2026-10-30): "Jupiter has so many moons because a lot of them are former asteroids from the belt; if gravitational bodies are near a belt, they should pick up asteroidal moons." Ties to `belt-evolution-and-body-ledger-2026-10-30.md` (the ledger, member IDs `<belt ID>.<rank hex>`, the cut rule E) and to the moon code (`generation/planets.generate_moons`, `physics/moon_composition.py`, `physics/tuning.MOON_*`).

Evidence tags: [S] read this session in the paper text (arXiv or the authors' PDF); [R] recalled, **not opened this session**, the building lane reads it before the number is used; [REPO] read in the repo; [C] computed here. Everything [R] or invented is a stand-in for GEN.725 with tier and tweak level.

Papers opened this session: Jewitt and Haghighipour 2007 (ARA&A 45, 261; arXiv astro-ph/0703059) [S]; Nesvorny, Vokrouhlicky and Deienno 2014 (ApJ 784, 22; arXiv 1401.0253) [S]; Vokrouhlicky, Nesvorny and Levison 2008 (AJ 136, 1463) [S]; Sharkey et al. 2025 (JWST, arXiv 2501.16484) abstract [S]; Bhatt et al. 2017 (arXiv 1710.06200) title and authors only.

## 1. The short answer to Boss's idea

The idea is half right, and the half that is wrong changes the design.

**Right:** most of the "many moons" of the giants are small, captured, irregular bodies, not formed in a disc round the planet. About 100 irregulars were known at the four giants by 2007, and they outnumber the regular moons (53 regular, 159 in all, Table 2 of Jewitt and Haghighipour) [S]. The irregulars are small-body-reservoir material: Jupiter's are spectrally Trojan-like or asteroid-like (Sharkey 2025: Himalia looks like an ammoniated, Ceres-like asteroid; Carme, Sinope and Themisto look like "red" Trojans) [S].

**Wrong, in its simple form:** a planet standing beside a belt does not pick up that belt's bodies. A body that comes into a planet's Hill sphere from outside leaves again unless something takes its energy away, and in a gas-free system nothing does. Nesvorny et al. 2014 give a capture probability of only (1.3 to 3.6) x 10^-8 per disc particle even in the best (instability) case, so a belt of a few hundred big bodies gives about 1e-5 captures [C]. What actually put the Solar System's irregulars where they are, in the leading model, was the **early instability**: Jupiter and Saturn met an ice giant in close encounters while a massive outer planetesimal disc was still there, and the three-body interaction (planet, ice giant, planetesimal) left some planetesimals bound [S, NVM07 and NDD14]. The source was the outer (transplanetary) disc, not the main belt.

**So the rule that fits both the physics and Boss's picture:** *a planet gets irregular moons from the small-body reservoirs it interacted with while it could dissipate energy.* The reservoirs are the system's belts and discs **at the time of a capture window**, and the windows are (1) the gas era (gas drag, pull-down), (2) an instability with giant-giant encounters (three-body capture), (3) an exchange with a binary planetesimal, (4) a collision with an existing moon. A belt near a planet matters through its **mass at the window** (a primordial main belt held 100 to 1000 times its present mass [R, Walsh 2011; Raymond and Izidoro 2017], so the belt-evolution design's early-depletion step is exactly the sweep that feeds captures), and through the **composition** it gives the captured bodies. Without a window, a planet beside a belt gets none, and that is also physically right (Earth, Venus, Mercury have none; Mars's Phobos and Deimos are now thought to be impact-disc moons, Rosenblatt 2016 [R]).

## 2. Physics

### 2.1 Where a moon can stay: the Hill sphere
* `r_H = a_p (m_p / 3 M_star)^(1/3)` [S, JH07 eq. 1]. Giants: 0.35, 0.43, 0.47, 0.77 AU for Jupiter, Saturn, Uranus, Neptune [S, Table 1].
* Stable limits (time-averaged), from numerical orbit surveys: prograde about **0.4 r_H** (the repo's `MOON_PROGRADE_STABLE_HILL_FRACTION = 0.4895` is the Hamilton and Burns 1991 figure [REPO, R]), retrograde about **0.5 to 0.6 r_H**, depending on the mean eccentricity (Nesvorny et al. 2003, quoted in VNL08) [S]. Retrograde orbits are stable out to larger distance, so after any loss process retrograde bodies outnumber prograde: **88 retrograde to 19 prograde, 4.5 to 1** [S, JH07 Table 2]. No known observation bias can produce that. Capture itself is symmetric in inclination [S].
* **Kozai band:** no irregular is found with inclination 55 to 130 degrees [S]. The Sun drives the eccentricity up on a timescale of 65 to 180 yr at Jupiter and the orbit hits the planet or the regular moons. The repo already encodes the band (`KOZAI_FORBIDDEN_DEG = (55, 130)`, `draw_inclination_deg`) [REPO].
* **Inner edge:** irregulars with pericentre `q = a(1-e)` inside the outermost regular moon are removed (collisions): q > 0.013 AU at Jupiter (Callisto), 0.023 at Saturn (Iapetus), 0.004 at Uranus (Oberon) [S, VNL08].
* Median `a / r_H` of the observed irregulars: **0.44, 0.29, 0.17, 0.19** for Jupiter, Saturn, Uranus, Neptune [S, JH07]; the outer half of the Hill sphere is empty. Captured orbits start at 0.03 to 0.14 AU prograde, 0.03 to 0.2 AU retrograde at Jupiter, eccentricity 0 to 0.7 [S, NDD14].
* **Scale for the generator [C]:** an Earth-mass planet at 1 AU has `r_H` = 1.5 million km, a Mars-like one 1.1 million km, a super-Earth 2.6 million km, Jupiter 53 million km, a hot Jupiter at 0.05 AU only 0.5 million km (7 Jupiter radii, so almost no irregular zone). A giant round a 0.3 solar-mass star at 1.5 AU: 23 million km.

### 2.2 Capture mechanisms and when each applies
| Mechanism | Needs | Efficiency | Applies to | Source |
|---|---|---|---|---|
| **Gas drag** in the planet's extended envelope or circumplanetary disc | gas giant, young, a dense extended envelope | high while it lasts, window of Myr | gas giants only; ends when the gas goes | Pollack 1979 [R]; JH07 sec. 6 [S] |
| **Pull-down** | planet mass grows fast (runaway), Hill sphere jumps, bodies in the old Hill sphere are bound | works for retrograde and prograde | gas giants only (needs runaway growth) | Heppenheimer and Porco 1977 [R]; JH07 sec. 6.1 [S] |
| **Three-body capture in planetary encounters** | an instability: two or more planets meeting while a planetesimal disc is present; the closest pass, the geometry | (1.3 to 3.6) x 10^-8 per disc particle at Jupiter | every giant; ice giants need it (gas drag and pull-down cannot reach them) | NVM07 [R], NDD14 [S] |
| **Exchange / binary dissociation** | a planetesimal binary flying past; one component keeps the energy | far too low for small irregulars, but the only route found for **Triton** (a binary of two to five Triton masses, early, with the gas still present) | one massive capture per ice giant at most | Agnor and Hamilton 2006 [R]; VNL08 [S] |
| **Collision with a regular moon** or tidal friction | an existing large moon, a pass through its orbit | low | sets the inner edge, destroys more than it makes | Colombo and Franklin 1971 [R]; JH07 [S] |
| **Collisional fragmentation after capture** | a captured body, time | not a capture, a multiplier: one body becomes a family | all | Bottke 2010 [R]; JH07 sec. 5 [S] |
| **Temporary capture** | a body passes through the Hill sphere at low speed (Earth always holds a metre-size mini-moon) | the body leaves again after months to years | every planet, continuous | Granvik 2012 [R] |

Jewitt and Haghighipour's observation that the four giants hold nearly the same irregular populations despite masses of 17 to 310 Earth masses (JH07 sec. 3) points at the mechanism whose efficiency does not depend on the planet's accretion history: three-body capture, helped by the larger Hill spheres of the distant planets [S].

### 2.3 How many (the numbers a generator can use)
* Disc: a transplanetary planetesimal disc of (3 to 4) x 10^7 bodies above R = 40 km (calibrated from the Jupiter Trojans), cumulative size law index **5 above R = 70 km and 1.8 below** [S, NDD14].
* Capture probability per disc particle at Jupiter in three instability histories: **(1.3 to 3.6) x 10^-8**; so about **one irregular above R = 40 km** is expected at Jupiter, against the observed Himalia (R about 85 km) and Elara (43 km) [S, NDD14]. Fraction of captured moons that survive on the orbits 0.05 to 0.1 AU: 0.01 to 0.3 depending on the history [S].
* The same process captured a **much larger population of small bodies** than the known R < 10 km irregulars: either the survey is incomplete, or collisions removed them, or the size law has a wave below 10 km [S, NDD14]. The observed law is differential q = 2 (cumulative index 1) at Jupiter, Saturn and Uranus above about 5 km and steeper (q about 3.5, cumulative 2.5) below 5 km at Jupiter [S, JH07]. A flatter law than the source's q of about 4 means the final population is **not** a direct copy of the source.
* NDD14 confirm NVM07's capture at Saturn, Uranus and Neptune, and find Jupiter's probability 1.5 to 4.2 times that of the one isolated NVM07 case [S]; the per-planet probabilities themselves are in NDD14 Table 1, which the building lane reads (item 1).
* Collision timescales: the irregulars of every giant would be destroyed by mutual collisions within the age of the Solar System if they were at Jupiter-like space density, so "the final distribution is shaped by capture, long-term instability, collisions with regular moons, and collisions among the irregulars" [S, NDD14 sec. 3.2; JH07 sec. 4 and 5]. The invariance of the final populations is, in my reading, the signature of **collisional saturation**: a population ground down as `N(t) = N0 / (1 + t / t_c)` ends up independent of `N0` for `t >> t_c`. That is the same law the repo already has for belts (`belt_mass_earth`, Wyatt 2007, [REPO]). This is an educated guess from physics, tier "Marginal educated guess from physics", to be tested in item 5 below.
* Rough observed head-counts for calibration [R, not opened]: irregulars with R > 5 km number about 10 at each giant (Jupiter 15, Saturn 9, Uranus 8, Neptune 7 within a factor of two). The calibration set the lane must read first: Jewitt and Haghighipour Table 2 [S: totals 63, 56, 27, 13; irregulars prograde/retrograde 6/49, 8/27, 1/8, 4/4 by 2007], Sheppard and Jewitt 2003 and Sheppard et al. 2005/2006 for the size laws.
* Population of the same capture family that are **not** moons: Jupiter Trojans (the same disc, capture probability about 6e-7 per particle, [C] from the 3 to 4e7 calibration, assuming 25 Trojans above R = 40 km [R]), so Trojans outnumber irregulars about 20 to 50 times above 40 km. Useful as a test target and as the generator's co-orbital population (belt-evolution doc sec. 6, "Trojan swarms").

### 2.4 Families
Irregulars cluster in `a/r_H` versus inclination (JH07 Fig. 3): Jupiter's major families are Himalia (Himalia, Elara, Lysithea, Leda), Carme, Ananke, Pasiphae, plus loners; Saturn has the Gallic, Inuit and Norse (Phoebe-linked) groups; the four giants hold "a handful of families" each [S, JH07 sec. 3; names recalled from Sheppard and Jewitt 2003 [R]]. A family is the debris of one disruptive collision after capture, with a velocity spread of tens to a few hundred m/s [R]. Himalia's group is spectrally non-uniform (Himalia, Elara and Lysithea differ in phyllosilicate content), so the parent was heterogeneous [S, Sharkey 2025]. Sources: Nesvorny et al. 2003 and 2004 [R]; Bottke et al. 2010 [R].

## 3. What this means for the existing generator [REPO]

* `generate_moons` places moons log-uniformly out to `0.4895 r_H`, each with a draw from `ORIGIN_TABLE` whose weights are **open** (the rows a moon supports share the probability equally): for a giant, disk and captured (and ring when a ring is near) get a third each regardless of the moon's mass. That makes a captured moon of 1e20 kg or more (R of 260 km and up) as common as a disc moon, but real captured bodies that large are Triton, Phoebe-size at most; N(> 100 km) is about 0.1 to 0.5 per giant, not one in three [S/C].
* `MOON_MIN_MASS_KG = 1e20` means **none of the real irregulars can exist today**: Himalia is about 7e18 kg [R], everything else less. The existing moon list cannot hold them, so Jupiter always comes out with a few regular and no irregular moons.
* `draw_inclination_deg("captured")` and `KOZAI_FORBIDDEN_DEG` are good and stay; the retrograde excess does not come out of them (inclinations are only half the story), it comes from the stability limit and collisional loss, which the code lacks.
* `planet_architecture` (Planets 1, GEN.728, merged or on its branch as the lane reported) already moves and scatters giants (hot Jupiters clear the inner chain). That **event is exactly the capture window**: a count of giant-giant close encounters is what NDD14 shows decides the yield (Jupiter needs a phase of encounters with an ice giant; "case 2 has the largest number of encounters and the most captures") [S]. The architecture must emit it.

## 4. Proposed model

### 4.1 Reservoirs and windows (what a planet can capture from)
A **reservoir** is a small-body population with an annulus, a mass at the window, a size law and a composition. Reservoirs come from the belt-evolution design: each belt or disc in the ledger, at its **birth mass** (before the early-depletion multiply) and again at its present mass.

A **window** is one of:
1. `gas`: only a body with a gas envelope (`envelope_mass_fraction` at or above `GIANT_ENVELOPE_FRACTION`), while the disc lasted. Stand-in strength 0.1 to 0.3 of the instability channel (JH07 says gas drag and pull-down cannot explain the ice giants; for gas giants they may add) [marginal educated guess].
2. `encounter`: the architecture's instability event, with `n_enc` giant-giant close encounters (one number per system, from the architecture). Efficiency `eps_enc` = 1e-8 to 4e-8 per reservoir particle that passes within the 3 AU encounter sphere (`f_3AU` of NDD14), scaled by the planet's Hill volume.
3. `exchange`: a binary planetesimal fly-by, only for bodies at least a few hundred km (Triton); probability tiny, one draw per giant against the binary fraction of the reservoir (the belt doc's L3 binary fraction), a massive moon at most.
4. `temporary`: no permanent body (a flag only, no object; a hook for a "mini-moon" description later).

### 4.2 The count: three tiers like the belt ledger
For each planet `p` and window `w`, the capture pool is an **expected number per size bin**:

    N0(>R) = sum over reservoirs  eps_w * f_geometry(p) * N_src(>R)

* `N_src(>R)`: from the belt/disc ledger (size law and top-end order statistic, `belt-evolution` sec. 3.1). For the outer disc use the 3 to 4e7 above R = 40 km, index 5/1.8 as a stand-in until the Kuiper-type disc of belt item 12 exists.
* `f_geometry(p) = (r_H(p) / r_H,Jup)^3`: the **stable volume** is the capture target, the Hill volume (stand-in, tier educated guess from physics, tweak L1). With this, an Earth-like planet gets 2e-5 of Jupiter's yield [C], which is why terrestrial planets have no irregular moons and Mars has none from this channel.
* After capture, **grind** with the belt law: `N(t) = N0 / (1 + t / t_c)` with `t_c` from the mean collision rate of the irregular zone (volume `~ r_H^3`, speeds from the orbit at `~0.4 r_H`), so heavy captures flatten onto one saturated population (test: the four Solar System giants must come out at the same order; see 4.5).
* Convert the expected numbers to **integers** with the draw module's order statistics: L1 tier is the top of the size law (at the lane's chosen cut), the remainder is L2 reproducible from the planet's seed and L3 a total, exactly as the belt ledger does. Poisson at the top end: the Jupiter-like case has mean about 1 body above 40 km, so some Jupiters have none and some three.

### 4.3 Orbit and survival
For each body: inclination from isotropic (`draw_inclination_deg`), a/r_H from a lognormal whose median is the planet's collision-set value (stand-in, interpolating 0.44 at 5 AU to 0.17 to 0.19 at 20 to 30 AU from the observed medians, [S] as data points), eccentricity 0 to 0.7 [S], q greater than the planet's outermost regular moon distance. Then **survival**: prograde bodies beyond `0.4 r_H` and retrograde beyond `0.55 r_H` (centre of the 0.5 to 0.6 range) are dropped, Kozai-band inclinations never arise. The retrograde to prograde ratio then falls out of the stability limits; check it against 88:19 as a calibration test rather than forcing it.

### 4.4 Composition and identity
* **Composition** follows the reservoir. A body captured from a belt keeps the belt's taxon share at its source distance (the ledger's gradient; belt item 8), from the outer disc the icy-body rule (`draw_icy_body`) as now ("captured": ice lost by U(0, 0.5), unchanged). A heterogeneous breakup (Himalia-style) draws sub-members with different volatile content from the parent [S, Sharkey 2025].
* **Identity.** Belt-evolution sec. 9 (decided by Boss): members print `<belt ID>.<rank hex>`, rank sorted on a mass drawn once at birth. For captured bodies:
  - A capture is a **flag and a parent link on the existing object**. A body that was an L1 member of a belt in the same system and is captured keeps its designation `<belt ID>.<rank hex>`; it gains `captured_by = <planet ID>` and `captured_epoch`. Its mass leaves the belt's ledger total (mass conservation). Which members are captured is a pure function of (belt seed, planet ID, window), never of insertion order.
  - A body from a reservoir that is not a belt of this system (the outer disc before it is a ledger belt, a passing population) takes `<planet ID>.<rank hex>`, the same sub-number pattern, 12-bit body field untouched, rank on a mass drawn once. Its source stays recorded as "outer disc".
  - Irregulars take **no names** (same rule as belt members); the moons with names are the regular, big moons and the few massive captures (Triton-class).
  - Planetary moons of the big kind keep today's issue (`store._UidIssuer`). The belt doc's warning that body numbers must be derived, not issued in insertion order, applies the same way: irregulars must not take the next insertion-order body number.

### 4.5 Counts to expect (generator targets; a Sun-like star, an instability occurred, gas giants then ice giants)
| Planet (present orbit) | `r_H` (million km) | stable limit pro / retro (million km) | Irregulars above R = 5 km | Above R = 40 km | Largest |
|---|---|---|---|---|---|
| Jupiter-like, 5.2 AU, 318 Earth masses | 53 | 26 / 36 to 29 | about 10 to 20 | about 1 | 50 to 150 km |
| Saturn-like, 9.6 AU, 95 | 66 | 32 / 36 to 40 | 5 to 15 | 0 to 2 | 50 to 150 km |
| Ice giant, 20 to 30 AU, 15 to 17 | 70 to 116 | 34 to 57 / 38 to 70 | 5 to 12 | 0 to 2 | 30 to 100 km |
| Warm Jupiter, 0.5 AU | 5.1 | 2.5 / 3 | 0 to 3 (the Hill sphere is small and stellar tide strips the outer zone) | about 0 | few km |
| Hot Jupiter, 0.05 AU | 0.5 | 0.25 / 0.3 (about 4 Jupiter radii) | 0 | 0 | none |
| Earth, Mars, super-Earth at 1 AU | 0.7 to 2.6 | 0.5 to 1.3 | 0 (expected 1e-4) | 0 | none by this channel |

The giant rows are the Solar System's observed head counts, so they are **calibration targets, not an independent prediction** [R]. The new and informative part is the other rows and the dependence on the instability (no instability: gas windows only, a fraction of these counts, and none for ice giants), on the Hill volume, and on the orbit.

An M-dwarf giant at 1.5 AU (23 million km Hill radius) gets the same scaling with its smaller disc mass: fewer.

Breakdown by origin for a Jupiter-like row, expected values [C from sec. 2.3]: instability capture dominant; belt-origin share = (mass of the belts inside the encounter region at the window / total reservoir mass), so a system with a heavy primordial belt next to the giant gives a larger asteroidal share, which is Boss's statement as a number.

### 4.6 What is not modelled (named so the table above is honest)
Temporary captures (flag only); capture at Trojan points (the belt doc's Trojan swarms are separate objects); regular-moon impacts that eject; captured binaries (a binary irregular is possible, the belt doc's L3 binary fraction can be given to it later); a planet that migrates after capture: **the Hill sphere shrinks with a smaller orbit and the moons beyond the new stable limit are lost** (`drop_unstable_moons` already does that for regular moons; irregulars need the same call with the new limits); a star that loses mass: `r_H` grows by `(M_formation / M_now)^(4/3)` for adiabatic widening [C], bodies that were just outside become bound candidates, a weird-world hook only.

## 5. Stand-ins (register in GEN.725 when built)
| Value | Stand-in | Tier | Tweak |
|---|---|---|---|
| Capture probability per reservoir particle | 1.3e-8 to 3.6e-8 (NDD14, Jupiter), scaled by Hill volume to other planets | Verified for Jupiter; educated guess from physics for scaling | L0 / L1 |
| Gas-window strength | 0.1 to 0.3 of the encounter channel | Marginal educated guess from physics | L3 |
| Retrograde stable limit | 0.5 to 0.6 `r_H`, 0.55 default | Verified range (VNL08) | L0 |
| Median `a/r_H` versus planet distance | 0.44 at 5 AU falling to 0.17 to 0.19 at 20 to 30 AU | from needed output (Table of medians) | L2 |
| Collisional timescale of the irregular zone | from Hill volume, orbit speed, grinding law | Marginal educated guess from physics | L3 |
| Disruption probability and family velocity spread | tens to hundreds of m/s | Unverified [R] | L3 |
| Binary-exchange probability for Triton-class | tiny; at most one per ice giant | Unverified | L4 until a source is read |
| Tier cut for irregulars (below) | see sec. 6 | from needed output | L1 |

## 6. Storage: where the irregulars go

Real irregulars are small: Himalia is about 7e18 kg, most are 1e12 to 1e17 kg [R]. The belt design's individual cut (1e19 kg) would store none, and the existing moon floor (1e20 kg) does not allow them. Options for the planet's irregular list (Jupiter-like, counts [C] from sec. 4.5 with the observed size law):

| Option | Cut | Stored per giant | Note |
|---|---|---|---|
| A | all above 1e19 kg (about R = 130 km) | 0 or 1 | too thin; the moons that are the point are missing |
| B | above 1e16 kg (about R = 11 km) | about 8 to 12 | the named-in-textbooks ones; cheap |
| C | above 1e14 kg (about R = 3 km) | 30 to 80 | the list a reader sees on a real page; paging needed on the system view |
| D | rule: from the top until 90% of the irregular mass is catalogued, cap 100 | 8 to 100 | consistent with the belt doc's rule E |

Recommended: **D, with floor B** (store every irregular above 1e16 kg; add the next heaviest until 90% of the population mass, up to 100). Smaller ones are L2 (drawn on demand from the planet seed, addressable by `<planet ID>.<rank hex>`) and L3 (a total mass and a size law on the planet), as the belt design has it. This is the belt doc's open question 3 again, applied to a second kind of reservoir; the decision for belts and irregulars should be the same rule so one function serves both.

The moon table gains a few columns (DB.70 style, `db/column_units.py`): `captured_by`, `capture_window` (gas, encounter, exchange), `source_reservoir` (belt ID or disc), `family_id`; new numeric columns need unit entries.

## 7. Items for TODO (numbers for TODO to assign)

Lock means it touches a lock (epoch, migration, golden, API version or a shared generator path).

**Non-lock (separate branch, one lane in sequence; none waits for a belt-ledger merge, they take simple inputs and later read the ledger):**
1. Read list and register: Nesvorny et al. 2003 and NVM07 (capture and stability), Hamilton and Burns 1991, Agnor and Hamilton 2006 and VNL08 (exchange, partly read), Bottke et al. 2010 (collisional history), Sheppard and Jewitt 2003 / Sheppard et al. 2005, 2006 (size laws, head counts), Granvik 2012 (temporary capture), Rosenblatt 2016 (Phobos, Deimos). Output: a short report with the head counts per giant read from tables and the register rows of section 5. Boss cannot access some papers; use second-hand citations if a source is closed.
2. `physics/irregular_zone.py`: `r_H`, stable limits (pro/retro), the inner edge from the outermost regular moon, Kozai exclusion (reuses the `KOZAI_FORBIDDEN_DEG`). Tests: Jupiter, Saturn, Uranus, Neptune values of section 2.1.
3. `physics/irregular_population.py`: windows, `N0(>R)`, grinding, order-statistic draw into L1/L2/L3 tiers. Tests: Jupiter-like gives about one body above 40 km, ten to twenty above 5 km; the four giants come out the same order; terrestrial planets get zero.
4. `physics/irregular_orbits.py`: a/r_H, e, i, survival, family drawing (parent breakup with a velocity spread). Tests: retrograde to prograde near 4.5:1 from stability alone, no body in the Kozai band, no q inside the outermost regular moon.
5. Collisional-saturation check (the guess of sec. 2.3): reuse `belt_mass_earth` for the irregular zone and show the invariance across Jupiter, Saturn, Uranus, Neptune, or replace the guess with a fitted scalar per planet and register it L3.
6. Composition from the reservoir (waits for or reads the GEN.473 taxon kernel: belt item 8) and the heterogeneous parent.
7. Designation for captured bodies: `<belt ID>.<rank hex>` kept with a parent link, and `<planet ID>.<rank hex>` for outside-reservoir bodies; a test that both are pure functions and round-trip through `parse_id` (shares belt item 6).
8. Make the captured-origin weight depend on moon mass (`ORIGIN_TABLE` row `captured`), so a 1e20 kg or larger moon is rarely captured. Counted as non-lock only if done as a pure helper; its wiring is lock (item 10).

**Lock:**
9. The architecture must record its capture windows (the instability event and the number of giant-giant close encounters, and the gas-disc lifetime already kept in `DiskSystem`): a field on the system, one epoch with the GEN.728 wiring (Planets 1).
10. Wire irregulars into `generate_moons` as a second tier (above the existing regular moons; `MOON_MIN_MASS_KG` is left alone for the regular list), and the captured-origin weight. New epoch, `--affects` plan, golden.
11. Store the L1 irregulars and the new moon columns (migration, API version), and the stability call for irregulars in `drop_unstable_moons` and in bring-to-epoch when a planet moves (GEN.1006). Replaces the open weights, no backward-compat shim.
12. Belt-ledger link: a captured member leaves the belt's mass total (belt items 9 to 11, the ledger storage); until those merge, captured bodies come from the outer-disc reservoir only. Epoch with the belt storage.
13. System map and pages: irregular orbits (eccentric, inclined, retrograde) on the moon-orbit view (`web/maps/moonorbits.py`), a body page for an irregular, and the clearance ledger (GEN.1014) taking the irregular zone as occupied space around the planet. Part of 8.2 object pages.

Order: 1; then 2 to 5 in parallel (2 and 3 first); 6 and 7 after the belt reads; 9 with the architecture train; 10 and 11 as one train; 12 and 13 after.

## 8. Decisions (Boss, 4:34 PM EDT 2026-10-30)
1. **Source of the asteroidal moons: decided, keep the windows** (gas era, instability, binary exchange); no continuous capture.
2. **Individual-cut rule: decided, the belt storage rule with a 1e16 kg floor** for irregular moons.

Items: GEN.1075 to GEN.1086 and TEST.495 (goal #4130).
