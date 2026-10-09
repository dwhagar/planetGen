# Planets and orbits in extreme environments: nebulae, pulsars, black holes and compact binaries

Boss's own research of 2026-10-09 (18:23Z, "Exotic Planetary and Stellar Dynamics"), digested against the design notes that already cover the same ground, with the numbers checked by calculation. It settles the open points of GEN.94 (can planets form in each nebula class) and GEN.130 (exotic star systems with a compact primary), and adds a few rules for GEN.95, GEN.129 and the later galaxy-scale items. Where Boss's text and an existing note agree, this note says so and adds nothing; where it adds, corrects or conflicts, it says which.

Informs: GEN.94, GEN.95, GEN.130, GEN.129, GEN.100, GEN.103, GEN.84 to GEN.89 (radiation inputs), GEN.9 and VIEW.2 (other galaxies)
Status: research, 2026-10-09; decisions marked "Boss" are his, everything else is a recommendation.

Evidence tags: [B] a claim from Boss's research text, with its cited paper not opened by the digest (the environment can read search snippets only); [C] computed here, scripts in `/mnt/project-files/research/scripts/exotic/`; [R] recalled and unconfirmed; [S] seen in an earlier note's search. The paper list Boss gave is reproduced under Sources.

Related notes: [nebula-and-asteroid-field-classes.md](nebula-and-asteroid-field-classes.md) (the GEN.94 rule table and GEN.95 application), [multistar-and-compact-systems.md](multistar-and-compact-systems.md) (stability criteria, the compact-primary rules and the GEN.130 go and no-go), [anomalies.md](anomalies.md) (pulsing and wind fractions), [galactic-potential.md](galactic-potential.md), [collisions-and-mergers.md](collisions-and-mergers.md).

## Summary

**GEN.94 is answered, and the existing rule table stands.** Boss's research reaches the same conclusion as the note: planet formation survives where a disc keeps its inner part, and almost every class is either too young or too irradiated for in-situ planets. Four things are new and worth adding to GEN.95:

1. **The photoevaporation cut should scale with the host's mass.** The existing table cuts outer discs at fixed radii (200, 50 and 10 AU for G0 above 1e2, 1e3 and 1e4). The physics Boss's text cites is the gravitational radius, `r_g = G M / c_s^2`, inside which the gas cannot escape: a 0.2 Msun star keeps a fifth of the disc a 1 Msun star keeps (section 1.1).
2. **Planetary-nebula second-generation planets need a binary.** The post-AGB discs seen are circumbinary, last only 1e4 to 1e5 yr, and only massive ones (about 0.1 Msun) form planetesimals; keep the rule "no in-situ planets in H to L" except a flagged rare case for a central binary (section 1.3).
3. **Engulfed giants become events, not losses.** A giant of more than about 5 Jupiter masses swallowed on the AGB is a luminous-red-nova candidate; the engulfment radius rule is unchanged (section 1.2).
4. **The formula in Boss's text for orbit expansion is garbled.** The correct adiabatic law is `a (M_star + M_p) = constant`; the existing note's `M_initial / M_final` factor is right. Boss's engine verified the law (Jeans 1924, Hadjidemetriou 1963, Veras et al. 2011) on 2026-10-09 18:45Z (section 1.2).

**GEN.130 stays a go in the three slices of [multistar-and-compact-systems.md](multistar-and-compact-systems.md) section 8, with these refinements:**

1. Pulsar planets: occurrence is about 0.1% of known pulsars (4 hosts in about 3,700 [R]) and about 0.7% of millisecond pulsars, which matches the note's "1% of millisecond pulsars" within the scatter; three origins are real and worth three presets (disrupted white dwarf fallback disc, captured circumbinary giant in a globular cluster, ablated white dwarf core), and the wind pressure makes every such planet a no-atmosphere, radiation-tier-3 world (section 2.1).
2. A wide tertiary around a black hole survives only if the black hole formed with almost no kick, because the outer orbit's escape speed is 2 km/s. V404 Cygni's 3,700 AU tertiary is the evidence (section 2.3).
3. Kozai-Lidov cannot be skipped for compact inner binaries, but it must be screened against relativistic precession, which quenches it in most triples, including V404 Cygni (section 2.3).
4. Black hole and neutron star "plunge versus disruption" is a lookup on mass ratio and spin, not a one-line tidal-radius test (section 2.3). It matters only for events, which stay out of scope.
5. Supermassive black hole planets ("blanets") and supermassive binaries stay no-go for star systems; the rule for the galactic nucleus phenomenon and for other galaxies is in section 2.4.

## 1. GEN.94 and GEN.95: planets in and around nebulae

### 1.1 Proplyds, the lifetime problem and the gravitational radius

Boss's text [B]: in the Orion Nebula Cluster, discs within 0.3 pc of the O6 star theta1 C Orionis show cometary ionisation fronts and lose `1e-7` to `1e-6 Msun/yr`; seven proplyds in NGC 1977 lie 0.04 to 0.27 pc from the B1V star 42 Orionis, so even a B star drives winds, at lower rates. A disc loses its mass in `t = M_disc / Mdot` [C]:

| Disc mass | Mass-loss rate | Lifetime |
|---|---|---|
| 0.04 Msun (HST 182-413, from the existing note) | 4.1e-7 Msun/yr | 9.8e4 yr |
| 0.01 Msun | 1e-7 | 1.0e5 yr |
| 0.01 Msun | 1e-6 | 1.0e4 yr |
| 0.1 Msun | 1e-7 | 1.0e6 yr |

This is the "proplyd lifetime problem" Boss's text names: core accretion needs several Myr, the disc is gone in 0.1 Myr. The existing note already gives the same lifetime. The resolutions Boss lists (very fast formation by gravitational instability or sequential inside-out growth; solids migrating inside the gravitational radius, where the gas is bound and shielded; a continuing supply of newly formed stars) are consistent with the disc surviving at small radius.

**New for GEN.95: scale the cut by host mass.** Wind launches from radii beyond about `0.2 r_g` [R, Liffman 2003 and later models], with `r_g = G M / c_s^2` and `c_s` the sound speed of the heated gas [C]:

| Host mass | Gas 100 K | 1,000 K (FUV heated) | 3,000 K | 10,000 K (EUV heated) |
|---|---|---|---|---|
| 0.2 Msun | 0.2 r_g = 56 AU | 5.6 AU | 1.9 AU | 0.3 AU |
| 0.5 Msun | 140 AU | 14 AU | 4.7 AU | 0.6 AU |
| 1.0 Msun | 280 AU | 28 AU | 9.3 AU | 1.3 AU |

Sound speed uses mu = 1.3 below 5,000 K and 0.6 above, so these are order-of-magnitude. The existing table's fixed cuts (200, 50 and 10 AU) match the 1 Msun column at roughly 140, 560 and 2,800 K, so they are a sensible calibration for a solar-mass star; the improvement is only to multiply the cut radius by `M / Msun`. A 0.3 Msun M dwarf at G0 = 1e3 then keeps about 15 AU instead of 50. A discs-survive test for the stars inside classes C, D, E and G is still the age rule (`PLANET_MIN_STAR_AGE_GY`) first and the G0 cut second.

ALMA survey [B]: 23 discs in the central ONC down to about 1.2 Jupiter masses of dust. This confirms solids remain in an irradiated cluster, which is why the rule keeps a reduced inner system rather than deleting planets.

### 1.2 The giant branches: expansion, engulfment, red novae

Boss's text [B]: stellar mass loss widens orbits while tides in the giant's envelope pull them in; engulfed planets feel drag `F = -0.5 C_d pi R_p^2 rho v^2`, spiralling in over 100 to 1,000 yr and depositing energy that can drive a "luminous red nova". Giants above about 5 Jupiter masses are named as the likely progenitors of low-luminosity red novae.

**Correction to the formula.** The text of 18:23Z gives `a_dot / a = -2 M_p_dot (M_star - M_p) / (M_star M_p)`, which depends on the planet's mass-loss rate and does not reduce to the known result. For mass lost from the star, slowly and isotropically, angular momentum conservation gives `a (M_star + M_p) = constant`, so

`a_dot / a = - M_star_dot / (M_star + M_p)`

(positive, an expansion, when the star loses mass). Integrated, a planet's orbit grows by `M_initial / M_final`.

**Verified by Boss's engine (2026-10-09 18:45Z).** The law `a (M_star + M_p) = constant` is the adiabatic result of Jeans 1924 (Monthly Notices of the Royal Astronomical Society 85, 2-11), formalised for variable-mass binaries by Hadjidemetriou 1963 (Icarus 2, 440-451) and treated for dying stars by Veras, Wyatt, Mustill, Bonsor and Eldridge 2011 (MNRAS 417, 2104-2123). Its conditions, which the generator should check before applying it:

- *Adiabatic*: the mass-loss timescale `M / Mdot` is much longer than the orbital period. Stellar winds on the giant branches and the planetary-nebula ejection of a white dwarf progenitor qualify.
- *Isotropic*: spherically symmetric loss. An asymmetric ejection kicks the star and moves the barycentre.
- *Eccentricity* stays constant in the adiabatic, isotropic limit.
- *Impulsive loss* (faster than an orbit, as in a supernova): position and velocity do not change, so the new orbit comes from vis-viva, `1/a_new = 2/r - v^2 / (G M_new)`, and the eccentricity depends on the orbital phase at the moment of loss. A previously circular orbit becomes unbound when `M_new < 0.5 M_old`. This is the Blaauw limit used in section 2.3 and in section 7.1 of the multi-star note.
- *Engulfment* overrides all of it: a planet inside the giant's envelope is dragged inward, not expanded.

Worked example (Boss's engine, corrected): a 1 Msun star becoming a 0.6 Msun white dwarf moves a planet from 1.0 AU to `1.0 / 0.6 = 1.67` AU (the engine printed 1.66). The generator uses the law for white dwarf hosts [C]:

| Progenitor to remnant | Outer orbit grows by |
|---|---|
| 1.0 to 0.54 Msun | 1.85 |
| 2.0 to 0.65 | 3.1 |
| 5.0 to 0.9 | 5.6 |
| 8.0 to 1.3 | 6.2 |

Keep the engulfment radius as the note has it (`1.5 + 0.8 (M_prog - 1)` AU, `WD_PROGENITOR_ENGULFMENT_AU = 1.5` for 1 to 2 Msun): planets inside it are lost; those just outside it are tightened by the envelope, not widened, so the survivor list near the radius needs an eccentricity excitation (a draw, e up to 0.3 [R]) rather than a pure expansion. New, optional: when the generator removes a planet of more than 5 Jupiter masses by engulfment, record a "swallowed giant" flag on the system; it gives the red-nova anomaly something to point at ([anomalies.md](anomalies.md)). The drag law is not needed in the generator.

### 1.3 Planetary nebulae and second-generation discs

The ejected envelope can form a circumbinary or post-AGB disc [B]. Observed discs look like young protoplanetary discs but last 1e4 to 1e5 yr, an order of magnitude shorter. Dust grows to millimetres; where the disc approaches 0.1 Msun the dust-to-gas ratio can trigger the streaming instability and pebble accretion, but at 0.01 Msun planetesimals do not reach isolation mass within the lifetime [B, Pourmand et al. 2025 and De Marco et al. 2025 in the list].

Rule for classes H to L in the rule table: keep "In situ: no", and add a flagged exception: a planetary nebula whose central star is a close binary (a post-common-envelope pair, section 5.3 of the multi-star note) may carry a young circumbinary disc, never planets. Planets in a planetary nebula are survivors beyond the engulfment radius, as the table says. There is no occurrence rate in Boss's text, so none is added; the rule is a flavour flag, not a population.

## 2. GEN.130: compact objects at the centre

### 2.1 Pulsar planets

Boss's text [B]: only a handful of pulsar planets are confirmed. The origins differ:

| System | Planets | Origin | Notes |
|---|---|---|---|
| PSR B1257+12 | 0.02, 4.3 and 3.9 Earth masses; 6.2 ms pulsar | Fallback disc of a tidally disrupted carbon-oxygen white dwarf | Second generation; pulsar speed above 326 km/s shows a large kick; planets would be carbon rich, with diamond interiors in the model of Margalit and Metzger 2016 |
| PSR B1620-26 b | 2.5 +- 1 Jupiter masses; orbits a pulsar and 0.34 Msun white dwarf at about 23 AU, almost 100 yr | Dynamical capture of a main-sequence star and planet in the globular cluster M4, about 12.7 Gyr old | Circumbinary; shows planets formed in a metal-poor early cluster |
| PSR J1719-1438 b | about 1 Jupiter mass, under 0.4 Jupiter radii | Stripped white dwarf core | "Chthonian" |
| PSR J2322-2650 b | about 1 Jupiter mass | Stripped white dwarf core; JWST sees a C2 and C3 rich atmosphere | Evaporating |

The existing note already lists B1257+12, B1620-26 b and J1719-1438 b with the same masses; J2322-2650 b is new. **Check on the stated density:** a Jupiter mass in 0.4 Jupiter radii has `1.33 g/cm3 / 0.4^3 = 21 g/cm3` [C], not the 11 g/cm3 in the text; the literature value is of that order [R]. Use the radius limit, not the density, when making the preset.

**Occurrence.** Four host systems among about 3,700 known pulsars [R, ATNF count] is 0.1%. Among millisecond pulsars (about 600 [R]) it is about 0.7% [C], consistent with the multi-star note's "about 1% of millisecond pulsars have planets, 0 of 151 young pulsars". Keep the 1% figure as the upper bound and the rule that young pulsars get none. Split by origin from four systems (low confidence): disrupted-white-dwarf disc 25%, captured circumbinary giant 25% (only where the system is in a globular cluster, so zero for the field unless clusters are modelled), ablated companion 50%.

**Wind pressure.** The pulsar wind's ram pressure on a planet is `P = xi * Edot / (4 pi r^2 c)` with `xi` of order 1. Computed with spin-down powers from the multi-star note [C]:

| Pulsar | Edot | Distance | Ram pressure | Times solar wind at 1 AU (2 nPa) | Earth-field standoff |
|---|---|---|---|---|---|
| PSR B1257+12 | 1.9e27 W | 0.19 AU | 6.2e-4 Pa | 3.1e5 | 0.9 planet radii |
| same | | 0.36 AU | 1.7e-4 | 8.7e4 | 1.1 |
| same | | 0.46 AU | 1.1e-4 | 5.3e4 | 1.2 |
| typical millisecond pulsar | 1.5e27 | 1 AU | 1.8e-5 | 8.9e3 | 1.7 |
| Crab | 4.5e31 | 1 AU | 0.53 | 2.7e8 | 0.3 |
| old 1 s pulsar | 4e24 | 1 AU | 4.7e-8 | 24 | 4.5 |

(Standoff is `R (B0^2 / 2 mu0 P)^(1/6)` for an Earth-strength dipole; it is the order of magnitude the model needs.) Close millisecond pulsar planets therefore have a magnetosphere at or below the planet's surface; they keep no atmosphere unless the wind is shielded or the wind is much weaker than the spin-down power suggests (a striped wind carries much of its power as Poynting flux [R]). Boss's text adds Alfven wings: the planet's motion through the magnetised wind drives currents that radiate low-frequency radio beamed along the wings [B]. That is a good phenomenon flavour for "pulsar planet" systems and needs no physics in the generator.

**Rules for the generator** (replacing the millisecond-pulsar row of section 7.6 of the multi-star note only where stated):

- A millisecond pulsar draws a planetary system with probability 0.7% (upper bound 1%).
- Type A (disc): 1 to 3 planets at 0.1 to 0.6 AU, masses log-uniform from 0.01 to 10 Earth masses, near-coplanar and near-resonant (B1257+12's outer pair is near a 3:2 resonance [R]); composition carbon rich ("carbon planet" flag); no atmosphere; radiation tier 3.
- Type B (captured giant): needs a cluster; mass 1 to 4 Jupiter masses at 20 to 50 AU; leave out until clusters exist.
- Type C (ablated core): a companion of 0.001 to 0.05 Msun (black widow, redback, chthonian) at an orbital period of hours to a day; stripped, carbon rich.
- Habitability: lowest tier always. This agrees with the existing no-go.

### 2.2 White dwarf and neutron star formation channels

Boss's text states the channels in prose: white dwarf mergers leave a carbon-dominated fallback disc whose high solid-to-gas ratio suppresses vertical shear instability and lets rocky planets form by gravitational instability before the disc dissipates [B]. This is the mechanism for type A above, and it also explains why young pulsars get none (no white dwarf to disrupt). No further generator rule.

### 2.3 Hierarchical compact systems

**Stability criterion.** Boss's text of 18:23Z wrote a generalised "Vynatheya" criterion; the formula as pasted does not match the published one and is not used. Boss's engine then identified the paper on 2026-10-09 18:45Z: Vynatheya, Hamers, Mardling and Bellinger 2022, "Algebraic and machine learning approach to hierarchical triple-star stability", Monthly Notices of the Royal Astronomical Society 516(3), 4146-4155 (arXiv 2207.03151, doi 10.1093/mnras/stac2540). The earlier recollection of the journal (Publications of the Astronomical Society of Australia) was wrong. Boss then supplied the paper (arXiv v2) at 18:51Z, and the equations below are read from it.

What is verified (by Boss's engine, from the paper's text) and what is not:

- The baseline is Mardling and Aarseth 2001 (MNRAS 321, 398-420), as in section 5.1 of the multi-star note: `R_p,out / a_in = 2.8 [(1 + m3/(m1+m2)) (1 + e_out) / sqrt(1 - e_out)]^(2/5) (1 - 0.3 i_mut / pi)`, with `R_p,out = a_out (1 - e_out)`.
- **Verified from the paper itself (2026-10-09 18:51Z).** Boss supplied the PDF (arXiv 2207.03151v2, preprint dated 7 September 2022); everything below is read from it, with the equations checked on the page image. The published MNRAS 516, 4146 text was not compared with this arXiv version. Definitions in the paper: `q_in = m2/m1 <= 1`, `q_out = m3/(m1+m2)`, `alpha = a_in/a_out`, `i_mut` in radians, `Y = a_out (1 - e_out) / [a_in (1 + e_in)]`; a triple is stable when `Y > Y_crit`.
  - **Equation 2 (Mardling and Aarseth 2001)**, quoted by the paper as `R_p,crit / a_in = 2.8 [(1 + q_out)(1 + e_out)/(1 - e_out)^(1/2)]^(2/5) (1 - 0.3 i_mut/pi)` with `R_p = a_out (1 - e_out)`. This matches the form in the multi-star note.
  - **Equation 3 (the effective inner eccentricity).** `e_in,max = sqrt(1 - (5/3) cos^2 i_mut)` (the quadrupole-order Kozai-Lidov maximum), `e_in,avg = 0.5 e_in,max^2` (the time-averaged separation is `a_in (1 + e_in,avg)`, Stein and Elsner 1977), and `e~_in = max(e_in, e_in,avg)`. The trick is Grishin et al. 2017's.
  - **Equation 4 (the updated criterion).** `Y~_crit = 2.4 [(1 + q_out) / ((1 + e~_in)(1 - e_out)^(1/2))]^(2/5) x [((1 - 0.2 e~_in + e_out)/8)(cos i_mut - 1) + 1]`, where `Y~` is `Y` with `e_in` replaced by `e~_in`. The system is stable when `Y~ = a_out (1 - e_out) / [a_in (1 + e~_in)] > Y~_crit`. The constant is 2.4, not 2.8, and the `(1 + e_out)` of Equation 2 does not appear; the `(1 + e~_in)` in the denominator and in `Y~` carries the inner eccentricity. The bracket in `cos i_mut` is monotonic; the near-polar "bowl" comes only from `e~_in`, which is largest near 90 degrees.
  - **Fit, not derivation.** Fitted to 10^6 MSTAR N-body triples, post-Newtonian terms off so the problem is scale-free. A triple counts as stable if it stays bound for 100 outer orbits and neither semi-major axis changes by more than 10%. Table 4 scores: Eggleton and Kiseleva 0.86, Mardling and Aarseth 0.90, Equation 4 0.93, the neural network 0.95. Equation 4 is weakest for retrograde orbits and high `e_in`. The paper does not model the argument of periapsis or the true anomaly, so the real boundary has a step structure the formula smooths over.
  - **Validity.** `0.01 <= q_in <= 1`, `0.01 <= q_out <= 100`, `1e-4 < alpha < 1`, `0 <= e_in, e_out < 1`, `0 <= i_mut <= pi`. Like Mardling and Aarseth it has no `q_in` term, which the paper says is fine except when both `q_in` and `q_out` are about 0.1 or lower. A planet has `q = m_p/m` near 1e-3, outside the fitted range, so for planets the formula is an extrapolation.
  - **Two gaps in the paper's text.** Equation 3 is real only for `cos^2 i <= 3/5` (39.2 to 140.8 degrees); the paper does not say what to do outside that window. The implementation in `vynatheya2022.py` takes `e_in,max = 0` there, so `e~_in = e_in` [C, implementation choice]. And the neural network's weights are on the authors' GitHub (Appendix A), which was not reachable; the six-layer, 50-neuron network cannot be rebuilt from the paper alone, so Equation 4 is the implementable form.
  - **Corrections to the engine's earlier descriptions.** Both passes were wrong about the form: the constant 2.8 is not replaced by a function of `e_in`, there is no piecewise polynomial `f(i)`, and the accuracies are 93% (Equation 4) and 95% (neural network) against 90% for Mardling and Aarseth, not 97% against 92%.
  - **Decided (Boss, 2026-10-09 18:56Z, "yes" to the coordinator's question).** Mardling and Aarseth for planets; Equation 4 for star-only triples inside its fitted mass-ratio range; Eggleton and Kiseleva as a cross-check. Put all three in one function so the choice is a switch. The code transcription is `/mnt/project-files/research/scripts/exotic/vynatheya2022.py`.
- Also unverified: a parallel 2022 boundary by Tory, Grishin and Mandel for low mass ratios, quoted as `a_in / R_p,out <= 10^(-0.6 + 0.04 q_out) q_out^(0.32 + 0.1 q_out)`. As pasted it gives `R_p,out / a_in >= 17` at `q_out = 0.01`, against 2.8 for the test-particle limit of Mardling and Aarseth [C], so the definition of `q_out` in that paper must differ from the one used here. Do not implement it from the pasted text.

**Unit tests.** The first two are Boss's engine's cases for Mardling and Aarseth, with the second corrected [C] (the engine multiplied the eccentricity term by 2.12 without the 2/5 power, giving 7.8 and 15.6; the whole bracket carries the exponent). The Equation 4 column is new, computed from the paper's equations [C, `vynatheya2022.py`]; "needed `a_out / a_in`" is `Y~_crit (1 + e~_in) / (1 - e_out)`.

| Case | Inputs | Mardling and Aarseth: critical `R_p,out / a_in` | Mardling and Aarseth: needed `a_out / a_in` | Equation 4: needed `a_out / a_in` |
|---|---|---|---|---|
| P-type, coplanar, circular | `m1 = m2 = 1`, `m3 -> 0`, `e = 0`, `i = 0` | 2.80 | 2.80 | 2.40 |
| S-type, coplanar, eccentric | `m1 = 1`, `m2 -> 0`, `m3 = 1`, `e_out = 0.5`, `i = 0` | `2.8 (2 x 2.121)^0.4 = 4.99` | 9.98 | 7.28 (`Y~_crit = 3.638`) |
| Equal masses, coplanar, circular | `q_out = 0.5`, all `e = 0`, `i = 0` | 3.29 | 3.29 | 2.82 |
| Polar | `q_out = 0.5`, all `e = 0`, `i = 90 deg` | 2.80 | 2.80 | 3.19 (`e~_in = 0.5`) |
| Retrograde | `q_out = 0.5`, all `e = 0`, `i = 180 deg` | 2.31 | 2.31 | 2.12 |
| Eccentric inner orbit | `q_out = 0.5`, `e_in = 0.6`, `i = 0` | 3.29 | 3.29 | 3.74 |

Equation 4 is looser than Mardling and Aarseth for prograde coplanar systems (2.40 against 2.80) and stricter for polar and eccentric-inner ones, which is the paper's stated gain. All of these belong in the unit tests of the criterion code (GEN.129). The V404 Cygni figures below use the same formulas.

**V404 Cygni** [B, Burdge et al. 2024]: a 9 Msun black hole in a 6.5-day orbit with a K star, and a tertiary at about 70,000 yr. Computed [C]:

| Quantity | Value |
|---|---|
| Inner semi-major axis (9.7 Msun, 6.5 d) | 0.145 AU |
| Outer semi-major axis (10.6 Msun, 70,000 yr) | 3,731 AU |
| Ratio | 25,700 |
| Mardling and Aarseth requirement | 2.9 (circular), 7.8 (e = 0.5), 59 (e = 0.9) |
| Equation 4 requirement (`q_out` = 0.093, `i` = 0) | 2.5 (circular), 5.7 (e = 0.5), 39 (e = 0.9) [C]; both mass ratios are about 0.1, the regime where the paper warns the formula is weakest |
| Outer circular and escape speed | 1.6 and 2.2 km/s |
| Kozai-Lidov timescale `(P_out^2 / P_in) (M/m3) (1-e^2)^1.5` | 3e11 to 3e12 yr |
| Inner GR apsidal precession period | 9,000 yr |

The triple is stable by a factor of 400 to 9,000 beyond Mardling and Aarseth (and about 10,000 beyond Equation 4 for a circular outer orbit). Kozai-Lidov cannot operate: its timescale is 20 to 200 times the age of the universe and GR precession is faster than it by more than seven orders of magnitude (in general, a KL cycle is quenched when the precession period is shorter than the cycle). The binary evidence, then, is only that the tertiary survived.

**Survival of a wide tertiary.** Mass loss alone leaves a circular outer orbit bound when the lost mass is under half the total before collapse (the Blaauw limit): for a 9 Msun black hole with companions of 0.7 and 0.9 Msun, a progenitor up to about 19 Msun keeps the tertiary and 20 Msun or more does not [C]. A kick has to be smaller than the outer escape speed, 2.2 km/s, on top of that. A hole formed by a normal supernova kick (tens to hundreds of km/s) loses the tertiary at once, so the survival is evidence for direct collapse with almost no mass or momentum loss [B], which matches the Monte Carlo of section 7.1 of the multi-star note (a 1,000 AU companion survives in 0.1% of cases even with a momentum-conserving kick). **Rule:** a bound companion beyond about 1,000 AU around a black hole or neutron star is allowed only for black holes flagged as direct collapse (progenitor above roughly 40 Msun at low metallicity [R]); otherwise the pre-collapse tree has to pass the kick and mass-loss check in 7.1, and fails for wide companions.

**Kozai-Lidov screen for compact inner binaries.** The multi-star generator already has a Kozai screen. Add the quenching terms: the screen rejects (or marks) an inner pair only if (1) the mutual inclination is between 39.2 and 140.8 degrees (test-particle limit, [B]), (2) the KL timescale is shorter than the system's age, and (3) the KL timescale is shorter than the GR precession period of the inner orbit. If all three hold, the pair's eccentricity cycles to near 1; whether it merges depends on the gravitational-wave time at the peak eccentricity (next paragraph). The existing note leaves mergers out of scope; the screen marks the system `kl_active` and leaves its handling to GEN.110 and the collisions note.

**Gravitational-wave merger time.** For a circular inner binary `t = (5/256) c^5 a^4 / (G^3 m1 m2 (m1 + m2))`, multiplied by `(1 - e^2)^3.5` for eccentric orbits (the approximate small-eccentricity-corrected form [R]). The orbit that merges in 13.8 Gyr [C]:

| Pair | Circular a | Eccentric a, e = 0.9 (about 300 times faster) |
|---|---|---|
| 1.4 + 1.4 Msun (two neutron stars) | 0.022 AU (4.7 Rsun) | about 0.09 AU [C, scaling a to the 1/4 power] |
| 1.4 + 8 Msun (neutron star and black hole) | 0.046 AU (9.9 Rsun) | about 0.19 AU |
| 10 + 10 Msun | 0.096 AU | about 0.4 AU |

So a compact binary closer than a few hundredths of an AU is a merger waiting to happen, and a generator that draws one older than its merger time has drawn an impossible system. **Rule:** when a compact-compact binary is drawn, reject it if the Peters time at its drawn separation and eccentricity is shorter than the system's age; this is a one-line test using the formula above.

**Neutron star and black hole plunge versus disruption.** Boss's text gives the physics [B]: the outcome depends on the mass ratio `q`, the black hole's spin through the innermost stable circular orbit, and the neutron star's equation of state (stiffer: disrupts; softer: plunges). A crude screen comparing the tidal radius `r_t = k R_NS q^(1/3)` with the ISCO radius [C] shows how sensitive it is (R_NS = 12 km, M_NS = 1.4 Msun):

| Black hole spin | r_ISCO | q at which they cross, k = 1 | k = 2 |
|---|---|---|---|
| 0 | 6.0 M | 0.95 (BH 1.3 Msun) | 2.7 (3.8 Msun) |
| 0.5 | 4.2 M | 1.6 (2.2 Msun) | 4.5 (6.4 Msun) |
| 0.9 | 2.3 M | 4.0 (5.5 Msun) | 11 (16 Msun) |

The answer moves by a factor of 3 with the unknown prefactor, so the screen is not usable as a rule. Numerical-relativity fits (Foucart 2012 and 2018 [R]) give the disrupted mass as a function of `q`, spin, compactness and the inclination of the spin; for a typical stellar black hole of 5 to 10 Msun and low spin the neutron star plunges and there is no counterpart. Mergers are events and stay out of scope; if a later item wants a "this pair merged with a kilonova" flavour line, take the fit from the paper then.

### 2.4 Supermassive black holes

**Blanets** [B, Wada, Tsukamoto and Kokubo 2019; Giang et al. 2021]: in the outer, snow-line-cooled part of a low-luminosity active galactic nucleus disc, ice-mantled grains coagulate (collisions between neighbours are slow because the disc moves as one fluid), growing to bodies of 20 to 3,000 Earth masses over the AGN's tens of Myr lifetime; Bondi accretion by the hole strips gas from the planet's Hill sphere, so they have almost no envelope. Radiative torque disruption spins grains up until they fracture and is a barrier except in dust-shielded zones.

Computed sanity checks for a 4e6 Msun hole [C]:

| Orbit | Speed | Period | Hill radius of a 1,000 Earth-mass body |
|---|---|---|---|
| 0.1 pc | 415 km/s | 1,480 yr | 13 AU |
| 1 pc | 131 km/s | 46,800 yr | 130 AU |
| 10 pc | 41 km/s | 1.5 Myr | 1,300 AU |

The objects would be dynamically ordinary (long periods, wide Hill spheres), so nothing in the orbital update fails for them. The issue is existence and observability: no blanet is confirmed. Boss's text names M51-ULS-1 b as an extragalactic microlensing candidate; the published candidate is an X-ray transit around an X-ray binary in M51 [R, Di Stefano et al. 2021]. Check the reference before citing it.

**Rule for the generator.** No stellar system in the galaxy has a supermassive hole as its primary, so blanets stay out of star systems and the sector fill. The Milky Way's Sgr A* is quiescent: there is no AGN disc and no blanets. A neighbouring-galaxy generator (GEN.9, VIEW.2, [multiple-galaxies.md](multiple-galaxies.md)) can flag the nucleus of an active neighbour as "blanet region, undetected" as text only.

**Supermassive binaries.** OJ 287 [B]: primary about 1.8e10 Msun, secondary about 1.5e8 Msun, e about 0.7, outbursts about every 12 yr as the secondary crosses the primary's disc twice an orbit; the timing drifts through relativistic apsidal precession, which constrains the primary's spin. The observed 12 yr is in the observer's frame; at redshift 0.306 [R] the rest-frame period is 9.2 yr, a semi-major axis of 11,500 AU (0.056 pc, 32 Schwarzschild radii of the primary) and a pericentre of about 10 Schwarzschild radii [C]. The 2023 NANOGrav 15-year data show a nanohertz gravitational-wave background consistent with a population of such binaries and with the final-parsec problem being solved in practice [B]. For the generator these matter only for the nucleus phenomenon: a galaxy produced by a merger (a neighbour, GEN.9) can have a binary nucleus with a separation in parsecs; the Milky Way has a single hole. No orbital-update work is needed.

## 3. Go and no-go, updated

| Item | Answer | Basis |
|---|---|---|
| GEN.94 rule table (A to W) | Stands. Add host-mass scaling to the G0 cut (1.1), the circumbinary exception for H to L (1.3), the swallowed-giant flag (1.2). | Boss's text agrees; the additions are small. |
| GEN.95 application | As in the nebula note, with the three additions above. | |
| GEN.130 slice (a) collapsed leaves with a kick check | Go. Add: a wide (over 1,000 AU) bound companion only for direct-collapse holes; reject compact-compact binaries whose Peters time is shorter than their age. | 2.3. |
| GEN.130 slice (b) presets | Go. Millisecond pulsars: 0.7% planet rate, three types (2.1). Black holes: unchanged. | 2.1. |
| GEN.130 slice (c) radiation fields | Go. Use the wind ram pressure for pulsar-planet atmosphere loss. | 2.1. |
| Young-pulsar planets | No-go. | Needs a disrupted white dwarf; none for a young pulsar. |
| Blanets in star systems | No-go. | 2.4. |
| Supermassive binaries | Text-only flag on the nucleus of merged neighbour galaxies. | 2.4. |
| Neutron star and black hole merger or kilonova events | Out of scope, as before. | 2.3. |

## 4. Corrections to the text and open verification

- The orbit-expansion law (section 1.2): verified by Boss's engine; the 18:23Z formula was wrong.
- The pasted stability criterion (section 2.3): the Vynatheya paper is real (MNRAS 2022, not PASA) and its Equation 4 is now read from the paper itself. The engine's descriptions of it (a function replacing the constant 2.8, a piecewise polynomial in the inclination, 97% accuracy) were wrong. The engine's S-type unit test should read 4.99 and 9.98, not 7.8 and 15.6.
- PSR J1719-1438 b's minimum density: 21 g/cm3 from the stated mass and radius limit, not 11 (section 2.1).
- M51-ULS-1 b is an X-ray transit candidate, not microlensing [R] (section 2.4).
- "Pulsar planets are about 1% of pulsars" in the nebula note is too high for normal pulsars: about 0.1% of all known pulsars and 0.7% of millisecond pulsars (section 2.1). The rule table's "about 1% of pulsars" for T and U should read "about 1% of millisecond pulsars".
- The proplyd numbers (mass-loss rates, distances, ALMA counts) and the V404 Cygni and pulsar planet parameters are [B] and agree with the earlier notes where the earlier notes cover them; none was opened at source.

## Evidence notes

[R] items to confirm when paper access is allowed: the `0.2 r_g` launch radius; the ATNF pulsar count (3,700) and the millisecond pulsar count (600); the J1719-1438 b density; the Tory, Grishin and Mandel boundary; the Foucart disruption fits; the direct-collapse progenitor threshold (40 Msun); the OJ 287 redshift (0.306); the M51-ULS-1 b detection method; the B1257+12 3:2 resonance and the approximate Peters-time eccentricity correction. All [C] values come from `calc.py` and `calc2.py` in `/mnt/project-files/research/scripts/exotic/`.

## Sources

Boss's research text of 2026-10-09 18:23Z cites the following (not opened by this digest): Agazie et al. 2023 (NANOGrav 15-year, ApJL); Amaro-Seoane et al. 2024 (SMBHs in hierarchical triples, MNRAS); Bhowmick et al. 2026 (arXiv 2606.12851); Burdge et al. 2024 (V404 Cygni triple, Nature); Chatterjee and Tan 2013 (inside-out planet formation); De Marco, Aleman and Akras 2025 (planetary nebulae, arXiv 2501.07869); Giang et al. 2021 (grain disruption barriers by radiative torques); Konacki and Wolszczan 2003 (PSR B1257+12 masses, ApJ 591 L147); Liu and Lai 2018 and Liu, Lai and Wang 2019 (NS and BH mergers in triples); Margalit and Metzger 2016 (white dwarf-neutron star merger and the pulsar planets, MNRAS); Pourmand et al. 2025 (second generation planet formation in post-AGB discs, arXiv 2509.03894); Ressler et al. 2025 (BH collisions with thin discs, OJ 287, arXiv 2509.18241); Toonen et al. 2016 (hierarchical triple evolution, arXiv 1612.06172); Vynatheya, Hamers, Mardling and Bellinger 2022 (algebraic and machine learning approach to hierarchical triple-star stability, MNRAS 516(3) 4146-4155, arXiv 2207.03151); Wada, Tsukamoto and Kokubo 2019 (planet formation around SMBHs); Winter et al. 2019 (a solution to the proplyd lifetime problem, MNRAS). Wikipedia pages on star formation, blanet, PSR B1620-26, pulsar planet and V404 Cygni are also cited.

Verified by Boss's engine on 2026-10-09 18:45Z: Jeans 1924, MNRAS 85, 2-11; Hadjidemetriou 1963, Icarus 2, 440-451, doi 10.1016/0019-1035(63)90074-6; Mardling and Aarseth 2001, MNRAS 321, 398-420, doi 10.1046/j.1365-8711.2001.03974.x; Veras et al. 2011, MNRAS 417, 2104-2123, doi 10.1111/j.1365-2966.2011.19393.x; Vynatheya et al. 2022 (doi 10.1093/mnras/stac2540; equations read from the arXiv v2 PDF Boss supplied on 2026-10-09 18:51Z). Also cited by that paper and used here: Grishin et al. 2017, MNRAS 466, 276 (the `e~_in` substitution); Stein and Elsner 1977 (time-averaged separation); Eggleton and Kiseleva 1995, ApJ 455, 640 (Equation 1 of the paper). The engine's reported citation counts are not used.

Formulas used: Peters 1964 gravitational-wave inspiral time; Mardling and Aarseth 2001 and Eggleton and Kiseleva 1995 (as in [multistar-and-compact-systems.md](multistar-and-compact-systems.md) 5.1); Kozai-Lidov timescale and test-particle inclination limits (standard, [R]); Blaauw 1961 mass-loss limit; Bardeen, Press and Teukolsky 1972 ISCO radius.
