# Group 1: sub-Neptune and super-Earth envelope and water-world occurrence

Scope: part 09 rows b4, e1, e2, e3, e4, e5 and the weight that links envelope loss to stellar mass, metallicity and age (called **e-link** below; part 09 has no explicit row for it, the dependence enters only through S and E_cum in e1/e2). Labels follow `_worker-rules.md`: READ (arXiv PDF text read by me), 2ND-HAND, DATA, COMPUTED, SNIPPET, RECALLED, NOT REACHED. Every arXiv id below was checked by downloading the PDF and reading its title and abstract. Part 09 section 2 table e and section 4 item R1 were read; the repo design doc `atmospheres-retention-and-classes.md` section 3.4 was read for the e1/e2 formula (the code in `src/` does not implement table e yet).

## 1. Summary and verdicts

The radius valley is measured well enough to tune an envelope-loss draw: the CKS sample (N = 970 planets, 703 hosts after cuts, [R3500]) gives the valley location, slope in period, slope with stellar mass and a null result for metallicity, each with a published 1 sigma, and the Kepler DR25 injection products and a fitted recovery curve are public, so an injection-recovery test is possible. The catch is that the two knobs named in part 09 (eps and f_env0) cannot be tuned alone: published fits ([R3797], N about 900 CKS planets) constrain them only jointly with a core-mass distribution (peak 3.7 Earth masses), and the data do not separate XUV photoevaporation (the mechanism of e1/e2) from core-powered mass loss ([R0354] favours core-powered at about 1 to 2 sigma). The water-world share (b4, e3, e4) is not constrained by any data set: published estimates run from under 3 percent water by mass for super-Earths [R3801] to 59 percent of mass-radius planets being water-rich [R0715], and every route that fixes the answer depends on the model. The envelope metal fraction Z (e5) has about eight JWST retrievals plus one 15-planet exo-Neptune retrieval; it is not a population distribution, and atmospheric Z is not bulk envelope Z.

| Row | What | Verdict | One-line reason | Status |
|---|---|---|---|---|
| e1/e2 | bare rock vs H/He envelope (eps, f_env0) | **a** with a proviso | CKS ensemble + public injections + published valley statistics; tune a three-knob set (eps or eta(v_esc), f_env0 distribution, core-mass distribution), not two | READ |
| e-link | stellar mass, [Fe/H], age weight | **b** | mass slope measured three ways (0.18 to 0.26); metallicity null (0.01 +0.05/-0.06); age only through one ratio (0.61 to 1.00) and a weak slope | READ |
| b4 | water-rich share | **c** | model-dependent spread 0 to 87 percent; no data-defined share; derive from formation physics, do not draw | READ |
| e3/e4 | water world vs hycean (eta) | **c** | eta has no population constraint; one contested planet (K2-18 b); hycean is a temperature and H2 gate, so derive by physics | READ |
| e5 | envelope Z, bins 0.01 to 0.6 | **b** | N about 8 planets plus 15 exo-Neptunes; upper edge 0.6 is below several retrievals; atmosphere Z differs from bulk Z | READ |

Where a verified source overrides a partly verified one (rule of 3:02 AM EDT): (1) the fixed eps range 0.05 to 0.2 (part 09, 2ND-HAND) is overridden by the read variable-efficiency law eta = 0.1 (v_esc/15 km/s)^-2 (Owen and Wu 2017, READ in part 04), which gives 0.35 at 8 km/s and 0.036 at 25 km/s; (2) "K2-18b near solar to tens of x" (part 09) is overridden by [R0340] (READ): water 0.6 to 162 times solar; (3) the "20 to 35 percent of M-dwarf bare planets water-rich" quoted second-hand through [R0978] is replaced by the read [R0714] text (about 10 to 30 percent G dwarfs, 20 to 35 percent M dwarfs, model output, not data); (4) [R1478]'s "flat-bottomed, band of forbidden sizes" reading is overridden by [R1477] (READ): the gap is partially filled and the planets in it are real.

No-fixed-content flags: none of the e1/e2, b4, e3/e4 knobs rests on the Solar System alone (it has no sub-Neptune). Two places lean on it: the Z-versus-mass trend that [R2802] extrapolates from Solar System giants (do not use it to set e5), and the Earth-like core density (32.5 percent iron) used as the rocky end member; the latter is also supported by exoplanet data ([R3797], rho_tilde = 0.26 +0.08/-0.09 against 0.33 for Earth). E_cum = 3.1e16 J/m2 comes from the project XUV model, which is not read here (see section 5).

## 2. Detail

### 2A. e1/e2 (bare rock vs retained H/He) and e-link (stellar mass, metallicity, age)

These share one data set (CKS and Kepler DR25), one fit (valley location, slope, depth) and one test, so one builder should do them together.

**1. Ensemble.**

| Data set | N | Measured | Accuracy and selection | Public | Source | Status |
|---|---|---|---|---|---|---|
| CKS DR1 [R1478] | 2025 planets; 900 after cuts; N_star = 36,075 for completeness | Rp, occurrence per radius bin | median Rp error 12 percent; Kp<14.2, b<0.7, P<100 d, 4700<Teff<6500 K; K-S p = 0.003 vs log-uniform, dip test p = 1.4e-3; peaks 1.3 and 2.4 Earth radii, gap 1.5 to 2.0 with deficit factor >= 2 | table 2 machine-readable at the journal; no licence text found | [R1478] | READ |
| CKS-Gaia [R1477] | 907 planets, 1189 stars | Rp to 5 percent; gap vs P, S, M* | IDEM occurrence with 25 percent completeness cut; ABC cross-check agrees within 1 sigma; gap partially filled; super-Earth and sub-Neptune frequencies nearly equal | code `github.com/California-Planet-Search/cksgaia` (README reachable, data on Dropbox links); table 4 machine-readable | [R1477] | READ |
| CKS DR2 [R3500] | 970 planets, 703 hosts, M* 0.5 to 1.4 | gap vs P, S, M*, [Fe/H], age | Rp error <= 20 percent; b<0.8; T/Tmax>0.6; detection model from [R0786] injections | VizieR J/AJ/163/179 (tables 1, 2; front page reachable, table download not tested); software `github.com/petigura/CKS-Cool` named in the paper, raw README returned 404 | [R3500] | READ |
| Asteroseismic sample [R4504] | 117 planets | gap slope | Rp median error 3.3 percent | not checked | [R4504] | READ |
| Short-cadence refit (Ho and [R1963]) | 431 planets | gap vs P, M*, age; depth | refit transits, deeper valley than earlier | not checked | Ho and [R1963] | READ |
| Gaia-Kepler ([R0353] [R0353, R0354] [R0354]) | sample size not extracted | gap vs M*, age, S | 3-D gap in (Rp, S, M*) | not checked | Berger | READ (abstracts) |
| Kepler DR25 injections | 146,294 on-target injections (INJ1) | pipeline recovery | five injection runs; detection efficiency and robovetter results as tables; last updated 19 Dec 2017 | NASA Exoplanet Archive "Simulated Data" page, reachable; licence text not on the page | NASA Exoplanet Archive | DATA (page read) |
| [R0786] | 146,295 injected targets | completeness of DR25 | 90 to 95 percent recovery for strong signals | archive page above | [R0786] | READ |
| Low-mass stars: [R0828] | 328 planets (275 Kepler + 53 K2), 126 for M* 0.08 to 0.65 | valley around Teff<4700 K | valley centre 1.54 +-0.16; flux slope -0.060 +-0.025; rocky:non-rocky 0.5 +-0.1 (mid-K) to 8.5 +-4.6 (mid-M); valley not visible in the 126-planet low-mass subset (counts and that remark are 2ND-HAND via [R3500]) | not checked | [R0828] | READ (abstract) / 2ND-HAND (counts) |
| Ho 2024 [R1962] | 72 refit planets + the 431-planet FGK set; bins N = 63, 141, 148, 92 | valley depth E vs M* | E_avg = 1.25 (+0.55/-0.17) for M*<0.7; 2.72 (+1.10/-0.74), 2.87 (+1.13/-0.76), 4.62 (+2.80/-1.62) for the higher bins; linear slope 5.36 (+3.19/-3.23) per solar mass | not checked | Ho 2024 | READ |
| TESS mid-to-late M ([R1568] preprint, [R1568]) | 8134 stars, 77 candidates | occurrence | 0.954 +-0.147 super-Earths and 0.148 +-0.045 sub-Neptunes per star, ratio 5.5:1; valley absent (unimodal, peak 1.25 +-0.05) | not checked; preprint, review status unknown | [R1568] | READ |

Occurrence per star of sub-Neptunes from [R3500] (READ, P = 1 to 300 d): 1.37 (+0.19/-0.18), 0.79 (+-0.08), 0.41 (+-0.04) for M* = 0.5 to 0.7, 0.7 to 1.0, 1.0 to 1.4. [R3797] give 0.72 +-0.03 planets per star for P = 1 to 100 d, 0.95 to 6 Earth radii (READ). [R3252] ([R3252], N = 1130 planets) find a mixture preferred by the data with an intrinsically rocky fraction Q = 0.20 (+0.06/-0.05) and 1.05 +-0.05 planets per star, and warn the occurrence is model dependent (READ).

**2. Physical bracket.**

| Quantity | Value | Source | Status |
|---|---|---|---|
| energy-limited loss | Mdot = eta pi R^3 F_XUV / (G M K_tide); eta 0.05 to 0.2 | [R2740] via part 04 | 2ND-HAND |
| variable efficiency | eta = 0.1 (v_esc/15 km/s)^-2; valley radius scales as eta^0.18 X_HE^0.19 rho^-0.24 | Owen and Wu 2017 via part 04 | READ (part 04) |
| initial envelope | "more than a few percent" of mass, cores near 3 Earth masses, Earth-like core density; formation inside the ice line | Owen and Wu 2017 abstract | READ |
| fitted birth population, model I (log-Gaussians) | core mass mu_M = 3.72 (+0.45/-0.33) Earth masses, sigma_M = 0.44 (+0.06/-0.03); initial envelope mu_X = 0.040 (+0.015/-0.016), sigma_X = 0.51 (+0.20/-0.12) (log units of its Eq. 1) | [R3797: Table 1] | READ |
| fitted birth population, model II (free shape) | envelope distribution "strongly peaked at about 1e-2"; abstract quotes a typical 4 percent; core peak about 4 Earth masses; core-accretion theory over-predicts H/He by about 5 | [R3797] | READ |
| photoevaporation efficiency index | alpha_eta < 1.36 at 2 sigma (data allow a variable eta) | [R3797] model III | READ |
| core-powered alternative | valley slope dlogR/dlogP = 1/(3(1-beta)) = -0.11; ice fraction of cores < about 20 percent | [R1734] | READ |
| core-powered stellar dependence | dlogR/dlogM* about 0.33-0.35; location independent of age and [Fe/H] to first order; sub-Neptune radius slopes about +0.1 with [Fe/H] and -0.1 with age | [R1735] | READ |
| repo formula, COMPUTED check | f_lost = eps pi R^3 E_cum S / (G M^2), R = M^0.27 Earth radii; at S = 1 and M = 5, f_lost = 1.56e-4. Insolation where f_lost = f_env0: M = 4: S = 49 (f 1 percent), 148 (3 percent), 492 (10 percent) at eps = 0.1; the eps bracket 0.05 to 0.2 moves each by x2 to /2 | COMPUTED from the repo doc inputs | COMPUTED |
| observed transition | planets at 1.3 Earth radii sit at S > 200, planets at 2.4 at S < 80 | [R1478] | READ |

The bracket f_env0 = 0.005 to 0.10 and eps 0.05 to 0.2 therefore contains the fitted values (4 percent, 1 percent, eta about 0.1) and spans S_crit from about 25 to 1000 for a 4 Earth-mass core, which covers the observed transition. The proposed draw is consistent with the data but wide enough to be uninformative until tuned.

**3. Verdict.** (a), conditional. A published ensemble and bracket exist. The proviso: the two part-09 knobs are degenerate with the core-mass distribution and with the mass-loss mechanism, so the builder must tune a joint set, and the repo's bare/retained binary has no intrinsic spread, so the gap would be empty where the data show it partially filled [R1477]. A radius model for retained envelopes (not in the repo) and a width knob are needed. The quantity should be derived by physics (energy-limited or Bondi-limited loss integrated over the star's XUV history) and only the birth population (core mass, f_env0) drawn.

**4. Knobs and targets.**

| Knob | Tuned against | Target and published uncertainty | Valid range | Source | Status |
|---|---|---|---|---|---|
| eps (or eta(v_esc) index) | valley slope in period | m = dlogRp/dlogP = -0.11 +-0.02 (Petigura); -0.09 (+0.02/-0.04) (Van Eylen); -0.096 (+0.023/-0.027) (Ho and Van Eylen) | P 1 to 100 d, M* 0.5 to 1.4 | as listed | READ |
| eps and f_env0 jointly | valley location | Rp,0 = 1.84 +-0.03 Earth radii at P = 10 d (Petigura) | same | [R3500] | READ |
| f_env0 (width and centre) | super-Earth to sub-Neptune ratio | about 1:1 over P < 100 d (Fulton and Petigura); 0.61 +-0.09 (age < 1 Gyr) and 1.00 +-0.10 (> 1 Gyr) [R0353] | FGK | as listed | READ |
| core-mass distribution (new knob) | peak and width | mu_M = 3.72 (+0.45/-0.33) Earth masses, sigma_M = 0.44 (+0.06/-0.03) | CKS FGK | [R3797] | READ |
| mechanism choice (XUV vs core-powered) | 3-D gap slopes | alpha = dlogRgap/dlogS = 0.069 (+0.019/-0.023), beta = dlogRgap/dlogM* at fixed S = -0.046 (+0.125/-0.117) [R0354]; predicted core-powered alpha 0.08, beta 0.00; photoevaporation alpha 0.12, beta -0.17 (as quoted by Berger from [R3796]) | needs about 5000 planets with Rp error <= 5 percent to separate [R3796] | [R0354, R3796] | READ |
| stellar-mass weight | valley vs M* | dlogRp/dlogM* = 0.18 (+0.08/-0.07) (Petigura); 0.231 (+0.053/-0.064) (Ho and Van Eylen); 0.26 (+0.21/-0.16) [R0353]; sub-Neptune size exponent alpha = 0.25 +-0.03, super-Earth 0.02 +-0.03 (Petigura) | M* 0.5 to 1.4 | as listed | READ |
| stellar-mass weight (depth) | depth E_avg | 1.25 (+0.55/-0.17), 2.72 (+1.10/-0.74), 2.87 (+1.13/-0.76), 4.62 (+2.80/-1.62) by M* bin | M* below 0.7 to above 1.1 | Ho 2024 | READ |
| metallicity weight | valley and sizes vs [Fe/H] | dlogRp/d[Fe/H] = 0.01 (+0.05/-0.06); sub-Neptune exponent -0.01 +-0.02, super-Earth -0.00 +-0.02 | [Fe/H] -0.5 to +0.5 | [R3500] | READ |
| age weight | valley vs age | dlogRp,valley/dlog(age) = 0.033 (+0.017/-0.025) (Ho and Van Eylen); sub-Neptune size gamma = -0.01 +-0.03 (Petigura, Teff > 5500 K, ages 1 to 10 Gyr) | 1 to 10 Gyr | as listed | READ |

Contradictions inside the verified set (no override, report the spread): Petigura finds no sub-Neptune size dependence on [Fe/H] (+-0.02) and none on age, while [R1735] predict +0.1 and -0.1; the data disfavour the core-powered metallicity and age slopes at about 5 sigma for metallicity and 3 sigma for age if the quoted errors are taken at face value (COMPUTED as 0.1/0.02 and 0.1/0.03). Petigura also notes the age uncertainty is 0.1 to 0.5 dex, so the age null is weak. The M-dwarf valley is disputed: Cloutier and Menou report a positive slope in P (0.058 +-0.022) while Petigura measures -0.12 (+0.06/-0.04) for M* = 0.5 to 0.7; Ho 2024 finds a shallower valley; [R1568] finds none.

**5. Consistency test.** Simulate planets for the generator's own host stars (M*, [Fe/H], age), apply e1/e2 and a retained-envelope radius law, then apply a selection model: geometric transit probability p_tr = 0.8 R*/a [R3500] times the recovery curve p_det(MES) = c * Gamma-CDF(MES; a, b) with (a, b, c) = (26.1, 0.320, 0.941), zero below MES = 10 (READ; COMPUTED values 0.41 at MES 8, 0.80 at 10, 0.92 at 12, 0.94 above 15). Then compare the simulated valley (KDE minima train, 1000 bootstraps as in Petigura) to the numbers above with the published 1 sigma as tolerance, in M* bins 0.5 to 0.7, 0.7 to 1.0, 1.0 to 1.4. Because the 3-D gap test is uninformative below about 5000 planets [R3796, R3797], the test cannot choose between mechanisms; it can only check that the tuned set reproduces m, Rp,0, the ratio and the depth. Systematic offsets between gap-finding methods are acknowledged by [R0354], so use one method on both sides. The CKS-Cool code could not be fetched (raw README 404), so the fit routine must be rewritten from the paper (NOT REACHED for the code).

**6. No-fixed-content flags.** None. E_cum (project XUV model, not read here) is calibrated on a Sun-like history; replace it by the generator's own star XUV track when e-link is wired.

**7. Group for a builder.** Builder A: e1, e2, e-link with the CKS fit, the injection-recovery step and the valley test. The same selection model and the same Petigura/Ho tables also serve table-d shoreline checks only loosely (not this group).

### 2B. b4, e3, e4 (water-rich share, water world, hycean)

**1. Ensemble.**

| Data set or model result | N | Result | Source | Status |
|---|---|---|---|---|
| [R3797] forward model on CKS | about 900 planets | cores "ice-poor, Earth-like": rho_tilde = 0.26 (+0.08/-0.09), sigma_rho_tilde < 0.16 (2 sigma), Earth 0.33; assumes H/He on every sub-Neptune and photoevaporation; two-species mass-radius tables [R1429], a 3-species check says ice content stays negligible | [R3797] | READ |
| [R3795, R3801] | super-Earths below the valley | population water mass fraction <= 3 percent at 1 sigma; end-members tested | [R3801] | READ |
| [R1734] | CKS valley | cores with water-ice fraction < about 20 percent | [R1734] | READ |
| [R4880] | growth model + Monte Carlo, no N | planets of 2 to 4 Earth radii are "water worlds" with > 1/4, possibly > 1/2 water by mass; no occurrence fraction stated | [R4880] | READ |
| [R0600, R0601] | population synthesis | sub-Neptunes water-rich; planets with some H/He are 13 to 16 percent; needs photoevaporation to fill the super-Earth peak | [R0600, R0601] | READ |
| [R4528], 2024 [R4529] | population synthesis | rocky cores near 3 (max 5) Earth masses, icy near 10; depletion of water mass fractions 0 to 20 percent; steam above 400 K gives +30 percent radius at 1 Earth mass; valley fades for M* 0.1 to 0.4 | Venturini | READ |
| [R0713, R0714] | Genesis library + loss models | about 10 to 30 percent (G) and 20 to 35 percent (M) of bare planets water-rich; most outside P = 10 d; bare water-rich and rocky-with-thin-air are degenerate | [R0714] | READ |
| [R0715] | Kepler P-Rp mixture (rocky + water worlds, no H/He) and TEPCat | water mass fraction peaks near 41 percent; TEPCat sample classified 26.2 percent rocky, 59 percent water-rich, 14.8 percent H/He-rich; at least 20 percent of sub-Neptunes need H/He; TEPCat N not extracted; a preprint | [R0715] | READ |
| [R2774] | 34 M-dwarf planets (of 43 listed) | three populations: rocky 0.94 +-0.13, water 0.47 +-0.05, puffy 0.24 +-0.04 of scaled Earth density; separation at 0.65 | LP22 | READ |
| [R3430] | 46 M, 61 K, 63 G + 9 late F planets | no distinct M-dwarf water-world population; no gap in density or radius; minimum sub-Neptune mass 1.9, 3.4, 4.3 Earth masses for M, K, FG | [R3430] | READ |
| [R0978] | 60 M-dwarf planets | two populations preferred: GMM AIC/BIC 9.37/19.76 (two) vs 14.05/30.67 (three) | [R0978] | READ |
| [R3798] | M-dwarf sample of LP22 | H/He sub-Neptunes with boil-off are degenerate with 1:1 water worlds below 10 Earth masses; mass and radius alone insufficient | [R3798] | READ |
| [R3443] (SPRIGHT) | not fetched | reported (via Dainese) to disfavour a water-world population | via [R0978] | 2ND-HAND |

**2. Physical bracket.**

| Constraint | Value | Source | Status |
|---|---|---|---|
| hycean definition | H2 <= 0.1 wt percent, water >= 10 wt percent, ocean bottom 300 to 400 K, radii up to 2.6 Earth radii at 10 Earth masses | [R2804] via part 05 | READ (part 05) |
| liquid vs steam | water is fully steam for Teq above 400 K | [R4529] abstract | READ |
| solids cap | refractories > 50 wt percent of solids caps ice plus water at 50 wt percent | [R4335] via part 09 | READ (part 09) |
| ice-line material | about 50 percent water by mass beyond the snow line | [R0300] review [R0300]; Chakrabarty | READ |
| TRAPPIST-1 densities | all seven fit one rocky relation with 21 wt percent Fe (vs 32 Earth) or water-layer alternative | [R0039] | READ |
| rocky/non-rocky threshold | R_thresh median 1.48, 95 percent upper 1.59 Earth radii (N = 27 RV planets) | [R3800] | READ |

**3. Verdict.** (c) for b4 and e3/e4. The spread of published water-rich shares is 0 to 87 percent depending on the model (< 3 percent mass, [R3801]; ice-poor cores, Rogers and Owen; 10 to 35 percent, [R0714]; 59 percent, [R0715]; 84 to 87 percent, [R0600, R0601]), and the two ends of the spread are both READ, so neither overrides the other. The part-09 value c[0.05, 0.4] s^3 sits inside the spread but cannot be narrowed from public data. The share is the output of a formation model (migration of ice-rich embryos), so it should be derived from the snow-line law and a migration rule, not drawn. The eta split in e3/e4 is not a probability: hycean is a gate on Teq (< about 400 K, steam otherwise), on H2 mass fraction (<= 0.1 percent) and on ocean mass (>= 10 percent), so e4 = p[b4] times a physical gate and e3 the remainder. The only hycean candidate, K2-18 b, is contested: [R2805] read it as hycean-compatible with about 1 percent CH4 and CO2, [R4731] ([R4731], READ) argues a lifeless hycean world makes < 1 ppm CH4 and a 100x solar mini-Neptune fits, and [R3299] ([R3299], READ) shows magma-ocean chemistry can raise sub-Neptune H2O, CH4, CO2 without icy accretion.

**What breaks the degeneracy, with N.** Mass and radius alone do not [R0714, R3798]. The discriminators that exist: (i) population inference on the radius distribution, which excludes ice-rich cores only if every sub-Neptune holds H/He (N about 900, Rogers and Owen); (ii) the dependence of the valley on period and stellar mass, where water-world models predict a flatter or filled valley at low M* ([R4529]; Ho 2024, N = 72 + 431; Cloutier and Menou, N = 328) and photoevaporation predicts a steeper one, but both families can be tuned; (iii) the density gap of M-dwarf planets, where three independent analyses disagree (LP22 N = 34 vs Parc N = 46 vs Dainese N = 60); (iv) transmission spectra giving mean molecular weight: GJ 9827 d is steam-rich (H2O volume mixing ratio > 31 percent, no escaping He, [R3517]), TOI-270 d has a metal-rich miscible envelope [R0339], GJ 3470 b (11.2 Earth masses) needs a deep H2 envelope, and He escape is detected in K2-18 b, GJ 3470 b, GJ 436 b and tentatively GJ 1214 b (caption of [R3798], READ). That is about 6 planets with an informative spectrum and the interpretation of each depends on the interior model. Monte Carlo power estimates (Sun 2024) say about 45 JWST planets are needed to see even a mass-abundance trend (Sun 2024, [R4258]); no equivalent number was found for the water-world fraction.

**4. Knobs.** Proposed: s^3 weight with amplitude c (b4) kept as a bracket 0 to 0.87 until a formation model is chosen; if one is chosen, tune the migration-efficiency knob against the density-valley depth and the sub-Neptune minimum mass vs M* (Parc: 1.9, 3.4, 4.3 Earth masses for M, K, FG dwarfs, READ, observational-bias caveat in the paper). Targets for a consistency test with their uncertainties: radius-valley depth by M* (Ho 2024 table above), super-Earth:sub-Neptune ratio 5.5:1 for mid-to-late M ([R1568], preprint), 0.61 to 1.00 with age [R0353]. For e4: gate Teq < 400 K, H2 fraction <= 0.001, water >= 0.1, no draw.

**5. Consistency test.** Mass-radius planets need a mass-radius sample, not the Kepler radius sample: compare the simulated density distribution (after TEPCat or PlanetS selection: mass error <= 25 percent, radius error <= 8 percent, as used by LP22, Parc and Dainese) with the observed one using the same classifier. No published tolerance exists because the three classifiers disagree; use the spread. A selection model is needed (RV and TTV follow-up bias towards bright hosts and large masses, not modelled in any of the cited papers).

**6. No-fixed-content flags.** None for b4/e3/e4. The 1:1 rock to water ratio (LP22, Zeng) is taken from Solar System ice-to-rock composition by LP22 (they cite the Solar System, 2ND-HAND via the paper text "if it has the same composition as the Solar System"); the generator should use its own snow-line and ice-fraction law, not that number.

**7. Group for a builder.** Builder B: b4, e3, e4 with the formation-model choice (these need the owner's decision first). Parc's PlanetS catalog and TEPCat are the shared data set.

### 2C. e5 (envelope metal mass fraction Z)

**1. Ensemble (JWST and related; per-planet, no population fit).**

| Planet | Mass, radius, Teq | Result | Source | Status |
|---|---|---|---|---|
| TOI-270 d | 4.78 Earth masses, 2.2 Earth radii, 350 to 380 K | mean molecular weight 5.47 (+1.25/-1.14), envelope metal mass fraction 58 (+8/-12) percent, C/O 0.47 (+0.16/-0.19), interior 90 (+3/-4) wt percent rock/iron | [R0339] | READ |
| K2-18 b | 8.63 +-1.35 Earth masses, 2.6 Earth radii | water 0.6 to 162 times solar (HST); JWST CH4 and CO2 near 1 percent each; 100x solar mini-Neptune fits | [R0340, R2805, R4731] | READ |
| GJ 9827 d | 1.98 Earth radii, 620 K | O/H about 4 by mass, H2O volume mixing ratio > 31 percent | [R3517] | READ |
| GJ 3470 b | 11.2 Earth masses, 600 K | metallicity about 100x solar from H2O, CO2, SO2; CH4 below equilibrium | [R0304] | READ |
| GJ 1214 b | n/a | metallicity above 100x solar (models 1 to 3000x) with thick haze | [R2316] | READ |
| GJ 3090 b | 2.13 Earth radii, about 700 K | > 100x solar at 3 sigma (best fit about 2300x); text cites 100 to 500x for five earlier planets and 1000x for GJ 1214 b | [R0051] | READ |
| V1298 Tau b | 9.85 Earth radii, 12 to 15 Earth masses (from scale height), 10 to 30 Myr | log Z = 0.6 (+0.4/-0.6) times solar; low by an order of magnitude vs mature sub-Neptunes | [R0255] | READ |
| 15 exo-Neptunes, 200 to 1000 K | HST transmission spectra | best described by low sedimentation efficiency (f_sed about 0.1) and 100x solar, intrinsic scatter about 0.5 scale height | [R0516] | READ |
| 19 planets, H2O, Na, K | HST | H2O abundance rises as mass falls; revisit with host-star abundances cannot tell a linear trend from independent draws; needs about 58 planets (HST) or 45 (JWST) for Bayes factor > 10 | [R4663]; Sun 2024 [R4258] | READ |
| carbon abundances | 2 sub-Neptunes (K2-18 b, TOI-270 d) | "too early to make a robust assessment" of a trend; possible carbon plateau | [R2802] review, [R2802] | READ |
| TOI-270 d interpretation | model | magma-atmosphere interaction fits the spectrum without icy accretion | [R3299] | READ |

[R4590] (HAT-P-26 b) and [R0932] ([R0932], READ abstract only: spectral-feature amplitude correlates with Teq or H/He fraction at 95 percent, "complex degeneracies" remain) add no population Z distribution; Wakeford was not fetched (RECALLED). The [R0710] bulk-Z result for 147 giants in part 05 is for giant planets and does not reach sub-Neptunes.

**2. Physical bracket.** Z (mass fraction) from "k times solar" with Z_sun = 0.0139 ([R0168] via part 05): Z = k Z0 / (1 - Z0 + k Z0). COMPUTED: k = 1 gives 0.014, 4 gives 0.053, 10 gives 0.124, 100 gives 0.585 (matches the 58 percent of TOI-270 d), 162 gives 0.70, 1000 gives 0.93, 2291 gives 0.97. The part-09 bin edges 0.01 to 0.6 therefore cover 0.7x to 106x solar; the upper edge excludes K2-18 b's 162x upper limit, GJ 1214 b (above 100x, up to about 1000x) and GJ 3090 b. Fully miscible envelopes can reach Z = 0.58 or more [R0339].

**3. Verdict.** (b). About eight sub-Neptune retrievals and one 15-planet exo-Neptune retrieval, with different methods and a stated lack of trend evidence. Not a population distribution. Also, atmospheric Z (what spectra measure) differs from bulk envelope Z (what e5's mean molecular weight formula uses) when the envelope is not well mixed or the magma ocean exchanges volatiles ([R0339, R3299] miscible case). Do not draw Z log-uniform as if it were occurrence-weighted.

**4. Knobs.** Centre, width and upper edge of the log-uniform draw: lower edge 0.01 (solar-like, V1298 Tau b at 4x gives 0.053 and is an outlier low; Z below 0.01 not excluded, K2-18 b lower limit 0.6x), upper edge raised to about 0.95 (1000x) or kept at 0.6 with a stated cut; width target: the Brande scatter of about 0.5 scale height is the only population-level spread (READ), but it is an aerosol-model quantity, not a Z spread. No target value with a published uncertainty exists for the mean of log Z.

**5. Consistency test.** Compare the simulated Z of the planets matching JWST targets (mass, radius, Teq) with the table above; tolerance is each retrieval's 1 sigma (for example TOI-270 d +8/-12 percent); a trend test needs about 45 JWST planets (Sun 2024) and is not yet possible.

**6. No-fixed-content flags.** The Solar System giants' carbon-versus-mass trend that [R2802] plots is the only population trend extrapolated to sub-Neptunes, and it is Solar System only (Uranus, Neptune, Saturn, Jupiter). Do not use it to set e5; the review itself finds the sub-Neptune carbon abundances do not follow it strongly.

**7. Group for a builder.** Builder C: e5 and e6 (CH4 vs CO2 carbon carrier), because they share the JWST table and the same retrievals; e6 is outside this group.

## 3. Sources verified (id, N, access)

All fetched with `curl` from `export.arxiv.org/pdf/<id>`, text via `pdftotext`; all arXiv content is under the arXiv non-exclusive licence unless the paper states CC (not individually checked; data licences were not stated in any paper read).

| Candidate in prompt | Verified id | Check | Status |
|---|---|---|---|
| [R1478] CKS III | 1703.10375 | correct; N 2025 / 900 | READ |
| [R1477] | 1805.01453 | correct; CKS VII; 907 planets | READ |
| [R3500] | 2201.10020 | CKS X; 970 planets (the id 2205.03385 first guessed was an unrelated paper) | READ |
| [R2020] | 1902.01417 | FGK ABC occurrence, DR25 + Gaia DR2; code on GitHub (ExoplanetsSysSim.jl, SysSimData) | READ |
| [R2021] (M dwarfs) | 2002.02573 | authors Hsu, Ford, Terrien (not Ragozzine, Ashby); 1,746 targets, 89 candidates; MIT-licensed code | READ |
| [R0572] | 2010.14812 | exists (rocky HZ occurrence); not used for the valley; code `github.com/stevepur/DR25-occurrence-public` (README read) | READ (abstract) |
| Owen and Wu 2017 | 1705.10810 | correct | READ |
| [R4504] | 1710.05398 | correct; 117 planets | READ |
| [R3797] | 2007.11006 | correct | READ |
| [R1734] | 1811.03202 (the id 1907.03745 is unrelated); 2020 stellar-dependence paper is 1907.03732 | read | READ |
| [R4880] PNAS | 1906.04253 | exists; gives no water-world occurrence number | READ |
| [R2774] | 2209.03871 | N = 34 | READ |
| [R3384] | 1911.04745 | two populations (rocky, volatile-rich); M-R laws M = 0.9 R^3.45 and 1.74 R^1.58; N not extracted | READ (abstract) |
| [R3800] | 1407.4457 | R_thresh 1.48 (95 percent < 1.59), 27 RV planets | READ |
| [R0300] | 2010.11867 (the id 2010.13219 is unrelated) | review | READ (abstract) |
| [R3252] | 1911.03582 | N = 1130; Q_rocky 0.20 | READ |
| [R0039] | 2010.01074 | TRAPPIST-1 | READ (abstract) |
| [R3430] | 2406.04311 (Parc, Bouchy, Venturini, Dorn, Helled) | 46 M-dwarf planets | READ |
| [R4528] / 2024 | 2008.05513 / 2404.01967 (the 2024 letter is by Venturini, not Parc) | pebble models | READ |
| [R0828] | 1912.02170 (the id 1912.02188 is unrelated) | read | READ (abstract) |
| Ho and [R1963]; Ho 2024 | 2301.04062; 2401.12378 | read | READ |
| [R0353, R0354] | 2005.14671; 2302.00009 | read | READ |
| [R3796] (3-D gap) | 2105.03443 | needs about 5000 planets | READ |
| [R3798, R3801] | 2301.04321; 2409.17394 | read | READ |
| [R0714]; 2026 | 2310.03593; 2602.11923 | read; code `github.com/GijsMulders/Genesis`, `github.com/arcunique/GPS` named in the papers | READ |
| [R0978] | 2503.02451 | N = 60 | READ |
| [R1568] | 2602.23364 | preprint | READ |
| [R0601] | 2401.04380 | read; KOBE code on GitHub | READ |
| JWST: [R0051, R0255, R0304, R0339, R0340, R2316, R2805, R3299, R3517, R4731] | 1909.04642; 2309.05566; 2401.11082; 2305.06240; 2406.04450; 2403.03325; 2410.03527; 2504.20428; 2507.08837; 2510.07367 | abstracts and key lines read | READ |
| [R4663]; Sun 2024 (Sun, Wang, Welbanks, Teske, Buchner); [R0516, R0932, R2802] | 1912.04904; 2402.08292; 1708.00016; 2310.07714; 2509.19247 | read | READ |
| NASA Exoplanet Archive Kepler "Simulated Data" page | URL `exoplanetarchive.ipac.caltech.edu/docs/KeplerSimulated.html` returned 200 | page read; no licence text on it | DATA |
| VizieR J/AJ/163/179 (CKS DR2) | front page returned 200; README download returned 500 | table sizes not read | DATA (front page) |
| CKS-Cool README | 404 on raw GitHub | NOT REACHED | NOT REACHED |
| [R3443, R4590] (SPRIGHT), [R4352] | not fetched | | NOT REACHED |

## 4. Boss must supply

1. A decision on the mechanism for e1/e2 (XUV energy-limited, as in the repo, versus core-powered) and a radius law for retained envelopes ([R2741] is the usual choice; I did not fetch it). Without a radius law no valley statistic can be computed.
2. For b4/e3/e4: the formation model (snow line plus migration rule) that sets the water-rich share; the public record gives only the model-dependent spread above. The Genesis library or the Burn/Venturini pipelines are the nearest public candidates; their licences and data tables were not checked.
3. The CKS DR2 and CKS-Gaia planet tables (VizieR J/AJ/163/179; Fulton and Petigura Table 4) and the DR25 injection efficiency tables, downloaded and checked in: the proxy here did not return the table files.
4. The XUV history (E_cum versus stellar mass and age) used by the repo, with its source; this lane did not read the activity document behind E_cum = 3.1e16 J/m2.
5. TEPCat or PlanetS tables of precise masses and radii for the mass-radius test (Parc catalogue; TEPCat), with a licence check.
6. For e5: an occurrence-weighted metal-fraction distribution does not exist; the owner must accept the bracket plus per-planet TEST points, or wait for JWST population studies (about 45 planets).

## 5. Recalled / to verify

- Z_sun = 0.0139 is READ in part 05 [R0168], taken second-hand here.
- E_cum = 3.1e16 J/m2 and the core radius law R = M^0.27 are from the repo doc (project XUV model [C]); not verified by me.
- [R2741] radius law; [R2740] eps range (2ND-HAND via part 04).
- [R4590] (HAT-P-26 b) details: RECALLED, not fetched.
- [R3384] sample size N: not extracted.
- [R0715] TEPCat N and its peer-review status; [R1568] review status.
- Cloutier and Menou sample counts (275 + 53): 2ND-HAND from [R3500]; verify in the paper's body.
- [R4335] cap and [R2804] hycean definition: taken from part 05 and part 09 (READ there, not re-read here).
- Licences of the CKS and DR25 tables: not stated in any paper read; verify at the archive before redistribution.
