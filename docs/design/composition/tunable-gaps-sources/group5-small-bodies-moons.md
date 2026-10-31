# Group 5: comets, small bodies, moons and surfaces (source check for tunable gaps)

Read-only research lane, 2026-10-17 (EDT clock times only). Labels as in _worker-rules.md: READ, 2ND-HAND, DATA, COMPUTED, SNIPPET, RECALLED, NOT REACHED. Nothing in the repo was touched.

## 1. Summary and verdict table

The public record is rich for one thing only: Solar System coma abundances of about 20 to 30 comets (k4, k5, k7 partly, T16). For those, a seeded log-normal with published centre, width and standard errors can be tuned and tested. Everything else is a handful of Solar System objects (17 KBO densities, the Solar System moons), one zone mix of one asteroid belt (T13), a modelling switch that is not a probability (T15, E2, E3), or has no ensemble (T14 portable weights, T17, T18, T21). Extrasolar data break the Solar System comet prior: two interstellar comets lie 2 to 4 sigma outside the Solar System CO/H2O and CO2/H2O distributions, so a per-system offset knob is needed, and two objects are all that constrain it. No exomoon is confirmed (Nature, 22 July 2026: "no exomoon has ever been confidently detected"). Every row has the Solar System as its only or main constraint and is flagged under the no-fixed-content rule.

| Row | Quantity | Verdict | Reason | Status |
|---|---|---|---|---|
| k1, k2, k3 | Coma class CO2-dominated / CO-dominated / comparable | (b) | N=25, counts 13/9/3, but class depends on r_h and on an upper-limit admission rule; better one continuous log(CO/CO2) with an r_h term | READ + COMPUTED |
| k4 | CO2/H2O log-normal | (a) Solar System; (b) with extrasolar offset | N=21: geo-mean 10.9 %, sigma 0.51 dex; 3I at 3.6 sigma | DATA + COMPUTED |
| k5 | CO/H2O log-normal | (a) Solar System; (b) with offset | N=16: 2.5 %, 0.65 dex; 2I and 3I at 1.8 to 2.8 sigma | DATA + COMPUTED |
| k7 | Single-comet species | (a) CH4, C2H6, H2CO, H2S, OCS; (c) N2, O2, Ar, Kr, Xe, Ne, SO2 | populations of 6 to 18 comets exist for the first group; the rest are 67P only | DATA + READ |
| T13 | Taxon shares in a/a_snow | (b) | one zone mix of one belt (34,503 asteroids); scatter is a model | READ |
| T14 | Within-taxon meteorite weights | (b) | fall shares and NEO spectral shares disagree (LL 10 % vs about half) | READ |
| T15 | Refractory organics p = 0.5 | not a probability | source runs two extreme models; real quantity is a continuous organic carbon fraction | READ |
| T16 | Comet volatile abundances | (a), 10 species | same ensemble as k4, k5; Willacy 2025 adds 12 curated comets | DATA + READ |
| T17 | Moon origin weights | (c) | rules exist, no population gives weights; 0 confirmed exomoons | READ |
| T18 | Moon ice fraction ramp width | (c) | two READ population models disagree (Sasaki vs Cilibrasi) | READ |
| T19 | Dwarf planet / KBO rock fraction | (b) | 17 binary KBO densities; rock fraction degenerate with porosity | READ |
| T21 | Grain state | (c), derive by physics | one extrasolar object (LHS 3844 b) rules out fresh powder | READ |
| E2 | Belt density | not a probability; derive from belt mass | debris-disc distribution exists for cold belts only | READ |
| E3 | Field family frequencies 30/20/20/10/8/5/4/3 % | (c) | no ensemble; design vocabulary | none |

## 2. Detail

### 2.1 Comets: k1 to k5, k7, T15, T16

**Ensemble (point 1).**

| Data set | Id | N | Measured and selection | Access, licence | Status |
|---|---|---|---|---|---|
| Harrington Pinto, Womack, Fernandez, Bauer 2022 | arXiv:2209.09985 | 25 comets with simultaneous CO and CO2 (Table 4); Table 3 at r_h <= 2.5 AU: 32 CO, 23 CO2, 18 CO+CO2 | production rates 0.71 to 6.18 AU; CO upper limits kept only if Q(CO)/Q(CO2) <= 2/3; Spitzer/NEOWISE photometry excluded; fluorescence efficiency ratio 11.6 | arXiv; licence not stated | READ |
| Robinson et al. 2024 database | arXiv:2403.15644; doi 10.7488/ds/7723, df_comp_all.csv | 909 rows, 96 designations (earlier lane counted 69, unresolved); 37 comets have a parent-volatile X/H2O, 21 have five or more, 13 have eight or more; at r_h <= 2.5 AU: 30, 18, 14 | X/H2O for CO, CO2, CH4, C2H6, CH3OH, HCN, NH3, H2CO, H2S, C2H2; only comets with a published nucleus radius; sources unhomogenised (Dello Russo 2016: 97 rows) | Edinburgh DataShare, CC BY 4.0 for the compilation; source papers may differ | DATA |
| Willacy et al. 2025 | arXiv:2509.20530 | 12 comets with >= 4 COH species, 12 with >= 4 N species, 7 with >= 4 S species; 9 of 12 are OCC | curated tables, all within 2 AU | arXiv | READ |
| Rubin et al. 2019 (67P) | arXiv:1907.11044 | 1 | in situ bulk ratios incl. N2, O2, Ar, Kr, Xe | arXiv | READ (prior lane) |
| Dello Russo 2016; Mumma and Charnley 2011; Bockelee-Morvan and Biver 2017; A'Hearn 1995 | doi 10.1016/j.icarus.2016.05.039; 10.1146/annurev-astro-081309-130811; 10.1098/rsta.2016.0252; 10.1006/icar.1995.1190 | 30 comets (title); more than 100 surveyed; review; 85 comets, optical only | Crossref metadata and abstracts only; texts refused | publisher | NOT REACHED (metadata READ) |
| A'Hearn 2012 | not found in Crossref | | | | RECALLED |

Answer to question 1 (count): a multi-species parent-volatile table exists for about 20 to 30 comets. Only 8 to 9 comets have CO and CO2 together with five other species, so the joint distribution rests on about one dozen objects.

**Class fractions (k1 to k3).** HP22 Table 2 (READ): CO2-dominated 13 of 25 (52 %), CO-dominated 9 (36 %), comparable 3 (12 %); Clopper-Pearson 68 % intervals 0.40 to 0.64, 0.25 to 0.48, 0.055 to 0.22 (COMPUTED). The repo weights 0.5, 1/3, 0.1 are the rounded abstract words and sum to 0.93; the printed counts replace them.

| Cut (COMPUTED from HP22 Table 4 unless noted) | N | CO2-dom | CO-dom | comparable |
|---|---|---|---|---|
| As published (READ) | 25 | 13 | 9 | 3 |
| Lowest-r_h row per comet, r_h <= 2.5 AU | 18 | 10 (0.56, 0.41 to 0.69) | 4 (0.22, 0.12 to 0.36) | 4 (0.22) |
| r_h <= 3.5 AU | 24 | 12 | 8 | 4 |
| JFC only (READ) | 9 | 8 (0.89, 0.68 to 0.98) | 1 | 0 |
| OCC only (READ) | 14 | 5 | 6 | 3 |

1. Distance: every coma beyond about 3.5 AU is CO-dominated (29P, C/1995 O1, C/2006 W3, C/2016 R2, 67P at 3.5 AU; READ). The CO-dominated share reflects where comets were observed as much as what they are made of.
2. Selection: 9 of the 13 CO2-dominated entries rest on an upper limit on CO (22P, 46P, 81P, 88P, 144P, C/2006 OF2, Q1, C/2007 N3, Q3; COMPUTED). Limits with ratio above 2/3 were not admitted, so the CO2 class is inflated by an amount that cannot be corrected without the excluded list.
3. The only class effect surviving a permutation test is CO (28 % of variance, p = 0.035; G1 lane, COMPUTED).

Recommendation: draw log10(CO/CO2) at the comet's own r_h and derive the three labels from HP22's thresholds (<= 2/3, >= 4/3).

**Correlated log-normal versus independent (question 1b).** Robinson database, r_h <= 2.5 AU, per-comet median of detections (DATA, COMPUTED):

| Species | N | Geo-mean % H2O | sigma dex (SE) | SE of mean | Range % |
|---|---|---|---|---|---|
| CO | 16 | 2.5 | 0.65 (0.12) | 0.16 | 0.30 to 26 |
| CO2 | 21 | 10.9 | 0.51 (0.08) | 0.11 | 0.14 to 37 |
| CH4 | 12 | 0.68 | 0.33 (0.07) | 0.09 | 0.15 to 2.2 |
| C2H6 | 18 | 0.41 | 0.33 (0.06) | 0.08 | 0.09 to 1.9 |
| CH3OH | 21 | 1.96 | 0.33 (0.05) | 0.07 | 0.21 to 4.5 |
| HCN | 22 | 0.14 | 0.28 (0.04) | 0.06 | 0.03 to 0.54 |
| NH3 | 15 | 0.57 | 0.23 (0.04) | 0.06 | 0.13 to 0.9 |
| H2CO | 15 | 0.25 | 0.35 (0.07) | 0.09 | 0.09 to 1.1 |
| H2S | 9 | 0.71 | 0.39 (0.10) | 0.13 | 0.09 to 1.5 |
| C2H2 | 12 | 0.10 | 0.32 (0.07) | 0.09 | 0.03 to 0.34 |

HP22 Table 3 medians, CO 3 +- 1 % (N=32) and CO2 12 +- 2 % (N=23), agree (READ). The CO2 width depends on the cut: 0.27 dex at r_h <= 2 AU (G1 lane, N=14), 0.51 dex at <= 2.5 AU.

Correlations (35 pairs with N between 6 and 22): mean r = 0.13 (23 positive, 12 negative); 9 pairs have p < 0.05 where 1.8 are expected, 1 has p < 0.01 (COMPUTED, Fisher z). CO versus CO2 is -0.87 at r_h <= 2 AU (N=8, p = 0.003) but -0.59 at <= 2.5 AU (N=9, p = 0.10), so it is not robust to the cut. A full covariance cannot be fitted: the best complete set is 9 comets for 7 species. Bartlett sphericity tests on 4 to 7 species (N = 9 to 17) give p from 0.002 to 0.8, four of seven below 0.07. A plausible common cause is the shared divisor Q(H2O). Verdict: correlated is better than independent, but only as a one-factor model (a shared per-comet offset on all X/H2O with loading about 0.13 in correlation, plus independent per-species scatter); a free covariance is not identifiable. Keep the G1 lane's three-step draw ((CO+CO2)/H2O median 18 +- 4 %, then a CO fraction) for the CO versus CO2 trade-off.

**Physical bracket (point 2).** Ice lines set the nucleus inventory; Eistrup, Walsh, van Dishoeck 2019 (arXiv:1907.11255, G1 lane READ): 15 comets lie on a trail tracing the CO ice line from about 30 AU early to 12 AU later. Coma is not nucleus: Cordiner 2025 shows coma CO2/H2O falls with r_h.

**Extrasolar constraint (question 2).**

| Object or system | Value | Source | Status |
|---|---|---|---|
| 2I/Borisov, ALMA | CO/H2O 35 to 105 %; HCN/H2O 0.06 to 0.16 % (Solar System-like) | Cordiner 2020, arXiv:2004.09586 | READ |
| 2I/Borisov, HST+Swift | CO/H2O 130 to 155 %; at least 173 % from HST; Solar System inner comets 0.2 to 23 %, typical 4 % | Bodewits 2020, arXiv:2004.08972 | READ |
| 3I/ATLAS, JWST NIRSpec, 3.32 AU | CO2/H2O 7.6 +- 0.3; CO/H2O 1.66 (3.7e26 over 2.23e26, COMPUTED); 4.5 sigma above the Solar System r_h trend | Cordiner 2025, arXiv:2508.18209 | READ |
| 3I/ATLAS other epochs, CO2/H2O | > 4.4 (3.3 to 3.1 AU); 0.21 +- 0.05 (2.0 to 2.2 AU); 3.16 and 7.28 (2.19, 2.54 AU); 0.31 to 2.39 (2.87 AU, [O I]) | Table 4 of arXiv:2603.25002 (cites arXiv:2512.07318, 2601.22034) | READ (table) |
| 1I/'Oumuamua | no coma | prior lane | READ |
| Fomalhaut / HD 181327 / beta Pic belts | CO+CO2 ice mass fraction 4.6 to 76 % / 0.8 to 3.7 % / below 32 %, "consistent with Solar System comets" | Matra 2017, arXiv:1705.05868 (id and title correct) | READ |
| beta Pic gas | HCN/(CO+CO2) below 2.5 %, low end of Solar System range | Matra 2018, arXiv:1801.03945 | READ |
| CO-bearing debris discs | 6 confirmed exocometary, 14 more gas detections of unconfirmed origin; 15 CO systems in 2018 | arXiv:1904.02715, 1801.03945 | READ |

COMPUTED z-scores against the Solar System fits: 3I CO2/H2O 7.6 is 3.6 sigma high, 3.16 is 2.9 sigma, 0.21 is 0.6 sigma; 3I CO/H2O 1.66 is 2.8 sigma; 2I CO/H2O is 1.8 to 2.8 sigma. A Solar System-only prior gives a joint one-sided probability of order 1e-5 to 1e-4 for both ISO CO/H2O values. The comparison is confounded by distance (3I at 3.3 AU against a fit at <= 2.5 AU) and the same object spans a factor 36 in CO2/H2O across epochs. Verdict: the Solar System population does not describe the two interstellar comets, but it is not cleanly refuted either. Exocomet belts do not break it: the Fomalhaut interval alone spans 1.2 dex, so belts test the prior only to about 1 dex.

Offset knob (PARAMETRISATION): per system draw delta_sys ~ N(0, sigma_sys) in log10 units, common to CO, CO2 and the other hypervolatile over H2O ratios (the shared H2O divisor makes them move together), then per-comet scatter. Observed offsets: +1.1 to +1.8 dex (2I CO/H2O), +1.8 (3I CO/H2O), +1.8 (3I CO2/H2O). With two objects, sigma_sys can only be bracketed: 0 is disfavoured at about 1e-4 if the two are typical, an upper edge near 1.5 dex; ISO detection is biased toward active, volatile-rich bodies.

**Verdict (point 3).** See table in section 1. Tolerances below come from the SE columns.

**Knobs (point 4).**

| Knob | Statistic and target (uncertainty) | Range | Label |
|---|---|---|---|
| mu, sigma CO2 | 10.9 % (0.11 dex), 0.51 dex (0.08); HP22 median 12 +- 2 % | r_h <= 2.5 AU, 21 comets | DATA + READ |
| mu, sigma CO | 2.5 % (0.16 dex), 0.65 dex (0.12); HP22 3 +- 1 % | same, 16 comets | DATA + READ |
| class shares | 13/9/3 of 25, binomial intervals above | per r_h bin | READ + COMPUTED |
| (CO+CO2)/H2O | median 18 +- 4 % (N=18) | r_h <= 2.5 AU | READ |
| other species sigma | 0.23 to 0.39 dex | same | DATA |
| common factor | mean pair correlation 0.13 | N=6 to 22 per pair | COMPUTED |
| delta_sys, sigma_sys | ISO offsets +1.1 to +1.8 dex | N=2 | READ + COMPUTED |
| carbon-chain-depleted flag | 0.37 JFC, 0.185 long-period, from about 250 optical comets (Cochran 2015, arXiv:1507.00761) | optical classes | 2ND-HAND |

**Consistency test (point 5).** Simulate comets at the observed r_h distribution, apply the HP22 admission rule, and compare class counts with 13/9/3 within binomial intervals and each species' geo-mean and sigma within 2 SE (CO2 +-0.22 dex on the mean, CO +-0.32). A selection model is needed (admission rule, r_h distribution, CO versus CO2 detection limits). For the ISOs, simulate systems with the offset knob and require that 2 of 2 draws above CO/H2O 0.3 are not improbable.

**Flags (point 6).** For k1 to k5, k7, T16 the only constraint on within-system scatter is the Solar System (N = 9 to 32); extrasolar objects constrain the offset only.

**k7.** Replace 67P-centred values by populations where they exist: CH4 0.68 % (N=12), C2H6 0.41 % (N=18), H2CO 0.25 % (N=15), H2S 0.71 % (N=9). Willacy 2025 Table 4 (READ): OCS 0.034 to 0.4 % in 6 comets; SO2 detected in 2 comets only (0.127, 0.23). N2 (0.089 %), O2 (3.1 %), Ar, Kr, Xe, Ne are 67P only (Rubin 2019); O2 detection probability cannot be set from one comet.

**T15.** Marboeuf et al. 2014 (arXiv:1407.7271, READ) "study two extremes cases": interstellar refractory organics survive (organic carbon 60 +- 15 % of total carbon in molecular clouds, C:O:N = 1:0.5:0.12, cited from Pollack 1994, 2ND-HAND) or are destroyed. There is no probability. Observed dust is never organic-free (67P dust about half organics; arXiv:2501.17864, SNIPPET). Replace the Bernoulli by a continuous f_org ~ N(0.60, 0.15) clipped to [0,1].

### 2.2 Asteroids and meteorites: T13, T14, E2, E3

| Data set | Id | N | Measured, selection | Status |
|---|---|---|---|---|
| DeMeo and Carry 2013 | arXiv:1307.2424 | 34,503 SDSS asteroids (17,575 with albedo), down to 5 km | bias-corrected taxon distribution by number, area, volume, mass; Table 6 zone mix; no uncertainties; cells with fewer than 20 objects flagged (A all zones; B Hildas; C Trojans; D, K, L, S Cybeles; V middle and outer; X Trojans) | READ |
| DeMeo and Carry 2014 | arXiv:1408.2787 | abstract | belt more diverse with size and distance than known; mixing | READ (abstract) |
| Carry 2012 | arXiv:1203.4336 | 287 bodies | densities by class | READ (abstract) |
| MITHNEOS (Binzel 2019) | arXiv:2004.05090 | NEOs | ordinary chondrite NEOs about half LL (194 objects); falls: L 47, H 43, LL 10 % of ordinary chondrite falls by argon ages; ordinary chondrites about 81 % of falls | READ |
| Meteorite falls | via arXiv:1310.2000 | about 1,100 falls | H 34.1, L 38.0, LL 7.9 % of all falls | 2ND-HAND |
| Mulders et al. 2015 | arXiv:1505.03516 | model | 0.4 dex accretion-rate spread gives about 0.2 dex snow-line spread at fixed mass and age | READ |

Candidate corrections: "DeMeo et al. 2019" does not match "Connecting asteroids and meteorites with visible and near-infrared spectroscopy", which is DeMeo, Burt, Marsset, Polishook, Burbine, Carry, Icarus 380, 114971 (2022) (Crossref). Mulders 2015 supports 0.2 dex only as a model consequence, not a measured scatter.

Extrasolar ensemble (question 3): none for taxon mixes. Polluted white dwarfs belong to Group 3; cross-reference only: Putirka and Xu 2021 (doi 10.1038/s41467-021-26403-8, abstract READ) find precisely measured Mg, Si, Ca, Fe for a few polluted WDs, "most are exotic in composition and mineralogy", beyond the spread of more than 4,000 nearby main-sequence stars. That is rock-type data, not Bus-DeMeo shares, so it cannot tune T13, but it says the mix should not be assumed Solar System-like. A white-dwarf study (arXiv:2410.06248, READ) confirms exoasteroids are inferred only from WD pollution and photometry.

**Bracket.** Taxon follows formation temperature (dry inside, hydrated C in the middle, anhydrous ice-rich P and D beyond; DeMeo and Carry 2013); x = a/a_snow is a PARAMETRISATION, untested for other stars; DeMeo and Carry 2014 show position alone does not set taxon.

**Verdict.** T13 (b): one system, no uncertainties, which is exactly the realisation R1 wants to replace. T14 (b): within-taxon weights exist only as fall shares, which disagree with NEO spectral shares ("an unsolved problem", MITHNEOS), so they depend on a delivery model.

**Knobs and test.** sigma_snow = 0.2 dex (Mulders, model); zone mix of Table 6 as target; H:L:LL 0.43:0.47:0.10 (falls, ordinary chondrites). Test: pooled simulated mix evaluated at the Solar System x equals Table 6 within simulated scatter; this is self-consistency, not an independent test. Flags: the only constraint is one belt.

**E2 (belt density).** Not a probability; derive from belt mass per area. Extrasolar ensembles bear on cold belts: Sibthorpe 2018 (arXiv:1803.00072, READ): 275 F to K stars, 47 discs, 17.1 +2.6/-2.3 %, fractional luminosity near 1e-5 at blackbody radii 7 to 40 AU. Warm dust: Kennedy and Wyatt 2013 (arXiv:1305.6607, READ): dustiest old systems 1 in 10,000, about 1 % for ages below 120 Myr. They give a distribution of belt dust amount, but no mapping to "dense, sparse, typical". Verdict (c) as a probability; derive as thresholds on drawn belt mass.

**E3 (families 30/20/20/10/8/5/4/3 %).** No public ensemble; the families are a design vocabulary and no survey defines standalone asteroid fields. The only physical correspondence is taxon mass shares in the Solar System belt (C 52.5, P 11, S 8.4, V 9.6, M 3.3 %; Table 5, READ). Verdict (c): keep as a flagged stand-in or derive the family from the parent belt's T13 shares.

### 2.3 Moons: T17, T18

| Data set | Id | N | Result | Status |
|---|---|---|---|---|
| Canup and Ward 2006 | doi 10.1038/nature04860 | 3 systems | satellite system mass about 1e-4 of the planet; "could limit the largest moons of extrasolar Jupiter-mass planets to Moon-to-Mars size" | READ (abstract) |
| Canup and Ward 2009 | arXiv:0812.4995 | scaling | M_T/M_p about 2.5e-4 (alpha/f / 3e-5)^(1/3) | READ (prior lane) |
| Sasaki, Stewart, Ida 2010 | arXiv:1003.5737 | 67 Saturn-type runs; 39 Jupiter-type runs with four large moons | inner two rocky, outer two icy (Jupiter); one dominant icy body (Saturn) | READ |
| Heller and Pudritz 2015 | arXiv:1410.5802 | simulations | icy moons larger than the smallest known exoplanet at about 15 to 30 R_J around super-Jovians | READ (abstract) |
| Cilibrasi et al. 2018 | arXiv:1801.06094 | population synthesis, Jupiter | 85 % of the moon population icy; Galilean-like systems in 20 % of cases | READ |
| Ronnet and Johansen 2020 | arXiv:1912.04765 (candidate "Ronnet et al.") | model | planetesimal capture and ablation feed pebble accretion | READ (abstract) |
| Quintana et al. 2016 | arXiv:1511.03663 | 280 N-body runs, 164 Earth-analogues | median 3 giant impacts of Moon-forming energy; last one at median 43 Myr; 15 of 164 stripped an atmosphere | READ |
| Malamud et al. 2020 | arXiv:1904.12854 | 28 SPH collisions | 10 of 28 made discs above 0.1 Earth masses, 2 above 0.2; 10 to 55 % of disc becomes moon | 2ND-HAND |
| Madeira et al. 2025 | arXiv:2507.04954 | N-body discs | disc above 0.03 planet masses gives one dominant satellite; large scatter | READ (abstract) |
| Kipping et al. 2022 | arXiv:2201.04643 | 70 cool giant planets | one candidate (Kepler-1708 b-i, 4.8 sigma, 1 % false-positive probability), occurrence rate unknown | READ |
| Teachey and Kipping 2018 | arXiv:1810.02362 | 1 candidate | Kepler-1625 b-i; contested (Heller 2019, Kreidberg 2019,) | READ |
| Hoy et al. 2026 | doi 10.1038/s41586-026-10751-w | 1 candidate | exosatellite of brown dwarf CD-35 2722 B, minimum mass about 0.9 Jupiter masses, period about 170 d; "no exomoon has ever been confidently detected" | READ (abstract) |

Question 4: no published ensemble defines origin-type weights. Giant planets have a mass scaling (3 systems) and Jupiter-only population syntheses. Terrestrial planets have giant-impact frequencies from N-body formation (Quintana), but whether an impact leaves a moon-forming disc depends on impactor mass, speed and angle (10 of 28 SPH cases), not on a frequency. Capture (Triton, doi 10.1038/nature04792) and ring origin (doi 10.1126/science.1226477) are single Solar System cases (Crossref titles READ). Candidates Jacobson and Morbidelli, Citron and Chambers were not located. Batygin and Morbidelli 2020 (ApJ 894, 143) is SNIPPET only.

Verdict. T17 (c): the gate rules are physics but no weight is published. Derive instead: P(impact moon) from the impact-energy frequency (Quintana-type) times a disc-yield function (SPH), and disk moons from the 1e-4 scaling. T18 (c): Sasaki (ice line r_ice = 20 R_J (M_p/MJ)^(2/3) (tau/5 Myr)^(-1/3), COMPUTED, fitted around Jupiter) and Cilibrasi (85 % icy) are both READ and disagree; neither overrides, so the spread 0.5 to 0.85 icy is the honest bracket. Neither paper gives a ramp width.

Knobs and test: M_T/M_p = 2.5e-4 (alpha/f / 3e-5)^(1/3), 0.8e-4 to 6.4e-4 over alpha/f 1e-6 to 5e-4 (COMPUTED from Eq. 4.15); simulated giant-planet systems must reproduce about 1e-4 within a factor 2, which uses only the Solar System (calibration point). P(disc > 0.1 Earth masses | giant impact) = 10/28 = 0.36 (CP68 0.24 to 0.50, COMPUTED), from a parameter-space suite rather than a population. Flags: the only constraint is the Solar System and simulations.

### 2.4 KBOs and dwarf planets: T19

| Data set | Id | N | Result | Status |
|---|---|---|---|---|
| Brown 2013 | arXiv:1311.0553 | 17 objects in Fig. 3 (binaries only; masses from orbits) | densities 0.5 to 2.6 g/cm3; "smaller objects have low densities, larger increasingly higher"; 2002 UX25 (650 km) 0.82 +- 0.11; small-KBO weighted mean about 0.6, porosity-free about 1.2 (one third rock) | READ |
| Barr and Schwamb 2016 | arXiv:1603.06224 | 7 bodies | radius (km), density: Triton 1353 / 2.061; Eris 1163 / 2.52; Haumea 620 / 2.6; Orcus 430 to 450 / 1.605 to 1.747; Pluto 1187 / 1.86; Charon 606 / 1.702; Quaoar 555 / 1.99 +- 0.46; "no clear correlation between size and density" among them | READ |
| Bierson and Nimmo 2019 | doi 10.1016/j.icarus.2019.01.027 | 17 KBOs | constant 70 % rock fits 14 of 17 within 2 sigma; porosity explains size trend; formation after 4 Myr; text behind a captcha | SNIPPET; NOT REACHED |

COMPUTED from Barr and Schwamb (ice 917, rock 2700 or 3500 kg/m3): rock fraction of the seven dwarf planets, mean 0.82 (SD 0.12) or 0.74 (SD 0.10), range 0.61 to 0.98; log density versus log radius slope 0.14, r = 0.35. For small bodies (Salacia 1.29, porosity-free 1.2) the rock fraction is 0.32 to 0.44. 2002 UX25 is below pure ice density, so density alone cannot give rock fraction without a porosity law. Candidate corrections: "Brown 2012" is arXiv:1112.2764 (not analysed); my search returned Grundy 2009 (arXiv:0812.3126, six binaries, assumed densities).

Bracket: x between 0 and 1 through 1/(x/rho_rock + (1-x)/rho_ice) before porosity; asteroid compaction near 350 km (Brown 2013). Verdict (b): both READ sources read the same 17 densities differently (rock fraction rising with size versus constant 70 % with porosity doing the work); the data cannot decide. Verified numbers win over "about 70 %, spread 0.58 to 0.82": dwarf planets 0.74 to 0.82 mean, small KBOs 0.3 to 0.45. Knobs: dwarf-planet rock fraction 0.78 +- 0.12 (N=7); porosity law with a transition near 350 km. Test: simulated density-diameter relation within the scatter of the 17 objects, with a binary-detection selection model. Flag: Solar System only, no extrasolar KBO density.

### 2.5 Grain state: T21 (question 6)

| Data set | Id | N | Result | Status |
|---|---|---|---|---|
| Paragas et al. 2025 | arXiv:2502.04433 | 11 samples, 3 textures | texture changes albedo up to a factor 7 for one sample; similar compositions differ by up to 2; slab versus powder changes a LHS 3844 b analogue dayside temperature by about 70 K; data at Zenodo 17986129 | READ |
| Zieba et al. 2026 | arXiv:2605.00100 | LHS 3844 b, JWST 5 to 12 micron | Bond albedo 0.14 +0.13/-0.14; slabs within about 2 sigma, crushed marginal, fresh powders excluded except hematite; old space-weathered powder also fits; no CO2 above 100 mbar | READ |
| Kreidberg et al. 2019 | arXiv:1908.06834 | LHS 3844 b | no atmosphere above 10 bar; Bond albedo below 0.2 | READ (abstract) |
| Hu, Ehlmann, Seager 2012; Mansfield 2019 | arXiv:1204.1544; 1907.13150 | library; method | Hu is the room-temperature library Paragas extends; Mansfield uses high albedo as an atmosphere test and says nothing on grain statistics | READ (titles) |
| Regolith model | arXiv:1802.03664 | Apollo 15 site | stochastic model; growth slower than crater equilibrium; thickness varies laterally | READ (abstract) |

Is there an extrasolar constraint? One object (LHS 3844 b) that picks between a young slab and an old powdered, space-weathered surface. Is it derivable? Largely: grain state follows from surface age, impactor flux, resurfacing, heat flux and (with an atmosphere) wind and weathering. Only a Solar System calibrated model exists (one lunar site, Solar System crater distribution). Verdict (c) as a weight; derive: P(powder) rises with age times flux, P(slab) with recent resurfacing, and the Paragas library maps texture to albedo. Test point: LHS 3844 b brightness-temperature ratio 0.96 +0.03/-0.04 (a calibration point, not a prior). Flag: the Moon and asteroids are the only other constraints.

### 2.6 Group for a builder (point 7)

| Group | Rows | Shared data and fit | Test |
|---|---|---|---|
| G5-A Comet coma sampler | k1 to k5, k7, T15, T16 | Robinson database, HP22 Table 4, Willacy 2025; log-normals, one-factor correlation, r_h-dependent CO/CO2, ISO offset | class counts, geo-means with SE, ISO exceedance |
| G5-B Belt composition | T13, T14, E2, E3 | DeMeo and Carry Table 6, fall shares, Mulders, DEBRIS; stand-ins | self-consistency |
| G5-C Moons and KBOs | T17, T18, T19 | Canup and Ward, Sasaki, Cilibrasi, Quintana, Brown, Barr and Schwamb | scaling and density table |
| G5-D Surface texture | T21 | Paragas library, Zieba and Kreidberg points; age and flux | LHS 3844 b analogue |

### 2.7 Where a verified source overrides a partly verified one

| Item | Earlier value | Winning verified value | Source |
|---|---|---|---|
| k1 to k3 | 0.5, 1/3, 0.1 (SNIPPET rounding) | 13/25, 9/25, 3/25 with intervals | HP22 Table 2, READ |
| k4 width | "by eye" 0.5 dex | 0.51 dex (SE 0.08), N=21 | DATA |
| k5 | median SNIPPET, 0.6 dex | 2.5 % (0.16 dex), 0.65 dex (0.12), N=16 | DATA + READ |
| k7 | 67P-centred | population geo-means (CH4, C2H6, H2CO, H2S, OCS) | DATA, READ |
| snow-line scatter 0.2 dex | RECALLED | READ, but a model consequence | arXiv:1505.03516 |
| T19 spread 0.58 to 0.82 | abstract only | 0.74 to 0.82 (dwarf planets), 0.3 to 0.45 (small) | COMPUTED from READ |
| ISO range | CO/H2O to 0.23 | CO/H2O to 1.66, CO2/H2O to 7.6 | READ |
| Sasaki vs Cilibrasi (T18) | both READ | no override; spread reported | READ |

## 3. Sources verified (id, N, access, licence)

All arXiv ids above were downloaded from export.arxiv.org and read (PDF text) unless marked abstract-only; none states a licence in the text read. DATA: Robinson database, doi 10.7488/ds/7723 (datashare.ed.ac.uk API; CC BY 4.0 for the compilation; 909 rows, 96 designations). Nature pages for Hoy 2026, Putirka and Xu 2021 and Canup and Ward 2006 were read through doi.org (abstracts, open metadata). Crossref metadata for Dello Russo 2016, Mumma and Charnley 2011, Bockelee-Morvan and Biver 2017, A'Hearn 1995, Crida and Charnoz 2012, Agnor and Hamilton 2006. Downloaded but not analysed: arXiv:2510.26308, 2512.07318 (Table 4 values via 2603.25002), 1112.2764, 1204.1544. Refused or captcha: PMC (Bierson and Nimmo), publisher PDFs.

## 4. Boss must supply

1. Dello Russo 2016 (Icarus 278, 301): taxonomy tiers and per-comet table; Mumma and Charnley 2011 and Bockelee-Morvan and Biver 2017 range tables.
2. Tabulated 3I/ATLAS abundances at fixed r_h, and any further interstellar comet with a coma, to constrain the system offset.
3. A population synthesis or survey of moon origin classes versus planet mass and orbit (none found).
4. Bierson and Nimmo 2019 full text and table of 17 KBO densities and radii.
5. Meteoritical Bulletin fall counts for CI, CM, CO, CV, enstatite and achondrites, plus a delivery model.
6. A published relation between surface age, impact flux and regolith texture for airless planets; Zenodo 17986129 licence.

## 5. Recalled / to verify

A'Hearn 2012 (ApJ 758, 29) title and venue; DeMeo 2019 A-type paper (Icarus 322); Kiefer 2014 Nature id (Ca II exocomets, only SNIPPET; Zieba 2019, arXiv:1903.11071, is TESS photometry with no composition); Jacobson and Morbidelli, Citron, Chambers on impact-moon frequency; Batygin and Morbidelli 2020 id; Halley CHON organic fraction; 67P organic share (Bardyn 2017, via snippets); H, L, LL fall shares original; Noll 2008 KBO binary fraction (relevant to a binary moon share); Eistrup 2019, Cochran 2015, Rubin 2019 (G1 lane READ, not re-read).
