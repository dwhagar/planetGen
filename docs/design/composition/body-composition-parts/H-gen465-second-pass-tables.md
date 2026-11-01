# Body composition part H: the GEN.465 second extraction pass and the water-share result (GEN.468)

Read-only research by Core 4 (about 7:03 AM EDT 2026-10-23, lane clock), merged here by the TODO thread. The lane's findings file and the extracted files are in `/mnt/project-files/design-drafts/gen465-second-pass/` (`README-gen465-gen468-findings.md` first). Only derived values are kept in the repo (the raw per-star file of 1,615 stars and the 37 KB text extract of [R2726] stay in the draft folder; they are source data, not repo content).

Labels: READ = number taken from the source text or data file; COMPUTED = arithmetic on READ numbers; 2ND-HAND = quoted by a source that was read, not the original; NOT REACHED = the host refused or served only a preview; RECALLED = from memory, to be verified.

## 1. Source status

| Source | Status |
|---|---|
| UNSCEAR 2008 Annex B (radon, soil Ra-226) | READ (full PDF) |
| VizieR J/ApJS/225/32 ([R0535], 1,615 stars, 15 elements) | READ |
| Zenodo 10527170 ([R4529], Bern model, 31,868 planets, nominal run) | READ |
| [R2726] (the rank 1 paper of the Boss list) | READ in full (28 pages, through the ADS link gateway; the IOP page is flaky and sometimes returns a captcha) |
| [R2729] (the 11-group chondrite table, K, Th, U) | NOT REACHED: Springer serves the abstract and a preview. A person has to fetch [R2729] |
| [R4658: Table 1] (metal volume percent per group) | READ; replaces the RECALLED metal list for the volume percent |
| Science, Royal Society, ScienceDirect, AGU, Oxford Academic, Springer | still refuse; [R0067, R0638, R1445, R3417] and the rest of the part files' lists stay NOT REACHED |
| [R2144] | opens as a 15-page scan with no text layer; needs reading by eye or OCR |

## 2. GEN.468: the water share beyond the snow line

[R2726: Table 11] (condensate distribution, mass percent of the total nebular mass, equilibrium condensation, noble gases excluded, solar-system composition), READ:

| Condensate | % of total mass | share of all condensate (COMPUTED) |
|---|---|---|
| Silicates and oxides | 0.303 | 0.204 |
| Metal + FeS | 0.186 | 0.125 |
| Total rock | 0.489 | 0.329 |
| H2O ice | 0.571 | 0.384 |
| CH4 ice | 0.330 | 0.222 |
| NH3 ice | 0.097 | 0.065 |
| Total ices | 0.998 | 0.671 |
| Total condensate | 1.487 | 1.000 |

Printed: H2O ice/rock 1.17; oxygen in rock 22.75 percent before water ice condenses at 182 K (10^-4 bar). The photospheric column gives the same ratios. The [R0100] column gives water/rock 2.09, and the text says the revision brought it to about unity.

Water as a share of the solids, by which ices have condensed (COMPUTED from the table; temperatures READ from Lodders Tables 7, 8 and 10):

| Where | Solids counted | Water share |
|---|---|---|
| Just outside the water line (182 K to 131 K) | rock + water | 0.539 |
| Past the NH3 hydrate line (131 K) | rock + water + NH3 | 0.494 |
| Past the CH4 line (41 K equilibrium; CH4 clathrate at 78 K takes about 9 percent of the carbon) | everything | 0.384 |

The Bern model ([R1268, R1269], NGPPS I, [R1268], sections 3.3.3 and 4, READ) uses the [R2865, R4335] condensation model, "H2O makes up about 60 % of all ices in mass", and maps all ice species to water ice for the planet structure. In the Zenodo run the water mass fraction of planets with a steam envelope tops out at 0.51 to 0.53; that cap is the model's total ice share of solids, labelled water. With a 60 percent water share of the ices it is about 0.30 of the solids in water proper (COMPUTED), inside the Marboeuf-derived 0.14 to 0.39.

Reading (the numbers are READ, the reading is the lane's):
1. Different denominators. "About 0.5" is the water share of rock plus water (or of rock plus all ices lumped as water). The Marboeuf-derived 0.14 to 0.39 is the water share of all solids. Same chemistry, two bookkeeping choices. Lodders gives 0.54, 0.49 and 0.38 for the three rows above.
2. Different total ice to rock. Lodders (equilibrium, 2003 abundances) has total ice/rock 0.998 / 0.489 = 2.04 with water about 57 percent of the ices; Marboeuf's comet-calibrated value is 1 +/- 0.5 with water about 60 percent of the ices (2ND-HAND via NGPPS I). So water/rock is 1.17 against 0.6 (0.3 to 0.9). A factor of about two is real and no source read settles it: it is equilibrium nebular chemistry with all carbon as CH4 against a composition read off comets (CO and CO2 kept in the gas or lost).

Result for GEN.468: keep the default as the Marboeuf-derived sequence (Boss, 1:21 AM EDT 2026-10-17, "keep the switch") and name the denominator in the switch: "water over all solids" 0.14 to 0.39 (default), "water over rock plus water" about 0.5 (0.49 to 0.54 from [R2726]). The physics-supported range for the water share of all solids beyond the line is 0.38 (equilibrium, all ices) down to 0.14 (comet-calibrated); the upper ends agree at 0.38 to 0.39. Open: whether planetesimals follow nebular equilibrium (ice/rock about 2.0) or comets (about 1.0). The snow-line temperature 182 K at 10^-4 bar is READ; Boss's decision keeps the 170 K centre with 0.2 dex scatter.

## 3. Extracted tables

### 3.1 Zenodo 10527170: water fraction of the planets, by final semi-major axis (READ, COMPUTED)

Nominal run, 5 Gyr, `nominal/ref_evo5e9.dat`, 31,868 rows. Water fraction W = Menv x Zhomo / total mass (water is carried in the envelope; the authors' `fice` column is 0 in this run). No formation distance is in the file, so W is tabulated against the final a; planets that migrated from outside the snow line sit inside it. Planets of 1 to 30 Earth masses:

| a (AU) | planets | with water envelope | share | W p05 | W p25 | W p50 | W p75 | W p95 |
|---|---|---|---|---|---|---|---|---|
| 0.05 to 0.2 | 1036 | 227 | 0.219 | 0.039 | 0.196 | 0.336 | 0.426 | 0.483 |
| 0.2 to 0.5 | 789 | 157 | 0.199 | 0.014 | 0.086 | 0.211 | 0.431 | 0.499 |
| 0.5 to 1 | 834 | 331 | 0.397 | 0.008 | 0.174 | 0.298 | 0.478 | 0.503 |
| 1 to 2 | 952 | 413 | 0.434 | 0.188 | 0.298 | 0.401 | 0.479 | 0.503 |
| 2 to 5 | 1159 | 403 | 0.348 | 0.285 | 0.395 | 0.461 | 0.491 | 0.504 |
| 5 to 10 | 898 | 454 | 0.506 | 0.366 | 0.494 | 0.504 | 0.506 | 0.511 |
| 10 to 30 | 618 | 416 | 0.673 | 0.475 | 0.505 | 0.509 | 0.514 | 0.529 |
| 30 to 100 | 36 | 12 | 0.333 | 0.509 | 0.523 | 0.544 | 0.548 | 0.564 |

By total mass (Earth masses), same population:

| M (Earth) | n | W p05 | W p25 | W p50 | W p75 | W p95 |
|---|---|---|---|---|---|---|
| 1 to 3 | 1457 | 0.065 | 0.443 | 0.504 | 0.508 | 0.521 |
| 3 to 6 | 727 | 0.069 | 0.286 | 0.441 | 0.499 | 0.509 |
| 6 to 10 | 407 | 0.083 | 0.275 | 0.390 | 0.468 | 0.504 |
| 10 to 20 | 247 | 0.078 | 0.219 | 0.346 | 0.441 | 0.498 |
| 20 to 30 | 43 | 0.125 | 0.182 | 0.266 | 0.390 | 0.488 |

Maximum over all rows 0.568; 7.2 percent of rows have no envelope; the median total mass of a row is 0.12 Earth masses (embryos are included). The cap of about 0.5 is a model input, not a result. This answers part C's request for planet-by-planet ice fractions for the migration population.

### 3.2 [R0535] (VizieR J/ApJS/225/32): 1,615 FGK stars, [X/H] against the Sun (READ, COMPUTED)

Statistics after dropping |[X/H]| >= 1 (one catalogue sentinel, [Na/H] = -21.41, and some outliers):

| Element | n kept | mean | std | median | robust sigma |
|---|---|---|---|---|---|
| C | 1592 | -0.025 | 0.180 | -0.02 | 0.178 |
| N | 1590 | 0.050 | 0.257 | 0.07 | 0.267 |
| O | 1592 | 0.079 | 0.167 | 0.08 | 0.178 |
| Na | 1611 | 0.010 | 0.238 | -0.01 | 0.222 |
| Mg | 1614 | 0.004 | 0.160 | 0.00 | 0.163 |
| Al | 1587 | 0.011 | 0.202 | 0.01 | 0.208 |
| Si | 1613 | 0.004 | 0.181 | 0.01 | 0.193 |
| Ca | 1613 | 0.053 | 0.183 | 0.07 | 0.178 |
| Ti | 1615 | 0.052 | 0.158 | 0.06 | 0.148 |
| V | 1592 | 0.008 | 0.173 | 0.02 | 0.178 |
| Cr | 1613 | 0.027 | 0.198 | 0.05 | 0.193 |
| Mn | 1610 | -0.016 | 0.274 | 0.00 | 0.282 |
| Fe | 1614 | 0.038 | 0.198 | 0.06 | 0.193 |
| Ni | 1614 | 0.016 | 0.211 | 0.01 | 0.215 |
| Y | 1609 | 0.064 | 0.221 | 0.09 | 0.193 |

- Joint [Mg/H], [Si/H], [Fe/H], n = 1613: correlations Mg-Si 0.936, Mg-Fe 0.924, Si-Fe 0.873; covariance (Mg, Si, Fe): 0.0257, 0.0271, 0.0293 / 0.0326, 0.0312 / 0.0390.
- Regression on [Fe/H]: [Mg/H] = 0.750 [Fe/H] - 0.024 (scatter 0.061 dex); [Si/H] = 0.799 [Fe/H] - 0.027 (scatter 0.088 dex).
- [Mg/Fe] mean -0.034, sigma 0.079; [Si/Fe] mean -0.034, sigma 0.097.
- Number ratios against the [R0169] solar values (Mg 7.60, Si 7.51, Fe 7.50; RECALLED reference): Mg/Si 5, 50, 95 percent = 1.00, 1.20, 1.62; Fe/Si 0.79, 1.02, 1.59.
- Ca and Al: [Ca/H] mean 0.052, sigma 0.188; [Al/H] mean 0.009, sigma 0.215; [Ca/Fe] mean 0.015, sigma 0.071; [Al/Fe] mean -0.032, sigma 0.137. This closes part B item 1 for the stellar side; the planet-minus-star scatter is still open.
- C/O number ratio (Asplund C 8.43, O 8.69): median 0.457, 5 and 95 percent 0.26 and 0.59.

### 3.3 UNSCEAR 2008 Annex B (READ)

- World soil averages (paragraph 76): 238U 33, 226Ra 32, 232Th 45, 40K 412 Bq/kg (the previous 40K value was 420); reported maxima 1000 (U), 360 (Th), 3200 (40K) Bq/kg. This replaces the RECALLED 33, 45, 420 of part E.
- Outdoor absorbed dose rate in air, worldwide 58 nGy/h (distribution mean 50 to 59).
- Indoor radon (paragraph 97): country averages 9 to 184 Bq/m3; log-normal with geometric mean 45 Bq/m3 (previous 40), geometric standard deviation 2.1; per caput dose 1.15 mSv/yr (radon) and 0.1 mSv/yr (thoron).
- Table 12 totals: external terrestrial 0.48 mSv/yr (outdoors 0.07, indoors 0.41; range 0.3 to 1.0); inhalation 1.26 (range 0.2 to 10); ingestion 0.29; cosmic 0.39; total 2.4 (range 1 to 13). This confirms the repo's `GROUND_MSV_YR` source of 0.48 mSv/yr (part E verify item 4).
- Not found: a radon emanation fraction from Ra-226 in soil, or soil Ra-226 ranges per country (those are in the annex tables of the 2000 report, not this one).

### 3.4 [R2726] (READ)

- 50 percent condensation temperatures at 10^-4 bar, solar-system composition: K 1006 K (feldspar), Na 958 K, Mg 1336 K, Si 1310 K, Th 1659 K (hibonite), U 1610 K (hibonite), S 704 K (FeS, troilite), Fe 1357 K (alloy), H2O ice 182 K, NH3 hydrate 131 K (50 percent at 123 K), CH4 hydrate 78 K (50 percent at 40 K). Photospheric composition: K 1001, Th 1647, U 1604, Fe 1351, H2O 181. This turns part E's RECALLED K condensation of 1006 K into READ.
- CI chondrite group mean (Table 3, ppm by mass): K 530 +/- 24, Th 0.0309, U 0.0084 (the printed uncertainties for Th and U are garbled in the text layer; check the page); Fe 182,150 +/- 1,910; O 458,200 +/- 5,750.
- Oxygen: about 23 percent of the O is bound to rock before water ice (26 percent after magnetite at 371 K).

### 3.5 Metal per chondrite group ([R4658: Table 1], READ, volume percent)

CI 0, CM 0.1, CO 1 to 5, CV 0 to 5, CK 0 to 5, CR 5 to 8, CH 20, CB 60 to 80, H 8, L 4, LL 2, EH 10, EL 10, R 0.1, K 7. COMPUTED weight percent for the ordinary groups, with metal 7.9 g/cm3 and the rest 3.3 g/cm3 (assumed grain densities, so approximate): H 17, L 9, LL 5, EH and EL 21. Part D's RECALLED values (H about 17, L about 8, LL about 3) hold for H and are within 2 points for L and LL. Not reached: the weight percent measured per group ([R2729], Rubin).

## 4. What stays Boss's

1. The [R2729] PDF (Springer, preview only): the 11-group chondrite table with K, Th and U per group.
2. Whether planetesimal ice/rock follows nebular equilibrium (about 2.0) or comets (about 1.0): the default of the water-share switch (GEN.468; no source read settles it).
3. The hosts that still refuse (Science, Royal Society, ScienceDirect, AGU, Oxford Academic).
