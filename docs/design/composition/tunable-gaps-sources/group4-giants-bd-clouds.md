# Group 4: giant planets, brown dwarfs, hot Jupiters, clouds and hazes

Rows: part 09 f5, f6, f7, g1, g4, h2, h3, h4, n6, n7; body table T12. Labels follow `_worker-rules.md`. Nothing in the repo was edited; tools run: curl, pdftotext, python3 -I, and a small Gibbs-energy script in the scratchpad.

Usability rule (Boss, 3:02 AM EDT): where two verified sources disagree, the later, larger, better-controlled one is used and both are named; where a verified source contradicts a partly verified one (SNIPPET, RECALLED, 2ND-HAND), the verified source wins. Every such place is listed in section 2.9. Solar System-only constraints are flagged per row (f = 0.2 in f5, T12 structure, the Titan lab rate in n7, the solar C/O 0.55 anchor in f7).

## 1. Summary and verdict table

The public record is richest for stellar C/O, the bulk heavy-element mass of giants and the brown-dwarf Teff-to-type map; it is thin for clouds, hazes and Kzz. Two rows are not probabilities: g1/g4 (type is a function of Teff with scatter, so draw mass and age, evolve, then map) and h2/h3 (regimes are outputs of equilibrium chemistry, computable with MIT-licensed code and NASA-format polynomials). A public data set found in this pass changes f5 to f7: LExACoM (Lothringer, Lowson and Fu 2025, arXiv 2510.26785; CSV public; 69 planets with retrieved C/O and metallicity, updated to 2026-08-28) has no Solar System member, so it can tune m and C/O within the no-fixed-content rule.

| Row | What it is | Verdict | One-line reason |
|---|---|---|---|
| f5, f6 | giant metallicity m | (a) for bulk Z; (b) for atmospheric m | Bulk: N=47, 147 and 403 fits with published scatter. Atmosphere: N=69 retrievals, scatter 0.4 to 0.9 dex, far wider than the design 0.3 dex; the bulk-to-atmosphere factor f has only Jupiter and Saturn behind it. |
| f7 | giant C/O core plus tail | (a) stars, (b) planets | Stellar C/O: N=849, 499, 355 agree (median 0.45 to 0.50, no star above 0.8). Planet C/O: N=69, median 0.57, 6 of 69 above 0.8, method-dependent. The design tail bracket [0.0013, 0.34] is wider than the data allow. |
| g1, g4 | brown-dwarf type vs Teff | not a probability; derive | Teff to spectral type is a fitted curve with RMS 79 to 134 K (Kirkpatrick 2021). The "T/Y gap 500 to 700 K" is late T (T8, T9). Draw mass from a power law (alpha 0.6 +-0.1) and age, cool with evolution tracks. |
| h2, h3 | hot-Jupiter regimes at 1500 and 2200 K | not a probability; derive by physics | Equilibrium chemistry computed here: CH4=CO at 966 to 1438 K (0.1 to 10 bar), H2O 10 percent depleted at 2331 to 3095 K (1 mbar to 1 bar). Edges are pressure-dependent outputs, not constants. |
| h4 | TiO/VO retention, Kzz vs grain size | (c) for Kzz prior; (b) for a Teq threshold | No population measures Kzz (WASP-107 b gives 10^11.6 in one paper and 10^8.4 to 10^9.0 in another). Ti-bearing species are retained in 1 of 6 ultra-hot Jupiters in one survey, at the hottest Teq. |
| n6 | silicate vs haze crossover near 950 K | (b) | A model output (Gao 2020) tuned to the N=34 WFC3 feature compilation; no measured width. Estrela 2022 (N=62) gives a different quantity, the cloud minimum at 1461 K. |
| n7 | haze yield LU[0.01, 1] | (c) | Lab rates span 250x but only for 100 to 10000x solar sub-Neptune gases at 300 to 600 K; yield is degenerate with metallicity, Kzz and cloud cover in every population statistic. |
| T12 | diluted vs compact core | (c) | Exoplanet bulk Z cannot tell the two apart (READ in two papers); the only structural evidence is Juno/Cassini gravity of Jupiter and Saturn, a Solar System-only constraint. |

Groups for builders: (A) f5, f6, f7, T12 share the T16/C25/TF19 bulk-Z fit and the LExACoM table; (B) g1, g4 share Kirkpatrick 2021 Tables 13 and 14 and the mass function; (C) h2, h3 share one equilibrium-chemistry module; (D) h4, n6, n7 share the WFC3/high-resolution feature compilations and are best built as one "aerosol amplitude vs Teq" test.

## 2. Per-row detail

### 2.1 f5 and f6: giant metallicity

**Question 1: is 0.26 dex a tunable width, and what does the retrieval trend add?**

1. Ensemble.

| Data set | N | What was measured | Status |
|---|---|---|---|
| Thorngren 2016, arXiv 1511.07854 | 47 giants, 20 Earth masses to 20 MJ, no inflation (F* < 2e8 erg/s/cm2) | Mz = (57.9 +-7.03) M^(0.61+-0.08) Earth masses, intrinsic spread factor 10^sigma = 1.82 +-0.09; Zpl/Zstar = (9.7 +-1.28) M^(-0.45+-0.09), Kendall tau 0.44, p 1.3e-5; no scatter printed for this second fit | READ |
| Chachan 2025, arXiv 2509.20428 | 147 warm giants (< 1000 K, age > 0.5 Gyr) | MZ = Mcore + fZ (Mp - Mcore), Mcore 14.7 (+1.8, -1.6) Earth masses, fZ 0.09 +-0.01, scatter 0.66 +-0.08 times MZ; Zp/Z* = a/Mp + b, a 2.6 +-0.3, b 3.3 +-0.5 (MJ); normalising by Z* does not cut the scatter (0.67) | READ |
| Thorngren and Fortney 2019, arXiv 1811.11859 | 403 giants | Zp and atmospheric upper limits; Jupiter and Saturn atmospheres about 20 percent of the bulk limit; table "available for download" | READ |
| Welbanks 2019, arXiv 1912.04904; Sun 2024, arXiv 2402.08292 | 19 | H2O in 14, Na 6, K 6; H2O rises with falling mass but sits below the Solar System carbon trend; Sun: linear trend not distinguishable from no trend, 45 to 58 planets needed | READ |
| LExACoM, Lothringer 2025, arXiv 2510.26785 | 66 (paper), 69 unique planets (CSV 2026-08-28) | retrieved [M/H], C/O, Teq, mass; CSV public at jlothringer.github.io (table licence not stated), code MIT | READ, DATA |

2. Physical bracket. Atmospheric Z cannot exceed bulk Zp, and Zp is at most 1 (TF19, READ). The design range [-0.5, 2.0] dex is exceeded by LExACoM ([M/H] from -1.99 to +2.63, COMPUTED, mixed C/H, O/H, (O+C)/H definitions).

3. Verdict. The 0.26 dex is the intrinsic scatter of Mz at fixed mass in the Mz-M fit (COMPUTED: log10 1.82 = 0.260 dex, band 0.25 to 0.27), not of the Zpl/Zstar fit that gives f5 its centre. Chachan says the scatter of Zp/Z* is no smaller than that of Zp (0.67 vs 0.66 multiplicative, 0.22 to 0.29 dex under a lognormal or small-sigma reading, COMPUTED), so 0.26 dex is a sound proxy for the bulk scatter. But the generator needs the atmosphere. Retrievals add two things: the atmospheric trend is steeper (ExoComp text: -1 to -2 dex per dex of mass; Solar System shallowest at -1.13 +-0.09; the O/H fits print -1.80 +-0.36 transit and -1.81 +-0.16 eclipse, column assignment inferred from a garbled table) against the bulk exponent -0.45 +-0.09; and the scatter about the trend is large. From the CSV (COMPUTED, my linear fit of [M/H] on log10 M, one row per planet): all 69 slope -0.32, residual 0.80 dex; transit 27 slope -0.64, 0.94 dex; eclipse 18 slope -1.30, 0.74 dex; direct imaging 24 slope -0.22, 0.41 dex; mean quoted error 0.2 to 0.5 dex, naive intrinsic scatter 0.36 to 0.85 dex. The design width 0.3 dex fits the bulk scatter but is 1.5 to 3 times narrower than the retrieved atmospheric scatter.

   Verdict: (a) if f5 is redefined as bulk Z (centre C25, width 0.22 to 0.29 dex); (b) for atmospheric m, with f (0.2, Jupiter and Saturn only) replaced by a distribution tuned to LExACoM.

   Design note corrected: part 09 f5 says "Chachan 2025 gives lower Z for sub-Saturns". The paper says sub-Saturns (< 0.3 MJ) have high bulk Zp, 0.2 or more, not a selection effect. C25 agrees with T16 at 0.1 MJ (Zp/Z* 29.3 vs 27.3, COMPUTED) and is 1.4 to 1.6 times lower from 0.3 to 3 MJ (12.0 vs 16.7 at 0.3; 5.9 vs 9.7 at 1; 4.2 vs 5.9 at 3), rejoining at 10 MJ (3.56 vs 3.4). C25 supersedes T16 for the centre above 0.3 MJ.

4. Knobs.

| Knob | Target | Value (uncertainty) | Source, status |
|---|---|---|---|
| bulk centre | Zp/Z* = 2.6/Mp + 3.3; MZ = 14.7 + 0.09 (Mp - 14.7) | a +-0.3, b +-0.5, Mcore +1.8/-1.6, fZ +-0.01; 0.03 to 12 MJ | C25 READ |
| bulk width | scatter of Mz | 0.260 dex (0.25 to 0.27) or 0.66 +-0.08 multiplicative | T16, C25 READ |
| factor f | atmosphere/bulk | 0.2 (Jupiter, Saturn only); upper bound 1 | TF19 READ; Solar System only |
| atmospheric slope | d[M/H]/d log M | -0.3 to -1.3 by geometry | COMPUTED |
| atmospheric width | residual sd | 0.4 to 0.9 dex | COMPUTED |

5. Consistency test. Compare simulated residual sd and slope with LExACoM (tolerance 0.2 dex on width, 0.3 on slope; the paper's fit errors are +-0.16 to 0.4); no completeness model exists for retrievals, so cut on mass and Teq only. For the bulk, compare median Zp/Z* in mass bins with C25 (29.3, 5.9, 3.6 at 0.1, 1, 10 MJ).

6. Flags. f = 0.2 and any use of Jupiter or Saturn atmospheric C/H are Solar System-only. The bulk fits and LExACoM contain no Solar System body in the comparison sample (LExACoM fits Jupiter to Neptune separately; ignore).

7. Group. Builder A.

### 2.2 f7: giant C/O, core and tail

**Question 2.**

1. Ensemble (stars).

| Data set | N | Result | Status |
|---|---|---|---|
| Brewer and Fischer 2016, arXiv 1608.06286 | 849 FGK dwarfs, SNR > 100 | median C/O 0.47, sd 0.07, long low tail, maximum 0.66 +-0.068, none carbon-rich; "frequency of carbon-rich dwarfs < 0.13 percent"; quadratic in [Fe/H], rms 0.06; Table 1 electronic; planet hosts (156) follow the same distribution | READ |
| Suarez-Andres 2018, arXiv 1801.09474 | 499 HARPS (99 hosts) | 100 percent of hosts below 0.8; mean 0.46 +-0.11 (low-mass hosts), 0.50 +-0.10 (high-mass), 0.45 +-0.13 (singles); median error 0.14; solar reference 0.61 | READ |
| Delgado Mena 2010, arXiv 1009.5224 | 61 HARPS plus 39 other hosts, 270 comparison | "34 percent of stars with known planets have C/O greater than 0.8"; oxygen from the single [O I] 6300 line | READ |
| Hypatia via ExoComp Table 5 | 355 FGKM | median 0.468, sd 0.143 (columns inferred) | 2ND-HAND |
| Costa-Almeida 2026, arXiv 2606.30722 | 290 HARPS hosts | stellar C/O 0.17 to 0.95, giant hosts median 0.43 (+0.02, -0.03), solar 0.53 | READ |
| Nissen 2014, arXiv 1406.5218 (candidate "Nissen 2013" is a proceeding, 1306.4836, SNIPPET); Teske 2013, arXiv 1309.6032 (candidate "Teske 2014" is arXiv 1403.6891, not read) | F and G dwarfs; 55 Cnc | thin-disc stars have higher C/O than thick-disc and halo at equal [O/H]; the [O I] 6300 line is blended with Ni I in metal-rich stars (55 Cnc 1.12 falls to 0.78 +-0.08 with other lines) | READ (abstracts) |

The 34 percent vs sub-percent conflict. Both numbers are READ but measure different things with different method quality. DM2010 used one oxygen line that Teske 2013 shows is blended at high metallicity; B&F use hundreds of lines and drop it; SA18, a larger HARPS sample, finds none above 0.8. B&F's "< 0.13 percent" is the carbon-rich fraction; their READ maximum of 0.66 gives zero stars above 0.8 in 849, so an upper limit of 3/849 = 0.35 percent at 95 percent (COMPUTED). DM2010's 34 percent is superseded on method grounds, and the design range [0.0013, 0.34] for f_tail takes its two ends from this disagreement; the upper end should be dropped as a stellar number.

Ensemble (planets). LExACoM, 69 unique planets (DATA, COMPUTED): C/O median 0.57, sd 0.27; transit 27 planets median 0.42 sd 0.38, eclipse 18 median 0.69 sd 0.18, direct imaging 24 median 0.58 sd 0.12. ExoComp text (66 planets): medians 0.350 (transit), 0.59 (direct), 0.7 (eclipse), sd about 0.2, K-S p 0.011 and 0.012 for transit against direct and eclipse; none above 1 at 1 sigma except TOI-5205 b (1.3 +-0.4). Design bins (<0.35, 0.35-0.7, 0.7-0.8, 0.8-1.0, >=1.0) hold 12, 43, 8, 3, 3 of 69 planets (0.174, 0.623, 0.116, 0.043, 0.043), against the design example 0.106, 0.747, 0.117, 0.009, 0.021.

Planet-minus-star offset (2ND-HAND via Costa-Almeida 2026): for six planets shared with Changeat 2022 (HST and Spitzer free retrievals) every planet exceeds its star, median +0.43 +-0.13; for directly imaged planets (Baburaj 2025) median -0.04 +-0.04. WASP-77A b alone has 0.59 +-0.08 (Line 2021, arXiv 2110.14821, READ), 0.36 +0.10/-0.09 (August 2023, 2ND-HAND) and a host star at 0.44 +0.07/-0.08, an offset of -0.08 to +0.15 for one planet. The offset has no single distribution (0.2 to 0.5 by geometry). The design core centre 0.55 rests on one planet (Line, N=1).

2. Physical bracket. C/O is bounded by 0 and about 1 for accretion including oxygen-bearing solids; 3 of 69 planets have a value at or above 1, all with errors reaching below 1.

3. Verdict. (a) for the stellar C/O distribution (three independent analyses agree on no tail above 0.8). (b) for planets: centre, width and tail can be tuned to LExACoM only to its pooled statistics, with a 0.2 to 0.3 geometry systematic in the median. The stellar catalogue does not fix the planet tail: the planet fraction above 0.8 is 6/69 = 0.087 (Wilson 95 percent 0.040 to 0.177), or 3/69 = 0.043 (0.015 to 0.120) if the 1-sigma lower bound must exceed 0.8 (COMPUTED). Both intervals exclude the design ends 0.0013 and 0.34. The C/O above 0.8 is carried by hot and ultra-hot planets (WASP-121 b, WASP-19 b, HD 149026 b, HIP 67522 b, HR 8799 e), so it mixes chemistry, retrieval systematics and formation.

4. Knobs and targets.

| Knob | Target | Value (uncertainty) | Status |
|---|---|---|---|
| stellar centre, width | median, sd | 0.47 and 0.07 (B&F); 0.45 to 0.50 and 0.10 to 0.13 (SA18); 0.468 and 0.143 (Hypatia) | READ, 2ND-HAND |
| stellar tail above 0.8 | fraction | < 0.35 percent (95 percent) | COMPUTED |
| planet minus star | median by geometry | +0.43 +-0.13 (6 planets, HST), -0.04 +-0.04 (direct) | 2ND-HAND |
| planet core centre | median of 69 | 0.57 (0.42 transit, 0.69 eclipse, 0.58 direct) | COMPUTED |
| planet core width | sd below 0.8 | about 0.2 | READ (ExoComp text) |
| f_tail (planets) | fraction above 0.8 | 0.043 to 0.087, 95 percent 0.015 to 0.177 | COMPUTED |

5. Consistency test. Compare simulated fractions in the five bins and the median with the counts above (binomial, N=69, tolerance from the Wilson intervals); on the stellar side compare host C/O with B&F (median 0.47, sd 0.07, none above 0.66). Hot-Jupiter hosts are metal rich and C/O rises with [Fe/H] (rms 0.06), so sample hosts by [Fe/H].

6. Flags. The core centre 0.55 and the solar 0.589 are Solar System anchors; LExACoM and the stellar catalogues replace them.

7. Group. Builder A with f5 (one table, one loader).

### 2.3 g1 and g4: brown-dwarf type vs Teff

**Question 3.** Brown-dwarf spectral type is a deterministic function of Teff with scatter, not a weighted draw; the mass and age are what should be drawn.

1. Ensemble.

| Data set | Source | N | Result | Status |
|---|---|---|---|---|
| Kirkpatrick et al. 2021, arXiv 2011.11616 (ApJS 253, 7) | READ (PDF) | 525 L, T, Y dwarfs within 20 pc; Table 14 lists dwarfs with Teff | piecewise Teff(SpT) fits in Table 13 with RMS; complete to 20 pc except types >= T8 and Teff < 600 K; mass function dN/dM ~ M^(-alpha), alpha = 0.6 +-0.1; simulation range 125 to 2800 K, 0.01 to 0.10 solar masses, constant birthrate over 10 Gyr | READ |
| Filippazzo et al. 2015, arXiv 1508.01767 | READ (abstract, sections 8 and 9) | 145 field-age and 53 young M6 to T9 dwarfs | Teff typical uncertainty 6 percent; Teff vs SpT dispersion minimum 200 K at late T; young low-gravity L dwarfs up to 300 K cooler than field dwarfs of the same type | READ |
| astro-ph/0402475 | READ (abstract) | 42 dwarfs | Teff falls steeply for M6 to L7 and T4 to T9 and is about constant, 1450 K, for L7 to T4 | READ |
| Beiler et al. 2024, arXiv 2407.08518 (the 2023 paper, arXiv 2306.11807, is one Y0 dwarf at 467 K) | READ (abstract, figure labels) | 23 late-T and Y dwarfs with JWST | Teff 1000 to 350 K, median error 20 K; per type: T6 897 and 987 K, T6.5 1044 K, T8 525 to 805 K, T9 565 to 585 K, Y0 371 to 510 K, Y0.5 406 and 412 K, Y1 362 to 526 K | READ |

2. Physical bracket. Teff of a brown dwarf is fixed by mass and age through cooling tracks (Saumon and Marley 2008, Baraffe 2003, as used by Kirkpatrick 2021, READ for their use); the hydrogen-burning limit near 0.075 to 0.08 solar masses, the deuterium limit near 13 MJ.

3. Verdict. Not a probability: derive. The empirical map (Kirkpatrick 2021 Table 13, READ via the layout text of the page; coefficients and RMS as printed):

| SpT range | Teff polynomial (SpT: L0 = 0, T0 = 10, Y0 = 20) | RMS |
|---|---|---|
| 0 to 8.75 (L0 to L8.75) | 2237.5 - 144.96 s + 4.0301 s^2 | 134 K |
| 8.75 to 14.75 (L8.75 to T4.75) | 1437.9 - 18.309 s | 79 K |
| 14.75 to 22 (T4.75 to Y2) | 5141.3 - 368.65 s + 6.7301 s^2 | 79 K |
| Teff from absolute H (9.5 to 25 mag) | quartic in MH | 88.1 K |

Evaluated (COMPUTED): L0 2238 K, L5 1613, L8 1336, T0 1255, T4 1182, T6 966, T8 686, T9 567, Y0 460, Y1 368, Y2 288 K. These agree with the individual JWST Teff above (Y0 at 468 to 510 K, T8 at 525 to 805 K), so the scatter at fixed type is about 80 to 130 K, with a plateau between L8.75 and T4.75 where the slope is only -18 K per subtype, so Teff cannot fix the type to better than about 4 subtypes there (79 K / 18 K).

Consequences for the rows. g2 (L, 1300 to 2200 K), g3 (T, 700 to 1300 K) and g5 (Y, below 500 K) are the same map and need no weights. g4, the "T/Y transition 500 to 700 K", is not a type: in the empirical map 500 to 700 K is T8 to T9 (686 K at T8, 567 K at T9), and the measured Y0 dwarfs start at 510 K and below. The gap came from rounding the T and Y rows; use the continuous curve above, with Y0 near 460 K (RMS 79 K; observed Y0 range 371 to 510 K). g1 (late-M, 2200 to 2800 K) is the extension of the map above L0 (2238 K): the READ table sample gives M7.5 at 2871 K and M9.5 at 2390 K (Kirkpatrick Table 14 excerpt), and the model grid stops at 2800 K, so the 2200 and 2800 K edges fit but the M-type fit is not in Table 13 (NOT REACHED). Low-gravity young L dwarfs are up to 300 K cooler than field dwarfs of the same type (Filippazzo, READ), so the map needs a gravity term below 500 Myr.

If a mass draw is needed: dN/dM ~ M^(-0.6 +-0.1), constant birthrate over 10 Gyr, masses 0.01 to 0.10 solar masses (READ for the fitted model and range; alpha applies to the field, not to young clusters). The 2019 paper (candidate, ApJS 240, 19; abstract says a power law with alpha about 0.6 and a cutoff no higher than about 5 MJ) is SNIPPET and superseded by the 2021 one.

4. Knobs: alpha 0.6 +-0.1 (target: space density per 150 K Teff bin, Kirkpatrick Table 15, numbers NOT REACHED); the Teff map with the RMS as Gaussian scatter at fixed type; a gravity offset up to -300 K below 500 Myr; low-mass cutoff (5 MJ or lower, SNIPPET).

5. Consistency test. Simulate a 20 pc volume with the mass function and constant birthrate, evolve with a track, map to type, and compare the counts per Teff bin with Kirkpatrick 2021 over 600 to 2250 K, where the census is complete (their completeness limit: all types except >= T8 and Teff < 600 K); tolerance from the Poisson error of each bin and the quoted alpha error. The observed space-density spike at 1200 to 1350 K (READ) is a test of the track.

6. Flags. None: all constraints are field brown dwarfs.

7. Group. Builder B.

### 2.4 h2 and h3: hot-Jupiter regimes at 1500 and 2200 K

**Question 4.**

1. Ensemble. The edges are not observed; the observational anchors are WASP-18 b (design cites Arcangeli 2018, not re-read here), the high-resolution detections of Fe, Ca, Mg and ions in WASP-76 b, WASP-121 b, WASP-189 b (Gandhi 2023, READ) and the JWST emission spectra (WASP-121 b, SNIPPET).

2. Physical bracket and what was computed. The regimes are outputs of Gibbs-energy minimisation. I wrote a minimiser (element potentials, scipy) and fed it the NASA-9 polynomials in `easychem/thermo_easy_chem_simp_own.inp` from the easyCHEM 2.1.0 sdist (DATA; NASA CEA thermo.inp format, 115 species parsed, including TiO, VO, SiO, Na, K, Fe, Mg, FeH, H-, e-, Na+, K+, Fe+; no Ti+, Mg+, AlO, CaO or condensates). Inputs: H=1, O/H 4.90e-4 times Z, C/O 0.55, He/H 0.0851 (Asplund 2009 values from memory, RECALLED), ideal gas, ten neutral species, no ions, no condensation. Sanity checks: H2O formation free energy at 1500 K reproduces the standard value (-13.2 RT, COMPUTED); CH4 dominates at 800 K, 1 bar. Results (COMPUTED, Z = 1; Z = 10 changes T by less than 70 K except where stated):

| P (bar) | T where CH4 = CO | T where H2O is 90 percent of the CO-limited value | T where H2O is 50 percent | T where H2 is half dissociated |
|---|---|---|---|---|
| 0.001 | not run | 2331 K | 2567 K | 2574 K |
| 0.01 | not run | 2542 K | 2843 K | 2884 K |
| 0.1 | 966 K (Z=10: 892 K) | 2792 K | 3182 K | 3278 K |
| 1 | 1156 K (1052 K) | 3095 K | 3609 K | 3794 K |
| 10 | 1438 K (1280 K) | not run | 4164 K | 4499 K |
| 100 | 1907 K (1635 K) | not run | not run | not run |

At 1500 K, CH4/CO is 7e-9 at 1 mbar and 7e-3 at 1 bar, so below about 1 bar the gas is CO-dominated at 1500 K. The 1500 K edge of h1/h2 is therefore the equilibrium CH4-to-CO line at about 5 to 10 bar; at photospheric pressures (0.01 to 0.1 bar) the line sits near 900 to 1000 K, and quenching (part 05 s7.3) decides what is seen. The 2200 K edge of h2/h3 corresponds to only 1 to 4 percent H2O depletion at 1 mbar to 0.01 bar; real thermal dissociation (half the water) lies at 2570 to 3600 K for 1 mbar to 1 bar, and observed depletion in ultra-hot Jupiters also involves H- opacity and ionisation, which this subset cannot do (Ti+ and Mg+ absent). The dayside is hotter than Teq. The edge is a pressure-dependent output; a fixed 2200 K is not supported by the computed table.

3. Verdict. Derive, do not draw. The "weights" Phi((T-1500)/100) and Phi((T-2200)/100) have no probabilistic content. What can be tuned are the physical inputs (P-T profile, Z, C/O).

   Licence-compatible tools (the repo is CC0, LICENSE.md, so GPL code cannot be linked):

| Code | Licence | Source, status |
|---|---|---|
| easyCHEM 2.1.0 (a from-scratch NASA CEA-style Gibbs solver, Fortran, plus the thermo file) | MIT | PyPI metadata (DATA); data file origin is the NASA CEA polynomial format, public-domain status of the NASA/TP-2002-211556 polynomials is RECALLED |
| chemicals 1.5.2, thermo 0.6.1 | MIT | PyPI metadata (DATA) |
| pyFastChem 4.0.3, exo_k 1.3.2 | GPL-3.0 | PyPI metadata (DATA); not usable as a linked dependency |
| Cantera 3.2.0, GGchem | PyPI field empty (Cantera); licences RECALLED | to verify |
| NASA/JANAF polynomials | NASA CEA thermo.inp is a US government product, RECALLED public domain; the JANAF tables are a NIST product whose licence the project doc says is not public domain (body master 15.2 item 5) | NIST SRD needs permission or re-derivation |

4. Knobs: the physical inputs only. Targets for a consistency test: LExACoM retrieved C/O and [M/H] for hot and ultra-hot Jupiters (8 planets in the paper's UHJ subsample: C/O median 0.497, [M/H] mean 0.86 +-0.52), and the dissociation-limited water depth in WASP-121 b-type planets (not tabulated in the public record I reached).

5. Consistency test: for a set of simulated hot Jupiters with computed photospheric (T, P), compare the fraction in the "H2O-rich", "CO + H2O", and "dissociated" regimes to the retrieved detection of H2O (Tsiaras 2018: water in every one of the 16 of 30 planets with significant atmospheres, READ). No selection model beyond transmission signal-to-noise.

6. Flags: none.

7. Group. Builder C (one chemistry module serves h1 to h5, f8, g6, m7 to m9 and n1 to n4).

### 2.5 h4: TiO/VO kept, Kzz vs grain size

**Question 5.**

1. Ensemble.

| Data set | Source | N | Result | Status |
|---|---|---|---|---|
| Spiegel, Silverio and Burrows 2009, arXiv 0902.3995 | READ (abstract, introduction) | model of 4 planets | to loft Ti so TiO can cause inversions Kzz must be at least a few times 10^7 cm2/s for 0.1 micron grains, rising to nearly 10^11 for 10 micron; HD 209458 b needs ~10^9 to ~10^11 for 0.1 to 10 micron; WASP-12 b ~2e7 | READ |
| Gandhi et al. 2023, arXiv 2305.17228 | READ (abstract, section 3) | 6 ultra-hot Jupiters (Teq 2210 to 2641 K) | Ca, Na, Ti and TiO are under-abundant; Ti-bearing species are constrained only in WASP-189 b | READ |
| Prinoth et al. 2022, arXiv 2111.12732 | READ | WASP-189 b (Teq 2641 +-34 K) | first unambiguous TiO detection in a transmission spectrum, 4.9 sigma; earlier WASP-33 b claims not reproduced | READ |
| Hoeijmakers et al. 2024, arXiv 2210.12847 | READ (abstract) | WASP-121 b | detections of Ca, V, Cr, Mn, Fe, Co, Ni; no Ti or TiO; titanium cold-trap | READ |
| Pelletier et al. 2023, arXiv 2306.08739 | READ (abstract, discussion) | WASP-76 b | VO detected; species with Tcond below ~1550 K near solar, above strongly depleted; Teq of WASP-121 b given as ~2350 K | READ |
| Tsiaras et al. 2018, arXiv 1704.05413 | READ (abstract) | 30 planets, 600 to 2400 K | TiO and/or VO at 4 sigma in WASP-76 b, probably in WASP-121 b (later high-resolution work finds no TiO in either) | READ |
| Sing et al. 2024, arXiv 2405.11027 and Welbanks et al. 2024, arXiv 2405.11018 | READ (abstracts) | WASP-107 b (about 750 K) | Kzz 10^11.6 +-0.1 (Sing) vs 10^8.4 to 10^9.0 (Welbanks) from partly the same data; metallicity 43 +-8 vs 10 to 18 solar | READ |

2. Physical bracket. A cold trap needs a nightside or deep layer cooler than the Ti condensation temperature (Tcond near 1550 K or above for the Ti-bearing perovskite in Pelletier's and Gandhi's text). Retention requires vertical mixing faster than grain settling, which is Spiegel's Kmin(a), a requirement not a measured Kzz.

3. Verdict. (c) for the Kzz prior LU[1e7, 1e11]: no population measures Kzz, and the two WASP-107 b analyses differ by 2.6 to 3.2 dex for the same planet (verified vs verified: both stand). The design's WASP-107 b 10^11.6 is right for Sing 2024 but not unique. (b) for a single empirical constraint that can replace the Kzz-grain machinery: the Teq above which Ti-bearing gas survives. In the verified data, WASP-76 b (Teq 2210 K) and WASP-121 b (about 2350 K) have no TiO/Ti, MASCARA-2 b has no TiO, and WASP-189 b (2641 K) has TiO and Ti: the threshold lies between about 2350 and 2640 K, from N = 3 to 6 planets. The design states "TiO above 2000 to 2300 K, VO above 1700 to 1900 K"; the verified record puts the Ti threshold 100 to 300 K or more higher, and shows VO retained at 2210 K (WASP-76 b). Here the verified record overrides the design's partly verified edges.

   Observed rate: TiO detected in 1 of 4 clearly tested ultra-hot Jupiters (WASP-189 b yes; WASP-76 b, WASP-121 b, MASCARA-2 b no; WASP-33 b disputed), no larger compilation found. Wilson 95 percent interval for 1 of 4 is 0.05 to 0.70 (COMPUTED), not a tuning constraint.

4. Knobs: T_trap (threshold Teq for Ti retention): target 2350 to 2640 K bracket from N=3; the Kzz-grain pair stays unconstrained and degenerate: for fixed Teq, retention probability depends on Kzz and grain radius only through Kmin(a) = 10^(7.48 + 1.76 (log10 a + 1)) (design interpolation; reproduces Spiegel's 3e7 at 0.1 micron and 1e11 at 10 micron), so a detection rate cannot separate them.

5. Consistency test: simulated fraction of ultra-hot Jupiters with gaseous TiO against Teq, compared with the list above (yes at 2641 K; no at 2210, 2250, 2260, 2350 K). A step at T_trap reproduces it; the Kzz-grain form (0.44 at 1 micron in the design example) cannot be told from 0.25 or 0.70 at N=4.

6. Flags: none.

7. Group. Builder D.

### 2.6 n6: silicate vs haze aerosol crossover near 950 K

1. Ensemble. Gao 2020, arXiv 2005.11939: an aerosol microphysics model on a grid Teq 700 to 2800 K, gravities 4, 10, 25 m/s2, metallicity 1 and 10 times solar, compared with the 1.4 micron water-band amplitude A_H (difference of maximum and minimum radius between 1.3 and 1.65 micron) of the WFC3 sample of Fu et al. 2017 (arXiv 1709.07385, N=34 hot Jupiters to hot Neptunes, with planets of radius below 0.5 RJ removed by Gao). Result (READ): silicates dominate aerosol opacity above about 950 K, hydrocarbon haze below; water-band amplitude is non-monotonic with Teq (increasing cloudiness with Teq above 2300 K and between 1100 and 1600 K, decreasing between 1600 and 2300 K and below 1100 K); robust to gravity and metallicity. Fu 2017 (READ): Spearman 0.43 between A_H and Teq (binned 0.85), slope 0.0008 +-0.00016 per K, degenerate with mass, hot planets under-sampled. Estrela et al. 2022, arXiv 2206.09934 (READ, abstract): 62 planets, uniformly processed WFC3; two populations, a cloud-free one over a wide Teq range and a "partial cloud" one that goes from mostly cloudy near 500 K to mostly clear near 1500 K; cloud/haze minimum at 1460.86 (+316, -405) K. Parmentier et al. 2016, arXiv 1602.03088 (READ, abstract): Kepler offsets suggest a silicate to MnS transition near 1600 K. Sing 2016 (arXiv 1512.04341, READ abstract): a clear-to-cloudy continuum in 10 hot Jupiters.

2. Physical bracket. Silicate clouds form where Teq is high enough for Mg2SiO4 and MgSiO3 to condense at observable pressures, and photochemical hydrocarbon haze needs methane in the upper atmosphere, which needs T below roughly 1000 K (the equilibrium CH4=CO line computed above is 892 to 1156 K at 0.1 to 1 bar). The 950 K of Gao is therefore the CH4 line plus the silicate-cloud sink depth, not an independent constant. So compute the crossover from the chemistry module (builder C) instead of drawing it.

3. Verdict. (b). A published ensemble (N=34 and N=62) exists, but it constrains cloud amplitude, not which aerosol dominates: the 950 K is a model inference, and no width is published. The 100 K width is a design choice. Estrela's 1461 K minimum and Gao's 950 K are different quantities and do not contradict. Gao predicts mid-infrared silicate features; JWST/MIRI quartz-cloud detections (Grant 2023, WASP-17 b, RECALLED) would test it and were not reached.

4. Knobs: crossover Teq (target 950 K, model; observed amplitude minimum 1461 +316/-405 K) and width (no published target). Test amplitude: A_H vs Teq slope 0.0008 +-0.00016 per K (Fu 2017) over 700 to 2400 K.

5. Consistency test: simulate A_H for the population with cloud and haze opacities and compare slope and binned Spearman (0.43 raw, 0.85 binned) with Fu 2017; with the Estrela cloud-free fraction as a second test. A transmission-signal selection model is required (only transiting planets with H*>1 ppm at WFC3 precision were observed).

6. Flags: none.

7. Group. Builder D.

### 2.7 n7: haze yield LU[0.01, 1]

1. Ensemble. Horst 2018, arXiv 1801.06512 (READ abstract): lab haze for super-Earth and mini-Neptune gases at 300 to 600 K and 100 to 10000 times solar, rates 0.04 to about 10 mg/hr (250x), Titan simulation 7.4 mg/hr. Crossfield and Kreidberg 2017, arXiv 1708.00016 (READ abstract): 6 warm Neptunes, feature amplitude correlates with Teq or with H/He mass fraction at 95 percent, with more haze at low Teq or higher metallicity (not separable). Libby-Roberts 2020, arXiv 1910.12988 (READ abstract): two Kepler-51 super-puffs flat to 0.6 scale heights, aerosols above 3 mbar; supports the Teq hypothesis, disfavours the H/He one. Brande 2024, arXiv 2310.07714 ("Clouds and clarity", READ abstract): 15 exo-Neptunes, 200 to 1000 K, best described by low sedimentation efficiency f_sed about 0.1 and 100 times solar metallicity, intrinsic scatter about 0.5 scale height, attenuation up to 6 scale heights at 500 to 800 K.

2. Physical bracket. Haze mass production is limited by the photolysis rate of methane (Gao 2020 parameterise it by comparing methane upwelling with photodissociation, READ), so a yield above one is impossible and zero for carbon-free gas.

3. Verdict. (c). The lab rates were measured on sub-Neptune gases, not on hot giants, in mg/hr, not a fraction of photolysed carbon, so "y = LU[0.01, 1]" is not the lab quantity. The population statistics constrain a combination of haze, cloud, metallicity and mean molecular weight, so the yield cannot be separated: Crossfield and Kreidberg say explicitly that temperature and H/He fraction cannot be separated at N=6, Brande's intrinsic scatter of about 0.5 scale height at N=15 absorbs any single knob, and Fu 2017 finds the temperature trend degenerate with mass. Boss would have to supply a haze production rate per unit methane photolysis for hot gases (none reached) or agree to treat y as a nuisance with the A_H scatter as its prior.

4. Knobs: y (no constraint), with the target amplitude statistic being the Brande 0.5 scale-height intrinsic scatter for 200 to 1000 K planets (N=15, READ abstract).

5. Consistency test: A_H vs Teq for 200 to 1000 K simulated Neptunes against the Brande attenuation (up to 6 scale heights at 500 to 800 K); selection model: transit signal-to-noise cut.

6. Flags: the only lab calibration is a Titan simulation (7.4 mg/hr); Titan is a Solar System body, so any scaling of y to the Titan rate is Solar System-only.

7. Group. Builder D.

### 2.8 T12: giant interior, diluted vs compact core

1. Ensemble. The core-plus-envelope population fit (C25: Mcore 14.7 +1.8/-1.6 Earth masses, fZ 0.09 +-0.01) and T16's 57.9 M^0.61 Earth masses fit. T16 tests homogeneous vs central-core extremes on four planets and finds the model effect "much smaller than observational uncertainties" (READ); C25 says the metal distribution (core vs uniform) "leads to differences that are much smaller than the astrophysical scatter" (READ) and adopts 10 Earth masses in the core and the rest in the envelope. WASP-107 b core 11.5 (+3.0, -3.6) Earth masses (Sing 2024, READ). The core-mass placeholder log-uniform 5 to 20 is consistent with C25's 14.7.

2. Verdict. (c). The population does not distinguish diluted from compact cores; C25's "Mcore" is a population-level parameter of the mass-metal relation, not a measured core. Structure is known only for Jupiter and Saturn (gravity fields) and weakly for WASP-107 b: Solar System-only for the structural question. The generator should keep the continuous Mz and draw the core fraction as a nuisance parameter between "all in envelope" and "10 Earth masses in the core".

3. Knobs: Mcore 14.7 (+1.8/-1.6) for the mass-metal relation; no knob for dilution.

4. Flag: structure evidence for Jupiter and Saturn is Solar System-only (not re-read here; RECALLED from the project doc).

5. Group. Builder A.

### 2.9 Where a verified source overrides a partly verified one (or disagrees with another verified one)

| Place | Design value (status) | Verified value | Winner and why |
|---|---|---|---|
| f7 tail | DM2010 34 percent (READ) | B&F 2016 and SA18 (both READ): none above 0.8 in 849 and 99 hosts | The two are verified; B&F and SA18 are used because DM2010 relied on one blended line (Teske 2013, READ) |
| f5 note | "C25 gives lower Z for sub-Saturns" | C25: sub-Saturns have high Z (>= 0.2); C25 Zp/Z* is within 7 percent of T16 at 0.1 MJ and 1.4 to 1.6 times lower at 0.3 to 3 MJ | READ C25 overrides the note |
| f5 width | 0.3 RECALLED, 0.26 READ for the Mz fit | 0.26 dex stands for Mz; no scatter is printed for the Zpl/Zstar fit; C25 gives 0.66 multiplicative for both | the 0.26 stays as the bulk scatter; atmospheric scatter is 0.4 to 0.9 dex (COMPUTED) |
| g4 | gap in part 05 (design) | 500 to 700 K is T8 to T9; Y0 at 371 to 510 K (Beiler 2024, READ) | continuous map replaces the gap row |
| h4 | TiO above 2000 to 2300 K (design) | Ti-bearing species retained only at Teq 2641 K in the six-planet sample; absent at 2210 to 2350 K | READ overrides the design edges |
| h4 | WASP-107 b Kzz 10^11.6 (READ, Sing 2024) | Welbanks 2024: 10^8.4 to 10^9.0 (READ) | both verified, disagree by 2.6+ dex; the disagreement is reported, not resolved |
| h2, h3 | edges 1500 and 2200 K (design) | CH4=CO at 966 to 1438 K, 50 percent H2O at 2567 to 3609 K for 1 mbar to 1 bar (COMPUTED) | the computed table replaces fixed edges |

## 3. Sources actually verified

Each item opened (PDF, or abstract only where stated); N in brackets.

- Bulk metallicity: Thorngren 2016 arXiv 1511.07854 [47]; Chachan 2025 arXiv 2509.20428 [147, no data statement found]; Thorngren and Fortney 2019 arXiv 1811.11859 [403, table "available for download", host not found]; Welbanks 2019 arXiv 1912.04904 [19]; Sun 2024 arXiv 2402.08292 [19].
- Population retrievals: LExACoM, Lothringer 2025 arXiv 2510.26785 [66 in paper, 69 unique planets in the 2026-08-28 CSV]; CSV and MRT fetched from jlothringer.github.io (DATA, table licence not stated); code MIT (raw GitHub LICENSE).
- Stellar C/O: Brewer and Fischer 2016 arXiv 1608.06286 [849]; Suarez-Andres 2018 arXiv 1801.09474 [499]; Delgado Mena 2010 arXiv 1009.5224 [100 hosts, 270 comparison]; Nissen 2014 arXiv 1406.5218 and Teske 2013 arXiv 1309.6032 (abstracts and one section); Costa-Almeida 2026 arXiv 2606.30722 [290]; Xie 2026 arXiv 2604.06809 (abstract, 598 hosts).
- Planet C/O and metals: Line 2021 arXiv 2110.14821 [1]; JWST WASP-39 b abstracts arXiv 2211.10488, 10487, 10489, 10493, 10490 (metallicity 1 to 100x; C/O sub-solar to solar); Sing 2024 arXiv 2405.11027 and Welbanks 2024 arXiv 2405.11018 (abstracts, WASP-107 b).
- Hot Jupiters: Gandhi 2023 arXiv 2305.17228 [6]; Prinoth 2022 arXiv 2111.12732 [1]; Pelletier 2023 arXiv 2306.08739 [1]; Hoeijmakers 2024 arXiv 2210.12847 [1]; Spiegel 2009 arXiv 0902.3995 [4 model planets].
- Clouds and hazes: Gao 2020 arXiv 2005.11939 [model]; Fu 2017 arXiv 1709.07385 [34]; Estrela 2022 arXiv 2206.09934 [62]; Tsiaras 2018 arXiv 1704.05413 [30]; Sing 2016 arXiv 1512.04341 [10]; Parmentier 2016 arXiv 1602.03088; Crossfield and Kreidberg 2017 arXiv 1708.00016 [6]; Brande 2024 arXiv 2310.07714 [15]; Libby-Roberts 2020 arXiv 1910.12988 [2]; Horst 2018 arXiv 1801.06512 [lab] (abstracts for the last four).
- Brown dwarfs: Kirkpatrick 2021 arXiv 2011.11616 [525]; Filippazzo 2015 arXiv 1508.01767 [198]; astro-ph/0402475 [42]; Beiler 2024 arXiv 2407.08518 [23]; Beiler 2023 arXiv 2306.11807 [1].
- PyPI metadata (DATA): pyfastchem 4.0.3 GPL-3.0; exo_k 1.3.2 GPLv3; easychem 2.1.0 MIT; chemicals 1.5.2 MIT; thermo 0.6.1 MIT; cantera 3.2.0 field blank; petitRADTRANS 3.4.0 MIT.

## 4. Boss must supply

1. The Thorngren and Fortney 2019 table (403 planets) and the Chachan 2025 table (147 planets) as machine-readable files, to re-fit the bulk scatter; neither is in the PDFs.
2. A TiO/VO/Ti detection and non-detection compilation for ultra-hot Jupiters with more than six planets (a high-resolution atlas or a JWST sample paper); none was found.
3. A published haze production rate per unit CH4 photolysis for 700 to 1500 K H2-dominated gas, or a decision that n7 is a nuisance parameter.
4. A population measurement of Kzz or a decision to drop the LU[1e7, 1e11] draw.
5. Kirkpatrick et al. 2021 Tables 11, 14 and 15 as machine-readable files (ApJS online) for the Teff scatter by bin and the space densities per 150 K.
6. The structural-interior evidence (Juno and Cassini gravity fits, diluted core) from primary papers, for T12.
7. Brewer and Fischer 2016 Table 1, to draw stellar C/O against [Fe/H] directly.

## 5. Recalled or to verify

1. Asplund 2009 log abundances O 8.69, C 8.43, He/H 0.0851 used in the equilibrium script (RECALLED).
2. Public-domain status of NASA CEA polynomials (NASA/TP-2002-211556) and the Burcat set; Cantera 3.2.0 and GGchem licences; ExoComp CSV licence (RECALLED).
3. The spectroscopic definition of the Y class (NH3 feature) is RECALLED; the Teff of Y0 (about 460 K from the fit, 371 to 510 K observed) is READ.
4. Candidates not resolved or not read: Teske 2014 (arXiv 1403.6891); August 2023, Changeat 2022, Baburaj 2025 (numbers via ExoComp or Costa-Almeida only); Grant 2023 quartz clouds.
5. Table 13 columns for the ExoComp mass-metallicity fits (garbled extraction).
6. The three Kirkpatrick polynomial branches join continuously (COMPUTED: 1277.6 and 1277.7 K at SpT 8.75; 1167.9 and 1168.1 K at 14.75), so the printed coefficients were transcribed correctly.
