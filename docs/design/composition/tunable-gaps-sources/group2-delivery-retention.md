# Group 2: volatile delivery, retention and lava worlds (tunable-gaps source check, 2026-10-17)

Scope: part 09 rows b1, b2, b3, b5, b6, d1 to d6, i6, and body-master tables T5, T9, T11 (with sections 14.2 and 8 of the body master, items R2, R4, R7 of part 09 section 4). Labels READ, 2ND-HAND, DATA, COMPUTED, SNIPPET, RECALLED, NOT REACHED are those of `_worker-rules.md`. Nothing in the repo was edited. All numbers below were read in the arXiv or Zenodo files named in section 3 unless a label says otherwise.

## 1. Summary and verdicts

The public record is thick for the shoreline (d rows) and thin for delivery. For the shoreline there is a published probabilistic fit with posteriors and a labelled planet sample [R0364] and an empirical line [R3037]. Both contradict the repo's slope-4 K window, which rests on 13 Solar System bodies only. For delivery, formation distance is not observable, so the only source of "water fraction as a function of formation distance" is population synthesis (Bern NGPPS). Its output is public (Zenodo 10527170, open, no licence field) and I extracted a mass-dependent migrant probability from it (COMPUTED, section 2.2). It can be checked against observed densities of 34 M-dwarf planets and seven TRAPPIST-1 planets, but not against any measurement of delivery itself. For lava worlds exactly two small planets have a volatile atmosphere detection (55 Cnc e, TOI-561 b), none has a verified bare-rock-vapour atmosphere. For super-Mercuries the published fraction is 0/20 to 5/22 depending on sample, and a mock test shows the outliers can arise from noise alone.

| Row | What it is | Verdict | One-line reason | Status |
|---|---|---|---|---|
| b1, b2, b3 | dry / low / high delivery weights vs x = a_form/a_snow | (b), with SS flag | weights can only come from NGPPS-type models; class payloads (0, 0.05, 0.5 wt% H2O) are Earth delivery scenarios | READ (Ortenzi), COMPUTED |
| b5 | continuous H2O fraction inside b2/b3 | (c) | no ensemble measures sub-percent water; TRAPPIST-1 gives 2 to 7 wt% with CMF degeneracy | READ (Agol) |
| b6 | ejected rogue water | (c), and not a separate draw | no composition of ejected planets exists; derive from the planet's own record | READ (Ma 2016 via part C) |
| d1 to d5 | shoreline regime weights | (b) | two-class posterior exists (N = 917 labelled, 15 small-planet labels); five-way split has no data | READ |
| d6 | pressure caps | (c); a derived cap, not a probability | only upper limits (about 1e-5 to 10 bar by composition) plus one 0.01 to 100 bar detection | READ |
| i6 | volatile-bearing vs bare vapour lava worlds | (b) | N = 2 of 2 volatile-bearing, 0 bare-vapour; Beta(3,1) posterior | READ |
| T5 | iron branch (super-Mercury, coreless) | (b) for p_SM as an upper limit; (c) for coreless | 0/20, 1/21, 1/11, 5/22 by sample; no coreless count | READ |
| T9 | born beyond snow line / migrant / in situ dry | (b) | NGPPS public data gives P(water-rich or migrant | M); observed check N = 34 | COMPUTED, READ |
| T11 | water delivered inside snow line | (c) | the only bracket is Earth's water budget; no exoplanet ensemble | READ |

Usability-rule overrides (verified over partly verified or snippet):
1. The repo doc says Meni-Gallardo and Palle keep slope 4 with a higher zero point (snippet). READ: they adopt slope 5.77 and intercept -4.35 (EECS). The slope-4 line through 55 Cnc e is shown and rejected because Solar System moons fall on the wrong side. Verified wins.
2. The repo doc says TRAPPIST-1 b and c are "bare rock or under 0.1 bar" (snippet). READ: b is airless-consistent [R1570, R1678] but a thick CO2 atmosphere with haze also fits; both planets: surface pressure at or above 1 bar strongly disfavoured [R1570]. The 0.1 bar figure is not in the papers I read. Verified wins.
3. Part 09 d1 cites [R4862] as the shoreline source (snippet). READ: the I proportional to v^4 line is "drawn by eye" in their Figure 1; the paper gives no K. K is the repo's own fit.
4. Part 09 R2 names "[R4335] data" as a water fraction against distance source. READ: that paper is "From stellar nebula to planets: the refractory components"; it tabulates the refractory mass fraction against distance, not water fractions, and mentions water five times. The 2015 paper (1507.01343) is a Research Note on Fe/Si, Mg/Si, C/O in stars versus planets.
5. The repo doc rates 55 Cnc e "likely CO2/CO atmosphere [S, Hu 2024]". Upgraded to READ: Hu et al. confirm a volatile atmosphere rich in CO2 or CO, 0.01 to 100 bar, and rule out a tenuous rock-vapour atmosphere.

## 2. Detail by subgroup

### 2.1 Answers to the four specific questions (short form)

1. Water fraction against formation distance outside population synthesis: no. Formation distance is not an observable. Observations give water fraction against mass, radius, period and stellar mass ([R2774], N = 34; [R0039], N = 7; [R3515], one planet). Population synthesis is the only source of the function; using it as the weight source is allowed under the no-fixed-content rule because it is a physics model with disc priors, not a copy of the Solar System, provided each weight is labelled MODEL and tuned only through the checks in 2.2. The statistics it can be checked against are the density gap and the water-world mass range of M-dwarf planets, the TRAPPIST-1 water fractions, and the radius-valley position. All of these are degenerate with photoevaporation (see the valley caveat in 2.2).
2. Rocky planets with an atmosphere detection or airless constraint: 15 with a direct label in [R0364] (2 positive, 13 negative; table in 2.3) plus 6 labelled ambiguous; Meni-Gallardo and Palle list 17 M-dwarf planets with reported non-detections (many from HST transmission, weak) and 4 tentative detections. K can be tuned to them, but the verified detections of 55 Cnc e and TOI-561 b need K above the repo window by factors of 12 to 190 (COMPUTED, 2.3), so a single K v^4 line cannot satisfy both the Solar System and the hot exoplanets. The test is a logistic classification of the labelled sample.
3. Super-Mercury share: 5 of 22 ([R0032], heterogeneous), 1 of 21 ([R3546], homogeneous), 1 of 11 [R3978], 0 distinct of 20 [R0541]. See 2.5.
4. Ultra-short-period or lava worlds with an atmosphere constraint: 7 USPs (P below 1 d) have a label (2 positive, 5 negative), see 2.3 and 2.4. Of these only the two positives are hot enough for a molten dayside (Teq above about 1700 K, the Berta-Thompson cut). That gives 2 of 2 volatile-bearing, not a fraction in any robust sense.

### 2.2 Delivery: b1, b2, b3, b5, b6, T9, T11

**Ensemble.**
- Bern NGPPS: I = [R1268, R1269], [R1268] (model, READ title); II = [R1268, R1269], [R1269] (solar-like stars, 1000-system populations, DACE identifiers; READ); III = [R3955] (warm super-Earths and cold Jupiters); IV = [R0602] (low-mass stars). The candidate list's "[R0602] NGPPS III" is NGPPS IV. NGPPS II sets the ice-poor/ice-rich boundary at 1 % volatiles by mass in the core and uses planetesimal densities 3.2 (rock) and 1.0 g/cm3 (icy) (READ).
- [R0600, R0601], [R0601] (Nature Astronomy), evolved NGPPS planets with a water-in-envelope equation of state; raw data "accessible at dace.unige.ch under identifier NG76" (READ). Its derived data is open on Zenodo 10527170 (code_upload_v2.zip, 136.7 MB, access "open", licence field empty; record created 2024-01-18). DACE front pages answer (HTTP 200) but the population tables are JavaScript-rendered and I did not test an API or download; DACE licence terms: NOT REACHED.
- Observational checks: [R2774] (Science): 43 transiting M-dwarf planets below 4 Re, 34 precisely characterised (mass better than 25 %, radius better than 8 %); three populations: rocky (mean density 0.94 +/- 0.13 of scaled Earth, radius 1.21 +/- 0.28 Re), "water worlds" taken as 1:1 rock to water by mass (0.47 +/- 0.05, 1.97 +/- 0.28 Re), puffy sub-Neptunes (0.24 +/- 0.04, 2.85 +/- 0.63 Re); density separation at 0.65; rocky planets exist from 0.3 to about 10 Me and water worlds start at 2 to 3 Me. [R0039]: TRAPPIST-1 water mass fractions with an Earth-like iron share: b 2.8, c 2.3, d 4.4, e 2.9, f 4.5, g 6.4, h 5.5 wt% (each with +2 to +4.5 / -1.2 to -3.1 errors), but 0.002 to 14 wt% when the core mass fraction is set to 50 %; only 88 % of linear fits to water against period have a negative slope. [R3515]: Kepler-138 d (1.51 Re, Teq about 350 K) volatile-rich, about 11 % volatiles by mass (SNIPPET value; abstract says "large fraction"). [R2111] (ApJL 939 L19): migration model reproduces valley with 1.4 Re planets mostly rocky and 2.4 Re mostly water-ice rich. [R4528]: rocky cores about 3 Me (max about 5), icy cores peak about 10 Me. [R1477]: valley near 1.6 Re at 0.8 Msun and 2.0 Re at 1.2 Msun (SNIPPET). [R2772]: no planets of 2.2 to 3.8 Re above 650 Earth fluxes (READ abstract).
- Delivery inside the snow line: [R3681] (READ): a forming giant scatters planetesimals inward, "whenever a giant planet forms it invariably pollutes its inner planetary system with water-rich bodies"; they match Earth's water, which is debated from under one to over ten oceans. The body master's 1.2e-4 to 5e-3 is the earlier worker's arithmetic (2.5e-3 to 25e-3 Me of C-type material times 5 to 20 % water), not a number printed by Raymond and Izidoro.
- [R3194] ([R3194], ApJL 920 L1) is a pebble-accretion explanation of why M dwarfs host more transiting super-Earths. It has no water-content result, so it does not say what the candidate list claims.
- [R3378] (Sci. Rep. 10, 10907, [R3378]; not on arXiv): Table 1 volatile contents are "based on three volatile delivery cases investigated for Earth" in their ref. 21 (READ in the Nature page text; table cells not rendered in what I fetched, values taken from the earlier part 09 reading). The dry/low/high payloads are therefore Earth scenarios.

**Computed from the public NGPPS population (NG76 nominal, 1000 systems, evolved to 5 Gyr; COMPUTED; column meaning taken from the authors' own plotting script `plot_biased_compo.py`, not from DACE documentation).** In the window 0.5 < P < 100 d and 0 < R < 20 Re the file has 2554 planets: 1727 rock, 502 water (envelope fully water, Zhomo = 1), 325 H/He. Fraction water among rock+water against total mass: 1 to 2 Me 0.02 (N = 518); 2 to 3 Me 0.09; 3 to 5 Me 0.23; 5 to 8 Me 0.49; 8 to 12 Me 0.75; 12 to 20 Me 0.91; 20 to 60 Me 0.92. Logistic fit over 0.5 < M < 30 Me (N = 2076, 490 water): p_mig = 1/(1 + exp(-6.28 (log10 M - 0.788))); centre 6.1 Me (bootstrap 5.9 to 6.3), slope 6.0 to 6.6 per dex. The bootstrap is sampling error inside one population, not model uncertainty. Water planets have water mass fraction 0.03 to 0.49 (5 to 95 %), median 0.28; 82 % have at least 10 wt%, none above 0.51. The repo's "migrant probability 0.14 to 0.67 above 3 Me, zero below" is a coarse step; the NGPPS shape is a smooth logistic with a centre about 6 Me. The formation-start distance is not in the evolved file, so P(class | x) itself cannot be extracted from it; DACE's formation-time tables or the six track directories in the Zenodo archive would be needed.

**Physical bracket.** Snow line: 170 K blackbody gives 2.69 AU at 1 Lsun; viscous-disc fits give 3.3, 2.6, 1.9, 1.2, 0.5 AU for 1.0, 0.8, 0.6, 0.4, 0.2 Msun with 0.2 dex scatter ([R3193], via part C, READ there); NGPPS median 2.81 AU (READ via part C). Ice-rich cores about 50 wt% water in pebble models; dry inside unless delivered (part C summary of Venturini, Bitsch, Ciesla 1 to 3 %, Raymond and Izidoro).

**Verdict.** b1 to b3: (b). A physics-model ensemble exists and is open, but the classes (0, 0.05, 0.5 wt%) are Earth scenarios and the only observational statistic that bears on them is the water-world fraction against mass, which the migration parameters of T9 control, not the b weights. b5: (c). b6: (c); rogue water should not be its own draw: the ejected planet keeps whatever T9/T11 gave it; Ma 2016 ([R2783], via part C) gives ejection statistics and "no composition". T9: (b). T11: (c).

**Knobs and targets.**
| Knob | Proposed form | Target statistic | Target value and uncertainty | Valid range |
|---|---|---|---|---|
| p_mig centre and slope (T9) | logistic in log10 M; centre 0.79 +/- 0.02, slope 6.3 +/- 0.3 per dex (NGPPS, 1 Msun) | NGPPS output itself, then Luque and Palle mass range of water worlds | water worlds from 2 to 3 Me; rocky up to about 10 Me (no error given) | 1 to 30 Me, 1 Msun host |
| M_ice,min stellar-mass scaling | 3 (M*/Msun)^0.48 from two NGPPS IV anchors (about 1 Me at 0.1, about 3 Me at 1 Msun; figure reads) | Luque and Palle 2 to 3 Me for M dwarfs | as left | 0.1 to 1 Msun |
| water share W of migrants | uniform 0.03 to 0.5, median 0.28 (NGPPS) | Luque and Palle density 0.47 +/- 0.05 (assumes 1:1) | density ratio 0.47 +/- 0.05; Agol b to h 2 to 7 wt% | rocky-water boundary |
| x_c, width 0.4 in s(x) (b1 to b3) | centre 1.0, width 0.4 | none; formation distance not observed | NO TARGET | n/a |
| delivered water inside the snow line (T11) | log-uniform 1.2e-4 to 5e-3 | none beyond Earth | NO TARGET | n/a |

**Consistency test.** Simulate a 1 Msun and an M-dwarf population with the generator's T9/b rules, apply the Luque and Palle selection (mass better than 25 %, radius better than 8 %, transiting, M-dwarf host) and compare the water-world and rocky fractions and the density modes with 0.94, 0.47, 0.24 (tolerance their stated errors). A transit-detection completeness and a radial-velocity mass-precision cut are needed on the simulated side. The valley caveat: the same observed gap is reproduced by evaporation of rocky cores ([R2772]; [R0600, R0601] combine both), so valley position alone cannot tune delivery.

**No-fixed-content flags.** Dry/low/high payloads (Ortenzi Table 1): Earth delivery scenarios, Solar System only. T_snow 170 K and the 2.69 AU: a blackbody value, not Solar System, but the nebula pressure behind the bracket is unsourced. T11 range: Earth's water budget and Solar System C-type asteroids are the only constraint, flagged. b6 defaults (0.35, 0.8, LU[1e-4, 0.1]): tuning.py defaults with no data, Earth ocean mass 2.3e-4 is the only anchor.

### 2.3 Retention: d1 to d6

**Ensemble.**
- [R4862] (ApJ 843, 122; READ): about 590 exoplanets with radius, mass and period from exoplanet.eu as of 26 August 2016 plus Solar System bodies; the I proportional to v_esc^4 line is "drawn by eye"; XUV scaling I_xuv = I (L*/Lsun)^-0.6 (their Eq. 27). An energy-limited shoreline should go as I_xuv v^3 sqrt(rho), distinct from the apparent v^4 (READ abstract). No K and no uncertainty are given.
- [R0364] (v2, 19 March 2026; READ): logistic shoreline in (bolometric flux f, v_esc, L*), fitted to the Solar System plus exoplanets, ages above 750 Myr, magma-ocean planets (Teq above 1700 K) excluded from the main fit. Samples: no-magma-ocean 869 atmosphere-labelled and 48 airless; any-temperature 1026 and 48; no-magma-ocean-no-CO2-freezeout 858 and 17. Labels: planets above 2.2 Re count as atmosphere-bearing by radius; small exoplanets are labelled only by direct measurement (list in the table below). Posterior (no magma ocean): log10 f0/f_earth = 2.62 (+0.45/-0.31); p = 5.9 (+0.61/-0.43); q = 1.17 (+0.28/-0.20); ln w = -1.40 (+0.32/-0.30), w = 0.247 dex (+0.092/-0.065), 5 to 95 % width w95 = 1.46 dex (a factor 28.7). Solar System alone: p = 3.76 (+0.73/-0.51); exoplanets alone: p = 6.37 (+3.1/-2.2); each consistent with 4. With TOI-561 b included: p = 6.67 (+0.77/-0.58), q = 1.54. Atmospheres are disfavoured for Earth-size planets in the habitable zone below log10 L* = -2.23 (+0.18/-0.21), about M4V.
- [R3037] (MNRAS, accepted 2026 June 15; READ): ExoAtmospheres database of 274 explored planets and 58 species; the empirical shoreline log10 I_xuv = 5.77 log10 v_esc - 4.35, anchored on Mars, GJ 9827 d, L 98-59 d, GJ 3090 b and Pi Men c (tentative detections); 55 Cnc e then falls on the line. Eris and Titan on it; Pluto is the exception. They state ultra-short-period planets that keep atmospheres "may follow a different retention criterion".
- Ji, [R2162] (ApJ 992, 198; READ abstract): hydrodynamic escape with initial volatile inventory 0.01 % to 1 % of planet mass; the inventory range alone prevents a unique shoreline, and stellar XUV and age uncertainty widen it further.
- [R3445] (READ abstract): historic XUV fluence of mid-to-late M-dwarf planets is 2.1 to 3.1 times the canonical scaling; only the largest terrestrial planets of such stars plausibly retain atmospheres.
- Rocky Worlds DDT: the [R3696] working-group report ([R3696]; READ) recommends about 500 hours of JWST, 15 to 20 rocky M-dwarf planets with Teq 200 to 600 K, 15 micron eclipses; it states the aim is to "define the cosmic shoreline". That is the report, not the program page (STScI page: NOT REACHED).

**Labelled small-planet sample (all READ in the paper named).**
| Planet | Host | Teq (K) | M (Me) | R (Re) | P (d) | Result | Source id | Label |
|---|---|---|---|---|---|---|---|---|
| 55 Cnc e | G8V | 1996 (MG25); about 2000 (Hu) | 8.8 (Hu); 8.58 (MG25) | 1.95 | 0.74 (MG25) | volatile atmosphere, CO2 or CO, 0.01 to 100 bar; rock-vapour-only ruled out; dayside depth varies (Patel: 4.5 micron brightness 873 to 2256 K) | 2405.04744, 2407.12898 | READ |
| TOI-561 b | thick-disk, 10 Gyr, Teff 5372 K | about 2300 | 2.24 (error printed as 2.0 in my extraction; probably 0.20, verify) | 1.4195 | 0.4466 | dayside inconsistent with bare rock; thick volatile envelope; density 4.3 +/- 0.4 | 2509.17231 | READ |
| LHS 3844 b | M dwarf | not read | not read | 1.3 | 0.46 (11 h) | dayside 1040 +/- 40 K; above 10 bar excluded at 3 sigma; bare rock fits | 1908.06834 | READ |
| GJ 367 b | M1V | 1367 | 0.63 | 0.70 | 0.32 | CO2 above about 1 bar, outgassed above 10 mbar (reducing) or 0.01 mbar (oxidising) excluded; 1 bar excluded for a wide range | 2401.01400 | READ |
| GJ 1252 b | M dwarf | 1089 +/- 69 | 1.32 | 1.213 | not read | Tday 1410 K; surface pressure below about 10 bar; no significant atmosphere | 2208.09479 | READ |
| TOI-1685 b | M2.5V | about 1000 | not read | not read | not read | featureless, 0.98 +/- 0.07 of blackbody; thick high-weight atmospheres disfavoured; thin ones not excluded | 2412.03411 | READ |
| GJ 486 b | M dwarf | 696.3 | 3.0 | 1.3 | not read | Tday 865 +/- 14 K, 0.97 of maximum; thin below 1 % H2O or 1 ppm CO2; Earth or Venus excluded | 2408.15123, 2305.00868 | READ |
| GJ 1132 b | M3.5 | 589 (MG25) | 1.65 (MG25) | 1.16 (MG25) | 1.63 (MG25) | Tday 709 +/- 31 K; 1 bar with 1 % H2O and Venus-like with 1 ppm CO2 or H2O excluded | 2408.13340 | READ (planet values 2ND-HAND) |
| GJ 3929 b | M dwarf | Tmax 737 +/- 14 | 1.13 +/- 0.09 | 1.09 +/- 0.04 | 2.616 | Tday 782 +/- 79 K; CO2-rich above 100 mbar excluded above 3 sigma (a four-visit follow-up, [R0869], leaves room for thin atmospheres: SNIPPET) | 2508.12516 | READ |
| TRAPPIST-1 b | M8 | S about 4.15 | 1.374 | 1.116 | 1.511 | Tday 490 +/- 17 K, no nightside flux; at or above 1 bar disfavoured; Ducrot: airless ultramafic or thick CO2 with haze | 2303.14849, 2412.11627, 2509.02128 | READ |
| TRAPPIST-1 c | M8 | S about 2.2 | 1.308 | 1.097 | not read | Tday 369 +/- 23 K; tenuous O2 or airless; no thick CO2 (flux ratio 421 +/- 94 ppm) | 2306.10150, 2509.02128 | READ |
| LHS 1140 c, LTT 3780 b, LTT 1445 A b, TOI-1468 b | M dwarfs | not read | not read | not read | not read | labelled airless by Berta-Thompson from eclipses (their refs) | 2507.02136 | 2ND-HAND |

Berta-Thompson leave Kepler-10 b, Kepler-78 b, K2-141 b, LHS 1478 b, TOI-431 b and L 98-59 b unlabeled (ambiguous), which are the other USP and hot-rock cases. Counting the rows above that have a label: 2 positive, 13 negative; USPs (P below 1 d) among them: 55 Cnc e, TOI-561 b, LHS 3844 b, GJ 367 b, TOI-1685 b (period RECALLED, about 0.67 d), GJ 1252 b (RECALLED about 0.52 d), LTT 3780 b (RECALLED about 0.77 d).

**Physical bracket.** Energy-limited escape gives I_xuv proportional to v^3 sqrt(rho) (Zahnle and Catling, READ); hydrodynamic escape broadens it with inventory (Ji, READ); runaway-greenhouse and CO2 freeze-out edges at 0.24 f_earth (194 K) and magma ocean at 1400 f_earth (Berta-Thompson, READ).

**Computed check of the repo's line against the labelled sample (COMPUTED; r = S/(K v^4), K = 3.85e-4 centre, S = (Teq/278.3)^4 with the Teq printed above, X = S so no M-dwarf boost; v = 11.186 sqrt(M/R)).**
| Planet | v (km/s) | S | r at K = 3.85e-4 | r range over K window | Observed |
|---|---|---|---|---|---|
| 55 Cnc e | 23.8 | 2646 | 21.6 | 12.6 to 27.7 | atmosphere |
| TOI-561 b | 14.1 | 4665 | 311 | 181 to 399 | atmosphere |
| GJ 367 b | 10.6 | 582 | 119 | 70 to 153 | airless |
| GJ 1252 b | 11.7 | 234 | 32.9 | 19 to 42 | airless |
| GJ 486 b | 17.0 | 40.3 | 1.26 | 0.73 to 1.61 | airless |
| GJ 3929 b | 11.4 | 18.4 | 2.85 | 1.7 to 3.7 | airless |
| GJ 1132 b | 13.3 | 20.1 | 1.65 | 0.96 to 2.1 | airless |
| TRAPPIST-1 b | 12.4 | 4.0 | 0.44 | 0.26 to 0.56 | airless-consistent |
The two atmosphere-bearing planets sit at r = 22 and 311, deep on the "airless" side; for the slope-4 line to put them on the atmosphere side K must exceed 8.3e-3 (55 Cnc e) and 0.118 (TOI-561 b), that is 12 to 190 times the top of the window. TRAPPIST-1 b is airless only if the XUV boost X/S exceeds 2.3 (range 1.8 to 3.9 over the K window); Pass 2025 gives 2.1 to 3.1, so this one is consistent at the edge. The same arithmetic with the Meni-Gallardo line (G-star hosts, S used for I_xuv) gives 55 Cnc e a ratio 0.68 to the line, matching their statement that it lies on it. The repo line is 10 to 200 times below the Berta-Thompson L = 1 shoreline for v = 4 to 20 km/s (at 11.186 km/s: 6.0 against 418 Earth fluxes), because the Solar System fit leaves out the exoplanet anchors.

**Verdict.** d1 to d5: (b). The tunable object is the two-class probability P(atmosphere) = 1/(1 + 10^(Delta/w)) with the Berta-Thompson posteriors; the five-way split of d2 to d5 is not constrained by any measurement. d6: (c), and it is a derived cap (min of supply, vapour pressure, and escape), not a probability; the only positive pressure is 55 Cnc e at 0.01 to 100 bar.

**Knobs and targets.**
| Knob | Form | Target statistic | Target (published uncertainty) | Range |
|---|---|---|---|---|
| log10 f0 | intercept at v_esc = Earth, L* = 1 | atmosphere labels, N = 917 | 2.62 (+0.45/-0.31) | 0.5 to 6.5 Earth fluxes in log; v_esc 0.1 to 30 Earth |
| p | slope in v_esc | same | 5.9 (+0.61/-0.43); SS-only 3.76 (+0.73/-0.51) | same |
| q | slope in L* | same | 1.17 (+0.28/-0.20) | log10 L* -3.5 to 0.5 |
| w | logistic width, dex | same | 0.247 (+0.092/-0.065) main; 0.090 (+0.086/-0.048, from ln w = -2.41 +0.67/-0.76, COMPUTED) in the zoomed sample | |
| Repo K window | LU[3.0e-4, 6.6e-4], slope 4, no width | 13 Solar System bodies | NO exoplanet target; width 0.34 dex against 1.46 dex (w95) | |
| Pressure caps (d6) | 0, 1e-5, 0.01, 0.1 bar | upper limits only | 0.1 bar (GJ 3929 b, CO2), 0.01 mbar (GJ 367 b, oxidising) | |
| Hot-planet branch | separate rule for Teq above 1700 K | 55 Cnc e, TOI-561 b | 2 of 2 atmosphere | |

**Consistency test.** Draw a planet population with the generator, compute (f, v_esc, L*) per planet, apply the labelled-sample rule (direct measurement or radius above 2.2 Re; exclude ages below 750 Myr and Teq above 1700 K), and compare log-likelihood and the fitted (f0, p, q, w) with the published posterior within its quoted intervals. A selection model is needed: JWST targets were chosen for high spectroscopy metrics and bright hosts, and "airless" only means no atmosphere above about 1 to 10 bar for hot tidally locked planets ([R2434], cited by Berta-Thompson; not read), so the simulated side must apply the same detection threshold. The published data cannot test the 0.01 and 0.1 bar caps of d3 and d4.

**No-fixed-content flags.** K window: the only constraint is 13 Solar System bodies (Titan forces K above 2.26e-4, Pluto and Triton raise it to 3.0e-4, Ganymede caps it at 6.56e-4); this is the largest single violation in the group. Pressure caps d4/d5 values (0.1 bar, none): Mars, Titan, Venus, Earth only. The Solar System sample is still legitimate as a test point inside the Berta-Thompson fit, which uses it jointly with exoplanets.

### 2.4 Lava worlds: i6 and the volatile-bearing fraction

**Ensemble.** The two positives above (Hu 2024; [R4329]) and the negatives that are USPs but not molten (GJ 367 b, GJ 1252 b, LHS 3844 b, TOI-1685 b). Hu 2024 ([R2023], Nature 630, 609): NIRCam and MIRI 4 to 12 micron emission spectrum; rules out a tenuous rock-vapour atmosphere; best fits a CO or N2 background with 0.01 to 100 bar and CO2 mixing ratio above 1e-5. [R3448]: five NIRCam occultations, depth from non-detection to 96 +/- 8 ppm (2.1 micron) and 119 (+34/-19) ppm (4.5 micron); different visits favour thin CO/CO2, silicate vapour, or a flat line, so the planet's atmosphere is variable or episodic; this does not contradict Hu on the presence of volatiles but shows a snapshot can mislead. [R4329]: four NIRSpec 3 to 5 micron eclipses; dayside well below the about 3000 K of bare rock; rejects "complete atmospheric desiccation". Berta-Thompson (READ): K2-141 b phase curve ambiguous ([R4108, R4905], not read), Kepler-10 b and Kepler-78 b optical only. Demographics: [R3904] (READ abstract-level: 106 USP candidates from about 200,000 Kepler targets); [R4720] (READ): occurrence (0.51 +/- 0.07) % for R above 0.84 Re and P below 1 d, 1.1 +/- 0.4 % for M dwarfs and 0.15 +/- 0.05 % for F dwarfs, USPs almost always below 2 Re (these occurrence numbers are 2ND-HAND, quoted by Winn from Sanchis-Ojeda). USP radii and a missing 2.2 to 3.8 Re population above 650 Earth fluxes: [R2772] (READ).

**Physical bracket.** [R2380] ([R2380]; the repo's name list "Kite, Fegley, Schaefer, Ford" is wrong: authors are Kite, Fegley, Schaefer, Gaidos): log10 P_eq = k1 exp(k2 Ts), k1 = -34.7, -60.1, -60.2 and k2 = -0.00112, -0.00122, -0.00101 for vapour fraction below 0.01, about 0.5, above 0.9, for 1400 to 2800 K, valid for an initial bulk-silicate-Earth composition (READ). Magma pools are well mixed below about 2400 K and variegated above (READ). The coefficients assume Earth's silicate, a Solar System anchor that applies to the vapour pressure only.

**Verdict.** i6: (b). Posterior for the fraction f_vol of molten-dayside rocky planets that keep volatiles, uniform prior: 2 of 2 gives Beta(3,1), median 0.79, 95 % interval 0.29 to 0.99 (COMPUTED). If K2-141 b and Kepler-10 b are counted as bare, 2 of 4 gives median 0.5, interval 0.15 to 0.85. The row's c[0.1, 0.9] is therefore not excluded but its lower end is poorly supported; the data say nothing about the 0.1 floor. The USP population is a mix of molten and non-molten planets, so the "volatile-bearing fraction" is only defined above the molten threshold.

**Knobs and test.** f_vol with the Beta posterior; target is the count 2 of N_molten, N_molten is 2 to 4 in the papers read; the test compares the simulated fraction of lava worlds with a CO2/CO signal against 2/2 (tolerance: the Beta interval). Selection: only bright, high-metric targets were observed, a model of JWST target choice is needed on the simulated side, and the molten threshold (1700 K) is itself a model choice.

**No-fixed-content flags.** The vapour-pressure fit assumes bulk-silicate Earth; Hu and Teske are exoplanet data, fine.

### 2.5 Iron branch: T5

**Ensemble.** [R0032] (Science; READ): 32 planets below 10 Me around 27 FGK stars, 22 below the radius gap treated as rocky, five "super-Mercuries" (K2-38 b, K2-106 b, K2-229 b, Kepler-107 c, Kepler-406 b), correlation slope above 4 (OLS-type), heterogeneous literature stellar data. [R3978] (READ): 11 probable rocky planets, 10 of 11 consistent with the host-star core mass fraction at 2 sigma; Kepler-107 c is the only 2 sigma super-Mercury, K2-229 b does not meet the 1 or 2 sigma threshold. [R3545] (READ): 33 planets in the sample, cmf of stars 0.32 (+0.14/-0.12) against 0.24 (+0.33/-0.18) for planets, maximum enrichment about cmf 0.8; also an iron-depleted group of highly irradiated planets (mean Fe/Si 0.86 against 1.69 for their stars) and "we lack a reliable theory that forms massive iron depleted planets". [R0541] (READ): 20 planets, 1 to 1.8 Re, homogeneous stellar data; OLS slope 1.3 +/- 1.0, orthogonal-distance slope 5.6 +/- 1.6; 75 % within 1 sigma of the host; no distinct super-Mercury population. [R3546]/26, [R3546] (v2; READ abstract and text; the candidate "[R3546]" is this paper): 21 planets, 20 stars, homogeneous; slope 0.94 (+1.02/-1.07), intercept -0.02 (+0.31/-0.29); super-Mercuries drop from four to one (Kepler-107 c); outliers recovered in mock populations about 15 % of the time; at least 150 planets needed. Individual: [R0456] (Kepler-107 c about 12.6 g/cm3, giant impact, SNIPPET value); [R3915] (K2-229 b, 2.59 +/- 0.43 Me, 1.165 +/- 0.066 Re, cmf Mercury-like: SNIPPET values, title READ); K2-106 b was reported dense (about 13 g/cm3, 80 % iron) but the reanalysis (11 USPs, Gaia DR2) gives 7.72 Me, 1.712 Re, 8.5 +/- 1.9 g/cm3 and a lower iron fraction (2ND-HAND, read in the reanalysis [R3791]).

**Counts and intervals (COMPUTED, exact binomial 95 %).** 1/21 = 0.048 (0.001 to 0.238); 5/22 = 0.227 (0.078 to 0.454); 1/11 = 0.091 (0.002 to 0.413); 0/20 distinct = 0 (0 to 0.168). The body-master range p_SM 0.05 to 0.23 is the two point estimates 1/21 and 5/22.

**Selection effects.** Only planets with mass precision better than 20 to 30 % and radius better than 20 % are used; radial-velocity mass measurement favours massive, dense, short-period planets around bright FGK stars (the sample is mostly USPs and hot planets); the radius cut at the valley is model-dependent; stellar abundances are heterogeneous in Adibekyan; [R3546] reports that mass-radius scatter alone reproduces the observed outliers about 15 % of the time, so even p_SM = 0 is compatible with the data.

**Verdict.** T5 super-Mercury branch: (b) as an upper limit, not a measured share: p_SM below about 0.24 (95 %, homogeneous sample). The coreless branch has no count (c): the only hints are the iron-poor outliers [R3545] and the TRAPPIST-1 fit to an iron-depleted or core-free structure ([R0039], READ). Derive the coreless branch by physics from the water fraction ([R1248], via part B) rather than drawing it.

**Knobs and test.** p_SM (prior Beta from 1/21), CMF range for the branch 0.5 to 0.8 ([R3545] cap about 0.8), slope m of planet against star iron share (0.94, +1.02/-1.07), intercept (-0.02, +0.31/-0.29), scatter. Test: draw rocky planets with host abundances, add the sample's mass and radius errors, recompute core mass fractions, and compare the number of 2 sigma outliers with 1 of 21 (and the recovered-outlier rate 15 % from Plotnykov's mock test). Needs injection of the measurement errors and the mass-precision selection.

**No-fixed-content flags.** Mercury is only a test label; the CMF cap 0.8 is from exoplanets. No flag.

## 3. Sources actually verified

| Source | Id | N or content | Access, licence | Label |
|---|---|---|---|---|
| [R1268, R1269], NGPPS I | [R1268] | model paper | arXiv; DACE populations | READ (title, scope) |
| [R1268, R1269], NGPPS II | [R1269] | 1000-system populations; ice boundary 1 % | arXiv; DACE | READ |
| [R3955], NGPPS III | [R3955] | warm super-Earths and cold Jupiters | arXiv | READ (title) |
| [R0602], NGPPS IV | [R0602] | low-mass stars | arXiv | READ (title); numbers via part C |
| [R0600, R0601] | [R0601]; Zenodo 10527170 | NG76 evolved population | Zenodo open, licence not stated; DACE NG76 | READ, DATA (analysed) |
| [R0364] | [R0364] | 917 labelled (Solar System plus exoplanets), code linked (not fetched) | arXiv | READ |
| [R3037] | [R3037] | ExoAtmospheres database, 274 planets | arXiv; database page NOT REACHED | READ |
| [R4862] | [R4862] | about 590 exoplanets plus SS | arXiv | READ |
| Ji 2025 | [R2162] | models | arXiv | READ (abstract, intro) |
| Pass 2025 | [R3445] | M-dwarf XUV | arXiv | READ (abstract) |
| [R3696] | [R3696] | working group report | arXiv | READ |
| JWST/Spitzer rocky-planet papers | 1908.06834, 2303.14849, 2306.10150, 2305.00868, 2408.15123, 2408.13340, 2508.10579, 2208.09479, 2405.04744, 2407.12898, 2509.17231, 2401.01400, 2412.03411, 2508.12516, 2412.11627, 2509.02128 | see 2.3 table | arXiv | READ (abstracts and parameter tables) |
| [R2774] | [R2774] | 34 M-dwarf planets | arXiv | READ |
| [R0039] | [R0039] | TRAPPIST-1 | arXiv | READ |
| [R3515] | [R3515] | Kepler-138 d | arXiv | READ (abstract) |
| [R2111, R2772, R4528] | 2210.05595, 2008.05513, 1604.05220 | valley and migration | arXiv | READ (abstracts) |
| [R3681] | [R3681] | scattered planetesimals | arXiv | READ |
| [R3194] | [R3194] | pebble model, no water result | arXiv | READ (abstract only; confirms the mismatch) |
| [R3378] | [R3378] | 28,080 evolutionary scenarios, Table 1 Earth cases | Nature page text; cells not rendered | READ (text), values 2ND-HAND via part 09 |
| [R4335], 2015 | 1312.3085, 1507.01343 | refractory components; stars vs planets | arXiv | READ (title, water mentions) |
| [R2380] | [R2380] | vapour pressure fit | arXiv | READ |
| [R3904, R4720] | 1403.2379; 1803.03303 | 106 USP candidates; occurrence | arXiv | READ / 2ND-HAND |
| [R0032, R0541, R3545, R3546, R3978] | 2102.12444, 2011.08893, 2010.06480, 2409.08361, 2511.17717 | 22, 11, 33, 20, 21 planets | arXiv | READ |
| [R0456, R3915] | 1902.01316; 1805.08405 | Kepler-107 c; K2-229 b | arXiv | READ (title); numbers SNIPPET |
| Rodriguez [R3791] | [R3791] | K2-106 b reanalysis | arXiv | READ |

## 4. Boss must supply

1. NGPPS formation-stage tables (NG76 or the NGPPS II populations on DACE): planet-by-planet formation-start distance, core ice and rock mass, final a and mass, to compute P(dry/low/high | a_form/a_snow, M). The evolved file I analysed has final a only. Also the DACE terms of use (not stated on the pages I could reach).
2. A published exoplanet-system N-body survey with giants and water mass fractions of inner planets (the T11 population); [R3681] is Solar System only and [R0803] (via part C) gives 1 to 3 % for low-mass hosts. I could not name a verified general survey.
3. The [R0364] code and label table (linked in the paper as the </> symbol; GitHub not reachable), and the IAC ExoAtmospheres database export, to rerun the fit with the generator's own population.
4. JWST results beyond the ones read: K2-141 b [R4108, R4905], LHS 1478 b (August 2025), TOI-431 b, the TOI-561 b phase curve [R0487], and the Rocky Worlds DDT results as they appear, to grow the lava-world count above N = 2.
5. The Kepler/K2 mass-radius sample with homogeneous stellar abundances beyond 21 planets (Plotnykov says 150 are needed) before p_SM or a coreless fraction can be tuned.
6. A decision on the nebula water partial pressure that fixes T_snow (161 to 213 K bracket) with a source.

## 5. Recalled / to verify

- Periods of GJ 1252 b (about 0.52 d), TOI-1685 b (about 0.67 d), LTT 3780 b (about 0.77 d) and P of TRAPPIST-1 c (about 2.42 d): RECALLED, used only to count USPs.
- Which reference number 21 in [R3378] is (probably [R1931]): RECALLED.
- TOI-561 b mass uncertainty (0.20 vs the printed 2.0); GJ 1252 b host type and period; LHS 3844 b mass and Teq; planet values for LTT 1445 A b, LHS 1140 c, TOI-1468 b, LTT 3780 b: not read.
- Numbers from snippets only: Kepler-138 d volatile share about 11 %, Fulton and Petigura valley positions 1.6 and 2.0 Re, Kepler-107 c 12.6 g/cm3, K2-229 b mass and radius, the four-visit GJ 3929 b follow-up [R0869].
- [R2434] (heat redistribution limit of 1 to 10 bar for hot daysides), [R0803, R1248], Ma 2016 and [R3193] were taken from part B/C or from Berta-Thompson's text, not re-read.
- Candidate "[R4335] and 2015" as a water source: does not say that (see section 1).
