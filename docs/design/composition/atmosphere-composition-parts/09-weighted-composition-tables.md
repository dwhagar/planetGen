# Part 09: weighted composition tables (draw tables for the generator)

Design requirement (project owner): for every composition category there is a table of options with weighted probabilities; each row has a class gate; the generator draws among the rows that pass the gate with weights normalised to 100 percent; a weight may be a function of star, formation distance, mass, temperature, pressure, redox and age; every weight carries its source. Rows without a sourced weight are marked `DESIGN Rn` with the nearest physical bracket so the generator can use a seeded range meanwhile. This part restates the recipes of part 05 (sections 2 to 13), part 06 (clouds, hazes, dust) and part 07 (ground radioisotopes) in that format, plus part 04 for escape and UV. Labels READ, 2ND-HAND, DATA, COMPUTED, SNIPPET, RECALLED, NOT REACHED are those of `_worker-rules.md`. Solar System bodies appear only as TEST calibration points. Nothing in the repo was edited.

## 1. Schema

**Columns.** `id | option | gate | weight w | weight source | status | notes`. The category is the section heading. Status is one of: **SOURCED** (the gate edges or weight numbers are taken directly from a read source), **DERIVED** (computed from a physical formula given in the row or a cited part), **DESIGN Rn** (no sourced weight; the bracket is in the weight column and `Rn` is the research item of section 4). A weight source ends in its verification label.

**Symbols.** M mass (Earth masses; MJ Jupiter masses), R radius (Earth radii), Teq, Ts, Tday (K), P (bar), x = a_form / a_snow, Phi standard normal CDF, N(m,s) normal, LN(med,s) log-normal with log10 x ~ N(log10 med, s), LU[a,b] log-uniform, c[lo,hi] a DESIGN constant, dIW log10 fO2 above the iron-wuestite buffer, X = xuv_index (insolation S times the star's cumulative-XUV ratio to the Sun's).

**Gate vocabulary (names the real repo classes; `src/planetgen/tuning.py` PLANET_CLASSES, read in full).** The repo has 20 classes: A B C D E F G H I J K L M N O P Q S T V. Type t (terrestrial): A B C D E F G H K L M N O P Q S V. Type g: I J T. Zone flags h (hot), e (ecosphere), c (cold), r (rogue-eligible): A,B hot only; C,D h e c r; E F G H K L M N O P Q V e only; I e c r; J h e c r; S h e c r; T c r. `HABITABLE_PLANET_CLASSES` = E F G H K L M O P Q V. `MOON_BLACKLIST` = Q S V. Shorthands used below: ROCKY = the 17 type-t letters; AIR = A B E F G H K L M N O P Q V (repo atmosphere text present); NOAIR = C D S (atmosphere None); GIANT = I J T; ROGUE = r flag (C D I J S T); BD = free-floating brown dwarf (rogue `mass_bin` 'brown-dwarf', 13 to 80 MJ, no class letter); COMET (bound or interstellar), AST (asteroid-field family: carbonaceous, stony, metallic, icy, basaltic, mixed, dust, collisional), WD and PSR hosts (white dwarf; `neutron-star:millisecond`). `[P]` after a letter (e.g. R[P]) marks a planned letter that is NOT in PLANET_CLASSES but is specified in `docs/design/habitability/atmospheres-retention-and-classes.md` (R sub-Neptune, U icy dwarf, W small ocean body, X ocean world, Y Titan-like, Z lifeless temperate). Physical clauses: M, R, Teq, P, T, host SpT, zone, age.

**Corrected (GEN.840, the class cut-over GEN.646).** The machine gates of the rows no longer name class letters: the class is derived from the finished body (`planet_classes.classify`), so it cannot be an input of the draw that finishes the body. A row names the kinds of body it supports (`terrestrial`, `giant`, `brown-dwarf`, `asteroid`, `comet`; a planet or moon is `terrestrial` or `giant` by its stored body type) and the stored conditions that must also hold. In the vocabulary above, read ROCKY as the kind terrestrial; AIR as terrestrial with a surface pressure of at least 1e-4 bar (10 Pa, the line under which `classify` calls a body airless); GIANT as the kind giant; LAVA as terrestrial with the zone and temperature conditions of the row; a sub-Neptune as a giant of 1 to 30 Earth masses (the envelope rows e and the flag d7 are giant rows, so a terrestrial body keeps no envelope until the formation module makes it a giant); the water-delivery rows b1 to b5 and b6 support terrestrial bodies and giants of at most 35 Earth masses (the classifier's ice-giant limit). A letter in the gate text of a row below is the old wording and is kept only as the history of the row.

**Usability rule (Boss, 3:02 AM EDT).** A row built on verified or partly verified material is usable if nothing verified contradicts it; where it is contradicted the verified source wins. Rows whose only support is RECALLED, SNIPPET, UNSUPPORTED or NOT REACHED stay DESIGN rows with a seeded bracket. See the master file section 9.

**Combining a conditional weight with the gate.** For a row i with input vector x: W_i = g_i(x) * w_i(x), where g_i is 1 if every clause of the gate holds, else 0. An edge clause written `T>=T0 (s)` is softened to Phi((T-T0)/s) instead of 0/1 (only rows that list s). A `Phi-window (s)` weight on a range lo to hi is Phi((hi-T)/s) - Phi((lo-T)/s) (a missing end drops its term). Rows whose gates overlap share probability in the overlap (examples c, j).

**Normalisation.** p_i = W_i / sum_j W_j over the rows of that table. If the sum is 0 the table's fallback is used (the airless / bare / absent row; stated per table). **Hierarchy:** a table marked `L2, parent=id` is drawn only after its parent row won; joint p = p(parent) * p(row | parent). A weight may call another table's probability, written p[id]; draw order: host star, a, b, d, c, e, f, g, h, i, j, k, l, then m, n, o, p (m to p read earlier draws).

**Continuous quantities.** A PARAM row gives dist, parameters and truncation, e.g. dIW ~ N(mu, s). The same quantity as a BINS table: p_k = [Phi((hi_k-mu)/s) - Phi((lo_k-mu)/s)] / [Phi((hi-mu)/s) - Phi((lo-mu)/s)] over the truncated range; for LN use log10 edges; for LU p_k = ln(hi_k/lo_k)/ln(hi/lo). Either form may be coded; both are given for metallicity (f6) and carbon-to-oxygen ratio (f7).

**DESIGN rows.** The weight is a constant c[lo,hi]; the seeded draw picks c log-uniform in the bracket (uniform if lo=0) before normalising. Worked examples use the geometric midpoint, stated.

**Seeded draw.** u from the project's single seeded random wrapper; take the first row whose cumulative p is at least u; a PARAM row uses further seeded draws. No other random source.

## 2. Category tables

### a. Mantle redox (magma-ocean stage, later upper mantle)

| id | option | gate | weight w (symbols above) | weight source | status | notes |
|---|---|---|---|---|---|---|
| a1 | Magma-ocean-stage fO2 (PARAM) | ROCKY, MOON, ROGUE rocky; any M (the clamp holds -2 below 0.01) | dIW ~ N(mu, sigma); mu = clamp(-2+1.25(log10 M+2), -2, 0.5) + clamp(log10 M, 0, 0.7324); sigma = 0.5 + 0.4 clamp(log10 M+1, 0, 1) (GEN.482b); floor on sigma 0.22 = 4 log10(1+0.005/0.037) | anchors Sossi 2020 via part 08 s2: +0.5 at about 1 Me, about -2 for core-formation-set small bodies (READ); floor COMPUTED; interpolation, 0.5 DESIGN | DERIVED | sets the first atmosphere; Fe3+/sumFe 0.037+-0.005; GEN.482 (docs/design/composition/magma-ocean-redox-against-mass.md): the mean is inside the first-hand ranges at the Moon, Mars and Earth, the width and the clamp above 1 Me are not sourced |
| a2 | Later upper mantle (PARAM) | ROCKY with air; tectonic or long-lived | dIW ~ N(3.5+2 clamp(log10 M, -3, 0), 1.0), clip [-6,6]; the draw takes max(a2, a1) (GEN.482b) | repo `atmosphere.py` REDOX_* ([D] fitted by eye to Moon, Mars, Earth) | DERIVED | far above meteoritic building blocks (part 05 s3.2); GEN.482 (docs/design/composition/magma-ocean-redox-against-mass.md): inside one standard deviation of the three first-hand anchors, no source above 1 Me |
| a3 | Building-block redox (PARAM) | AST, C, D, small bodies | dIW ~ U[-1.4, 0] | Schaefer and Fegley 2017 Table 2 via part 05 s3.2: EH about -1, eucrite -1.4 to -0.8, ranking CI>CV>H>EH (READ); uniform shape DESIGN | DESIGN R3 | bracket sourced, shape not |
| a4 | White-dwarf debris (PARAM) | host WD, rocky | dIW ~ N(-1, 0.5), bracket [-1.37, -0.6] | Doyle 2019 arXiv:1910.12989 data (READ); centre and width RECALLED design (part 05 s10) | DESIGN R10 | |
| a5 | Redox class R (reduced) | ROCKY, any dIW | Phi((0-dIW)/0.5) | cut 0: Liggins 2022 classes (READ) and part 05 s3.4 crossovers (COMPUTED); width 0.5 DESIGN | DERIVED | parent of c1, c2 |
| a6 | Redox class I (intermediate) | same | 1 - w[a5] - w[a7] | same | DERIVED | parent of c3 |
| a7 | Redox class O (oxidised) | same | Phi((dIW-2.5)/0.5) | cut 2.5: H2S/SO2 crossover near +3 at 1200 K (COMPUTED, part 05 s3.4); equals repo REDOX_OXIDIZED_DELTA_IW | DERIVED | parent of c4 |

**Example a (COMPUTED, Python).** Planet M=1.0: a1 mu=+0.50; integrating a5 to a7 over the draw N(0.5,0.5) gives a5 R 0.240, a6 I 0.758, a7 O 0.002; the later stage a2 (mu 3.5, s 1.0) gives a5 R 0.001, a6 I 0.185, a7 O 0.814. Planet M=0.05: a1 mu=-1.13, classes a5 R 0.944, a6 I 0.056; later stage mu=0.90: a5 R 0.211, a6 I 0.713, a7 O 0.076. For a single drawn dIW=+0.5: R 0.159, I 0.841, O 0.000.

### b. Delivered volatile inventory (dry / low / high) and water mass fraction, with snow-line dependence

| id | option | gate | weight w (symbols above) | weight source | status | notes |
|---|---|---|---|---|---|---|
| b1 | Dry delivery (H2O about 0, CO2 0.018 to 0.1 wt%) | ROCKY, MOON, ROGUE; x known | 1 - P(not dry); x = a_form/a_snow; a_snow = sqrt(L_form/Lsun) (278.6/T_snow)^2 AU, T_snow=170 K gives 2.69 AU; P(not dry) a logit in ln x, log10 M and a hinge at x = 0.6 (the three nested events of `ngpps_formation_water.json`) | classes and wt%: Ortenzi 2020 Table 1 (READ), except the dry CO2: the delivered 0.6 wt% is replaced by a log bracket from the stored 0.018 to the High case's 0.1 wt%, tuned to the bulk silicate Earth's 100 ppm carbon (GEN.813, `volatile-inventory-against-planet-totals.md`); the weight is weight: the Bern NGPPS population against the formation distance and the mass (10506 planets of 0.5 to 30 Earth masses; COMPUTED, GEN.802, `ngpps_formation_water.json`; the class thresholds that map the model's water fraction to Ortenzi's cases are DESIGN); T_snow bracket 161 to 213 K COMPUTED from the water-ice p_sat fit (part 05 s12.1) at nebula partial pressure 1e-9 to 1e-5 bar (that pressure unsourced) | DERIVED | a_snow ranges 1.8 to 3.4 AU over the bracket |
| b2 | Low delivery (H2O 0.05, CO2 0.01 wt%) | same | P(not dry) - P(high or rich) | Ortenzi 2020 Table 1 (READ); weight: the Bern NGPPS population against the formation distance and the mass (10506 planets of 0.5 to 30 Earth masses; COMPUTED, GEN.802, `ngpps_formation_water.json`; the class thresholds that map the model's water fraction to Ortenzi's cases are DESIGN) | DERIVED | stored in mantle 0.015 and 0.0022 wt% |
| b3 | High delivery (H2O 0.5, CO2 0.1 wt%) | same | P(high or rich) - P(rich) + (1 - c) P(rich) | Ortenzi 2020 Table 1 (READ); weight: the Bern NGPPS population against the formation distance and the mass (10506 planets of 0.5 to 30 Earth masses; COMPUTED, GEN.802, `ngpps_formation_water.json`; the class thresholds that map the model's water fraction to Ortenzi's cases are DESIGN); c is the plateau of b4 | DERIVED | stored 0.045 and 0.005 wt%; D_H2O 0.01 |
| b4 | Water-rich, ocean world (>=10 wt% water) | ROCKY, S, V, or sub-Neptune (I, R[P]) | c P(rich), c in [0,1], 1 as fitted | hycean needs >=10 wt% water (Madhusudhan 2021 READ); refractories >50 wt% of solids caps water+ices at 50 wt% (Thiabaud 2014 READ); the water fraction f >= 0.1 share is weight: the Bern NGPPS population against the formation distance and the mass (10506 planets of 0.5 to 30 Earth masses; COMPUTED, GEN.802, `ngpps_formation_water.json`; the class thresholds that map the model's water fraction to Ortenzi's cases are DESIGN); c bracketed 0 to 1 because published shares run from 3 percent (Rogers 2025) to 87 percent | DERIVED | feeds e3, e4 |
| b5 | Continuous H2O fraction given b2 or b3 (PARAM) | b2 or b3 won | log10 W_H2O ~ U[log10 5e-4, log10 5e-3] | part 05 s3.1 rule between Ortenzi cases (READ values, rule DESIGN) | DESIGN R2 | bins low/high 50/50 |
| b6 | Ejected rogue, water-rich | ROGUE | chance 0.35 (terrestrial), 0.8 (sub-Neptune); fraction LU[1e-4, 0.1] | repo `tuning.py` ROGUE_WATER_* ([D] defaults) | DESIGN R2 | Earth oceans 2.3e-4 |

**Example b (COMPUTED).** L=1 Lsun, a_snow 2.686 AU, c = sqrt(0.05*0.4) = 0.141. Formed at 1.2 AU (x=0.447, s=0.118): b1 dry 0.7782, b2 low 0.2077, b3 high 0.0139, b4 water-rich 0.0002. Formed at 4.03 AU (x=1.5, s=0.734): b1 dry 0.0671, b2 low 0.3701, b3 high 0.5099, b4 water-rich 0.0529.

### c. Outgassed atmosphere regime for rocky bodies (L2 under the redox class)

| id | option | gate | weight w (symbols above) | weight source | status | notes |
|---|---|---|---|---|---|---|
| c1 | R, hot: H2-H2O-CO (CO2, H2S) | L2, parent=a5; T>=950 K (s 50) | 1 | speciation at 1200 K, 1 bar, dIW 0: H2 .48 H2O .30 CO .11 CO2 .051 H2S .033 (part 05 s3.4, COMPUTED) | DERIVED | CH4 1.7e-5 there |
| c2 | R, cool and dense: CH4-H2O-NH3-H2 | L2, parent=a5; T<=1100 K (s 50); P>=10 | 1 | CH4 about 0.24 at 1000 K 100 bar; NH3 .015 to .02 at 800 K 100 bar (part 05 s3.4, COMPUTED) | DERIVED | Ar_40 rule in c5 |
| c3 | I: H2O-CO2 (CO, N2) | L2, parent=a6 | 1 | Schaefer 2012 Table 3: H2O 49 to 80 %, CO2 12 to 42 %, 1000 K, 100 bar and 2000 K (part 05 s3.4, READ) | DERIVED | crust gives HF, HCl (m5, m6) |
| c4 | O: H2O-CO2-SO2 | L2, parent=a7 | 1 | BSE 62 % H2O, 20 % SO2, 12 % CO2 (Schaefer 2012, READ); dIW+4 at 1200 K 1 bar: H2O .78 CO2 .16 SO2 .021 (COMPUTED) | DERIVED | |
| c5 | After cooling: H2O condensed, CO2-N2 residual | after c1 to c4; Ts<647 K (RECALLED T_crit of water) and ocean mass | 1 | Sossi 2020: CO-rich gas then CO2-N2 near 97:3 once H2O condenses; about 99 % of H2O and 35 % of CO2 dissolve in the melt, N2 about 0 % (part 08 s2, READ) | DERIVED | Ar_40 rule: n = f_degas f_Ar K40 (exp(lam t)-1), f_degas LU[0.05,0.7] (part 05 s3.6; DESIGN R3) |
| c6 | Steam atmosphere kept | Ts>=647 K | 1 - w[c5] | same | DERIVED | fallback row |

**Example c (COMPUTED).** dIW=+0.5, Ts=1000 K, P=100 bar. L1 (a5/a6/a7) = 0.159 / 0.841 / 0.000. L2 for R: c1 gate Phi((1000-950)/50) = 0.841, c2 gate Phi((1100-1000)/50) = 0.977 (P>=10 passes), so p(c1|R) = 0.463, p(c2|R) = 0.537; joint c1 0.073, c2 0.085, c3 0.841.

### d. Atmospheric retention regime (airless, exosphere, thin, moderate, dense) and primary H/He

| id | option | gate | weight w (symbols above) | weight source | status | notes |
|---|---|---|---|---|---|---|
| d1 | Airless | ROCKY, MOON, AST | P_K(r>30) | shoreline S_cs = K v^4, K LU[3.0e-4, 6.6e-4] fitted to 13 bodies (docs/design/habitability/atmospheres-retention-and-classes.md s3.1, [C]); underlying Zahnle and Catling 2017 (SNIPPET) | DERIVED | fallback row |
| d2 | Trace or exosphere (outgassed only) | same | P_K(2.5<r<=30) | same | DERIVED | |
| d3 | Tenuous | same | P_K(1<r<=2.5) | same | DERIVED | |
| d4 | Moderate | same | P_K(0.7<r<=1) | same | DERIVED | |
| d5 | Dense | same | P_K(r<=0.7) | same | DERIVED | |
| d6 | Pressure cap of the winning row (PARAM) | d1 to d5 | P_max = 0, 1e-5, 0.01, 0.1 bar for d1, d2, d3, d4; none for d5; P = min(supply, P_max, vapour) | smooth caps suggested in the project doc s3.6 (values DESIGN) | DESIGN R4 | |
| d7 | Primary H/He kept (flag outside the d1 to d5 draw) | M 1 to 30, lambda_H2=2.016 amu v_esc^2/(2 k T_exo)>=30 (T_exo = Teq + min(1500, 750 sqrt(F_xuv_rel)) K) | goes to table e (e1 or e2) | same doc s3.2 and s3.4 ([C]) | DERIVED | gate only |

**Example d (COMPUTED).** r = X/(K v^4), v = 11.186 sqrt(M/R) km/s, K log-uniform, so P_K(r in [r0,r1]) is the log-length of the K interval. Planet M=0.3, R=0.68 gives v=7.43 km/s. Sun-like host with X=1.5: d5 dense 0.000, d4 moderate 0.372, d3 tenuous 0.628, d2 trace 0.000, d1 none 0.000. Same planet, M-dwarf host with X=15 (the doc's factor of 10 in XUV exposure): d5 dense 0.000, d4 moderate 0.000, d3 tenuous 0.000, d2 trace 1.000, d1 none 0.000.

### e. Sub-Neptune / super-Earth envelope (L1 regime, L2 metallicity)

| id | option | gate | weight w (symbols above) | weight source | status | notes |
|---|---|---|---|---|---|---|
| e1 | Bare rock (envelope stripped) | S, V, I, T, R[P]; M 1 to 30; not d1-only | (1-p[b4]) (share lost); share lost = 1 - envelope left / f_env0 M, the envelope left from the hydro-based loss rate of Kubyshkina et al. 2018 integrated over the star's XUV history from the onset age (GEN.822); f_env0=c[0.05,0.2] | hydro-based rate: Kubyshkina et al. 2018 (READ); XUV history: Johnstone-type decay of `activity` (DESIGN); onset age and f_env0 tuned to the radius valley (GEN.509, GEN.822) | DERIVED | envelope fraction 0.5 to 10 % of mass |
| e2 | H/He envelope retained | same | (1-p[b4]) (1 - share lost) | same | DERIVED | goes to e5, e6 |
| e3 | Water world, steam/ocean, no H2 envelope | same, p[b4]>0 | p[b4] (1-eta), eta=c[0.1,0.5] | TEST: GJ 9827 d H2O>31 % (READ) | DESIGN R1 | |
| e4 | Hycean (H2 <=0.1 wt%, water >=10 wt%, bottom 300 to 400 K) | same | p[b4] eta | Madhusudhan 2021 (READ) | DESIGN R1 | radii up to 2.6 R_E at 10 M_E |
| e5 | Envelope metal mass fraction Z (BINS) | L2, parent=e2 | LU[0.01, 0.6]: bins 0.01-0.03, 0.03-0.1, 0.1-0.3, 0.3-0.6 | bracket from TEST planets: GJ 1214b >100x solar, TOI-270 d 58 wt%, K2-18b near solar to tens of x (READ); occurrence-weighted distribution NOT REACHED | DESIGN R1 | mean molecular weight mu = 2.4 + 3.1 log10(Z/0.01)/log10(58) (matches 2.4 and 5.5, part 05 s5) |
| e6 | L2 carbon carrier CH4 vs CO2/CO | parent=e2 | CH4: Phi((500-Teq)/100) | CH4 to CO near 500 K (400 to 600) for planets, Zahnle and Marley 2014 (READ); width = quoted half-range | DERIVED | N as N2; NH3 not seen (K2-18b) |

**Example e (COMPUTED).** M=5, core R=5^0.27=1.54, S=100: f_lost = 0.0156 of the planet mass. With f_env0=0.03 the ratio is 0.520. Formed at x=1.5 (p[b4] = 0.053), eta = 0.224 (bracket midpoint): e1 0.492, e2 0.455, e3 0.041, e4 0.012. L2 (e5) bin probabilities 0.268, 0.294, 0.268, 0.169 with mean mu 2.82, 3.70, 4.58, 5.26.

### f. Giant planets: temperature regime, metallicity bins, C/O bins, carbon carrier

Sudarsky class boundaries were not re-read in parts 05 to 07 (part 06 s2.6 marks them 2ND-HAND at best), so the regimes below follow part 05 s7.4.

| id | option | gate | weight w (symbols above) | weight source | status | notes |
|---|---|---|---|---|---|---|
| f1 | Cold giant, ice-giant envelope (H2O, NH3, H2S cloud-depleted; CH4 5e-4 m) | GIANT; Teq<200 K | Phi-window (s 100 K) | boundaries 200, 700, 1300 from part 05 s7.4 (COMPUTED/design) and Burrows and Sharp 1999: NH3 clouds about 200 K, water cloud about 350 K, N2 to NH3 about 700 K (READ); width DESIGN | DERIVED | |
| f2 | Warm giant, 200 to 700 K (CO quench-set 1e-9 to 1e-5) | GIANT | Phi-window | same | DERIVED | |
| f3 | T-like, 700 to 1300 K (CH4, NH3 1.2e-4 m, Na, K gas above clouds) | GIANT, BD | Phi-window | same | DERIVED | |
| f4 | Hot, Teq>=1300 K (go to table h; BD to table g) | GIANT zone h | Phi-window | same | DERIVED | |
| f5 | Metallicity m (PARAM) | GIANT | log10 m ~ N(log10[9.7 (M/MJ)^-0.45], 0.3), range [-0.5, 2.0] | centre: Thorngren 2016 Zpl/Zstar relation (READ); 0.26 dex scatter of the Mz fit (READ) is the bracket for 0.3 (0.3 RECALLED design) | DERIVED | Chachan 2025 gives lower Z for sub-Saturns |
| f6 | Metallicity bins of f5 (BINS) | same | edges log10 m = -0.5, 0, 0.5, 1.0, 1.5, 2.0 | same | DERIVED | |
| f7 | C/O (PARAM, mixture) | GIANT, BD | (1-f_tail) core + f_tail tail; core log10(C/O) ~ N(log10 0.55, 0.15) truncated to [0.2, 0.8]; tail U[0.8, 1.5]; f_tail = LU[0.0013, 0.34], default 0.03 | stellar C/O median 0.47 with no star above 0.8 (N=849; recounted on the full 1,592-star Brewer and Fischer table on the Asplund 2021 scale: median 0.49, 0.9 percent above 0.8, 97.5th percentile 0.68); planets (LExACoM, 69): median 0.57, tail 0.043 to 0.087 above 0.8 (tunable-gaps 3.4, COMPUTED from READ tables); core width 0.21 dex reproduces the 12 of 69 planets below 0.35 (GEN.516; the single planet WASP-77Ab 0.59 +/- 0.08, Line 2021, READ, gave 0.15) (READ, COMPUTED) | DERIVED R5 | bins <0.35, 0.35-0.7, 0.7-0.8, 0.8-1.0, >1.0 |
| f8 | Carbon carrier CH4 vs CO | GIANT, BD | CH4: Phi((T_trans-T)/100); T_trans 500 (planets), 1100 (BD) | Zahnle and Marley 2014 (READ); generator should evaluate the quench formula t_q with its own g, m, Kzz instead | DERIVED | HCN, C2H2 appear at C/O>=1 (Madhusudhan 2012 READ) |

**Example f (COMPUTED).** Giant M=0.3 MJ, Teq=950 K. L1: f2 warm 0.0062, f3 T-like 0.9936, f4 hot (to h) 0.0002. Zpl/Zstar = 16.7 (log10 = 1.222). Metallicity bins (sigma 0.3): 0.000, 0.008, 0.223, 0.596, 0.173; with sigma 0.26: 0.000, 0.003, 0.194, 0.662, 0.141. C/O bins with f_tail 0.03: 0.2-0.35: 0.106, 0.35-0.7: 0.747, 0.7-0.8: 0.117, 0.8-1.0: 0.009, 1.0-1.5: 0.021. An unmodified log-normal (0.55, 0.15 dex) would put 0.139 above 0.8; see conflicts.

### g. Brown dwarf spectral-type sequence (BD; free-floating, 13 to 80 MJ)

| id | option | gate | weight w (symbols above) | weight source | status | notes |
|---|---|---|---|---|---|---|
| g1 | Young hot, late-M-like | BD; Teff 2200 to 2950 K | Phi-window (s 100 K) | cap 2950 K = the BHAC15 temperature of the hydrogen-burning limit at its peak (`substellar_ceiling.HOTTEST_K`, GEN.494); lower edge DESIGN | DESIGN R15 | |
| g2 | L dwarf: H2O, CO (5e-4 m), Na, K, FeH, CrH; Fe and silicate clouds | BD; 1300 to 2200 K | Phi-window | part 05 s7.4 (COMPUTED); Burrows and Sharp 1999 (READ) | DERIVED | |
| g3 | T dwarf: CH4, H2O, NH3, Na, K; sulfide and salt clouds | BD; 700 to 1300 K | Phi-window | same | DERIVED | |
| g4 | Late T, NH3-rich (T8 to Y0 at 460 to 686 K), water cloud near 350 K deeper | BD; 460 to 686 K | Phi-window | the edges are the type-temperature map itself: T8 686 K and Y0 460 K, scatter 79 K (Kirkpatrick et al. 2021 Table 13, 525 dwarfs within 20 pc; READ); 500 to 700 K is T8 to T9, not a type of its own (GEN.494) | DERIVED R15 | |
| g5 | Y dwarf: CH4, water cloud | BD; Teff<500 K | Phi-window | part 05 s7.4 | DERIVED | |
| g6 | Carbon carrier CH4 vs CO | BD | CH4: Phi((1100-Teff)/100) | Zahnle and Marley 2014 1100+-100 K (READ) | DERIVED | |

**Example g (COMPUTED).** M=40 MJ, age 3 Gyr: Teff from L = 4e-5 Lsun (t/Gyr)^-1.3 (M/0.05 Msun)^2.64 (repo ROGUE_BROWN_DWARF_*, citing Burrows and Liebert; 2ND-HAND), radius 1 RJ, gives Teff = 848 K. Table: g3 T 0.9310, g4 T/Y 0.0687, g5 Y 0.0002. P(CH4 is the carbon carrier) = 0.994.

### h. Hot and ultra-hot Jupiters

| id | option | gate | weight w (symbols above) | weight source | status | notes |
|---|---|---|---|---|---|---|
| h1 | Warm-like, quench-enhanced CH4, NH3, HCN | GIANT zone h; Tday<1500 K | Phi-window (s 100 K) | part 05 s8 recipe; Moses 2011 (READ) | DERIVED | |
| h2 | CO + H2O + H2; CO2 rises with m | Tday 1500 to 2200 K | Phi-window | the edges are computed equilibrium lines at the photospheric pressure (methane equal to carbon monoxide; water 10 percent destroyed) with the NASA polynomials at solar metals and C/O 0.55, gas only and no ions; the softness of 100 K is DESIGN (GEN.518) (COMPUTED, DESIGN softness) | DERIVED R6 | |
| h3 | Thermally dissociated: H, OH, H-, Na, K, Fe, Mg, SiO, ions | Tday>2200 K | Phi-window | the edge is the computed temperature at which a tenth of the water is destroyed, at the photospheric pressure, gas only and no ions; the softness of 100 K is DESIGN (GEN.518) (COMPUTED, DESIGN softness) | DERIVED R6 | no mole fractions without FastChem or GGchem |
| h4 | TiO/VO gas kept (leaky cold trap) | h2 or h3; grain radius a (um) | P(log Kzz >= Kmin(a)); Kzz ~ LU[1e7, 1e11]; log Kmin = 7.48 + 3.52 (log10 a + 1)/2 | a Ti threshold between 2350 and 2641 K from N=5 planets (TiO seen at 2641 K, not at 2210, 2250, 2260 and 2350 K); no population measures Kzz, so none is drawn (GEN.519) (READ, COMPUTED) | DERIVED R6 | TiO above 2000 to 2300 K, VO above 1700 to 1900 K |
| h5 | C/O>=1 branch (HCN, C2H2; TiO, VO depleted) | all h | p[f7 bins >=1.0] | Madhusudhan 2012 (READ) | DERIVED | |

**Example h (COMPUTED).** Tday=2500 K: h2 CO-H2O 0.0013, h3 dissociated 0.9987. Tday=2150 K: h2 CO-H2O 0.6915, h3 dissociated 0.3085. TiO/VO kept for 1 um grains: log Kmin = 9.24, p = 0.440.

### i. Lava-world vapour regimes

| id | option | gate | weight w (symbols above) | weight source | status | notes |
|---|---|---|---|---|---|---|
| i1 | Type I: Na, O2, O, Fe | ROCKY zone h (A, B, S, C); Tp<2033 K (s 60) | Phi-window | Miguel 2011 text boundaries (READ); s covers the caption list, which differs by 2 to 155 K | SOURCED | prior if Tp undrawn: 189 of 193 Kepler candidates (0.979) |
| i2 | Type II: + SiO, Mg | 2033 to 2588 K | Phi-window | same | SOURCED | 3 of 193 (0.016) |
| i3 | Type III: O2, O, SiO columns equal | 2588 to 2890 K | Phi-window | same | SOURCED | 0 of 193 |
| i4 | Type IV: SiO > O2 | 2890 to 3168 K | Phi-window | same | SOURCED | 1 of 193 (0.005) |
| i5 | Type V: SiO, O, O2, Na | Tp>3168 K | Phi-window | same | SOURCED | 0 of 193 |
| i6 | Bare silicate vapour (volatile-free melt) vs volatile-bearing CO2/CO | ROCKY zone h, 1<=M<=10 | volatile-bearing c[0.1,0.9], bare 1-c | 2 of 2 molten-dayside ultra-short-period planets keep volatiles (N=2): the bracket is the 95 percent interval 0.29 to 0.99 of Beta(3, 1) (GEN.511) (READ, COMPUTED) | DERIVED R7 | P_eq: log10 P = k1 exp(k2 Ts), small vapour fraction k1=-34.7, k2=-0.00112 (Kite 2016, 1400 to 2800 K, READ) |

**Example i (COMPUTED).** Tp=2060 K: i1 I 0.326, i2 II 0.674; P_eq = 3.5e-04 bar. Tp=2500 K: i2 II 0.929, i3 III 0.071; P_eq = 7.8e-03 bar.

### j. Frost-equilibrium and evolved atmospheres (Titan-, Pluto-, Venus-, Mars-like)

| id | option | gate | weight w (symbols above) | weight source | status | notes |
|---|---|---|---|---|---|---|
| j1 | Venus-like dense CO2-N2 (97:3), SO2 trace | N, V, A, S[with air]; S_eff above the runaway limit (1.04 to 1.11 for a Sun-like star); dIW>=1; water lost | c[0.1,1] | CO2-N2 ratio: Sossi 2020 (part 08 s2, READ); runaway limit Kopparapu 2013/2014 ([R] in project doc); prior DESIGN | DESIGN R8 | S_eff limit depends on Teff (Kopparapu polynomial) |
| j2 | Mars-like thin CO2-N2-Ar | K, L, U[P]; M<0.5; 1<r<=2.5 (d3); 190<=Ts<=270 K | c[0.1,1] | same CO2-N2 ratio at lower pressure (Sossi, READ); TEST Mars CO2 .96 N2 .019 Ar .019 (Mahaffy 2013, READ) | DESIGN R8 | |
| j3 | Titan-like N2-CH4 (N2 0.5 to 2 bar) | Y[P], MOON; 57<=Ts<=125 K (s 5), r<0.7, S<=0.02, R 1500 to 4000 km | c[0.1,1] | T window: N2 ice p_sat 30.6 to 62.2 K for 1e-7 to 0.1 bar (part 06 s2.4, 2ND-HAND); TEST Titan CH4 1.5 to 5.7 % (Niemann 2010, READ) | DESIGN R8 | |
| j4 | Pluto/Triton-like frost equilibrium: N2 >99 %, CH4 0.3 %, CO 0.05 % | D, C, U[P]; 30<=Ts<=67 K (s 5), r<1 | c[0.1,1] | p_i = chi_i p_sat,i(Ts) (Raoult); TEST values (Mousis and Anderson 2024 2ND-HAND); N2/CO/CH4 vapour-pressure relations not reached | DESIGN R8 | CO column vs near-surface ratio differ 8x |
| j5 | Abiotic O2 world (10 to 100+ bar O2, or CO build-up) | Z[P]; host M dwarf, age>1 Gyr, HZ, 0.3<=M<=5 | c[0.05,0.5] | Luger and Barnes 2015 via Meadows 2018 (2ND-HAND): 240 bar O2 per ocean; Gao 2015 cap 50 % CO2 loss (READ) | DESIGN R8 | |

**Example j (COMPUTED).** R=2500 km, Ts=61 K, r=0.5, S=0.012, all priors at the midpoint 0.316: gates j3 and j4 both pass (edge overlap 57 to 67 K): j3 Titan-like 0.471, j4 frost-eq 0.529. (This reproduces the transition region where N2 is neither thick nor frozen out.)

### k. Comet and KBO coma volatile classes (ratios to water)

| id | option | gate | weight w (symbols above) | weight source | status | notes |
|---|---|---|---|---|---|---|
| k1 | CO2-dominated coma | COMET active (perihelion<3 AU, repo threshold), icy AST | 0.5 | 13 of 25 comets (HP22 Table 2); the class depends on the heliocentric distance (every coma beyond about 3.5 AU is CO-dominated) and 9 of the 13 rest on an upper limit on CO, so the share is inflated by an amount that cannot be corrected (GEN.520) (READ) | DERIVED R9 | |
| k2 | CO-dominated coma | same | 1/3 | 9 of 25 comets (HP22 Table 2); the heliocentric distance and the CO admission rule bias the share (GEN.520) (READ) | DERIVED R9 | Hale-Bopp CO 23 %, Hyakutake 18.2 % of H2O (READ) |
| k3 | Comparable CO2 and CO | same | 0.1 | 3 of 25 comets (HP22 Table 2); the heliocentric distance and the CO admission rule bias the share (GEN.520) (READ) | DERIVED R9 | the three weights sum to 0.93; normalised |
| k4 | CO2/H2O (PARAM) | all k | LN(0.12, 0.5 dex) | geometric mean 10.9 percent and 0.51 dex of 21 comets at r_h up to 2.5 AU (Robinson database); one common offset per system is not drawn (N=2 interstellar objects) (GEN.520) (READ, COMPUTED) | DERIVED R9 | |
| k5 | CO/H2O (PARAM) | all k | LN(0.03, 0.6 dex), range 0.005 to 0.23 | geometric mean 2.5 percent and 0.65 dex of 16 comets at r_h up to 2.5 AU (Robinson database); one common offset per system is not drawn (N=2 interstellar objects) (GEN.520) (READ, COMPUTED) | DERIVED R9 | |
| k6 | CH3OH, NH3, HCN to H2O (PARAM) | all k | LU[0.002, 0.026], LU[0.001, 0.0073], LU[0.0009, 0.0029] | Willacy 2025 cross-comet ranges (READ) | SOURCED | log-uniform shape DESIGN |
| k7 | Single-comet species: CH4 3.4e-3, C2H6 2.9e-3, H2CO 3.2e-3, H2S 1.1e-2, N2 8.9e-4, O2 3.1e-2, OCS 4.1e-4, SO2 1.3e-3, Ar 5.8e-6, Kr 4.9e-7, Xe 2.4e-7, Ne <5e-8 | all k | LN(67P value, 0.6 dex) | CH4, C2H6, H2CO and H2S take the geometric mean and spread of 9 to 18 comets at r_h up to 2.5 AU (Robinson database); N2, O2, OCS, SO2, Ar, Kr and Xe are from one comet (67P, Rubin 2019) with a design width of 0.6 dex and stay stand-in values (GEN.520) (READ, DESIGN width) | DERIVED R9 | O2 presence probability unsourced |

**Example k (COMPUTED).** Class table normalised: k1 CO2-dominated 0.536, k2 CO-dominated 0.357, k3 comparable 0.107. CO2/H2O bins (<0.03, 0.03-0.1, 0.1-0.3, >0.3): 0.114, 0.323, 0.350, 0.213.

### l. White-dwarf and pulsar planet options

| id | option | gate | weight w (symbols above) | weight source | status | notes |
|---|---|---|---|---|---|---|
| l1 | Rocky body with silicate-vapour exosphere (Na, SiO, Mg, Fe) | host WD or PSR; rocky | Phi((log10 P_eq(Ts) + 6)/0.5), P_eq from i6 formula | Kite 2016 (READ); threshold 1e-6 bar and width DESIGN | DESIGN R10 | composition via a4 |
| l2 | Airless | same | 1 - w[l1] | complement of l1 | DESIGN R10 | fallback row |
| l3 | Surviving H/He giant | host WD, GIANT | c[0,0.1] | no more than 1 percent of white dwarfs should host close-in giants (Mustill, Veras and Villaver 2014, three giants, read in the abstract and discussion, and as Veras 2016 quotes it; READ); the bracket 0 to 0.1 holds the wider orbits that survive | DESIGN R10 | composition via tables f and g |

**Example l (COMPUTED).** WD host, Ts=1800 K: P_eq = 2.39e-05 bar; w[l1] = 0.997. No sourced composition exists for pulsar-planet atmospheres.

### m. Trace and minor gas presence (binary rows: present with probability w, absent with 1-w)

| id | option | gate | weight w (symbols above) | weight source | status | notes |
|---|---|---|---|---|---|---|
| m1 | O3, O2-supported | O2>=1e-3 PAL (any host) | 1 | column LU[192, 279] DU at 100 % PAL, Sun (Kozakis fit via project doc 2ND-HAND; Cooke 2022, READ); 18 DU at 0.1 % PAL | SOURCED | M-dwarf columns larger (Segura 2005) |
| m2 | O3, abiotic without O2 | O2<1e-3 PAL; CO2-rich, H2-poor, M or K host | c[0.1,0.6] | O3 about 10 % of Earth column (Domagal-Goldman 2014, 2ND-HAND); O3>1 ppm (Gao 2015, READ) | DESIGN R11 | |
| m3 | N2O | HAB with life_chemical, or SEP-driven weakly reducing air | c[0.05,0.5] | TEST ~400 ppbv (Airapetian 2020 READ); presence unsourced | DESIGN R11 | longer lifetime around M dwarfs |
| m4 | PH3 | ROCKY with air | c[0,0.05] | TEST 20 ppb Venus claim (Greaves 2020 READ); lifetime hours (Sousa-Silva 2020, READ) | DESIGN R11 | |
| m5 | HCl | c1 to c4; crust-like inventory | c[0.1,0.5] | 1.4 % at 1000 K 100 bar, crust composition (Schaefer 2012 READ) | DESIGN R11 | |
| m6 | HF | same | c[0.1,0.5] | 4.5 % at 1000 K 1e-6 bar, crust composition (READ) | DESIGN R11 | |
| m7 | Na gas | f2, f3, h, i, BD; T in the 0.1-1 bar layer | Phi((T-962)/33) | Na2S T_cond 929 to 996 K at 0.1 to 1 bar (part 06 s2.6, COMPUTED from READ fit); Na 1.66e-6 per H times m | DERIVED | |
| m8 | K gas | same | Phi((T-775)/26) | KCl T_cond 749 to 801 K (part 06 s2.6); K 1.17e-7 per H times m | DERIVED | |
| m9 | OCS (COS) | gas has CO and H2S | log10 K = -1.477 + 4/T (CO+H2S=COS+H2) | part 05 s3.3 (COMPUTED from NASA polynomials) | DERIVED | formation equilibrium, not a probability |
| m10 | Radon-222 | ROCKY with air and solid ground | 1 | flux J = C_Ra f rho_b sqrt(lam D_e), f 0.05 to 0.7 (UNSCEAR 2008, READ); D_e unsourced | DERIVED | |
| m11 | NO2, CH3Cl | HAB with life, or lightning | c[0,0.3] | NOT REACHED (Harman 2018) | DESIGN R11 | |

**Example m (COMPUTED).** Warm giant, T(0.1 to 1 bar)=900 K, m=10: m7 Na present with p=0.030 (Na is in the cloud), m8 K present p=1.000, m9 K_eq(COS) = 0.034 at 900 K.

### n. Clouds, hazes and dust by temperature window and class

| id | option | gate | weight w (symbols above) | weight source | status | notes |
|---|---|---|---|---|---|---|
| n1 | Oxide, silicate, iron clouds: Al-Ti oxides 1850-2000 K; Fe 1690-1840; forsterite 1597-1698; enstatite 1513-1597; Cr 1416-1521; MnS 1271-1343 | GIANT, BD, i-types; profile crosses window | 1 | T_cond fits Visscher 2010 (READ), part 06 s2.6 (COMPUTED) | DERIVED | 10^4/T_cond = 5.44-0.48 logP-0.48[Fe/H] (Fe) |
| n2 | Na2S 929-996 K; KCl/ZnS 749-801 K | same | 1 | same | DERIVED | |
| n3 | NH4SH 189-223 K (0.1-10 bar); H2S ice 105-128 K only if NH3 absent | GIANT, BD, f1, f2 | 1 | same | DERIVED | |
| n4 | Water cloud or ice 161-253 K (pp 1e-9 to 1e-3 bar) | any with H2O | 1 | same | DERIVED | |
| n5 | CO2 ice 114-195 K; NH3 ice 135-179 K; CH4 ice 53-90 K; N2 ice 31-62 K, CO ice 34-66 K | cold GIANT, K, Y[P], D, comets | 1 | part 06 s2.4, s2.6 (2ND-HAND ice curves) | DERIVED | |
| n6 | Silicate-dominated aerosol | GIANT, BD, hot rocky | Phi((Teq-950)/100) | crossover Teq about 950 K, Gao 2020 (READ); width DESIGN | DESIGN R12 | |
| n7 | Photochemical haze (tholin, soot) | carbon-bearing air | (1 - w[n6]) * y; yield y = LU[0.01, 1]; monomer 1 to 50 nm, 1.0 g/cm3 | lab rates span 250x (Horst 2018 READ); haze threshold CH4/CO2 0.1 to 0.2 (Arney 2016 READ) | DESIGN R12 | near zero for carbon-free gas, lower in H2 gas |
| n8 | Mineral dust (Mars-type, r_eff 1 to 2 um) | ROCKY with air and surface wind u*>u*_t | 1 if u*>u*_t | threshold from Shao-Lu and Iversen-White formulas (READ, COMPUTED to 1.6x of Kok 2012) | DERIVED | |
| n9 | S8 / H2SO4 aerosol | outgassing S flux >=100x Earth volcanic (N2, CO2 air); Earth flux in H2 air | 1 | Hu 2013 (READ); H2SO4-H2O vapour curve NOT REACHED | DERIVED | |

**Example n (COMPUTED).** Giant with profile 300 to 2000 K and Teq=800 K: windows crossed (n1, n2): Al-Ti oxides, Fe, forsterite, enstatite, Cr, MnS, Na2S, KCl/ZnS; n3 to n5 are below 300 K and absent. Aerosol dominance: n6 silicate 0.067, n7 haze 0.933 (before the yield factor).

### o. Ground radioisotope abundance regimes by host (part 07)

| id | option | gate | weight w (symbols above) | weight source | status | notes |
|---|---|---|---|---|---|---|
| o1 | Thorium (and U) scale, thin-disc Sun-like host (PARAM) | host thin disc, FGK; rocky | log10 f_Th ~ N(0.08, 0.058), f relative to BSE | Botelho 2019: [Th/H] mean +0.080, scatter 0.058 dex, 53 solar twins (READ) | SOURCED | BSE: U 0.020, Th 0.077, K 280 ppm (McDonough 2020 READ) |
| o2 | Other populations (thick disc, halo, young, M host) | not o1 | log10 f ~ U[-0.23, +0.40] | Unterborn 2015: 59 to 251 % of solar (READ); shape DESIGN | DESIGN R13 | |
| o3 | K/U (mass) | ROCKY | LU[8.9e3, 1.4e4] | crust 8900 to 11800, BSE 1.4e4 (McDonough 2020, Wipperfurth 2020, READ); volatile-depleted bodies lower, NOT REACHED | DESIGN R13 | |
| o4 | Th/U | ROCKY | U[3.8, 4.5] | BSE 3.77, crust 4.0 to 4.5 (READ); stellar Th/U vs cosmic age NOT REACHED | DESIGN R13 | |
| o5 | Heating vs age | ROCKY | h(t)/h(4.5 Gyr) = sum s_i exp(lam_i (4.5 - t)), shares U238 .384, U235 .016, Th .402, K40 .193 | McDonough 2020 shares and half-lives (READ); COMPUTED | DERIVED | Al-26 adds 4.2e-7 W/kg for the first 20 Myr (RECALLED Al) |
| o6 | Al-26 early heating present | planetesimals formed within a few Myr | c[0.1,1] | 26Al/27Al 5.25e-5 (McDonough 2020 READ); presence in a given system unsourced; `birth_radioisotopes.py` (GEN.476) gives the share of stars born above the melting ratio: 12 percent (6 to 20) | DESIGN R13 | |

**Example o (COMPUTED).** Thin-disc Sun-like host. Bins f<0.7, 0.7 to 1.4, >1.4: 0.0000, 0.8729, 0.1271. Heating at age 3 Gyr relative to 4.5 Gyr: 1.432 (part 07 table: 1.43); at 8 Gyr 0.589. Heat per kg of rock = 98.29 U + 26.18 Th + 3.387e-3 K (uW/kg per kg/kg, McDonough 2020).

### p. Host-star UV regime (ozone and photochemistry)

| id | option | gate | weight w (symbols above) | weight source | status | notes |
|---|---|---|---|---|---|---|
| p1 | F/G, quiescent: strong NUV, FUV/NUV(280-400 nm) about 1e-3 | host F, G; age>tau_sat | 1 | band fluxes at S=1 (part 04 s2.3, COMPUTED from MUSCLES-type DATA); Sun FUV/NUV 9.8e-4; F-star FUV NOT REACHED | DERIVED | |
| p2 | K, quiescent | host K | 1 | K2V 1.9e-3 (same) | DERIVED | |
| p3 | M, post-saturation: FUV and Ly-alpha dominated | host M; age>tau_sat | 1 | M1.5 1.2e-2, M2.5 1.6e-2, M4 2.7e-2 to 3.6e-2, TRAPPIST-1 (M8) 0.60 (same); O3 lifetime 0.8 to 2 h vs 124 s for the Sun | DERIVED | long CH4, N2O, CH3Cl lifetimes (Segura 2005 READ) |
| p4 | M, saturated and flaring (L_X/L_bol 10^-3.14 for early M) | host M; age<tau_sat; tau_sat = LU[0.162, 2] Gyr | 1 | tau_sat 162 Myr early M (Jackson 2012 Table 2, READ) vs 2 Gyr project stand-in | DERIVED | proton ozone loss only if no magnetosphere (Segura 2010 READ) |
| p5 | Young star (any type), age<0.1 Gyr | age<0.1 Gyr | 1 | 1-1180 A flux 110x the 4.57 Gyr value at 0.1 Gyr (Ribas 2005, READ) | DERIVED | |
| p6 | Rotation track (L2 on p4, p5, solar mass): slow / median / fast | any | percentile of the initial rotation distribution uniform 0.02 to 0.98 | an observed distribution has uniform percentiles by construction, so the draw carries no weights of its own; the saturation time follows the percentile and the mass: Johnstone et al. 2021 Fig. 10 read at the 5th, 50th and 95th percentiles for 0.1 to 1.2 solar masses (`saturation_time.py`; matches Tu 2015 at one solar mass, 5.2, 22.8 and 297 Myr against 5.7, 23 and 226) (READ) (GEN.493) | DERIVED R14 | L_X at 100 Myr differs 200x |

**Example p (COMPUTED).** M1.5 host, planet age 1 Gyr. Jackson 2012: L_X/L_bol = 10^-3.14 (t/tau_sat)^-1.18 with tau_sat = 0.162 Gyr gives 8.47e-05 (the paper's table value is 8.5e-5). With tau_sat ~ LU[0.162, 2] Gyr, the chance that age 1 Gyr is still saturated is ln(2/1)/ln(2/0.162) = 0.276: p4 0.276, p3 0.724.

## 3. Formula quick list for the coder

dIW_class cuts: a5, a7. Snow line: a_snow = sqrt(L/Lsun) (278.6/T_snow)^2 AU. Shoreline: r = X/(K v^4), v = 11.186 sqrt(M/R) km/s, K ~ LU[3.0e-4, 6.6e-4]. Envelope loss: f_lost = eps pi R^3 E_cum S/(G M^2). Metallicity: log10 m ~ N(log10[9.7 (M/MJ)^-0.45], 0.3). Brown dwarf: L = 4e-5 Lsun (t/Gyr)^-1.3 (M/0.05 Msun)^2.64, Teff = (L/(4 pi R^2 sigma))^(1/4), R = 1 RJ. Lava vapour: log10 P_eq = k1 exp(k2 Ts) (k1 -34.7, -60.1, -60.2; k2 -0.00112, -0.00122, -0.00101 for vapour fraction <0.01, about 0.5, >0.9). Quench: part 05 s7.3. Condensation: part 05 s7.2, part 06 s2.1.


**Row counts (parsed from the tables above):** 104 rows: SOURCED 8, DERIVED 48, DESIGN 48.

## 4. No sourced weights yet

Every DESIGN row, grouped by the research that would source it. Until then the generator draws the bracket in the row's weight column with the seeded draw.

- **R1** (b4, e3, e4, e5): Occurrence-weighted distribution of bare rock / H/He envelope / water-rich worlds and of envelope metal fraction by radius, mass and insolation: radius-valley occurrence (Fulton 2017), Owen and Wu 2017 (READ in part 04), mass-radius population fits, JWST sub-Neptune retrievals.
- **R2** (b1, b2, b3, b4, b5, b6): Volatile delivery versus formation distance: population-synthesis tables of water mass fraction against a_form (Thiabaud 2014 arXiv:1312.3085 data; 2015 arXiv:1507.01343 body unread), nebula water partial pressure and T_snow versus time (sets the snow-line law), mapping to Ortenzi's Dry/Low/High classes.
- **R3** (a1, a2, a3): Magma-ocean fO2 prior against mass and the later-mantle oxidation law: Sossi 2020 data (part 08), Gaillard 2014 and 2021, Hirschmann 2012 (graphite, C solubility), degassed 40Ar fraction by tectonic regime, redox of small-body building blocks.
- **R4** (d6): Cosmic-shoreline constants and pressure caps by regime: Zahnle and Catling 2017 arXiv:1702.03386 (fetch via export.arxiv.org), Meni-Gallardo and Palle 2025, XUV histories (Johnstone 2021, part 04).
- **R5** (f7): Giant-planet C/O and C-rich host fraction: Brewer and Fischer 2016 arXiv:1608.06286 (PDF downloaded, body unread), Delgado Mena 2010 sample selection, Welbanks 2019 retrieval spread, Thorngren 2016 scatter of Zpl/Zstar.
- **R6** (h2, h3, h4): Hot-Jupiter thermochemistry above 2200 K (FastChem or GGchem tables), Venot 2012 arXiv:1208.0560 body, a Kzz distribution, edges of the 1500 and 2200 K regimes.
- **R7** (i6): Lava-world vapour: Miguel 2011 numeric Type I to V tables, Schaefer and Fegley 2009 arXiv:0906.1204, LavAtmos arXiv:2210.10463, and the fraction of lava worlds that keep volatiles (JWST 55 Cnc e and other targets).
- **R8** (j1, j2, j3, j4, j5): Evolved and frost atmospheres: Fray and Schmitt 2009 and Grundy 2024 vapour pressures (N2, CO, CH4, CO2), Kopparapu c and d coefficients and the 5 and 0.1 Earth-mass sets, occurrence of Titan-, Pluto-, Venus-like bodies, abiotic-O2 frequency (Luger and Barnes 2015 arXiv:1411.7412 body).
- **R9** (k1, k2, k3, k4, k5, k7): Comet coma classes and ratios: Harrington Pinto 2022 full per-comet table (the class fractions here are a SNIPPET), Mumma and Charnley 2011, Bockelee-Morvan 2004, Dello Russo 2016, A'Hearn 2012, and a second comet with O2, noble gases, N2.
- **R10** (a4, l1, l2, l3): White-dwarf and pulsar planets: Veras 2016 arXiv:1601.05419 body, polluted-WD debris compositions beyond Doyle's six, atmospheres of pulsar planets (none found).
- **R11** (m2, m3, m4, m5, m6, m11): Trace gases: Harman 2018 (ApJ 866, 56) NOx and N2O, Domagal-Goldman 2014 numeric abiotic O3, Bains 2021 PH3, other Schaefer 2012 inventories (HCl, HF by composition), lightning yields.
- **R12** (n6, n7): Clouds and hazes: Gao 2020 width of the 950 K crossover, Horst 2018 yield against T and C/O, H2SO4-H2O vapour curve, perchlorate and carbonate decomposition equilibria.
- **R13** (o2, o3, o4, o6): Radioisotopes: Th, U, K against [Fe/H] and Galactic population (Frebel and Hill cosmochronometers), K volatile depletion by body class, White and Klein 2014 oceanic crust, 26Al abundance spread between systems.
- **R14** (p6): Initial-rotation distribution to weight Tu 2015 tracks, and tau_sat against stellar mass for M dwarfs (MUSCLES-type samples).
- **R15** (g1, g4): Brown-dwarf type-Teff relation and the T/Y boundary near 500 to 700 K from an observational compilation; hot young end above 2200 K.

## 5. Conflicts between the repo and part 05 (and inside part 05)

1. **Class letters.** `PLANET_CLASSES` has 20 classes (A B C D E F G H I J K L M N O P Q S T V). R and W were removed (comment at `tuning.py` line 1601, which also names four removed classes without letters); U, X, Y, Z are absent. The design doc plans R, U, W, X, Y, Z and retiring T. Part 05's columns are regimes, not letters; mapping used here: RR and RO = redox rows a5 to a7 on ROCKY; LV = A, B, S, C in zone h (table i); SN = I, T, V and R[P]; WW = O and X[P]; GG = I, J, T; HJ = J zone h; CM = COMET and icy AST.
2. **I and T overlap.** I (15,000 to 50,000 km, 2.4 to 7.8 R_E) and T (15,000 to 55,000 km) cover the same sizes; part 05's sub-Neptune regime (1.6 to 3.2 R_E) is only the low end of I, and V (10,000 to 15,000 km, 1.6 to 2.4 R_E) straddles the radius valley while the repo calls it rocky. The design doc's rocky test is R<=1.12 M^0.27.
3. **NOAIR classes.** C, D and S have `atmosphere` None in every zone; `CLASS_MIXES` has no key for them, so a pressure would fall to DEFAULT_MIX (N2 .78, O2 .2). Table d can give them trace or tenuous air (part 04 and the design doc say a 2 to 10 Earth-mass rock in e or c is far inside the shoreline).
4. **Free O2 without life.** CLASS_MIXES has O2 0.10 (G), 0.25 (H), 0.30 (L), 0.22 (O, P), 0.15 (Q), 0.08 (V). Abiotic O2 needs water loss around M or K stars (part 04 s3.6); PLANET_CLASS_MAX_LIFE_STAGE caps G at photosynthesis but no rule ties O2 to life stage. L has Ar 0.60, far above the radiogenic budget (Earth: Ar 0.93 % with about 43 % of 40Ar degassed, part 07 s2.10).
5. **E and F.** Text and mixes (E: H2O .55, CH4 .25, NH3 balance; F: CO2 .60, CH4 .20, NH3 balance) put NH3 and CH4 in 340 to 400 K air; NH3 appears in part 05 only at 800 K and 100 bar for dIW<=1 and is photolysed fast (part 04); redox rows a5 to a7 should select these mixes.
6. **N.** Text says "dense, reducing mix of CO2 and sulfides"; CLASS_MIXES is oxidised (CO2 .965, SO2 1.5e-4, H2S 1e-4). Sulfur is H2S only below dIW about +2.5 (part 05 s3.4). There is no hot-zone Venus-like letter (N is e only; A, B are small).
7. **I, J, T mixes.** One CH4 3e-3, H2O 1e-3, H2S 1e-4 (T: CH4 .02) for every zone. Part 05 s7.4: CH4 holds carbon only below about 500 K; J in zone h is a hot Jupiter (CO, tables h); H2O and H2S are cloud-depleted in cold giants.
8. **Redox.** Repo mean dIW 3.5 + 2 log10 M (a2) against magma-ocean +0.5 at 1 Earth mass (a1, Sossi); REDOX_OXIDIZED_DELTA_IW=2.5 equals the a7 cut (consistent). Meteoritic dIW is -1.4 to 0.
9. **Brown dwarfs** have no class letter (planet_class NULL; rogue-eligible letters are C D I J S T only); their type must come from Teff (table g). Part 05 s7.4 leaves 500 to 700 K uncovered (row g4 added).
10. **Part 05 C/O recipe.** C/O = 0.55 * 10^N(0, 0.15) truncated to [0.2, 1.5] puts 13.9 % above 0.8 (COMPUTED), not the intended 2 to 5 %; row f7 truncates the core at 0.8 and adds an explicit tail.
11. **Part 05 s4 lava boundaries.** It states the I/II boundary differs by "up to 130 K (2033 vs 2096)"; the difference is 63 K, and the caption list has five numbers for four boundaries (2735 vs 2890 differ by 155 K). Rows i1 to i5 use the text edges with width 60 K.
12. **HABITABLE_PLANET_CLASSES** (E F G H K L M O P Q V) is a life flag, not an atmosphere regime; N has no life_chemical and is excluded; Q is zone e only. Moons draw any letter except Q, S, V, so Titan-like moons have no letter before Y.
13. **Shoreline status.** The project's K window and r thresholds are fitted to 13 Solar System bodies ([C] in the design doc) and the cited paper is SNIPPET only, so table d is DERIVED from a calibration, not a published fit.

## 6. Recalled / to verify

RECALLED or SNIPPET values used: T_crit of water 647 K (c5, c6); the K window 3.0e-4 to 6.6e-4 and r thresholds (design doc fit, Zahnle and Catling 2017 SNIPPET); 278.6 K blackbody constant (design doc; 278.3 K for S0=1361); Kopparapu runaway limit 1.04 to 1.11 (design doc [R]); Harrington Pinto class fractions and median CO2/H2O 0.12 (SNIPPET); width 0.3 dex on metallicity, C/O 0.55 and 0.15 dex, dIW N(-1, 0.5) for WD debris (part 05 s15); all edge widths (100 K, 60 K, 50 K, 0.5 dIW) and priors c[.] are DESIGN choices, not measurements; Al 2.35 wt% for 26Al heating (part 07).

## 7. Cross-links

Part 05 (recipes restated here), part 06 s2.1 to 2.11 (rows n1 to n9, i), part 07 s2.9 to 2.10 and F10 to F12 (rows o), part 04 s2.3, s3, s4.5 (rows p, m1, m2, d). Repo: `src/planetgen/physics/atmosphere.py` CLASS_MIXES and REDOX_*, `tuning.py` PLANET_CLASSES, ROGUE_* constants, COMET_* constants.
