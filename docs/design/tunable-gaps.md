# Tunable gaps: can the unsourced composition values be seeded random draws that stay tunable and scientifically accurate?

**Status (2026-10-17):** result of the Tunable gaps investigation (a read-only thread, finished 3:36 AM EDT as the coordinator reported; Boss's question at 3:00 AM EDT), merged here by the TODO thread. It feeds the atmosphere tables GEN.478 and the loader GEN.472, the research items GEN.480 to GEN.494, and the new harness, selection models and builders GEN.503 to GEN.525. The six defaults of section 6 (D1 to D6) stand as Decided unless Boss says otherwise (composition-model.md section 4a item 11). The worker reports are in [tunable-gaps-sources/](tunable-gaps-sources/group1-subneptunes.md) and keep their own labels.


Date: 2026-10-17, written after Boss's 3:02 AM EDT usability rule (lane clocks drift). Read-only investigation, no repo or docs edits. Requested by Boss (question typed at 3:00 AM EDT): "We may be reaching the limits of what we can find in the public domain. Do we have enough information to draw the remaining composition values from random ranges and still have them tunable and scientifically accurate to what we know?"

Method: the two composition masters, part 09 and `docs/design/composition-model.md` were read by this lane. Six workers, one per data group, then checked the public record paper by paper (about 330 sources opened; arXiv ids downloaded and read, Zenodo and Crossref reached, ADS, Science, OUP, ScienceDirect and A&A refused). Their six reports are the detail behind this file and are listed in section 8. This lane read all six in full and re-ran three of their statistics (the 22-value white-dwarf redox statistics, the Wilson intervals for the giant C/O tail, the shoreline K needed for 55 Cnc e); all three reproduced. It did not re-open the papers behind every number, so each number keeps its worker label (READ, 2ND-HAND, DATA, COMPUTED, SNIPPET, RECALLED, NOT REACHED).

## 1. Answer

**Short answer: no, we are not at the limit of the public domain for most rows. For a minority we are, and for those the limit is not locked papers but data that does not exist.** Of the 48 atmosphere rows with no sourced weight:

| Class | Meaning | Rows | Count |
|---|---|---|---|
| **T** tunable draw | a published ensemble plus a physical bracket defines a seeded draw whose knobs can be tuned to a stated population statistic | a3, a4, p6, k4, k5 | 5 |
| **F** few constraints | tunable, but only to a few stated constraints (N small, or method-dependent) | b1, b2, b3, e5, f7, h4, i6, k1, k2, k3, n6, o2, o3, o6 | 14 |
| **P** derive | not a probability: a deterministic map or a physics outcome. Compute it; the only randomness is upstream (inventory, redox, mass, age, temperature) | b6, d6, e3, e4, g1, g4, h2, h3, j1, j5, l1, l2, l3, m2, m5, m6, o4 | 17 |
| **S** needs a source | the public record has no data set that defines the draw; the exact paper or data set is named in section 6 | a1, a2, b4, b5, n7, k7 | 6 |
| **X** design flag | no ensemble can exist (no exoplanet measurement, or the quantity is whether life makes a gas); keep a labelled bracket, default off or flagged, never tune | j2, j3, j4, m3, m4, m11 | 6 |

So 36 of 48 (T, F, P) can be handled from the public record, 17 of them by removing the draw. The 12 rows in S and X are the real limit. Counts are this lane's classification of the workers' verdicts (COMPUTED tally; the per-row reasoning is in section 3 and the worker files).

Three findings matter more than the counts:

1. **Tuning would make the generator more accurate, not less.** The public record contradicts four current brackets: the cosmic-shoreline K window (fitted to 13 Solar System bodies; two verified exoplanet papers disagree, and 55 Cnc e and TOI-561 b need K 12 to 190 times above it), the white-dwarf redox centre (a4 uses -1.0 from six stars; 22 values give -1.25 +- 0.10), the giant C/O tail (f_tail upper end 0.34 comes from a result no larger sample reproduces), and the M-dwarf saturation time (p4 mixes masses in one range).
2. **Many "weights" are really physics.** Brown-dwarf type, the hot-Jupiter regime edges, the lava-world and white-dwarf exosphere gates, Venus-like and abiotic-O2 outcomes, ozone presence and halogen gases are outputs of Teff, chemistry, escape or inventory. Drawing them as weights hides the physics and cannot be tuned; deriving them removes 17 unsourced rows at once.
3. **Where tuning is possible it needs a selection model on the simulated side.** Five shared modules are enough: Kepler transit-plus-recovery, a mass-radius sample cut, a JWST atmosphere target selection, the comet admission rule by heliocentric distance, and a white-dwarf pollution detection model (section 4). Without them a tuned width is only an upper bound.

**Where the Solar System alone constrains a knob** (these break NO FIXED CONTENT if tuned to; section 5): the mass slope of mantle redox, the later-mantle law, the delivery payloads, K depletion and K/U, the Th/U formation duration, the Titan, Pluto, Mars and Venus-PH3 rows, the within-system comet scatter, the belt taxon mix, KBO densities, moon origin and the giant bulk-to-atmosphere factor.

## 2. Boss's rule applied (3:02 AM EDT: verified and partly verified usable; verified wins on contradiction)

Partly verified material was counted as available input. Places where a verified source overrides a partly verified, snippet or recalled value (all from the workers' reports, each READ unless stated):

| Row | Partly verified or recalled value | Verified value that wins | Source |
|---|---|---|---|
| d1 to d5 | K window LU[3.0e-4, 6.6e-4], slope 4, "Meni-Gallardo keep slope 4" (snippet) | Meni-Gallardo and Palle slope 5.77, intercept -4.35; Berta-Thompson logistic: log10 f0 = 2.62 (+0.45/-0.31), p = 5.9 (+0.61/-0.43), q = 1.17 (+0.28/-0.20), width 0.247 dex; Zahnle and Catling draw the line "by eye" and give no K | arXiv:2508.12865, 2507.02136, 1702.03386 |
| d (TRAPPIST-1) | "b and c bare or under 0.1 bar" (snippet) | at or above 1 bar disfavoured; 0.1 bar is in no paper read | arXiv:2303.14849, 2306.10150 |
| a4, a1 | WD debris centre -1.0 (6 stars, RECALLED); a1 floor -2 (Solar System Moon theory) | 22 values from Doyle 2020 and 2023 Tables: mean -1.25 (SE 0.10), median -1.26, sd 0.45, intrinsic about 0.22 (COMPUTED by worker, re-run by this lane) | arXiv:2008.05992, 2303.00063 |
| e1/e2 | eps = 0.05 to 0.2 fixed (2ND-HAND) | Owen and Wu variable efficiency eta = 0.1 (v_esc/15 km/s)^-2 | arXiv:1705.10810 |
| e5 | K2-18 b "near solar to tens of x" | water 0.6 to 162 times solar | Benneke 2019 (arXiv:1909.04642) |
| f7 | Delgado Mena 34 percent of hosts above C/O 0.8 | Brewer and Fischer (N=849) and Suarez-Andres (N=499): no star above 0.8; one blended oxygen line explains the 34 percent (Teske 2013). Both ends of the design tail bracket [0.0013, 0.34] are outside the planet data (0.043 to 0.087) | arXiv:1608.06286, 1801.09474, 1309.6032 |
| g4 | "T/Y gap 500 to 700 K" | that range is T8 to T9; Y0 at 371 to 510 K (Kirkpatrick 2021 Table 13; Beiler 2024) | arXiv:2011.11616, 2407.08518 |
| h2/h3 | edges at 1500 K and 2200 K | computed: CH4=CO at 966 to 1907 K (0.1 to 100 bar), half the water dissociated at 2570 to 3600 K (1 mbar to 1 bar) | worker Gibbs calculation, NASA-format polynomials |
| h4 | TiO above 2000 to 2300 K | Ti-bearing gas kept only at Teq 2641 K (WASP-189 b); absent at 2210 to 2350 K; VO kept at 2210 K | arXiv:2305.17228, 2111.12732, 2210.12847, 2306.08739 |
| k1 to k3 | weights 0.5, 1/3, 0.1 (rounded abstract words) | counts 13, 9, 3 of 25 | Harrington Pinto 2022 Table 2 (arXiv:2209.09985) |
| k4, k5, k7 | log-normal widths "by eye", 67P-centred values | N=21 (10.9 percent, 0.51 dex), N=16 (2.5 percent, 0.65 dex); population geo-means for CH4, C2H6, H2CO, H2S, OCS | Robinson 2024 database (CC BY 4.0) |
| m2 | O3 "about 10 percent of Earth" | 0.01 to 0.1 percent for G, K and M hosts, 15 percent for an F host | Domagal-Goldman 2014 Table 3 (arXiv:1407.2622) |
| p4 | tau_sat = LU[0.162, 2] Gyr | 0.16 Gyr (0.5 to 0.6 Msun), 0.6 +- 0.3 (0.2 to 0.3), 2.3 +- 1.3 (0.1 to 0.2) | arXiv:1111.0031, 2205.02331 |
| o2 | Hypatia thin-disc Th scatter 0.14 to 0.19 dex | Botelho solar twins 0.056 dex | arXiv:1810.10413 |
| T19 | rock fraction 0.58 to 0.82 | dwarf planets 0.74 to 0.82, small KBOs 0.3 to 0.45 | arXiv:1603.06224, 1311.0553 |
| R2 (part 09) | Thiabaud 2014 as a water-against-distance source | that paper is about refractories; mentions water five times | arXiv:1312.3085 |

Verified against verified, no override (spread reported, not resolved): Putirka and Xu (12 of 23 white-dwarf bodies "exotic") against Trierweiler (most chondritic); Sing 2024 (WASP-107 b Kzz 10^11.6) against Welbanks 2024 (10^8.4 to 10^9.0); Sasaki 2010 against Cilibrasi 2018 (icy-moon fraction); Rogers 2025 (water mass fraction at most 3 percent) against Chakrabarty 2026 (59 percent water-rich); Luque and Palle (N=34, three populations) against Parc 2024 (N=46) and Dainese 2025 (N=60) (no water-world population).

## 3. Row by row (atmosphere, part 09)

Notation: knob -> statistic -> target (uncertainty, source label). "Test" is what the consistency test compares. Details, N, access and licences are in the worker files.

### 3.1 Delivery, water worlds, envelopes (rows b, e; body T5, T9, T11)

| Row | Class | Finding | Knobs -> statistic -> target | Test | Solar System flag |
|---|---|---|---|---|---|
| e1/e2 bare rock vs retained H/He | T (conditional) | CKS (970 planets), public Kepler DR25 injections and a fitted recovery curve make an injection-recovery test possible. The two part-09 knobs are degenerate with the core-mass distribution and with the mechanism (photoevaporation against core-powered loss, favoured only at 1 to 2 sigma, Berta... Berger 2023). A retained-envelope radius law is missing; a bare/retained binary cannot give the partly filled gap | eta(v_esc) index, f_env0 centre and width, core-mass peak -> valley slope in period -> -0.11 +- 0.02 (READ, Petigura 2022, arXiv:2201.10020); valley location 1.84 +- 0.03 Earth radii at 10 d; core mass peak 3.72 (+0.45/-0.33) Earth masses (Rogers and Owen 2021, arXiv:2007.11006); valley depth by stellar-mass bin (Ho 2024) | simulate the generator's hosts, apply transit probability and the recovery curve, compare the valley with the published 1 sigma in three stellar-mass bins; one gap-finder method on both sides. Cannot choose a mechanism (needs about 5000 planets) | none |
| e-link (stellar mass, metallicity, age) | F | mass slope 0.18 to 0.26 (three papers); metallicity null (0.01 +0.05/-0.06); age only through one ratio (0.61 to 1.00) | slope -> dlogRp/dlogM* -> 0.18 (+0.08/-0.07) (Petigura) | same harness | none |
| b4 water-rich share | S | published shares run from at most 3 percent to 87 percent depending on the formation model; mass and radius alone cannot break the degeneracy with H/He (about 6 informative JWST spectra) | none until a formation model is chosen | density distribution of precise M-R planets (LP22, Parc and Dainese disagree) | LP22 uses a Solar System 1:1 ice-rock ratio; use the generator's own snow-line law |
| e3/e4 water world vs hycean | P | hycean is a gate: Teq below about 400 K, H2 at most 0.1 wt percent, water at least 10 wt percent; K2-18 b is contested (Madhusudhan 2023, Wogan 2024, Nixon 2025) | none (e4 = p[b4] x gate, e3 the remainder) | gate truth table | none |
| e5 envelope Z | F | about 8 JWST retrievals plus one 15-planet exo-Neptune fit near 100x solar; atmospheric Z is not bulk envelope Z; the 0.6 upper edge (about 106x solar) excludes K2-18 b (162x), GJ 1214 b and GJ 3090 b | bracket edges only; no population mean exists | per-planet retrievals at their own 1 sigma; a trend test needs about 45 JWST planets | the only mass trend (Madhusudhan 2025) is Solar System giants; do not use |
| b1, b2, b3 dry/low/high | F (MODEL label) | the Ortenzi payloads (0, 0.05, 0.5 wt percent H2O) are Earth delivery scenarios. Formation distance is not observable. The only function of water against formation distance is population synthesis; its output is public (Bern NGPPS, Zenodo 10527170, open, no licence field). The worker extracted a migrant probability: logistic in mass, centre 6.1 Earth masses, slope 6.3 per dex (COMPUTED) | p_mig centre and slope -> fraction of water worlds by mass -> NGPPS itself, checked against Luque and Palle (water worlds start at 2 to 3 Earth masses; rocky to about 10) | density modes 0.94, 0.47, 0.24 with the LP22 selection (mass better than 25 percent, radius better than 8 percent) | payloads and T11 range are Earth-only; flagged |
| b5, b6 | S, P | no ensemble measures sub-percent water; rogue water is not a separate draw (the ejected planet keeps what it formed with) | n/a | n/a | b6 defaults have no data |
| T9 migrant | F | NGPPS logistic above; minimum ice-bearing mass 3 (M*/Msun)^0.48 Earth masses from two NGPPS-IV anchors (figure reads) | see b1 | see b1 | none |
| T11 delivery inside the snow line | S | range is arithmetic from Earth's water budget; Raymond and Izidoro 2017 is a Solar System model | none | none | yes, Earth only |
| T5 iron branch | F (upper limit) | super-Mercury share: 0/20, 1/21, 1/11, 5/22 by sample (95 percent intervals 0.001 to 0.24 for 1/21); a mock test reproduces the outliers about 15 percent of the time from noise alone; no coreless count | p_SM prior from 1/21 -> number of 2-sigma outliers -> 1 of 21; slope planet-to-star iron 0.94 (+1.02/-1.07) (Plotnykov 2025, arXiv:2511.17717) | mock population with the sample's mass and radius errors | none (Mercury is only a label) |

### 3.2 Redox, radioisotopes, ground (rows a, o; body T8, T20, T22, T23, R2, R4)

The best new find of this investigation is a real extrasolar ensemble for the small-body end of mantle redox: **24 polluted white dwarfs** with O, Mg, Si, Fe measured (Doyle 2020 Table 1, Doyle 2023 Table 4), 22 values and 2 upper limits.

| Row | Class | Finding | Knobs -> statistic -> target | Test | Solar System flag |
|---|---|---|---|---|---|
| a3 building blocks | T | the uniform U[-1.4, 0] has mean -0.7, 0.55 dex above the WD mean | centre and width -> median, sd of 22 WD values -> -1.26 (SE 0.10), sd 0.45 (intrinsic about 0.22 from 13 well-constrained stars) | median, sd and a KS test against the 22 values, with the per-star errors (0.07 to 0.43), the detection truncation near -0.9, no lower bound for 7 of 24 and the DB settling bias (Teff 16000 to 19000 K) | none; the WD scale (core-formation dIW, Earth about -1.5) differs from the surface magma-ocean scale (Earth +0.5 on Sossi's); the project must say which each row draws (consistent with GEN.498) |
| a4 WD debris | T | same ensemble; N(-1, 0.5) from six stars | same | same | none |
| a1 magma-ocean fO2 against mass | S | low-mass end tunable to the WDs; the Earth end (+0.5) is Earth only; the mass slope has no observational anchor beyond theory (Deng 2020 uses three Solar System bodies, base depths 55, 14, 5 GPa; Guimond 2023: at least 2 dex of mineralogical spread, "real fO2 unconstrained for exoplanets"). The 1.25 dex-per-dex slope was interpolated between a floor and Earth | derive surface fO2 from ocean base pressure; draw only the scatter (Deng 0.5, Guimond 0.3 to 0.4 dex) | none possible | yes, Earth end and slope |
| a2 later mantle | S | law fitted to Moon, Mars, Earth | none | none | yes, only the Solar System |
| T8 | S | same as a1 and a2; no per-planet distribution exists | as a1 | none | yes |
| o2 Th scale, non-thin-disc | F | Hypatia (14,085 stars; Th II in 284, Eu in 4570, K in 1509, no U): [Th/Mg] thin -0.044 (sd 0.153, N=215), thick -0.063 (sd 0.193, N=32), no offset; halo, young and M hosts have no Th data | log f_Th ~ N(m_pop, s) with s = 0.056 (Botelho, intrinsic) to 0.19 (Hypatia, mixed catalogues) | [Th/Mg] by disc population | none |
| o3 K/U | F via Na proxy, else S | stellar K exists (N=1509, [K/Mg] sd 0.22 dex) but planetary K is volatile-depleted; only Solar System bodies give depletion (Moon 0.032, Earth 0.18, Mercury 0.18, Mars 0.28). Extrasolar: Na in 202 DZ white dwarfs spans formation temperatures below 1000 K to above 1400 K | d_K = f(T_form) with T_form tuned to the Na/Mg spread (Harrison 2021, arXiv:2103.05713) if Boss accepts Na as a proxy | WD Na/Mg histogram | yes: Solar System K/U must not be the tuning target |
| o4 Th/U | P | derivable: Th/U = PR x F(T) x exp[(l238 - l232) tau]; PR = 1.67 to 2.0 (READ via Frebel 2007); table at tau = 4.5 Gyr gives 3.19 to 4.72 for T = 3 to 13 Gyr. T is calibrated on CI chondrites (3.76 needs T = 6.5 Gyr); no U in any disc star | PR range; T by birth epoch (needs a chemical-evolution model); tau = planet age | halo-star chronometers (about 8 stars with U) and one Solar System datum | the T calibration is Solar System only |
| o6 Al-26 present | F | range 5.25e-6 to 5.25e-4 is 0.1 to 10 times the Solar System value. Extrasolar: ISM 3.0e-6; clouds 2 to 3e-5; cluster models enrich 9 to 28 percent of disks to 3.5e-5 to 2.4e-4 (Lichtenberg 2016); WD bodies need at least 3e-5 in at least 4 percent, likely tens of percent (Curry 2022). Evidence points to a mixture, not a log-uniform | enriched flag p_enr 0.04 to 0.66 (centre 0.25 to 0.45), level of enriched about 3e-5 to 1e-4, un-enriched 3e-6 -> fraction of differentiated WD pollutants -> 66 +4/-6 percent (Bonsor 2020) | WD differentiation fraction with a body-size model | only 5.2e-5 itself is Solar System |
| T20 ground type, T22 end-members, R2 | S (derive) | WD debris is disrupted bulk, not surfaces: it gives bulk mantle mineralogy, core-rich about 4 percent of 202 DZ systems (3 sigma), crust-bearing about 0.5 to 1 percent. Putirka and Xu (23 WDs, 12 "exotic") and Trierweiler (31 WDs, 15 chondritic on raw abundances) disagree; the first one's table is "on request" only | derive ground types by petrology from the stored mantle, volatile inventory, heat flux and tectonic state; test only the bulk-mantle assemblage and the crust-bearing frequency | normative mineralogy of simulated bodies against the 23 and 31 WDs (Ca selection, Teff, settling) | the Earth rock-type library stays a flagged stand-in |

### 3.3 Retention, lava worlds, evolved atmospheres (rows d, i, j, l, m, p)

| Row | Class | Finding | Knobs -> statistic -> target | Test | Solar System flag |
|---|---|---|---|---|---|
| d1 to d5 shoreline | F | the tunable object is a two-class probability of an atmosphere, not a five-way split (nothing measures d2 to d5). Berta-Thompson: 917 labelled planets (869 atmosphere-labelled and 48 airless), 15 small planets with a direct label (2 positive, 13 negative). Atmospheres disfavoured for Earth-size HZ planets below log10 L* = -2.23 (about M4V) | log10 f0, p, q, w -> logistic classification -> 2.62, 5.9, 1.17, 0.247 dex with the intervals above | log-likelihood and fitted (f0, p, q, w) within the published intervals, with the same labelling rule, ages above 750 Myr, Teq below 1700 K, the same detection threshold | the repo K window rests on 13 Solar System bodies; Titan forces K above 2.26e-4, Ganymede caps it at 6.56e-4; this is the largest single violation. The Solar System stays legitimate as test points inside the joint fit |
| d6 caps | P | only upper limits (about 1e-5 to 10 bar by composition) and one 0.01 to 100 bar detection (55 Cnc e); derive P = min(supply, vapour, escape) | none | none | caps d4 and d5 values are Mars, Titan, Venus, Earth only |
| i6 lava worlds | F | 2 of 2 molten-dayside USPs with an atmosphere constraint are volatile-bearing (55 Cnc e CO2 or CO, 0.01 to 100 bar, rock-vapour-only ruled out; TOI-561 b); none has a bare rock-vapour atmosphere. Beta(3,1): median 0.79, 95 percent 0.29 to 0.99 | f_vol -> fraction of lava worlds with a CO2 or CO signal -> 2 of N_molten (2 to 4) | the Beta interval | the Kite vapour-pressure fit assumes bulk silicate Earth (vapour pressure only) |
| j1 Venus-like | P (F check) | gate from Kopparapu: S_eff 1.0512 (Sun), 0.99, 1.107, 1.188 at 0.1, 1, 5 Earth masses; 3-D value about 1.10 (Leconte). TRAPPIST-1 b and c (N=2) allow bare rock | w_retain residual in [0, 1] | must not forbid 0 of 2 | none |
| j5 abiotic O2 | P | outcome of water loss, EUV and oxygen sinks (Luger and Barnes: 120 to 2400 bar per ocean; Harman 2018 about zero for wet HZ planets). The O2-sink fraction is unsourced | sink fraction in [0, 1] | none | none |
| j2, j3, j4 | X | no exoplanet ensemble (Mars-, Titan-, Pluto-like); no confirmed exomoon (Nature, 22 July 2026) | none; derive from inventory (b), escape (d) and vapour pressure, leave the leftover prior [0, 1] flagged | none | yes, Solar System only |
| l1, l2 | P | the polluted-WD fraction (27 to about 50 percent, N=85, Koester 2014, arXiv:1404.2617) gates "rocky debris exists", not "keeps a silicate-vapour exosphere"; l1 is P_eq(Ts) against an exosphere criterion | threshold log10 P_eq = -6 +- 0.5 is a design constant with no data | none | none |
| l3 surviving giant | F | upper limits (under 1.5 percent hot Jupiters in 1148 K2 WDs; 0.8 +0.6/-0.4 percent transit fraction in 313 polluted WDs) and two detections (WD 1856+534 b; MOA-2010-BLG-477Lb); a prediction of over half of WDs hosting Jovian planets is 2ND-HAND. Derive from the engulfment rule (Veras 2016) x the generator's own giant occurrence | engulfment radius 1.2 x R_max(M_prog) | transit limits (period below 27 d only) | none. Pulsar-planet atmospheres: no data, derive from escape (Patruno and Kama) |
| m2 abiotic O3 | P | a photochemistry look-up by star type, CO2, H2 flux, desiccation and lightning | none | none | none |
| m3, m4, m11 | X | no measurement of N2O, PH3, NO2 or CH3Cl on any exoplanet; the PH3 Venus claim is contested by two READ papers; N2O is a biosignature flag | none; default bracket [0, 1] flagged (m4 [0, 0.05]); NO2 from lightning | none | yes: tuning to Earth N2O 400 ppbv or Venus 20 ppb breaks the rule |
| m5, m6 HCl, HF | P | HCl and HF are speciation of the crust halogen budget (percent only for continental crust; below 0.1 percent for bulk silicate Earth); crust enrichment 26.5 (F) and 11.8 (Cl) is Earth-only (COMPUTED by worker from Schaefer 2012 Table 2). Stellar F and Cl: solar log eps 4.40 and 5.31; 16 stars with Cl, 49 K giants with F | depletion exponent (unsourced), crust enrichment | none | yes, planet-side depletion |
| p4 tau_sat | T for early M, F for mid and late M | Jackson 2012: 162 Myr, systematic range 79 to 332 Myr (N=55 of 717); Medina 2022: 0.6 +- 0.3 and 2.3 +- 1.3 Gyr (N=219, kinematic ages that assume constant star formation). No X-ray tau_sat in the record for 0.3 to 0.5 Msun | tau_sat(M) piecewise in log mass, log-normal with the published uncertainty; or compute from Rossby saturation (Ro_sat = 0.13 to 0.16 with Wright's times; 0.0605 with Spada's: never mix) | fraction saturated against age and mass: 2.4 +- 0.3 Gyr, 2.0 +- 1.2 Gyr for Prot below 10 d, active fractions 74/19/8 percent | none |
| p6 rotation tracks | T | replace 0.25/0.5/0.25 by a uniform percentile u in [0.02, 0.98] on the Johnstone 2021 tracks (MIT code, masses 0.1 to 1.25 Msun) | percentile u -> cluster rotation percentiles at 150 Myr | simulated X-ray distribution against cluster data within the published scatter | none |

### 3.4 Giants, brown dwarfs, hot Jupiters, clouds (rows f, g, h, n; body T12)

| Row | Class | Finding | Knobs -> statistic -> target | Test | Solar System flag |
|---|---|---|---|---|---|
| f5/f6 giant metallicity | T (bulk), F (atmosphere) | 0.26 dex is the intrinsic scatter of the Mz fit (47 giants, factor 1.82, Thorngren 2016, arXiv:1511.07854). Chachan 2025 (147 warm giants): Mcore 14.7 (+1.8/-1.6) Earth masses, fZ 0.09 +- 0.01, scatter 0.66 +- 0.08 times MZ; it is 1.4 to 1.6 times lower than Thorngren from 0.3 to 3 MJ, and sub-Saturns have high Z (the part-09 note says the opposite). Atmospheric m (LExACoM, 69 planets, public CSV, no Solar System member): residual 0.4 to 0.9 dex, steeper than the bulk exponent | bulk centre and width as published; atmospheric slope -0.3 to -1.3 and width 0.4 to 0.9 dex from LExACoM | median Zp/Z* in mass bins (29.3, 5.9, 3.6 at 0.1, 1, 10 MJ); atmospheric residual sd within 0.2 dex | the bulk-to-atmosphere factor f = 0.2 rests on Jupiter and Saturn only |
| f7 giant C/O | T (stars), F (planets) | stars: median 0.47, sd 0.07, maximum 0.66 (N=849), 0.45 to 0.50 (N=499); no star above 0.8 (upper limit 0.35 percent at 95 percent, COMPUTED). Planets (LExACoM): median 0.57; 6 of 69 above 0.8 (Wilson 0.040 to 0.177, re-run by this lane), 3 of 69 if the 1-sigma lower bound must exceed 0.8; method-dependent with a 0.2 to 0.3 geometry systematic; the offset to the host star has no single distribution | stellar centre and width as published; planet core centre 0.57 and tail 0.043 to 0.087 | binomial counts in the five bins, N=69; host C/O against the stellar catalogue by [Fe/H] | the core centre 0.55 and the solar 0.589 are Solar System anchors; LExACoM replaces them |
| T12 diluted or compact core | S | exoplanet bulk Z cannot tell the two apart (two papers READ); structure evidence is Juno and Cassini gravity of Jupiter and Saturn | keep Mz continuous; core fraction a nuisance between all-envelope and 10 Earth masses in the core | none | yes |
| g1, g4 brown-dwarf type | P | type is a function of Teff with scatter (Kirkpatrick 2021, N=525, Table 13: RMS 134, 79, 79 K; plateau of only -18 K per subtype from L8.75 to T4.75). Draw mass from dN/dM ~ M^-0.6 +- 0.1 and age, evolve, then map | alpha 0.6 +- 0.1; Teff scatter as Gaussian at fixed type; gravity offset up to -300 K below 500 Myr | 20 pc counts per Teff bin over 600 to 2250 K (census complete there), Poisson tolerance | none |
| h2, h3 hot-Jupiter regimes | P | edges are pressure-dependent outputs of equilibrium chemistry (computed table in 2.4 of the worker file). Licence-compatible: easyCHEM (MIT), chemicals and thermo (MIT); pyFastChem and exo_k are GPL-3 and the repo is CC0, so they cannot be linked | physical inputs only (P-T, Z, C/O) | fractions in H2O-rich, CO+H2O and dissociated regimes against retrieved water detections (Tsiaras 2018: water in 16 of 30 with significant atmospheres) | none |
| h4 TiO and VO | F (a Teq threshold), S (Kzz) | no population measures Kzz, and two papers on one planet differ by 2.6 to 3.2 dex. The empirical constraint is a Ti threshold between 2350 and 2640 K from N=3 to 6 | T_trap -> TiO detection against Teq -> yes at 2641 K, no at 2210, 2250, 2260, 2350 K | a step at T_trap against that list | none |
| n6 silicate against haze | F | Gao 2020 is a model tuned to the Fu 2017 WFC3 amplitudes (N=34); the 950 K is the CH4 line plus the silicate sink depth, so compute it from the chemistry module; Estrela 2022 (N=62) gives a different quantity (cloud minimum at 1461 +316/-405 K) | crossover Teq; no published width | A_H against Teq: slope 0.0008 +- 0.00016 per K, binned Spearman 0.85 | none |
| n7 haze yield | S | lab rates span 250 times but were measured on sub-Neptune gases at 300 to 600 K, in mg/hr, not as a fraction; every population statistic is degenerate with metallicity, Kzz and cloud cover | none; either a nuisance parameter with the Brande scatter (about 0.5 scale height, N=15) as prior, or a haze production rate per unit CH4 photolysis for hot gas | none | the only calibration is a Titan simulation |

### 3.5 Comets, belts, moons, surfaces (rows k; body T13 to T19, T21, E2, E3, R1, R3)

| Row | Class | Finding | Knobs -> statistic -> target | Test | Solar System flag |
|---|---|---|---|---|---|
| k1 to k3 coma class | F | counts 13/9/3 of 25 (Clopper-Pearson 68 percent: 0.40 to 0.64, 0.25 to 0.48, 0.055 to 0.22). Class depends on r_h (every coma beyond about 3.5 AU is CO-dominated) and 9 of 13 CO2-dominated entries rest on a CO upper limit. Better: draw log10(CO/CO2) at the comet's own r_h and derive the class | (CO+CO2)/H2O median 18 +- 4 percent (N=18) then a CO fraction | class counts within the binomial intervals under the HP22 admission rule | within-system scatter is Solar System only |
| k4, k5 | T (Solar System), F with an offset | CO2/H2O geo-mean 10.9 percent (0.11 dex SE), sigma 0.51 dex (N=21); CO/H2O 2.5 percent (0.16 dex), 0.65 dex (N=16); HP22 medians agree (12 +- 2, 3 +- 1). A free covariance is not identifiable; a one-factor model fits (mean pair correlation 0.13). 2I/Borisov and 3I/ATLAS lie 1.8 to 3.6 sigma outside the Solar System fits (3I CO2/H2O varies 36-fold across epochs) | mu, sigma as published; per-system offset delta_sys ~ N(0, sigma_sys), upper edge near 1.5 dex, from N=2 | geo-means within 2 SE; 2 of 2 interstellar draws above CO/H2O 0.3 must not be improbable | yes: the Solar System population is not a prior for other systems |
| k7 single-comet species | S (partly T) | populations exist for CH4 (0.68 percent, N=12), C2H6 (0.41, N=18), H2CO (0.25, N=15), H2S (0.71, N=9), OCS (6 comets); N2, O2, Ar, Kr, Xe, Ne, SO2 are 67P only | as published | as above | O2 presence cannot be set from one comet |
| T15 organics p=0.5 | P | Marboeuf runs two extreme models; there is no probability. Use a continuous organic fraction N(0.60, 0.15) clipped to [0, 1] | f_org | none | partly |
| T13 taxon mix, R1 | F | one zone mix of one belt (34,503 asteroids; DeMeo and Carry Table 6, no uncertainties); no extrasolar taxon ensemble; white-dwarf pollution says exoasteroid rock types should not be assumed Solar System-like. The 0.2 dex snow-line scatter is a model consequence of an assumed 0.4 dex accretion-rate spread | sigma_snow 0.2 dex | self-consistency: the pooled mix at the Solar System x equals Table 6 | yes, one belt |
| T14 within-taxon weights | F | fall shares (H 34, L 38, LL 8 percent) and NEO spectral shares (about half LL) disagree | fall shares as a flagged placeholder | none | yes |
| E2 belt density | P | derive from belt mass per area; debris-disc distribution exists for cold belts only (Sibthorpe 2018: 17.1 percent of 275 stars) | thresholds on belt mass | none | none |
| E3 field families | S | no ensemble; a design vocabulary. Derive the family from the parent belt's T13 shares or keep as a flagged stand-in | none | none | yes |
| T16 comet volatiles | T | same ensemble as k4, k5 plus Willacy 2025 (12 comets) | as published | as above | as above |
| T17 moon origin | S | rules exist, no population gives weights; zero confirmed exomoons. Disk moons: M_T/M_p about 1e-4 (Canup and Ward); impact moons: 10 of 28 SPH collisions made discs above 0.1 Earth masses (0.36, CP68 0.24 to 0.50) | derive P(impact moon) = impact-energy frequency (Quintana: median 3 Moon-forming impacts, N=164) x disc yield | scaling within a factor 2 | yes |
| T18 moon ice fraction ramp | S | Sasaki (r_ice = 20 R_J (M_p/MJ)^(2/3) (tau/5 Myr)^(-1/3)) and Cilibrasi (85 percent icy, Galilean-like in 20 percent) are both READ and disagree; neither gives a ramp width | bracket 0.5 to 0.85 icy | none | yes |
| T19 dwarf planet and KBO rock fraction | F | 17 binary KBO densities (0.5 to 2.6 g/cm3); both READ sources read them differently (rock fraction rising with size against constant 70 percent with porosity doing the work); dwarf planets 0.74 to 0.82, small KBOs 0.3 to 0.45 (COMPUTED); Bierson and Nimmo not reachable | dwarf-planet rock fraction 0.78 +- 0.12 (N=7); porosity transition near 350 km | density-diameter relation within the 17-object scatter, with a binary-detection selection | yes, no extrasolar KBO density |
| T21 grain state | S (derive) | one extrasolar object: LHS 3844 b (Bond albedo 0.14 +0.13/-0.14; fresh powders excluded; slabs or old weathered powder fit). Texture changes albedo by up to a factor 7 (Paragas 2025, 11 samples) | P(powder) rises with surface age x impact flux; P(slab) with resurfacing; the Paragas library maps texture to albedo | LHS 3844 b brightness-temperature ratio 0.96 +0.03/-0.04 (a calibration point, not a prior) | yes |
| R3 meteorite groups as outputs, R5 fence | S | no ensemble beyond the Solar System; R5 needs no data | n/a | n/a | yes |

### 3.6 Rows outside the 48 that the workers also checked

d1 to d5 and p4 (above), e1/e2 (above), f5/f6 (above), f8 (carbon carrier, derived). T3 (white-dwarf type against temperature) and T7 (core light elements) were **not checked in this pass**. Lane judgement, not verified: T3 is very likely class T (white-dwarf spectral evolution is a classic large ensemble; RECALLED candidates Kleinman 2013, Kepler 2019, Genest-Beaulieu and Bergeron 2019); T7 is class X or P (planet cores are unobservable; derive light-element content from redox and pressure partitioning, Hakim 2018 and Unterborn 2023 ranges). Section 6 item 11 asks for a short follow-up on both.

The body-composition "water share switch" (Marboeuf-derived 0.14 to 0.39 against about 0.5, GEN.468): the planet data cannot discriminate (NGPPS water planets 0.03 to 0.49, median 0.28; Chakrabarty 2026 peak near 41 percent; ice-line material about 50 percent in Bean 2021). The discriminating source stays Lodders 2003 (Boss's rank 1 research item), so the switch stays.

## 4. Builder groups (one worker each) and the shared harness

| Builder | Rows | Shared data and fit | Needs |
|---|---|---|---|
| B1 Radius valley and envelope | e1, e2, e-link | CKS DR1/DR2/Gaia tables, Kepler DR25 injections, Petigura detection model, Rogers and Owen and Ho 2024 targets | injection-recovery selection module (shared with B2); decision D2; Lopez and Fortney 2014 radius law |
| B2 Delivery and water worlds | b1 to b6, e3, e4, T9, T11, switch | NGPPS NG76 extraction, Luque and Palle, Parc, Dainese, TRAPPIST-1 (Agol) | decision D1; DACE formation-stage tables |
| B3 Retention and lava | d1 to d6, j1, j5, i6 | Berta-Thompson label table, Meni-Gallardo database, JWST rocky-planet list | JWST target-selection module; Berta-Thompson code and labels |
| B4 Iron branch | T5 | Plotnykov 2025, Brinkman 2024, Adibekyan 2021, Schulze 2021 | mass and radius error injection |
| B5 Redox and WD debris | a1 to a4, T8, l1 to l3 | Doyle 2020 and 2023 tables (transcribe 24 rows), Koester 2014, Buchan 2022, Harrison 2021 | WD detection and settling model; decision D3 |
| B6 Radioisotopes | o2, o3, o4, o6, T23 | Hypatia API (Th, Eu, K), Botelho, Mors-style cluster results (Lichtenberg 2016), Th/U formula | Hypatia licence; a chemical-evolution r-process history |
| B7 Ground and petrology | T20, T22, T21, R2 | petrology on the stored bulk; Putirka and Xu, Trierweiler, Paragas | Putirka and Xu Table A1 |
| B8 Giants | f5 to f7, T12 | Thorngren 2016, Chachan 2025, LExACoM CSV, Brewer and Fischer, Suarez-Andres | machine-readable T&F 2019 and Chachan tables |
| B9 Brown dwarfs | g1, g4 | Kirkpatrick 2021 Tables 13 to 15, Beiler 2024, Filippazzo 2015 | Tables 11, 14, 15 machine-readable |
| B10 Chemistry module | h2, h3, n6 crossover, f8, m7 to m9, n1 to n4, m2 | easyCHEM or own Gibbs solver with NASA-format polynomials | licence check on the polynomial source |
| B11 Aerosol amplitude | h4, n6, n7 | Fu 2017 (N=34), Estrela 2022 (N=62), Brande 2024 (N=15), UHJ TiO/VO list | decision D5 |
| B12 Comet sampler | k1 to k5, k7, T15, T16 | Robinson database (CC BY 4.0), HP22 Table 4, Willacy 2025, interstellar objects | decision D6; the three review papers |
| B13 Belts | T13, T14, E2, E3, R1, R3 | DeMeo and Carry Table 6, fall shares, Mulders | stand-in flags only |
| B14 Moons and KBOs | T17, T18, T19 | Canup and Ward, Sasaki, Cilibrasi, Quintana, Brown, Barr and Schwamb | Bierson and Nimmo |
| B15 Host-star activity | p4, p6 (and p3, p5) | Jackson 2012, Medina 2022, Wright 2011 and 2018, Johnstone 2021 tracks | the Mors percentile data |
| B16 Halogens and trace gases | m5, m6, m2 look-up | Schaefer 2012, Asplund 2021, Maas, Jonsson | Cl and F depletion source |
| B17 Design-flag registry | j2 to j4, m3, m4, m11, k7 rest, E3, T17, T18 | none | decision D4 |

**The shared tuning harness** (one build, used by B1 to B16): a target table with columns row id, statistic, value, sigma, source, label, validity range and selection-model id; a simulate-and-compare function per target (chi-square within 2 sigma, or binomial/Wilson for counts); and five selection models (Kepler transit plus recovery; mass-radius sample cut of mass better than 25 percent and radius better than 8 percent; JWST atmosphere target selection; comet admission by heliocentric distance; white-dwarf pollution detection with settling). The loader refuses a tunable row with no target-table entry (extends TEST.214). Knob scope is the validity range of the data (for example [Fe/H] -0.5 to +0.5 for the CKS metallicity null; 20 pc for brown dwarfs; FGK for the radius valley); outside it the draw is a flagged extrapolation.

**Built (GEN.503, lane M, 2026-10-23):** `planetgen/generation/tuning_harness.py` and the target table `planetgen/physics/data/tuning_targets.json`. A target has the row it tunes (`<table>:<row>`), the statistic (`mean`, `sd`, `median` or `proportion`), the published value and sigma (a proportion states its count and takes the Wilson interval), the source and its label, the validity range, the selection-model id, the fill id and the falsifier. `simulate_and_compare` draws one value per seed on the stream `tuning:<target id>:<i>` and passes within 2 sigma, where sigma joins the published uncertainty and the sampling error of the simulation. A possibility-table `Row` carries `tunable`; `TargetTable.require` (run when the atmosphere recipes load, and by `Registry.validate(targets=...)`) refuses a tunable row without a target that may be tuned to, and a target for a row that is not tunable. A value measured only in the Solar System is a `test-point` with the `solar_system` flag and is never tuned to; the loader refuses a `target` that names a Solar System body. A SNIPPET, RECALLED or NOT REACHED value is listed and never tuned to. A target that names a selection model refuses to simulate until GEN.504 to GEN.508 register it. Worked example: rows a3 and a4 (the 22 white-dwarf redox values) against mean -1.25 (sigma 0.10) and sd 0.45 (sigma 0.069). A retune is a table version change and a generator epoch bump; adding a target changes no draw. The other builders (GEN.509 to GEN.525) mark their rows `tunable` in `scripts/atmosphere_recipes_data.py` and add their targets to the table.

## 5. NO FIXED CONTENT: where the Solar System alone constrains a knob

Tuning a knob to these would break the rule; keep them as test points or flagged stand-ins.

| Where | What is Solar System only |
|---|---|
| d table | repo K window (13 bodies); d4 and d5 pressure caps (Mars, Titan, Venus, Earth) |
| a1, a2, T8 | the Earth end (+0.5 on Sossi's scale), the 1.25 slope, the later-mantle law fitted to Moon, Mars, Earth |
| b1 to b3, T11, b6 | Ortenzi payloads (Earth scenarios), Earth's water budget, rogue defaults |
| o3, o4, T23 | K/U and K depletion (Moon, Earth, Mercury, Mars), Th/U formation duration (CI chondrites) |
| j2, j3, j4, m3, m4, m11 | Mars, Titan, Pluto, Triton, Earth N2O, Venus PH3, modern-Earth lightning |
| m5, m6 | crust enrichment 26.5 and 11.8, Venus HCl and HF test values |
| k4, k5, k7, T16 | within-system scatter (N = 9 to 32 comets, all Solar System); N2, O2 and noble gases from 67P alone |
| T13, T14, E3 | one belt's taxon mix and meteorite fall shares |
| T17, T18, T19, T21 | Solar System moons, KBO densities (17), lunar regolith |
| f5, T12, n7 | f = 0.2 (Jupiter, Saturn), the interior structure of the two giants, the Titan lab haze rate |

## 6. What Boss must decide or supply, ranked by leverage

**Decisions (defaults picked so work proceeds; none parks the thread):**
- **D1 water-rich share (b4, e3, e4).** Default: keep the formation model as a switch with the NGPPS logistic as the MODEL default and the flagged bracket 0 to 0.87; the alternative is to wait for a formation model.
- **D2 envelope loss mechanism (e1/e2).** Default: XUV energy-limited with variable efficiency (Owen and Wu, the repo's choice), a core-mass distribution as a third knob, and Lopez and Fortney 2014 as the retained-envelope radius law (not yet fetched).
- **D3 which IW scale each redox row draws.** Already decided as one scale (GEN.498, Hirschmann 2021); the decision left is to state that WD dIW is the core-formation scale, a1 and a2 the surface magma-ocean scale, and to add the conversion test (G-4).
- **D4 life-linked gases (m3, m4, m11).** Default: a flagged design bracket, default off unless the life-stage flag allows, never tuned.
- **D5 Kzz and haze yield (h4, n7).** Default: nuisance parameters with the Brande scatter as prior; drop the LU[1e7, 1e11] Kzz draw in favour of a Teq threshold.
- **D6 comets.** Default: accept the Solar System population as the within-system shape plus a per-system offset knob bracketed by the two interstellar objects.

**Papers and data to obtain (named exactly):**
1. Putirka and Xu 2021 Supplementary Table A1 (23 white dwarfs; "on request" from the lead author), Hollands 2018 compositions (MNRAS 477:93) and Xu 2019 (AJ 158:242) tables, and the Doyle 2019, 2020, 2023 tables as machine-readable files.
2. NGPPS formation-stage tables (DACE, NG76 or the NGPPS II populations): planet-by-planet formation start distance, core ice and rock mass, final a and mass; DACE terms of use (not stated on the pages reached).
3. Thorngren and Fortney 2019 (403 giants) and Chachan 2025 (147 giants) as machine-readable tables; Brewer and Fischer 2016 Table 1; Kirkpatrick 2021 Tables 11, 14 and 15.
4. The Berta-Thompson 2025 code and label table and the ExoAtmospheres (IAC) export.
5. CKS DR2 tables (VizieR J/AJ/163/179), Fulton and Petigura Table 4, the DR25 injection efficiency tables; Lopez and Fortney 2014; the source behind E_cum = 3.1e16 J/m2 (the XUV history).
6. Dello Russo 2016 (Icarus 278:301), Mumma and Charnley 2011, Bockelee-Morvan and Biver 2017 (comet ranges); Bierson and Nimmo 2019; Schaller and Brown 2007; Harman 2018 full text.
7. Lodders 2003 (already Boss's rank 1), Hirschmann 2021 full text, Deng 2020 Fig. 3 data, Frank, Meyer and Mojzsis 2014 (r-process history for Th/U), a primary Na and K condensation-temperature source.
8. Licences: Hypatia (copyright notice, no licence stated), the Zenodo NGPPS record (open, no licence field), Zenodo 17986129 (Paragas texture data), LExACoM CSV (no licence stated), CKS and DR25 tables.
9. Schreiber 2019 and Wilson 2019 (the l3 upper bracket and the volume-limited polluted fraction 45 +- 4 percent), JWST rocky-planet papers that postdate the workers' reach.
10. A follow-up on T3 and T7 (section 3.6).

**Absent from the public record, not merely locked** (nothing to supply; build the derivation and flag the row): the mass dependence of mantle redox beyond theory; water as a function of formation distance outside population synthesis; Kzz and haze yield as populations; moon origin weights; surface (not bulk) rock-type frequencies of any exoplanet; K depletion beyond the Solar System; envelope Z occurrence; planet-side halogen depletion; coreless-planet counts; and any measurement of N2O, PH3, NO2 or CH3Cl on an exoplanet.

## 7. Findings for the TODO thread

Each is a failure or mismatch this lane does not own. Check each against current main before filing (the repo doc claims are from `docs/design/atmospheres-retention-and-classes.md` as the workers read it).

| # | Kind | Object | Finding | Evidence |
|---|---|---|---|---|
| 1 | design | d table K window `LU[3.0e-4, 6.6e-4]`, slope 4 | contradicted by exoplanet fits and by 55 Cnc e and TOI-561 b (need K 12 to 190 times higher); replace by the Berta-Thompson logistic | arXiv:2507.02136, 2508.12865 |
| 2 | doc | retention doc | "Meni-Gallardo and Palle keep slope 4" is wrong (slope 5.77); "TRAPPIST-1 b and c under 0.1 bar" is in no paper read; Zahnle and Catling give no K | arXiv:2508.12865, 2303.14849, 2306.10150, 1702.03386 |
| 3 | design | a4 N(-1, 0.5), a1 floor -2, a3 U[-1.4, 0] | WD data centre -1.25 +- 0.10, sd 0.45 | arXiv:2008.05992, 2303.00063 |
| 4 | design | f7 `f_tail` LU[0.0013, 0.34], core 0.55 | stellar tail is zero above 0.8; planet tail 0.043 to 0.087; core centre 0.57 (LExACoM) | arXiv:1608.06286, 2510.26785 |
| 5 | design | f5 note "Chachan gives lower Z for sub-Saturns" | opposite; Chachan is 1.4 to 1.6 times lower than Thorngren from 0.3 to 3 MJ | arXiv:2509.20428 |
| 6 | design | g4 row | the 500 to 700 K "gap" is T8 to T9; g1 to g5 are one Teff map | Kirkpatrick 2021 Table 13 |
| 7 | design | h2/h3 edges 1500 and 2200 K | pressure-dependent outputs; fixed edges not supported | worker Gibbs table |
| 8 | design | h4 TiO and VO edges | Ti retained only above about 2350 to 2640 K; VO kept at 2210 K | arXiv:2305.17228 and others |
| 9 | design | k1 to k3 weights | replace by printed counts and an r_h term; add per-system offset knob | arXiv:2209.09985, 2004.09586, 2508.18209 |
| 10 | design | p4 `LU[0.162, 2]` Gyr, p6 weights | tau_sat is mass dependent; p6 becomes a percentile draw | arXiv:1111.0031, 2205.02331, 2009.07695 |
| 11 | doc | m2 "O3 about 10 percent of Earth" | 0.01 to 0.1 percent for G, K, M hosts | arXiv:1407.2622 |
| 12 | doc | part 09 R2 source | Thiabaud 2014 is about refractories, not water | arXiv:1312.3085 |
| 13 | design | `physics` ozone and halogen rows | HCl and HF percent-level only for crust composition; bulk silicate Earth below 0.1 percent | arXiv:1108.4660 |
| 14 | design | o6 Al-26 range | extrasolar evidence points to a mixture (enriched flag), not log-uniform around the Solar System value | arXiv:2206.09675, 1608.01435 |
| 15 | design | o4 Th/U U[3.8, 4.5] | an age-4.5 Gyr Solar System value; derive from production ratio, birth epoch and decay | worker formula |
| 16 | doc | candidate-list errors in earlier notes | Thorngren id, Kane 2014 = 1409.2886, Blackman 2021 = 2110.07934, "Burn 2021 NGPPS III" = NGPPS IV, Petigura 2022 = 2201.10020, Gupta and Schlichting 2019 = 1811.03202, Zieba 2019 has no composition, DeMeo 2019 vs 2022 | worker files, section 3 each |
| 17 | process | selection models and target table | build once (section 4) | this report |
| 18 | licence | Hypatia, NGPPS Zenodo, LExACoM CSV, Paragas Zenodo | no licence stated; confirm before distributing derived tables | worker files |

## 8. Worker files (detail, N, access, licences, recalled lists)

All under `docs/design/tunable-gaps-sources/` (copied unchanged from `/mnt/project-files/reports/tunable-gaps-2026-10-17-sources/`):
- [brief.md](tunable-gaps-sources/brief.md) (the common brief)
- [group1-subneptunes.md](tunable-gaps-sources/group1-subneptunes.md) (e rows, b4)
- [group2-delivery-retention.md](tunable-gaps-sources/group2-delivery-retention.md) (b, d, i rows, T5, T9, T11)
- [group3-redox-wd-radio.md](tunable-gaps-sources/group3-redox-wd-radio.md) (a, o rows, T8, T20, T22, T23, R2, R4)
- [group4-giants-bd-clouds.md](tunable-gaps-sources/group4-giants-bd-clouds.md) (f, g, h, n rows, T12)
- [group5-small-bodies-moons.md](tunable-gaps-sources/group5-small-bodies-moons.md) (k rows, T13 to T19, T21, E2, E3)
- [group6-evolved-trace-uv.md](tunable-gaps-sources/group6-evolved-trace-uv.md) (j, l, m, p rows)

## 9. Limits of this report

- Worker numbers keep their own labels; this lane re-ran three statistics and read every report but did not reopen the papers behind every figure.
- Several figures are SNIPPET or RECALLED in the worker files (Kepler-138 d volatile share, the GJ 3929 b follow-up, Wilson 2019, Hot Rocks eclipse depths, Venus HCl and HF); none is used as an input here.
- Statistics from retrievals (LExACoM) are method dependent (0.2 to 0.3 geometry systematic in C/O); they are tuning targets with a wide tolerance, not truth.
- Population-level tuning uses what astronomers could observe (FGK dwarfs in the Kepler field, nearby M dwarfs, JWST targets chosen for signal); the generated galaxy covers other populations, so each knob carries a validity range and extrapolation is flagged.
- Classes T, F, P, S, X are this lane's own grouping of the workers' (a), (b), (c) verdicts and "derive" recommendations, written so a builder can see at a glance whether the work is a fit, a formula or a decision.
