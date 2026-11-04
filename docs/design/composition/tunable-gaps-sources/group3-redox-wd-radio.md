# Group 3: mantle redox, ground types, radioisotopes (white-dwarf debris as the extrasolar ensemble)

Read-only source check for the tunable-gaps question. Labels: READ (own reading of the paper text or abstract on arXiv/Nature), 2ND-HAND (read in a reachable paper or project part that cites the original), DATA (fetched data, named), COMPUTED (own Python from stated inputs), SNIPPET (search summary only), RECALLED (memory), NOT REACHED. Fetch date for all DATA: 2026-10-17. Repo and project text were not edited.

## 1. Summary and verdicts

The best new finding is a real extrasolar redox ensemble: polluted white dwarfs (WD) with O, Mg, Si, Fe measured give an intrinsic oxygen fugacity (dIW at metal-silicate equilibrium) for 24 parent bodies (22 values, 2 upper limits): median -1.26, mean -1.25 (standard error 0.10), sample sd 0.45, crude intrinsic sd 0.22 for the 13 best-constrained stars (all COMPUTED from READ Tables of [R1152, R1153]). That anchors the small-body end of the redox draw without the Solar System, and it contradicts two project values: the a4 centre -1.0 (from Doyle's first six stars) and the a1 floor -2. There is no published observational constraint on the mass dependence of mantle redox beyond theory ([R1075] with three Solar System bodies, [R1717] with no calibration). Ground types (T20) and crust/mantle end-members (T22) cannot be anchored to WD debris except for bulk-mantle mineralogy and a crust-bearing frequency of about 1 percent; surface ground types are not observable in disrupted bodies. For stellar Th, U, K, Eu there is a public catalogue (Hypatia, open API, no stated licence) with Th N=284, Eu N=4570, K N=1509, U none; Th/U in a new planet is derivable from a nucleosynthesis production ratio, the r-process duration before birth and radioactive decay, but the duration is calibrated on the Solar System today. For Al-26 there are extrasolar constraints (WD differentiation, gamma-ray star-forming-region levels, cluster simulations) that support a bimodal presence-plus-level draw, not a log-uniform one.

| Row | Quantity | Verdict | One-line reason | Status of anchor |
|---|---|---|---|---|
| a1 | Magma-ocean fO2 vs mass | (c) for the mass slope; (b) for the low-mass end | low-mass end: WD median -1.26, N=22 (extrasolar); high-mass end: Earth only (Solar System); slope: theory only. Derive by physics (base pressure to surface fO2), draw only the scatter | READ (WD), READ (Deng, 3 Solar System bodies) |
| a2 | Later upper mantle | (c) | only Solar System bodies plus theory (Guimond: at least 2 dex mineralogical spread at fixed Fe3+/sumFe); repo law is "fitted by eye" to Moon, Mars, Earth | flag: Solar System only |
| a3 | Building-block redox U[-1.4, 0] | (a) after replacing the uniform by the WD distribution | N=22 values; U[-1.4,0] mean -0.7 sits 0.55 dex above the WD mean | READ |
| a4 | WD debris N(-1, 0.5) | (a) | same ensemble; centre should be -1.25 (SE 0.10), width 0.45 total or about 0.22 to 0.35 intrinsic | READ, COMPUTED |
| o2 | Th scale for non-thin-disc hosts | (b) | Hypatia Th: thin N=217, thick N=33; Th/Mg shows no thick-thin offset; halo, young, M hosts have no Th data | DATA |
| o3 | K/U | (b) via Na as proxy, otherwise (c) | stellar K N=1509 but planetary K is volatile-depleted; extrasolar depletion only via WD Na (N=202 DZ) or one WD with K; Solar System K/U (4 bodies) must not be the tuning target | READ, flag |
| o4 | Th/U | derive by physics; (a) for the formula | formula below; production ratio spread READ; U never measured in a disc star | READ, COMPUTED, flag |
| o6 | Al-26 present c[0.1,1] | (b) | WD-based lower bound on enriched fraction at least 4 percent, tens of percent likely; Lichtenberg 9 to 28 percent; not a pure probability if tied to birth cluster | READ |
| T8 | Mantle oxidation state | same as a1 to a4; largest hole is mass and age dependence | no per-planet distribution exists | READ |
| T20 | Ground type weights | (c) for surface types; (b) for 3 sub-statistics | WD gives bulk mantle, core-rich 8/202, crust-bearing about 1 percent; no surface rocks | READ |
| T22 / R2 | Crust and mantle end-members | (c) for a drawn table; derive by petrology | Putirka and Xu data table is "on request" only; WD gives test targets (normative mineralogy of 23 PWDs) | READ |
| T23 / R4 | K depletion, Al-26 range | K: (b)/(c); Al-26: (b) | see o3, o6 | READ |

Usability-rule overrides named (verified beats partly verified or recalled): (1) [R1152]/2023 WD dIW (READ) overrides the RECALLED a4 centre -1 and width 0.5; (2) the same data overrides the a1 floor -2 for small bodies, which came from Solar System Moon theory; (3) [R0479] thin-disc Th scatter 0.056 dex (READ, N=53) overrides Hypatia's 0.14 to 0.18 dex as an intrinsic width (Hypatia mixes catalogues); (4) Putirka and Xu 2021 and [R4425] both READ and disagree on "exotic" versus "chondritic" mantles, so no override: both are kept as a spread.

Answers to the five specific questions: (1) 24 WD systems, centre -1.25, sd 0.45, bodies are asteroid to Vesta/Moon-mass for 13 systems and unknown for the rest, so dIW says nothing about mass dependence (section 2.1). (2) No published constraint on the mass dependence beyond theory (2.2). (3) Only partly (2.4). (4) Hypatia, with numbers in 2.5; Th/U is derivable, formula in 2.6. (5) Yes for presence and level, no for log-uniform shape (2.7).

## 2. Detail

### 2.1 Redox ensemble: polluted white dwarfs (rows a3, a4, part of a1, T8)

**Ensemble.**

| Item | Value | Source | Status |
|---|---|---|---|
| [R1154]:356 | 6 WDs (O, Fe, Mg, Si, Ca, Al measured, IR-excess disks required): dIW -1.17, -1.37, -1.14, -0.59, < -1.21, -0.63 | [R1154] via project part G3 table; text of arXiv version read for method | 2ND-HAND (values), READ (method) |
| [R1152], ApJ, "Where are the extrasolar Mercuries?" | 16 WDs, Table 1 (median dIW, 1 sigma, Teff, Oxs/Fe). Values: J1043+0855 -2.13 (no lower bound), 1536+520 -1.37, GD 40 -1.26, J0738+1835 -1.14, 1226+110 -0.63, 1145+017 -0.70, 1929+012 -0.49, Ton 345 -2.10 (no lower bound), G 29-38 < +0.95 (upper limit only), HS 2253+8023 -0.94, 1350-162 -1.01, 1232+563 -1.21, 2207+121 -1.33, 0446-255 -1.47, G 241-6 -1.54, 1551+175 -1.40. Per-star errors 0.07 to 0.43 where bounded | [R1152] | READ |
| [R1153], "New chondritic bodies in eight oxygen-bearing WDs" | 8 more WDs (He-atmosphere): GaiaJ0218+3625 -1.29, WD1244+498 -0.54, J1248+1005 -1.09, WD1415+234 < -0.87 (limit), J1734+6052 -1.04, GaiaJ1922+4709 -1.78 (no lower bound), EC22211-2525 -1.25, J2248+2632 -1.74 (no lower bound) | [R1153: Table 4] | READ |
| Combined | N=24 systems; 22 values, 2 upper limits; 7 of 24 have no lower bound | COMPUTED from the two tables | COMPUTED |
| Statistics (22 values) | mean -1.248, median -1.256, sd 0.447, range -2.13 to -0.49 | COMPUTED | COMPUTED |
| Statistics (13 stars with both bounds) | mean -1.115, sd 0.338, rms per-star error 0.259, intrinsic sd about 0.22 by quadrature subtraction (crude) | COMPUTED | COMPUTED |
| Fraction reduced | [R1152]: 11 of 16 have dIW > -3 (oxidised relative to solar gas), 5 of 16 cannot exclude reduced bodies (about 1/4 to 1/3). Mercury-like (dIW < -3) not confirmed in any | [R1152] sec. 4.2 | READ |
| Selection | stars need O, Mg, Si, Fe (and Ca) with errors; in 2020, 63 percent have confirmed IR excess; settling in DBs at Teff 16000 to 19000 K can fake dIW < -3 in a declining phase; method validated for dIW about -8 to about -0.9 on Solar System objects | [R1152] sec. 2.1, 4.1; [R1154] SI | READ |
| Measurement scale | Eq. 3: dIW = 2 log(x_FeO,rock / x_Fe,metal), activity coefficients set to 1, x_Fe,metal = 0.85 (varying 0.8 to 1.0 shifts dIW by about 0.2); plausible activity coefficients shift all values by at most about 1 log unit. [R1154] rounds dIW = dQFM - 4 | [R1152] sec. 1, 3.1 | READ |
| Public table? | Tables are in the papers (ApJ, open arXiv). No separate machine-readable file found; transcribing 24 rows is trivial. Licence: arXiv default, not stated in the PDFs | | READ |
| Scope | Pollution fraction: 25 to 50 percent of WDs polluted ([R2415, R4918]: 48 of 85 DAs with metals, 23 (27 percent) must be accreting, up to about 50 percent) | [R2415] abstract | READ |
| Size of the bodies | [R1832]: 13 of 202 DZ systems have finished accretion; accreted masses Moon to half of Vesta. [R0573]: of 42 systems, 3 favour small parent bodies (differentiation pressure 4 to 8 GPa), 2 favour Earth-sized (GD 61, WD 0446-255), 4 degenerate with fO2, 5 unconstrained, 26 no evidence of differentiation, 1 unphysical. [R1154]: bodies are "asteroids ... or fragments of rocky exoplanets" | [R1832] abstract; [R0573] sec. 3 | READ |

**Physical bracket.** dIW of a body with metal is set by FeO in rock against Fe in metal (Doyle Eq. 3). Solar gas is dIW -6 to -7 (via [R1152], 2ND-HAND); Solar System rocks span about 6 dex [R1154]. Rocks with no metal give only upper limits ([R1152] sec. 2.3).

**Verdict.** (a) for the small-body, core-formation quantity (rows a3, a4). Be blunt about two limits. First, Doyle's dIW is the fO2 at metal-silicate equilibrium on a scale with Earth at about -1.5 (Earth -1 to -2 depending on the FeO activity-coefficient ratio, G3 table, 2ND-HAND), whereas a1 uses the surface magma-ocean scale with Earth at +0.5 (Sossi) or +2 (Deng). These are different quantities and the project must say which one each row draws (G3 finding F2: the generator should fix one IW scale, [R1951], metadata verified only). Second, WD bodies have no known mass for 11 or more of the 24 systems, so the ensemble cannot supply any mass slope.

**Knobs.**

| Knob | Target statistic | Target value (uncertainty) | Range of validity |
|---|---|---|---|
| centre mu_small | median and mean of 22 WD dIW | -1.26 / -1.25 (SE 0.10; scale systematic up to 1) | small bodies (Moon to Vesta mass, 13 systems; others unknown), host progenitors 1 to 3 Msun |
| width sigma | sd of 22 values; sd of 13 well-constrained | 0.45 (SE about 0.07) ; intrinsic about 0.22 | same |
| tails | lower tail to -2.1 (3 of 22 below -1.7, all with weak lower bounds); upper edge -0.49 (method cannot resolve above about -0.9) | one-sided truncation at -0.5 not physical, it is a detection limit | |
| reduced-body fraction (dIW < -3) | 0 of 24 confirmed | upper limit about 0.12 at 95 percent (rule of three, COMPUTED) | |

**Consistency test.** Compare simulated median, sd and a KS test against the 22 values with tolerance the standard errors above; the simulated side needs a detection model (per-star error 0.07 to 0.43, truncation near -0.9 for oxidised bodies, no lower bound for 7 of 24, DB settling bias for Teff 16000 to 19000 K as in [R1152: Figs. 8 to 11]). Without a selection model the width is an upper bound.

**No-fixed-content.** None for the low-mass anchor (extrasolar). The a1 high-mass anchor (+0.5 at 1 Earth mass) is Earth only.

**Repo comparison (usability rule).** a4 N(-1, 0.5) used the first six stars (their mean of five values is -0.98, COMPUTED). The 22-value mean -1.25 differs by 2.6 standard errors; keep sd 0.45 or use the intrinsic 0.22 plus measurement noise. a3 U[-1.4, 0] has mean -0.7; 6 of 22 WD values lie below -1.4 but only 2 (-1.47, -1.54) are well constrained and both sit within 1 sigma of -1.4, so the uniform's lower edge is not hard-contradicted, its centre is. [R3944] (ApJ 843:120, [R3944]) meteorite values stay as Solar System test points.

### 2.2 Mass dependence of mantle redox (rows a1, a2, T8)

| Source | What it gives | Calibration | Status |
|---|---|---|---|
| Deng, Du, [R2614], Nat. Commun. 11:2007, [R1075] (the candidate "[R1075]" is right; authors Deng, Du, Karki, Ghosh, Lee) | First-principles plus thermodynamic model: magma-ocean surface fO2 rises with ocean depth. Base depths 55, 14, 5 GPa give base dIW -2, -1.5, -2 and surface about IW+2 (Earth), IW-0.3 (Mars), IW-2 (Moon). 1 sigma about 0.5 log. Earth surface falls about 1.5 log if base is 25 not 55 GPa. Mentions super-Earths only qualitatively | three Solar System bodies; no exoplanet datum | READ (Nature page) |
| [R4154], Sci. Adv. 6:eabd1387, [R4154] | Earth magma ocean about IW+0.5 on their scale at 2173 K, Fe3+/sumFe 0.037 +- 0.005 | Earth only | 2ND-HAND (part G3), metadata verified |
| [R1717], MNRAS, [R1717] | At fixed Fe3+/sumFe, mineralogy alone gives at least 2 dex fO2 spread across exoplanet mantles from Hypatia host-star refractory abundances; about 3 dex from Mg/Si; "Fe3+/sumFe and thus the real fO2 is unconstrained for exoplanets" (fixed at 3 percent for figures) | no calibration | READ |
| Gaillard, Bernadou, Roskosz, EPSL 577:117255 | Outgassed composition vs fO2 regimes (above IW+2 dry C-N-S atmosphere) | theory plus experiments | 2ND-HAND (part G3), metadata verified |
| [R1154]/2020/2023, [R0573, R1832] | Overlap check: WD 1350-162 (small parent, Buchan) dIW -1.01 +0.18/-0.21; WD 0446-255 (Earth-sized, Buchan) dIW -1.47 +0.20/-0.22. Difference 0.46 dex at about 1.5 sigma, N=2, differentiation pressure degenerate with fO2 in 4 systems | | COMPUTED from READ |

**Answer to question 2.** No published constraint on the mass (or pressure) dependence of mantle redox beyond theory exists in what I could reach. The only observational tie is the two-star overlap above, which has no power. a1's slope (1.25 dex per dex of mass) was interpolated between a WD-like floor and Earth, so it is a design interpolation.

**Verdict.** a1 slope and a2: (c). Recommend deriving the surface dIW from a physical model (base pressure of the ocean from M, Fe3+ budget, mineralogy) and drawing only a Gaussian scatter (Deng 0.5, Guimond 0.3 to 0.4 dex 2ND-HAND via G3). Boss would have to supply Deng's Fig. 3 data (G3 item 6) and a source for Fe3+/sumFe in exoplanet mantles (none exists). The a2 repo law (3.5 + 2 log10 M, sigma 1.0) is fitted to Moon, Mars, Earth only: Solar System only, flag.

### 2.3 Row a1 specifics under the usability rule

a1 as written: mu = clamp(-2 + 1.25 (log10 M + 2), -2, 0.5), sigma 0.5. WD data say small bodies centre at -1.25, not -2. Replace the floor by -1.25 (SE 0.10) and keep the Earth end at +0.5 (Sossi scale, Solar System only) until a second exoplanet anchor exists. Scale mismatch: Doyle's scale and Sossi's scale differ; do not mix without a conversion test (G3 G-4).

### 2.4 Ground types and crust/mantle end-members (T20, T22, R2)

| Ensemble | N | Finding | Source | Status |
|---|---|---|---|---|
| Putirka and Xu 2021, Nat. Commun. 12:6168 | 23 PWDs within 200 pc with Mg, Si, Ca, Fe | 11 fall in or next to the peridotite/pyroxenite ternary (Earth-like mantle assemblages); the other 12 are "exotic" (quartz pyroxenites, periclase dunites, etc.); only 1 matches bulk silicate Earth; no crust types found; sample is Ca-selected; Mg may be underestimated for cool WDs | Nature page (CC BY 4.0); data "available from the lead author upon request", Supplementary Table A1 | READ |
| Trierweiler, Doyle, PSJ 4:136 | 31 WDs with O plus other elements | chi-square test: 15 of 31 pass as chondritic on raw abundances (21 of 31 with steady-state correction); no WD better fit by MORB or continental crust; Hypatia stars: 75 percent chondritic (F 66.5, G 75.8, K 76.6, M 56.1 percent) | [R4425] | READ |
| Harrison, Bonsor, Kama, Buchan, Blouin | 202 DZ WDs | more than 60 percent consistent with primitive material; 8 systems core-rich at 3 sigma (4.0 percent); 1 crust-rich (0.5 percent); Na implies formation temperatures below 1000 K to above 1400 K | [R1832] abstract | READ |
| [R1833]:3814 | 17 systems | 1 to 9 core- or mantle-dominated; no crust-rich pollutants | [R1833] abstract | READ |
| [R0457]:2683 | Ca/Fe of WD sample | model fit: 66 +4/-6 percent accreted differentiated-planetesimal remnants, 31 +5/-5 percent altered by heating | [R0457] abstract | READ |
| [R1985], Nat. Astron. 5:451 | 524 northern WDs within 40 pc searched | Li in 4 cool WDs (K in one); Li, Na, Ca ratios imply crust fragments | [R1985] abstract | READ |
| [R4276]:202 | 9 WDs | volatile-poor, differentiated bodies; one mantle-like with small Fe-Ni core; some Ca- or Na-anomalous | [R4276] abstract | READ |
| [R3801] | 8 WDs with Spitzer dust | olivine/pyroxene dominance correlates tentatively with Mg/Si; GD 56 core mass fraction 0.59 | [R3801] abstract | READ |

**What can and cannot be anchored.** WD debris is disrupted bulk material, so it constrains: (i) mantle mineral assemblage, which the project already treats as a deterministic map of Mg/Si and C/O (T6); (ii) core-rich fraction about 4 percent of DZs at 3 sigma [R1832] and 66 percent differentiated in a model [R0457]; (iii) a crust-bearing frequency of order 0.5 to 1 percent (1/202; 4/524 of the 40 pc sample have Li crust signatures, but the denominator of polluted WDs in that sample is not given); (iv) volatile content (excess oxygen for a handful, 2ND-HAND). It does not constrain surface rock types (basaltic, felsic, anorthositic, evaporite, clay, sulfur, carbon, ice) of an intact planet, nor their weights.

**Verdict.** T20: (c) for the surface types as a drawn table; (b) only for three sub-statistics (core-rich/metal about 4 percent, crust-bearing about 1 percent, water-rich a "handful" with no count). Be blunt: T20's weights are a design table that cannot be tuned to any observed population of surfaces. Recommend deriving ground type by decision rules from the stored mantle, volatile, heat flux and tectonic state (as T20 already states), and testing the output against WD statistics only at the level of bulk mantle assemblage and the crust-bearing frequency. T22/R2: (c) for a drawn proportion table. Putirka and Xu and Trierweiler disagree (12 of 23 exotic versus most chondritic), so the spread of the "exotic mantle" fraction, 0 to about 50 percent, is the honest tolerance; no verified source overrides the other. Boss must request Putirka and Xu Supplementary Table A1 and the [R1984] compositions (no arXiv id found; MNRAS 477:93, [R1984], Crossref only).

**Test.** Simulated bulk-mantle normative mineralogy of bodies built from Hypatia-like host abundances versus the 23 PWDs (Putirka) and 31 WDs (Trierweiler chi-square on Si, Fe, Mg, Al, Ca, Ni, Cr); needs a WD selection model (Ca selection, Teff, settling) and a differentiation model.

### 2.5 Stellar Th, U, K, Eu (rows o1 to o4, T23)

| Ensemble | N | Measured | Scatter and trends | Access, licence | Status |
|---|---|---|---|---|---|
| Hypatia Catalog ([R1943]:54, [R1943]; database paper [R1942]) | front page: 14,085 stars, 104 elements and species, 383 catalogues, 452,599 measurements | element list includes K, Eu, Eu_II, Th, Th_II; no U | see below | open web API (v2.2.0, no key) at hypatiacatalog.com; page shows "Copyright 2026", no open licence stated; 2014 reduced table on VizieR J/AJ/148/54 (2685 rows) | READ (papers), DATA (API) |
| Hypatia counts by element via API scatter query (vs [Fe/H]) | Th I: 1 star; Th_II: 284 stars; Eu: 903; Eu_II: 4570; K: 1509; U: 0 | DATA, COMPUTED | | | DATA |
| Hypatia [Th_II/Fe] | 284 | median +0.04, sd 0.22; for |[Fe/H]|<0.3 (N=228) slope vs [Fe/H] -0.79, residual sd 0.14 | catalogue heterogeneity inflates sd; Botelho slope is -0.338 | | DATA, COMPUTED |
| Hypatia Th by disc (star-query disk flag) | thin 217, thick 33, unclassified 34 | thin [Th/H] mean -0.015 sd 0.142; thick -0.254 sd 0.306; [Th/Mg] thin -0.044 sd 0.153 (N=215), thick -0.063 sd 0.193 (N=32); [Th/Eu] thin -0.065 sd 0.197, thick -0.077 sd 0.292 | no thick-thin offset in Th/Mg, the planet-relevant ratio | | DATA, COMPUTED |
| Hypatia [K/Mg] | 1419 with both | mean +0.031, sd 0.224, 95 percent -0.44 to +0.46 dex ([R4468] with 179 stars: x1.13 solar, 95 percent x0.35 to x3.63, consistent) | NLTE-sensitive K I lines | | DATA, COMPUTED; READ (Unterborn) |
| Hypatia [Eu_II/Fe] | 4570; |[Fe/H]|<0.3: 3550 | median +0.10; slope -0.23; residual sd 0.17 | | | DATA, COMPUTED |
| [R0478, R0479], [R0479] | 53 thin-disc solar twins | Th II 4019 A, precision 0.025 dex; pristine [Th/H] super-solar, uniform dispersion 0.056 dex, +0.037 to +0.138 vs age; Th follows Eu, not Nd (60 percent s-process) | | arXiv | READ |
| [R4469] | 14 solar twins and analogues | Th 59 to 251 percent of solar (Sun low); solar log eps(Th) 0.09 | | arXiv | READ |
| [R3096]:3786, [R3096] | 276 disc stars, Th in 170 | first extended disc Th set; errors 0.07 to 0.16 dex; Th declines more steeply than Eu with [Fe/H]; GCE needs metallicity-dependent Th yield or a Th-rich source whose rate falls | | arXiv | READ (abstract) |
| [R0664] | Eu in 577 (HARPS, Delgado [R1062]) + 377 (Bensby/Battistini) stars | intrinsic [Eu/H] scatter 0.025 dex at fixed [alpha/H] and Teff; [Eu/alpha] anticorrelates with [alpha/H]; Th/Eu assumed 0 (Botelho) but Mishenina finds larger dispersion | | arXiv | READ |
| [R0309] | 79 solar twins, 30 elements to about 2 percent | C/O and Mg/Si homogeneous within 10 percent at solar metallicity; no Th or K reported in the summary | | | SNIPPET |
| [R4468] | Hypatia: 2040 FGK with Eu and Mg; 179 with K and Mg | no 40K/K data outside the Solar System; Earth 40K/K 0.1419 percent unknown elsewhere; doubling 40K/K lengthens degassing by about 1 Gyr | | arXiv | READ |

**Answer to question 4.** Catalogue: Hypatia (open API) for Th II, Eu, K; U is not measured in any disc star ([R1447] and the eight halo r-process-enhanced stars are the only U detections, per G2; HE 1523-0901 log eps(U) = -2.06, log eps(U/Th) = -0.86, READ). Scatter: thin-disc Th intrinsic 0.056 dex (Botelho, READ) up to 0.14 to 0.19 dex heterogeneous (Hypatia, DATA); Eu intrinsic 0.025 dex, observed 0.07 to 0.17; K/Mg 0.22 dex. Trends: [Th/Fe] falls with [Fe/H] (slope -0.34 in solar twins over +-0.1 dex, Botelho, via G2; -0.79 over +-0.3 dex in Hypatia, DATA), so rock-normalised Th/Mg shows no population offset (thin -0.044, thick -0.063). Age: Botelho gives +0.037 to +0.138 pristine [Th/H] across ages (READ); G2 reports -0.0121 dex per Gyr of stellar age for [Th/Si] (2ND-HAND). Halo, young (below 1 Gyr) and M-dwarf hosts: no Th data in what I reached.

**Verdicts.** o2: (b). A draw log f_Th ~ N(m_pop, s_pop) with thin-disc s from Botelho 0.056 and an upper bound s of 0.15 to 0.19 (Hypatia thin/thick Th/Mg) is tunable to those statistics (N=215 and 32). The repo bracket U[-0.23, +0.40] (59 to 251 percent) is a 14-star range; Hypatia's thick-disc [Th/Mg] -0.063 +- 0.19 sits inside it. Halo, young, M: Boss must supply Th data (none found). o1 stays (a) as stated.

### 2.6 Th/U derivation (row o4) and K/U (row o3)

**Th/U formula (nucleosynthesis, birth epoch, decay).** Let PR = (232Th/238U) production ratio (number). Published inputs: [R1447: Table 1] (READ, original sources 2ND-HAND) lists log(U/Th)_initial from -0.301 through -0.29, -0.256, -0.243 ([R1001], U/Th = 0.571) to -0.22, i.e. PR = 1.67 to 2.0; [R1001] abstract (Nature 435:1203, [R1001]) states prior estimates of the U/Th production ratio 0.4 to 0.7 and derives a precise value from meteorites plus halo stars (value 0.571 +0.037/-0.031 from a search summary, SNIPPET; consistent with Frebel's table). Chronometer: delta t = 21.8 [log(U/Th)_initial - log eps(U/Th)_now] ([R1447], READ).

For continuous r-process enrichment of the gas at constant rate for T Gyr before the planet's host star forms, then radioactive decay for the planet's age tau:

Th/U (number) = PR x F(T) x exp[(l238 - l232) tau], with F(T) = [(1 - exp(-l232 T))/l232] / [(1 - exp(-l238 T))/l238]; mass ratio = 0.9748 x number ratio. Decay constants l = ln2 / half-life: 232Th 14.05 Gyr (RECALLED; [R4468] lists 14), 238U 4.468 Gyr ([R4468] lists 4.47, READ), so l238 - l232 = 0.106 per Gyr. Real r-process input is episodic (neutron-star mergers), so T is an effective duration and the constant-rate form is a simplification (COMPUTED model, not a published fit).

COMPUTED table (PR = 1.751; mass Th/U):

| T (Gyr) | tau = 0 | 2 | 4.5 | 8 |
|---|---|---|---|---|
| 3 | 1.98 | 2.45 | 3.19 | 4.63 |
| 6.5 | 2.32 | 2.86 | 3.73 | 5.40 |
| 10 | 2.65 | 3.28 | 4.27 | 6.18 |
| 13 | 2.93 | 3.62 | 4.72 | 6.83 |

PR = 1.667 or 2.0 at T = 6.5 gives 3.55 or 4.26 at tau = 4.5. **What is READ and what is calibrated:** PR range (READ via Frebel Table 1) and decay constants are not Solar System. T is calibrated here on one datum, CI Th/U = 3.76 ([R2986] via G2): it requires T = 6.5 Gyr for PR 1.751 (7.7 for PR 1.667, 3.5 for PR 2.0), COMPUTED. So the degeneracy PR versus T is closed only by the Solar System: **flag.** An independent route would fit T(birth epoch) with a Galactic chemical-evolution mod[R3096] GCE with OMEGA+, READ abstract).

**Verdict o4.** Replace U[3.8, 4.5] by this formula; the repo range is an age-4.5 Gyr Solar System value and ignores that Th/U grows with age (about x1.16 over 1.5 Gyr, COMPUTED). Knobs: PR in [1.67, 2.0] (READ), T by birth epoch (calibrated), tau from the planet age. Test: halo-star chronometer ages at [Fe/H] about -3 (N=2 with U in the same pattern, 8 stars in total per G2) and the Solar System point only; no disc star has U. Conclusion: derive, do not draw; verdict (b) for tuning because the only comparison data are one Solar System datum and about 8 halo stars.

**o3 K/U.** K is moderately volatile; planetary K/U = (K/U)_star x d_K. Solar System d_K values (Moon 0.032, Earth 0.18, Mercury surface 0.18, Mars 0.28; COMPUTED in G2 from READ ratios) are the only direct constraints on d_K; flag: Solar System only. Extrasolar information: (i) [R1832]: Na in 202 DZ pollutants spans formation temperatures below 1000 K to above 1400 K (READ abstract); K and Na are both moderately volatile (50 percent condensation temperatures about 1006 K for K and 958 K for Na, [R2726], RECALLED), so a draw d_K = f(T_form) with T_form tuned to the Na/Mg spread of [R1832] is possible if Boss accepts the proxy; (ii) [R1985] detected K in 1 of 4 crust-bearing cool WDs (N=1); (iii) [R0046] finds moderate refractories (Mg, Si, Fe) depleted relative to Ca, Ti in one WD against its F-star companion (N=1). Verdict (b) via Na, otherwise (c). 40K/K is unconstrained outside the Solar System ([R4468], READ). [R0521] (Nature Geosci. 12:564, [R0521]; authors Braukmuller, Wombacher, Funk, Munker, not "Siebert et al."), Palme and O'Neill, Lodders: not read (NOT REACHED beyond metadata).

### 2.7 Al-26 (rows o6, T23, R4)

| Ensemble or constraint | Value | Source | Status |
|---|---|---|---|
| Canonical early-Solar-System 26Al/27Al | (5.23 +- 0.13)e-5 ([R2768] table); 5.2e-5 ([R2128] via [R2209]); 5.85e-5 ([R4353] via [R2679]). Disk deviations up to a factor 2 ([R2569, R2717] via [R2209]) | [R2209, R2679, R2768] | READ (tables), 2ND-HAND (originals); Solar System |
| Present-day interstellar average | 3.0e-6 [R4307] | via [R2209] | 2ND-HAND |
| Galactic 26Al mass | 2.8 +- 0.8 solar masses (INTEGRAL/SPI, 16 sigma line) | [R1108, R1109] A&A, [R1109] | READ |
| Star-forming regions today | molecular-cloud 26Al/27Al about 2 to 3e-5 if all interstellar 26Al sits in clouds; Orion about 3e-5 (5.8e-4 Msun 26Al in 2e5 Msun H2) | Jura, Xu, [R2209, R4837] | READ (COMPUTED by authors) |
| Extrasolar melting | WD parent bodies with Fe/Al varying over a factor above 100 require 26Al/27Al at least 3e-5 for melting of small rocks (formation time 1.5 Myr) | same | READ |
| WD differentiation to enriched fraction | at least 4 percent of systems enriched (conservative), most likely tens of percent; example 66 percent differentiated with 40 percent of un-enriched bodies differentiating gives 43 percent enriched: 25 percent of systems at Solar System level; [R1648]: about 1 percent above Solar System | [R1716], MNRAS, [R0967] | READ (Reiter, Gounelle 2ND-HAND) |
| Fraction differentiated pollutants | 66 +4/-6 percent [R0457]; 16 of 42 WDs (40 percent) core- or mantle-rich ([R0573] via Curry) | [R0457, R0967] | READ |
| Cluster simulations | fraction of disks enriched f_enr 9 to 28 percent; mean 26Al/27Al of enriched disks 3.5e-5 to 2.4e-4 (sd comparable to or above the mean); 3 to 10 percent of disks have 26Al above the Solar System; 7 to 16 percent have heating above Solar System | [R3050]:3979, [R2679: Table 4] | READ |
| Inheritance picture | Solar System 26Al explained by Wolf-Rayet-wind input to molecular clouds with cloud residence about 1e8 yr; no unusual nearby source needed. Title is "unexceptional", not "unique" | [R4837]:16, [R4837] | READ (abstract) |
| Consequence | high 26Al gives small dry planets; no 26Al gives ocean worlds, radii differ by up to about 10 percent | [R2678], Nat. Astron., [R2678] | READ (abstract) |
| [R0030] (ARAA birth environment) | | NOT REACHED (ARAA refused) | NOT REACHED |

**Answer to question 5.** The range 5.25e-6 to 5.25e-4 is 0.1 to 10 times the canonical Solar System value, so it is built on the Solar System. Without it: the low end is supported by the interstellar mean 3.0e-6 (2ND-HAND) and cloud values 2 to 3e-5 (READ); the high end by cluster-simulation means up to 2.4e-4 with large scatter (READ), with 3 to 10 percent of disks above Solar System; the WD data need at least 3e-5 in a tens-of-percent fraction of systems. The shape is not log-uniform: the evidence points to a mixture (un-enriched near the ISM level, enriched at tens of 1e-5 with a tail), and the distribution of levels for enriched disks is only a figure in [R2679] (numbers not in text). R4 replacement (derive from nucleosynthesis yields and birth cluster) is feasible: yields are in [R2679] (Mej,Al = 3.4e-6 M*, Mej,Fe = 5.6e-7 M*, READ) plus the Diehl galactic mass as a consistency check.

**Verdict o6: (b).** Draw an enriched flag with probability p_enr and a level from a mixture. Knobs and targets: p_enr 0.04 to about 0.66 (Curry; central "tens of percent" about 0.25 to 0.45); level-enriched centre about 3e-5 to 1e-4 (Jura, Lichtenberg); un-enriched level 3e-6 (ISM). The repo c[0.1,1] is inside that bracket; 0.1 is not contradicted by the lower bound 0.04. Test: simulated fraction of differentiated WD pollutants ([R0457], 66 +4/-6 percent) with a body-size model (Curry: differentiation by gravity needs bodies above about 1000 km if p_enr is below 50 percent); needs the WD selection and the polluter-size distribution. Only the Solar System fixes the exact 5.2e-5 point; extrasolar data fix the order of magnitude and the enriched fraction.

## 3. Verification of candidate sources

| Candidate (from prompt) | Result | Id / DOI | N | Access and licence | Status |
|---|---|---|---|---|---|
| [R1154] Science | exists | [R1154] | 6 WDs | arXiv; licence not stated | READ (method), values 2ND-HAND |
| [R1152] | right title "Where are the extrasolar Mercuries?", ApJ | [R1152] | 16 | arXiv | READ |
| (extra) [R1153] | new | [R1153] | 8 | arXiv | READ |
| [R1833] | exists, MNRAS 479:3814 | [R1833] | 17 | arXiv | READ |
| [R1832] | MNRAS 504:2853 is Harrison, Bonsor, Kama, Buchan, Blouin, Koester | [R1832] | 202 DZ | arXiv | READ. A second [R0458]:116694 exists per a search summary, not read | SNIPPET |
| [R0457] | exists (MNRAS 492:2683) | [R0457] | WD Ca/Fe sample | arXiv | READ |
| [R0458] | pilot study of a WD in a wide binary: not found; follow-up found | [R0046] | 1 | arXiv | READ (abstract) |
| [R4276] | exists | [R4276] | 9 | arXiv | READ (abstract) |
| Putirka and Xu 2021 | exists, Nat. Commun. 12:6168 | [R3627] | 23 | Nature, CC BY 4.0; data on request | READ |
| [R4425] | no paper found | | | | NOT REACHED |
| [R4425] | "A Chondritic Solar Neighborhood", PSJ 4:136 | [R4425] | 31 WDs, Hypatia stars | arXiv | READ |
| Xu et al. | Xu 2019 AJ 158:242 | [R4799] | about 10 | Crossref metadata only | NOT READ |
| [R1984] | MNRAS 477:93; no arXiv id found | [R1984] | | Crossref only | NOT REACHED |
| [R1579] | exists | [R1579] | 4 | arXiv | READ (abstract) |
| [R2210] | ARE&PS 42:45 | [R2210] | | not reachable | NOT REACHED |
| [R2415] | exists | [R2415] | 85 DAs | arXiv | READ (abstract) |
| [R1951] | GCA 313, "Iron-wustite revisited" | [R1951] | | metadata only | NOT READ |
| [R4154] | Sci. Adv. 6:eabd1387 | [R4154] | Earth | metadata only (ANU repository CC BY-NC 4.0 per search) | 2ND-HAND (G3) |
| [R1075] | right id; authors Deng, Du, Karki, Ghosh, Lee | [R1075] | 3 bodies | Nature | READ |
| [R1717] | exists | [R1717] | Hypatia planet hosts | arXiv | READ |
| [R1488] | EPSL 577:117255, authors Gaillard, Bernadou, Roskosz, Bouhifd | [R1488] | | metadata only | 2ND-HAND (G3) |
| [R0501] | PSJ 3, "Retention of water in terrestrial magma oceans and carbon-rich early atmospheres" | [R0501] | | metadata only | NOT READ |
| [R3944] | ApJ 843:120 | [R3944] | | values via part 05 | 2ND-HAND |
| [R0479] | exists | [R0479] | 53 | arXiv | READ |
| [R4469] | exists | [R4469] | 14 | arXiv | READ |
| [R0309] | exists | [R0309] | 79 | arXiv | SNIPPET |
| [R1943] | exists | [R1943]; API docs hypatiacatalog.com/api | 14,085 stars (site) | open API; copyright 2026 notice, no licence | READ, DATA |
| Delgado Mena (Eu, Th) | Eu in 2017 paper (A&A 606:A94); no Th paper by this author found | [R1062] | 577 Eu stars (via Carrasco) | arXiv | 2ND-HAND |
| Frebel and Hill cosmochronometry | right paper is Frebel, Christlieb, Norris, Thom, Beers | [R1447] | 1 star (+CS 31082-001 remeasured) | arXiv | READ |
| [R1001] | Nature 435:1203 | [R1001] | | abstract READ; value SNIPPET and via Frebel Table 1 | READ (abstract) |
| [R2679] | exists | [R2679] | N-body 1e3 to 1e4 stars, 30 runs each | arXiv | READ |
| [R2678] | exists (Nat. Astron.) | [R2678] | model | arXiv | READ (abstract) |
| [R2768] | Prog. Part. Nucl. Phys. 102:1 | [R2768] | review | arXiv | READ (table) |
| Diehl et al. | A&A paper | [R1109] | | arXiv | READ |
| [R0030] | | | | ARAA refused | NOT REACHED |
| [R4837] | title is "unexceptional" | [R4837] | | arXiv | READ (abstract, table) |
| (extra) Jura, Xu, [R4837] | key extrasolar 26Al paper | [R2209] | | arXiv | READ |
| (extra) [R0967] | WD-based SLR prevalence | [R0967] | 202-WD framework | arXiv | READ |
| (extra) [R0573] | masses of WD pollutants | [R0573] | 42 systems | arXiv | READ |
| (extra) [R0664, R3096, R4468] | stellar Th/Eu | [R3096]; 2411.10711; 2204.04243 | 276; 577+377; Hypatia | arXiv | READ |

## 4. Boss must supply

1. Putirka and Xu 2021 Supplementary Table A1 (23 PWD abundances, "on request" from the lead author). Without it the exotic-mantle fraction cannot be re-fit.
2. White-dwarf compositions, and Xu et al. 2019 (AJ 158:242) tables, for a WD sample with Ca, Fe, Mg, Si, Na, K.
3. Machine-readable tables for [R1154], 2020, 2023 (24 systems), plus [R3795] oxygen budgets, so the redox test can be run with per-star errors.
4. Hypatia licence and permission to redistribute derived tables (copyright notice only; K, Eu, Th are mixtures of literature catalogues).
5. A decision on which IW scale each row draws (core-formation dIW from WDs versus surface magma-ocean dIW) and a conversion test (G3 G-4, [R1951] full text).
6. [R1075: Fig. 3] data (surface dIW vs ocean-base pressure) and any observational source on Fe3+/sumFe in exoplanet mantles (none found).
7. A Th (or Th/Eu) survey of halo, young and M-dwarf hosts, and any U in disc stars (none found).
8. A Galactic chemical-evolution r-process history ([R1438] or equivalent) to replace the Solar-System-calibrated duration T in the Th/U formula.
9. Distribution of enriched 26Al/27Al in [R2679: Fig. 7] as numbers, [R0030, R1648].
10. Na, K condensation temperatures from a primary source ([R2726] or Wasson) to map WD Na depletion to K.

## 5. Recalled / to verify

| Item | Where used | Why unverified |
|---|---|---|
| half-life 232Th 14.05 Gyr | Th/U table ([R4468] lists 14, 238U 4.47, READ) | 14.05 from memory |
| atomic masses 232.038, 238.051 (mass/number factor 0.9748) | Th/U table | standard values, not looked up |
| 50 percent condensation temperatures K 1006 K, Na 958 K | o3 proxy | memory [R2726] |
| [R1001] value 0.571 +0.037/-0.031 | PR range | search summary; matches Frebel Table 1 (-0.243) |
| [R1648] (1 percent) | o6 | read only insi[R0967] |
| [R4307] ISM 3.0e-6 | o6 | read only insi[R2209] |
| [R0458]:116694 | candidate | exists per search summary, not read |
| [R1488, R4154] specifics | a1 | via part G3, not read by me |
| Venera 9, 10 K/U | not my group | RECALLED in G2 |
| [R1984] crust claim cited by Putirka and Xu as "27 of 29 PWDs" (their ref 11) | T22 | I could not open the study; the attribution is from reading their numbered reference list |
