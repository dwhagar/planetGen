# Gap-filling plan: how we fill what the public record does not hold

Status: research plan. Written 2026-10-17 by the gap-filling research and planning thread (read-only: no code, no TODO.md, no docs/plan edits; it finished at about 11:12 AM EDT as it reported) and filed by the TODO thread on 2026-10-17 at about 11:16 AM EDT. The stand-in registry of section 6 and the work order of section 9 go to the build items GEN.401, GEN.402, GEN.403, GEN.405, GEN.369, GEN.386, GEN.447, GEN.472 and GEN.503 (GEN.280 is already merged, PR #1855). Boss's request list is GEN.531. Boss decided the seven open points of section 8 between 11:19 AM and 11:23 AM EDT as the coordinator reported; section 8 records each as Decided and the rows they touch say so. Everything attributed to Gemini is unverified until a checking thread verifies each claim (section 7); nothing in this note is a default until a builder registers it.


Written 2026-10-17, about 11:12 AM EDT as this thread's clock reported, by the gap-filling research and planning thread. Read-only: no code, no PRs, no edits to TODO.md or docs/plan. TODO files it. Inputs: `gen-406-galaxy-research-2026-10-17.md` and `docs/design/galaxy-laws-research.md` (the 15 requests, GEN.531), `tunable-gaps-2026-10-17.md` (extended here, not redone), section 6 of the atmosphere research master, the body composition master, and sections 4a and 4b of `docs/design/composition-model.md`. Boss's first batch of Gemini results (pasted in this thread at 10:52 AM EDT) is checked in section 7 and labelled "Gemini, unverified" wherever it appears.

Scripts, data and logs behind every number marked COMPUTED or DATA: `/mnt/project-files/reports/gap-filling-plan-2026-10-17-files/`.

## 0. Summary

1. **Of the open gaps, most can be filled by reasoning from verified inputs, and several of Boss's 15 research requests do not need him or Gemini at all.** I list 57 gaps (section 2). By their best fill (section 5): 15 get a derivation (D) from verified relations, 28 a bracket from verified neighbours (B), 11 a theory prior (T), 3 stay labelled placeholders (P). Seven items needed a person or a data set no reasoning can replace; Boss decided all seven (section 8).
2. **Eight fills were tested in this thread, not just proposed** (section 4); the four that matter most: (a) the water share of solids beyond the snow line follows from the oxygen budget of Asplund 2021 abundances and gives 0.17 to 0.32 for a solar mix (0.40 at C/O 0.35), which reproduces the 0.14 to 0.39 bracket of the research draft and puts the 0.5 convention at the all-methane limit; across 1592 real stars (Brewer 2016) it ranges 0.21 to 0.51 (16 to 84 percent, all carbon as CO), so it is a per-star quantity and the "switch" can be replaced by a derived function (BC1); (b) an accreted-halo integral built from Moster 2013 and Giocoli 2008 lands inside the measured stellar-halo range at 5e10 with no tuning (0.071 against 0.056), but it over-accretes by 2 to 3 times over the full range, exceeds 1 above 1e11, and its mass slope disagrees with the D'Souza stack (1.4 to 2.2 against 0.79), which is the first test to run (GR1); (c) the age-velocity normalisation (R4) is answered for the Milky Way by Sharma and Hayden 2020 Table 2 (21.1 km/s at 10 Gyr, exponent 0.441) and scales to other galaxies through the verified Martinsson relation with a 20 percent closure check (GR4); (d) the Comerón 2021 hot-component mass gives a hard upper bound on the thick-to-thin mass ratio that settles the Yoachim against Comerón 2018 conflict (GR3).
3. **Reachability changed.** By 10:59 AM EDT VizieR returned table data, ExoMol and HITRAN front pages answered, and the arXiv search interface returned ids for papers the GEN.406 thread "could not find by title" (Muñoz 2018 III is arXiv:1806.06891, Cooper and Cole 2010 is arXiv:0910.3211, Jordan 2005 is astro-ph/0508219, Pfeffer 2018 is arXiv:1712.00019). So requests R9a (read, answered below), R10, R13 and the Brewer and Fischer tables are reading jobs for us, not research for Boss (section 3 classifies each request A, B or C).
4. **Gemini's first batch is mixed.** I checked 20 claims: 8 agree with verified sources, 7 are contradicted by a verified figure (the occupation mass, the occupation edges, the thick-disc formula, three age-velocity numbers, the nuclear cluster slope), and 5 are unsupported, misattributed or internally inconsistent. I count about 110 distinct claims, about 45 checkable, so a separate checking thread is worth starting (section 7). It must enter through the intake rules of section 7, never as defaults.
5. **One convention makes every fill replaceable:** a stand-in registry with a fill id, a kind (D, B, T, P), an evidence label, assumptions, knobs, a validity range, a falsifier test and the research item that replaces it (section 6). It extends the `stand-in` and `open` markers already Decided for the GEN.472 loader and the target table of GEN.503, so no redesign is needed when Gemini or a paper delivers a value.
6. **What blocks what.** The galaxy gaps feed the 8.1 phase 5 "Galaxy laws build" items (GEN.401, 402, 403, 405, 369, 386); two of them also feed GEN.301's own prerequisites (GR4 through GEN.447, GR9c through GEN.280). The body and atmosphere gaps belong to the composition headline (9.0) and the generation changes of 8.4, apart from the building blocks already back-scattered to 8.1 (section 9).

## 1. Labels used in this plan

**Fill kind** (what the replacement value is):

| Kind | Meaning | Ranked how |
|---|---|---|
| D | derived: a physical or scaling relation computed from verified inputs; the inputs are named | first choice |
| B | bracket: bounded by neighbouring verified quantities; the draw spans the bracket | second |
| T | theory prior: a published or textbook theory supplies the form; the numbers are not measured | third |
| P | placeholder: my own choice or a recalled figure, labelled so, with a bracket | last |

**Evidence label** (how sure each input is), unchanged from the research lanes: READ, DATA, COMPUTED, 2ND-HAND, RECALLED, SNIPPET, NOT REACHED. New for this plan: **LEAD** (a Gemini claim or any claim whose paper has not been opened by us; it can steer where to look and widen a bracket, nothing else) and **POINTER** (a paper id found by search and not yet read).

**Standing rule applied everywhere** (Boss, 3:02 AM EDT): verified and partly verified material is usable unless contradicted; a verified source wins; speculation never beats a verified figure; every theory-based value says so.

**Ids.** GR1 to GR15 are R1 to R15 of GEN.531 (galaxy; "GR" avoids a clash with the atmosphere research items AR1 to AR15 = GEN.480 to GEN.494 and the body items R1 to R5 = GEN.473 to GEN.476). GX1 to GX5 are galaxy gaps the research found that are not among the 15. BC1 to BC19 are body composition gaps. AT.x are atmosphere rows by their tunable-gaps id (for example AT.a1); AD1 to AD8 are atmosphere data and input gaps from section 6 of the atmosphere master. Mass and time conventions: 200c masses with h = 0.671 inside the generator; the papers use other definitions, so each fill that crosses a definition says so.

## 2. Every open gap in one table

"Gates" says what the gap holds up: **8.1** = an item of the 8.1 phase 5 "Galaxy laws build" (in the lane report `core1-5-lanes-to-8-1-2026-10-17.md`; they gate the 8.1 stamp); **301** = also named in GEN.301's own prerequisite line; **9.0** = the composition and atmosphere-colour headline; **8.4** = a generation change that alters stored values.

### 2.1 Galaxy laws (GEN.531 and the open points of 5.13.6)

| Id | Quantity | What is missing | Why the public record lacks it | Feeds | Gates |
|---|---|---|---|---|---|
| GR1 | Stellar halo mass against M*, 3e9 to 3e10 | an observed halo mass or accreted fraction | halos fainter than about 30 mag/arcsec² are below survey limits; only simulations reach down | GEN.401, 403, 386 | 8.1 |
| GR2 | Numbers behind three figures (Rodriguez-Gomez 2016 Fig 3, D'Souza 2014 Fig 16, Benitez-Llambay and Frenk 2020 Fig 10) | digitised curves | printed as plots only | GEN.401, 369 | 8.1 |
| GR3 | Thick-disc mass, length and height ratios as formulas | one formula that survives the Yoachim against Comerón conflict | Comerón 2018 prints ranges and Spearman coefficients, no fit equation | GEN.401 | 8.1 |
| GR4 | Normalisation of the age-velocity relation for other galaxies | sigma at a reference age against Vc | stellar ages exist for the Milky Way only | GEN.447, 401 | 8.1, 301 |
| GR5 | Scatter of log(Vc/sigma0) (Courteau 2007) | the RMS | not printed in the text; the Courteau data table is not on VizieR and not reachable (the Letter's samples are 154, 24, 54 and 11 galaxies, no 64; section 4.9) | GEN.401, 402 | 8.1 |
| GR6 | Oldest thick-disc age and the first-disc epoch | start of thick-disc formation against mass | few external spectroscopic thick-disc histories (S0s) | GEN.386 | 8.1 |
| GR7 | Radius ratio of metal-poor to metal-rich cluster systems | ratio against galaxy mass | papers disagree (more extended against equal within errors); no ratio printed | GEN.403 | 8.1 |
| GR8 | Cluster system radius for log M* 8.5 to 10.3 | measured radii | spirals and intermediate masses are poorly sampled | GEN.403 | 8.1 |
| GR9a | Cluster half-light radius against mass and environment | the law | (answered by Jordan 2005, section 4.6) | GEN.297, 403 | 8.1 |
| GR9b | Cluster ages against host mass | the law | only Milky Way two-track data; simulations (E-MOSAICS) hold it | GEN.403 | 8.1 |
| GR9c | Nuclear cluster and nuclear disc ages and dispersions | the law | review exists, ages unread | GEN.280 | 8.1, 301 |
| GR10 | Stellar halo age against halo mass, more than six haloes | the relation | six measured haloes; simulation papers print it in figures | GEN.401, 403 | 8.1 |
| GR11 | Dwarf stellar-to-halo relation and scatter for 1e8.5 to 1e10.5 | median and intrinsic scatter | semi-empirical fits stop at 1e10.5 (checked here); Fitts 0.5 dex against Nadler under 0.19 dex | GEN.369 | 8.1 |
| GR12 | Spread of the reionization redshift | per-galaxy z_re | patchy-reionization papers print maps, not a law | GEN.369, 403 | 8.1 |
| GR13 | Ultrafaint size relation and scatter | r_h against M_V below log M* 5.5 | survey papers, not yet read | GEN.369 | 8.1 |
| GR14 | Accreted outer envelope above 1e11 | a halo law | the Harmsen method is unsuited to bulge-dominated galaxies | GEN.401, 386 | 8.1 |
| GR15 | Width of the halo occupation transition | fitted width and curve | Fig 10 of Benitez-Llambay and Frenk is a plot; Graus 2019 not located | GEN.369 | 8.1 |
| GX1 | Halo mass given M* (conditional, not the inverse of the median) | a halo mass function | Sheth and Tormen or Tinker not read | GEN.400 follow-on, 401 | 8.1 |
| GX2 | Arm pitch against bulge fraction | the sign | Hart 2017 and Yu and Ho 2018 disagree | GEN.401 | 8.1 |
| GX3 | Nuclear cluster mass slope | choice between 0.48 (Neumayer) and 1.0 to 1.36 (Georgiev) | genuine disagreement between verified papers | GEN.402 | 8.1 |
| GX4 | Environment factor of the cluster count (0.3 dex, Local Volume against Virgo) | an input | the generator holds no environment | GEN.403 | 8.1 |
| GX5 | Halo metallicity baseline (-0.6 against -0.9) and the law outside 5e8 to 1e10 of halo mass | one defensible baseline | the only printed equation covers 5e8 to 1e10 | GEN.405, 369 | 8.1 |

### 2.2 Body composition (body master section 8, 14.4, 15.4, tunable-gaps 3.5)

| Id | Quantity | What is missing | Why the record lacks it | Feeds | Gates |
|---|---|---|---|---|---|
| BC1 | Water share of solids beyond the snow line; condensation temperatures | Lodders 2003 (switch 0.14 to 0.39 against about 0.5) | IOP full text refused | GEN.428, 468, 410 | 9.0 |
| BC2 | Stellar Mg, Si, O, C against [Fe/H] (the project stores only [Fe/H]) | trends and covariance | tables sat behind VizieR, reachable since 10:59 AM EDT | GEN.439, 440 | 9.0 |
| BC3 | Mantle mineral equations of state (originals) | Stixrude and Lithgow-Bertelloni | publisher refused; a 2ND-HAND copy exists and is usable | GEN.425 | 9.0 (not a blocker) |
| BC4 | Oxide and K/Th/U compilations | Palme and O'Neill, Rudnick and Gao and others | publisher refused | GEN.428 | 9.0 |
| BC5 | Giant interiors: core mass, diluted core, inflation efficiency | structure of exoplanet giants | not observable; Juno and Cassini are two planets | GEN.416, T12 | 9.0 |
| BC6 | Moon composition against distance and origin | weights and ramps | zero confirmed exomoons | T17, T18, GEN.457 | 9.0 |
| BC7 | Ice densities and the ice compaction law | Fray and Schmitt (density), Bierson and Nimmo | publisher refused | GEN.442, 457 | 9.0 |
| BC8 | Planet minus star scatter in Fe/Mg and the FeO against star and water | the joint distribution | few precise planet compositions (Plotnykov slope 0.94, +1.02 to -1.07) | GEN.428 | 9.0 |
| BC9 | Water equation of state at high pressure (superionic) | one table | French 2009 and AQUA not extracted | GEN.425 | 9.0 |
| BC10 | Albedo of weathered regolith and frost | lab to Bond conversion | Hapke 2012 not read | GEN.441, 412 | 9.0 |
| BC11 | Rogue-planet volatile stripping at ejection | any measurement | none exists | GEN.428 | 9.0 |
| BC12 | Vapour-pressure validity limits and the cold end | Fray and Schmitt 2009 per-fit windows | publisher refused | GEN.470, 497 | 9.0 |
| BC13 | Radiation weighting and tissue factors | ICRP 103 tables | publisher paywall | GEN.436, 444 | 9.0 |
| BC14 | Belt composition kernel (R1, GEN.473) | derived general distribution | measured only in one belt | GEN.473, 477 | 9.0 |
| BC15 | Rock types from bulk oxides (R2, GEN.474) | petrology from the stored mantle | Earth library only | GEN.474 | 9.0 |
| BC16 | Meteorite groups as outputs (R3, GEN.475) | composition-space map | Solar System groups only | GEN.475 | 9.0 |
| BC17 | Al-26 and Fe-60 from the birth environment (R4, GEN.476) | yields and enrichment odds | one measured system | GEN.476 | 9.0 |
| BC18 | Core light elements (T7) | not checked in tunable-gaps | cores are unobservable | GEN.428 | 9.0 |
| BC19 | White-dwarf spectral type against temperature (T3) | not checked in tunable-gaps | large ensembles exist but were not read | GEN.463 family | 9.0 |

### 2.3 Atmosphere (tunable-gaps classes S and X, the thin F rows, and the data inputs)

| Id | Quantity | What is missing | Why the record lacks it | Feeds | Gates |
|---|---|---|---|---|---|
| AT.a1, a2 (T8) | Mantle redox against mass; later-mantle law | an ensemble beyond 22 white dwarfs | exoplanet mantles are not measured; fits use three Solar System bodies | GEN.482, 498 | 9.0 |
| AT.b4 | Water-rich share | a formation model | published shares run 0 to 87 percent by model | GEN.481, 510 | 9.0 |
| AT.b5 | Sub-percent water fraction | any ensemble | none measures it | GEN.481 | 9.0 |
| AT.n7 | Haze yield | a population | degenerate with Kzz, metallicity and cloud cover | GEN.491, 519 | 9.0 |
| AT.k7 | Single-comet species (N2, O2, Ar, Kr, Xe, SO2) | a population | 67P is the only comet with them | GEN.488, 495 | 9.0 |
| AT.j2, j3, j4 | Mars-, Titan-, Pluto-like occurrence | any ensemble | no exoplanet analogue measured | GEN.487 | 9.0 |
| AT.m3, m4, m11 | N2O, PH3, NO2 and CH3Cl presence | any measurement | none on any exoplanet | GEN.490, 525 | 9.0 |
| AT.F | Thin F rows: b1 to b3, d1 to d5, e5, f7, h4, i6, k1 to k3, n6, o2, o3, o6 | more data | few constraints each (tunable-gaps 3) | GEN.509 to 524 | 9.0 |
| AD1 | Visible absorption of CH4, NH3, H2O in giants | band strengths | Karkoschka tables not extracted | GEN.422, atmosphere colour | 9.0 |
| AD2 | Mars dust and Titan haze optical spectra | n and k tables | tables unread | GEN.423 | 9.0 |
| AD3 | Infrared opacity tables | correlated-k grids | line data reachable now, not extracted | GEN.413, 414 | 9.0 |
| AD4 | Rayleigh originals (King factors) | depolarization of SO2, H2S | publisher refused | GEN.422 | 9.0 |
| AD5 | Energetic particle spectra and shower runs | measured spectra | several originals refused | GEN.436, 444 | 9.0 |
| AD6 | Hot-giant line lists and vapour pressures | TiO, VO, FeH lists | ExoMol front page reachable, data not pulled | GEN.485, 499 | 9.0 |
| AD7 | H2SO4-H2O and brine vapour pressure | Tabazadeh, Zeleznik | publishers refused | GEN.471, 491 | 9.0 |
| AD8 | Photochemistry rates (JPL 19-5) | falloff parameters | JPL pages refused | GEN.496, 501 | 9.0 |

Counts: 22 galaxy rows, 19 body rows, 16 atmosphere rows (AT.F stands for 16 thin rows, so it is one line here and 16 in the registry), 57 lines.

## 3. Which of Boss's 15 requests need Boss, and which are reading jobs for us (17 lines because R9 splits in three)

Class **A**: the paper id is found (POINTER) and the thing is reachable, so a builder or research thread reads it; no one needs to supply it. **B**: the paper was searched for and not located, or is behind a refused publisher; try a wider search first. **C**: no source exists in the public record; fill by derivation or bracket and flag. **Done** means this thread read the answer. Ids found by the arXiv search interface at 10:55 to 10:59 AM EDT (titles checked, papers not yet read except where stated).

| Request | Class | What to read or do | Pointer |
|---|---|---|---|
| GR1 observed halo, 3e9 to 3e10 | C | no survey exists; fill D (section 4.2) | none |
| GR2 three figure curves | A | extract vector paths from the PDFs (Rodriguez-Gomez 2016 is arXiv:1511.08804, downloaded and its abstract read; D'Souza 2014 is arXiv:1404.2123; Benitez-Llambay and Frenk is arXiv:2004.06124) | as listed |
| GR3 thick-disc formulas | A then C | the Comerón 2018 text (arXiv:1709.00457, downloaded) prints no fit equation; Comerón 2011, 2012, 2014 tables remain to try; the bound of section 4.4 settles the choice either way | arXiv:1709.00457 read |
| GR4 age-velocity normalisation | Done for the Milky Way, C for other galaxies | Sharma and Hayden 2020, arXiv:2004.06556, Table 2 read (section 4.3) | read |
| GR5 Courteau scatter | B | the paper text has no RMS; searched 2026-10-17 (section 4.9): the table is not on VizieR and the authors' file is refused; the Letter is ApJ 655, L21 | astro-ph/0610721 read |
| GR6 oldest thick disc, first-disc epoch | B | Pinna 2019 not located by author and topic search; Haywood 2013 and Bensby 2014 are Milky Way only; try a title search by words in the abstract | none found |
| GR7 blue against red radius ratio | A | Kartha 2014 (arXiv:1310.1979, downloaded, ratio not in the text lines searched); Forbes and Ponman 2012 (arXiv:1205.5315); Hargis and Rhode 2014 (arXiv:1409.6306) | listed |
| GR8 intermediate-mass cluster systems | A | Hargis and Rhode 2014 and Georgiev 2010 (arXiv:1004.2039) | 1409.6306 |
| GR9a cluster half-light radius | Done | Jordan 2005, astro-ph/0508219, abstract read (section 4.6) | read |
| GR9b cluster ages against host | A | Pfeffer 2018 (arXiv:1712.00019) and the E-MOSAICS age papers | listed |
| GR9c nuclear populations | A | Neumayer, Seth and Böker 2020, arXiv:2001.03626 | listed |
| GR10 halo age against halo mass | A | Cooper and Cole 2010, arXiv:0910.3211 | listed |
| GR11 dwarf stellar-to-halo relation | A | Garrison-Kimmel 2017 (arXiv:1701.03792); Behroozi 2019 fit range read (section 4.5) | listed |
| GR12 reionization spread | A | Busha and Alvarez 2010 (arXiv:0901.3553), Alvarez and Busha (arXiv:0812.3405), Katz 2020 (arXiv:1905.11414); the effect is second order (section 4.5) | listed |
| GR13 ultrafaint sizes | A | Muñoz 2018 III, arXiv:1806.06891 | listed |
| GR14 envelope above 1e11 | C | no method isolates an accreted halo there; bracket | none |
| GR15 occupation width | B | Graus 2019 not located; the width is derived from the verified edges (section 4.5) | none found |

So of 17 lines, 2 are done (GR4 for the Milky Way, GR9a), 10 are class A reading jobs, 3 are B searches (GR5, GR6, GR15), 2 are C (derive or bracket: GR1, GR14). No galaxy request needs Boss's own research; two need a data table that may not exist (GR5, GR6), and Boss decided the fills (section 8).

The same reachability change applies to the body and atmosphere lists: VizieR returned the Brewer and Fischer table (done, section 4.7); ExoMol and HITRAN front pages answered (AD1, AD3, AD6 become download jobs; Boss decided they are downloaded at build time, section 8 item 4); NIST and PubChem are reachable for BC7 and AD4; the publisher pages in BC3, BC4, BC12, BC13 and AD7 stay refused.

## 4. Fills tested in this thread

Each subsection says what was computed, from which inputs, and what would show it wrong. Nothing here is a measurement of a galaxy or a planet. Every input marked READ was opened in this thread; others carry the label they had in the research lanes.

### 4.1 Water share of solids beyond the snow line from the oxygen budget (BC1)

**The derivation (D).** Take the stored star's abundances. Oxygen goes first into the rock: MgO, SiO2, Al2O3, CaO, Na2O, and FeO for the oxidised fraction of the iron; the rest of the iron is metal or FeS (sulfur taken from the same pool); nickel is metal. Carbon is split into refractory organics (counted as solid carbon only, hydrogen ignored, so the organic mass is a lower bound), gas-phase CO (one oxygen per carbon), CH4 (no oxygen) and CO2 (two). Whatever oxygen is left is water. The water share is that water mass over water plus rock plus metal plus refractory carbon. This is bookkeeping, not disk chemistry, and it has no kinetics, hydrated silicates, nitrogen ices or hydrogen in the organics.

**Inputs.** Asplund, Amarsi and Grevesse 2021 (arXiv:2105.01661) Table 2, photospheric log eps: C 8.46, O 8.69, Na 6.22, Mg 7.55, Al 6.43, Si 7.51, S 7.12, Ca 6.30, Fe 7.46, Ni 6.20 (READ).

| Case (Asplund 2021 mix, C/O 0.589, Mg/Si 1.10) | Water share of solids |
|---|---|
| Fe metal and FeS, all carbon as CO gas | 0.227 |
| Fe half FeO, all carbon as CO | 0.193 |
| Fe all FeO, all carbon as CO | 0.168 |
| Fe metal, 30 percent of carbon refractory, rest CO | 0.323 |
| Fe metal, all carbon as CH4 | 0.543 |
| all carbon as CO2 | 0 (more oxygen demanded than exists; not a physical option) |

The research draft's bracket of 0.14 to 0.39 (Marboeuf-derived) is reproduced by CO-dominated carbon with iron metal to oxidised and 0 to about 30 percent refractory carbon (0.17 to 0.32; 0.40 at C/O 0.35). The "about 0.5" convention is the all-methane limit (0.543), which a solar C/O does not reach in a hot inner disc but can in a cold, high-pressure, reducing one. So the two values of the Decided "keep the switch" are the two ends of one derived function of the star and the carbon partition.

**Dependence on the star** (same bookkeeping, one variable changed): C/O 0.35 gives 0.396, 0.589 gives 0.227, 0.80 gives 0.000 (all CO; no oxygen is left for water, the carbon-rich regime); Mg 20 percent low or high moves it to 0.249 or 0.205.

**Across real stars (DATA).** The 1592 stars of Brewer et al. 2016 with C, O, Na, Mg, Al, Si, Ca, Fe, Ni (VizieR J/ApJS/225/32, downloaded 10:59 AM EDT; elements scaled by their own [X/H], sulfur following iron, Asplund solar zero points): all-CO case median 0.341, 16 to 84 percent 0.213 to 0.514, 2.6 percent of stars below 0.05; 30-percent-refractory case median 0.394, 16 to 84 percent 0.305 to 0.528, 0.8 percent below 0.05. So both switch values sit inside the spread among ordinary stars; a fixed share is wrong for any one of them by tens of percent.

**Tunable draw.** The carbon partition (refractory share, CO, CH4) as a Dirichlet with concentration chosen so the 16 to 84 percent range reproduces 0.14 to 0.39 for a solar star; the iron oxidation fraction uniform on 0 to 1 until the redox rows (AT.a1) tie it to the mantle. These two numbers are the same oxidation variables the redox rows need, so BC1, BC8 and AT.a1 should share them (section 9, builder X1).

**Wrong if.** The derived range misses 0.14 to 0.39 for a solar mix (it does not); the Lodders 2003 condensation sequence, once obtained, puts water below 150 K at a different mass fraction by more than 0.1; measured rock-to-ice ratios of outer-disc bodies (comet nuclei, small KBOs, rest of BC7) fall outside the star-by-star spread.

### 4.2 Accreted stellar halo from the subhalo mass function (GR1, GR10, GX5)

**The derivation (D).** The stars of a galaxy's accreted component come from satellites that merged. Let mu = m_v/M0 be a satellite's virial mass at accretion over the host's present mass. Giocoli, Tormen and van den Bosch 2008 (arXiv:0712.1563) Eq. 2 give the universal "unevolved" subhalo mass function: dN/dln(mu) = N0 x^-a exp(-6.283 x^3) with x = mu/a, a = 0.8, N0 = 0.21, independent of host mass (READ, section 3 and Fig. 2). The satellite's stars come from the stellar-to-halo relation: Moster, Naab and White 2013 at z = 0 (READ in the GEN.368 paper check; the parameters are in `intergalactic-and-halo-contents.md` row 2), used at the accretion mass (a stand-in for the relation at the accretion epoch). A satellite below the reionization gate is dark: the Benitez-Llambay and Frenk occupation erf (centre 10^9.4 Msun/h, width 0.25 dex placeholder). The accreted stellar mass is the sum over mu below a ceiling mu_max (Rodriguez-Gomez 2016, abstract: ex situ stars are laid down at larger radii in order of decreasing merger mass ratio, so the halo takes the small mergers).

**Check on my implementation.** The host halo masses it returns (1.65e11, 3.15e11, 8.27e11, 2.00e12, 1.07e13 Msun for M* = 3e9, 1e10, 3e10, 5e10, 1e11) reproduce the worked table of row 2 of the GEN.368 laws table to its printed digits.

**Result (COMPUTED, `acc.py`).** Accreted stars over the host's stellar mass:

| Host M* | mu_max 0.1 | mu_max 0.3 | mu_max 1.0 | GEN.406 placeholder f = 0.05 (M*/5e10)^0.79 |
|---|---|---|---|---|
| 3e9 | 0.003 | 0.017 | 0.030 | 0.0054 |
| 1e10 | 0.005 | 0.024 | 0.041 | 0.014 |
| 3e10 | 0.016 | 0.071 | 0.106 | 0.033 |
| 5e10 | 0.071 | 0.230 | 0.289 | 0.050 |
| 1e11 | 0.71 | 1.05 | 1.12 | 0.087 |
| 2e11 | 2.57 | 2.93 | 3.00 | 0.15 |

Reading it: (1) at 5e10 with mu_max = 0.1 the integral gives 0.071 against the measured 0.056 (range 0.020 to 0.142, six discs of 3.8e10 to 8e10, Harmsen 2017, READ in the GEN.406 note), with no tuning; (2) the all-mu total at 5e10 (0.29) is about three times the Illustris ex situ fraction for Milky Way-sized galaxies (about 0.10, 1-sigma range 0.05 to 0.30; Rodriguez-Gomez 2016 abstract, READ), so the integral over-accretes by a factor of about 2 to 3 even before the halo ceiling; (3) above 1e11 it exceeds 1, which is impossible; the Moster z = 0 relation applied at accretion is wrong there, so the derivation is valid for M* up to about 5e10 and GR14 takes over above; (4) below 3e10 the derived halo is two to three times lighter than the placeholder; (5) the log slope of the derived fraction over 1e10 to 1e11 is 1.4 to 2.2 against 0.79 from the D'Souza stack (which measures an outer light fraction of disc-dominated galaxies, not the same quantity). Items 2 and 5 are the real risk: the derived mass dependence is too steep or the ceiling is wrong.

**Proposed use.** Median by the placeholder power law anchored to the measured six (B); the integral as the shape source below 3e10 and the consistency test above (D); the draw width 0.33 dex (Harmsen scatter, D) widened by half the gap between the two medians below 3e10. The same integral, with each progenitor's infall time drawn from the host's assembly history, gives the halo's age (GR10) and the progenitor-mass-weighted metallicity rule of D'Souza and Bell 2018 (GX5): [Fe/H]_halo = relation(M*_progenitor) - 0.15 (READ in the GEN.406 note).

**Wrong if.** The six Harmsen fractions fall outside 2 sigma of the median at their masses; the integral's slope cannot be brought within 0.3 of 0.79 by a physically motivated ceiling or deposit fraction; Dragonfly's mean 0.009 with 3 of 8 undetected (Merritt 2016, partly) sits outside the 16 to 84 percent range at 1e10 to 1e11.

### 4.3 Age-velocity normalisation for any disc (GR4)

**Milky Way anchor (READ).** Sharma and Hayden 2020 (arXiv:2004.06556, GALAH with Gaia DR2, LAMOST, APOGEE, Kepler and K2) Table 2: vertical dispersion sigma_0,z = 21.1 ± 0.2 km/s for midplane stars of solar metallicity, solar angular momentum and age 10 Gyr, exponent beta_z = 0.441 ± 0.007; radial sigma_0,R = 39.4 ± 0.3, beta_R = 0.251 ± 0.006. Dispersions rise linearly with height |z|; low-alpha and high-alpha stars follow the same relations; valid for 3 to 20 kpc (abstract). So sigma_z/sigma_R = 0.54 at 10 Gyr, inside the 0.5 to 0.6 of Aumer, Binney and Schönrich 2016 (READ in the GEN.406 note).

**Rule for other galaxies (D).** Do not copy the Milky Way values. Use the verified Martinsson et al. 2013 relation, sigma_z,0 = (0.248 ± 0.038) V(2.2 h_R), falling as exp(-R / 2 h_R) (READ), as the age-averaged dispersion at radius R, and the verified exponent beta_z for the shape in age. For a thin disc whose star formation was flat over 0 to tau Gyr, the age average of (t/tau)^beta is 1/(1+beta), so sigma_z(t, R) = sigma_bar(R) (1+beta) (t/tau)^beta. The Sharma and Hayden normalisation then becomes the independent check instead of an input.

**Closure check (COMPUTED from READ values; Milky Way numbers are check targets only).** For V = 220 km/s, h_R = 2.6 kpc, R = 8 kpc: Martinsson gives sigma_bar = 0.248 x 220 x exp(-8/5.2) = 11.7 km/s, so sigma_z(10 Gyr) = 16.9 km/s. Sharma and Hayden give 21.1. The rule is 20 percent low, inside the 15 percent scatter of the ratio plus the 0.441 against 0.5 exponent choice.

**Tunable draw.** beta_z drawn from N(0.441, 0.03); the Martinsson ratio drawn with its 15 percent scatter; sigma_R from beta_R = 0.251 and the ratio 0.54 ± 0.05.

**Wrong if.** The Milky Way closure misses by more than 35 percent after the flat-beyond-1.5-h_R correction (Martinsson section 8.3) is applied; a second galaxy with stellar ages disagrees by more than the scatter (none is known; this is what class C means for other galaxies).

### 4.4 Thick-to-thin mass ratio bounded by the hot component (GR3)

Comerón et al. 2021 (arXiv:2101.04478) Eq. 1 gives the mass of the hot component (thick disc plus classical bulge) as log M_h = 0.489 + 0.9161 log M* (READ in the GEN.406 note), so M_h/M* = 0.49, 0.39, 0.35 at M* = 3e9, 5e10, 2e11. With the thin disc taken as M* - M_h at most and the thick disc as M_h minus the classical bulge B (as a share of M*):

M_thick / M_thin = (M_h/M* - B/T_classical) / (1 - M_h/M*), bounded above by M_h / (M* - M_h).

| M* | Upper bound (no classical bulge) | with B/T 0.1 | with B/T 0.3 | Yoachim and Dalcanton formula | Gemini formula (LEAD) |
|---|---|---|---|---|---|
| 3e9 (Vc about 89) | 0.96 | 0.76 | 0.37 | 0.69 | 2.54 |
| 5e10 (Vc about 181) | 0.64 | 0.48 | 0.15 | 0.14 | 1.17 |
| 2e11 (Vc about 255) | 0.54 | 0.38 | 0.08 | 0.06 | 0.80 |

Both verified sources are consistent with the bound: Yoachim's low values need a classical bulge share of about 0.3 (early-type spirals, B/T 0.35 to 0.38 in row 6 of the laws table) and Comerón 2018's 0.2 to 1 above 120 km/s needs B/T up to about 0.2. The conflict is therefore morphology: the same bound with the drawn B/T of the galaxy gives both. The Gemini formula exceeds the upper bound at all three masses.

**Tunable draw.** B/T_classical from the existing row 6 law; ratio scatter 0.2 dex until a fit with scatter is read. Length ratio 1.15 (Vc/100)^-0.2 and height ratio 2.5 with 30 percent (Yoachim and Dalcanton, READ).

**Wrong if.** The drawn ratios violate the bound in any seed (this is a test, not a tolerance); for galaxies with a known B/T the ratio falls outside Comerón 2018's range.

### 4.5 Reionization and the occupation transition (GR12, GR15)

Benitez-Llambay and Frenk 2020 (READ in the GEN.406 note): half occupation at 10^9.4 Msun/h for early reionization (z_re above 10), a factor of 4 lower for z_re = 4, no galaxy below 3e8 and every halo hosting one above 5e9 (Msun/h); insensitive for z_re above 10 and below 2. The shift is 0.6 dex over 6 units of z_re, 0.1 dex per unit. The global mean is z_re = 7.68 ± 0.79 (Planck 2018 Eq. 18, TT,TE,EE+lowE, mid-point of a tanh model, and z_re below 9.0 at 95 percent; READ, arXiv:1807.06209). That puts the half-occupation mass at about 10^9.17 Msun/h (2.2e9 Msun with h = 0.671) for the mean, and a patchy spread of one unit in z_re moves it by 0.1 dex, which adds in quadrature to the transition width 0.25 dex as sqrt(0.25^2 + 0.1^2) = 0.27. The width itself follows from the two verified edges: from 3e8 to 5e9 is 1.22 dex, and if the edges are 2.4 to 3 sigma from the centre the width is 0.20 to 0.25 dex. So R12 is second order and R15's width is a derived bracket; the centre is the real unknown, and the centre depends on z_re, which the Planck mean fixes.

### 4.6 Cluster half-light radius (GR9a is answered)

Jordan et al. 2005 (astro-ph/0508219, abstract READ): half-light radii measured for thousands of clusters in 100 early-type Virgo galaxies; the average increases with galactocentric distance or with decreasing galaxy surface brightness, decreases with the host galaxy's colour, shows no variation with the cluster's luminosity, and red clusters are about 17 percent smaller than blue; the mean is 2.7 ± 0.35 pc for clusters of colour (g-z) = 1.2 in a host of colour 1.5 at surface brightness mu_z = 21. The Harris file (157 clusters, median 2.9 pc, GEN.406) agrees as a check. Fill: log-normal about 2.7 pc with the three dependences as small multipliers (coefficients to be extracted from the paper by the cluster builder), the red/blue offset 0.83, and no dependence on cluster mass.

### 4.7 Star abundances: Mg, Si, O, C against [Fe/H] (BC2 is answered)

From the 1592 stars of Brewer et al. 2016 with all four elements (DATA, VizieR J/ApJS/225/32, script `brewer_fit.py`; [X/Fe] = a + b [Fe/H], residual sd, zero points on the survey's own scale, number ratios with Asplund 2021 solar values):

| Relation | Slope b | Intercept a | Residual sd (dex) |
|---|---|---|---|
| [C/Fe] | -0.275 | -0.054 | 0.111 |
| [O/Fe] | -0.446 | 0.057 | 0.127 |
| [Mg/Fe] | -0.247 | -0.024 | 0.060 |
| [Si/Fe] | -0.196 | -0.027 | 0.087 |

Residual correlations: C-O 0.52, C-Mg 0.61, C-Si 0.66, O-Mg 0.32, O-Si 0.47, Mg-Si 0.69. Number ratios: C/O median 0.49, 16 to 84 percent 0.37 to 0.58; Mg/Si median 1.07, 16 to 84 percent 0.98 to 1.29. Caution: with no quality cut 0.9 percent of the catalogue lies above C/O 0.8 (maximum 1.00), whereas the paper's cut sample (N = 849 in the tunable-gaps note) has none above 0.66; the cut (Teff and line quality) must be applied before use, which the builder takes from the paper. Validity: [Fe/H] -1.2 to +0.6, F, G and K stars; outside it the draw is a flagged extrapolation, and thick-disc and halo stars need another sample (Hypatia, licence open).

**Fill (V, DATA).** Draw [Fe/H] from the stored star, then a four-vector of residuals from a multivariate normal with these sds and correlations. This replaces the need to store Mg, Si, O and C and keeps them reproducible from the seed. It also gives BC1 and BC4 their inputs.

### 4.8 Mantle redox from an iron-wuestite identity (AT.a1, a2)

For equilibrium between metal and silicate, dIW = 2 log10(a_FeO / a_Fe); in the ideal-solution limit a is the mole fraction. With FeO at 5 to 15 mole percent of the silicate and iron at 85 to 95 percent of the metal, dIW runs from -2.5 to -1.6; real melts add an activity-coefficient correction of a few tenths of a dex. The white-dwarf ensemble is -1.25 ± 0.10 (22 values, COMPUTED and re-run in the tunable-gaps note), so the identity gets the sign and the order (a more iron-oxidised mantle sits higher) and the ensemble sets the zero point within the correction. This is the derivation the "surface from ocean base pressure" row (AT.a1) can use in place of the 1.25 dex-per-dex slope interpolated between a floor and Earth. The iron oxidation fraction is the same variable as in section 4.1.

### 4.9 VizieR search for GR5 and GR3 (2026-10-17)

Boss's decision 6 (section 8) ran the quick VizieR search first. Full report: [vizier-courteau-comeron.md](vizier-courteau-comeron.md). Short:

- GR5: the Courteau data table is not on VizieR and the authors' file is refused, so the residual scatter about Eq. 5 cannot be computed. The raw scatter of log10(Vc/sigma0) of the spirals in two source samples is 0.114 dex (Prugniel, Maubon and Simien 2001, J/A+A/366/68, 69 spirals), 0.119 dex (Pizzella et al. 2005 high-surface-brightness, 40) and 0.123 dex (113 de-duplicated), and 0.104 dex about Hubble-type means. That is inside or within 0.003 dex of the top of the plan's 0.06 to 0.12 bracket, so the plan's replacement rule is not triggered and the fill stays; the true scatter about Eq. 5 is smaller. The paper is ApJ 655, L21 (not ApJ 671, 203).
- GR3 rank 2: Comeron, Salo and Knapen 2018 (A&A 610, A5, VizieR J/A+A/610/A5, 124 galaxies with at least two discs) gives log10(M_thick/M_thin) = 2.113 - 1.086 log10(Vc) with scatter 0.27 dex (0.32 with the single-disc galaxies; the plan had 0.2), a median thick-to-thin height ratio of 4.2 (the plan had 2.5) and a median length ratio of 1.46 (the plan's Yoachim and Dalcanton ratio falls with Vc, this one rises). At Vc 89, 181 and 255 km/s the fit gives 0.99, 0.46 and 0.32 against the plan's hot-component bound 0.98, 0.64 and 0.53, so it lies near the bound and well above the Yoachim and Dalcanton formula at high mass. Taken galaxy by galaxy, 48 of 124 (39 percent) exceed the bound, so the bound test applies to the mean relation or carries a margin.
- References: A&A 562, A121 is the ARRAKIS ring atlas; the 2014 thick-disc paper is A&A 571, A58 (arXiv:1409.0466). The Comeron 2011 and 2012 tables are not on VizieR and are superseded by 2018.

## 5. Ranked fill for every gap

Reading the tables: **Kind** is the kind of the first-ranked fill (D, B, T, P, section 1). **Fill** gives the ranked options in order; the first one is used unless its inputs are missing. **Draw** is how the uncertainty enters generation as a tunable draw (distribution and width); a width marked (P) is my own choice, labelled so, and a width marked (V) comes from a verified scatter. **Wrong if** is the test that would show the fill false; each one becomes the falsifier test of the registry (section 6). **Result lands as** says what a Gemini or paper result changes. Where the draw says "none", the quantity is computed, not drawn, and the scatter comes from its inputs.

### 5.1 Galaxy laws

| Id | Kind | Fill (ranked) | Draw | Wrong if | Result lands as |
|---|---|---|---|---|---|
| GR1 | B | 1 Placeholder median of GEN.406 anchored on the six measured haloes (Harmsen 2017). 2 The accreted-halo integral of 4.2 as shape source below 3e10 and as a consistency test. | log-normal about the median, 0.33 dex (V, Harmsen scatter), widened by half the gap between the two medians below 3e10 | a measured halo falls outside 2 sigma; the integral slope cannot get within 0.3 of 0.79 | a simulation or survey fraction replaces rank 1; Gemini's "ex situ below 0.10 under 1e10.5" (LEAD) is only an upper cap |
| GR2 | B | 1 The anchor values already read from the three figures (GEN.406). 2 Digitise the vector paths of arXiv:1511.08804, arXiv:1404.2123 and arXiv:2004.06124 (reading job). | interpolation error 0.1 dex (P) until digitised | digitised curve differs from rank 1 by more than 0.1 dex | digitised table replaces; Gemini's 8-row table (LEAD) is compared with it, not used |
| GR3 | D | 1 The hot-component bound of 4.4 with the drawn classical bulge fraction. 2 The Comeron, Salo and Knapen 2018 fit (VizieR J/A+A/610/A5, section 4.9): log10(M_thick/M_thin) = 2.113 - 1.086 log10(Vc), scatter 0.27 dex, median height ratio 4.2, median length ratio 1.46 (outermost section); Yoachim and Dalcanton length 1.15 (Vc/100)^-0.2 and height 2.5 ± 30 percent (READ) stay as the length and height fills until the section choice of the length ratio is decided. | ratio scatter 0.27 dex (V, Comeron 2018); B/T from the existing row 6 law; Decided (Boss, 2026-10-17, between 11:19 AM and 11:23 AM EDT, card answers relayed in words by the coordinator): ship this fill and run the quick VizieR search for the Comeron 2011, 2012 and 2014 tables first | 39 percent of the real galaxies exceed the bound galaxy by galaxy (section 4.9), so the test applies to the mean relation or carries a margin of 0.13 dex, and does not reject each seed | a published fit equation with scatter replaces the scatter only; Gemini's formula (LEAD) is rejected by the bound |
| GR4 | D | 1 Martinsson 2013 age-averaged sigma_z times (1+beta)(t/tau)^beta, beta from Sharma and Hayden 2020 (4.3). 2 The Sharma and Hayden Table 2 values as a closure target for the Milky Way only. | beta_z N(0.441, 0.03); Martinsson ratio 15 percent (V); sigma_R/sigma_z 0.54 ± 0.05 | the Milky Way closure misses by more than 35 percent | a second galaxy with ages tests the rule; Gemini's 38.5 and 62.0 km/s (LEAD) are contradicted and not used |
| GR5 | B | 1 Bracket from Martinsson's 15 percent and Gemini's 0.12 dex (LEAD may only widen), so 0.06 to 0.12 dex centred on 0.09. 2 The Courteau data file if found. Raw scatter of the source samples is 0.104 to 0.123 dex (VizieR search, section 4.9), inside or at the top of the bracket, so the fill is unchanged. | log-normal 0.09 dex, bracket as width (P); Decided (Boss, 2026-10-17, between 11:19 AM and 11:23 AM EDT, card answers relayed in words by the coordinator): ship this fill and run the quick VizieR search for Courteau 2007 first | the table gives an RMS outside 0.06 to 0.12 | the RMS from the table or the authors replaces the centre |
| GR6 | T | 1 Thick-disc start set by the gas-rich turbulent collapse picture: the host's halo formation epoch (GR10 integral) plus a delay. 2 A bracket between the Milky Way thick-disc age (Haywood, Bensby, Milky Way only) and the S0 histories if located. | start lookback: GR10 epoch, delay uniform 0 to 1 Gyr (P); Decided (Boss, 2026-10-17, between 11:19 AM and 11:23 AM EDT, card answers relayed in words by the coordinator): halo formation epoch plus a 0 to 1 Gyr delay, a labelled placeholder, replaced when S0 histories are found | S0 spectroscopy shows a start independent of mass while GR10 gives a mass dependence | external-galaxy histories replace the delay; Gemini's 12 to 13 Gyr start (LEAD, section 7 claim 19) cannot set the epoch because its redshifts do not match its lookback times |
| GR7 | B | 1 Lower end 1.0 (the papers that find equal extent within errors, READ in GEN.406); upper end the largest printed ratio among Kartha 2014, Forbes and Ponman 2012, Hargis and Rhode 2014 (to read). 2 Until read, upper end 1.5 (P). | log-uniform between the ends; Decided (Boss, 2026-10-17, between 11:19 AM and 11:23 AM EDT, card answers relayed in words by the coordinator): widen and keep both, neither source dropped | a fourth survey puts the ratio outside the bracket | printed ratios replace the ends; Gemini's 1.20 ± 0.08 (LEAD, attribution unverified) may centre the draw after one of the three papers shows it |
| GR8 | B | 1 Interpolate the radius of the cluster system in log against log M* between the verified low-mass and massive ends. 2 Hargis and Rhode 2014 and Georgiev 2010 measurements (reading job). | 0.3 dex (P) | measured radii sit outside the interpolation by more than 0.3 dex | measured radii replace the interpolation; Gemini's ranges (2 to 8, 10 to 25, over 100 kpc) are LEAD |
| GR9a | D | 1 Jordan 2005: mean 2.7 ± 0.35 pc, red 17 percent smaller than blue, small dependences on distance, surface brightness and host colour, none on cluster luminosity (abstract READ). 2 The Harris file median 2.9 pc as a check. | log-normal about 2.7 pc; width from the paper's table (to extract); red/blue 0.83 | the 157 Harris clusters disagree with the median by more than the paper's scatter | coefficients of the three dependences replace the multipliers |
| GR9b | B | 1 Old track ages from the infall times of the accreted satellites (GR10 integral), young track from the disc epoch (GR6), weights from the blue/red split. 2 E-MOSAICS ages (reading job, arXiv:1712.00019 and companions). | 1 Gyr (P) per track | simulated ages differ from the rule by more than 2 Gyr | simulation table replaces the rule |
| GR9c | B | 1 Dispersion from the virial theorem with nuclear cluster mass and half-light radius (D); ages as an old base from the host epoch plus recurrent young populations for M* above 1e9. 2 Neumayer, Seth and Böker 2020 (arXiv:2001.03626, reading job). | young mass fraction uniform 0 to 0.3 above 1e9, zero below (P) | the review shows a young fraction outside 0 to 0.3 or dispersions far from the virial value | review values replace the knobs; Gemini's 15 to 50 km/s and young below 100 Myr (LEAD) are consistent and widen nothing |
| GR10 | D | 1 Infall times from the accreted-halo integral (4.2) with the host's own assembly history. 2 Cooper and Cole 2010 (arXiv:0910.3211, reading job) as the check. | width from the spread of infall times in the integral | Cooper and Cole's ages differ by more than 1.5 Gyr at the same halo mass | their ages replace the integral; Gemini's z50 formula (LEAD, claim 20) is compared |
| GR11 | B | 1 Moster 2013 extended below 1e10.5 and flagged as an extrapolation (Behroozi 2019 fit range READ: only 10^10.5 to 10^15). 2 Garrison-Kimmel 2017 (arXiv:1701.03792, reading job). | scatter uniform 0.19 to 0.5 dex (bracket of Nadler and Fitts, READ in GEN.406) | the median at 1e9.5 differs from the dwarf-galaxy data by more than the scatter | a dwarf fit replaces rank 1; Gemini's slope 1.90 and zero point (LEAD, claim 18) are not used as a median |
| GR12 | D | 1 Occupation centre shift 0.1 dex per unit z_re from Benitez-Llambay and Frenk and z_re = 7.68 ± 0.79 (Planck 2018) (4.5). | z_re from N(7.68, 0.79) | the patchy-reionization maps give a spread more than three times 0.79 | Busha and Alvarez or Katz widths replace 0.79; Gemini's void and overdense ranges (LEAD) widen only after the maps are read |
| GR13 | B | 1 Continue the verified dwarf size law downward to M_V near -2, floor at the cluster size 2.7 pc. 2 Muñoz 2018 (arXiv:1806.06891, reading job). | log-normal 0.3 dex (P) | measured r_h sit outside the continuation by more than 0.3 dex | the paper's relation replaces; Gemini's mean 80 pc, 20 to 300 pc, 0.35 dex (LEAD) is a first comparison |
| GR14 | B | 1 Logistic rise of the ex situ fraction in log M* from the verified 0.1 at 5e10 to a ceiling of 1 at the top, centre fixed by the digitised Illustris curve (GR2). 2 Until then, centre between 5e10 and 1e12. | centre 0.25 dex (P) | the digitised curve crosses 0.5 outside the bracket | digitised curve replaces; Gemini's crossover 1.5e11 ± 0.2e11 (LEAD) tests the result |
| GR15 | D | 1 Width from the verified edges: 3e8 to 5e9 is 1.22 dex, so 0.20 to 0.25 dex (4.5). | 0.25 dex plus the 0.1 dex z_re term in quadrature = 0.27 | a fit to the Benitez-Llambay and Frenk simulations gives a width outside 0.15 to 0.35 | the fitted width replaces; Gemini's 0.25 ± 0.04 agrees (claim not counted: its centre is contradicted) |
| GX1 | D | 1 Bayes: the probability of a halo mass given the stellar mass is proportional to the probability of the stellar mass given the halo mass times the halo mass function, with the scatter of the Moster relation. 2 Sheth and Tormen or Tinker function (reading job). | none (the conditional itself) | the sampled halo masses at fixed M* have a median shifted from the inverse of the median relation by something other than the Eddington-bias term | a published conditional replaces |
| GX2 | B | 1 Mixture of the two published signs (Hart 2017, Yu and Ho 2018) with equal weights; neither wins because both are verified. | sign 50/50, magnitude from each paper; Decided (Boss, 2026-10-17, between 11:19 AM and 11:23 AM EDT, card answers relayed in words by the coordinator): widen and keep both, neither source dropped | a third survey fixes the sign | the surviving sign replaces the mixture |
| GX3 | B | 1 Slope uniform on 0.48 to 1.36 (Neumayer against Georgiev). Gemini's 0.65 ± 0.10 sits inside and changes nothing. | uniform 0.48 to 1.36; Decided (Boss, 2026-10-17, between 11:19 AM and 11:23 AM EDT, card answers relayed in words by the coordinator): widen and keep both | the review shows one fit to be a selection artefact | a reconciling paper replaces |
| GX4 | T | 1 Environment drawn from the halo mass (more massive halos are in groups and clusters); the 0.3 dex factor applies with that probability. | probability from the halo mass function; Decided (Boss, 2026-10-17, between 11:19 AM and 11:23 AM EDT, card answers relayed in words by the coordinator): the environment is drawn from the halo mass and no environment concept is added now, so the 0.3 dex cluster-count factor applies with the group probability (also GR11's scatter choice and AT.j) | the generator gains an environment concept that disagrees | a later environment concept replaces the draw |
| GX5 | D | 1 [Fe/H] of the halo = relation(M* of the dominant progenitor) minus 0.15 (D'Souza and Bell 2018, READ), progenitors from the 4.2 integral. 2 Outside 5e8 to 1e10 of halo mass, extend flagged. | 0.2 dex (P) around the progenitor relation | the six measured haloes lie outside 2 sigma | a metallicity-mass paper for dwarfs replaces the extension |

### 5.2 Body composition

| Id | Kind | Fill (ranked) | Draw | Wrong if | Result lands as |
|---|---|---|---|---|---|
| BC1 | D | 1 The oxygen budget of 4.1 per star. 2 The two end values 0.14 to 0.39 and 0.5 as the bracket check. | Dirichlet over carbon partition; iron oxidation uniform 0 to 1 until tied to the mantle | Lodders 2003 condensation sequence gives water below 150 K more than 0.1 away | Lodders temperatures add condensation limits; no replacement of the share formula |
| BC2 | D | 1 The Brewer 2016 regressions and residual correlations of 4.7. 2 Hypatia for thick-disc and halo stars (Decided (Boss, 2026-10-17, between 11:19 AM and 11:23 AM EDT, card answers relayed in words by the coordinator): downloaded at build time, no derived table stored in the repository, tests offline, 6.5). | multivariate normal with the printed sds | a second survey (Hypatia or APOGEE) gives slopes outside the errors | a larger sample replaces the regression |
| BC3 | B | 1 The 2ND-HAND Stixrude values (usable under the standing rule). 2 The original paper when fetched. | none (equation of state) | densities at 100 GPa differ by more than 2 percent against the original | the original replaces the copy |
| BC4 | B | 1 Major oxides from the star by the 4.1 bookkeeping (D part); radioactive K, Th, U scaled with [Mg/H] and [Fe/H]. 2 Rocky-body compilations when fetched. | 0.2 dex on K, Th, U (P) | meteoritic ratios in a compilation sit outside the draw | compilations replace the scaling |
| BC5 | T | 1 Structure-model prior: core mass fraction Dirichlet bounded by the two measured planets as outer limits; inflation efficiency as in the table of tunable-gaps. | Dirichlet; inflation lognormal (P) | a statistical ensemble of radii and masses (or a tidal Love number set) lies outside | giant interior fits replace; Juno and Cassini are anchors only |
| BC6 | T | 1 Moons from the circumplanetary disc: composition ramp with distance in snow-line units. | ramp width 0.3 of the snow-line distance (P) | an exomoon or the regular-satellite census contradicts the ramp | any exomoon measurement replaces |
| BC7 | B | 1 NIST and PubChem solid densities where printed; others bracketed by the neighbours in the same molecular family. 2 Fray and Schmitt 2009 when fetched. | 3 percent (P) | a measured density falls outside the bracket | the originals replace |
| BC8 | B | 1 Plotnykov slope 0.94 of planet against star Fe/Mg with the printed band (body master). 2 Giant-impact stripping theory as a tail. | the printed band as width | larger samples change the slope outside its errors | a precise sample replaces; GEN.468-type changes ride the epoch |
| BC9 | B | 1 Bracket between the verified ice VII and ice X tables and the superionic range. 2 French 2009 and AQUA when extracted. | 5 percent on density (P) | an ab initio table differs by more than 5 percent | the table replaces |
| BC10 | B | 1 Lab reflectances as bracket for weathered regolith and frost, Bond conversion by a phase integral. 2 Hapke 2012 when read. | 0.05 in albedo (P) | an observed Bond albedo is outside | Hapke replaces the conversion |
| BC11 | P | 1 The loss fraction of volatiles at ejection is drawn uniform 0 to 1 and the draw is the fill; there is no "retain all volatiles" default. A labelled placeholder (Decided (Boss, 2026-10-17, between 11:19 AM and 11:23 AM EDT, card answers relayed in words by the coordinator): ships as a placeholder). | uniform 0 to 1 loss fraction (P) | any measurement of a rogue planet atmosphere | any such measurement |
| BC12 | B | 1 Each fit inside its stated window, Clausius-Clapeyron beyond, flagged and doubled in width. 2 Fray and Schmitt originals. | width doubles outside the window | a lab point outside the window disagrees | per-fit windows replace |
| BC13 | B | 1 Tissue factors from secondary sources (NASA and UNSCEAR reports on NTRS). 2 ICRP 103 when fetched. | none | the ICRP table differs | ICRP replaces |
| BC14 | T | 1 Belt composition by condensation zoning from BC1: rocky inside the snow line, water-rich beyond, drawing 4 of 84 names as already Decided. | zoning boundary at the star's snow line, width 0.3 AU in log (P) | belt spectra in other systems contradict | any measured belt replaces |
| BC15 | D | 1 CIPW-style norm computed from the bulk oxides and the IUGS total-alkali-silica classes. | none | a computed norm for a known rock composition fails to reproduce its mineral list | none (deterministic) |
| BC16 | D | 1 Groups as regions of Fe/Si, Mg/Si and oxidation state computed from BC1, BC2 and BC15. | none | grouping in the computed space does not reproduce the measured groups' separation | a larger meteorite sample refines region edges |
| BC17 | T | 1 Al-26 and Fe-60 from supernova and Wolf-Rayet yields with the odds of a nearby massive star from the embedded cluster mass function. | Bernoulli with cluster-size probability | enrichment ratios in a second system fall outside | measured ratios replace the odds |
| BC18 | T | 1 Core light elements from partitioning theory with the mantle redox of 4.8: more oxidised, less oxygen in metal. | Dirichlet over S, Si, O, H (P) | seismic or density deficits for a second body contradict | any measured second body |
| BC19 | B | 1 DA fraction against effective temperature from the large catalogues reachable through VizieR (DATA job; Decided (Boss, 2026-10-17, between 11:19 AM and 11:23 AM EDT, card answers relayed in words by the coordinator): downloaded at build time, tests offline, 6.5). | binomial about the fitted fraction | a larger catalogue changes the fraction outside its errors | the catalogue replaces the bracket |

### 5.3 Atmosphere

| Id | Kind | Fill (ranked) | Draw | Wrong if | Result lands as |
|---|---|---|---|---|---|
| AT.a1, a2 | D | 1 dIW = 2 log10(x_FeO / x_Fe) (4.8) with the 22-value white-dwarf ensemble (-1.25 ± 0.10) as zero point. 2 The existing 1.25 dex-per-dex slope as a cross-check. | iron oxidation fraction shared with BC1; activity correction 0.3 dex (P) | the white-dwarf ensemble disagrees with the identity by more than the correction | a larger ensemble replaces the zero point |
| AT.b4 | T | 1 Water-rich share from formation theory (the published 0 to 87 percent runs by model). | draw spans the range uniformly (P) | a population study measures the share | a population study |
| AT.b5 | P | 1 Log-uniform sub-percent water, labelled (Decided (Boss, 2026-10-17, between 11:19 AM and 11:23 AM EDT, card answers relayed in words by the coordinator): ships as a placeholder). | log-uniform (P) | any measurement | any ensemble |
| AT.n7 | P | 1 Haze yield log-uniform within the degenerate range, labelled (Decided (Boss, 2026-10-17, between 11:19 AM and 11:23 AM EDT, card answers relayed in words by the coordinator): ships as a placeholder). | log-uniform (P) | degeneracy resolved by spectra | spectra |
| AT.k7 | T | 1 Single-comet species ranked by trapping in amorphous ice: ordering by volatility. | Dirichlet in volatility order (P) | a second comet shows the opposite order | any second comet |
| AT.j2, j3, j4 | T | 1 Mars-, Titan-, Pluto-like occurrence from the escape and retention physics of the generator (shoreline argument), not a separate frequency. | none (an output) | an observed analogue population contradicts | any measured analogue population |
| AT.m3, m4, m11 | T | 1 Abiotic trace yields from lightning and photochemistry; presence follows the chemistry rows. | per-species lognormal (P) | a detection or strong limit on a temperate planet | any measurement |
| AT.F | B | 1 Each thin row: bracket from its neighbours as listed in the tunable-gaps note (16 rows). | per row as listed there | each row's falsifier | per-row data |
| AD1 | B | 1 Band strengths from neighbouring molecule bands; 2 Karkoschka tables extracted. | 10 percent (P); data downloaded at build time if it comes from a table, tests offline (6.5) | a measured giant colour disagrees | tables replace |
| AD2 | B | 1 Dust and haze spectra from similar materials. | 10 percent (P); data downloaded at build time if it comes from a table, tests offline (6.5) | measured spectra | tables replace |
| AD3 | D | 1 Compute correlated-k tables from the HITRAN and ExoMol line lists (reachable by 10:59 AM EDT). | none (computed); Decided (Boss, 2026-10-17, between 11:19 AM and 11:23 AM EDT, card answers relayed in words by the coordinator): downloaded at build time, no derived table stored in the repository, tests run offline (6.5) | the grid disagrees with a published spectrum | none |
| AD4 | B | 1 King factors from NIST for molecules with printed values; others by family. 2 Originals when fetched (Boss fetches them, section 8 item 5). | 5 percent (P); tests offline (6.5) | depolarization of SO2 or H2S measured outside | originals replace |
| AD5 | B | 1 Spectra from reachable NTRS reports; 2 originals. | 10 percent (P); data downloaded at build time, tests offline (6.5) | a shower run disagrees | originals replace |
| AD6 | D | 1 Download the hot-giant line lists (TiO, VO, FeH) and vapour pressures from ExoMol. | none; Decided (Boss, 2026-10-17, between 11:19 AM and 11:23 AM EDT, card answers relayed in words by the coordinator): downloaded at build time, no derived table stored in the repository, tests run offline (6.5) | tables disagree with lab data | none |
| AD7 | B | 1 H2SO4-H2O and brine vapour pressure from the nearest reachable fits; 2 Tabazadeh and Zeleznik originals (Boss fetches them, section 8 item 5). | 10 percent in log pressure (P); tests offline (6.5) | lab data differ | originals replace |
| AD8 | B | 1 Photochemistry rates from the NIST kinetics database where reachable; 2 JPL 19-5 (Boss fetches it, section 8 item 5). | factor 2 (P); data downloaded at build time, tests offline (6.5) | JPL rates differ by more than factor 2 | JPL replaces |

### 5.4 Count of the best fills

57 lines. By the kind of the first-ranked fill: **D 15** (GR3, GR4, GR9a, GR10, GR12, GR15, GX1, GX5, BC1, BC2, BC15, BC16, AT.a1/a2, AD3, AD6), **B 28**, **T 11**, **P 3** (BC11, AT.b5, AT.n7).

## 6. One convention for every filled value: the stand-in registry

Goal: any value that is not a verified measurement says so in the code, in the docs and in the loader report, and can be replaced by a result without a redesign. This extends what is Decided: the GEN.472 loader already carries source states (verified, partly verified, open, stand-in), Dirichlet where only a range is known, and a table version that bumps the generator epoch; GEN.503 already holds a target table (row id, statistic, value, sigma, source, label, validity range, selection-model id). The registry adds the fill itself and the link between them. Names below are proposals for the builder; nothing here exists in the repository.

### 6.1 The record

One record per fill, in a single module (proposed `planetgen/fills.py`), data not behaviour:

| Field | Meaning |
|---|---|
| `fill_id` | the id of this plan (GR4, BC1, AT.a1); permanent, never reused |
| `kind` | D, B, T or P (section 1) |
| `evidence` | READ, DATA, COMPUTED, 2ND-HAND, RECALLED, SNIPPET, LEAD; one label per input, worst one wins for the record |
| `source_state` | the GEN.472 state it maps to: D from READ or DATA inputs is "partly verified" (verified once the closure check passes); B is "partly verified"; T is "stand-in (theory)"; P is "stand-in (placeholder)"; anything resting on a LEAD is "open" |
| `form` | formula or table reference |
| `assumptions` | one line each, as written in this report |
| `validity` | the range inside which the fill holds, and what the code does outside it: draw with a doubled width and count it as an extrapolation |
| `knobs` | distribution, centre, width, with (V) or (P) |
| `target_id` | the GEN.503 target row it is tuned or checked against, if any |
| `falsifier` | the test id that fails when "wrong if" in section 5 is true |
| `replaced_by` | the research item that supersedes it (for example GEN.531 request id, or the reading job) |
| `sources` | citations; each marked READ only if opened by a thread, with the date and the thread |
| `licence` | for data tables, the licence or "public record" |
| `epoch` | the generator epoch in which this version took effect |

### 6.2 In code and in docs

Constants keep the style of `tuning.py`: a docstring with the source and the label. A filled value adds one line: `Fill GR4 (D, COMPUTED): ...` and reads its numbers from the registry rather than repeating them. In documents the tag is `[fill GR4 · D · COMPUTED]` right after the figure, and the figure is a range or a distribution, never a bare number. A sentence that rests on a LEAD says "Gemini, unverified" in the same sentence. The loader's start-up report prints counts by kind and by state, and the number of draws made outside a validity range; a build check fails if a fill with `source_state` "open" feeds a stored value (so unverified material can never reach a seed silently).

### 6.3 Replacing a stand-in

1. The result arrives (a paper read by a thread, a table downloaded, a Gemini claim).
2. The owner classifies it: verified (opened and read), partly verified (second-hand but consistent), or LEAD.
3. A LEAD may only change the registry's `knobs` width (widen a bracket) or add a comparison to the falsifier; the centre and the kind do not change.
4. A verified result replaces the fill: new record version, `replaced_by` filled, old record kept as superseded with its epoch, falsifier test updated to the new target.
5. If the replacement changes any stored value, it is a table version change and therefore a generator epoch bump, riding the pending epoch (Decided); if it only widens a draw's width in a way that changes stored values, same.
6. The docs tag changes from `[fill ...]` to a plain citation.

### 6.4 Conflicts

| Case | Rule |
|---|---|
| Verified against verified | neither wins; mixture or widened bracket; both stay in the record (GX2, GX3, GR11 scatter are examples) |
| Verified against Gemini or any LEAD | verified wins; the LEAD is logged as contradicted with the verified figure beside it |
| Partly verified against LEAD | partly verified wins unless the LEAD points to a primary source that is then read |
| LEAD with no contradiction | usable only to widen a bracket; stays "open" |
| Theory against data | data wins as soon as it is READ; the theory fill is kept as the label of what the draw was |
| Two theory fills | the simpler one; the other is recorded as an alternative in `assumptions` |
| Derived against measured | a derived fill that disagrees with a READ measurement fails its falsifier and is demoted to P until corrected |

### 6.5 Downloaded data and offline tests

Decided (Boss, 2026-10-17, between 11:19 AM and 11:23 AM EDT, relayed by the coordinator): line lists and star catalogues (HITRAN, ExoMol, Brewer 2016, Hypatia) are downloaded at build time. No derived table of them is stored in the repository. The cost is that a build needs the network, so the design must also say how tests run offline: tests never reach the network; they read either fixtures the test builds itself from synthetic values or a cache of earlier downloads kept outside the repository, and the record's `licence` field names the dataset licence. The data work of AD1 to AD8, BC2 and BC19 (group X9 and the BC2 and BC19 parts of X6 and X7) is gated on that design (the download step, the cache location, the fixtures and the offline test rule), not on reading the licences.

## 7. Gemini's first batch

### 7.1 Counts and recommendation

The paste of 10:52 AM EDT is about 34,000 characters. I count about 110 distinct factual claims (numbers, relations, attributions) and about 45 that can be checked against a paper reachable by us. I checked 20 against papers or computations already in hand: **8 agree, 7 are contradicted by a verified figure, 5 are unsupported, misattributed or internally inconsistent.** A full claim-by-claim check is its own job: I recommend one separate checking thread that works the other 25 checkable claims by opening the cited papers, and a second only if Boss wants the unverifiable theory content (the theory frontiers section, about 30 claims) screened for factual errors.

The paste is a copy of the document the thread received and sits unchanged at `/mnt/project-files/reports/gap-filling-plan-2026-10-17-files/gemini-batch-1-unverified.txt`.

### 7.2 The 20 claims

| # | Gemini, unverified | Check | Verdict |
|---|---|---|---|
| 1 | Halo star counts reach no lower than M* about 3.8e10 | the six measured haloes span 3.8e10 to 8e10 (Harmsen 2017, READ in GEN.406) | agrees |
| 2 | Ex situ fraction rises steeply with M*, about 0.1 at Milky Way mass | Rodriguez-Gomez 2016 abstract (READ): about 0.10 at that mass, rising | agrees |
| 3 | Scale height ratio thick/thin 2.5 to 3.8 | Yoachim and Dalcanton: 2.5 ± 30 percent (READ); ranges overlap | agrees |
| 4 | Scale length ratio 1.0 to 1.3 | Yoachim and Dalcanton formula gives 0.95 to 1.18 over 89 to 255 km/s | agrees |
| 5 | Vertical heating index 0.46 ± 0.03 | Sharma and Hayden 2020 Table 2: 0.441 ± 0.007 | agrees |
| 6 | Dispersions rise linearly with height; low- and high-alpha follow the same law; 3 to 20 kpc | Sharma and Hayden 2020 abstract, same wording as Gemini's sentence (Gemini cites "Sharma et al. 2021 GALAH DR3"; the paper read is Sharma and Hayden 2020) | agrees (attribution loose) |
| 7 | Mean cluster half-light radius 2.7 ± 0.3 pc, no dependence on cluster luminosity | Jordan 2005 abstract: 2.7 ± 0.35 pc, none on luminosity; but the paper finds dependence on host distance, surface brightness and colour, so "invariant across host luminosity" overstates | agrees (partly) |
| 8 | Blue clusters about 20 percent larger than red | Jordan 2005: red about 17 percent smaller, i.e. blue/red = 1.20 | agrees |
| 9 | Occupation half-point 10^8.6 Msun | Benitez-Llambay and Frenk 2020: 10^9.4 Msun/h for early reionization, only 4 times lower for z_re = 4, i.e. 10^8.8 at the extreme | contradicted |
| 10 | Occupation reaches 1 at 10^9.3 and 0 at 10^7.8 | the same paper: none below 3e8 (10^8.5), all above 5e9 (10^9.7) | contradicted |
| 11 | Thick/thin mass ratio log = 0.35 - 1.10 log(Vc/100) | gives 2.5, 1.2, 0.8 at 89, 181, 255 km/s, above the hard bound of 0.96, 0.64, 0.54 (4.4); also gives 0.94 at 220 km/s while the same passage says 0.2 to 0.45 above 200 | contradicted |
| 12 | sigma_z(10 Gyr) = 38.5 ± 1.8 km/s | Sharma and Hayden 2020 Table 2: 21.1 ± 0.2 | contradicted |
| 13 | sigma_R(10 Gyr) = 62.0 ± 2.5 km/s | same table: 39.4 ± 0.3 | contradicted |
| 14 | Radial heating index 0.42 ± 0.03 | same table: 0.251 ± 0.006 | contradicted |
| 15 | Nuclear cluster mass slope 0.65 ± 0.10, attributed to Neumayer 2020 | Neumayer 2020 gives 0.48; Georgiev gives 1.0 to 1.36 (both READ in GEN.406); 0.65 matches neither | contradicted |
| 16 | RMS of log(Vc/sigma0) = 0.12 dex, 64 galaxies | the Courteau 2007 text (read) prints no RMS; the passage also says both Vc = sqrt(2) sigma0 and Vc at least 2 sigma0 | unsupported, inconsistent |
| 17 | Blue/red cluster system radius ratio 1.20 ± 0.08 from Forbes 2012, Rhode and Zepf 2004, Kartha 2014 | the Kartha 2014 text searched, ratio not found; GEN.406 found papers that disagree (equal within errors against more extended) | unsupported |
| 18 | Low-mass stellar-halo relation slope 1.90, zero point 6.80, from UniverseMachine | Behroozi 2019 Appendix J fits only 10^10.5 to 10^15 (READ), so the paper has no low-mass parameterisation to cite | misattributed |
| 19 | Thick-disc start 12 to 13 Gyr lookback at z = 2.5 to 4.0 (Pinna 2019) | with Planck-like cosmology 12 Gyr is z = 3.5 and 13 Gyr is z = 6.8 (computed), and 2.5 to 4 is 11.2 to 12.3 Gyr; Pinna 2019 not located by search | inconsistent, unsupported |
| 20 | Halo formation redshift z50 = 1.20 - 0.35 log(M/1e12), from Aquarius, Cooper 2010, Deason 2016, D'Souza and Bell 2018 | Cooper 2010 haloes are Milky Way mass only; no relation of this form located in the sources we hold | unsupported |

Not counted: the occupation width 0.25 ± 0.04 dex (agrees with the 0.20 to 0.25 derived from the verified edges, 4.5) is in the same row as the contradicted centre and is kept in the bracket for GR15 as a LEAD that agrees.

### 7.2a Result of the checking thread (2026-10-17)

The checking thread finished at 11:49 AM EDT as it reported ([gemini-batch1-check.md](gemini-batch1-check.md)). It checked 51 claims (the 20 above plus 31 more): 10 agree, 17 partly, 16 contradicted by a read source, 4 unsupported, 4 misattributed. Five of the 20 above moved after more reading. The corrections that matter for the build items: the occupation centre is 10^9.33 Msun/h with width 0.21 dex (Gemini's 10^8.6 is the paper's no-reionization limit); the 1.9 low-mass stellar-to-halo slope is Garrison-Kimmel 2014 (1.92), not UniverseMachine; the halo formation redshift z50 is 1.02 at 1e10 and 0.74 at 1e12 (Gemini 1.9 and 1.2); the age-velocity normalisations are 21.1 and 39.4 km/s with indices 0.441 and 0.251 (Gemini 38.5, 62.0, 0.46, 0.42); the Courteau scatter is 0.097 dex about the paper's own line (79 markers digitised from its Fig 2; the 64 galaxies were a misreading). Requests now answered by a usable table: R2 (all three curves), R5, R12 (the z_re table with Busha's spread), R13 (Munoz ultra-faint sizes) and parts of R3, R6 and R10. Still open for Boss's list: R1, R4 beyond the Milky Way, R7, R8, R14 and parts of R6, R9, R10, R11. Three conflicts between two verified sources (thick-to-thin height ratio 2.5 Yoachim and Dalcanton against 4.2 Comeron 2018 median; thick-to-thin mass ratio against Vc; cluster system radius about 3.5 against about 6 times the galaxy's) are open and wait on Boss: his decision 2 (widen and keep both, 11:19 AM EDT) covered arm pitch against bulge fraction, the nuclear cluster mass slope and the blue to red cluster radius ratio only, and the coordinator posted a card for each of the three at about 11:56 AM EDT.

### 7.3 How results enter

Every Gemini figure is stored as a LEAD record next to the fill it concerns, with the claim number from 7.2, the verified figure it was checked against and a verdict. Agreeing claims raise the fill's evidence only if a verified source is also cited; contradicted claims are kept and show in the loader report as "LEAD rejected". The checking thread of 7.1 turns the remaining checkable claims into verdicts of the same shape; the second batch Boss expects uses the same intake. The theory frontier content (WDM, SIDM, fuzzy dark matter, thick-disc origin paradigms, assembly bias) is not a measurement and is never a default; it can only supply the form of a T-kind fill, for example the gas-rich turbulent collapse picture behind GR6, and it is labelled theory wherever it appears.

## 8. What needed Boss or a data set: Decided

Seven items needed a person or a data set; nothing else in the 57 lines waits on a person. Boss answered all seven between 11:19 AM and 11:23 AM EDT on 2026-10-17 (card answers relayed in words by the coordinator). Each row is Decided.

| # | Item | Decided | Gaps |
|---|---|---|---|
| 1 | Does the generator have an environment (field, group, cluster)? | Draw the environment from the halo mass; no new environment concept now. The 0.3 dex cluster-count factor applies with the group probability. | GX4, GR11 (scatter choice), AT.j |
| 2 | Policy for verified-against-verified disagreements | Widen and keep both; neither source is dropped. GX2: sign 50/50. GX3: uniform 0.48 to 1.36. GR7: log-uniform from 1.0 to the largest printed ratio, 1.5 until the ratios are read. | GX2, GX3, GR7 |
| 3 | Three labelled placeholders | Ship as labelled placeholders: BC11 a uniform 0 to 1 loss fraction (the draw is the fill, there is no "retain all volatiles" default), AT.b5 and AT.n7 log-uniform. | BC11, AT.b5, AT.n7 |
| 4 | Licence for downloaded line lists and star catalogues (HITRAN, ExoMol, Brewer 2016, Hypatia) | Download at build time; do not store derived tables. Builds need the network, so the design states how tests run offline (6.5); group X9 is gated on that design, not on the licence read. | AD3, AD6, BC2, BC19 (and AD1 to AD8 where they use a table) |
| 5 | Papers behind refused publishers | Boss fetches them himself. The labelled stand-ins stay in the items meanwhile; the body reading batch X7 and the items that need the originals wait for his files. The papers are listed below and on GEN.531. | BC3, BC4, BC7, BC12, BC13, AD4, AD7, AD8 |
| 6 | Courteau 2007 table and Comeron fit equations | Ship the plan's fills and also run the quick VizieR search for Courteau 2007 and the Comeron 2011, 2012 and 2014 tables first. | GR5, GR3 (rank 2) |
| 7 | Thick-disc start in other galaxies | The halo formation epoch plus a 0 to 1 Gyr delay, a labelled placeholder, replaced when S0 histories are found. | GR6 |

### 8.1 The papers Boss fetches (item 5)

The full list is [papers-to-fetch.md](papers-to-fetch.md): the nine papers Boss asked for first (section A there), the other papers the files could not read (B) and what is not needed (D), each with a Crossref citation and the places to try a free copy. Corrections to the first version of this table, from that list:

- BC3: Stixrude and Lithgow-Bertelloni has two parts, Part I (GJI 162, 610, 2005, doi 10.1111/j.1365-246x.2005.02642.x, bronze open access) and Part II (GJI 184, 1180, 2011, doi 10.1111/j.1365-246x.2010.04890.x); the 2024 paper III named earlier is what the second-hand copy carries.
- BC4: Palme and O'Neill and Rudnick and Gao each have a 2014 second edition of the Treatise on Geochemistry (prefer it) and an earlier first edition.
- BC7, BC12: Fray and Schmitt 2009 is unchanged; a HAL record exists. Bierson and Nimmo 2019 has a free author copy in PubMed Central (PMC7058130), which a thread can fetch itself.
- BC13: the free extract of ICRP 103 stops before the wR and wT tables (section 4, pages 61 to 78), so those pages are the ones to fetch.
- AD4: the numbers needed (depolarization of SO2 and H2S) are in Bogaard, Buckingham, Pierens and White 1978 (J. Chem. Soc., Faraday Trans. 1, 74, 3008, doi 10.1039/f19787403008), not in King 1923; King (Proc. R. Soc. A 104, 333, doi 10.1098/rspa.1923.0113) is the origin of the factor, so the item names both. Sneep and Ubachs 2005, Bates 1984 and Bucholtz 1995 (Appl. Opt. 34, 2765) stay on the list.
- AD7: Tabazadeh, Toon, Clegg and Hamill 1997 (GRL 24, 1931, bronze open access at Wiley) and Zeleznik 1991 (JPCRD 20, 1157).
- AD8: JPL Publication 19-5 is Burkholder et al. 2020, Evaluation No. 19.
- Free copies that make some "refused" rows readable by a thread: Lodders 2003 (IOP bronze PDF, for BC1), Sossi 2020 (open access, already read).
- On arXiv, so a research thread reads them and Boss need not: Aumer and Binney 2009 (0905.2512, GR4), Haywood 2013 (1305.4663) and Bensby 2014 (1309.2631) (GR6), Pinna 2019a and b (1901.04310, 1904.01260) (GR6), Graus 2019 (1808.03654, an inference) and Graus's FIRE age-gradient paper (1901.05487) (GR15). GR4, GR6 and GR15 therefore move from class B to class A in section 3.

## 9. Order of work

Grouped for builders; groups with no arrow between them run in parallel. "Gate" says what the group holds up.

| Group | Contents | Needs | Gate |
|---|---|---|---|
| X1 Registry | `fills.py` record, loader report, build check on "open" fills, docs tag; hooks into GEN.472 source states and the GEN.503 target table | nothing | first, because every other group registers into it |
| X2 Galaxy derived fills | GR3 bound, GR4 rule, GR12 and GR15, GR9a, GX1 | X1 | 8.1 phase 5 Galaxy laws build (GEN.401, 369, 403); GR4 feeds GEN.447, which GEN.301 names |
| X3 Accreted halo | integral of 4.2: GR1, GR10, GX5, GR14 with the slope test against 0.79 | X1, X2 (GX1 for halo mass given M*) | GEN.401, 403, 405, 386 |
| X4 Galaxy reading jobs | GR2 digitising, GR5, GR7, GR8, GR9b, GR9c, GR11, GR13 | none; research threads | GEN.402, 403, 369, and GR9c feeds GEN.280, which GEN.301 names |
| X5 Galaxy priors | GR6, GX2, GX3, GX4 | X3 (GR6 uses the GR10 epoch); Boss items 1, 2 and 7 Decided | GEN.386, 402, 403 |
| X6 Star and oxidation package | BC2 regression, BC1 water share, shared iron-oxidation variable, AT.a1 and a2, BC4 major oxides, BC15, BC16 | X1 | 9.0; building blocks already back-scattered to 8.1 (GEN.472 phase 2) |
| X7 Body reading jobs | BC3, BC4 originals, BC7, BC9, BC10, BC12, BC13, BC19 | the originals wait for the files Boss fetches (item 5, 8.1); BC19 downloads at build time (6.5) | 9.0 |
| X8 Theory priors | BC5, BC6, BC8, BC14, BC17, BC18, AT.b4, AT.k7, AT.j, AT.m | X6 for the shared variables | 9.0 |
| X9 Atmosphere data | AD1 to AD8 (downloads and extractions) | the download and offline-test design (6.5), not the licence read; AD4, AD7 and AD8 originals wait for Boss's files | 9.0 |
| X10 Placeholders | BC11, AT.b5, AT.n7 | none (Boss item 3 Decided) | 9.0 |
| X11 Thin F rows | the 16 AT.F rows | the tunable-gaps builders | 9.0 |

Which gaps block what, as far as I traced:

- **Block the 8.1 headline chain**: the two gaps named in GEN.301's own prerequisite line, GR4 (through GEN.447) and GR9c (through GEN.280). GR4 is derived and ready (X2). GR9c is a bracket with a virial-theorem part and one reading job (X4); it can ship as a labelled bracket without waiting for the review. I did not trace the other prerequisites of GEN.301 to gaps, so there may be more.
- **Block the 8.1 stamp** (the six-phase "Galaxy laws build" items listed in the lane report): all galaxy rows of 5.1, in the order X2, X3, X4, X5. Boss decided items 1, 2, 6 and 7 of section 8, and each has a fill that ships without further answers (item 6 also runs the quick VizieR search first).
- **Do not block 8.1**: every body and atmosphere row. They belong to 9.0, with the generation changes that alter stored values staying with 8.4 as already Decided; X6 is the exception worth pulling early because its outputs (the star regression and the water-share function) are inputs to GEN.472 phase 2 and the 8.1 back-scattered building blocks.
- **Epoch**: a fill that replaces a value later changes stored output and rides the pending epoch; the X1 registry records the epoch per fill so the builder of GEN.322-324 style rescatters can see which stand-ins sat inside it.

## 10. Limits of this report

- Nothing in the repository was changed. No papers were taken from refused publishers. Every number marked COMPUTED can be reproduced with the scripts in the files folder.
- READ means a thread opened the paper or table; abstracts only are marked "abstract READ". The Jordan 2005 and Rodriguez-Gomez 2016 statements rest on abstracts. Nothing in 7.2 is more than one check deep.
- The accreted-halo integral is a sketch: it applies the Moster z = 0 relation at the accretion mass, ignores satellite stripping, and over-predicts by a factor of 2 to 3 against Illustris. It is a shape source and a test, not a law.
- The oxygen-budget water share has no disk chemistry. The age-velocity rule has one closure target (the Milky Way) and so cannot yet be tested on a second galaxy.
- The Brewer 2016 regressions use the whole catalogue without its quality cut; the cut must be applied by the builder (4.7).
- Gemini claim counts (110 distinct, 45 checkable) are my estimates from reading the paste; the 20 checked are exactly those in 7.2.

Files: `/mnt/project-files/reports/gap-filling-plan-2026-10-17-files/` holds `water_share.py`, `acc.py`, `ws_stars.py`, `arx.py`, `brewer2016_vizier_J_ApJS_225_32.tsv` and `gemini-batch-1-unverified.txt`.
