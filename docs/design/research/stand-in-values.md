# Every stand-in value, with its source tier and its tweak level (GEN.725)

Status: generated from the stand-in registry of [gap-filling-plan.md](gap-filling-plan.md) by `scripts/build_stand_in_values.py`; do not edit it by hand (`python scripts/build_stand_in_values.py --write` rewrites it, and a test fails when it differs from the registry). The tier and the tweak level of each row are the registry's; a row's owner confirms them.

A stand-in is any value that is not a verified measurement of the thing it stands for. Every one carries two labels (Boss, 3:40 PM EDT on 2026-10-24; the tweak scale decided 4:01 PM EDT): **where its number came from** (the source tier) and **how far it was bent to give the outcome wanted** (the tweak level). Any change that adds a stand-in adds its row to the registry with both labels in the same pull request.

There are 142 rows below. A row may be a family of values (a recipe table, a set of bands); its tier cell then names the tier of each part, and its tweak is the highest level any part has.

## The source tiers

| Tier | Name | Meaning |
|---|---|---|
| T1 | Verified | read in the source and the value is the source's |
| T2 | Confident unverified | from a source we trust but have not opened or cannot check |
| T3 | Marginal unverified | from a source we doubt, or a bracket from sources that disagree |
| T4 | Unverified | no source at all for the value, only a reason for its range |
| T5 | Confident educated guess based on physics | a law of physics fixes the form and the size to within a factor of a few |
| T6 | Marginal educated guess based on physics | physics fixes the form, but the size is argued, not derived |
| T7 | Educated guess based on inferred physics | the physics is itself inferred from other bodies or other regimes |
| T8 | Educated guess based on reputable fringe or new physics | a published but not yet accepted mechanism |
| T9 | Educated guess based on needed output (full tweak) | chosen so that the generated sky matches what has to come out |

Unverified and guessed tiers may be tweaked to hit target probabilities (Boss, 3:40 PM EDT on 2026-10-24).

## The tweak levels

| Level | Name | Meaning |
|---|---|---|
| L0 | As envisioned | the value is used as its source gives it |
| L1 | Nudged inside the source's uncertainty | moved within the error bar the source prints |
| L2 | Tuned to a published statistic | moved until a published count or distribution is reproduced (a tuning target) |
| L3 | Hand-tuned against a sweep | chosen by running the generator over many draws and looking at the outcome |
| L4 | Fitted to the outcome with no independent source | the number is whatever makes the output right |

## How many rows are at each tier and tweak level

| Tier | L0 | L1 | L2 | L3 | L4 | All |
|---|---|---|---|---|---|---|
| T1 | 21 | 0 | 0 | 0 | 0 | 21 |
| T2 | 38 | 0 | 0 | 0 | 0 | 38 |
| T3 | 30 | 1 | 1 | 1 | 0 | 33 |
| T4 | 5 | 0 | 0 | 0 | 0 | 5 |
| T5 | 6 | 2 | 0 | 0 | 0 | 8 |
| T6 | 14 | 4 | 1 | 4 | 0 | 23 |
| T7 | 10 | 2 | 0 | 1 | 0 | 13 |
| T8 | 0 | 0 | 0 | 0 | 1 | 1 |
| T9 | 0 | 0 | 0 | 0 | 0 | 0 |
| All | 124 | 9 | 2 | 6 | 1 | 142 |

## How the labels were given

- **Tier** follows the kind of the fill and the evidence label of its inputs: a derived (D) fill from a source we read is T1 or T2; a bracket (B) is T2 or T3 by how many sources agree and whether we opened them; a theory (T) fill is T5 to T7 by how much of the physics is derived and how much inferred; a placeholder (P) with no source is T4 or, when physics argues its size, T6 or T7. Rows whose parts differ list each part.
- **Tweak** is L0 unless the value was moved after the fact. The L3 rows (hand-tuned against a sweep) are the mantle redox laws AT.a1 and a2, the atmosphere draws AT.b4 and b5, BC6 and AM19 (the exobase ratios); the L2 row is AM10 (the alkali activity, tuned to Kite et al. 2016); L1 rows are nudged inside the source's error bar. No row is L4: a row at L4 needs the item that justified it.
- The recipe rows of the atmosphere draw (`atmosphere_recipes_data.py`) are not rows here one by one. Those that carry a tuning target in `tuning_targets.json` (a3, a4, e1, f5, f7, k1 to k5, k7, o2) are L2: their weights were tuned to a published statistic.

## 5.1 Galaxy laws

| Id | Kind | Tier | Tweak | Fill (ranked) | Wrong if | Result lands as |
|---|---|---|---|---|---|---|
| GR1 | B | T3 | L0 | Placeholder median of GEN.406 anchored on the six measured haloes (Harmsen 2017). | a measured halo falls outside 2 sigma; the integral slope cannot get within 0.3 of 0.79 | a simulation or survey fraction replaces rank 1; Gemini's "ex situ below 0.10 under 1e10.5" (LEAD) is only an upper cap |
| GR2 | B | T2 | L0 | The anchor values already read from the three figures (GEN.406). | digitised curve differs from rank 1 by more than 0.1 dex | digitised table replaces; Gemini's 8-row table (LEAD) is compared with it, not used |
| GR3 | D | T1 | L0 | The hot-component bound of 4.4 with the drawn classical bulge fraction. | 39 percent of the real galaxies exceed the bound galaxy by galaxy (section 4.9), so the test applies to the mean relation or carries a margin of 0.13 dex, and does not reject each seed | a published fit equation with scatter replaces the scatter only; Gemini's formula (LEAD) is rejected by the bound |
| GR4 | D | T2 | L0 | Martinsson 2013 age-averaged sigma_z times (1+beta)(t/tau)^beta, beta from Sharma and Hayden 2020 (4.3). | the Milky Way closure misses by more than 35 percent | a second galaxy with ages tests the rule; Gemini's 38.5 and 62.0 km/s (LEAD) are contradicted and not used |
| GR5 | B | T3 | L0 | Bracket from Martinsson's 15 percent and Gemini's 0.12 dex (LEAD may only widen), so 0.06 to 0.12 dex centred on 0.09. | the table gives an RMS outside 0.06 to 0.12 | the RMS from the table or the authors replaces the centre |
| GR6 | T | T6 | L0 | Thick-disc start set by the gas-rich turbulent collapse picture: the host's halo formation epoch (GR10 integral) plus a delay. | S0 spectroscopy shows a start independent of mass while GR10 gives a mass dependence | external-galaxy histories replace the delay; Gemini's 12 to 13 Gyr start (LEAD, section 7 claim 19) cannot set the epoch because its redshifts do not match its lookback times |
| GR7 | B | T3 | L0 | Lower end 1.0 (the papers that find equal extent within errors, READ in GEN.406); upper end the largest printed ratio among Kartha 2014, Forbes and Ponman 2012, Hargis and Rhode 2014 (to read). | a fourth survey puts the ratio outside the bracket | printed ratios replace the ends; Gemini's 1.20 ± 0.08 (LEAD, attribution unverified) may centre the draw after one of the three papers shows it |
| GR8 | B | T3 | L0 | Interpolate the radius of the cluster system in log against log M* between the verified low-mass and massive ends. | measured radii sit outside the interpolation by more than 0.3 dex | measured radii replace the interpolation; Gemini's ranges (2 to 8, 10 to 25, over 100 kpc) are LEAD |
| GR9a | D | T1 | L0 | Jordan 2005: mean 2.7 ± 0.35 pc, red 17 percent smaller than blue, small dependences on distance, surface brightness and host colour, none on cluster luminosity (abstract READ). | the 157 Harris clusters disagree with the median by more than the paper's scatter | coefficients of the three dependences replace the multipliers |
| GR9b | B | T3 | L0 | Old track ages from the infall times of the accreted satellites (GR10 integral), young track from the disc epoch (GR6), weights from the blue/red split. | simulated ages differ from the rule by more than 2 Gyr | simulation table replaces the rule |
| GR9c | B | T6 | L0 | Dispersion from the virial theorem with nuclear cluster mass and half-light radius (D); ages as an old base from the host epoch plus recurrent young populations for M* above 1e9. | the review shows a young fraction outside 0 to 0.3 or dispersions far from the virial value | review values replace the knobs; Gemini's 15 to 50 km/s and young below 100 Myr (LEAD) are consistent and widen nothing |
| GR10 | D | T2 | L0 | Infall times from the accreted-halo integral (4.2) with the host's own assembly history. | Cooper and Cole's ages differ by more than 1.5 Gyr at the same halo mass | their ages replace the integral; Gemini's z50 formula (LEAD, claim 20) is compared |
| GR11 | B | T2 | L0 | Moster 2013 extended below 1e10.5 and flagged as an extrapolation (Behroozi 2019 fit range READ: only 10^10.5 to 10^15). | the median at 1e9.5 differs from the dwarf-galaxy data by more than the scatter | a dwarf fit replaces rank 1; Gemini's slope 1.90 and zero point (LEAD, claim 18) are not used as a median |
| GR12 | D | T2 | L0 | Occupation centre shift 0.1 dex per unit z_re from Benitez-Llambay and Frenk and z_re = 7.68 ± 0.79 (Planck 2018) (4.5). | the patchy-reionization maps give a spread more than three times 0.79 | Busha and Alvarez or Katz widths replace 0.79; Gemini's void and overdense ranges (LEAD) widen only after the maps are read |
| GR13 | B | T3 | L0 | Continue the verified dwarf size law downward to M_V near -2, floor at the cluster size 2.7 pc. | measured r_h sit outside the continuation by more than 0.3 dex | the paper's relation replaces; Gemini's mean 80 pc, 20 to 300 pc, 0.35 dex (LEAD) is a first comparison |
| GR14 | B | T3 | L0 | Logistic rise of the ex situ fraction in log M* from the verified 0.1 at 5e10 to a ceiling of 1 at the top, centre fixed by the digitised Illustris curve (GR2). | the digitised curve crosses 0.5 outside the bracket | digitised curve replaces; Gemini's crossover 1.5e11 ± 0.2e11 (LEAD) tests the result |
| GR15 | D | T2 | L0 | Width from the verified edges: 3e8 to 5e9 is 1.22 dex, so 0.20 to 0.25 dex (4.5). | a fit to the Benitez-Llambay and Frenk simulations gives a width outside 0.15 to 0.35 | the fitted width replaces; Gemini's 0.25 ± 0.04 agrees (claim not counted: its centre is contradicted) |
| GX1 | D | T1 | L0 | Bayes: the probability of a halo mass given the stellar mass is proportional to the probability of the stellar mass given the halo mass times the halo mass function, with the scatter of the Moster relation. | the sampled halo masses at fixed M* have a median shifted from the inverse of the median relation by something other than the Eddington-bias term | a published conditional replaces |
| GX2 | B | T1 | L0 | Mixture of the two published signs (Hart 2017, Yu and Ho 2019) with equal weights; neither wins because both are verified. | a third survey fixes the sign | the surviving sign replaces the mixture |
| GX3 | B | T2 | L0 | Slope uniform on 0.48 to 1.36 (Neumayer against Georgiev). Gemini's 0.65 ± 0.10 sits inside and changes nothing. | the review shows one fit to be a selection artefact | a reconciling paper replaces |
| GX4 | T | T5 | L0 | Environment drawn from the halo mass (more massive halos are in groups and clusters); the 0.3 dex factor applies with that probability. | the generator gains an environment concept that disagrees | a later environment concept replaces the draw |
| GX5 | D | T2 | L0 | [Fe/H] of the halo = relation(M* of the dominant progenitor) minus 0.15 (D'Souza and Bell 2018, READ), progenitors from the 4.2 integral. | the six measured haloes lie outside 2 sigma | a metallicity-mass paper for dwarfs replaces the extension |

## 5.2 Body composition

| Id | Kind | Tier | Tweak | Fill (ranked) | Wrong if | Result lands as |
|---|---|---|---|---|---|---|
| BC1 | D | T2 | L0 | The oxygen budget of 4.1 per star. | Lodders 2003 condensation sequence gives water below 150 K more than 0.1 away | Lodders temperatures add condensation limits; no replacement of the share formula |
| BC2 | D | T1 | L0 | The Brewer 2016 regressions and residual correlations of 4.7. | a second survey (Hypatia or APOGEE) gives slopes outside the errors | a larger sample replaces the regression |
| BC3 | B | T2 | L0 | The 2ND-HAND Stixrude values (usable under the standing rule). | densities at 100 GPa differ by more than 2 percent against the original | the original replaces the copy |
| BC4 | B | T3 | L0 | Major oxides from the star by the 4.1 bookkeeping (D part); radioactive K, Th, U scaled with [Mg/H] and [Fe/H]. | meteoritic ratios in a compilation sit outside the draw | compilations replace the scaling |
| BC5 | T | T7 (core), T2 (inflation peak) | L1 | Structure-model prior: core mass fraction Dirichlet bounded by the two measured planets as outer limits; inflation efficiency as in the table of tunable-gaps. | a statistical ensemble of radii and masses (or a tidal Love number set) lies outside | giant interior fits replace; Juno and Cassini are anchors only |
| BC6 | T | T6 | L3 | Moons from the circumplanetary disc: composition ramp with distance in snow-line units. | an exomoon or the regular-satellite census contradicts the ramp | any exomoon measurement replaces |
| BC7 | B | T2 | L0 | NIST and PubChem solid densities where printed; others bracketed by the neighbours in the same molecular family. | a measured density falls outside the bracket | the originals replace |
| BC8 | B | T6 (tail), T1 (slope, one paper) | L1 | Plotnykov slope 0.94 of planet against star Fe/Mg with the printed band (body master). | larger samples change the slope outside its errors | a precise sample replaces; GEN.468-type changes ride the epoch |
| BC9 | B | T3 | L0 | Bracket between the verified ice VII and ice X tables and the superionic range. | an ab initio table differs by more than 5 percent | the table replaces |
| BC10 | B | T3 | L0 | Lab reflectances as bracket for weathered regolith and frost, Bond conversion by a phase integral. | an observed Bond albedo is outside | Hapke replaces the conversion |
| BC11 | P | T6 | L0 | The loss fraction of volatiles at ejection is drawn uniform 0 to 1 and the draw is the fill; there is no "retain all volatiles" default. A labelled placeholder (Decided (Boss, 2026-10-17, between 11:19 AM and 11:23 AM EDT, card answers relayed in words by the coordinator): ships as a placeholder). | any measurement of a rogue planet atmosphere | any such measurement |
| BC12 | B | T2 | L0 | Each fit inside its stated window, Clausius-Clapeyron beyond, flagged and doubled in width. | a lab point outside the window disagrees | per-fit windows replace |
| BC13 | B | T2 | L0 | Tissue factors from secondary sources (NASA and UNSCEAR reports on NTRS). | the ICRP table differs | ICRP replaces |
| BC14 | T | T6 | L0 | Belt composition by condensation zoning from BC1: rocky inside the snow line, water-rich beyond, drawing 4 of 84 names as already Decided. | belt spectra in other systems contradict | any measured belt replaces |
| BC15 | D | T1 | L0 | CIPW-style norm computed from the bulk oxides and the IUGS total-alkali-silica classes. | a computed norm for a known rock composition fails to reproduce its mineral list | none (deterministic) |
| BC16 | D | T2 | L0 | Groups as regions of Fe/Si, Mg/Si and oxidation state computed from BC1, BC2 and BC15. | grouping in the computed space does not reproduce the measured groups' separation | a larger meteorite sample refines region edges |
| BC17 | T | T6 | L0 | Al-26 and Fe-60 from supernova and Wolf-Rayet yields with the odds of a nearby massive star from the embedded cluster mass function. | enrichment ratios in a second system fall outside | measured ratios replace the odds |
| BC18 | T | T7 | L1 | Core light elements from partitioning theory with the mantle redox of 4.8: more oxidised, less oxygen in metal. | seismic or density deficits for a second body contradict | any measured second body |
| BC19 | B | T2 | L0 | DA fraction against effective temperature from the large catalogues reachable through VizieR (DATA job; Decided (Boss, 2026-10-17, between 11:19 AM and 11:23 AM EDT, card answers relayed in words by the coordinator): downloaded at build time, tests offline, 6.5). | a larger catalogue changes the fraction outside its errors | the catalogue replaces the bracket |

## 5.3 Atmosphere

| Id | Kind | Tier | Tweak | Fill (ranked) | Wrong if | Result lands as |
|---|---|---|---|---|---|---|
| AT.a1, a2 | D | T7 | L3 | dIW = 2 log10(x_FeO / x_Fe) (4.8) with the 22-value white-dwarf ensemble (-1.25 ± 0.10) as zero point. | the white-dwarf ensemble disagrees with the identity by more than the correction | a larger ensemble replaces the zero point |
| AT.b4 | T | T6 | L3 | Water-rich share from formation theory (the published 0 to 87 percent runs by model). | a population study measures the share | a population study |
| AT.b5 | P | T6 | L3 | Log-uniform sub-percent water, labelled (Decided (Boss, 2026-10-17, between 11:19 AM and 11:23 AM EDT, card answers relayed in words by the coordinator): ships as a placeholder). | any measurement | any ensemble |
| AT.n7 | P | T6 | L1 | Haze yield log-uniform within the degenerate range, labelled (Decided (Boss, 2026-10-17, between 11:19 AM and 11:23 AM EDT, card answers relayed in words by the coordinator): ships as a placeholder). | degeneracy resolved by spectra | spectra |
| AT.k7 | T | T1 (O2, 67P noble gases), T2 (N2 over CO) | L0 | Single-comet species ranked by trapping in amorphous ice: ordering by volatility. | a second comet shows the opposite order | any second comet |
| AT.j2, j3, j4 | T | T6 | L0 | Mars-, Titan-, Pluto-like occurrence from the escape and retention physics of the generator (shoreline argument), not a separate frequency. | an observed analogue population contradicts | any measured analogue population |
| AT.m3, m4, m11 | T | T5 (N2O, NO2), T7 (PH3, CH3Cl) | L0 | Abiotic trace yields from lightning and photochemistry; presence follows the chemistry rows. | a detection or strong limit on a temperate planet | any measurement |
| AT.F | B | T3 | L0 | Each thin row: bracket from its neighbours as listed in the tunable-gaps note (16 rows). | each row's falsifier | per-row data |
| AD1 | B | T2 | L0 | Band strengths from neighbouring molecule bands; | a measured giant colour disagrees | tables replace |
| AD2 | B | T3 (dust and haze spectra), T1 (Wolff 2009 solar-band single-scattering albedo) | L0 | Dust and haze spectra from similar materials. | measured spectra | tables replace |
| AD3 | D | T1 | L0 | Compute correlated-k tables from the HITRAN and ExoMol line lists (reachable by 10:59 AM EDT). | the grid disagrees with a published spectrum | none |
| AD4 | B | T2 (air gases), T3 (SO2 and H2S King factors via PICASO) | L0 | King factors from NIST for molecules with printed values; others by family. | depolarization of SO2 or H2S measured outside | originals replace |
| AD5 | B | T3 | L0 | Spectra from reachable NTRS reports; 2 originals. | a shower run disagrees | originals replace |
| AD6 | D | T1 | L0 | Download the hot-giant line lists (TiO, VO, FeH) and vapour pressures from ExoMol. | tables disagree with lab data | none |
| AD7 | B | T1 (H2SO4-H2O, Zeleznik 1991), T2 (acid vapour pressure via Gao 2014), T4 (brines) | L0 | H2SO4-H2O and brine vapour pressure from the nearest reachable fits; | lab data differ | originals replace |
| AD8 | B | T2 | L0 | Photochemistry rates from the NIST kinetics database where reachable; | JPL rates differ by more than factor 2 | JPL replaces |

## 6.3a Ranges of the galaxy shape dials

| Lever | Range | Label | Tier | Tweak | Source |
|---|---|---|---|---|---|
| `--shape` | spiral | theory | T5 | L0 | the galaxy families; only the spiral exists, later families fill the list |
| `--mass` | 281 to 1.55e+11 solar masses | derived | T2 | L0 | Boss, 10-21 (7:04 to 7:31 PM EDT): floor half the stellar mass of Segue 1, 562 solar masses (Local Volume Database v1.1.1; Simon et al. 2011, arXiv:1007.4198); ceiling 1.5 x Andromeda's 1.03e11 (Sick et al. 2014, arXiv:1410.0017), read by the galaxy size research thread |
| `--compactness` | -3 to 3 sigma | derived | T2 | L0 | scatter of the late-type size-mass law, Shen et al. 2003 (read in the repository); three sigma is the draw's own cap |
| `--density` | 0.04 to 100000 solar masses per pc^2 | derived | T2 | L0 | density-limits report (10-21): the smallest density with a planned sector at a 4 pc scale height is 0.043; the largest found in dense stellar systems is 1e5 (Hopkins et al. 2010, arXiv:0908.4088); the floor also grows with the scale height (check_outline) |
| `--disk-scale-length` | 28 to 8250 pc | derived | T2 | L0 | Boss, 10-21 (7:04 PM to 8:08 PM EDT): floor 28 pc (Segue 1's half, 7.1 pc, is under the 4 pc cell grid's reach: a variable sector edge, a later item, would bring it back); ceiling 1.5 x Andromeda's 5.5 kpc (Courteau et al. 2011, arXiv:1106.3564; Dorman et al. 2013, arXiv:1310.4179), read by the galaxy size research thread |
| `--flatness` | 3 to 20 | bracket | T3 | L0 | median 7.3 and 0.13 dex read in Kregel, van der Kruit and de Grijs 2002 (shape_draw); the range 3 to 20 is the proposal's suggestion (galaxy size levers, 10-21), not read in a source |
| `--bulge-share` | 0 to 0.9 | bracket | T3 | L0 | stage medians from Graham and Worley 2008 Table 4 (shape_draw); the range 0 to 0.9 is the proposal's suggestion (galaxy size levers, 10-21), not read in a source |
| `--bulge-size` | 0.02 to 1 | bracket | T3 | L0 | bar length from Erwin 2005 Table 5 and Gadotti 2011 (shape_draw), which draws 0.067 to 1; the floor of 0.02 is the proposal's suggestion (galaxy size levers, 10-21) |
| `--nuclear-amplitude` | 0 to 10 | placeholder | T4 | L0 | hand option of the first plan, no source |
| `--edge` | 4 to 262140 pc, or 3 to 20 L | derived | T2 | L0 | relative: Boss, 10-21 (7:01 PM EDT), 3 to 20 scale lengths; absolute: one 4 pc shell up to the grid's reach of 65,535 shells, 262 kpc (spherical.MAX_ADDRESSABLE_SHELL) |
| `--halo-ratio` | 10 to 3000 | derived | T2 | L0 | the relation's own range: a halo holds 26 times its stars at the peak of star formation efficiency (Moster et al. 2013, read in the repository) and 2200 times at the halo threshold of 1e10 solar masses; the ends are those, rounded |
| `--stellar-halo-share` | 0.001 to 0.3 | bracket | T3 | L0 | the draw's own cut: 0.1 percent (Dragonfly galaxies with no detected halo) to 30 percent (D'Souza and Bell 2014's disc-dominated galaxies), read in structure_draw |
| `--richness` | 1 to 3 systems per solar mass | bracket | T3 | L0 | the derived mean system mass runs 0.47 (old stars) to 0.79 solar masses (stars just formed), 1.3 to 2.1 systems per solar mass (the generator's own laws); the band 1 to 3 widens it by about a factor of 1.4 |
| `--arm-count` | 1 to 10 | placeholder | T4 | L0 | Hart et al. 2017 counts 1 to 5 or more arms; the ceiling of 10 is a placeholder |
| `--pitch-angle-deg` | -90 to 90 degrees, at least 1 away from 0 | placeholder | T4 | L0 | Hart 2017 draws 10 to 40; the full quarter turn is a placeholder; under 1 degree from 0 is refused because a flat pitch divides by zero (density-limits report, 10-21) |
| `--arm-contrast` | 0 to 1 | placeholder | T4 | L0 | 0 to 1 because the density cannot go negative (theory); the 0.4 default is the first plan's hand value, no source |

## 6.6 Fills of the shared download step

| fill_id | kind | tier | tweak | stand-in | replaced_by | sources |
|---|---|---|---|---|---|---|
| SC1 | D | T1 | L0 | Baraffe et al. 2015 track at solar metallicity and fitted relations | the MIST download (`isochrones.py`) | Choi et al. 2016; Dotter 2016 |
| SC2 | D | T1 | L0 | the D13 radius law and a fixed layer split | the sequence download (`white_dwarf_models.py`) | Bédard et al. 2020 |
| SC3 | D | T1 | L0 | the nickel-to-krypton crust table and the D13 radius floor | the equation-of-state and crust download (`neutron_star_matter.py`) | Malik et al. 2024; Servais and Chamel 2026 |
| AD9 | D | T1 | L0 | the 35 band means of the gas table for ozone (20 nm bins; the UV index they give is off by a large factor in the UV-B) | the ozone download (`ozone_uv.py`, read by `physics/uv_cross_sections.py`, GEN.443) | Keller-Rudek et al. 2013 (ESSD 5, 365); Daumont et al. 1992; Malicet et al. 1995; Brion et al. 1998 |

## 6.7 Stand-ins of the atmosphere mixture draw

| fill_id | kind | tier | tweak | constant | knobs | validity | falsifier | replaced_by | evidence |
|---|---|---|---|---|---|---|---|---|---|
| AM1 | B | T3 | L0 | mixture_evolved.WATER_RESIDUAL | 1e-5 to 1e-3, log | water-lost bodies (class N, V rows) | test_atmosphere_mixture_evolved::test_a_venus_like_body_is_carbon_dioxide_with_a_few_percent_nitrogen | GEN.481 volcanic delivery and loss | READ (Venus's 20 ppm of water, a TEST point) |
| AM2 | B | T3 | L0 | mixture_evolved.SULFUR_SINK | 1e-4 to 1e-2, log | Venus-like and Mars-like air | test_atmosphere_mixture_evolved::test_a_venus_like_body_is_carbon_dioxide_with_a_few_percent_nitrogen | GEN.481 | READ (Venus's 150 ppm of SO2, a TEST point) |
| AM3 | P | T6 | L0 | mixture_evolved.WATER_PHOTOLYSED | 0.3 to 1.0, uniform | abiotic-oxygen worlds (row j5) | test_atmosphere_mixture_evolved::test_an_abiotic_oxygen_world_keeps_oxygen_in_its_air | GEN.511 (row weights) and the water-loss model | design |
| AM4 | B | T2 | L0 | mixture_evolved.OCEAN_WATER_WT_PERCENT | 0.0233 to 0.233 weight percent, log | abiotic-oxygen worlds | test_atmosphere_mixture_evolved::test_an_abiotic_oxygen_world_keeps_oxygen_in_its_air | the body water budget of lane J | READ (Luger and Barnes 2015: up to ten oceans) |
| AM5 | B | T3 | L0 | mixture_evolved.CARBON_PER_OCEAN_WATER | 0.01 to 0.2, log | abiotic-oxygen worlds | test_atmosphere_mixture_evolved::test_an_abiotic_oxygen_world_keeps_oxygen_in_its_air | GEN.481 | COMPUTED (Earth's carbon over ocean water 0.054) |
| AM6 | P | T7 | L0 | mixture_evolved.OXYGEN_SINK | 0 to 0.9, uniform | abiotic-oxygen worlds | none: no measurement of the share the crust oxidises away; the 95 percent ceiling test bounds the outcome | the redox and crust model (GEN.481, lane K) | design |
| AM7 | B | T2 | L0 | mixture_evolved.TITAN_HUMIDITY | 0.3 to 1.0, uniform | methane-nitrogen air of a cold body | test_atmosphere_mixture_evolved::test_a_titan_like_air_is_nitrogen_with_methane_at_the_saturation_and_a_hydrogen_trace | a methane cycle model | READ (Niemann et al. 2010: Titan methane humidity about 0.5) |
| AM8 | B | T2 | L0 | mixture_evolved.TITAN_HYDROGEN | 3e-4 to 3e-3, log | methane-nitrogen air | test_atmosphere_mixture_evolved::test_a_titan_like_air_is_nitrogen_with_methane_at_the_saturation_and_a_hydrogen_trace | photochemistry (GEN.496) | READ (Niemann et al. 2010: H2 1.0e-3 +/- 0.16e-3) |
| AM9 | P | T7 | L0 | mixture_evolved.FROST_ICE | CH4 3e-3 to 0.3, CO 1e-3 to 0.1, both log | surface ice of a cold nitrogen body | test_atmosphere_mixture_evolved::test_a_pluto_like_air_follows_raoults_law_over_the_ice | the ice composition item of lane K | design |
| AM10 | D | T6 | L2 | mixture_lava.ALKALI_ACTIVITY_COEFFICIENT | 3e-9 to 3e-8, log | silicate melt, 1500 to 3000 K | test_atmosphere_mixture_lava::test_the_vapour_pressure_over_the_melt_is_kites_within_a_factor_of_four_from_2000_to_2500_k | GEN.499 melt activity model | READ (Kite et al. 2016: vapour over a BSE melt, 2e-4 bar at 2000 K, 8e-3 bar at 2500 K) |
| AM11 | P | T7 | L0 | mixture_lava.COMPOSITION_SCATTER | 0.15 dex per oxide | lava worlds with no stored melt | test_atmosphere_mixture_lava::test_a_stored_melt_composition_replaces_the_drawn_one | the bulk composition of lane J | design |
| AM12 | P | T3 | L0 | mixture_coma.DOMINANCE | rows k1 (CO2 over CO), k2 (CO over CO2), k3 (comparable) | comets and icy asteroids | test_atmosphere_mixture_evolved::test_a_coma_is_water_with_the_drawn_ratios_and_the_class_ordering_of_carbon_dioxide_and_monoxide | the comet composition survey | design |
| AM13 | D | T1 | L0 | mixture_rocky.INVENTORY | wet class water 0.05 to 0.5 wt%, CO2 over H2O 0.2; dry class CO2 0.6 wt% | rocky bodies | test_atmosphere_mixture_rocky::test_the_inventory_follows_the_delivered_water_class | GEN.481 | READ (Ortenzi et al. 2020 Table 1) |
| AM14 | B | T3 | L0 | mixture_rocky.INVENTORY | N/C 1e-3 to 0.1 and S/C 0.1 to 3, log; dry water 1e-4 to 5e-3 wt%; rich water 10 to 50 wt%, CO2 over H2O 5e-3 to 0.2; potassium 0.005 to 0.1 wt%; degassed 40Ar share 0.05 to 0.7 | rocky bodies | test_atmosphere_mixture_rocky::test_argon_grows_with_the_age_of_the_body | GEN.481 (delivery), GEN.492 (radioisotopes) | READ (bulk silicate Earth and CI chondrite ratios) |
| AM15 | B | T2 | L0 | mixture_rocky.MAGMATIC_TEMPERATURE_K | 1100 to 1500 K, uniform | volcanic gas | test_atmosphere_mixture_rocky::test_an_earth_like_oxidised_volcano_gas_is_water_carbon_dioxide_and_sulfur_dioxide_at_the_vent | GEN.482 | READ (part 05 section 3.4: basalt erupts at 1450 to 1550 K) |
| AM16 | B | T2 | L0 | mixture_rocky.QUENCH_TEMPERATURE_K | 400 to 900 K, log | cool bodies | test_atmosphere_mixture_rocky::test_the_equilibrium_temperature_is_the_surface_one_above_the_quench_and_the_quench_below_it | GEN.482 | READ (Zahnle and Marley 2014: 400 to 600 K; Liggins 2022: 800 K) |
| AM17 | B | T3 | L0 | mixture_rocky.OXIDATION_TIME_GY | 0.1 to 1.5 Gy, log | bodies with long-lived geology | test_atmosphere_mixture_rocky::test_the_redox_of_a_young_or_inert_body_is_the_magma_ocean_stage_and_of_an_old_tectonic_one_the_later_mantle | GEN.482 | READ (Aulbach and Stagno 2016: present upper-mantle redox by about 3.5 Ga) |
| AM18 | P | T7 | L0 | mixture_trace.ABUNDANCES | N2O 1e-8 to 1e-6, PH3 1e-10 to 2e-8, HCl 1e-6 to 1e-2, HF 1e-7 to 1e-3, NO2 1e-11 to 1e-8, CH3Cl 1e-11 to 1e-9, all log | rows m3 to m6 and m11 when drawn present | test_atmosphere_mixture_evolved::test_the_trace_gases_follow_the_rows_drawn_after_the_life_stage | the photochemistry and biosignature items | design |
| AM19 | B | T6 | L3 | mixture_escape.EXOBASE_RATIO | cool (CO2 at least half) 1.3, hot 4.0 | the Jeans filter of GEN.693 | test_atmosphere_mixture_escape::test_a_carbon_dioxide_air_has_a_cold_exobase_and_keeps_its_gas_on_a_mars_sized_body | a thermosphere model | fitted to the Solar System as ratios (Mars 1.2, Venus 1.4, Earth 4) |
| AM20 | P | T6 | L0 | mixture_escape.XUV_EXPONENT | 0.25, clamped to 0.5 to 3 times | the Jeans filter | none: no measurement; the clamp keeps the exobase inside the table of AM19 | a thermosphere model | design |
| AM21 | P | T7 | L0 | mixture_escape.ESCAPE_WIDTH | 4.0 in lambda | the Jeans filter | test_atmosphere_mixture_escape::test_the_share_kept_is_smooth_falls_with_the_parameter_and_older_bodies_need_more | hydrodynamic escape model | design |
| AM22 | P | T7 | L0 | mixture_escape.LOST_BELOW | 1e-3 of a species | the Jeans filter | test_atmosphere_mixture_escape::test_a_small_warm_body_loses_everything_it_would_outgas_and_is_airless | the same | design |
| AM23 | P | T7 | L0 | mixture_escape.LOST_TOTAL | 1e-2 of the gas | the Jeans filter | test_atmosphere_mixture_escape::test_a_small_warm_body_loses_everything_it_would_outgas_and_is_airless | the same | design |
| AM24 | P | T7 | L0 | mixture_evolved.TITAN_MAX_METHANE | 0.5 of the gas | nitrogen-methane air when the saturation exceeds the pressure | test_atmosphere_mixture_evolved::test_a_titan_like_air_is_nitrogen_with_methane_at_the_saturation_and_a_hydrogen_trace | a methane cycle model | design |
| AM25 | P | T7 | L0 | mixture_oxygen.INERT_FLOOR | 0.5 percent nitrogen plus argon | any air with oxygen above a few percent | test_atmosphere_mixture_escape::test_a_lifeless_oxygen_world_is_capped_at_95_percent_with_a_nitrogen_argon_floor | the volcanic delivery and noble-gas inventory (GEN.481) | design |
| AM26 | B | T3 | L0 | mixture_rocky.INVENTORY | weathering_sink 1e-5 to 1e-2, log, on the carbon and sulfur of a body with an ocean and a long-lived geology | rocky bodies with `ocean_present` and `long_lived_geology` | test_atmosphere_recipe_gaps::test_an_ocean_and_a_long_lived_geology_leave_little_carbon_and_sulfur_and_all_the_nitrogen_and_argon | GEN.481 (weathering feedback) | READ from memory, not re-checked (Earth keeps about 1.5e-5 of its near-surface carbon in the air: 0.9e15 kg C in the CO2 against about 6e19 kg C in the carbonates) |
| AM27 | D | T5 | L1 | surface_pressure.PROCESSING_ENERGY_J_PER_KG | 9.17e6 J per kg of mantle, fixed | rocky bodies with a stored interior (`interior_heat`) | test_surface_pressure::test_venus_degasses_about_what_the_earth_does_and_its_air_asks_for | the melt-production model of Tosi et al. 2017 (equations 15 to 19) | CALIBRATED on the Earth and Venus surface carbon, 1e22 mol each (source tier 5, tweak level 1; Tosi et al. 2017 picture) |
| AM28 | B | T6 | L1 | surface_pressure.DEGAS_EFFICIENCY_H | 0.05 to 0.5, log | rocky bodies | test_surface_pressure::test_the_inventory_degassed_is_each_element_times_its_share_and_its_efficiency | the water solubility of the melt (Sossi et al. 2022) in a melt-production model | design (source tier 6, tweak level 1; Marty 2012 up to 90 percent of the water in the mantle, Bower et al. 2022 magma-ocean retention) |
| AM29 | P | T6 | L0 | surface_pressure.DEGAS_EFFICIENCY | C, N, Ar 1; S 0.5 | rocky bodies | test_surface_pressure::test_the_inventory_degassed_is_each_element_times_its_share_and_its_efficiency | sulfur solubility in basalt (a melt model) | design (source tier 6, tweak level 0; sulfide stays in the melt) |
| AM30 | B | T1 | L0 | surface_pressure.WEATHERING_ALPHA | 0.1 to 0.5, uniform | bodies with liquid water at the surface | test_surface_pressure::test_the_cycle_draws_stay_in_their_brackets_and_are_recorded | a weathering model with runoff and land (GEN.726 follow-up) | READ (Lehmer et al. 2020 Table 1: 0.1 to 0.5; source tier 1, tweak level 0) |
| AM31 | B | T1 | L0 | surface_pressure.WEATHERING_T_E_K | 10 to 40 K, uniform | bodies with liquid water at the surface | test_surface_pressure::test_the_cycle_draws_stay_in_their_brackets_and_are_recorded | the same | READ (Lehmer et al. 2020 Table 1: 10 to 40 K; source tier 1, tweak level 0) |
| AM32 | B | T1 | L0 | surface_pressure.OUTGASSING_EXPONENT | 1 to 2, uniform | bodies with liquid water at the surface | test_surface_pressure::test_the_cycle_draws_stay_in_their_brackets_and_are_recorded | the same | READ (Lehmer et al. 2020 Table 1: 1 to 2; source tier 1, tweak level 0) |
| AM33 | D | T5 | L1 | surface_pressure.HEAT_FLUX_EARTH_W_M2 | 0.0614 W per m2, fixed | the reference of the outgassing flux | test_surface_pressure::test_the_earth_is_the_balance_at_its_modern_flux_land_and_temperature | a mobile-lid heat loss in `interior_heat` | CALIBRATED on the interior model's own Earth twin, mid rheology (source tier 5, tweak level 1) |
| AM34 | P | T6 | L0 | mixture_subject.PHOTOSPHERE_BAR | 0.1 bar, fixed | the regime of hot giants and brown dwarfs (rows h and g read the edge at this pressure) | test_mixture_subject::test_a_giant_gets_its_layer_temperatures_and_the_photosphere_pressure | the cooling-track and atmosphere model of the giant lane | design (the recipe quantity's own definition: about 0.1 bar for the dayside emission of a hot giant) |
| AM35 | D | T5 | L0 | mixture_subject.DAYSIDE_OVER_EQUILIBRIUM | 1.278, fixed | every body the draw reads a dayside temperature of (hot giants, lava worlds) | test_mixture_subject::test_the_earth_gives_the_quantities_the_world_gives | a heat-redistribution model | COMPUTED (a day hemisphere alone catches the starlight: f = 2/3 against 1/4, so T rises by (8/3)^(1/4); the upper bound of the dayside, no redistribution) |

## 6.8 Stand-ins of the orbit elements

| fill_id | kind | tier | tweak | constant | knobs | validity | replaced_by | evidence |
|---|---|---|---|---|---|---|---|---|
| OE1 | B | T1 | L0 | orbit_elements.MULTI_SIGMA | 0.083, fixed | planets of systems with two or more planets | a law fitted to planets of every period and size | READ (Van Eylen et al. 2019 abstract: multi-planet systems consistent with a half-Gaussian of 0.083, +0.015 -0.020) |
| OE2 | B | T1 | L0 | orbit_elements.SINGLE_MODERATE_FRACTION | 0.76, fixed | a system with one planet | the same | READ (the same abstract: 0.76, +0.21 -0.12, of single transiting systems have a moderate eccentricity) |
| OE3 | B | T1 | L0 | orbit_elements.SINGLE_RAYLEIGH_MODE | 0.26, fixed | the moderate share of single planets | the same | READ (the same abstract: a Rayleigh law that peaks at 0.26, +0.04 -0.06) |
| OE4 | T | T2 | L0 | orbit_elements.draw_planet_eccentricity | the three above | every planet | radial-velocity and direct-imaging eccentricities by period and mass | the laws of OE1 to OE3 used outside the planets they were fitted to (Kepler's, mostly under 100 days) and with the planets the generator counts, not the ones a transit survey finds |
| OE5 | P | T4 | L0 | orbit_elements.MAX_ECCENTRICITY | 0.9, fixed | every planet and moon | none needed | design (the laws put almost nothing beyond it; the periapsis cut comes first for a close orbit) |
| OE6 | B | T2 | L0 | orbit_elements.MUTUAL_INCLINATION_SIGMA_DEG | 1.6 degrees, fixed | planets | a law of the plane scatter of giants and of systems beyond the Kepler range | Fabrycky et al. 2014: Rayleigh mutual inclination of 1.0 to 2.2 degrees (read in an earlier session, not again for this item); a planet's own angle to the plane is that over the square root of two |
| OE7 | P | T3 | L0 | orbit_elements.REGULAR_MOON_ECCENTRICITY_RANGE | 0.001 to 0.01, log | moons of disk, impact, ring or binary origin | a forced eccentricity from the resonances and neighbours of the moon system | the forced eccentricities of the regular satellites of the giants run from 0.001 to 0.03 and the Moon 0.055; the range is the one the tide of GEN.700 already drew from |
| OE8 | T | T5 | L0 | orbit_elements.draw_moon_eccentricity | the Hill cut of `tuning.MOON_PROGRADE_STABLE_HILL_FRACTION` | captured moons | capture and evolution of an irregular moon | the thermal distribution of eccentricities, f(e) = 2 e, of a body captured on a wide orbit (Jeans 1919, recalled), cut at the stable part of the parent's Hill sphere |
| OE9 | P | T3 | L0 | orbit_elements.GIANT_TIDAL_Q, GIANT_TIDAL_LOVE_NUMBER_K2 | 1e5 and 0.3, fixed | the tide of a gas giant | per-body Q and k2 (GEN.731) | Jupiter: k2 about 0.38 and Q from 3.5e4 to about 1e6 from the Galilean moons' orbits (recalled; GEN.731 gives every body its own) |
| OE10 | T | T5 | L0 | orbit_elements.circularisation_time_yr | Q 100, k2 0.03 for rock and ice | every planet and moon | per-body Q and k2 (GEN.731) | the circularisation time of a synchronously rotating body, tau = (2/21)(Q/k2)(m/M)(a/R)^5/n (Goldreich and Soter 1966), with the rocky Q and k2 of `constants.MOON_TIDAL_DISSIPATION_Q` and `MOON_TIDAL_LOVE_NUMBER_K2` |
| OE11 | T | T6 | L0 | insolation.companion_swing | in phase, fixed | planets of a binary | a time series of the sum | the swing of the planet's own orbit and the companion's are taken in phase, the largest the pair can give; their periods differ and the extremes meet only now and then |

## 6.9 Stand-ins of the tidal heating

| fill_id | kind | tier | tweak | constant | knobs | validity | replaced_by | evidence |
|---|---|---|---|---|---|---|---|---|
| TH1 | P | T2 | L0 | tidal_heating.ROCK_RIGIDITY_PA | 7e10, fixed | the solid mantle of a rocky body | a rigidity from the body's composition and temperature | the shear modulus of the mantle in PREM runs from 65 to 80 GPa (recalled) |
| TH2 | T | T2 | L0 | tidal_heating.rigid_k2 | rigidity only | a rocky body without a large liquid core | a layered Love number | the Love number of a homogeneous, incompressible elastic sphere, k2 = (3/2) / (1 + 19 mu / (2 rho g R)) (Love 1911, recalled); it gives the Moon 0.021 against the observed 0.0242 |
| TH3 | P | T3 | L0 | tidal_heating.ROCK_TIDAL_Q | 100, fixed | solid rock | a Q from the body's temperature and frequency | Q of 100 for solid rock (Goldreich and Soter 1966, recalled), the same value that circularises the orbit (OE10); the Moon's own is 38 |
| TH4 | P | T3 | L2 | tidal_heating.ICY_K2_OVER_Q | 3e-4, fixed | an icy body (a water mass fraction of 0.01 or more) | k2/Q from the shell thickness and the ocean beneath it | tuned to the published figure: Europa radiates 1.4e11 W of tide at e = 0.0094 with k2/Q = 3e-4 (Peale et al. 1979 applied to the ice by Ojakangas and Stevenson 1989, recalled) |
| TH5 | P | T8 | L4 | tidal_heating.MOLTEN_K2_OVER_Q | 0.015, fixed | a melt-limited rocky body | the Maxwell rheology of the partly molten mantle in the box | fitted to the one body that shows it: Io radiates 0.9 to 1.5e14 W (Veeder et al. 2004, recalled) with k2 = 1.5 and Q = 100 |
| TH6 | P | T2 | L0 | tidal_heating.SOLIDUS_K, LIQUIDUS_K | 1700 and 2200, fixed | the one-box mantle temperature | a pressure-dependent melting curve | the peridotite solidus runs from 1394 K at the surface to 1750 K at 3 GPa (Hirschmann 2000, recalled); the liquidus lies 400 to 500 K above |
| TH7 | P | T3 | L3 | tidal_heating.HEAT_PIPE_MELT_FRACTION | 0.3, fixed | a rocky body whose tide warms its mantle | a resurfacing-rate model | tuned against Io: melt begins to segregate at a melt fraction of 0.2 to 0.4 (the rheological transition, Arzi 1978, recalled; Moore 2003 puts Io's mantle near 0.2 to 0.3) |
| TH8 | T | T6 | L0 | interior_heat.heat | constant, fixed | every body with a stored eccentricity | the eccentricity history | the tide is taken constant over the age at the stored eccentricity, though a planet that has circularised had more of it in the past |
| TH9 | T | T6 | L1 | interior_heat.heat | the two conditions | rocky bodies | a heat-pipe loss law in the box | a body is melt-limited when its solid-state tide is at least its own heat and warms the mantle to the heat-pipe temperature; the real transition is gradual and is bistable in a Maxwell box, so the state is chosen, not integrated |

## 6.10 Stand-ins of the spin and its despinning

| fill_id | kind | tier | tweak | constant | knobs | validity | replaced_by | evidence |
|---|---|---|---|---|---|---|---|---|
| SP1 | P | T2 | L0 | spin.ACCRETION_SPIN_SIGMA_PER_HOUR, ACCRETION_CRITICAL_PER_HOUR | 1.26 to 1.49 (their three components); 2.15 for their perfect-accretion model | a rocky planet of 0.1 to 1.2 Earth masses built by giant impacts | a spin law that includes the later impacts and the planetesimal stage | READ through a summary of Kokubo and Genda 2010 (ApJ 714:L21, arXiv:1003.4384): each component of a planet's birth spin is Gaussian, 1.32, 1.26 and 1.49 rad/h with their realistic accretion model, rms 2.35 rad/h, 0.71 of the breakup rate at their 3000 kg/m3; paper not opened in full |
| SP2 | T | T2 | L0 | spin.accretion_spin | isotropic, fixed | a rocky planet | a measured obliquity distribution | the same paper: the obliquity is isotropic, n(e) = sin(e)/2, which comes from components of one width; it replaces the impact-tilt law with gamma 3, which was planetGen's choice |
| SP3 | T | T7 | L0 | spin.accretion_spin | scale by breakup, fixed | a rocky planet of any density | the paper's mass and density scan | the paper finds the mean spin independent of mass at one density; a planet of another density is scaled by its own breakup rate sqrt(G M / R^3), which is assumed, not read |
| SP4 | P | T6 | L0 | spin.BREAKUP_FRACTION | 0.7 to 1 | a rocky planet | a source on the spin of a planet at the end of accretion | a birth spin is held at 0.85 of the breakup rate, the cap the hot stars use (planetGen's choice); about 14 percent of the Gaussian vectors lie above it |
| SP5 | P | T3 | L1 | constants.ROTATION_PERIOD_RANGE_HOURS | 8 to 48 and 8 to 20 | moons; gas giants | a circumplanetary-disk spin law for moons; a giant spin-up law | a moon's birth spin is log-uniform over 8 to 48 hours (the range this generator has used for rocky bodies since GEN.104, whose source was the impact papers), a gas giant's uniform over 8 to 20 hours (Jupiter 9.9, Saturn 10.7, Neptune 16.1, Uranus 17.2 h, recalled); the tide then moves both |
| SP6 | P | T2 | L0 | tidal_spin.ROCKY_MOMENT_FACTOR, GIANT_MOMENT_FACTOR | 0.25 to 0.4 | rocky bodies; gas giants | a factor from the body's core fraction | the moment-of-inertia factor I / (m R^2): Earth 0.331, Venus 0.33, Mars 0.365, Mercury 0.346, the Moon 0.393 (rock, 0.33); Jupiter 0.254, Saturn 0.22 (giants, 0.25); recalled |
| SP7 | T | T3 | L0 | tidal_spin.despin_rate_rad_s2 | none | a body spun down to the orbital rate | a frequency-dependent rheology (Efroimsky's) | the despinning torque is the constant-Q one of Gladman et al. 1996 (Icarus 122:166, recalled), independent of how fast the body turns; a real mantle's response depends on the frequency, which only the resonance windows (SP8) carry |
| SP8 | P | T2 | L0 | tidal_spin.RESONANCES | the window edges | a rocky planet brought down from a fast spin | the same windows for a layered planet | READ through a summary of Makarov 2012 (ApJ 752:73, arXiv:1110.2658): with a frequency-dependent rheology capture is certain in 0.34 < e < 0.49 for 2:1, 0.2 < e < 0.41 for 3:2 and e < 0.29 for 1:1; the first window met coming down takes the spin |
| SP9 | P | T3 | L0 | tidal_spin.RESONANCE_CAPTURE_PROBABILITY | 0.07 to 1 | Mercury-like captures | a draw once one model is chosen | capture probability inside a window is 1 (Makarov 2012: 40 of 40); Goldreich and Peale 1966 give 0.07 at Mercury's eccentricity with a constant-torque tide and Correia and Laskar 2004 give 0.55 with the eccentricity's chaotic history; the sources split, the bracket is kept |
| SP10 | T | T6 | L0 | planets.generate_spin | none | every planet and moon | a coupled spin, obliquity and eccentricity history | approximations of the spin state: a locked body turns exactly once an orbit (not 1 + 6 e^2 of a constant-time-lag tide); the eccentricity is the stored one, not the larger past one; the obliquity stays as drawn until the spin has come all the way down, then is zero; no atmospheric tide, core-mantle coupling or the tide of the body's own moons |
