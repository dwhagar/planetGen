# Atmosphere composition research: compounds, gases, condensates, radiation and the stats to compute with them

**Status (2026-10-17):** research result of GEN.431 (the Atmosphere research lane), merged here by the TODO thread so that every lane has it in the repository. It feeds GEN.409 (the checked atmosphere colour design), GEN.410 (the gas and compound table), GEN.436 (the radiation design) and the generation chain of [atmosphere-colour.md](atmosphere-colour.md); [composition-model.md](composition-model.md) is the master note. Times in it are as the lane's clock reported. Where it says a source was not reached, the input is open and is Boss's to research (section 6); each such input carries a seeded physical range meanwhile, never a fixed value.

Date: 2026-10-17 (written between 12:03 AM and 1:00 AM EDT as the session clock reported; lane clocks drift). Research lane, read-only on code and docs (the TODO thread merged the result into the repository). Request: Boss, 12:00 AM EDT (the message the coordinator relayed), atmosphere composition as a joint headline of the next major release with body composition; every body gets composition figures from real science, drawn from a bulk density and pressure first and then filled in, or built up from basics. Chain: goal issue #1944, plan note `docs/design/atmosphere-colour.md`, baseline reports `sky-color-analysis-2026-10-17.md`, `atmosphere-sources-check-2026-10-17.md`, `composition-density-link-2026-10-17.md` under `/mnt/project-files/reports/`.

This file is the master section. The detail (about 440 KB of tables, formulas and source lists) sits in seven part files in `docs/design/atmosphere-composition-parts/`, copied unchanged from the lane's draft (`/mnt/project-files/design-drafts/atmosphere-parts/`).

## 1. What was delivered

| Part | File | Size | What it holds |
|---|---|---|---|
| 01 | [01-gas-table-thermodynamics.md](atmosphere-composition-parts/01-gas-table-thermodynamics.md) | 61 KB | 72 gases and volatile compounds plus metal hydrides and NH4SH: formula, CAS, molar mass, flags, photolysis threshold, triple and critical points, heat capacity at three temperatures, gamma, latent heats, saturation vapour pressure curves with errors at the boiling and triple points, Henry constants for 35 species, viscosity, real-gas Z, mixture mean molar mass and density, a checked replacement for the code's two-point condensation table, a 10-row test table |
| 02 | [02-gas-optics-uv-visible.md](atmosphere-composition-parts/02-gas-optics-uv-visible.md) | 64 KB | Rayleigh cross sections for 28 gases (13 verified, 15 from one polarizability formula good to -9 to +6 percent), depolarization ratios with real spreads, UV/visible absorption of 22 gases in 35 bins with temperature trends, alkali and metal-oxide band data for hot giants, surface UV from the ozone column, UV index and erythemal and DNA action spectra, air mass and Chapman function, ground-bounce formula, blackbody chromaticity 2000 to 30000 K, and the re-check of the unreached sources |
| 03 | [03-infrared-greenhouse-thermal.md](atmosphere-composition-parts/03-infrared-greenhouse-thermal.md) | 65 KB | Greenhouse gas bands and forcing laws, collision-induced and water continuum absorption, heat capacity and lapse rates, analytic grey and two-band temperature structure, the inverse for albedo and greenhouse multiplier, Rayleigh and cloud albedo, heat redistribution, runaway and moist greenhouse limits, a split between what can be code and what needs full radiative transfer |
| 04 | [04-photochemistry-lifetimes-escape.md](atmosphere-composition-parts/04-photochemistry-lifetimes-escape.md) | 59 KB | A starter reaction network with rate constants, photolysis lifetimes of about 35 species for F, G, K and M hosts, UV shielding, ozone column and layer depth, escape thresholds and formulas, stellar XUV laws, exospheres of airless bodies, flares, SEPs and lightning |
| 05 | [05-mixing-ratios-by-class.md](atmosphere-composition-parts/05-mixing-ratios-by-class.md) | 54 KB | The solar element budget, mass-metallicity draw rule, outgassed rocky-planet atmospheres by mantle redox, lava worlds, sub-Neptunes and hycean worlds, Venus-, Mars-, Titan-, Pluto-like frost-equilibrium atmospheres, giants and brown dwarfs, hot Jupiters, comets, white-dwarf and pulsar planets, host-star dependence, and a species-by-class table |
| 06 | [06-condensates-clouds-hazes-dust-minerals.md](atmosphere-composition-parts/06-condensates-clouds-hazes-dust-minerals.md) | 66 KB | Condensation temperature fits for oxides, silicates, iron, sulfides and salts, sublimation of about 19 ices, cloud-base rule, optical table for 47 materials (density, refractive index, single-scattering albedo, colour slope from Mie runs), haze yields and sizes, dust lofting and settling, cloud albedo against optical depth, lava-world vapours and mineral densities |
| 07 | [07-radiation-shielding-and-ground-sources.md](atmosphere-composition-parts/07-radiation-shielding-and-ground-sources.md) | 62 KB | EUV and X-ray absorption by composition, photon attenuation for 13 materials, proton ranges (NIST PSTAR) for 10 materials, muon and secondary data, dose against atmospheric column, cutoff rigidity, flare and SEP scaling, radioisotopes in the ground (K, U, Th, Rb, Al-26, Fe-60), gamma dose conversion, radon, Ar-40 and He-4 budgets, radiolysis, W values, dose limits, a dose recipe with code status |

Row counts by verification label, summed from the seven parts (each worker's own count; a row with two labels counts under each): READ 399, DATA 479, COMPUTED 533, 2ND-HAND 118, NOT REACHED 46, RECALLED 34 (all also listed in the parts' to-verify sections), SNIPPET 9. Part 02 adds about 770 computed or data cells for the absorption bins that are not counted here.

Labels: READ is the original paper's own text; DATA is a data file in a package or repository (name and version given); COMPUTED is re-derived in Python from stated inputs; 2ND-HAND is a number read in a paper or code that cites the original; RECALLED is memory and unverified; SNIPPET is a search summary and not relied on; NOT REACHED is needed but no reachable source gave it.

## 2. Checks by this lane on the worker output

The parts are worker output; these are the checks the lane made itself, not a full audit.
- All 72 molar masses in the part 01 species table were recomputed from the `periodictable` package: 70 agree to 0.01 g/mol; the two exceptions are 3He (an isotope mass, my parser read the formula as a carbon-12 compound) and Rn (the isotope-222 mass 222.0176, as labelled).
- The part 07 finding that the dipole cutoff rigidity scales as M/r^2 (not r^3) agrees with the Stormer vertical cutoff Rc = 14.9 GV / L^2 for Earth, which depends on the moment over the square of the distance.
- The part 01 finding that the code's H2S condensation row is a copy of the SO2 row is true in `src/planetgen/physics/atmosphere.py` (`CONDENSATION_K["h2s"] = CONDENSATION_K["so2"]`), and H2S boils at about 213 K at 1 atm, 50 K below SO2.
- The earlier verified values (Rayleigh sigma of N2, CO2/N2 ratio 2.51, Earth tau_R 0.0969) reappear in part 02 unchanged.
- The vapour-pressure, optical-constant and Mie tables were not independently recomputed. TODO should have each item's checked design (GEN.409 and the items it spawns) sample a few rows against the cited file before coding a table.

## 3. Findings that change the existing design or code

| # | Where | Finding | Source part |
|---|---|---|---|
| 1 | `atmosphere.py` `CONDENSATION_K` | Two-point Clausius-Clapeyron rows are off: the H2S row is a copy of SO2 (20 to 50 K wrong), CO2 up to 7 K warm, SO2 up to 16 K cold, N2, CO, Ar, O2 and CH4 up to 5 K cold. Rows stay acceptable only for H2, N2, CO, Ar, O2, CH4, H2S (own curve needed) and C2H6. A checked replacement with temperatures at seven pressures is given for 14 entries including NH3, H2S, NH4SH, HCN, C2H6 and H2SO4. Fray and Schmitt sublimation fits exceed the lab curves by up to 6 times below 40 K | 01, 06 |
| 2 | `atmosphere.py` `REDOX_EARTH_DELTA_IW = 3.5` | The mantle redox prior is far more oxidised than meteoritic building blocks (about 1 below to 0 above the iron-wustite buffer) and than a magma-ocean stage (near IW to IW+1); a young rocky planet should start reduced and oxidise with time and mass, not start at Earth's present mantle. White dwarf debris shows bodies near IW-1 | 05 |
| 3 | Rayleigh table (GEN.410, GEN.417) | The older coded CH4 refractivity form is about 30 percent above the measured value (He 2021 gives the replacement); the H2/N2 ratio spreads from 0.202 to 0.217; N2 depolarization spreads from 0.020 to 0.031 (2 percent in sigma); the water-vapour refractivity has no known reference temperature, so use the generic polarizability route; 15 further gases are covered by one formula with a stated error of -9 to +6 percent | 02 |
| 4 | `radiation.py` and the radiation design | The dipole cutoff scales as M/r^2 where the code uses (R_E/R_p)^3; the GCR dose-against-column curve is 1.5 to 2.7 times below proton-primary shower simulations; the flat 0.6 magnetic factor is wrong in thin atmospheres; the range law `R = 0.0022 E^1.77` holds below 300 MeV but is 5 times too long at 10 GeV; the Europa 5.4 Sv/day figure has no primary source; the 0.64 mSv/day Mars figure is a snippet in the design file and the replacement is Guo 2015 (read) | 07 |
| 5 | Greenhouse multiplier | The existing multiplier is exactly (T_s/T_e)^4 - 1, so every composition-derived result can be reported in the repo's own variable (Earth 0.625, Venus 99.6, Titan 0.69, computed). The inverse for albedo and infrared optical depth tau0 works for those three test sets. The mean molar mass of the mixture replaces `atm_molar_density` in the scale height and in tau0. Infrared optical depth from line absorbers (CO2, H2O, CH4) needs correlated-k grids that are not reachable, so it stays a parametrisation until sourced | 03 |
| 6 | Cloud-species constants | The Mg2SiO4 saturation constant printed in a standard source is likely a misprint (4.88 printed, 14.88 fits); the SiO2 constant is labelled "made up" in its own code; the FeS optical row looks suspect | 06 |
| 7 | Atmosphere document references | "Guzman-Marquez 2021" in the imported Atmospheric Toxicity document was not found; Harman 2018 (ApJ 866, 56) is real | 05 |
| 8 | UV index | The repo's `uv_index` is not the CIE index; part 02 gives the CIE erythemal and Setlow DNA formulas for a proper one | 02 |
| 9 | Ozone | Ozone column against O2 and host star: Kozakis fit stays (earlier report); Chapman layer depth and peak pressure against O2, and abiotic O2/O3 channels, are in part 04; the fit is valid only for 0.01 to 150 percent of the present level | 04 |
| 10 | The ten stored species | The stored set (O2, CO2, CO, N2, Ar, H2, H2O, CH4, H2S, SO2) misses He, Ne, Kr, Xe, NH3, NH4SH, HCN, C2H6, C2H2, OCS, HCl, HF, N2O, NO2, O3, Na, K, SiO, Mg, Fe, TiO, VO and PH3, which the class mixes and descriptions already name or which the physics needs; the species-by-class table in part 05 section 13 lists the add list with the class that needs each | 05 |

## 4. Design synthesis: one registry, many consumers

A single gas and compound registry (the GEN.410 table) should hold one row per species with these fields, each filled from a part; fields with no reachable source carry a NOT REACHED flag and a seeded range rather than a fixed value.

| Field group | Fields | Part |
|---|---|---|
| Identity | name, formula, CAS, molar mass, isotopologue flags, flags (greenhouse, condensible, photochemically fragile, toxic) | 01 |
| Phase | triple point, critical point, normal boiling or sublimation point, saturation vapour pressure curve with validity range, latent heats | 01, 06 |
| Thermal | heat capacity at three temperatures (or polynomial), gamma, viscosity, Henry constant, real-gas Z | 01, 03 |
| Scattering | Rayleigh cross section (polarizability, depolarization, reference temperature of the formula) | 02 |
| UV and visible absorption | coarse bins with temperature trend, band positions for hot-giant absorbers | 02 |
| Infrared | band centres, strength class, forcing law, collision-induced pairs | 03 |
| Photochemistry | photolysis threshold, lifetime by host star class, reaction rows | 04 |
| Escape | mass for the thermal escape test, diffusion-limit data | 04 |
| High energy | EUV and X-ray photoabsorption, mass attenuation, proton stopping, ion-pair energy W | 07 |
| Distribution | median mole fraction and log range by body class, redox, temperature, host star | 05 |

Condensed species (the second registry, GEN.414's table): name, formula, density, condensation or sublimation curve and metallicity term, refractive index n and k by wavelength, single-scattering albedo, Angstrom exponent, particle size range by process (part 06). Radioactive ground species: half-life, decay energy, abundance ranges in crust and mantle, gamma dose coefficient, host-star scaling (part 07).

How the parts feed the chain items of #1944:

| Chain item | Reads from |
|---|---|
| GEN.410 gas and compound table | 01 (all), 02 (Rayleigh, UV bins), 03 (bands), 05 (species-by-class), 06 (condensates) |
| GEN.411 molar density from partial pressures, GEN.424 balance gas | 01 (mean molar mass, scale height, real gas), 05 |
| GEN.412 albedo and greenhouse recorded with an inverse; GEN.423 | 03 (inverse, tau0, Rayleigh and cloud albedo), 01 |
| GEN.422 freeze-out ordering | 01 and 06 (checked curves), 05 |
| GEN.413 ozone | 04 section 3, 02 (cross sections, Chapman), 05 |
| GEN.414 haze, dust, cloud | 06, 04 (haze yields), 03 (anti-greenhouse) |
| GEN.417 Rayleigh of the mixture | 02 |
| GEN.418, GEN.419 colour of solid worlds and giants | 02 (chromaticity), 06 (cloud and haze optics), 05 (giant chemistry), 03 |
| Existing radiation and escape items (GEN.87, GEN.91 shoreline) | 07, 04 |

New work these parts support, and where the TODO thread placed it (2026-10-17; version words are Boss's, owners are proposals):
1. Atmospheric radiation: EUV and X-ray column absorption by composition, ion production and aurora (07 section 2), escape-relevant XUV heating (04): GEN.443 (9.0 phase 2).
2. Surface radiation from composition: GCR and SEP dose against the column and gas, ground gamma from K, U, Th drawn from the host star, radon, UV index (07, 02): GEN.444 and GEN.445 (9.0 phase 2); the corrections of the existing constants are GEN.453.
3. Ground radioisotope abundances and radiogenic heating drawn from host-star composition (07, linked to the bulk composition chain): GEN.444 with GEN.439.
4. Outgassing speciation from mantle redox, volatile budget and temperature, with a reduced start (05): GEN.451 and GEN.411 (8.4 phase 2).
5. Giant and brown-dwarf chemistry sequences by temperature, quench and metallicity (05), hot-Jupiter thermochemistry (05): GEN.454 (8.4 phase 2).
6. Exospheres of airless bodies, comet coma abundances and production rates (04, 05): GEN.455 (8.4 phase 2).
7. Dust lofting and settling by pressure and gravity (06): GEN.414 (8.4 phase 3).
8. Lava-world and hot-rocky vapours (05, 06): GEN.456 (8.4 phase 3).
9. A real UV index and biological UV doses (02): GEN.445 (9.0 phase 2).
10. The checked saturation curves replacing the two-point table (01, 06): GEN.452 (8.4 phase 2); the added species and their storage: GEN.410 and DB.71 (8.4 phase 2).

## 5. Sources and domains

At the time the lane started, only export.arxiv.org, ar5iv, raw.githubusercontent.com, PyPI and search snippets worked. At about 12:03 AM EDT the network was widened. Tested by this lane and by the workers:

| Reachable | Refused |
|---|---|
| arxiv.org (abstracts and PDFs), export.arxiv.org, NIST WebBook, NIST XCOM, PSTAR and ASD, HITRAN (pages), ExoMol (pages), PubChem, NASA NTRS, Wikipedia, Copernicus (AMT, ACP), Optica (page and citation, PDF behind a captcha), IOPscience, Nature, Springer, Crossref, Zenodo, OSTI, UNSCEAR, PDG, NNDC, ICRP, CDS and VizieR, IPCC, PDS, archive.org (catalogue, lending 401), Mainz UV/VIS atlas, IUP Bremen ozone page, webmineral.com, gitlab.com, par.nsf.gov, europepmc | doi.org, science.org, sciencedirect.com, agupubs and Wiley, aanda.org, academic.oup.com, cambridge.org, pubs.acs.org, researchgate, USGS, Mindat, lpi.usra.edu, github.com web pages (raw files work), arXiv API (rate limited), JPL data evaluation (502), HITRAN line-by-line (login), several publisher captchas |

Re-check of the sources that the first report could not reach (part 02 section 2.10):

| Source | Result now |
|---|---|
| Peck and Huang 1977 (H2 dispersion) | Optica gives only the citation. The dispersion constants (0.0148956, 180.7, 0.0049037, 92 at 273.15 K, 0.168 to 1.6945 um) were read in the refractiveindex.info database file `H2_nk_Peck.yml` (CC0). DATA for the constants; paper not reached |
| Burkholder and Talukdar 1994 | Mainz atlas files read for 220, 240, 260, 280 and 298 K, 407 to 762 nm; peak 5.163e-21 cm2 at 602 nm. The atlas says the cross sections were rebuilt from ratios the authors supplied in 2005. DATA; paper not reached |
| Serdyuchenko 2014 | Bremen data read: 11 temperatures from 193 to 293 K, 213 to 1100 nm, 0.01 nm grid, systematic uncertainty under 0.5 percent. The page states no licence. DATA; licence not reached |
| Chandrasekhar 1960 | Not reachable (lending 401). The ground-bounce factor R_tot = R + T A (1 - s) / (1 - A s) was re-derived and checked by Monte Carlo. Formula COMPUTED |
| Sneep and Ubachs 2005, Bates 1984, Kasting and Donahue 1980, Segura 2003, Kasten and Young 1989 | Still second-hand (ScienceDirect and Cloudflare pages); He 2021 supports the Sneep and Ubachs form and measured cross sections for O2 and CO2 at 532 nm agree with it (CO2 13.3e-27 cm2 measured against 13.32e-27 n-based) |
| Bucholtz 1995 Table 3, Penndorf 1957, Hulburt 1953 | Citations confirmed; text behind paywall or captcha; still not reached. The Earth tau_R test value stays computed |

The first report's open point on CO2 at 532 nm (computed 13.3e-27 cm2 against a snippet 12.4e-27) is closed in favour of 13.3e-27 by the measured value in He 2021.

## 6. What Boss must research (consolidated)

Ranked by how much the chain depends on the item. The parts hold the full lists (part 01 section 13, 02 section 4, 03 section 11, 04 section 8, 05 section 14, 06 section 4, 07 section 5).

1. **Fray and Schmitt 2009** (Planet. Space Sci. 57, 2053), with Grundy et al. 2024 and Jager and Span 2012: validity ranges and coefficients of every cold sublimation curve (N2, CO, CH4, CO2, SO2, H2S, NH3). The ScienceDirect copy is refused and no arXiv copy exists; the current numbers come from a supplementary table of arXiv:2009.02277 with two typos. Needed by the freeze-out ordering and by every cold world.
2. **Magma-ocean redox and carbon solubility:** Sossi et al. 2020 (Sci. Adv.), Gaillard et al. 2014 and 2021, Hirschmann 2012. Needed to set the redox prior (finding 2 above) and the CO against CO2 split.
3. **Methane, ammonia and water absorption strengths in the visible for the giants** (Karkoschka 1994; Karkoschka and Tomasko 2009, 2010), and the visible NH3, SO2, NO2 tails.
4. **Mars dust and Titan haze optical spectra:** single-scattering albedo and imaginary index of Mars dust (Wolff 2009, Ockert-Bell 1997, Lemmon 2004, Tomasko 1999) and Titan haze absolute depth (Tomasko 2008, Doose 2016). Part 06 gives a flagged parametrisation meanwhile.
5. **Line and continuum data for composition-derived infrared optical depth:** correlated-k grids or Rosseland and Planck means for CO2-dominated, N2-O2-H2O and steam atmospheres (HITRAN, HITEMP or ExoMol with CIA; Pierrehumbert 2010 for the water side; Myhre 1998 original and Parmentier and Guillot 2014 coefficient formulas). Without them the greenhouse stays a parametrisation (finding 5).
6. **Rayleigh originals:** Bogaard et al. 1978 depolarization of SO2, H2S and others; Sneep and Ubachs 2005; Bates 1984; Bucholtz 1995 Table 3. They remove the plus or minus 5 to 9 percent generic error.
7. **Cosmic-ray and energetic-particle inputs:** Tylka and Dietrich 2009 and Xapsos 2000 SEP spectra; Hassler et al. 2014 Mars dose; Nordheim 2018 and Cooper 2001 (Europa); AP9/AE9 belts; shower runs for H2, He, CH4, NH3 and CO2 columns; ICRP 103 numeric limits; W values and EUV outer-shell cross sections for Ar, Ne, Kr, Xe, HCN, H2CO, HCl, HF, CH3Cl, SiO, TiO, VO and the metals.
8. **Hot-giant and lava-world line lists and vapour pressures:** TiO, VO, FeH, CrH, Na, K line lists and pressure-broadened wings; partial pressures of Na, K, SiO, Mg, Fe, TiO, VO over magma (Schaefer, Lodders and Fegley 2012; LavAtmos arXiv:2210.10463; Miguel 2011 numeric tables).
9. **H2SO4-H2O and brines:** Tabazadeh 1997, Zeleznik 1991, Ayers 1980 for the acid vapour pressure by weight percent; perchlorate eutectics and deliquescence (Chevrier, Hanley, Gough); NH4Cl and HNO3 and HCl hydrates; clathrate dissociation (Sloan and Koh, Choukroun 2010).
10. **Comet per-species abundances and production-rate law:** Mumma and Charnley 2011, Bockelee-Morvan 2004, Dello Russo 2016, A'Hearn 2012 with usable OCR, Combi and Meech production rates; Harrington Pinto 2022.
11. **Photochemistry gaps:** JPL 19-5 evaluation (ClO + O, ClO + NO, ClO + ClO, Troe parameters for all falloff reactions), phosphorus chemistry rates and PH3 cross sections beyond Sousa-Silva 2020, F-star far-UV and extreme-UV spectra, Zerkle 2012 and Trainer 2006 originals (the haze threshold is read through Arney 2016), Harman 2015, Zahnle 2006, Hunten 1973.
12. **Stellar abundance trends:** Th, U, Eu and K against [Fe/H] and [alpha/Fe] for radioisotope scaling; measured K/U and Th/U for oceanic crust, Moon, Mars, Venus, Mercury and chondrite classes; sub-Neptune envelope metallicity occurrence distribution; Thorngren 2016 scatter of the mass-metallicity fit.
13. **Licences:** per-dataset licences of the Heays 2017 (Leiden), Huebner 2015 and Burkholder 2020 cross sections and the Serdyuchenko and Mainz files before any table built from them is distributed.
14. **Proxy widening still useful:** doi.org, science.org, sciencedirect.com, agupubs (Wiley), aanda.org, academic.oup.com, cambridge.org, pubs.acs.org, usgs.gov, mindat.org, lpi.usra.edu, github.com pages, jpldataeval.jpl.nasa.gov (Project settings, environment network). With those most of the "second-hand" and "not reached" rows above can be read from the papers.

## 7. Recalled values

Every RECALLED value is flagged in its row and listed in the "Recalled / to verify" section of its part (about 34 values, mostly band centres, lab constants and the Peng-Robinson constants, which matched NIST compressibility within a few percent). None of them is used as a verified input. Parts 01, 02, 03 and 06 hold nearly all of them.

## 8. Method and limits

Seven workers each took one domain under the same rules (scientifically accurate, nothing fixed or copied from Solar System presets, every number labelled, no unmarked recall). Solar System bodies appear only as rows marked TEST or CALIBRATION. The parts were written between 12:05 AM and 1:00 AM EDT. The workers' row counts are their own. Several widths in the generator recipes (for example the seeded spread on ozone, scatter on mixing ratios) are design choices and are flagged as such in the parts. This master and the parts say what the science supports; where the science is missing the part says so and offers a seeded physical range for the generator to use meanwhile, in line with the rule that nothing is fixed per class.

## 9. Boss-supplied research awaiting verification (2026-10-17, 1:22 AM EDT)

Boss posted a research text (tiers 1 to 5 and a fourteen-row priority matrix) in the body composition thread; the coordinator routed it to the Atmosphere research lane to verify and fold in. Until that is done it is not a verified input. It covers the cold sublimation curves (Fray and Schmitt 2009 with Grundy 2024 and Jager and Span 2012, the Feistel and Wagner water-ice formulation), magma-ocean redox (Sossi 2020, Gaillard 2014 and 2021, Hirschmann 2012), the visible absorption of CH4, NH3 and H2O for the giants, Mars dust and Titan haze optics, infrared optical depth and the Parmentier and Guillot profile, Rayleigh King factors, radiation belts and the ICRP weights, silicate vaporisation, brines and clathrates, comet abundances, photochemical rate gaps, stellar abundance systematics and radioisotopes, dataset licences and the network domains. Where it overlaps the lane's own reads the design compares them first (for example the sublimation validity ranges against the lane's finding that the Fray and Schmitt fits exceed the laboratory curves by up to 6 times below 40 K, and the Thorngren scatter of 0.3 to 0.5 dex against the 0.26 dex the Body composition lane read). The notes are on GEN.452, GEN.451, GEN.419, GEN.414, GEN.417, GEN.423, GEN.456, GEN.455, GEN.459, GEN.460, GEN.461, GEN.436 and GEN.410.

