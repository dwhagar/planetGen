# Data sources and their terms

The program's own code and text are released under CC0 (LICENSE.md). This file records where the data files in the repository came
from and what is known of their terms. The data tables carry no licence column (Boss, 2026-10-22: "Remove a license column and
move license information to its own file in root"; the file is this one, DATA-LICENSE.md). Boss manages licensing.

Rule (Boss, 2026-10-17, 2:35 PM EDT): every paper and dataset acquired is usable. The raw papers and downloaded datasets are
never in the repository and are never redistributed or sold; the values read from them and derived from them are published, and
the program stays free forever. Where no licence statement was found, the entry says so.

## The materials table (`src/planetgen/physics/data/materials.json`, GEN.427)

The keys are the source keys of the file.

| Source key | Source | Terms known |
|---|---|---|
| `research-e` | Part E of the body composition research (this repository, `docs/design/composition/body-composition-parts/`) | The project's own text; its values are read from the sources below |
| `nist-webbook` | NIST Chemistry WebBook phase-change and critical data | NIST Standard Reference Data are copyrighted by the Secretary of Commerce; no permission was sought |
| `crc-chemicals` | CRC Handbook tables as carried in the `chemicals` 1.5.2 package | Copyright of the handbook's publisher; the package's own licence is open source; read second-hand |
| `pdg-nuclear` | Particle Data Group atomic and nuclear properties | PDG data are free to use with citation |
| `slb2011-burnman` | Mineral parameters of Stixrude and Lithgow-Bertelloni (2011), as shipped in `burnman` | `burnman` is open source (GPL); read second-hand, only values are kept |
| `slb2024-burnman` | Stixrude and Lithgow-Bertelloni (2024), as parameters read in `burnman` 2.1.0 | As above |
| `e-table-2.6` | Band reflectance computed in part E from spectral library values | Computed values; the libraries (USGS splib07a and others) are public or open data |
| `computed-e` | Values computed in part E from its other columns | None |
| `handbook-mineralogy` | Handbook of Mineralogy, Mineral Data Publishing and the Mineralogical Society of America (2001-2005) | "All rights reserved" is printed on each page; only the density values read from the pages are kept, with the citation |
| `webmineral` | Webmineral Mineralogy Database (David Barthelmy) | No licence statement found; only the listed density values are kept |
| `crc-density-wiki` | Wikipedia, "Densities of the elements (data page)" | CC BY-SA 4.0 for the page text; the values are CRC and WebElements figures; only the values are kept |
| `pubchem` | PubChem compound pages, section Density | PubChem data are open; each annotation quotes its own source (CRC, USCG, NTP, EPA, HSDB) |
| `seager2007` | Seager, Kuchner, Hier-Majumder and Militzer 2007, ApJ 669, 1279 (arXiv:0707.2895) | Published paper; Tables 1 to 3 were read and their values are kept with the citation |
| `hakim2018` | Hakim et al. 2018, Icarus 313, 61, as read in part B | Published paper; values only |
| `part-b-fit` | Part B's own fit (this repository) | The project's own work |
| `exoplex` | ExoPlex 1.0.1a1 source distribution | Open source; read second-hand |
| `hp2011-burnman` | Holland and Powell (2011) data set as shipped in `burnman` | As for `burnman` |
| `journaux2020` | Journaux et al. 2020 Gibbs splines as distributed in SeaFreeze 1.1.3, evaluated in part C | SeaFreeze is GPL-3; only the spline data were used and the evaluator is the project's own; values only |
| `blum-bischoff` | Blum and Bischoff, arXiv:2207.12731, as read in part C | Preprint; values only |
| `ideal-gas` | The ideal gas law | None (exact law) |

## The gas table (`src/planetgen/physics/data/gases.json`, GEN.410)

Every row was recorded as `open` (no licence found). NIST Standard Reference Data are copyrighted (permission is noted in body
composition part G4); PubChem and USGS data are public domain. The sources are listed in the file's `sources`.

## The proton range file (`src/planetgen/physics/data/proton_range.json`, GEN.453)

The CSDA range and total stopping power of protons at 26 energies from 1 MeV to 10 GeV in 17 materials, as printed by NIST
PSTAR (NIST Standard Reference Database 124, ICRU Report 49), thinned from PSTAR's 133 energies by
`scripts/build_proton_range.py`. NIST Standard Reference Data are copyrighted by the Secretary of Commerce; no permission was
sought. The raw pages are never in the repository.

## The element and isotope files (`elements.json`, `isotopes.json`, `decay_series.json`, GEN.438)

- `elements.json`: the tables are CIAAW's (c) 2007-2024; no licence statement was found.
- `isotopes.json`: NUBASE2020 and AME2020 are the AMDC's, NIST's listing carries no licence statement, and the two papers (Lugaro et
  al. 2018, Dye 2012) are published.
- `decay_series.json`: computed from the NUBASE2020 file.

## Datasets fetched at setup (`src/planetgen/datasets/`)

| Dataset | Terms recorded |
|---|---|
| `isochrones` (MIST v1.2) | Open; cite Choi et al. 2016 and Dotter 2016 |
| `neutron_star_matter` | CC BY 4.0 |
| `white_dwarf_models` | Open; cite Bedard et al. 2020 |
| `white_dwarf_mass_radius.json` | Computed here from physical constants; no third-party data |

## The comet volatile population (`src/planetgen/physics/data/comet_population.json`, GEN.495)

The keys are the source keys of the file. Only values computed or read from the papers are kept; the database is not in the repository.

| Source key | Source | Terms known |
|---|---|---|
| `rob24` | Robinson, Malamud, Opitom, Perets and Blum 2024 (arXiv:2403.15644) and its database (Edinburgh DataShare, doi 10.7488/ds/7723) | The database was downloaded for the research; per-species statistics were computed and only those are kept; no licence statement was read |
| `hp22` | Harrington Pinto et al. 2022 (arXiv:2209.09985) | Published paper; values only |
| `coch15` | Cochran et al. 2015, Space Science Reviews 197, 9 (arXiv:1507.00761) | Published paper; values only |
| `cor25` | Cordiner et al. 2025 (arXiv:2508.18209) | Preprint; two measured ratios only |

## The polluted white dwarfs (`src/planetgen/physics/data/wd_pollution_sample.json`, GEN.508)

The keys are the source keys of the file. Values read from the tables of the papers (the median oxygen fugacity with its bounds, the type and the effective temperature of each star); the papers' text is not in the repository.

| Source key | Source | Terms known |
|---|---|---|
| `doyle20` | Doyle, Klein, Schlichting and Young 2020, ApJ (arXiv:2008.05992), Table 1 | Published paper; values only |
| `doyle23` | Doyle, Klein, Dufour et al. 2023 (arXiv:2303.00063), Tables 2 and 4 | Published paper; values only |

## The precisely measured planets (`src/planetgen/physics/data/mass_radius_sample.json`, GEN.505)

Fitted numbers only (two regression lines and their scatter per sub-sample, and counts); the archive's rows are not kept.

| Source | Terms known |
|---|---|
| NASA Exoplanet Archive, table `pscomppars` (https://exoplanetarchive.ipac.caltech.edu), retrieved 2026-10-23 | Public data of the archive, operated by the California Institute of Technology under contract with NASA; the archive asks that its use be acknowledged, as done here |
| Luque and Palle 2022, Science 377, 1211 | Published paper; the two cuts (mass better than 25 percent, radius better than 8 percent) only |

## The planets with a JWST spectrum (`src/planetgen/physics/data/jwst_selection.json`, GEN.506)

Fitted numbers only (logistic coefficients, the host magnitude model and per-class counts); the archive's rows are not kept.

| Source | Terms known |
|---|---|
| NASA Exoplanet Archive, tables `pscomppars` and `spectra` (https://exoplanetarchive.ipac.caltech.edu), retrieved 2026-10-23 | Public data of the archive, operated by the California Institute of Technology under contract with NASA; the archive asks that its use be acknowledged, as done here |
| Kempton et al. 2018, PASP 130, 114401 (arXiv:1805.03671) | Published paper; the two formulas and the scale factors of its Table 1 only |

## The comets with a measured composition (`src/planetgen/physics/data/comet_admission.json`, GEN.507)

Fitted numbers only (counts of known and measured comets by group and perihelion bin, and two regression lines and their scatter); the database's rows are not kept.

| Source | Terms known |
|---|---|
| Robinson, Malamud, Opitom, Perets and Blum 2024, A&A (arXiv:2403.15644), database doi 10.7488/ds/7723 | CC BY 4.0; attribution given here |
| JPL Small-Body Database query API (https://ssd-api.jpl.nasa.gov/sbdb_query.api), retrieved 2026-10-23 | Public data of NASA/JPL; orbit class and perihelion distance only |

## The habitable-zone fits (`src/planetgen/physics/orbits.py`, GEN.146)

Eight coefficients of Kopparapu et al. 2014, ApJL 787, L29 (arXiv:1404.5292), Table 1, read from the paper (1 Earth-mass runaway greenhouse and
maximum greenhouse rows). Published paper; values only. The paper's raw text is not in the repository.
