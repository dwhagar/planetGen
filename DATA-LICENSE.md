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
| `cod-part-d` | Densities computed in part D (section 2.1) from the Crystallography Open Database cell of the same formula | COD entries are public domain (CC0); only the computed densities are kept |
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
| `ozone_uv` (AD9, GEN.443; MPI-Mainz UV/VIS Spectral Atlas snapshot of 2013-06-24, Zenodo record 6951) | CC BY 4.0; cite Keller-Rudek et al. 2013 (Earth Syst. Sci. Data 5, 365); the cross sections read are those of Daumont et al. 1992, Malicet et al. 1995 and Brion et al. 1998; the archive is downloaded at setup and removed after the four tables are kept |

## The Chappuis cross sections (`src/planetgen/physics/data/ozone_chappuis.json`, GEN.418)

The ozone cross section of the Chappuis band, 380 to 780 nm, as 5 nm cell means of the 295 K, 1 nm table of the Brion, Daumont, Malicet and Jenouvrier group (Daumont et al. 1992, Malicet et al. 1995, Brion et al. 1998) in the MPI-Mainz UV/VIS Spectral Atlas snapshot of 2013-06-24 (Zenodo record 6951, CC BY 4.0; cite Keller-Rudek et al. 2013, Earth Syst. Sci. Data 5, 365). Built by `scripts/build_ozone_chappuis.py` from the `ozone_uv` download; only these 81 derived values are in the repository.

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

## The Kepler recovery curve (`src/planetgen/physics/data/kepler_recovery.json`, GEN.504)

Fitted numbers only (the recovery curve, the effective-depth law, the noise lines, the quantiles of the Kepler magnitude and a table of observed against predicted recovery); the injection and stellar rows are not kept.

| Source | Terms known |
|---|---|
| NASA Exoplanet Archive, Kepler DR25 flux-level transit injections, group INJ1 (https://exoplanetarchive.ipac.caltech.edu/docs/KeplerSimulated.html), retrieved 2026-10-23 | Public data of the archive, operated by the California Institute of Technology under contract with NASA; the archive asks that its use be acknowledged, as done here |
| NASA Exoplanet Archive, table `keplerstellar` (delivery q1_q17_dr25_stellar), retrieved 2026-10-23 | Same |
| Burke and Catanzarite 2017, Kepler Science Document KSCI-19083 | Published document; the form of the recovery curve only |

## The thorium of thick-disc stars (`src/planetgen/physics/data/tuning_targets.json` targets o2-mean and o2-sd, GEN.514)

Two summary statistics only (the mean and the standard deviation of [Th/Mg] of 32 stars); the catalogue's rows are not kept.

| Source | Terms known |
|---|---|
| Hypatia Catalog (Hinkel et al. 2014, AJ 148, 54, arXiv:1405.6719), open web API, queried 2026-10 | The site states no open licence; the two statistics are computed from it and no rows are kept |

## The comets with a measured composition (`src/planetgen/physics/data/comet_admission.json`, GEN.507)

Fitted numbers only (counts of known and measured comets by group and perihelion bin, and two regression lines and their scatter); the database's rows are not kept.

| Source | Terms known |
|---|---|
| Robinson, Malamud, Opitom, Perets and Blum 2024, A&A (arXiv:2403.15644), database doi 10.7488/ds/7723 | CC BY 4.0; attribution given here |
| JPL Small-Body Database query API (https://ssd-api.jpl.nasa.gov/sbdb_query.api), retrieved 2026-10-23 | Public data of NASA/JPL; orbit class and perihelion distance only |

## The compaction table (`src/planetgen/physics/data/compaction_table.json`, GEN.442)

Computed here by `scripts/build_compaction_table.py` from the porosity laws and grain densities of Bierson and Nimmo 2019 (Icarus 326, 10, Table 1 and Eqs. 1 and 2; published paper, values only) and third-order Birch-Murnaghan constants that are the design note's stand-ins. No third-party data file is in the repository.

## The habitable-zone fits (`src/planetgen/physics/orbits.py`, GEN.146)

Eight coefficients of Kopparapu et al. 2014, ApJL 787, L29 (arXiv:1404.5292), Table 1, read from the paper (1 Earth-mass runaway greenhouse and
maximum greenhouse rows). Published paper; values only. The paper's raw text is not in the repository.

Ten more coefficients of the same table (GEN.639): the runaway greenhouse sets for 0.1 and 5 Earth masses (`KOPPARAPU_RUNAWAY_GREENHOUSE_BY_MASS`). Published paper; values only.

## The water layer of the structure solver (`src/planetgen/physics/bulk_structure.py`, GEN.429)

| Source | What was kept | Terms known |
|---|---|---|
| NIST Chemistry WebBook fluid properties for water (IAPWS-95), isotherm at 300 K from 0.1 to 950 MPa, read 2026-10-23 | Twenty density values | NIST Standard Reference Data are copyrighted by the Secretary of Commerce; no permission was sought; values only |
| Ice VI, VII and X reference points | Read from the materials table (`journaux2020`, `seager2007`) | As the materials table |

## Clathrate hydrates and eutectic brines (`src/planetgen/physics/data/clathrate_brine.json`, GEN.471)

The keys are the source keys of the file. Published values only; no third-party data file is in the repository.

| Source key | Source | Terms known |
|---|---|---|
| `marboeuf-2014` | Marboeuf et al. 2014, Table 5 (from Mohammadi and Richon 2010, Fray 2010, Hersant 2004) | Published paper; the six fit constants only |
| `semenov-2023` | Semenov, Stoporev et al. 2023, Int. J. Mol. Sci. 24, 9321 | Open access (MDPI); two quadruple points only |
| `chevrier-2009` | Chevrier, Hanley and Altheide 2009, GRL 36, L10202 | Published abstract; two eutectic points only |
| `pestova-2005` | Pestova et al. 2005, Marion et al. 2010, as quoted by arXiv:1811.00862 and arXiv:2012.00100 | Second-hand; one eutectic point |

## Papers, datasets and code read for the second-hand citation pass (10-24, `reports/secondhand-citations-2026-10-24.md`)

No code or data from these sources is in the repository. The pass produced a report; any value a lane later copies into a data
file takes its terms from the row here, and the lane adds its own section above. Raw papers and downloaded files stay out of the
repository. The "Basis" column says how the terms were learned: **file** (printed in the PDF Boss uploaded), **arXiv** (the licence
shown on the paper's arXiv abstract page), **repo** (the LICENSE file of the code repository, read 10-24), or **publisher** (the
publisher's usual terms; no statement was seen, so treat as unconfirmed).

### Uploaded by Boss

| Source | Terms | Basis |
|---|---|---|
| Curtis et al. 2020, ApJ 904, 140 | CC BY 4.0 | file |
| Anenburg 2026, "Oxygen Fugacity Buffers and an Interactive Calculator", G-cubed 10.1029/2026GC013032 | Open access, Creative Commons Attribution (version not printed) | file |
| Summa, Monaco, Lazzeretti and Zanasi 2021, J. Phys. Chem. Lett. 12, 8855 | (c) 2021 The Authors; ACS open licence, type not printed | file |
| Pahlevan et al., Science Advances (Redox state of Earth's magma ocean), 2020 | CC BY-NC (NonCommercial); only values are read | file |
| JPCA perspective on atmospheric photochemistry (2024) | CC BY 4.0 | file |
| Wolff et al. 2009, JGR 114, E00D04 | (c) AGU 2009; the Wiley page states open-access articles fall under their Creative Commons licence, type not printed | file |
| Lodders 2003, ApJ 591, 1220 | (c) American Astronomical Society, all rights reserved; values only | file |
| Palme and O'Neill 2003, Treatise on Geochemistry 2.01 | (c) Elsevier; values only | publisher |
| Pahlevan, Schaefer and Hirschmann 2019, magma ocean outgassing and hydrogen isotopes (EPSL) | (c) Elsevier; values only | publisher |
| Rudnick and Gao 2003, Treatise on Geochemistry 3.01 | (c) 2003 Elsevier Ltd; values only | file |
| Hirschmann 2021, Geochim. Cosmochim. Acta 313, 74 | (c) Elsevier; values only (the IW buffer is already used by the repository) | publisher |
| Bierson and Nimmo 2019, Icarus (density of Kuiper Belt objects) | (c) 2019 Elsevier, accepted manuscript | file |
| Gaillard and Scaillet 2014, EPSL 403, 307 | HAL open archive copy; Elsevier holds the copyright | file |
| Li, Cheng and Yang 2019, Geo-fO2, G-cubed 10.1029/2019GC008273 | (c) 2019 AGU; type of any open licence not printed | file |
| Stixrude and Lithgow-Bertelloni 2011, Geophys. J. Int. 184, 1180 | (c) Oxford University Press and the Royal Astronomical Society; values only | publisher |
| Zeleznik 1991, J. Phys. Chem. Ref. Data 20, 1157 | (c) 1991 US Secretary of Commerce, assigned to the American Institute of Physics and ACS; reprint read from the NIST reference data pages; values only | file |
| Tabazadeh et al. 1997, GRL 24, 1931 | (c) American Geophysical Union; values only | file |
| Myhre et al. 1998, GRL 25, 2715 | (c) American Geophysical Union; values only | file |
| Bucholtz 1995, Appl. Opt. 34, 2765 | (c) Optica (formerly OSA); values only | publisher |
| Bates 1984, Planet. Space Sci. 32, 785 | (c) Elsevier; values only | publisher |
| King 1923, Proc. R. Soc. A 104, 333 | Royal Society archive copy; a 1923 paper, no licence statement; values only | publisher |
| Sneep and Ubachs 2005, JQSRT 92, 293 | (c) Elsevier; values only | publisher |
| Irwin et al. 2018, detection of H2S above the clouds of Uranus | Author manuscript; Nature Astronomy terms apply | publisher |
| Hassler et al. 2014, Science 343, 1244797 (Mars RAD) | (c) AAAS; values only | publisher |
| Mumma and Charnley 2011, Annu. Rev. Astron. Astrophys. 49, 471 | (c) 2011 Annual Reviews | file |
| Paranicas, Cooper and Garrett 2009, Europa (Univ. of Arizona Press), 529 | (c) University of Arizona Press; no statement in the chapter text; values only | publisher |
| Gronoff et al. 2011, A&A 529, A143 | (c) ESO; no Creative Commons statement printed in the file | publisher |
| Poyatos et al. 2026 and Postolec et al. 2026, A&A manuscripts | (c) ESO 2026; author manuscripts, no Creative Commons statement printed | file |
| Gunther et al. 2020, TESS flares, AJ | DTU Orbit copy; copyright and moral rights remain with the authors or other copyright owners | file |
| McBride, Zehe and Gordon 2002, NASA Glenn coefficients (NASA/TP-2002-211556) | US Government work, no copyright asserted | publisher |
| ICRP Publication 103 free extract and preface | (c) ICRP; the full report was not read | file |

### Open papers read through arXiv, with the licence shown on the abstract page

| Licence | arXiv papers |
|---|---|
| CC BY 4.0 | Hammond et al. 2025 (2409.04386); He et al. 2022, Titan haze optical constants (2112.08521); Bower et al. 2022 (2110.08029); Drozdovskaya et al. 2023 (2307.11486); Beth et al., ROSINA ion zoo (2008.08430); Merg et al. 2025 (2505.09490); Blum and Bischoff (2207.12731) |
| CC BY-NC-SA 4.0 | Seignovert et al. 2017 (1704.00842); Maurice et al. 2023 (2301.07505) |
| CC0 1.0 | Chabrier 2021, H/He equation of state II (2107.04434) |
| Public domain dedication | Davenport 2014 (1410.7779) |
| arXiv non-exclusive licence (copyright stays with the authors; read and cite, no redistribution) | Vainio et al. 2017 (1707.00485); Koldobskiy et al. 2021 (2101.10234); Asvestari et al. 2017 (1612.02446); Wang, Lineweaver and Ireland 2018 (1708.08718); Lisse et al. 2021 (2009.02277); Mazevet et al. 2019 (1810.05658); Haldemann et al. 2020, AQUA (2009.10098); Gao et al. 2014 (1312.3750); the Mauna Kea extinction paper (1210.2619); Irwin et al. 2019 (1812.05383); Bailey et al. 2011 (1102.2554); Corlies et al. 2021 (2012.02247); Zahnle et al. 2020 (2001.00095); Ranjan and Sasselov 2017 (1610.06223); Heays, Bosman and van Dishoeck 2017 (1701.04459); Biver et al. 2012 (1201.4318); Chabrier, Mazevet and Soubiran 2019 (1902.01852); Schlichting, Sari and Yalinewich 2015 (1406.6435); Foley 2015 (1509.00427); Charnoz et al. 2009 (0912.3017); the stagnant-lid carbon cycle paper (1712.03614); Yoshizaki and McDonough (1909.07560) |
| arXiv assumed-licence (pre-2003 paper) | Crida, Morbidelli and Masset 2006 (astro-ph/0511082) |

Only formulas, numbers and tables of facts were read from these papers, each with its citation in the report.

### Datasets, tables and code

| Source | Terms | Basis |
|---|---|---|
| PetroBuffer (P. Liggins; Frost 1991 and Campbell 2009 buffer constants as coded) | MIT | repo |
| GEOS-Chem (JPL 19-5 and 15-10 rate lines) | MIT | repo |
| PICASO (King factors from Bogaard 1978) | GNU GPL v3 | repo |
| BurnMan (Stixrude and Lithgow-Bertelloni parameters) | GNU GPL (v2 or later per its package metadata) | publisher |
| VolFe (one benchmark value only) | Not confirmed; only a number quoted on its documentation page was used | n/a |
| Chabrier DirEOS2019 and DirEOS2021 tables (`perso.ens-lyon.fr/gilles.chabrier/DirEOS/`) | No licence in the README; the 2021 paper is CC0 | file |
| Leiden photodissociation and photoionisation database (home.strw.leidenuniv.nl/~ewine/photo/) | Free download; cite Heays et al. 2017 | publisher |
| IPCC Third Assessment Report, WG1 Chapter 6 (Myhre expressions, Table 6.2) | (c) IPCC; reproduction with attribution allowed | publisher |
| Wikipedia "Equivalent dose" and "Effective dose" | CC BY-SA 4.0 for the text; values are ICRP's | publisher |
| Canadian Nuclear Safety Commission radiation protection regulations schedule | Crown copyright, reproduction allowed with attribution | publisher |
| Putirka and Xu 2021, Nature Communications, supplement and peer review file | CC BY 4.0 (Nature Communications open access) | publisher |
| Vainio 2017 event table typed into `sep-band-parameters-vainio2017.csv` (report folder only) | Facts from the paper; stays out of the repository | n/a |

### Not read at all

Fray and Schmitt 2009, ICRP 103 (full report), Karkoschka 1994, Karkoschka and Tomasko 2009/2010, Tomasko et al. 2008, Ayers,
Gillett and Gras 1980, Frost 1991, Hirschmann 2012, Palme and O'Neill 2014, JPL Publication 19-5, Bogaard et al. 1978, Huebner et
al. 1992, Huebner and Mukherjee 2015, Cooper et al. 2001, Patterson et al. 2012, Nordheim et al. 2018, Gronoff et al. 2009,
Dartnell. Where a lane uses a number carried by another work, the report cites that work and its tier. Their terms are the
publishers' and are not recorded here because nothing was taken from the originals.

## Water against formation distance (`docs/design/composition/volatile-delivery-formation-distance.md`, GEN.481)

Derived counts and logistic fits from the Bern NGPPS NG76 population: Burn et al. 2024 (Nature Astronomy 8, 463), data and code of Zenodo record 10527170 (no licence field), and the starting distances of the same planets from the public DACE population ng76. Emsenhuber et al. 2021 (A&A 656, A69; arXiv:2007.05561) was read for the disc model. The raw files are never in the repository and are downloaded only by the reproduction script kept in the project files.
