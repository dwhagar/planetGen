# Data catalog: what data and tables the project already holds

Written by the TODO thread, 2026-10-29 (Boss, 8:51 PM EDT: every TODO item and every thread must know what data we already have, in case it helps what they are doing). Built from the data-use audit (`/mnt/project-files/reports/data-use-audit-2026-10-29.md`, which holds the measurements and the import graph), `/mnt/project-files/data`, the setup fetch registry (`src/planetgen/datasets/catalogue.py`), `src/planetgen/physics/data/`, the `*-files` folders under `/mnt/project-files/reports/`, and the tables the code carries in Python.

**Required reading for every lane before hand-tuning a number or searching for a source.** If a dataset below could serve your item, use it first (Boss's source rules of 8:16 PM EDT 2026-10-29: physics wins, secondary sources are accepted, a registered stand-in is fine). If you use a held dataset that a fresh server does not have, add its setup catalogue entry in the same change (raw datasets never go in the repo; the licence line in `DATA-LICENSE.md` is Boss's).

Where things live: **S** = downloaded by `update.sh` into `/var/cache/planetgen/datasets` (exists on every fresh server); **P** = held only in `/mnt/project-files` (NOT on a fresh server: a code reader needs a setup entry first); **R** = a table in the repo (always there).

## 1. Ranked by usefulness (how many wrong or hand-tuned things it can fix, then how ready it is)

| Rank | Data | Where | Status | Serves |
|---|---|---|---|---|
| 1 | ExoMol correlated-k tables (H2O, CO2, CH4) and HITRAN collision-induced absorption | S, 1,453 MB downloaded, 6 MB kept | built, merged, **read by nothing** | GEN.850, GEN.851, GEN.852, GEN.972 to GEN.975, GEN.752 |
| 2 | MIST v1.2 isochrones, all 25 columns and every phase (setup keeps 6 columns of the main sequence only) | S, 221.5 MB | post-main-sequence and pre-main-sequence rows are read and discarded at setup | GEN.935, GEN.940, GEN.947, GEN.939, GEN.941, GEN.943, GEN.744, GEN.922 and GEN.927 to GEN.931 (checks), GEN.942 |
| 3 | [R0254] (BHAC15) isochrone table | S, setup entry and reader `planetgen.datasets.baraffe_isochrones` (#3808) | reader built, not yet read by the generator; three constants still fitted by hand | GEN.936, GEN.948, GEN.947 |
| 4 | Exoplanet occurrence and mass-radius data (Kepler DR25 and Exoplanet Archive fits, the 2,000-system class distribution check) | R (4 harness JSONs), P (`reports/data/planet-class-distribution`) | harness only | GEN.728, GEN.729, GEN.730, GEN.722, GEN.828, GEN.829, GEN.920, TEST.395, GEN.753 |
| 5 | NOAA SWPC proton events above 10 MeV 1976 to 2017 (265 events by Sources 1's count) and 47 Tylka and Dietrich events, SEP band parameters [R4478] | S, setup entry and reader `planetgen.datasets.noaa_proton_events` (#3808); the SEP band file stays in `reports/secondhand-citations-2026-10-24-files/` | reader built, not yet read by the generator (GEN.747) | GEN.747, GEN.746, GEN.749 |
| 6 | Neutron star matter: 21 equations of state [R2834] and outer-crust compositions [R4028] | S, 2.0 MB | read by tests only | GEN.743, GEN.968, GEN.969, GEN.937 |
| 7 | [R1827] (2010 edition) globular clusters, 157 clusters | S, setup entry and reader `planetgen.datasets.harris_clusters` (#3808) | reader built, not yet read by the generator (GEN.297) | GEN.297, GEN.209, GEN.531 family, cluster laws in `tuning.py` |
| 8 | AD1 methane and ammonia tables: [R2265] k-tables, PDS gbat_0001, NMSU planetary datasets | P, `data/ad1-methane/` (about 19 MB) | not read; licence entry is Boss's | GEN.949, GEN.594 |
| 9 | [R0535] (1,615 stars, 15 elements) | P, `reports/gap-filling-plan-2026-10-17-files/brewer2016_vizier_J_ApJS_225_32.tsv` | fitted to trends and a 4 by 4 matrix only | TEST.456, GEN.465 family, star abundance items |
| 10 | [R2726] condensation tables (Table 11), [R4658] chondrite groups, UNSCEAR 2008 extracts | P, `design-drafts/gen465-second-pass/` | read once by a research pass | GEN.950, GEN.468, GEN.473 to GEN.477 |
| 11 | MIST phase durations at 4 [Fe/H] (derived) | P, `reports/best-guess-standins-2026-10-24-files/mist_post_main_sequence_durations.csv` | not read | GEN.935, GEN.940 |
| 12 | Second-hand citation files: [R1464] buffers, JPL 19-5 rate lines, Leiden solar photo rates | P, `reports/secondhand-citations-2026-10-24-files/` | photolysis lifetimes JSON exists, its module has no caller | GEN.794, GEN.748, GEN.746 to GEN.749, GEN.854, atmosphere chemistry items |
| 13 | Courteau, Comeron, Pizzella and VizieR galaxy tables (c18 tableb and tablec, p05 tables, digitised figures of [R0338] occupation, the [R3209] cluster table, [R1169]) | P, `reports/vizier-courteau-comeron-2026-10-17-files/` | used for the galaxy-laws research | GEN.612 to GEN.618, GEN.531 family, GEN.380 to GEN.387 (galaxy types), GEN.535, GEN.536 |
| 14 | Bern population synthesis ([R4529], Zenodo 10527170, 31,868 planets), as a binned fit | R (`ngpps_formation_water.json`) | used well, as a test target | GEN.849, GEN.510, water delivery rows b2 and b3 |
| 15 | [R0308] white dwarf cooling sequences | S, 2.5 MB | used (cooling, radius); 9 of 16 numbers | GEN.742 (done), GEN.743 (helium-core and ONe layers) |
| 16 | Ozone cross sections, MPI-Mainz atlas | S, 126.3 MB; 4 of 5,756 files kept | used | GEN.748, GEN.746 |
| 17 | [R1429, R2741, R4881] tables | P (`reports/bulk-composition-sources-check-2026-10-17/`), R (`giant_radius_tables.py`) | Fortney transcribed and used; the others as laws or checks | GEN.733, GEN.964, GEN.965, GEN.809 to GEN.811, GEN.816, GEN.784 |

## 2. Setup datasets (S): the six `update.sh` fetches

All are in `src/planetgen/datasets/catalogue.py`, fetched by `planetgen datasets fetch`, with the exact size and sha256 in the entry. A "transient" piece is deleted after setup reduces it. Readers are in `src/planetgen/datasets/`.

| Dataset | Downloaded | Kept | Reader module | Read by (generator) |
|---|---|---|---|---|
| MIST v1.2 isochrones [R0773, R1150] | 221.5 MB | `main_sequence.npz` about 4 MB, 340,000 rows, 6 of 25 columns, phase 0 only | `datasets/isochrones.py` | `physics/main_sequence.py`, `stellar_composition.py` |
| [R0308] white dwarf sequences | 2.5 MB | the archive, 46 files | `datasets/white_dwarf_models.py` | `stellar_evolution.py`, `stellar_composition.py` |
| Neutron star matter (21 EOS, 20 crust models) | 2.0 MB | the files | `datasets/neutron_star_matter.py` | crust only, on demand; EOS by tests |
| Ozone cross sections (MPI-Mainz) | 126.3 MB | `ozone_uv.npz` | `datasets/ozone_uv.py` | `uv_cross_sections.py`, `body_radiation.py`, `atmospheric_radiation.py` |
| ExoMolOP correlated-k tables | 1,257.8 MB | `ir_ktables.npz` about 6 MB | `datasets/ir_ktables.py` | `line_opacity.py`, only by `radiative_convective.py`, which nothing calls |
| HITRAN CIA, 11 pair files | 195.4 MB | `ir_cia.npz` | `datasets/hitran_cia.py` | same chain, not wired |

Setup entries built by Sources 1 (#3808, 9:18 PM EDT 10-29, as the lane reported): `planetgen.datasets.harris_clusters` (GEN.297), `noaa_proton_events` (GEN.747), `baraffe_isochrones` (GEN.936 and GEN.948); registry rows SC4 to SC6. Section 3 rows for these three are now S, not P.

Setup entry built by Star 1 (GEN.943): `planetgen.datasets.mamajek_sequence`, the Pecaut and Mamajek mean dwarf spectral-type sequence (version 2022.04.16, 56 kB, 67 whole-number types O3V to M9V), read by `stellar_evolution.spectral_letter_and_subclass`; registry row PS12.

Wiring rule (Sources 1, GEN.754): a new lock-free module, dataset reader, data file or registry row that nothing reads yet needs one line in `docs/design/research/wiring-register.md` naming its open wiring item (`python scripts/wiring_guard.py --check`); a new L3 or L4 registry row needs a line in `docs/design/research/burn-down-register.md`.

## 3. Held in the project files (P)

| Data | Path under `/mnt/project-files` | What it is | Could feed |
|---|---|---|---|
| [R1827] | `data/harris/mwgc.dat` | 157 Milky Way globular clusters, about 50 columns | cluster mass, radius and metallicity laws (GEN.297, GEN.209) |
| [R0254] | `data/gen549-baraffe2015/BHAC15_iso.2mass` | pre-main-sequence and low-mass isochrones by age and mass, down into the brown dwarfs | GEN.936, GEN.948, GEN.947 |
| AD1 methane | `data/ad1-methane/` (README, `oxford-karkoschka-ktables`, `pds-gbat_0001`, `nmsu-planetary-datasets`) | 3,689-wavelength methane k-tables at 100, 198, 296 K; giant-planet spectra and absorption tables | giant visible colour (GEN.949), the AD1 fallback (GEN.594) |
| NOAA proton events | `reports/best-guess-standins-2026-10-24-files/noaa_swpc_proton_events_gt10MeV_1976_2017.csv` | 255 events | SEP rates and Weibull parameters (GEN.747) |
| MIST phase durations | `reports/best-guess-standins-2026-10-24-files/mist_post_main_sequence_durations.csv` | phase lengths at 4 [Fe/H] | GEN.935, GEN.940 |
| [R1464], JPL 19-5, [R4478] | `reports/secondhand-citations-2026-10-24-files/` | oxygen buffers, photo-reaction rate lines, solar photo rates at 1 AU, SEP band parameters | photolysis, redox, SEP (GEN.794, GEN.748, GEN.747) |
| [R0535] | `reports/gap-filling-plan-2026-10-17-files/brewer2016_vizier_J_ApJS_225_32.tsv` | 1,615 stars, 15 element abundances | star composition draw and TEST.456 |
| Gap plan scripts | `reports/gap-filling-plan-2026-10-17-files/` (`acc.py`, `arx.py`, `water_share.py`, `ws_stars.py`, `brewer_fit.py`) | the fills tested in the 10-17 plan | GEN.950, GEN.468 and the registry |
| VizieR galaxy tables | `reports/vizier-courteau-comeron-2026-10-17-files/` | Courteau and Comeron tables (c18), [R3540] (p05), digitised occupation and size figures (pms01: blf20 occupation, [R1169, R3209] cluster table) | galaxy laws, halo and dwarf laws |
| Gemini batch check | `reports/gemini-batch1-check-2026-10-17-files/` | the verification outputs for the first unverified batch | tier checks only (LEAD material) |
| Bern analysis scripts | `reports/gen481-files/`, `reports/gen802-files/` | analysis and fit of the formation-water statistics | GEN.849, GEN.510 |
| Planet class distribution check | `reports/data/planet-class-distribution/` (`classes.py`, JSON, summary) | 2,000 generated systems by class and mass bin | GEN.828, GEN.829, TEST.395 |
| Life chemistry review | `reports/life-chemistry-review-2026-10-28-files/` | population analyses behind the life classes | GEN.862 to GEN.894, TEST.435, TEST.436 |
| Eccentricity-first probes | `reports/gen842-eccentricity-first-2026-10-29-files/`, `reports/gen842/`, `reports/gen842-tides-2026-10-29/` | prototype patch and measurement scripts for the giant radius, spin and tide censuses | GEN.961, GEN.976, GEN.953 to GEN.960, GEN.978, GEN.979 |
| Boss's papers | `reports/boss-papers-2026-10-17-files/` and `docs/design/research/boss-papers-2026-10-17.md` | optics, sulfuric acid forcing, radiation and comets, mantle and crust, galaxy disc and halo | atmosphere, radiation, comet and galaxy items |
| Papers Boss cannot access | list in `reports/secondhand-citations-2026-10-24.md` | second-hand values for Fray and Schmitt, ICRP 103, [R4392, R4394], [R0210, R1464, R1953, R3417], JPL 19-5, [R0443], Huebner | any item that would cite them: use the second-hand value (a secondary source is accepted) |

## 4. Tables in the repo (R)

`src/planetgen/physics/data/*.json` (30 files, all small). Read by the live generator: `age_needs`, `atmosphere_recipes` (101 KB), `brown_dwarf_teff_map`, `clathrate_brine`, `comet_population`, `compaction_table`, `decay_series`, `elements`, `gases` (232 KB), `isotopes` (129 KB), `life_classes`, `life_class_text`, `life_ledger`, `life_supply`, `materials` (160 KB), `nasa_glenn_species`, `ozone_chappuis`, `planet_classes`, `proton_range`, `saturation_pressure`, `tuning_targets`, `white_dwarf_mass_radius`. Read only by the tuning harness (by design): `comet_admission`, `jwst_selection`, `kepler_recovery`, `mass_radius_sample`, `wd_pollution_sample`. Read by modules nothing calls: `design_flags` (`mixture_trace.py`, GEN.479), `photolysis_lifetimes` (`photolysis.py`). Read by tests only, as a validation target: `ngpps_formation_water`.

Tables in Python: `physics/giant_radius_tables.py` ([R1429] radii, all three tables), `tuning.py` (4,275 lines: the anchors, ranges and laws, most with their source in the docstring), `tests/data/prevalence_expectations.json` (the researched prevalence ranges the census checks), the ozone Chappuis table in `ozone_chappuis.json`, the unit registry (`GEN.568`).

The **stand-in registry** (`docs/design/research/stand-in-values.md`, 316 rows, GEN.725) lists every registered stand-in with its source, tier, tweak level and the item that replaces it: 216 at L0, 18 at L1, 10 at L2, 55 at L3, 17 at L4. Check it before adding a constant: the number may already be there with a source.

Design research with data in its text: `docs/design/research/gap-filling-plan.md` (every filled gap with its source and fill), `docs/design/research/open-gaps-closest-laws-2026-10-17.md`, `docs/design/research/ad1-methane-sources-check-2026-10-17.md`, `docs/design/research/imported/` (Boss's imported research notes on habitability, radiation, toxicity, orbital algorithms).

## 5. What we do NOT hold (fetch or second-hand needed)

Occurrence data as a setup dataset (the Exoplanet Archive tables are fitted into JSONs, not held), hydrogen-helium equation of state tables (GEN.964), BT-Settl and MUSCLES spectra (GEN.970), AP-8 and AE-8 radiation belt models, measured rotation periods of any rocky exoplanet (none exist), radiation shower simulations for GEN.746 (second-hand only). For the Solar System, statistics are allowed as a check but nothing fixed is copied (NO FIXED CONTENT).

## 6. Items that need no data

UX, MAP, NAV, API, DB, OPS, SEC, PERF, ADM, POP, VIEW, USR and DOC items (interface, storage, speed, security, docs) do not use these datasets and carry no "Data we hold" note. Each open GEN and TEST item that a dataset above could serve carries a one-line "Data we hold" note in docs/TODO.md.
