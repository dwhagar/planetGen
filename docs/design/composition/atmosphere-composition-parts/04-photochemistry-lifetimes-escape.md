# Part 04: Photochemistry, ozone, lifetimes, escape, exospheres and secondary effects

Read-only research note for the planetGen atmosphere effort. Nothing in the repository was edited. No Solar System body is used as a fixed input; Solar System bodies appear only as rows marked TEST or CALIBRATION.

Status labels: READ (original paper text), 2ND-HAND (read in a paper or code that cites the original), DATA (data file in a package or data repository), COMPUTED (re-derived by me with Python from stated inputs), SNIPPET, RECALLED (memory, unverified, also listed in section 9), NOT REACHED.

## 0. Provenance and network report

Domains that worked in this session: arxiv.org (PDF, abs, search HTML), zenodo.org (API and files), raw.githubusercontent.com, archive.stsci.edu (MAST MUSCLES spectra), lasp.colorado.edu/lisird (solar reference spectra), acp.copernicus.org (PDF), link.springer.com (open-access PDFs), uv-vis-spectral-atlas-mainz.org (host without "www"), nicholaswogan.github.io (Photochem docs), kinetics.nist.gov (front page only; no reaction query was made), iupac-aeris.ipsl.fr and cdsarc/vizier (front pages only).

Domains that failed: jpldataeval.jpl.nasa.gov (JPL 19-5 evaluation, unreachable), diva-portal (connection reset), nrel.gov and rredc.nrel.gov, www2.iap.tuwien.ac.at, one NTRS record (404), api.github.com (403), PyPI `photochem` (no such package), export.arxiv.org API (rate limited), the "www" Mainz host (reset).

Primary data sources used:

| Item | Where | Licence note | Status |
|---|---|---|---|
| Reaction mechanism `zahnle_earth.yaml` (706 reaction entries) | Zenodo `photochem_clima_data` v0.3.0, path `photochem_clima_data/data/reaction_mechanisms/`; same file on raw GitHub `Nicholaswogan/photochem_clima_data/main/...` | Zenodo metadata says gpl-3.0 | DATA |
| Photolysis cross sections and per-branch quantum yields `xsections/*.h5` | same dataset; its README says data came from Kevin Zahnle and Phidrates | gpl-3.0 (Zenodo metadata) | DATA |
| Solar spectrum (solar-minimum reference week) | LISIRD `whi_ref_spectra.csv`, column irradiance_apr10_apr16, Woods et al. 2009 (Whole Heliosphere Interval) | not checked | DATA |
| Stellar spectra eps Eri (K2V), GJ 832 (M1.5), GJ 176 (M2.5), GJ 876 (M4), GJ 699 (M4), TRAPPIST-1 (M8) | MAST MUSCLES adapt-const-res-sed FITS (v22 for the first four, v25 for GJ 699 and TRAPPIST-1) | not checked | DATA |
| Photochem code description: 98 neutral gas species, 610 forward reactions, 95 photolysis reactions, H/O/C/N/S/Cl/He only | Wogan 2025, arXiv:2509.25578 | n/a | READ |

The shipped mechanism has no phosphorus, so there is no PH3 network in it. It also lacks ClO+O -> Cl+O2, ClO+NO and ClO+ClO, so the chlorine catalytic cycle is incomplete (see section 2.2).

## 1. Scope and how the numbers feed a generator

Inputs from earlier generator stages: host spectral class and age (UV/XUV history), orbital distance, planet mass and radius, surface pressure and temperature, bulk composition (N2, CO2, H2, O2, H2O ...), volcanic and biological source fluxes, magnetic field flag.

Outputs of this part:

1. Photochemical lifetime of every trace species against UV, scaled by host star (sections 2.3 and 2.4). The generator multiplies the table value by the orbital flux scaling S and a geometry factor, then compares it with sink and source terms to decide which species can accumulate.
2. A rule set for the CO2/CO/O2 triple, methane and haze, sulfur and nitrogen families (section 2.6).
3. An ozone column, its host-star and O2 dependence, and the false-positive channels for O2/O3 (section 3).
4. Escape rates and which gases are lost over a given age, using XUV evolution laws and a choice of escape formula (section 4).
5. Exosphere density and source/sink scalings for airless bodies (section 5).
6. Flare, particle and lightning perturbations (section 6).

The cosmic shoreline, Jeans lambda thresholds with the project convention, the exobase temperature model and the `retention()` code already live in `/home/user/planetGen/docs/design/atmospheres-retention-and-classes.md`. They are cited here and not redone. The Kozakis 2022 ozone coefficients and the Cooke, Segura and Rugheimer cross-checks live in `/home/user/planetGen/docs/design/atmosphere-sources-check.md` and are also cited, not redone.

## 2. Photochemistry

### 2.1 Starter network (rate formats and values)

Format used by the Photochem mechanism file (DATA, checked by COMPUTED reproduction of one value): k = A T^b exp(-Ea/T) with T in K and Ea in K, in cm3 molecule-1 s-1 for bimolecular reactions and cm6 molecule-2 s-1 for low-pressure limits. Falloff reactions give only a low-pressure pair (k0) and a high-pressure pair (kinf) and no Troe keys, so the form is the Lindemann expression k = k0 [M] / (1 + k0 [M] / kinf). That reading of the file is an assumption (RECALLED convention, section 9). Check: O + O2 + M with k0 = 2.989e-28 T^-2.3 gives 6.0e-34 cm6 s-1 at 298 K, which matches the commonly tabulated JPL value (the JPL value itself was not reached; COMPUTED).

Column "k" gives COMPUTED values at the two temperatures from the file coefficients. For falloff rows, M = 2.5e19 cm-3 at 298 K and M = 1e18 cm-3 at 220 K (an upper-troposphere-to-stratosphere illustration only).

Chapman and oxygen family:

| Reaction | A, b, Ea (K) | k at 220 K | k at 298 K | File ref | Status |
|---|---|---|---|---|---|
| O + O2 + M -> O3 + M | k0: 2.989e-28, -2.3, 0; kinf: 2.8e-11, 0, 0 | 1.2e-15 (M=1e18) | 1.5e-14 (M=2.5e19) | De97 | DATA, COMPUTED |
| O + O + M -> O2 + M | k0: 5.2e-35, 0, -900; kinf 1e-12 | n/a | n/a | Ts86 | DATA |
| O + O3 -> 2 O2 | 8.0e-12, 0, 2060 | 6.9e-16 | 8.0e-15 | At04 | DATA, COMPUTED |
| O1D + N2 -> O + N2 | 1.8e-11, 0, -110 | 3.0e-11 | 2.6e-11 | De97, At97 | DATA, COMPUTED |
| O1D + O2 -> O + O2 | 3.2e-11, 0, -70 | 4.4e-11 | 4.0e-11 | At07 | DATA, COMPUTED |
| O1D + CO2 -> O + CO2 | 7.4e-11, 0, -120 | 1.3e-10 | 1.1e-10 | De97 | DATA, COMPUTED |
| O1D + H2O -> 2 OH | 2.2e-10, 0, 0 | 2.2e-10 | 2.2e-10 | Sa03 | DATA |
| O1D + CH4 -> CH3 + OH | present in file (second channel H2COH + H) | not extracted | not extracted | not extracted | DATA |

HOx family:

| Reaction | A, b, Ea (K) | k at 220 K | k at 298 K | File ref | Status |
|---|---|---|---|---|---|
| OH + O3 -> HO2 + O2 | 1.6e-12, 0, 940 | 2.2e-14 | 6.8e-14 | De97 | DATA, COMPUTED |
| H + O3 -> OH + O2 | 1.4e-10, 0, 470 | 1.7e-11 | 2.9e-11 | De97 | DATA, COMPUTED |
| OH + HO2 -> H2O + O2 | 4.8e-11, 0, -250 | 1.5e-10 | 1.1e-10 | At04 | DATA, COMPUTED |
| HO2 + O -> OH + O2 | 2.7e-11, 0, -240 | 8.0e-11 | 6.0e-11 | At04 | DATA, COMPUTED |
| HO2 + HO2 -> H2O2 + O2 | 3.19e-19, 2.01, -1274.3 | 5.3e-12 | 2.2e-12 | Ki84, Kl22 | DATA, COMPUTED |
| OH + OH + M -> H2O2 + M | k0: 6.615e-29, -0.8, 0; kinf 2.6e-11 | n/a | n/a | At97 | DATA |
| H + O2 + M -> HO2 + M | k0: 1.582e-27, -1.8, 0; kinf 7.5e-11 | 9.6e-14 (M=1e18) | 1.4e-12 (M=2.5e19) | At97 | DATA, COMPUTED |
| HO2 + O3 | not found as a separate entry in the file | n/a | n/a | n/a | NOT REACHED |

NOx family:

| Reaction | A, b, Ea (K) | k at 220 K | k at 298 K | File ref | Status |
|---|---|---|---|---|---|
| NO + O3 -> NO2 + O2 | 1.4e-12, 0, 1310 | 3.6e-15 | 1.7e-14 | At04 | DATA, COMPUTED |
| NO2 + O -> NO + O2 | 5.5e-12, 0, -188 | 1.3e-11 | 1.0e-11 | At04 | DATA, COMPUTED |
| NO + HO2 -> NO2 + OH | 3.6e-12, 0, -270 | 1.2e-11 | 8.9e-12 | At04 | DATA, COMPUTED |
| NO2 + OH + M -> HNO3 + M | k0: 3.968e-23, -2.9, 0; kinf 2.298e-9, -0.6, 0 | 6.0e-12 (M=1e18) | 3.5e-11 (M=2.5e19) | At97 | DATA, COMPUTED |
| NO + OH + M -> HNO2 + M | k0: 1.930e-24, -2.6, 0; kinf 6.368e-11, -0.1, 0 | n/a | n/a | De97 | DATA |
| N2O formation: NO + HNO -> N2O + OH; O1D + N2 + M -> N2O + M | in file | not extracted | not extracted | n/a | DATA |

ClOx family (partly missing):

| Reaction | A, b, Ea (K) | k at 298 K | File ref | Status |
|---|---|---|---|---|
| Cl + O3 -> ClO + O2 | 2.8e-11, 0, 250 | 1.2e-11 (COMPUTED) | At07 | DATA, COMPUTED |
| ClO + O + M -> OClO + M | present (type falloff) | not extracted | n/a | DATA |
| ClO + hv -> Cl + O, Cl + O1D | photolysis entries present | see section 2.4 | n/a | DATA |
| ClO + O -> Cl + O2 (the catalytic closing step) | absent from the file | n/a | n/a | NOT REACHED |
| ClO + NO -> Cl + NO2 | absent from the file | n/a | n/a | NOT REACHED |
| ClO + ClO channels, ClONO2 formation | absent from the file | n/a | n/a | NOT REACHED |

Carbon, CO2 and CH4 family:

| Reaction | A, b, Ea (K) | k at 220 K | k at 298 K | File ref | Status |
|---|---|---|---|---|---|
| CH4 + OH -> CH3 + H2O | 2.578e-17, 1.83, 1396 | 8.75e-16 | 8.0e-15 | Ba92 | DATA, COMPUTED |
| CO + OH -> CO2 + H | 3.662e-19, 1.89, -583 | 1.4e-13 | 1.2e-13 | Li07 | DATA, COMPUTED |
| CO + O + M -> CO2 + M | k0: 1.7e-33, 0, 1510; kinf 1e-14, 0, 1630 | 1.4e-18 (M=1e18) | 3.6e-17 (M=2.5e19) | Ts86, To84 | DATA, COMPUTED |
| O + CO2 -> CO + O2 | 2.8e-11, 0, 26500 | negligible | negligible | Ts86 | DATA |
| Cl + CH4 -> HCl + CH3 | 1.79e-15, 1.27, 940 (as extracted) | not computed | not computed | n/a | DATA |
| CH3 + O2 + M -> CH3O2 + M | falloff entry present | not extracted | not extracted | n/a | DATA |
| HCN + OH -> CO + NH2 | 1.07e-13, 0, 5890 (as extracted) | not computed | not computed | n/a | DATA |

Sulfur, nitrogen and other rows from the same file:

| Reaction | Parameters | File ref | Status |
|---|---|---|---|
| SO2 + OH + M -> HSO3 + M | k0: 5.978e-23, -3.3, 0; kinf 7.046e-11, -0.7, 0; k = 6.6e-13 at 220 K with M=1e18 and 1.16e-12 at 298 K with M=2.5e19 | At97, At04 | DATA, COMPUTED |
| HSO3 + O2 -> SO3 + HO2; SO3 + H2O + H2O -> H2SO4 | present | n/a | DATA |
| N2H4 + OH | 4.34e-8 T^-1.33 (as extracted) | n/a | DATA |
| N2O + OH | 2.03e-26 T^4.72 exp(-18400/T) (as extracted) | n/a | DATA |
| N + O2, N + NO | present | n/a | DATA |
| H2O2 + OH | present | n/a | DATA |

The file holds 95 photolysis entries and 4 Mie-type entries, plus 3 three-body and 103 falloff entries. Cross sections and quantum yields come from the per-species HDF5 files. Table 4 of Wogan 2025 lists the cross-section sources (Heays 2017, Huebner and Mukherjee 2015, Keller-Rudek 2013 and others) (READ).

### 2.2 What the reachable mechanism cannot give

The JPL 19-5 evaluation was unreachable, so Troe falloff parameters (F_c, N) and the ClOx reactions listed above have no verified source. The generator should treat the Lindemann form as an upper bound on the transition-region rate and flag ClOx rates NOT REACHED in the "Boss must research" list.

### 2.3 Stellar UV bands used for lifetimes

Method (COMPUTED). Each spectrum was scaled to a bolometric flux of 1.361e6 erg cm-2 s-1 at the planet (S = 1, 1361 W m-2 treated as the unit). The Sun uses the solar-minimum reference week; its tail beyond 2400 nm (3.8% of the flux) was subtracted before scaling. The MUSCLES tails beyond 5.5 or 13 micrometres were ignored (about 10% bolometric uncertainty). The F-star row is a 6500 K blackbody valid only above 200 nm (flux set to zero below 200 nm), because no F-star FUV spectrum was reachable. FUV-driven numbers for F stars are NOT REACHED; the F column below is useful only for wavelengths above 200 nm.

Band fluxes at S = 1, in W m-2 (COMPUTED from the files in section 0; DATA inputs):

| Band | F (BB 6500 K) | Sun (solar min) | K2V eps Eri | M1.5 GJ 832 | M2.5 GJ 176 | M4 GJ 876 | M4 GJ 699 | M8 TRAPPIST-1 |
|---|---|---|---|---|---|---|---|---|
| XUV 1-120 nm | not meaningful | 2.76e-3 | 7.78e-2 | 6.72e-2 | 4.54e-2 | 3.20e-2 | 3.22e-2 | 2.25 |
| Ly-alpha 121-122 nm | not meaningful | 5.83e-3 | 8.69e-2 | 6.11e-2 | 4.34e-2 | 2.98e-2 | 4.32e-2 | 0.185 |
| FUV 120-200 nm | not meaningful | 9.73e-2 | 0.113 | 6.69e-2 | 5.36e-2 | 4.01e-2 | 4.91e-2 | 0.902 |
| MUV 200-280 nm | 50.9 | 6.19 | 1.54 | 5.74e-2 | 5.68e-2 | 2.13e-2 | 9.2e-3 | 1.30 |
| NUV 280-320 nm | 51.6 | 20.0 | 12.9 | 0.289 | 0.162 | 4.22e-2 | 6.8e-2 | 0.85 |
| NUV 320-400 nm | 141 | 79.5 | 48.3 | 5.37 | 3.22 | 1.46 | 1.31 | 0.643 |

Note: the M-dwarf MUV and NUV entries below about 300 nm in the MUSCLES files are reconstruction or model values and are the least certain part of the table. The quiescent Sun is a solar-minimum week, so the Sun values are low for an average solar cycle. Stellar activity scaling with age is in section 4.5.

### 2.4 Photolysis rates and lifetimes versus host star

Unattenuated top-of-atmosphere photolysis rate J (s-1, Sun overhead, S = 1; COMPUTED from DATA cross sections, branches and the spectra above). Lifetime = 1/J, unattenuated. Real lifetimes in a column are longer (self-shielding and overlying absorbers, section 2.5) and the surface-level value is far longer for FUV-driven species. Diurnal averaging multiplies lifetimes by about 2 to 4 (a day-side mean of cos(zenith) at the equator on an equinox gives a factor pi, a whole-sphere mean gives 4; COMPUTED by geometry; valid for optically thin absorption only). Scale J linearly with the stellar constant S.

| Reaction | F BB6500 (>200 nm only) | Sun | K2V | M1.5 | M2.5 | M4 GJ 876 | M4 GJ 699 | M8 T-1 |
|---|---|---|---|---|---|---|---|---|
| O2 -> O + O | 6.9e-9 | 9.9e-7 | 6.3e-6 | 1.5e-6 | 2.0e-6 | 2.3e-6 | 1.8e-6 | 5.4e-5 |
| O2 -> O + O1D | 0 | 2.0e-6 | 7.3e-6 | 1.8e-6 | 2.6e-6 | 2.7e-6 | 2.0e-6 | 9.4e-5 |
| O3 -> O + O2 | 5.6e-3 | 1.2e-3 | 6.0e-4 | 1.8e-4 | 1.1e-4 | 5.5e-5 | 6.1e-5 | 1.2e-4 |
| O3 -> O1D + O2 | 4.7e-2 | 6.9e-3 | 2.3e-3 | 1.8e-4 | 1.4e-4 | 8.1e-5 | 9.1e-5 | 1.5e-3 |
| H2O2 -> 2 OH | 8.0e-4 | 1.0e-4 | 8.6e-5 | 4.0e-5 | 3.0e-5 | 2.3e-5 | 2.8e-5 | 2.6e-4 |
| NO2 -> NO + O | 1.4e-2 | 7.4e-3 | 4.5e-3 | 5.4e-4 | 3.4e-4 | 1.7e-4 | 1.6e-4 | 7.1e-4 |
| CO2 -> CO + O | 1.3e-10 | 8.6e-9 | 3.3e-9 | 6.2e-10 | 9.8e-10 | 6.6e-10 | 1.2e-9 | 1.8e-7 |
| CH4 -> CH3 + H | 0 | 3.2e-6 | 4.6e-5 | 3.1e-5 | 2.2e-5 | 1.6e-5 | 2.2e-5 | 1.2e-4 |
| H2O -> OH + H | 4.7e-9 | 9.9e-6 | 8.0e-5 | 4.9e-5 | 3.7e-5 | 2.6e-5 | 3.4e-5 | 3.0e-4 |

Lifetimes 1/J of the total photolysis (sum over branches), unattenuated (COMPUTED). Entries "inf" mean zero J in the supplied spectrum (the F column has no flux below 200 nm; for species whose absorption lies in the FUV this "inf" is an artefact and must not be used). Unit abbreviations: s seconds, h hours, d days, y years.

| Species | F (>200 nm) | Sun | K2V | M1.5 | M2.5 | M4 GJ 876 | M4 GJ 699 | M8 T-1 |
|---|---|---|---|---|---|---|---|---|
| O2 | artefact | 3.9 d | 20 h | 3.6 d | 2.6 d | 2.3 d | 3.0 d | 1.9 h |
| O3 | 19 s | 124 s | 345 s | 0.8 h | 1.1 h | 2.0 h | 1.8 h | 0.2 h |
| H2O | artefact | 25 h | 2.8 h | 4.5 h | 6.0 h | 8.5 h | 6.5 h | 0.8 h |
| H2O2 | 0.3 h | 2.7 h | 3.2 h | 7.0 h | 9.3 h | 12 h | 9.9 h | 1.1 h |
| CO2 | artefact | 8.6 d | 8.2 h | 22 h | 30 h | 2.1 d | 26 h | 4.5 h |
| CO | artefact | 18 d | 12 h | 2.4 d | 3.1 d | 5.8 d | 1.8 d | 7.6 h |
| H2 | artefact | 66 d | 35 h | 9.2 d | 12 d | 25 d | 6.0 d | 22 h |
| N2 | artefact | 22 d | 13 h | 2.9 d | 3.7 d | 7.4 d | 3.0 d | 12 h |
| CH4 | artefact | 32 h | 2.1 h | 3.5 h | 4.8 h | 6.6 h | 4.6 h | 0.8 h |
| C2H6 | artefact | 26 h | 1.8 h | 2.8 h | 3.8 h | 5.2 h | 3.7 h | 0.6 h |
| C2H2 | 27 h | 13.5 h | 1.1 h | 1.8 h | 2.4 h | 3.3 h | 3.3 h | 0.4 h |
| C2H4 | artefact | 8.1 h | 1.6 h | 2.7 h | 3.4 h | 4.6 h | 3.6 h | 0.2 h |
| CH3OH | 2.0 d | 26 h | 2.7 h | 4.6 h | 6.1 h | 8.3 h | 6.0 h | 0.8 h |
| H2CO | 0.7 h | 1.5 h | 1.5 h | 6.0 h | 7.8 h | 11 h | 8.9 h | 0.7 h |
| NH3 | 571 s | 1.9 h | 3.1 h | 7.1 h | 8.8 h | 12 h | 9.4 h | 0.5 h |
| N2H4 | 418 s | 1.3 h | 2.5 h | 5.0 h | 6.4 h | 8.5 h | 6.9 h | 0.4 h |
| HCN | artefact | 20 h | 1.4 h | 2.0 h | 2.8 h | 4.0 h | 2.8 h | 0.5 h |
| N2O | 2.8 d | 2.7 d | 7.4 h | 12.5 h | 15 h | 18 h | 15 h | 1.4 h |
| NO2 | 72 s | 135 s | 222 s | 0.5 h | 0.8 h | 1.7 h | 1.8 h | 0.4 h |
| NO | artefact | 33 h | 14 h | 2.3 d | 2.6 d | 3.5 d | 2.2 d | 2.6 h |
| HNO3 | 0.2 h | 1.5 h | 1.4 h | 2.9 h | 3.7 h | 5.3 h | 3.6 h | 573 s |
| SO2 | 390 s | 1.4 h | 1.1 h | 1.8 h | 2.4 h | 3.3 h | 2.5 h | 0.3 h |
| H2S | 319 s | 1.0 h | 0.8 h | 1.4 h | 1.8 h | 2.4 h | 2.1 h | 0.2 h |
| SO | 181 s | 0.6 h | 0.4 h | 0.7 h | 0.9 h | 1.2 h | 0.9 h | 493 s |
| SO3 | 0.6 h | 7.1 h | 1.9 d | 44 d | 33 d | 57 d | 73 d | 6.7 h |
| H2SO4 (gas) | artefact | 11 h | 0.8 h | 1.1 h | 1.5 h | 2.2 h | 1.5 h | 0.3 h |
| OCS | 0.6 h | 4.0 h | 2.5 h | 5.3 h | 4.7 h | 5.3 h | 7.1 h | 0.3 h |
| CS2 | 46 s | 547 s | 0.7 h | 2.7 h | 3.1 h | 4.4 h | 3.2 h | 203 s |
| HCl | 3.1 d | 2.3 d | 9.8 h | 19 h | 23 h | 30 h | 19 h | 1.3 h |
| ClO | 30 s | 152 s | 393 s | 3.8 h | 4.3 h | 15 h | 26 h | 0.5 h |
| HOCl | 0.2 h | 0.6 h | 0.8 h | 7.3 h | 9.8 h | 15 h | 9.2 h | 0.5 h |
| NO3 | 4 s | 4 s | 4 s | 8 s | 14 s | 31 s | 28 s | 225 s |
| HO2 | 112 s | 0.4 h | 2.3 h | 2.4 d | 2.0 d | 3.0 d | 5.6 d | 0.6 h |
| OH | 7.0 h | 36 h | 8.9 h | 14 h | 19 h | 26 h | 21 h | 2.3 h |
| S8 | no photolysis entry in this mechanism (treated as aerosol, section 2.6) | | | | | | | |

Reading the table. (1) The Sun column is low because the WHI week is solar minimum. (2) For K2V and M dwarfs the FUV and Ly-alpha dominate, so CH4, H2O, CO2, HCN, NH3 and hydrocarbons photolyse several to tens of times faster than for the Sun, while the NUV-driven species (O3, NO2, H2O2, HOCl) are slower by factors of 2 to 10. (3) TRAPPIST-1 (M8) has strong XUV and Ly-alpha, so every FUV-driven species has lifetimes of an hour or less unattenuated. (4) The N2, CO and H2 lifetimes of days come from EUV below 100 nm, which is absorbed above about 100 km altitude; they are not column photolysis lifetimes at the surface. (5) The O2 rows sum all O2 photolysis branches in the file at the top of the atmosphere; they set the Chapman source (section 3) before any self-shielding.

### 2.5 UV shielding data and the critical-pressure rule

Band-mean absorption cross sections from the Photochem HDF5 files (DATA) and the mixing ratio needed for a vertical optical depth of 1 through a 1 bar N2 column, N_tot = 2.22e25 cm-2 (COMPUTED). Entries "none" mean no data in that band.

| Species | 160-200 nm (cm2) | 200-242 nm | 242-290 nm | 290-340 nm | Mixing ratio for tau = 1 through 1 bar | Status |
|---|---|---|---|---|---|---|
| O3 | 7.9e-19 | 3.0e-18 | 7.8e-18 | 2.0e-19 | about 6e-9 to 2e-8 in the first three bands | DATA, COMPUTED |
| SO2 | 3.2e-18 | 2.1e-18 | 4.3e-19 | 2.6e-19 | about 1.4e-8 to 1.1e-7 in the first three bands | DATA, COMPUTED |
| H2S | 3.1e-18 | 2.0e-18 | 2.1e-20 | 8e-23 | not computed | DATA |
| O2 | 8.5e-19 | 4.2e-24 | none | none | not applicable | DATA |
| CO2 | 1.9e-20 | 1.2e-25 | none | none | not applicable | DATA |
| H2O | 2.3e-18 | 4.2e-24 | none | none | not applicable | DATA |
| NH3 | 5.5e-18 | 1.4e-18 | none | none | not computed | DATA |
| S8 (gas or aerosol cross section) | 2.0e-17 | 2.8e-17 | 1.4e-17 | 4.2e-18 | about 2e-9 | DATA, COMPUTED |

Consequence: parts per billion of SO2 or O3 (or S8 aerosol) block the 200-290 nm window for a 1 bar atmosphere, so an anoxic planet with volcanic SO2 or a sulfur haze gets a UV shield without any biology. Whether the shield blocks the surface 200-300 nm flux must be computed for the actual column; the band means above ignore line structure and temperature dependence.

Critical-pressure rule (Hu, Seager and Bains 2012, arXiv:1210.6885, Appendix D, eq. D5, READ): a gas of mixing ratio f is photolysed efficiently only above the level where P < P* = mu0 g M_atm / (sigma_UV F), where sigma_UV is a mean UV cross section and F the stellar photon flux factor defined in the paper. Table 6 of the paper gives mean cross sections and 1-ppm critical pressures P* in Pa (READ): H2 1.10e-17 cm2, 207 Pa; O2 7.18e-21, 3.17e5; H2O 1.72e-18, 1330; NH3 1.31e-18, 1730; N2O 1.38e-20, 1.65e5; CH4 1.81e-17, 126; CO2 7.44e-22, 3.06e6; H2S 1.11e-18, 2050; OCS 2.74e-20, 8.33e4; SO2 3.74e-18, 609. A gas with P* larger than the surface pressure is effectively shielded and is not photolysed in the lower atmosphere.

### 2.6 Cycles and family rules

CO2, CO and O2. The hydrogen-oxide cycle (OH, HO2, H2O2) recombines CO and O into CO2 catalytically. Gao et al. 2015 (arXiv:1501.06876, READ) show this cycle breaks down when H is below about 1 ppm and the direct CO + O -> CO2 path takes over at below 0.1 ppm. In their desiccated M-dwarf CO2 atmosphere (GJ 436 spectrum, Mars-like 1 bar, 240 K surface, diffusion-limited H escape) they found abiotic O2 up to about 20%, O3 above 1 ppm, CO of 10 to 40%, and a cap of about 50% on CO2 destruction. Hu et al. 2012 (READ) found abiotic O2 in a 1 bar CO2 atmosphere of 6.4e-7, 3.8e-6 or 1.3e-3 for H2 emission of 3e10, 3e9 or 0 molecules cm-2 s-1, with O3 columns up to a third of the Earth column, and CO lifetimes of 0.1 to 700 years across the benchmark atmospheres. Schwieterman et al. 2019 (arXiv:1903.07661, READ) find that biospheres can maintain about 100 ppmv CO even at low H2 flux, and around M dwarfs CO of hundreds of ppm to several percent, with a CO runaway for Proxima at the highest fluxes; the model fluxes were CH4 1.5e11 and CO 3.0e11 molecules cm-2 s-1 with a maximum CO deposition velocity of 1.2e-4 cm s-1. Earth-like modern CO of about 0.1 ppmv is a TEST value from the same paper. Harman et al. 2015 (the original paper, not on arXiv) is represented only through Meadows et al. 2018 and Gao 2015 (2ND-HAND).

Methane. Segura et al. 2005 (astro-ph/0510224, READ): methane lifetime of 10 to 12 years on an Earth-like planet versus about 200 years for an AD Leo planet; OH is about 1e6 cm-3 on the Earth case versus about 1 cm-3 for active M stars and about 1e-2 cm-3 for quiescent M stars. Hu et al. 2012 benchmark CH4 chemical lifetimes of 8e4 years (H2 atmosphere) and 6e4 years (N2 or CO2) at low OH (READ). Arney et al. 2016 (arXiv:1610.04515, READ): haze thickness scales with CH4/CO2; CH4/CO2 near 0.1 gives negligible haze and a thick haze needs about 0.2; the CH4 flux required to start haze is about 1e11 to 3e11 molecules cm-2 s-1; the pCO2 range was 0.0036 to 0.018 bar for Psurf 0.5 to 1 bar; monomer density 0.64 g cm-3 (Trainer 2006 value, 2ND-HAND). Serpentinization CH4 flux maxima are 6.8e8 (1 Earth mass) and 1.3e9 (5 Earth masses) molecules cm-2 s-1 (Guzman-Marmolejo 2013, 2ND-HAND via Arney 2016), far below the haze-starting flux, so abiotic methane alone does not trigger haze. The haze threshold "CH4/CO2 about 0.1" quoted in the project brief traces to Trainer 2006 and Zerkle 2012; those originals were not reached, so the generator should draw the threshold between 0.1 and 0.2 (Arney 2016, READ) and carry the uncertainty. Arney 2016 also gives an NH3 lifetime in the haze case of about 1e4 years versus 7e7 years from Wolf and Toon (tau about 11 case; READ in Arney 2016, original not reached).

Sulfur. Hu, Seager and Bains 2013 (arXiv:1302.6603, READ): tropospheric lifetime of H2S and SO2 on Earth is about 2 days; Earth's volcanic sulfur flux is about 3e9 S cm-2 s-1 and the volcanic H2S/SO2 ratio is about 0.1 (TEST values). Sulfur haze (S8 or H2SO4) becomes optically thick at 100 times the volcanic flux or more in N2 and CO2 atmospheres, and at an Earth-like flux in H2 atmospheres. H2S and SO2 accumulate to detectable levels only above about 1000 times the Earth flux. Dry deposition velocities used there: SO2 1 cm s-1, H2S 0.015, O3 0.4, H2O2 0.5, S8 and H2SO4 aerosols 0.2 (READ; use as defaults). Zahnle et al. 2006 was not reached.

Nitrogen. NO and NO2 cycle with ozone (reactions in 2.1); NOx is converted to HNO3 and rained out. N2O arises from N(4S) + NO2 -> N2O + O and NO + NH -> N2O + H in SEP-driven early-Earth models (Airapetian et al. 2020, arXiv:1905.05093, READ; the text states a production more than 300 times the earlier model, an Earth-like biological level of about 400 ppbv is quoted there as a TEST value). HCN forms in neutral reactions NO + CH -> HCN + O, CH2 + N(4S) -> HCN + H, CH3 + N(4S) -> HCN + 2 H and CH + CN -> HCN + H (same source, READ), and the same paper states that SEP-driven HCN production exceeds lightning production in a weakly reducing N2-CO/CO2-CH4-H2O atmosphere (READ, laboratory results cited there). The reaction HCN + OH has Ea about 5890 K (DATA), so HCN is long-lived against OH and is lost chiefly by FUV photolysis (lifetime table above) and rain-out.

Hydrocarbons and tholins. Horst 2017 (arXiv:1702.08611, READ): Titan haze production about 3.5e-14 g cm-2 s-1 (Podolak and Bar-Nun) and about 1e-14 in later work; ionospheric mass flux 3.2e-15; CH4 destroyed in about 30 Myr; N2/CH4/CO atmospheres favour particle production; "tholin" means the laboratory analogue only. Wogan 2025 (READ) reports a Titan haze production of about 7e-15 g cm-2 s-1 and an N2 photolysis that is about 50% UV and 50% galactic cosmic rays. These are TEST points, not generator constants. Condensate and haze optics are Part 06.

Phosphine. Venus data are TEST values only: Greaves et al. 2020 (arXiv:2009.06593) report about 20 ppb PH3 (READ). Bains et al. 2021 (arXiv:2009.06499, READ): PH3 lifetime below 1 s above 78 to 98 km, up to 1e11 s in the UV-shielded deep atmosphere, with a transport limit near 400 to 700 years. Sousa-Silva et al. 2020 (arXiv:1910.05224, READ): PH3 + OH k = 2e-11 cm3 s-1 at 288 K (A = 2.71e-11, E = 1.29 kJ mol-1), lifetime 28 h at night and 5 h by day, PH3 cross sections from Chen 1991 with quantum yield 1 below 230 nm. The shipped Photochem mechanism does not contain phosphorus, so there is no lifetime-versus-star row for PH3 (NOT REACHED).

Hydrogen peroxide. H2O2 forms from HO2 + HO2 and OH + OH + M (2.1) and is lost by photolysis (lifetimes 0.3 to 12 h unattenuated, table 2.4) and deposition (0.5 cm s-1, Hu 2013). In dry CO2 atmospheres the H2O2 route is the main catalytic CO2 recombination channel (Gao 2015).

## 3. Ozone

### 3.1 Mechanism and what the generator needs

Chapman scheme (textbook, RECALLED, listed in section 9): O2 + hv -> 2 O; O + O2 + M -> O3 + M; O3 + hv -> O + O2; O + O3 -> 2 O2. The O + O2 + M and O + O3 rate coefficients are DATA in section 2.1. Odd-oxygen catalysts are NOx (NO + O3, NO2 + O), HOx (OH + O3, HO2 + O, H + O3) and ClOx (Cl + O3); the first two have verified rows in 2.1, the closing ClOx step is NOT REACHED.

The Chapman steady-state scaling, the Kozakis 2022 fit coefficients and the Cooke, Segura and Rugheimer cross-checks, and the generator guidance ("Kozakis fit for the median, interpolate T_eff, seeded spread of about 30 percent") are in `atmosphere-sources-check.md`; they are cited, not repeated.

### 3.2 Published column values (calibration and test points)

| Item | Value | Source | Status |
|---|---|---|---|
| Maximum O3 column and the O2 level at which it occurs, by host (Kozakis 2022 Table 2): G0V | 7.74e18 cm-2 at 25% PAL | Kozakis et al. 2022, arXiv:2208.09415 | READ |
| Same, Sun | 5.64e18 cm-2 at 25% PAL | same | READ |
| Same, K2V | 4.37e18 cm-2 at 55% PAL | same | READ |
| Same, K5V | 4.96e18 cm-2 at 150% PAL | same | READ |
| Same, M5V | 4.65e18 cm-2 at 150% PAL | same | READ |
| O3 layer moves downward as O2 decreases (more strongly for hotter hosts); a secondary layer sits in the mesosphere; the pressure dependence of O + O2 + M puts the maximum at lower O2 for hotter stars | qualitative | same | READ |
| 3-D model O3 column at 100, 10, 1 and 0.1% PAL O2 | 279, 169, 66 and 18 DU | Cooke et al. 2022, arXiv:2102.11675 | READ |
| Ratio of other models to Cooke at 10% PAL (Way, Segura, Kasting and Donahue, Levine) | 1.46, 1.57, 1.76, 2.43 | same | READ |
| Same ratios at 1% PAL | 1.83, 1.87, 2.24, 2.89 | same | READ |
| TEST range for a modern Earth ozone column (O2 = 100% PAL, Sun) | 279 to 300 DU | Cooke 2022 (279); the 300 DU upper bound is the generator test range in the task brief | READ (279); the 300 DU edge is not independently verified |
| Kozakis fit for the Sun at 100% PAL | about 192 DU | `atmosphere-sources-check.md` (project doc), fit coefficients from Kozakis 2022 Table 3 | 2ND-HAND (project doc) |
| AD Leo ozone column in the Segura et al. runs | 1.67e19 molecules cm-2; Earth case 8.25e18 | Segura et al. 2010 (arXiv:1006.0022) and Segura 2005 | READ |
| O3 column in the strongest depletion case | 1.1e18 cm-2 | Segura et al. 2010 | READ |
| 1 DU conversion | 2.687e16 molecules cm-2 | standard constant (Loschmidt number times 1 cm3 per 1 cm2 column at STP, 2.687e19 cm-3 x 1e-3 cm) | COMPUTED |

The spread between the Kozakis 1D result (about 192 DU for the Sun at 100% PAL) and the 3-D Cooke value (279 DU) is a factor 1.45 and should be treated as the model-structure uncertainty. The generator should draw the column from a distribution spanning both.

### 3.3 Vertical profile and layer position

COMPUTED indicative model (O2 photolysis only, no O3 self-shielding, isothermal atmosphere, overhead Sun, solar-minimum spectrum, Chapman layer production per unit log pressure peaks near one scale height times ln(column/unit)). The pressure of the peak O2 photolysis rate is:

| O2 level (PAL) | Peak production pressure (hPa) |
|---|---|
| 1.5 | 42 |
| 1.0 | 63 |
| 0.25 | 252 |
| 0.1 | 630 |
| 0.01 | 6300 |

The K2V case sits about 1.25 times deeper. This is an upper bound on the depth because O3 self-shielding is omitted and temperature structure is ignored; it only shows the monotonic downward shift as O2 falls, which agrees in direction with Kozakis 2022 (READ). A Chapman production layer in log pressure has a full width at half maximum of 2.45 scale heights (COMPUTED). The M-dwarf results of this toy model (peaks at 0.01 to 11 hPa) are dominated by the FUV and Ly-alpha and are not reliable. The observed ozone-layer peak near 20 to 25 km on Earth is RECALLED (section 9); it is a calibration point only.

### 3.4 Host-star dependence

From the J-value table (2.4): the O3 photolysis lifetime (unattenuated) rises from 124 s (Sun) to 0.8 to 2.0 h for M dwarfs, and the O2 photolysis lifetime is 3.9 d (Sun), 20 h (K2V), 2 to 4 d (mid M dwarfs). The ratio of O3 destruction to O2 photolysis therefore falls by a factor of about 10 to 100 around M dwarfs, which is the origin of the large M-dwarf ozone columns in Segura 2005 and Kozakis 2022. The generator rule: ozone column = Kozakis fit at the host T_eff, then multiplied by a factor from a seeded spread; the table here gives the physical direction and the sizes of the J ratios.

### 3.5 Ozone shielding by other gases on anoxic worlds

Section 2.5 gives the mixing ratios at which SO2 (1.4e-8 to 1.1e-7) and O3 (6e-9 to 2e-8) become opaque at 200 to 290 nm for a 1 bar N2 column. H2O and CO2 do not shield beyond about 200 nm and 220 nm respectively (cross sections 4.2e-24 and 1.2e-25 cm2 at 200-242 nm, DATA). Hu 2013 gives the sulfur fluxes needed for S8 or H2SO4 haze (section 2.6, READ).

### 3.6 Abiotic O2 and O3 false positives

From Meadows et al. 2018 (arXiv:1705.07560, READ):

| Channel | Magnitude | Notes | Source cited there | Status |
|---|---|---|---|---|
| Low non-condensable inventory (water vapour photolysis, Wordsworth and Pierrehumbert 2014) | about 0.15 bar O2 | needs little N2 and little Ar | Wordsworth and Pierrehumbert 2014 | 2ND-HAND |
| Pre-main-sequence ocean loss around M dwarfs | hundreds to thousands of bar O2 | one Earth ocean gives about 240 bar O2; Proxima b gives 1 to 3 ocean masses (below 250 up to 750 to 2500 bar O2) | Luger and Barnes 2015 | 2ND-HAND (numbers read in Meadows 2018; the original arXiv:1411.7412 was downloaded but not extracted) |
| CO2 photolysis without sinks | O2 from 0.2% to 6%, reduced by orders of magnitude with sinks; Gao 2015 about 15% | | Tian 2014, Harman 2015, Gao 2015 | 2ND-HAND; Gao READ |
| O3 without O2 buildup | O3 at about 10% of the Earth column | | Domagal-Goldman 2014 | 2ND-HAND |
| TiO photocatalysis | | | no source | 2ND-HAND |

The numbers differ between sources by orders of magnitude because they depend on sinks (volcanic reductants, ocean oxidation), so the generator must treat them as scenario flags (set by M-dwarf pre-main-sequence phase length, ocean inventory and N2 inventory) and not as fixed levels.

## 4. Escape

The shoreline and the project's Jeans convention (30 retained, 15 lost) are in `atmospheres-retention-and-classes.md` and are only cited here. This section adds thresholds from the literature, energy-limited and hydro-fit formulas, the diffusion limit, nonthermal channels and stellar XUV laws.

### 4.1 Jeans parameter, thermal escape and thresholds

Jeans parameter at the exobase: lambda = G M m / (k T r), equal to X = v_esc^2 / v0^2 = R/H in Wurz et al. 2022 (Space Sci. Rev. 218, 10, DOI 10.1007/s11214-022-00875-6, READ). Jeans flux (Gronoff et al. 2020, arXiv:2003.03231, READ): Phi = N u / (2 sqrt(pi)) (1 + lambda) exp(-lambda), with N the exobase density and u the most probable speed. Wurz et al. 2022 write the escape fraction of the thermal flux as f_esc = (1 + X) exp(-X) / (2 sqrt(pi)) and the thermal release flux as Phi_th = n_i sqrt(kT / (2 pi m)) (READ).

| Threshold | Value | Source | Status |
|---|---|---|---|
| Hydrodynamic (blow-off) escape for the exobase parameter | lambda_ex below 2.5 | Gronoff 2020 | READ |
| Blow-off regime in X | X below 3/2; purely hydrodynamic for X_crit up to 2 to 3; purely Jeans for X_crit of 6 or more | Wurz 2022 | READ |
| A species is lost over geological time | X below 15 | Wurz 2022 | READ |
| Exobase non-Maxwellian corrections to Jeans | minus 25% (Chamberlain and Smith) and plus 30% (Merryfield and Shizgal) | Gronoff 2020 | READ |
| Project convention | 30 retained, 15 lost | project doc | cited, not redone |

The values f_esc(X) from the Wurz formula (COMPUTED): X = 1: 2.1e-1; 1.5: 1.6e-1; 3: 5.6e-2; 6: 4.9e-3; 10: 1.4e-4; 15: 1.4e-6; 20: 1.2e-8; 30: 8.2e-13. Wurz gives v_esc = 4250 m/s for Mercury and 2375 m/s for the Moon (READ).

Exobase parameters lambda_ex in Gronoff 2020 Table 2 (READ; these are TEST or calibration rows, not generator constants): Mercury 2.2, Venus 22.3, Earth 9.4 to 5.0, Moon 0.8, Mars 6.3 to 5.0, Titan 2.3 (CH4: 37.3), Pluto 15.1 (CH4: 8.5), Jupiter 311 to 218, Saturn 157 to 98, Uranus 34, Neptune 48. Exobase temperatures in the same table: Titan 180 K, Pluto 68 K, Uranus 800 K, Neptune 700 K, Saturn 500 to 800 K (READ).

### 4.2 Energy-limited and hydro-fit XUV escape

| Formula | Statement | Validity and notes | Source | Status |
|---|---|---|---|---|
| Energy-limited | Mdot = eta pi R_p^3 F_XUV / (G M_p K_tide) (Lopez and Fortney form with R_XUV in place of R_p in the numerator when the XUV photosphere is above the radius) | eta = 0.1 default | Lopez, Fortney and Miller 2012 (arXiv:1205.0010); Lopez and Fortney 2013 (arXiv:1305.0269) | READ |
| Tidal factor | K_tide = 1 - 3/(2 xi) + 1/(2 xi^3), xi = R_Hill / R_XUV | K_tide near 1 unless the planet is close to its Roche lobe | same | READ |
| R_XUV | radius where the atmosphere becomes optically thick to XUV (about 1 nbar) | | Lopez and Fortney 2013 | READ |
| XUV flux used | F_XUV (1 to 1200 A) from Ribas 2005 | mass loss starts at 10 Myr | Lopez and Fortney 2013 | READ |
| Efficiency range | eta from 0.05 to 0.2 | Owen and Jackson 2012 values quoted | Lopez and Fortney 2013 | 2ND-HAND (Owen and Jackson 2012 arXiv:1206.2367 was downloaded but not read) |
| Variable efficiency | eta = 0.1 (v_esc / 15 km s-1)^-2 | fits the radius-valley placement | Owen and Wu 2017 (arXiv:1705.10810) | READ |
| Owen and Wu mass loss | Mdot = eta pi R_p^3 L_HE / (4 pi a^2 G M_p) | L_HE = L_sat (t / t_sat)^(-1-a0) after saturation, form following Jackson 2012 | same | READ |
| Radius-valley scaling | R_valley proportional to eta^0.18 X_HE^0.19 rho^-0.24 | | same | READ |

Values of eta(v_esc) from the Owen and Wu law (COMPUTED): v_esc = 8 km s-1 gives 0.35; 15 km s-1 gives 0.10; 25 km s-1 gives 0.036.

Hydro-based fit of Kubyshkina et al. 2018 (arXiv:1810.06920, READ):

Mdot = exp(beta) F_XUV^alpha1 (d0/au)^alpha2 (R_pl/R_E)^alpha3 Lambda^K, with K = zeta + theta ln(d0/au), Lambda = G M_pl mH / (kB T_eq R_pl), Mdot in g s-1, F_XUV in erg cm-2 s-1.

| Branch | beta | alpha1 | alpha2 | alpha3 | zeta | theta |
|---|---|---|---|---|---|---|
| Lambda below exp(Sigma) | 32.0199 | 0.4222 | -1.7489 | 3.7679 | -6.8618 | 0.0095 |
| Lambda at or above exp(Sigma) | 16.4084 | 1.0000 | -3.2861 | 2.7500 | -1.2978 | 0.8846 |

Sigma = (15.611 - 0.578 ln F_XUV + 1.537 ln(d0/au) + 1.018 ln(R_pl/R_E)) / (5.564 + 0.894 ln(d0/au)).

Validity (READ): masses 1 to 39 Earth masses, radii 1 to 10 Earth radii, T_eq 300 to 2000 K, host mass 0.4 to 1.3 solar masses, orbital distance 0.002 to 1.3 au, L_XUV 1e26 to 5e30 erg s-1, built with eta = 15%, Lambda below 80, hydrogen-dominated atmospheres only. Do not apply to CO2- or N2-dominated atmospheres.

Worked checks (COMPUTED by me with the table above; not published values):

| M_p (M_E) | R_p (R_E) | T_eq (K) | d0 (au) | F_XUV (erg cm-2 s-1) | Lambda | exp(Sigma) | Kubyshkina Mdot (g s-1) | Energy-limited, eta = 0.15, R_XUV = R_p, K_tide = 1 (g s-1) |
|---|---|---|---|---|---|---|---|---|
| 5 | 1.8 | 700 | 0.1 | 460 | 30.1 | 13.5 | 7.0e8 | 1.6e8 |
| 10 | 2.5 | 900 | 0.05 | 1800 | 33.7 | 14.0 | 5.2e9 | 8.6e8 |
| 3 | 1.5 | 500 | 0.2 | 115 | 30.3 | 13.7 | 8.6e7 | 4.0e7 |

All three fall on the high-Lambda branch, and the hydro fit exceeds the simple energy-limited estimate by factors of 2 to 6 (because the energy-limited estimate used R_XUV = R_p; the larger R_XUV in real envelopes closes the gap). The generator should draw eta between 0.05 and 0.2 (or use the variable-efficiency law) and keep the hydro fit as a cross-check inside its validity box. Boil-off after disc dispersal is not covered by a source I could read (section 8).

### 4.3 Diffusion-limited escape (Hunten)

Hunten limit for a light trace gas i of total mixing ratio f in a heavier background (Gronoff 2020 form, READ): Phi = b_i,dom / H x n_i,homopause / n; the crossover mass with a trace flux is mc = m_esc + kT Phi / (b g X) (Hunten 1987 form as quoted by Gronoff). For m_i much smaller than the background, Phi_lim is approximately b f / H_a with H_a = kT / (m_a g) (the standard limiting form; COMPUTED usage below).

Binary diffusion coefficients (Hu, Seager and Bains 2012, arXiv:1210.6885, READ), D in cm2 s-1 with N the total number density (cm-3) and T in K, so b = D N:

| Pair | b = D N (cm-1 s-1) |
|---|---|
| H in N2 | 4.87e17 T^0.698 |
| H2 in N2 | 2.80e17 T^0.740 |
| H in CO2 | 3.87e17 T^0.750 |
| H2 in CO2 | 2.15e17 T^0.750 |
| H in H2 | 8.16e17 T^0.728 |

The same paper gives a thermal diffusion factor of -0.38 (READ).

Worked values (COMPUTED, illustrative inputs): H2 in N2, T = 250 K, f(H2) = 1e-4, g = 980 cm s-2, H_a = 7.6e5 cm gives Phi_lim = 2.2e9 cm-2 s-1; H in N2 with f = 1e-4 gives 3.0e9; H2 in CO2 at T = 200 K, g = 370, f = 1e-4 gives 1.1e9. The limit scales with the total hydrogen mixing ratio at the homopause (including H2O and CH4 H content) and is independent of the XUV flux, so thick atmospheres with a cold trap (low stratospheric H2O) lose hydrogen more slowly than XUV-limited values.

### 4.4 Nonthermal channels

Gronoff 2020 Table 1 lists the escape processes (READ): Jeans, hydrodynamic, charge exchange and trapped-ion escape, ion pickup, sputtering, photochemical escape, polar-wind ion outflow, unmagnetised ion outflow (Kelvin-Helmholtz and clouds), and bulk ion escape. Table 3 (READ) gives present escape for the three terrestrial planets in the format "total s-1 then flux cm-2 s-1" (TEST or calibration points only):

| Process | Venus | Earth | Mars |
|---|---|---|---|
| Jeans | 2.5e19 then 5.1 | H: 6e26 then 1e8 (solar maximum) | H: 1.6e26 to 1.1e27 then 1e8 to 6.9e8 |
| Charge exchange or trapped | H: 5e24 to 5e25 then 1e6 to 1e7 | H: 6e26 then 1e8 (solar minimum) | 1e22 to 1e23 then 1e4 to 1e5 |
| Ion pickup | H+: 1e25 then 2e6; O+: 1.5e25 then 3e6 | small | O+: 1e24 then 1e6; C+: 1.6e23 then 1e5 |
| Sputtering | O: 5e23 to 5e24 then 1e5 to 1e6 | small | O: 3e24 then 1.8e6 |
| Photochemical | 3.8e25 then 7.7e6 | small | O: 5e25 then 3e7; C: 1e24 then 1e6 |
| Magnetised ion outflow (polar wind) | not applicable | H+: 8e25 then 1.3e7; O+: 3e25 then 5e6 | not applicable |
| Unmagnetised ion outflow | 5e24 to 1e25 then 1e6 to 2e6 | O+: 6e24 then 1e6 | 1e25 then 1e7 |

Notes. The Earth Jeans flux at solar maximum equals the charge-exchange flux at solar minimum because hydrogen escape is diffusion limited (stated in the table caption, READ). Flux areas used by Gronoff: Venus exobase at 200 km, 4.9e18 cm2; Earth 500 km, 5.9e18 cm2; Mars 200 km, 1.6e18 cm2. The Moon and Mercury release channels are in section 5. Airapetian et al. 2020 (arXiv:1905.05093) and Lammer reviews describe the active-star ion escape; I did not extract transferable rate formulas from them (section 8).

### 4.5 Stellar XUV evolution laws

Ribas et al. 2005 (astro-ph/0412253, Table 5, READ): flux at 1 au in erg cm-2 s-1 equals alpha tau^beta with tau the age in Gyr.

| Interval | alpha | beta | F at 0.1 Gyr | F at 1 Gyr | F at 4.57 Gyr (COMPUTED) |
|---|---|---|---|---|---|
| 1-20 A | 2.40 | -1.92 | 200 | 2.40 | 0.13 |
| 20-100 A | 4.45 | -1.27 | 83 | 4.45 | 0.65 |
| 100-360 A | 13.5 | -1.20 | 214 | 13.5 | 2.18 |
| 360-920 A | 4.56 | (-1.0, value in parentheses in the table) | 46 | 4.56 | 1.00 |
| 920-1180 A | 2.53 | -0.85 | not computed | 2.53 | not computed |
| 1-360 A plus 920-1180 A | 24.8 | -1.27 | not computed | 24.8 | not computed |
| 1-1180 A | 29.7 | -1.23 | 504 | 29.7 | 4.58 |
| Ly-alpha (Eq. 2) | 19.2 | -0.72 | 101 | 19.2 | 6.4 |

Validity: solar-type (G) stars at 1 au. Ratio of the 1-1180 A flux to its 4.57 Gyr value (COMPUTED from the law): 2.7 at an age of 2.07 Gyr, 6.0 at 1.07 Gyr and 110 at 0.1 Gyr; Ribas quotes about 2.5, 6 and 100 for 2.5, 3.5 and 4.5 Gyr before the present (READ).

Jackson et al. 2012 (arXiv:1111.0031, Table 2, READ): L_X / L_bol = 10^sat for t at or below tau_sat and 10^sat (t / tau_sat)^-alpha afterwards.

| (B-V)0 bin | log(L_X/L_bol) sat | log tau_sat (yr) | alpha | L_X/L_bol at 1 Gyr and 4.5 Gyr (COMPUTED) |
|---|---|---|---|---|
| 0.29 to 0.45 | -4.28 | 7.87 | 1.22 | not computed |
| 0.45 to 0.565 | -4.24 | 8.35 | 1.24 | 9.0e-6, 1.4e-6 |
| 0.565 to 0.675 | -3.67 | 7.84 | 1.13 | 1.0e-5, 1.9e-6 |
| 0.675 to 0.79 | -3.71 | 8.03 | 1.28 | not computed |
| 0.79 to 0.935 | -3.36 | 7.90 | 1.40 | 1.3e-5, 1.5e-6 |
| 0.935 to 1.275 | -3.35 | 8.28 | 1.09 | not computed |
| 1.275 to 1.41 (early M) | -3.14 | 8.21 | 1.18 | 8.5e-5, 1.4e-5 |

Tu et al. 2015 (arXiv:1504.04546, READ): R_X = C Ro^beta with beta = -2.7, saturation Ro_sat = 0.13 and R_X,sat = 10^-3.13; t_sat = 2.9 Omega0^1.14 Myr. The 10th, 50th and 90th percentile initial rotation rates are 1.8, 6.2 and 45.6 times solar, with t_sat of 5.7, 23 and 226 Myr. For a solar-mass star, L_X (5-100 A, erg s-1, t in Myr) is 2.0e31 t^-1.12, 2.6e32 t^-1.42 and 2.3e36 t^-2.50 for the three percentiles; L_EUV (100-920 A) is 7.4e31 t^-0.96, 4.8e32 t^-1.22 and 1.2e36 t^-2.15. Values at 100, 1000 and 4570 Myr (COMPUTED): L_X 10th: 1.2e29, 8.7e27, 1.6e27; 50th: 3.8e29, 1.4e28, 1.7e27; 90th: 2.3e31, 7.3e28, 1.6e27 erg s-1. L_EUV 10th: 8.9e29, 9.8e28, 2.3e28; 50th: 1.7e30, 1.1e29, 1.6e28; 90th: 6.0e31, 4.3e29, 1.6e28. Check: L_X,sun = 10^27.2 at 4570 Myr (READ), consistent with the 1.6e27 above. The relation log L_EUV = 4.8 + 0.86 log L_X is attributed to Sanz-Forcada (2ND-HAND via Tu 2015). The range of rotation histories at 100 Myr is a factor of about 200 in L_X; the generator should draw the rotation percentile.

Johnstone et al. 2021 (arXiv:2009.07695, READ): R_X = C1 Ro^beta1 and C2 Ro^beta2 in two Rossby-number branches with slopes -0.135 +/- 0.030 and -1.889 +/- 0.079, Ro_sat = 0.0605 +/- 0.0033 and R_X,sat = 5.135e-4 +/- 3.3e-5. Which slope belongs to which branch must be rechecked (physically the shallow slope is the saturated branch; section 9). Surface-flux relations in erg cm-2 s-1: log F_EUV1 (10-36 nm) = 2.04 + 0.681 log F_X (eq. 19); log F_EUV2 (36-92 nm) = -0.341 + 0.920 log F_EUV1 (eq. 21); log F_Lya = 3.97 + 0.375 log F_X (eq. 23); F_EUV2 / F_EUV1 = 0.924 F_EUV1^-0.0798. Flare frequency scales as N(>E1) / N(>E2) = a power law in E (eqs. 25 to 27; exponent not transcribed).

The three families disagree at 4.5 Gyr by a few tens of percent and at 100 Myr by up to two orders of magnitude (Tu percentiles); that spread is real and should be drawn, not averaged.

## 5. Exospheres of airless bodies

All rows below are READ from Wurz et al. 2022 (Space Sci. Rev. 218, 10), Oza et al. 2019 (arXiv:1908.10732), Lammer et al. 2022 (arXiv:2203.01656), Lierle and Schmidt 2025 (arXiv:2505.23915), Vorburger et al. 2019 (arXiv:1909.01014) unless stated; Solar System bodies are calibration or test points.

### 5.1 Source and loss processes

| Process | Parameters | Source | Status |
|---|---|---|---|
| Thermal release | Phi_th = n_i sqrt(kT / (2 pi m)); subsolar temperature formula Tmax in the paper | Wurz 2022 | READ |
| Micrometeoroid impact vaporisation (MIV) | impact speeds about 20 km s-1 (Mercury) and 14 km s-1 (Moon); plume temperature 2500 to 5000 K (use 3000 to 4000 K); vapour volume follows; micrometeoroid flux Mercury (1.66 to 3.56)e-16 up to 2.4e-14 g cm-2 s-1, Moon 5.1e-17 to 8.4e-17 g cm-2 s-1 | Wurz 2022 Table 2 | READ |
| Photon-stimulated desorption (PSD) | cross section Q = 3 +/- 1 e-20 cm2 (Na, 250 to 400 nm) and 2e-20 cm2 (K); surface density N_S = 7.5e14 cm-2; photon flux 2e15 cm-2 s-1; Weibull velocity distribution kappa = 1.7, v0 = 575 m s-1; escape fraction below 1e-3 (Mercury) and below 1e-2 (Moon); photoionisation dominates the loss | Wurz 2022 | READ |
| Ion sputtering | yield 0.01 to 0.1 atoms per ion for H+ and 0.1 to 1.0 for He+ at solar-wind energies; Sigmund analytic yield (Wurz eqs. 29 to 32); Na mineral binding energy 2.6 to 8.4 eV; escaping fraction of sputtered atoms 31% (Fe) to 65% (O) at Mercury and 90% to about 100% at the Moon; CME enhancement up to about 50 times | Wurz 2022 | READ |
| Thermal escape of MIV species at T = 3500 K | for the Moon all MIV species (Na, K, Ca, Mg, Fe and others) escape promptly; for Mercury only species up to O escape immediately | Wurz 2022 Table 3 | READ |
| Thermal escape of volatiles from Mercury and the Moon | Table 1 of the paper lists n_i, v0, X, f_esc and flux for H, H2, He, OH, H2O, Ne, N2, O2, Ar, CO2 (not transcribed here) | Wurz 2022 Table 1 | READ (table not transcribed) |

### 5.2 Sodium and potassium data

| Item | Value | Source | Status |
|---|---|---|---|
| Na photoionisation rate at 1 au (symbol k_NaD2 in my notes of Oza Table 1; symbol to be rechecked) | 5.92e-6 s-1 for G stars, 1.3e-5 for F, 9.5e-7 for K (Huebner and Mukherjee 2015) | Oza 2019 | READ |
| Lifetime 1/rate (G, F, K) | 1.7e5 s (1.95 d), 7.7e4 s (0.89 d), 1.05e6 s (12.2 d), scaling with 1/d^2 | derived | COMPUTED |
| Io Na lifetime | about 4 h | Oza 2019 | 2ND-HAND |
| Na column and Na/K | comets about 8e11 cm-2, Na/K 54; Mercury about 8e11, 80 to 190; Io 5e9 to 1e12, 10; Europa about 2e10, 25; Moon about 5e10, 6; Enceladus at most 1e9, about 100; Earth mesosphere about 8e11, about 100; Jupiter about 7e9; solar Na/K about 16 | Oza 2019 Table 1 | READ (TEST points) |
| Silicate vapour pressure | P_vap = exp(-A/T + B): enstatite A = 68908 +/- 8773 K, B = 38.1 +/- 5.0; fayalite A = 60377 +/- 1082 K, B = 37.7 +/- 0.7 | Oza 2019 | 2ND-HAND |
| Evaporation source rate | Mdot_evap = x_i 4 pi R^2 P_vap sqrt(m_i / (2 pi kB T)); exo-Io Na mass fraction 0.05; T0 of 2000 K or more implies disintegration (Perez-Becker and Chiang 2013) | Oza 2019 | READ, 2ND-HAND (disintegration limit) |
| Mercury Na effective temperatures | about 1100 to 1400 K; Ca exosphere about 12000 to 20000 K; Na e-folding scale height about 200 km at the poles (about 1500 K) | Lierle and Schmidt 2025 | READ |

### 5.3 Moon, Callisto, Europa, exo-airless worlds

| Item | Value | Source | Status |
|---|---|---|---|
| Moon noble-gas loss fractions r_esc / r_ion (solar-wind implanted) | He 0.15 / 1.3e-3; Ar 0 / 7.8e-5; Ne 5e-6 / 8.8e-5; Kr 0 / 1.26e-4; Xe 0 / 1.6e-4 | Lammer 2022 Table 3 | READ |
| Moon model surface densities from that table (cm-3) | He 544; 36Ar+38Ar 374; Ne 4220; Kr 0.19; Xe 0.0182 | Lammer 2022 | READ |
| Moon observations quoted in the same text (cm-3) | He 2000 to 40000; Ar 40000 at sunrise or 8000 (LADEE); Ne about 2000 to 3000 | Lammer 2022 | READ |
| Solar-wind flux used; Moon H2 | 3.2e12 m-2 s-1; H2 about 2100 cm-3 | Lammer 2022 | READ |
| Callisto CO2 column | about 8e14 cm-2 | Vorburger 2019 | READ |
| Callisto O2 column | 3.4e15 (Cunningham) to (3 to 4)e16 (Kliore) cm-2 | Vorburger 2019 | READ |
| Callisto atomic O at exobase; sputtered H2O | up to 1e4 cm-3; 110 cm-3 at the surface (32 cm-3 at 1000 km); an ionosphere lowers densities about 2.5 times | Vorburger 2019 | READ |
| Europa O2 | radiolysis plus thermal release (Johnson 2018, arXiv:1804.10589); column not extracted | Johnson 2018 | NOT REACHED (value) |
| Exo-Mercury and lava-world exospheres | Na, Ca and Ca+ tail searches for CoRoT-7b (Guenther 2011, arXiv:1009.5500; Mura 2011); HD 219134b and 55 Cnc e discussed in Lammer 2022 | Lammer 2022 | READ (discussion), values not extracted |
| Silicate-vapour atmospheres of evaporating super-Earths | Schaefer and Fegley 2009 (arXiv:0906.1204) | | downloaded, NOT extracted |

The generator rule for an airless body: release flux from the dominant source (thermal for X below 15; PSD, MIV and sputtering for refractories), loss by photoionisation with the Na rates above scaled by 1/d^2 and by ballistic return. Typical Na columns span 1e9 to 1e12 cm-2 across the TEST points in the table (Io, Europa, Moon, Mercury, comets), so a log-uniform draw over about three decades is supported by data but the dependence on source strength is not calibrated here (section 8).

## 6. Secondary effects

### 6.1 Flares and energetic particles

| Item | Value | Source | Status |
|---|---|---|---|
| AD Leo flare proton flux (above 10 MeV) | 5.9e8 cm-2 sr-1 s-1 | Segura et al. 2010 (arXiv:1006.0022) | READ |
| Ozone depletion by flare protons for an unmagnetised planet | maximum 94% two years after the flare; the UV flare alone does little | Segura 2010 | READ |
| NOx yield | 1.27 N atoms per ion pair; 1.3 to 1.6 NOx per pair (Rusch 1981) | Segura 2010 | 2ND-HAND |
| Flare energy range of the Kepler GJ 1243 sample used | 10^30.5 to 10^34 erg | Segura 2010 / Tilley | READ |
| Repeated flares | 94% ozone loss within 10 years of repeated proton events, leaving UV-C of about 0.1 to 1 W m-2 at the surface | Tilley et al. 2019 (arXiv:1711.08484) | READ |
| Particle-driven chemistry | N2O and HCN production by solar energetic particles in a weakly reducing atmosphere (formation reactions in 2.6) | Airapetian 2020 | READ |
| Flare frequency scaling | power law N(>E1) / N(>E2) in energy | Johnstone 2021 | READ (exponent not transcribed) |

A magnetic field removes the proton channel in the Segura calculation; the generator should therefore make the proton-driven ozone loss conditional on the magnetosphere flag and flare rate (Part on radiation shielding).

### 6.2 Lightning NOx and N2O

| Item | Value | Source | Status |
|---|---|---|---|
| NO per flash | 15 (range 2 to 40) e25 molecules, equal to 250 mol or 3.5 kg N | Schumann and Huntrieser 2007 (ACP 7, 3823) | READ |
| Global source and flash rate | 5 +/- 3 Tg N yr-1 at 44 +/- 5 flashes s-1 | same | READ |
| Laboratory yields | 1e17 NO per J (corona discharge) and other values 4e12 to 1.5e13 molecules J-1 | same | READ |
| Consistency check | 15e25 x 44 x 3.156e7 s / 6.022e23 x 14 g gives 4.8 Tg N yr-1; per unit area 1.3e9 NO cm-2 s-1 over 5.1e18 cm2 | derived | COMPUTED |

The Earth values are TEST points. For other atmospheres the flash energy and the NOx yield per joule depend on composition (N2/O2 ratio, CO2, H2O); Airapetian 2020 reports SEP-driven HCN exceeding lightning HCN in weakly reducing atmospheres (READ). No composition-dependent lightning yield table was reached (section 8).

### 6.3 Volcanic gases: photochemical lifetimes

Use the lifetimes in section 2.4 divided by an attenuation factor and multiplied by the diurnal factor (2 to 4). Sun-column unattenuated lifetimes: SO2 1.4 h, H2S 1.0 h, OCS 4.0 h, CS2 9 min, HCl 2.3 d, H2SO4 gas 11 h, NH3 1.9 h, HCN 20 h. Observed tropospheric lifetimes are much longer because of shielding, OH reaction rates and deposition: Earth tropospheric H2S and SO2 about 2 days (Hu 2013, READ; TEST). SO2 + OH + M gives k = 6.6e-13 cm3 s-1 at 220 K and 1.2e-12 at 298 K (2.1), so an OH abundance of 1e6 cm-3 gives an OH-limited lifetime of 10 to 18 days (COMPUTED); around active M dwarfs OH is about 1 cm-3 or less (Segura 2005, READ) so this channel is shut and photolysis and deposition dominate. For a haze-forming sulfur flux see Hu 2013 thresholds in 2.6.

## 7. Generator recipe (formulas and validity)

1. Photochemical lifetime of species i against the host star: tau_i = 1 / (J_i,table x S x f_geom x f_shield). J_i,table from section 2.4 (S = 1); S in units of the stellar constant; f_geom 0.25 to 0.5 (diurnal or global mean); f_shield from the column (tau = 1 mixing ratios in 2.5, critical pressure in 2.5). Valid for the six host spectra tabulated; interpolate between hosts by T_eff. Not valid for F stars below 200 nm.
2. Ozone: Kozakis fit from the project doc, with the draw between the 1-D fit and the 3-D Cooke value (3.2); host-dependence direction from 3.4.
3. Escape: Jeans parameter with the project convention; energy-limited form with eta 0.05 to 0.2 or the Owen and Wu variable law (4.2); the hydro fit only inside its validity box (4.2); the diffusion limit as a cap for hydrogen (4.3); XUV history from Ribas for G stars, Jackson for other types, Tu for rotation spread (4.5).
4. Exospheres: release by thermal, PSD, MIV and sputtering (5.1), loss by photoionisation scaled as 1/d^2 (5.2).
5. Perturbations: flare protons conditional on the magnetosphere (6.1), lightning NOx scaling by flash density (6.2).

## 8. Boss must research

| Missing item | Why it matters | Where to get it |
|---|---|---|
| JPL 19-5 rate evaluation (ClO + O, ClO + NO, ClO + ClO, ClONO2 formation; Troe F_c and N for all falloff reactions; HO2 + O3) | closes the ClOx cycle and the falloff form | jpldataeval.jpl.nasa.gov (unreachable) or NIST Kinetics Database queries |
| Phosphorus chemistry rates and PH3 cross sections beyond Sousa-Silva 2020 | no PH3 lifetime by star | Bains 2021 supplement; Chen 1991 |
| F-star FUV and EUV spectra | the F column in 2.3 and 2.4 is blank below 200 nm | MUSCLES-type survey of F stars, or a stellar model library via MAST |
| Zerkle 2012, Trainer 2006 originals (methane haze threshold) | the CH4/CO2 threshold is read only through Arney 2016 | original journal papers |
| Harman 2015 (CO2 photolysis O2 and O3) | abiotic O2 numbers are 2ND-HAND | Astrophys. J. 812, 137 |
| Zahnle 2006 (sulfur) and Kasting and Donahue 1980 | sulfur and hydrogen-escape cycle details | original papers |
| Hunten 1973 and Hunten 1987 originals | the diffusion-limit form was read through Gronoff 2020 | original papers |
| Boil-off and hydro escape rates for non-hydrogen envelopes | the Kubyshkina fit is hydrogen-only | Owen and Wu 2016; later hydrodynamic grids |
| Owen and Jackson 2012 efficiencies; Lopez and Fortney tables | eta spread is 2ND-HAND | arXiv:1206.2367; 1305.0269 tables |
| Wurz 2022 Tables 1 and 3 values; Mercury and Moon Ca, Mg, K column densities | only formulas transcribed | Space Sci. Rev. 218, 10; Mura 2011 |
| Io and Europa original exosphere data; Schaefer and Fegley 2009 table values | silicate-vapour exospheres not quantified | original papers; arXiv:0906.1204 |
| Composition-dependent lightning NOx yield and N2O yield from SEPs | lightning scaling for non-Earth air | Hodosan 2016 style models; Airapetian 2016 |
| Troe or ion chemistry for N2/CH4 haze and ion escape rate formulas | nonthermal escape and haze at Titan-like conditions | Lammer and Gronoff reviews, individual papers |
| A temperature and pressure dependent UV cross-section set for anoxic shielding (SO2, S8, H2SO4) beyond band means | shield estimates are band means | Mainz UV/VIS atlas (reachable, bulk download not done) |

## 9. Recalled / to verify

| Item | Why flagged |
|---|---|
| Chapman reaction scheme (four reactions in 3.1) | textbook, not read in a source in this session |
| Earth ozone-layer peak near 20 to 25 km | standard figure, not read |
| Lindemann form k0 [M] / (1 + k0 [M] / kinf) for the falloff entries in `zahnle_earth.yaml` | the file shows no Troe keys; the code that reads it was not inspected |
| Which of the two Johnstone 2021 slopes belongs to the saturated branch | the shallow slope is probably the saturated branch; the read summary may have swapped the labels |
| 300 DU upper edge of the modern ozone test range | the 279 DU value is read; the upper edge comes from the task brief |
| Symbol of the Na photoionisation rate in Oza 2019 Table 1 (k_NaD2) | rate versus lifetime wording to be rechecked |
| Hu 2012 definitions of F and mu0 in eq. D5 | the formula was read, the variable definitions were not copied here |
| Standard limiting form Phi = b f / H_a for light trace gases | derived form, not read in that exact form |
| The statement that the shipped mechanism's O + O2 + M value matches JPL | the JPL value is recalled, not reached |

## 10. Cross-links

Part 03 (infrared greenhouse and thermal) for lapse rates, temperature profiles and the exobase temperature inputs; Part 05 (mixing ratios by planet class) for the bulk compositions that set the photochemical regime; Part 06 (condensates, clouds, hazes, dust) for haze production rates and S8 and H2SO4 aerosols; the gas-optics part for cross sections and UV opacity beyond the band means used here; the radiation-shielding part for magnetosphere and particle-flux inputs; the gas-table and thermodynamics part for molar masses, diffusion data and species lists. Project docs: `atmospheres-retention-and-classes.md` (shoreline, Jeans convention, `retention()`), `atmosphere-sources-check.md` (Kozakis, Cooke, Segura, Rugheimer ozone).

## 11. Row counts by status label

Counted from the Status columns by a script over the tables above (rows with two labels count once under each label).

Rows with a Status column: 120. READ 52, 2ND-HAND 11, DATA 49, COMPUTED 28, NOT REACHED 5, SNIPPET 0, RECALLED 0 (recalled items are in section 9 and flagged inline). The large lifetime and J-value tables in 2.3 and 2.4 and the Kubyshkina and Ribas tables carry their label in the heading text: COMPUTED from DATA inputs (2.3, 2.4), READ coefficients with COMPUTED evaluations (4.2, 4.5).
