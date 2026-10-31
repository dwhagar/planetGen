# Part 02: Gas optics from the extreme UV to the near IR

Research lane for the atmosphere-composition work (v9.0). Read-only research; nothing in the repo was edited. Solar System bodies appear only as labelled TEST points. Status labels follow `_worker-rules.md`: READ, 2ND-HAND, DATA, COMPUTED, SNIPPET, RECALLED, NOT REACHED. A row may carry two labels (for example DATA then COMPUTED) when a data file was read and a number was derived from it.

This part extends, and does not redo, sections 2 to 4 of `/mnt/project-files/reports/atmosphere-sources-check-2026-10-17.md` (verified Rayleigh law, per-gas refractivity of N2, O2, Ar, CO2, CO, CH4, H2O, NH3, H2, He; Earth test tau_R(550 nm) = 0.0969; ozone fit of Kozakis 2022) and `docs/design/composition/atmosphere-colour.md`.

Reachability record: worked were export.arxiv.org PDFs, PyPI packages (photochem_clima_data, musica/TUV-x, colour-science), the Mainz UV/VIS atlas data files (under `/uvvis_data/`), the IUP Bremen ozone page, physics.nist.gov ASD, archive.stsci.edu (BT-Settl), sbuv.gsfc.nasa.gov (solar spectrum), Wikipedia and Optica abstract pages (citation only). Failed or paywalled: ScienceDirect, Cloudflare and Radware challenge pages, archive.org lending texts (HTTP 401), the NTRS search for Kasting and Donahue.

---

## 1. Scope and how the numbers feed the generator

| Generator input | Output of this part | Section |
|---|---|---|
| Gas columns from the stored mixture (part 05), pressure, gravity | tau_R(lambda) of the mixture from per-gas sigma, incl. Ne, Kr, Xe, N2O, NH3, H2S, SO2, C2H6, C2H4, HCN, O3 (GEN.417) | 2.1, 3.1 |
| Columns of O3, NO2, SO2 and other absorbers | absorption optical depth in coarse bins, for colour and UV | 2.3 to 2.5 |
| Star spectrum (Planck, BT-Settl) | chromaticity, absolute sRGB, UV at the top of the air | 2.9, 3.6 |
| O2 and O3 columns, zenith, stellar UV | UV index, erythemal and DNA doses, UV-C at the ground | 2.6, 2.7, 3.4 |
| Cross sections times quantum yield times flux | photolysis rates j (GEN.413 and the photochemistry lane) | 2.3, 3.3 |
| Hot-giant temperature class | alkali, TiO, VO, FeH, CrH opacity regime | 2.5, 3.7 |
| tau_R, aerosol tau, albedo, zenith | single-scatter sky, two-stream, Chapman air mass | 3.5 |

New compared with the earlier report: a generic Rayleigh formula (3.1) for Ne, Kr, Xe, N2O, NH3, H2S, SO2, C2H6, C2H4, HCN, O3 and C2H2, OCS, HCl, H2CO within -9 to +6 percent of 13 verified gases; DFT-computed depolarization ratios (2.2); real spreads for H2 (0.202 to 0.217 of N2), CH4 (30 percent below the older coded formula) and N2 depolarization (0.020 to 0.031); 22 gases in 35 coarse absorption bins (2.3) with temperature trends (2.4); a computed surface-UV demonstration (2.7); and the Source re-check (2.10).

Rules kept: nothing is fixed per class; Solar System bodies appear only as test points (Earth tau_R, Earth ozone columns, UV index of Earth).

---

## 2. Tables

### 2.1 Rayleigh scattering cross sections, extension of the verified table

Definitions: sigma in cm2 per molecule; F_K = (6 + 3 rho)/(6 - 7 rho) with rho the natural-light depolarization ratio; ratio is sigma(550 nm)/sigma_N2(550 nm) = 4.635e-27. "n-based" means the Sneep and Ubachs formula of the earlier report with the refractivity at its own reference density (288.15 K or 273.15 K, see report section 2.1).

**2.1a Gases whose cross section comes from a measured refractivity (n-based), including the new noble gases, N2O and the re-checks.**

| Gas | sigma 450 nm | sigma 550 nm | sigma 650 nm | Ratio to N2 | F_K (550) | Source | Status |
|---|---|---|---|---|---|---|---|
| N2 | 1.054e-26 | 4.635e-27 | 2.351e-27 | 1.000 | 1.035 | report 2.2 formula, re-evaluated | COMPUTED |
| O2 | 9.32e-27 | 4.072e-27 | 2.062e-27 | 0.879 | 1.102 | report 2.2, re-evaluated; He 2021 Table 2 measured 4.55e-27 at 532 nm against n-based 4.642e-27 | COMPUTED, READ |
| Ar | 9.07e-27 | 3.991e-27 | 2.025e-27 | 0.861 | 1.000 | report 2.2, re-evaluated | COMPUTED |
| CO2 | 2.663e-26 | 1.164e-26 | 5.88e-27 | 2.51 | 1.145 | report 2.2; He 2021 Table 2: 13.32 n-based and 13.3 measured (1e-27 cm2) at 532 nm, so the 550 nm value 1.164e-26 is consistent | COMPUTED, READ |
| CH4 (He 2021 dispersion) | 2.270e-26 | 9.905e-27 | 5.00e-27 | 2.14 | 1.000 | He et al. 2021, ACP 21, 14927, Table 1: (n-1)e8 = 4869.8 + 4.1023e14/(1.133e10 - nu^2) at 288.15 K; Table 2: 11.34 n-based and 11.3 measured (1e-27 cm2) at 532 nm | READ, COMPUTED |
| CH4 (older coded form, n-1 = 4.6662e-4 + 4.02e-6/lam^2) | 2.942e-26 | 1.283e-26 | 6.47e-27 | 2.77 | 1.000 | report 2.2; **30 percent above the measured value; replace it with the He 2021 form** | COMPUTED |
| CO (Allen/Cox form) | 1.313e-26 | 5.734e-27 | 2.896e-27 | 1.24 | 1.016 | report 2.2 (2ND-HAND form, TauREx); no measurement reached | 2ND-HAND, COMPUTED |
| H2, route A: Peck and Huang n-based | 2.31e-27 | 1.008e-27 | 5.09e-28 | 0.217 | 1.038 | (n-1) = 0.0148956/(180.7 - s^2) + 0.0049037/(92 - s^2), s = 1/lam in um^-1, at 273.15 K, valid 0.168 to 1.6945 um; refractiveindex.info `H2_nk_Peck.yml` (CC0) citing Peck and Huang, JOSA 67, 1550 (1977) | DATA, COMPUTED |
| H2, route B: Dalgarno and Williams direct | 2.15e-27 | 9.38e-28 | 4.74e-28 | 0.202 | n/a | report 2.2 (2ND-HAND); evaluated here | 2ND-HAND, COMPUTED |
| He | 1.375e-28 | 6.11e-29 | 3.12e-29 | 0.0132 | 1.000 | n-based from Mansfield and Peck (refractiveindex.info `He_nk_Mansfield.yml`); the report formula gives 6.04e-29 (ratio 0.0130) | DATA, COMPUTED |
| Ne | 4.93e-28 | 2.187e-28 | 1.115e-28 | 0.0472 | 1.000 | Bideau-Mehu et al. fit, `Ne_nk_Bideau-Mehu.yml` (fitted 0.1404 to 0.5462 um, 650 nm is extrapolated); Cuthbertson gives 2.263e-28, Borzsonyi 2.260e-28. Spread 3.5 percent | DATA, COMPUTED |
| Kr | 2.106e-26 | 9.208e-27 | 4.656e-27 | 1.99 | 1.000 | Bideau-Mehu `Kr_nk_Bideau-Mehu.yml` (fitted to 0.6234 um); Cuthbertson 9.195e-27, Borzsonyi 9.271e-27. Spread 0.8 percent | DATA, COMPUTED |
| Xe | 5.567e-26 | 2.402e-26 | 1.205e-26 | 5.18 | 1.000 | Bideau-Mehu `Xe_nk_Bideau-Mehu.yml`; Cuthbertson 2.489e-26, Borzsonyi 2.498e-26 (the database notes a typo in the paper's coefficient, 12.75e-6 should be 12.75e-3). Spread 4.0 percent | DATA, COMPUTED |
| N2O | 3.664e-26 | 1.587e-26 | 7.978e-27 | 3.42 | 1.2237 | He et al. 2021 Table 1: (n-1)e8 = 46890 + 4.12e-6 nu^2 at 288.15 K, valid 15453 to 21838 cm^-1 (458 to 647 nm); F_K = (3.3462 + 70.8e-12 nu^2)/(2.7692 - 47.2e-12 nu^2); Table 2: 18.19 n-based and 18.3 measured (1e-27 cm2) at 532 nm | READ, COMPUTED |
| NH3, Cuthbertson | 1.661e-26 | 7.168e-27 | 3.597e-27 | 1.55 | 1.00 (assumed) | (n-1) = 0.032953/(90.392 - s^2) at 273.15 K, fitted 0.48 to 0.67 um; `NH3_nk_Cuthbertson.yml` | DATA, COMPUTED |
| NH3, TauREx/Allen | 1.716e-26 | 7.408e-27 | 3.716e-27 | 1.60 | 1.00 (assumed) | (n-1) = 37.0e-5 (1 + 12.0e-3/lam^2) at 273.15 K (report 2.2) | 2ND-HAND, COMPUTED |
| H2O vapour (generic check only) | 7.72e-27 | 3.38e-27 | 1.71e-27 | 0.73 | 1.0010 | report gives n - 1 = 2.362e-4 at 550 nm with unknown reference temperature: n-based sigma is 3.11e-27 (if 288.15 K) or 2.80e-27 (if 273.15 K); the generic formula with alpha = 1.501 A^3 gives 3.38e-27, 9 to 21 percent above. Use the report value | COMPUTED |

Findings (COMPUTED unless stated):

- H2/N2 is 0.202 by Dalgarno and Williams but 0.217 by the Peck and Huang refractivity with F_K 1.038 (0.215 with the DFT rho); both give n - 1 = 1.37e-4 at 1 um, so the difference is dispersion shape. Draw 0.20 to 0.22.
- N2 depolarization: F_K 1.034 (rho 0.0199, report), rho 0.0305 (Penndorf 1957 via Ranjan and Sasselov, 2ND-HAND, F_K 1.057), rho 0.0221 (my DFT); a 2.2 percent spread in sigma. Ranjan and Sasselov (arXiv:1610.06223, READ, appendix A) use rho(O2) = 0.054 and rho(CO2) = 0.0774, consistent with the report (0.057, 0.0747).
- He 2021 Table 2 (BBCES, READ) confirms CO2, N2O, CH4 and O2 n-based values to within about 2 percent at 330 to 710 nm. Ne, Kr, Xe, CO and NH3 have no direct measurement reached.

**2.1b Polarizability route for the gases without a reachable refractivity.** alpha is the CCCBDB orientation-averaged static polarizability (NIST CCCBDB table, public domain, file `web/pol.json`); rho is from my DFT anisotropy (2.2); sigma uses section 3.1 with B = 0.0075 um^2; the "B range" column shows sigma(550) at B = 0.0045 and 0.0105.

| Gas | alpha (A^3) | rho_n | F_K | sigma 450 | sigma 550 | sigma 650 | Ratio to N2 | B range at 550 nm | Status |
|---|---|---|---|---|---|---|---|---|---|
| H2S | 3.631 | 0.00004 | 1.000 | 4.52e-26 | 1.98e-26 | 9.99e-27 | 4.26 | 1.94e-26 to 2.01e-26 | DATA, COMPUTED |
| SO2 | 3.882 | 0.0388 | 1.068 | 5.51e-26 | 2.41e-26 | 1.22e-26 | 5.20 | 2.37e-26 to 2.46e-26 | DATA, COMPUTED |
| C2H6 | 4.226 | 0.0025 | 1.004 | 6.14e-26 | 2.69e-26 | 1.36e-26 | 5.80 | 2.64e-26 to 2.74e-26 | DATA, COMPUTED |
| C2H4 | 4.188 | 0.0245 | 1.042 | 6.26e-26 | 2.74e-26 | 1.39e-26 | 5.91 | 2.69e-26 to 2.79e-26 | DATA, COMPUTED |
| HCN | 2.593 | 0.0331 | 1.058 | 2.44e-26 | 1.07e-26 | 5.39e-27 | 2.30 | 1.05e-26 to 1.09e-26 | DATA, COMPUTED |
| C2H2 | 3.487 | 0.0356 | 1.062 | 4.42e-26 | 1.94e-26 | 9.78e-27 | 4.17 | 1.90e-26 to 1.97e-26 | DATA, COMPUTED |
| OCS | 5.090 | 0.0675 | 1.122 | 9.96e-26 | 4.36e-26 | 2.20e-26 | 9.40 | 4.27e-26 to 4.44e-26 | DATA, COMPUTED |
| HCl | 2.515 | 0.0016 | 1.003 | 2.17e-26 | 9.50e-27 | 4.81e-27 | 2.05 | 9.32e-27 to 9.69e-27 | DATA, COMPUTED |
| H2CO | 2.770 | 0.0275 | 1.047 | 2.75e-26 | 1.20e-26 | 6.09e-27 | 2.60 | 1.18e-26 to 1.23e-26 | DATA, COMPUTED |
| NH3 (generic) | 2.103 | 0.0016 | 1.003 | 1.52e-26 | 6.65e-27 | 3.36e-27 | 1.43 | 6.52e-27 to 6.77e-27 | DATA, COMPUTED |
| O3, CCCBDB alpha with DFT rho | 3.079 | 0.1097 | 1.210 | 3.93e-26 | 1.72e-26 | 8.69e-27 | 3.71 | 1.69e-26 to 1.75e-26 | DATA, COMPUTED |
| O3, Ranjan and Sasselov values | 3.21 | (KCF 1.06) | 1.060 | 3.74e-26 | 1.64e-26 | 8.28e-27 | 3.53 | n/a | 2ND-HAND, COMPUTED |
| Ne, Kr, Xe via generic (check of the formula) | 0.381, 2.498, 4.005 | 0 | 1.000 | 4.97e-28, 2.14e-26, 5.49e-26 | 2.18e-28, 9.35e-27, 2.40e-26 | 1.10e-28, 4.73e-27, 1.22e-26 | 0.047, 2.02, 5.19 | n/a | DATA, COMPUTED |

Notes on 2.1b:

- O3: Ranjan and Sasselov (READ, appendix A.8) use KCF = 1.06 and alpha = 3.21 A^3; my DFT anisotropy gives F_K 1.21 (a multi-reference molecule, where PBE0 is least trustworthy). The two O3 results differ by 5 percent in sigma; use 1.5e-26 to 1.75e-26 at 550 nm; O3 Rayleigh is minor because O3 is a trace gas.
- For SO2 and H2S Ranjan and Sasselov (appendix A.5, A.6, READ) use the same formula with alpha and delta from Bogaard et al. 1978 at 488, 514.5 and 632.8 nm but print no numbers, so Bogaard stays NOT REACHED; my values are independent of it.
- The generic formula errs by -9 to +6 percent on the 13 verified gases (3.1), so each value here carries about +-8 percent. CH4, CO2 and NH3 come out 7 to 9 percent low, so similar molecules (C2H6, HCl, H2S) may too.

### 2.2 Depolarization ratios and polarizability anisotropy (DFT-computed)

Method (COMPUTED, my own calculation): PySCF finite field, PBE0, geometry optimized with geomeTRIC at aug-cc-pVDZ, polarizability at aug-cc-pVTZ by central difference with field 0.002 a.u.; static, not frequency dependent. Mean alpha = (a1 + a2 + a3)/3; anisotropy gamma^2 = 1/2 sum (ai - aj)^2; natural-light rho_n = 6 gamma^2/(45 alpha^2 + 7 gamma^2) (relation checked by Monte Carlo, 2.8). Starting geometries were approximate (RECALLED) and optimized away.

| Gas | alpha DFT (A^3) | alpha CCCBDB (A^3) | rho_n DFT | rho_n or F_K literature in this file | DFT versus literature | Status |
|---|---|---|---|---|---|---|
| N2 | 1.762 | 1.710 | 0.0221 | F_K 1.034 gives 0.0199 (report); 0.0305 (Penndorf via Ranjan) | +11 percent, -28 percent | COMPUTED |
| H2 | 0.840 | 0.787 | 0.0158 | 0.0221 (assumed in table 2.1a) | -29 percent | COMPUTED |
| CO | 1.941 | 1.953 | 0.0093 | F_K 1.016 gives 0.0095 | -2 percent | COMPUTED |
| CO2 | 2.524 | 2.507 | 0.0732 | F_K 1.1364 gives 0.0747; Shemansky 0.0774 (Ranjan) | -2 percent | COMPUTED |
| N2O | 2.860 | 2.998 | 0.1107 | He 2021 F_K 1.2237 gives 0.116 | -5 percent | COMPUTED |
| H2O | 1.422 | 1.501 | 0.00083 | 0.000299 (Murphy 1977 via report) | 2.8 times larger, both tiny | COMPUTED |
| CH4 | 2.534 | 2.448 | 0 | 0.0002 assumed by Ranjan (as CCl4) | consistent | COMPUTED |
| NH3 | 2.122 | 2.103 | 0.0016 | none reached | unverified | COMPUTED |
| H2S | 3.677 | 3.631 | 0.00004 | none reached (Bogaard NOT REACHED) | unverified | COMPUTED |
| SO2 | 3.874 | 3.882 | 0.0388 | none reached | unverified | COMPUTED |
| HCN | 2.536 | 2.593 | 0.0331 | none reached | unverified | COMPUTED |
| C2H2 | 3.495 | 3.487 | 0.0356 | none reached | unverified | COMPUTED |
| C2H4 | 4.157 | 4.188 | 0.0245 | none reached | unverified | COMPUTED |
| C2H6 | 4.319 | 4.226 | 0.0025 | none reached | unverified | COMPUTED |
| O3 | 2.702 | 3.079 | 0.1097 | KCF 1.06 gives 0.035 (Ranjan) | factor 3, unreliable | COMPUTED |
| OCS | 5.044 | 5.090 | 0.0675 | none reached | unverified | COMPUTED |
| HCl | 2.565 | 2.515 | 0.0016 | none reached | unverified | COMPUTED |
| H2CO | 2.660 | 2.770 | 0.0275 | none reached | unverified | COMPUTED |

Validation reading: for gases whose F_K is known from n-based work the DFT F_K is within 1.1 percent of the n-based value for N2, CO2, N2O, CO and H2, so a wrong rho changes sigma by at most about 1 percent unless rho exceeds 0.05. Only O3 is unreliable. The alpha(DFT)/alpha(CCCBDB) ratio is 0.88 to 1.07 (median 1.00).

### 2.3 UV and visible absorption cross sections as coarse band tables (embeddable)

Source of every number in the three tables: `photochem_clima_data` 0.3.4 (PyPI; GNU GPL v3 for the package, while each dataset keeps the licence of its origin, which I did not find stated per file; the data origins are Heays et al. 2017 (arXiv:1701.04459; commonly called the Leiden VUV database, RECALLED name), Huebner 2015, Schmidt 2013 for CO2, Karkoschka 1994 for CH4, Keller-Rudek 2013 for NH3 and others as listed below), file `data/xsections/<species>.h5`, dataset `photoabsorption` (total photoabsorption, cm2 per molecule). I averaged the file onto 35 bins (COMPUTED): each cell is the mean over the part of the bin the data cover; `-` means the bin lies outside the data (the package treats cross sections as zero there); `*` means under 50 percent of the bin is covered; `0` means the file holds zero. Bins below 100 nm are photoabsorption in the X-ray and EUV (smooth, 1e-18 to 1e-16 cm2). The file also has `photodissociation` and `photoionisation` datasets; use photodissociation times the quantum yield (also in the file) for a photolysis rate j.

Provenance (DATA; `metadata.yaml`): coverage in nm is O2 4.2 to 239.9; O3 0.06 to 1099.9 (Huebner 2015 below 52.7 nm, Heays 2017 above); SO2 1.7 to 416.6; H2O 0.1 to 230.4 (Ranjan 2020 for 192 to 230 nm); CO2 0.2 to 212.8 (Schmidt 2013 at 270 K, reduced resolution, above 138.9 nm); CO 0.1 to 166.5; CH4 5.7 to 1000 (**Karkoschka 1994 for 152 to 1000 nm, so all visible CH4 bands rest on a paper that was NOT REACHED**); NH3 0.1 to 230; H2S 0.06 to 316; N2O 91.3 to 209.9; NO2 0.06 to 660; HCN 0.1 to 187.9; C2H2 6.3 to 299.4; C2H6 15.4 to 149.5; H2CO 6.3 to 375.9; HCl 4.5 to 229.9; OCS 3.5 to 279.4; H2O2 65.1 to 349.9; N2 11.6 to 99.9 (Ranjan and Sasselov, READ, quote no N2 extinction detected longward of 108 nm); NO 0.1 to 229.3; C2H4 0.1 to 194.9; H2 0.1 to 119.9. Unlisted origins are Heays 2017 or Huebner 2015.

Cross sections in cm2 per molecule (COMPUTED from DATA):

| bin nm | O2 | O3 | SO2 | H2O | CO2 | CO | CH4 |
|---|---|---|---|---|---|---|---|
| 10-30 | 9.5e-18 | 1.2e-17 | 9.8e-18 | 5.0e-18 | 1.1e-17 | 6.5e-18 | 2.3e-18 |
| 30-60 | 2.2e-17 | 2.6e-17 | 3.5e-17 | 1.7e-17 | 2.9e-17 | 1.8e-17 | 1.8e-17 |
| 60-80 | 2.5e-17 | 3.4e-17 | 5.7e-17 | 2.1e-17 | 3.5e-17 | 2.4e-17 | 3.9e-17 |
| 80-91.2 | 1.5e-17 | 2.0e-17 | 4.3e-17 | 1.9e-17 | 3.5e-17 | 3.3e-17 | 4.7e-17 |
| 91.2-100 | 1.2e-17 | 1.1e-17 | 4.6e-17 | 1.8e-17 | 5.8e-17 | 3.0e-17 | 4.4e-17 |
| 100-110 | 1.5e-18 | 9.3e-18 | 2.4e-17 | 8.1e-18 | 2.8e-17 | 2.0e-17 | 3.0e-17 |
| 110-120 | 6.2e-19 | 9.8e-18 | 1.9e-17 | 6.2e-18 | 1.7e-17 | 1.0e-18 | 1.9e-17 |
| 120-123 | 2.9e-18 | 2.4e-17 | 3.6e-17 | 1.0e-17 | 6.4e-20 | 7.0e-20 | 1.8e-17 |
| 123-130 | 1.7e-18 | 8.0e-18 | 6.1e-17 | 8.1e-18 | 3.3e-19 | 1.5e-19 | 1.8e-17 |
| 130-140 | 6.8e-18 | 1.1e-17 | 6.4e-18 | 4.0e-18 | 6.4e-19 | 4.8e-18 | 7.9e-18 |
| 140-150 | 1.4e-17 | 3.4e-18 | 4.5e-18 | 7.6e-19 | 3.7e-19 | 1.9e-17 | 9.5e-20 |
| 150-160 | 8.3e-18 | 1.7e-18 | 3.5e-18 | 2.6e-18 | 3.1e-19 | 1.0e-17 | 2.9e-24 |
| 160-170 | 2.8e-18 | 1.2e-18 | 5.1e-19 | 5.3e-18 | 6.7e-20 | 5.0e-20 | 0 |
| 170-180 | 5.2e-19 | 8.7e-19 | 1.1e-18 | 3.5e-18 | 6.9e-21 |  -  | 0 |
| 180-190 | 1.7e-20 | 6.6e-19 | 3.8e-18 | 2.0e-19 | 4.1e-22 |  -  | 0 |
| 190-200 | 3.9e-22 | 4.1e-19 | 7.3e-18 | 1.7e-21 | 1.7e-23 |  -  | 0 |
| 200-210 | 7.6e-24 | 3.8e-19 | 5.9e-18 | 1.7e-23 | 4.8e-25 |  -  | 0 |
| 210-220 | 5.1e-24 | 1.1e-18 | 2.6e-18 | 3.3e-25 | 3.8e-26* |  -  | 0 |
| 220-230 | 3.2e-24 | 3.0e-18 | 4.0e-19 | 5.0e-27 |  -  |  -  | 0 |
| 230-240 | 1.6e-24 | 6.4e-18 | 3.9e-20 | 2.9e-28* |  -  |  -  | 0 |
| 240-260 |  -  | 1.1e-17 | 1.1e-19 |  -  |  -  |  -  | 0 |
| 260-280 |  -  | 7.7e-18 | 5.3e-19 |  -  |  -  |  -  | 0 |
| 280-300 |  -  | 1.6e-18 | 7.5e-19 |  -  |  -  |  -  | 0 |
| 300-320 |  -  | 1.2e-19 | 2.9e-19 |  -  |  -  |  -  | 0 |
| 320-340 |  -  | 6.2e-21 | 7.3e-21 |  -  |  -  |  -  | 0 |
| 340-360 |  -  | 2.1e-22 | 1.1e-22 |  -  |  -  |  -  | 0 |
| 360-400 |  -  | 2.9e-24 | 7.3e-23 |  -  |  -  |  -  | 0 |
| 400-450 |  -  | 6.2e-23 | 4.7e-25* |  -  |  -  |  -  | 2.1e-28 |
| 450-500 |  -  | 5.4e-22 |  -  |  -  |  -  |  -  | 1.4e-27 |
| 500-550 |  -  | 2.2e-21 |  -  |  -  |  -  |  -  | 8.3e-27 |
| 550-600 |  -  | 4.3e-21 |  -  |  -  |  -  |  -  | 5.1e-27 |
| 600-650 |  -  | 3.8e-21 |  -  |  -  |  -  |  -  | 4.4e-26 |
| 650-700 |  -  | 1.5e-21 |  -  |  -  |  -  |  -  | 2.7e-26 |
| 700-800 |  -  | 3.9e-22 |  -  |  -  |  -  |  -  | 2.4e-25 |
| 800-1100 |  -  | 3.8e-23 |  -  |  -  |  -  |  -  | 1.6e-24 |

| bin nm | NH3 | H2S | N2O | NO2 | HCN | C2H2 | C2H6 |
|---|---|---|---|---|---|---|---|
| 10-30 | 3.6e-18 | 1.1e-18 |  -  | 1.1e-17 | 3.9e-18 | 3.9e-18 | 3.8e-18 |
| 30-60 | 1.8e-17 | 5.5e-18 |  -  | 2.8e-17 | 1.4e-17 | 2.1e-17 | 2.7e-17 |
| 60-80 | 3.1e-17 | 1.6e-17 |  -  | 4.0e-17 | 2.7e-17 | 4.4e-17 | 6.1e-17 |
| 80-91.2 | 2.6e-17 | 2.5e-17 |  -  | 4.4e-17 | 2.7e-17 | 4.6e-17 | 6.7e-17 |
| 91.2-100 | 1.8e-17 | 3.3e-17 | 4.8e-17 | 3.0e-17 | 2.6e-17 | 3.9e-17 | 5.3e-17 |
| 100-110 | 2.4e-17 | 4.2e-17 | 2.3e-17 | 1.4e-17 | 3.1e-17 | 3.0e-17 | 4.4e-17 |
| 110-120 | 1.9e-17 | 3.9e-17 | 1.5e-17 | 9.7e-18 | 4.0e-17 | 2.2e-17 | 2.9e-17 |
| 120-123 | 1.2e-17 | 4.3e-17 | 3.8e-18 | 9.8e-18 | 2.5e-17 | 3.2e-17 | 2.3e-17 |
| 123-130 | 1.1e-17 | 4.4e-17 | 5.1e-17 | 2.4e-17 | 3.6e-18 | 2.4e-17 | 2.4e-17 |
| 130-140 | 1.1e-17 | 3.3e-17 | 1.8e-17 | 9.8e-18 | 4.2e-18 | 3.8e-17 | 2.0e-17 |
| 140-150 | 2.3e-18 | 4.0e-17 | 3.8e-18 | 1.3e-17 | 2.2e-18 | 9.2e-18 | 3.1e-18 |
| 150-160 | 1.1e-18 | 1.6e-17 | 4.1e-19 | 1.5e-17 | 3.5e-19 | 1.8e-17 |  -  |
| 160-170 | 9.9e-19 | 2.8e-19 | 5.0e-20 | 1.6e-17 | 1.4e-19 | 9.3e-19 |  -  |
| 170-180 | 3.2e-18 | 1.2e-18 | 8.0e-20 | 1.1e-17 | 3.9e-20 | 1.0e-18 |  -  |
| 180-190 | 8.0e-18 | 4.0e-18 | 1.0e-19 | 1.8e-18 | 3.7e-21 | 3.4e-19 |  -  |
| 190-200 | 9.8e-18 | 6.7e-18 | 6.5e-20 | 2.5e-19 |  -  | 9.7e-20 |  -  |
| 200-210 | 5.3e-18 | 4.6e-18 | 1.5e-20 | 3.4e-19 |  -  | 2.9e-20 |  -  |
| 210-220 | 6.4e-19 | 2.3e-18 |  -  | 4.3e-19 |  -  | 7.7e-21 |  -  |
| 220-230 | 2.0e-20 | 1.0e-18 |  -  | 3.7e-19 |  -  | 9.4e-21 |  -  |
| 230-240 |  -  | 4.0e-19 |  -  | 1.8e-19 |  -  | 5.1e-21 |  -  |
| 240-260 |  -  | 7.0e-20 |  -  | 2.6e-20 |  -  | 3.6e-21 |  -  |
| 260-280 |  -  | 1.3e-21 |  -  | 3.1e-20 |  -  | 2.3e-21 |  -  |
| 280-300 |  -  | 3.1e-22 |  -  | 8.6e-20 |  -  | 6.4e-22 |  -  |
| 300-320 |  -  | 9.4e-23 |  -  | 1.8e-19 |  -  |  -  |  -  |
| 320-340 |  -  |  -  |  -  | 3.1e-19 |  -  |  -  |  -  |
| 340-360 |  -  |  -  |  -  | 4.4e-19 |  -  |  -  |  -  |
| 360-400 |  -  |  -  |  -  | 5.6e-19 |  -  |  -  |  -  |
| 400-450 |  -  |  -  |  -  | 5.6e-19 |  -  |  -  |  -  |
| 450-500 |  -  |  -  |  -  | 3.5e-19 |  -  |  -  |  -  |
| 500-550 |  -  |  -  |  -  | 1.6e-19 |  -  |  -  |  -  |
| 550-600 |  -  |  -  |  -  | 6.2e-20 |  -  |  -  |  -  |
| 600-650 |  -  |  -  |  -  | 2.1e-20 |  -  |  -  |  -  |
| 650-700 |  -  |  -  |  -  | 9.2e-21* |  -  |  -  |  -  |
| 700-800 |  -  |  -  |  -  |  -  |  -  |  -  |  -  |
| 800-1100 |  -  |  -  |  -  |  -  |  -  |  -  |  -  |

| bin nm | H2CO | HCl | OCS | H2O2 | N2 | NO | C2H4 | H2 |
|---|---|---|---|---|---|---|---|---|
| 10-30 | 6.5e-18 | 1.1e-18 | 7.3e-18 |  -  | 6.4e-18 | 7.0e-18 | 5.5e-18 | 3.5e-19 |
| 30-60 | 2.5e-17 | 1.8e-17 | 2.6e-17 |  -  | 2.1e-17 | 2.1e-17 | 3.0e-17 | 3.5e-18 |
| 60-80 | 3.6e-17 | 4.4e-17 | 5.1e-17 | 1.2e-17 | 2.7e-17 | 2.1e-17 | 5.4e-17 | 1.0e-17 |
| 80-91.2 | 3.2e-17 | 3.9e-17 | 7.0e-17 | 2.5e-17 | 2.2e-17 | 2.6e-17 | 5.3e-17 | 1.6e-17 |
| 91.2-100 | 3.9e-17 | 3.5e-17 | 7.2e-17 | 2.3e-17 | 4.3e-17 | 2.3e-17 | 5.4e-17 | 3.6e-17 |
| 100-110 | 1.6e-17 | 3.7e-17 | 5.7e-17 | 1.0e-17 |  -  | 1.3e-17 | 3.3e-17 | 2.0e-17 |
| 110-120 | 1.6e-17 | 1.8e-17 | 2.5e-17 | 1.0e-17 |  -  | 3.4e-18 | 1.9e-17 | 1.9e-19 |
| 120-123 | 1.3e-17 | 2.3e-18 | 1.7e-17 | 9.7e-18 |  -  | 2.4e-18 | 2.3e-17 |  -  |
| 123-130 | 1.6e-17 | 3.6e-17 | 2.3e-17 | 8.3e-18 |  -  | 2.3e-18 | 2.5e-17 |  -  |
| 130-140 | 1.1e-17 | 2.0e-18 | 1.8e-17 | 8.2e-18 |  -  | 2.6e-18 | 1.8e-17 |  -  |
| 140-150 | 5.8e-18 | 2.7e-18 | 3.2e-17 | 6.2e-18 |  -  | 3.1e-18 | 1.0e-17 |  -  |
| 150-160 | 1.4e-17 | 3.3e-18 | 9.2e-17 | 2.2e-18 |  -  | 1.7e-18 | 2.3e-17 |  -  |
| 160-170 | 6.0e-19 | 2.4e-18 | 2.0e-17 | 1.2e-18 |  -  | 1.2e-18 | 3.1e-17 |  -  |
| 170-180 | 8.9e-18 | 1.1e-18 | 1.1e-17 | 1.2e-18 |  -  | 3.1e-18 | 1.5e-17 |  -  |
| 180-190 | 8.6e-20 | 3.2e-19 | 5.6e-19 | 9.5e-19 |  -  | 3.0e-18 | 4.8e-19 |  -  |
| 190-200 | 0 | 7.0e-20 | 1.8e-20 | 5.7e-19 |  -  | 9.6e-19 | 9.4e-20* |  -  |
| 200-210 | 0 | 1.1e-20 | 8.1e-20 | 4.2e-19 |  -  | 2.3e-19 |  -  |  -  |
| 210-220 | 0 | 1.4e-21 | 1.9e-19 | 3.1e-19 |  -  | 2.9e-19 |  -  |  -  |
| 220-230 | 9.8e-23 | 1.3e-22 | 2.8e-19 | 2.2e-19 |  -  | 1.7e-19 |  -  |  -  |
| 230-240 | 4.7e-22 |  -  | 1.8e-19 | 1.6e-19 |  -  |  -  |  -  |  -  |
| 240-260 | 2.9e-21 |  -  | 3.1e-20 | 9.0e-20 |  -  |  -  |  -  |  -  |
| 260-280 | 1.3e-20 |  -  | 5.5e-21 | 3.9e-20 |  -  |  -  |  -  |  -  |
| 280-300 | 2.6e-20 |  -  |  -  | 1.2e-20 |  -  |  -  |  -  |  -  |
| 300-320 | 2.8e-20 |  -  |  -  | 4.0e-21 |  -  |  -  |  -  |  -  |
| 320-340 | 1.9e-20 |  -  |  -  | 1.4e-21 |  -  |  -  |  -  |  -  |
| 340-360 | 6.1e-21 |  -  |  -  | 6.0e-22* |  -  |  -  |  -  |  -  |
| 360-400 | 1.9e-22* |  -  |  -  |  -  |  -  |  -  |  -  |  -  |
| 400-450 |  -  |  -  |  -  |  -  |  -  |  -  |  -  |  -  |
| 450-500 |  -  |  -  |  -  |  -  |  -  |  -  |  -  |  -  |
| 500-550 |  -  |  -  |  -  |  -  |  -  |  -  |  -  |  -  |
| 550-600 |  -  |  -  |  -  |  -  |  -  |  -  |  -  |  -  |
| 600-650 |  -  |  -  |  -  |  -  |  -  |  -  |  -  |  -  |
| 650-700 |  -  |  -  |  -  |  -  |  -  |  -  |  -  |  -  |
| 700-800 |  -  |  -  |  -  |  -  |  -  |  -  |  -  |  -  |
| 800-1100 |  -  |  -  |  -  |  -  |  -  |  -  |  -  |  -  |


*Status for every cell of the three tables above: DATA, COMPUTED.*

Reading the tables (all COMPUTED from the file unless stated):

- **O2.** Schumann-Runge continuum peak: bin mean 1.4e-17 at 140 to 150 nm. Schumann-Runge bands: 1.9e-19 (175 to 180 nm), 1.7e-20 (180 to 190), 4.0e-22 (190 to 200). Herzberg continuum, 2 nm means: 1.5e-23 at 200 nm, 6.5e-24 at 205, 4.2e-24 at 220, 2.3e-24 at 230, 1.0e-24 at 240. Lyman alpha (121.6 nm): 1.1e-20, the window in the continuum. For 1 bar of 20.95 percent O2 (N_O2 = 4.5e24 cm^-2, TEST point) tau_O2 = 29 at 205 nm, 7 at 235 nm and 4.7 at 239 nm, so O2 alone blocks nearly everything below 240 nm. **Band-mean caveat (COMPUTED on the file's 0.1 nm grid, already smoother than the real lines):** transmission exp(-<sigma>N) against the true mean of exp(-sigma N) is 0.152 against 0.253 at N = 1e19 cm^-2 in the 175 to 180 nm bin, and 6.8e-9 against 2.6e-3 at 1e20; in the 180 to 190 nm bin 0.174 against 0.476 at 1e20. Band means understate transmission of thick columns of structured absorbers by orders of magnitude; the generator needs a k-distribution or effective cross section for O2 at 175 to 200 nm. The error matters above 1e18 cm^-2 and only for structured bands; the Hartley band is smooth.
- **O3.** Hartley peak 1.15e-17 cm2 at 255.4 nm (Serdyuchenko 293 K file, 1.1507e-17; Malicet/Brion 295 K in TUV-x, 1.149e-17), Huggins 310 to 360 nm falling by about 3 decades, Chappuis peak 5.1e-21 to 5.2e-21 cm2 at 602 nm. The Chappuis peak from three independent files agrees within 2.3 percent (5.11e-21 Serdyuchenko, 5.16e-21 Burkholder and Talukdar through Mainz, 5.23e-21 Brion/Malicet through TUV-x); the Hartley peak from two files within 0.2 percent.
- **SO2.** Two strong systems, 190 to 230 nm (bin means to 7.3e-18) and 260 to 330 nm (7.5e-19 at 280 to 300), weak forbidden absorption to 400 nm; SO2 shields 190 to 330 nm in a volcanic atmosphere.
- **H2O** bands near 120 to 140 and 150 to 185 nm (1.0e-17 at 120 to 123 nm, 5.3e-18 at 160 to 170); below 2e-21 beyond 190 nm, so water vapour does not shield the ground UV. **CO2** is below 1e-24 beyond 200 nm at 270 K but rises by decades when hot (2.4).
- **NO2** is the only listed species with a broad visible band: 5.6e-19 cm2 at 360 to 450 nm, 6.2e-20 at 550 to 600, 9e-21 at 650 to 700; a column of 1e17 cm^-2 gives tau = 0.056 at 400 to 450 nm, a yellow-brown tint.
- **CH4** visible and near-IR: 4.4e-26 cm2 at 600 to 650 nm, 2.4e-25 at 700 to 800, 1.6e-24 at 800 to 1100 (Karkoschka 1994). Giant-planet methane bands (near 619, 725, 890 nm, RECALLED) are one to two decades stronger than these bin means and need a band model (part 03).
- **H2S, NH3, HCN, C2H2, H2CO, HCl, OCS, H2O2** absorb below about 230 to 350 nm only; H2CO has a weak band 250 to 360 nm (2.8e-20 at 300 to 320 nm).
- **N2, H2, He, Ne, Ar, Kr, Xe, CO, C2H6** are transparent above 150 to 170 nm and never shield the surface.
- Compounds the project text mentions that have no table here: HF (Mainz `HF.spc` exists, file not read), CH3Cl, SiO, Na, K, perchlorates, carbonates, radon, NOx other than NO and NO2. They are listed in section 4.

### 2.4 Temperature dependence of cross sections

| Species and band | Data | Result | Source and status |
|---|---|---|---|
| O3 Hartley, 215 to 285 nm | Serdyuchenko, 11 temperatures 193 to 293 K | ratio 193 K to 293 K within 0.97 to 1.04 (Hartley peak +1.6 percent at 193 K) | DATA (Bremen file, no licence stated on the page; citation required), COMPUTED |
| O3 Huggins, 295 to 400 nm | same | ratio 193/293 K: 0.89 (295 to 305), 0.82 (305 to 315), 0.75 (315 to 325), 0.66 (325 to 335), 0.54 (335 to 345), 0.34 (345 to 355), 0.17 (355 to 375), 0.35 (375 to 400) | DATA, COMPUTED |
| O3 Chappuis, Serdyuchenko | same | 193/293 K: 0.85 (400 to 450), 0.95 (450 to 500), 0.99 (500 to 550), 1.01 (550 to 600), 1.00 (600 to 650), 0.97 (650 to 700), 0.95 (700 to 800); above 900 nm the file is noise-dominated | DATA, COMPUTED |
| O3 Chappuis, Burkholder and Talukdar 1994 | Mainz atlas files `O3_BurkholderTalukdar(1994)_{220,240,260,280,298}K_407-762nm.txt` (the atlas notes the cross sections were rebuilt from sigma(T)/sigma(298 K) ratios received from J. B. Burkholder in January 2005) | 220 K to 298 K ratio: 0.682 at 420 nm, 0.843 at 450, 0.974 at 500, 0.994 at 550, 1.007 at 602, 0.986 at 650, 0.935 at 700, 0.944 at 750; peak 5.163e-21 cm2 at 602 nm (298 K) | DATA (the 1994 paper itself NOT REACHED), COMPUTED |
| N2O, 200 to 230 nm | Mainz, 240 K and 296 K, 2 nm means | 240 K to 296 K ratio 0.81 (200 nm), 0.68 (210), 0.54 (220), 0.44 (230) | DATA, COMPUTED |
| SO2, 250 to 320 nm | Mainz, 203 K and 293 K, 2 nm means | no change at 250 to 300 nm (1.06e-19 against 1.05e-19 at 250; 7.69e-19 against 7.66e-19 at 280); 0.91 at 320 nm; above 340 nm the 203 K file is noise dominated | DATA, COMPUTED |
| SO2 hot | Mainz, 800 K (a few points) | 1.9e-19 at 250 nm against 1.1e-19 in the 293 K 240 to 260 nm bin (about 1.7 times); 6.9e-19 at 300 nm | DATA, coarse |
| CO2 hot | Mainz, 880, 1630, 3050 K, 2 nm means | 200 nm: 4.7e-21 (880 K), 1.0e-19 (1630 K), 4.8e-19 (3050 K), against 5e-25 in the 200 to 210 nm bin at 270 K (rise of 4 to 6 decades); 220 nm: 2.7e-20 (1630 K), 2.5e-19 (3050 K); 240 nm: 7.7e-21 (1630 K), 1.3e-19 (3050 K) | DATA, COMPUTED |
| H2S hot | Mainz, 423 K and 773 K | 210 nm: 3.2e-18 to 3.3e-18; 250 nm: 2.5e-20 to 4.9e-20; 300 nm: 2.4e-21 to 1.1e-20 (factor 2 to 4.6 per 350 K) | DATA, COMPUTED |
| NH3, 175 and 195 K | Mainz | 175 K and 195 K identical within 1 percent at 180 to 210 nm (peak 1.14e-17 at 190 nm) | DATA, COMPUTED |
| Inventory only | Mainz atlas set lists | CO2: 150 to 800 K plus shock-tube to 3050 K; N2O: 148 to 2000 K; H2S: 423, 573, 773 K; O2: at 203, 243 K plus Vattulainen and Zuev at high temperature; NH3: 70, 175, 195 K besides room temperature | READ (page listings), values not extracted |

Rule of thumb (design choice): O3 and NO2 need no temperature correction above 450 nm; the O3 Huggins band must be scaled when the stratosphere is below 250 K (the Bremen page offers a quadratic fit file for 300 to 370 nm, not downloaded); hot CO2, SO2 and H2S atmospheres absorb in the UV at 2 to 6 decades higher cross section than at 300 K.

### 2.5 Alkali metals, metal oxides and hydrides on giants and hot planets

**2.5a Resonance lines (neutral alkali and a few refractory atoms).** Wavelengths (air), oscillator strengths f and Einstein A are from NIST ASD (physics.nist.gov, public domain; DATA). Integrated cross section S = pi r_e lambda^2 f (cm2 nm) and the Doppler line-centre cross section sigma0 = sqrt(pi) r_e c f lambda/b with b = sqrt(2kT/m) at 1000 K and 2000 K are COMPUTED by me. Natural widths are 1e-5 to 2e-5 nm, so lines are pressure-broadened and the wings dominate.

| Line | lambda (nm) | f | A (s^-1) | S (cm2 nm) | sigma0 at 1000 K | sigma0 at 2000 K | Strength rank | Status |
|---|---|---|---|---|---|---|---|---|
| Li I D2 / D1 | 670.776 / 670.791 | 0.498 / 0.249 | 3.69e7 | 1.98e-14 / 0.99e-14 | 3.2e-12 / 1.6e-12 | 2.3e-12 / 1.1e-12 | 3 (mid) | DATA, COMPUTED |
| Na I D2 / D1 | 588.995 / 589.592 | 0.641 / 0.320 | 6.16e7 / 6.14e7 | 1.97e-14 / 0.99e-14 | 6.7e-12 / 3.3e-12 | 4.7e-12 / 2.4e-12 | 3 (mid) | DATA, COMPUTED |
| K I | 766.490 / 769.896 | 0.666 / 0.332 | 3.78e7 / 3.73e7 | 3.46e-14 / 1.74e-14 | 1.2e-11 / 5.9e-12 | 8.3e-12 / 4.2e-12 | 2 (strong) | DATA, COMPUTED |
| Rb I | 780.027 / 794.760 | 0.696 / 0.342 | 3.81e7 / 3.61e7 | 3.75e-14 / 1.91e-14 | 1.8e-11 / 9.2e-12 | 1.3e-11 / 6.5e-12 | 2 (strong) | DATA, COMPUTED |
| Cs I | 852.113 / 894.347 | 0.714 / 0.344 | 3.28e7 / 2.86e7 | 4.59e-14 / 2.43e-14 | 2.6e-11 / 1.3e-11 | 1.8e-11 / 9.2e-12 | 1 (strongest) | DATA, COMPUTED |
| Ca I | 422.673 | 1.75 | 2.18e8 | 2.77e-14 | 1.7e-11 | 1.2e-11 | 2 | DATA, COMPUTED |
| Ca II K | 393.366 | 0.682 | 1.47e8 | 9.3e-15 | 6.2e-12 | 4.4e-12 | 3 | DATA (NIST accuracy C), COMPUTED |
| Mg I | 285.213 | 1.80 | 4.91e8 | 1.30e-14 | 9.3e-12 | 6.6e-12 | 3 | DATA, COMPUTED |

Profile physics (Burrows and Volobuyev 2003, arXiv:astro-ph/0210086, READ): a Lorentzian core holds only to about 20 A detuning; beyond it the wings of Na D (5890 A) and K I (7700 A) pressure-broadened by H2 and He form the "continuum" of T dwarfs between 0.4 and 1.0 um, with natural red and blue cutoffs, and "rainbow" satellite features in the blue wings predicted at 0.515 um (Na) and 0.69 um (K) against the measured 0.55 um and 0.72 um (Na+He and K+He). Numerical widths are NOT REACHED here (section 4). Li I 670.8 nm is a substellarity test line; the L-dwarf optical-spectra paper astro-ph/0001060 (READ) finds its equivalent width is changed more by the added opacity of other absorbers than by Li depletion. Rb I and Cs I appear in L dwarfs (Cushing et al. 2005, astro-ph/0412313, READ).

**2.5b Molecular band systems (positions only; no cross sections reached).**

| Molecule | System and band heads (um) | Where it matters | Source | Status |
|---|---|---|---|---|
| TiO | gamma (A 3Phi-X 3Delta), gamma' (B 3Pi-X 3Delta), delta, epsilon, phi systems; red-optical heads near 0.705, 0.76 and 0.84 um; phi 0-0 head predicted near 1.104 um but not seen in M or L dwarfs | M dwarfs and hot giants (pM class) | Cushing et al. 2005 (system list and the 1.104 um non-detection, READ); the optical head positions RECALLED | READ, RECALLED |
| VO | B 4Pi-X 4Sigma red-optical bands; A 4Pi-X 4Sigma: 1-0 at 0.96, 0-0 at 1.06, 0-1 at 1.18 um | late M dwarfs, hot giants | Cushing et al. 2005 (READ) | READ |
| FeH | F 4Delta-X 4Delta Wing-Ford band, 0-0 head at 0.9969 um with a Q branch at 1.0060 um (Cushing 2005); the same head is quoted as 989.6 nm by Lodders and Fegley 2006 quoting Burgasser: spread 0.9896 to 0.9969 um, unresolved; 0-1 head 1.1939 um, Q branch 1.2221 um, 1-2 head 1.2389 um | L and T dwarfs, cool giants | Cushing et al. 2005 (READ); Lodders and Fegley 2006, astro-ph/0601381 (READ, Fig. 2 and spectral-type table) | READ, 2ND-HAND |
| CrH | A 6Sigma+-X 6Sigma+: 0-0 head 0.8611 um (defines the L class), 0-1 head 0.9969 um (but many features there are FeH), 0-2 near 1.18 um, 0-3 near 1.4 um; Burrows et al. 2002 computed line lists for 12 bands | L dwarfs, brown dwarfs | Cushing et al. 2005 (READ); Burrows, IAU 211 (astro-ph/0207399, READ, only the list of the bands) | READ |

**2.5c Which regime applies (hot-Jupiter classes and condensation).** Fortney et al. 2008 (arXiv:0710.2558, READ): planets warm enough to keep Ti and V in the gas, so TiO and VO bands dominate the optical and drive a thermal inversion, are "pM class"; cooler ones, where Ti and V are condensed, are "pL class", and their visible spectrum is set by alkali lines and water. At 1 mbar and 1700 K (before Ti/V condensation) the V-band opacity is about log kappa = -0.5 in cm2 g^-1; at 1600 K, only 33 K below the 90 percent Ti-condensation curve, it is still an order of magnitude above the value at 1400 K, where the strong TiO and VO bands have given way to alkali lines and water. The first Ti condensate (solar metallicity) depends on pressure: TiN above about 30 bar, Ca3Ti2O7 for 5 to 30 bar, Ca4Ti3O10 for 0.03 to 5 bar, CaTiO3 below 0.03 bar (Fortney quoting Lodders and Fegley 2006); my own reading of Lodders and Fegley 2006 gives different boundaries (CaTiO3 below 0.03 bar, Ca4Ti3O10 from 0.03 to 3.2 bar, Ca3Ti2O7 3.2 to about 10 bar, TiN above about 10 bar). The spread in the pressure boundaries is a factor of 2 to 3 and should be drawn. Na and K condense as Na2S and KCl at lower temperature (Lodders and Fegley 2006 Fig. 2a; curves in part 06). All READ.

### 2.6 Action spectra for surface UV

| Item | Definition | Source | Status |
|---|---|---|---|
| Erythemal action spectrum S_er, CIE 1998 / ISO 17166 | S = 1 for 250 to 298 nm; S = 10^(0.094 (298 - lambda)) for 298 < lambda <= 328 nm; S = 10^(0.015 (140 - lambda)) for 328 < lambda <= 400 nm; 0 outside | the reviewer exchange in the GMD discussion (gmd.pdf, READ) states the CIE (2014) piecewise function with the exponent (140 - lambda) and says the NOAA text has (139 - lambda) as a typo; same text gives the 328 nm switch | 2ND-HAND |
| Values of S_er | 0.649 at 300 nm, 0.220 at 305, 0.0745 at 310, 0.0086 at 320, 1.41e-3 at 330, 7.1e-4 at 350, 3.1e-4 at 380 | evaluated from the formula | COMPUTED |
| UV index | UVI = (integral of E(lambda) S_er(lambda) dlambda) / (25 mW m^-2) = 40 m2 W^-1 times the erythemally weighted irradiance; designed so that typical midday summer sunlight was about 250 mW m^-2, hence 0 to about 11 | Wikipedia UV Index article (reachable page, cites EPA/WHO) | 2ND-HAND |
| UV band names | UV-A 315 to 400 nm, UV-B 280 to 315 nm, UV-C 100 to 280 nm | same page | 2ND-HAND |
| Setlow 1974 DNA-damage action spectrum | digitised from the hand-drawn Setlow plot (received from R. Setlow in 1995), per quantum (multiply by lambda/300 for an irradiance in energy units). Weights: 1.07 (256 nm), 1.12 (260), 1.23 (266 to 268), 1.13 (270), 0.933 (276), 0.791 (280), 0.530 (286), 0.309 (290), 0.163 (294), 0.0552 (298), 0.0333 (300), 0.00901 (304), 0.00115 (310), 2.5e-4 (314), 4.3e-5 (320), 1.1e-5 (328), 4.6e-6 (340), 2.6e-6 (350), 1.8e-6 (360). The file header says the curve is normalised to unity at 300 nm but the stored weight at 300 nm is 0.0333; this inconsistency is reported, not resolved (ratio 260 to 300 nm is 33.6 in the file) | TUV-x data `DNA_damage_in_vitro(Setlow,1974)_spectral_wght_1.nc`, in `musica` 0.17.1 (Apache-2.0 for the package; data licence not found) | DATA |
| Other action spectra shipped with TUV-x | Anders 1995 erythema (minimal doses, a different quantity), CIE 2006 non-melanoma skin cancer, previtamin D3 (CIE 2006), SCUP-human, plant, phytoplankton and cataract weights | same package, folder `spectral_weights` | DATA (names read, values not extracted) |
| Solar reference | ATLAS-1 (Thuillier et al. 2004, GM 141) spectrum at 1 au, 26150 points, moderately high activity | `photochem_clima_data` `stellar_spectra/ATLAS1_2004.h5`, source sbuv.gsfc.nasa.gov | DATA |

### 2.7 Surface UV demonstration (Beer-Lambert, computed)

Setup (COMPUTED): ATLAS-1 solar spectrum on a 0.5 nm grid, 200 to 400 nm; Rayleigh tau of a 1 bar air column (0.0969 at 550 nm scaled by the air refractivity formula of the earlier report); O3 at 223 K (Serdyuchenko above 215 nm, Photochem below); O2 band means from the Photochem file (20.95 percent O2, N = 4.5e24 cm^-2); plane-parallel air mass 1/cos(z). Direct beam on a horizontal surface = F0 cos z exp(-m tau). "Crude global" adds a diffuse term equal to half of the Rayleigh-scattered-out direct light, attenuated by the O3 path; it is an order-of-magnitude estimate (+-50 percent), not a radiative-transfer result. No aerosol, cloud, ground albedo or altitude. TOA values: 200 to 280 nm 6.85 W m^-2, 280 to 315 nm 17.44, 315 to 400 nm 86.00; erythemally weighted UVI-equivalent at the top of the atmosphere 666; Setlow-weighted 8.66 W m^-2.

| O3 column (DU) | Zenith | UVI direct | UVI crude global | UV-C 200 to 280 nm direct (W m^-2) | UV-B 280 to 315 nm direct (W m^-2) | Setlow-weighted direct (W m^-2) |
|---|---|---|---|---|---|---|
| 0 | 0 deg | 133.6 | 370 | 0.704 | 5.24 | 1.23 |
| 0 | 60 deg | 16.5 | 160 | 0.0506 | 0.809 | 0.115 |
| 100 | 0 deg | 23.4 | 48.9 | 5.4e-7 | 2.47 | 0.0339 |
| 100 | 60 deg | 1.98 | 9.5 | 4e-13 | 0.296 | 1.5e-3 |
| 200 | 0 deg | 11.5 | 22.8 | 4e-12 | 1.70 | 9.7e-3 |
| 300 (TEST) | 0 deg | 7.27 | 13.9 | 1e-16 | 1.27 | 4.4e-3 |
| 300 (TEST) | 60 deg | 0.69 | 2.5 | 0 | 0.112 | 1.9e-4 |
| 500 | 0 deg | 4.11 | 7.4 | 6e-19 | 0.776 | 1.6e-3 |

Readings (COMPUTED): the direct UVI falls as DU^-1.08 from 100 to 500 DU at overhead sun (the usual radiation amplification factor for erythema is quoted near 1.1, RECALLED); the 300 DU, overhead-Sun crude global UVI of 13.9 is the right order for a clear tropical sky (the Wikipedia page says values in the tens occur in the tropics and mountains), which is a plausibility check, not a validation; **with no ozone at all, O2 alone still removes light below about 240 nm, but 0.70 W m^-2 of 240 to 280 nm UV-C reaches the ground at overhead Sun, with a Setlow-weighted dose about 280 times that of the 300 DU case**, so a planet with O2 but no ozone has lethal UV-C at the surface unless it has another absorber (SO2, thick haze, or ocean cover). With 100 DU the UV-C is already 1e-7 W m^-2. The generator therefore needs ozone column, O2 column and the stellar UV, in that order of importance.

### 2.8 Air mass and scattering checks (COMPUTED)

**Chapman function Ch(X, chi)** by direct numerical integration of exp(-(r - R)/H) along the slant ray for an exponential atmosphere, X = R/H, normalised to 1 at the zenith (Earth-like X is 757 for H = 8.4 km, a TEST value; gas giants 3000 or more; small moons 50 or less). Compare plane-parallel sec(chi) and the Kasten and Young 1989 fit (Earth only; formula AM = 1/(cos z + 0.50572 (96.07995 - z)^-1.6364), z in degrees; the same formula is in pvlib and on the Wikipedia air-mass page, both 2ND-HAND).

| chi (deg) | sec chi | X = 50 | X = 200 | X = 757 | X = 3000 | Kasten and Young |
|---|---|---|---|---|---|---|
| 0 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| 60 | 2.000 | 1.901 | 1.972 | 1.992 | 1.998 | 1.994 |
| 75 | 3.864 | 3.228 | 3.640 | 3.796 | 3.846 | 3.813 |
| 80 | 5.759 | 4.189 | 5.104 | 5.541 | 5.699 | 5.586 |
| 85 | 11.47 | 5.824 | 8.277 | 10.09 | 11.03 | 10.31 |
| 88 | 28.65 | 7.421 | 12.48 | 18.44 | 23.87 | 19.43 |
| 90 | infinite | 8.928 | 17.76 | 34.50 | 68.66 | 37.92 |

At 90 degrees the limit is sqrt(pi X/2) = 34.48 for X = 757, agreeing with the integral (34.50). The Kasten and Young fit is 10 percent above the exponential-atmosphere value at the horizon and cannot be reused for other worlds. Rule: use the Chapman integral (or the Chapman 1931 approximation) with X = R/H, never sec(chi) beyond chi = 75 degrees. Status: COMPUTED.

**Rayleigh phase function with depolarization** (COMPUTED by orientation-averaged Monte Carlo of random molecules; matches the formula to 4 digits): p(theta) = 3/(4 (1 + 2 g)) [(1 + 3 g) + (1 - g) cos^2 theta], g = rho_n/(2 - rho_n), normalised so that (1/2) integral of p over mu = 1. Forward-to-side ratio I(0)/I(90) = 1.8006 for rho = 0.022 (Monte Carlo 1.8009), 1.9569 for rho = 0.0004 (1.9566) and 2.0 for rho = 0; the 90-degree depolarization of the simulation equalled the input.

**Two-stream test** (COMPUTED, scalar Monte Carlo with 3e5 photons per case, conservative Rayleigh scattering, no polarization): Eddington transmission T(tau, mu0) = [(2/3 + mu0) + (2/3 - mu0) exp(-tau/mu0)]/(4/3 + tau), reflectance R = 1 - T, with ground albedo A and slab spherical albedo s the top-of-atmosphere reflectance is R_tot = R + T A (1 - s)/(1 - A s). Against the Monte Carlo (36 cases, A of 0, 0.3, 0.8; mu0 of 1, 0.5, 0.2; tau of 0.1 to 3), R_tot is within 1.3 percent for every case with mu0 of 1 or 0.5; at mu0 = 0.2 it is 0.1 to 3.5 percent low for A = 0 and 0.3 (worst: Monte Carlo 0.650 against Eddington 0.627 at A = 0, tau = 1) and within 0.5 percent for A = 0.8. The closed form s = 0.75 tau/(1 + 0.75 tau) is good only above tau = 2 (exact integral of Eddington R against the closed form: 0.084 against 0.070 at tau 0.1; 0.292 against 0.273 at 0.5; 0.440 against 0.429 at 1; 0.603 against 0.600 at 2); below tau = 2 integrate s = integral 2 mu (1 - T) dmu numerically.

### 2.9 Stellar light to colour

**Blackbody chromaticity (CIE 1931 2-degree), COMPUTED.** Inputs: the CIE 1931 colour matching functions from colour-science (BSD-3-Clause; the package's CIE 1931 2-degree table), Planck law with c2 = 1.438776877e-2 m K, trapezoid integration; sRGB (D65) matrix. "Linear" is linear sRGB scaled to the largest channel (a chromaticity, not a brightness); 8-bit is gamma-encoded. dxy is the distance to the Planckian-locus polynomial (as coded in colour-science, valid 1667 to 25000 K); the largest dxy is 0.0005.

| T (K) | x | y | locus x | locus y | dxy | linear sRGB / max (R, G, B) | sRGB 8-bit |
|---|---|---|---|---|---|---|---|
| 2000 | 0.5267 | 0.4133 | 0.5269 | 0.4133 | 0.0002 | 1.000, 0.257, 0.008 | 255, 139, 22 |
| 2500 | 0.4770 | 0.4137 | 0.4765 | 0.4137 | 0.0005 | 1.000, 0.372, 0.068 | 255, 164, 73 |
| 3000 | 0.4369 | 0.4041 | 0.4366 | 0.4042 | 0.0004 | 1.000, 0.477, 0.154 | 255, 184, 109 |
| 3500 | 0.4053 | 0.3907 | 0.4053 | 0.3909 | 0.0001 | 1.000, 0.571, 0.259 | 255, 199, 139 |
| 4000 | 0.3804 | 0.3767 | 0.3805 | 0.3767 | 0.0001 | 1.000, 0.654, 0.377 | 255, 211, 165 |
| 4500 | 0.3608 | 0.3635 | 0.3607 | 0.3635 | 0.0001 | 1.000, 0.727, 0.502 | 255, 222, 188 |
| 5000 | 0.3451 | 0.3516 | 0.3450 | 0.3516 | 0.0001 | 1.000, 0.792, 0.628 | 255, 230, 208 |
| 5500 | 0.3324 | 0.3410 | 0.3323 | 0.3410 | 0.0001 | 1.000, 0.849, 0.754 | 255, 237, 225 |
| 5772 | 0.3266 | 0.3358 | 0.3265 | 0.3358 | 0.0001 | 1.000, 0.877, 0.820 | 255, 241, 234 |
| 6000 | 0.3221 | 0.3318 | 0.3220 | 0.3318 | 0.0001 | 1.000, 0.899, 0.875 | 255, 243, 240 |
| 6504 | 0.3135 | 0.3236 | 0.3134 | 0.3236 | 0.0000 | 1.000, 0.943, 0.993 | 255, 249, 254 |
| 7000 | 0.3064 | 0.3165 | 0.3064 | 0.3166 | 0.0000 | 0.906, 0.890, 1.000 | 244, 242, 255 |
| 8000 | 0.2952 | 0.3048 | 0.2952 | 0.3048 | 0.0001 | 0.765, 0.802, 1.000 | 227, 231, 255 |
| 10000 | 0.2806 | 0.2883 | 0.2807 | 0.2883 | 0.0001 | 0.609, 0.695, 1.000 | 205, 217, 255 |
| 12000 | 0.2718 | 0.2776 | 0.2718 | 0.2776 | 0.0000 | 0.527, 0.635, 1.000 | 192, 209, 255 |
| 15000 | 0.2637 | 0.2673 | 0.2637 | 0.2672 | 0.0000 | 0.460, 0.582, 1.000 | 181, 201, 255 |
| 20000 | 0.2565 | 0.2576 | 0.2564 | 0.2576 | 0.0000 | 0.406, 0.536, 1.000 | 171, 193, 255 |
| 25000 | 0.2525 | 0.2522 | 0.2525 | 0.2523 | 0.0001 | 0.378, 0.512, 1.000 | 165, 190, 255 |
| 30000 | 0.2501 | 0.2488 | - | - | - | 0.361, 0.497, 1.000 | 162, 187, 255 |

Out-of-gamut note: 2000 K nearly clips the sRGB blue channel; chromaticity (x, y) is the portable result. White-balance caveat: the project decided on absolute colour, so no adaptation to the star is applied.

**Non-blackbody spectra, M dwarfs and cooler (COMPUTED from DATA).** Spectra: STScI reference-atlas BT-Settl grid `phoenixm00_<Teff>.fits` (solar metallicity, log g = 5.0, fluxes in erg cm^-2 s^-1 A^-1) from archive.stsci.edu/hlsps/reference-atlases/cdbs/grid/phoenix/phoenixm00/. Same colorimetry. The offsets against a blackbody at the same Teff:

| Teff (K) | log g | BT-Settl x, y | blackbody x, y | dx, dy | BT-Settl sRGB 8-bit | blackbody sRGB 8-bit | BT-Settl flux ratio 450/650 nm vs blackbody |
|---|---|---|---|---|---|---|---|
| 2300 | 5.0 | 0.4857, 0.3930 | 0.4957, 0.4152 | -0.0100, -0.0222 | 255, 149, 81 | 255, 155, 56 | 0.89 |
| 2500 | 5.0 | 0.4275, 0.3983 | 0.4770, 0.4137 | -0.0495, -0.0154 | 255, 187, 119 | 255, 164, 73 | 1.60 |
| 2800 | 5.0 | 0.4246, 0.3983 | 0.4519, 0.4086 | -0.0273, -0.0104 | 255, 189, 121 | 255, 176, 96 | 2.14 |
| 3000 | 5.0 | 0.4162, 0.3962 | 0.4369, 0.4041 | -0.0208, -0.0079 | 255, 194, 129 | 255, 184, 109 | 0.77 |
| 3300 | 5.0 | 0.4152, 0.3933 | 0.4170, 0.3963 | -0.0018, -0.0029 | 255, 193, 130 | 255, 193, 128 | 1.19 |
| 3500 | 5.0 | 0.4079, 0.3868 | 0.4053, 0.3907 | +0.0026, -0.0039 | 255, 194, 139 | 255, 199, 139 | 1.33 |
| 4000 | 5.0 | 0.3985, 0.3744 | 0.3804, 0.3767 | +0.0180, -0.0023 | 255, 194, 151 | 255, 211, 165 | 0.99 |
| 5000 | 5.0 | 0.3532, 0.3457 | 0.3451, 0.3516 | +0.0081, -0.0059 | 255, 216, 200 | 255, 230, 208 | 0.56 |
| 5800 | 5.0 | 0.3266, 0.3308 | 0.3260, 0.3353 | +0.0007, -0.0046 | 255, 236, 234 | 255, 241, 235 | 0.94 |

Reading: the offset reaches 0.05 in x at 2500 K and 0.018 at 4000 K and is irregular because TiO and VO bands cut the red-optical flux (the 450/650 nm ratio changes by factors 0.56 to 2.1 without a trend). At 5000 K and above the blackbody is within 0.008 in x and 0.006 in y. Rule: blackbody above about 4500 K; for cooler stars use the BT-Settl chromaticities (or integrate the file) with a seeded spread of 0.02 in x and 0.01 in y.

### 2.10 Source re-check (extra task, after the network was widened)

| Source | Result | Label |
|---|---|---|
| Peck and Huang 1977, JOSA 67, 1550 | Optica page gives only the citation (Edson R. Peck and Shiezen Huang, pp. 1550-1554). Dispersion constants READ in refractiveindex.info `H2_nk_Peck.yml` (CC0; spells the second author "Hung"): 0.0148956, 180.7, 0.0049037, 92 at 273.15 K, 0.168 to 1.6945 um | DATA (constants); paper NOT REACHED |
| Burkholder and Talukdar 1994, GRL 21, 581 | Mainz files for 220, 240, 260, 280, 298 K, 407 to 762 nm READ (356 points each); peak 5.163e-21 cm2 at 602 nm; the atlas says the cross sections were rebuilt from ratios supplied by the authors in 2005 | DATA; paper NOT REACHED |
| Serdyuchenko 2014 (AMT 7, 625) | Bremen data file READ (11 temperatures 193 to 293 K, 213 to 1100 nm, 0.01 nm grid, systematic uncertainty under 0.5 percent); the page asks users to check for updates and states no licence | DATA; licence NOT REACHED |
| Sneep and Ubachs 2005; Bates 1984 | ScienceDirect error page for both | still 2ND-HAND (He 2021 eq. 1 and Table 1; five codes in the report) |
| Bucholtz 1995 Table 3 | Optica page gives citation and tags only; PDF behind a captcha | NOT REACHED |
| Penndorf 1957; Hulburt 1953 | Optica pages confirm the citations (Penndorf JOSA 47, 176-182; Hulburt JOSA 43, 113-118, "Explanation of the Brightness and Color of the Sky, Particularly the Twilight Sky"); full text paywalled | NOT REACHED (text); Penndorf formula stays 2ND-HAND |
| Kasting and Donahue 1980 | Cloudflare pages; NTRS search zero hits | still 2ND-HAND (through Cooke 2022) |
| no source | archive.org finds three items, lending copies return HTTP 401. The ground-bounce factor R_tot = R + T A (1 - s)/(1 - A s) was instead re-derived and checked by Monte Carlo (2.8) | NOT REACHED; formula COMPUTED |
| Kasten and Young 1989 | Optica captcha; the Wikipedia air-mass page prints the formula and the citation | 2ND-HAND (formula) |

---

## 3. Formulas with validity ranges

### 3.1 Rayleigh cross section of any gas (GEN.417)

    sigma_i(lam) = 1.3057e-28 * alpha_i^2 * F_K,i / lam^4 * (1 + B/lam^2)^2        [cm2; alpha in A^3; lam in um]
    F_K = (6 + 3 rho)/(6 - 7 rho),   rho = 6 g2 / (45 a^2 + 7 g2)   (a = mean polarizability, g2 = anisotropy squared)
    tau_R(lam) = [ sum_i x_i sigma_i(lam) ] * P / (m_bar g)          [P column pressure, m_bar mean molecular mass in kg, g gravity; molecules per m2, convert to cm2]

The coefficient is 128 pi^5/3 in cgs. B = 0.0075 um^2 with a drawn spread of +-0.003. The fitted B of the 13 verified gases runs from 0.0024 (He), 0.0030 (Ne) and 0.0054 (N2, O2, Ar) through 0.0065 (CO2), 0.0074 (Kr), 0.0079 to 0.0088 (H2, CH4, CO, N2O) to 0.0114 (Xe) and 0.0115 (NH3) (COMPUTED); the B range changes sigma(550 nm) by only 2 percent. alpha is the CCCBDB static value; at long wavelength it agrees with the optical polarizability implied by the refractivity formulas to 0.98 to 1.05.

**Validation (COMPUTED, 13 gases).** sigma_generic/sigma_n-based at 550 nm with F_K taken from the n-based source: N2 0.979, O2 0.990, Ar 1.040, CO2 0.926, CH4 0.907, N2O 1.038, CO 1.013, H2 0.956, He 1.060, Ne 0.995, Kr 1.015, Xe 1.001, NH3 0.925 (Cuthbertson). Mean 0.988, standard deviation 0.046, range 0.907 to 1.060. At 450 and 650 nm the ratios change by under 0.01 for N2, CO2, CH4 and Xe. With F_K from the DFT rho instead of the source value the sigma changes by under 1.2 percent (section 2.2). **Draw the error as a log-normal with sigma 5 percent about the generic value when no refractivity exists.**

Validity: 300 to 1000 nm and any gas whose first electronic absorption lies below 200 nm (all gases in 2.1); molecules with strong visible absorption (NO2, O3 Chappuis, I2, S2) need the real refractivity near the band; below 300 nm the (1 + B/lam^2)^2 factor grows and the formula should not be used within about 80 nm of the first strong absorption. Mixtures: add cross sections linearly by number fraction; the verified Earth test (air, 550 nm, 288.15 K reference) is tau_R = 0.0969; the generic formula applied to N2, O2, Ar and CO2 (CCCBDB alpha, source F_K) gives 0.0951, 1.9 percent low (COMPUTED, TEST point).

### 3.2 Phase function, depolarization from anisotropy

p(Theta) = 3/(4 (1 + 2 g)) [(1 + 3 g) + (1 - g) cos^2 Theta] with g = rho/(2 - rho) (checked by Monte Carlo, section 2.8); for a mixture weight rho_eff by the scattering cross sections, g_eff = sum x_i sigma_i g_i / sum x_i sigma_i.

### 3.3 Beer-Lambert photolysis rate (feeds GEN.413 and the photochemistry lane)

    j_k(z) = integral over lam of sigma_k(lam, T) phi_k(lam) F(lam, z) dlam,    F(lam, z) = F_TOA(lam) exp(-m(chi) tau_slant(lam, z))
    tau_slant = sum over absorbers of N_col,j(z) sigma_abs,j(lam) + tau_R(lam, z)

with sigma_abs the photoabsorption (2.3) and sigma_k the photodissociation cross section times quantum yield from the same file. Use column amounts above the level z. For the O2 Schumann-Runge bands (175 to 200 nm) and any structured absorber do not use bin means for N above 1e18 cm^-2 (section 2.3 test); the stored Photochem 0.1 nm data are the minimum resolution. Add a scattered-light term where Rayleigh tau exceeds 0.3 (wavelengths below about 410 nm for 1 bar of air); the crude estimate of 2.7 puts the diffuse UV at or above the direct beam there.

### 3.4 Surface UV

    E_dir(lam) = F_TOA(lam) cos(chi) exp(-m [tau_R + N_O3 sigma_O3(lam, T) + N_O2 sigma_O2(lam) + sum_j N_j sigma_j(lam)])
    UVI = 40 m2/W * integral S_er(lam) E(lam) dlam      (2.6)
    D_DNA = integral S_Setlow(lam) (lam/300 nm) E(lam) dlam      (energy basis, 2.6)

with m(chi) from the Chapman function (section 3.5), O3 column from part 05 or the Kozakis fit in the earlier report (1 DU = 2.687e16 cm^-2), diffuse sky from a two-stream or the crude half of the Rayleigh-scattered light of 2.7 (+-50 percent). Direct UVI follows DU^-1.08 at overhead sun, 100 to 500 DU (2.7).

### 3.5 Sky colour chain

1. Air mass: Ch(X, chi) from the integral of exp(-(r - R)/H) along the ray (2.8), X = R/H = R m_bar g/(k T); sec(chi) only below 75 degrees.
2. Single-scatter sky radiance seen from the surface at view cosine mu for a sun cosine mu0 (derived here, checked by Monte Carlo): I(mu, Theta) = F0 p(Theta) mu0 [exp(-tau/mu0) - exp(-tau/mu)]/(4 pi (mu0 - mu)); and from space (upward) I_up(mu) = F0 p(Theta) mu0 [1 - exp(-tau (1/mu0 + 1/mu))]/(4 pi (mu0 + mu)). Integrated over the hemisphere these give the first-order diffuse flux. Against the Monte Carlo (COMPUTED, unpolarised, rho = 0) the single-scatter fluxes fall short by 10 percent (tau 0.05, mu0 1), 17 percent (tau 0.1, mu0 0.5) and 21 percent (tau 0.1, mu0 0.2), the contribution of the second and later scatterings, which grows about as tau squared. Earth at 443 nm has tau_R = 0.23 and at 350 nm 0.62, so the blue sky needs the correction.
3. Thick air: Eddington transmission T(tau, mu0) = [(2/3 + mu0) + (2/3 - mu0) exp(-tau/mu0)]/(4/3 + tau); R = 1 - T; spherical albedo s = integral 2 mu (1 - T) dmu; with ground albedo A: R_tot = R + T A (1 - s)/(1 - A s); ground flux = T/(1 - A s) (units of F0 mu0). For the planet seen from space at view cosine mu add the ground-reflected radiance A E_down exp(-tau/mu)/pi with E_down = F0 mu0 T/(1 - A s) (the ground term of R_tot in the 36-case check above); for a sky seen from the ground the albedo term enters only through the extra downward flux T/(1 - A s). Accuracy against Monte Carlo: within 1.3 percent for mu0 of 0.5 and 1, up to 3.5 percent low at mu0 = 0.2 (2.8). Use for tau above 0.3; for tau below 0.3 use the single-scatter radiance plus the second-order excess above, or use two-stream everywhere for the planet disk colour (GEN.418, 419).
4. Transmission of the direct beam at twilight: exp(-m tau_total(lam)) with the Chapman m, including the O3 Chappuis band (peak 5.1e-21 to 5.2e-21 cm2; 300 DU gives tau = 0.042 vertically at 600 nm, COMPUTED, so tau = 0.8 at m = 20, removing orange-red light).

### 3.6 Star to sRGB (absolute colour, GEN.416)

1. Spectrum S(lam): Planck B(lam, T) = 1/[lam^5 (exp(c2/(lam T)) - 1)], c2 = 1.438776877e-2 m K, for Teff above 4500 K; BT-Settl chromaticity from the table in 2.9 below that.
2. Multiply by the atmosphere's transmission or reflectance spectrum (for the sky: single-scatter or two-stream radiance per wavelength).
3. X, Y, Z = integral S(lam) xbar, ybar, zbar dlam (CIE 1931 2-degree tables at 1 nm; trapezoid); x = X/(X + Y + Z), y = Y/(X + Y + Z).
4. Linear sRGB = M [X, Y, Z] with M = [[3.2406, -1.5372, -0.4986], [-0.9689, 1.8758, 0.0415], [0.0557, -0.2040, 1.0570]] (D65); scale so the brightest channel is 1 for a hue-only swatch, clip negatives; encode with 12.92 v below 0.0031308 else 1.055 v^(1/2.4) - 0.055.
5. Table in 2.9 gives checks at 5772 K (x 0.3266, y 0.3358, 8-bit 255, 241, 234) and 6504 K (0.3135, 0.3236).

### 3.7 Rule for alkali and oxides on hot giants (GEN.419)

If T at 1 mbar to 0.1 bar exceeds the Ti condensation curve (about 1700 K at 1 mbar per Fortney 2008), use pM class: strong TiO/VO bands in 0.4 to 1.1 um (colour: dark, red-brown); between roughly 1400 and 1700 K use mixed; below, Na D and K I wings (S = 2e-14 and 3.5e-14 cm2 nm, section 2.5) set the 0.4 to 1.0 um absorption, with cloud and haze optical depths from part 06. Numerical cross sections for the molecules are NOT REACHED (section 4).

---

## 4. Boss must research

Ranked by how much the sky/colour/UV chain depends on it.

1. **Bogaard, Buckingham, Pierens and White 1978, J. Chem. Soc. Faraday Trans. 1, 74, 3008** (via Ranjan and Sasselov 2017 reference list): measured polarizability and depolarization ratio of SO2, H2S and other gases at 488, 514.5 and 632.8 nm. It would replace my DFT depolarization (unverified for NH3, H2S, SO2, HCN, C2H2, C2H4, C2H6, OCS, HCl, H2CO; O3 off by a factor 3). Pays off: the 2.1b table.
2. **Direct Rayleigh or refractivity measurements** of NH3 (King factor), CO, H2S, SO2, Ne, Kr, Xe: Sneep and Ubachs 2005 (JQSRT 92, 293), Bates 1984 (Planet. Space Sci. 32, 785); both behind ScienceDirect. Pays off: removes the +-5 percent generic error and the CH4/CO/NH3 spreads.
3. **Absorption by hot giants:** numerical line lists for TiO, VO, FeH, CrH, Na, K, Li, Rb, Cs (candidate names Toto, VOMYT, MoLLIST, VALD, RECALLED) (the ExoMol site is reachable but the large line lists were not downloaded), plus the Allard et al. Na-H2 and K-H2 pressure-broadened line profiles and Burrows and Volobuyev 2003 tables of widths; resolve the FeH/CrH 0.9896 to 0.9969 um head. Pays off: numeric opacity for the pM/pL colour, and the titanium-cloud pressure boundaries (factor 2 to 3 spread, 2.5c).
4. **Schumann-Runge band transmission:** a published parameterisation (opacity distribution or effective cross section for 175 to 205 nm; a Minschwaner et al. type scheme is a candidate, RECALLED) because band means fail for thick O2 columns (2.3). Pays off: correct UV-C at the surface and mesospheric O3 production in the ozone lane (GEN.413).
5. **Licences and provenance** for the embedded tables: per-dataset licence of Heays 2017 (Leiden database), Huebner 2015 and of the Serdyuchenko and Mainz files (Mainz asks for a citation of Keller-Rudek et al. 2013, ESSD 5, 365; Bremen gives no licence text). The embedded 35-bin table should be regenerated from a source with a clear licence or the licences confirmed.
6. **Visible absorbers for giants and haze worlds** (earlier report item 3): Karkoschka 1994 (Icarus 111, 174; the Photochem CH4 data above 152 nm rest on it); NH3 visible bands; SO2 and NO2 visible tails; S2, S3, S4, S8, SO3, H2SO4, HNO3 and CS2 cross sections exist as Photochem files and were not binned.
7. **Species with no table here:** HF, CH3Cl, SiO, SO, CS2 (Mainz pages exist for HF, CS2, SO, not read), perchlorates, carbonates, radon, atomic Na and K below 100 nm.
8. **Stellar UV for non-solar stars:** measured UV spectra of M, K and F dwarfs (MUSCLES or similar) and flare-time UV; the only measured spectrum used here is the ATLAS-1 Sun.
9. **O3 Huggins temperature fit** (quadratic coefficients 300 to 370 nm, a file on the Bremen page) and O2 and H2O high-temperature UV cross sections (Mainz inventory lists them, values not extracted), needed for hot and lava worlds.
10. **CIE 1998 / ISO 17166 erythemal text and Setlow 1974 (PNAS 71, 3363)** themselves, to confirm the exponent and the normalisation puzzle in 2.6.
11. **Bucholtz 1995 Table 3, Penndorf 1957, Hulburt 1953, Kasten and Young 1989, Kasting and Donahue 1980**: still behind paywalls or captchas (2.10).
12. A measured **depolarization of N2** that settles 0.0199 against 0.0305 (a 2 percent sigma issue).

---

## 5. Recalled / to verify

Every value below is RECALLED (not read in a source) and appears in the text with that label.

| Item | Where used | What to verify |
|---|---|---|
| Starting bond lengths and angles for the DFT geometries | 2.2 method | irrelevant to the result because the geometries were optimized, but not independently checked |
| TiO band heads near 0.705, 0.76 and 0.84 um | 2.5b | spectroscopic constants from a TiO line list |
| Methane bands near 619, 725 and 890 nm | 2.3 notes | Karkoschka 1994 |
| Radiation amplification factor for erythema near 1.1 | 2.7 | UNEP or WMO ozone assessment |
| Minschwaner et al. type Schumann-Runge parameterisation as a candidate | section 4 item 4 | exact reference |
| Names of the line lists (Toto, VOMYT, MoLLIST) and the Allard profile papers | section 4 item 3 | ExoMol site |
| Leiden VUV database as the common name of the Heays 2017 data | 2.3 | A&A 602, A105 |

SNIPPET-level items were not used.

---

## 6. Cross-links to neighbouring parts

- **Part 01 (`01-gas-table-thermodynamics.md`):** add alpha (CCCBDB) and rho_n (2.1b, 2.2) as gas-table columns; its photolysis-threshold column should agree with the cut-offs in 2.3.
- **Part 03 (`03-infrared-greenhouse-thermal.md`, section 7.1 and 7.2):** it shares the Rayleigh table and the conservative-layer albedo with a two-stream closure. Collision to resolve: part 03 uses F_K = 1.0 for H2 (this part: 1.027 to 1.038, a 3 to 4 percent change in sigma) and its Peck-Khanna N2 and Peck H2 files at 273.15 K give 4.6e-27 and 9.7e-28 at 550 nm, consistent with 2.1a. The Eddington accuracy figures of 2.8 apply to its albedo layer. IR bands and collision-induced absorption are there, not here.
- **Part 04 (`04-photochemistry-lifetimes-escape.md`, sections 2.3, 2.5 and 3):** photolysis rates, stellar UV bands, the UV-shielding band means (its 2.5 table agrees in kind with 2.3: O3 7.8e-18 at 242 to 290 nm against my 1.1e-17 and 7.7e-18 in the 240 to 280 nm bins) and the Chapman balance use the cross sections of 2.3 and the formulas of 3.3. The band-mean caveat of 2.3 applies to its tau = 1 mixing ratios.
- **Part 05 (`05-mixing-ratios-by-class.md`):** gas columns for tau_R and the UV absorbers; the Kozakis ozone fit of the earlier report.
- **Part 06 (`06-condensates-clouds-hazes-dust-minerals.md`):** aerosol optical depth, albedo and Angstrom exponent add to tau_R (3.5); Na2S and KCl condensation and the Ti-bearing condensate pressure boundaries (2.5c) belong with its condensation curves.
- **Part 07 (`07-radiation-shielding-and-ground-sources.md`):** it points to this part for the code's UV index and for ozone cross sections. Read-only check of the code (`physics/radiation.py`): `uv_index` = (blackbody 200 to 300 nm fraction relative to the Sun, floor 0.1) x flux ratio x ozone_ratio^-1.6, with ozone_ratio = sqrt(pO2/21.2 kPa) capped at 1.2 and a cap of 1000 for no ozone. Against 2.7 (COMPUTED, fixed Sun): erythemal UVI rises 18 (direct) to 27 (crude global) times from 300 DU to zero ozone, not 1000, while the Setlow-weighted dose rises about 280 times and the UV-C from 1e-16 to 0.7 W m^-2; the direct UVI exponent from 100 to 500 DU is -1.08. The code number is therefore neither the CIE UV index nor a DNA dose; GEN.413 and the radiation lane should define which quantity is stored.
- **Chain items:** GEN.413 (ozone), GEN.416 (blackbody colour: 2.9, 3.6), GEN.417 (Rayleigh optical depth of the mixture: 2.1, 3.1), GEN.418 and GEN.419 (solid and giant atmosphere colour: 3.5, 3.7), GEN.414 (aerosol depths), DB.68 (ozone and aerosol columns).
