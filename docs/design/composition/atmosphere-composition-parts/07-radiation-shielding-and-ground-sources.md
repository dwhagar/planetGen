# 07 Radiation: atmospheric shielding and ground sources

Part 07 of the atmosphere series. Scope: derive the ATMOSPHERIC radiation environment (what the gas column stops and makes) and the SURFACE radiation environment (atmosphere plus ground composition) from composition, pressure, gravity, magnetic field and stellar activity. Earth, Moon, Mars and Europa appear only as labelled calibration points; nothing is taken from presets.

Status labels: READ (original paper text), 2ND-HAND (read in a reachable document that cites the original), DATA (data file or static table, source named), COMPUTED (derived by me in Python from stated inputs), SNIPPET (search summary only), RECALLED (memory, also listed in section 7), NOT REACHED. A combined label such as "COMPUTED (DATA)" counts under its first word. Abbreviations: X = column mass P0/g in g/cm2; M = dipole moment in units of Earth's (7.768e22 A m2, project constant); mu/rho = mass attenuation coefficient in cm2/g; CSDA = continuous slowing down approximation; lambda_I = PDG nuclear interaction length; X0 = radiation length.

## 1. Scope and how the numbers feed a generator

Inputs the generator already has per body: surface pressure P0, gravity g, star spectrum class and age, flare rate, magnetic dipole moment and radius, land fraction, O2 partial pressure. New inputs this part needs from the composition work (parts 01, 05): mole fractions of the 12 listed gases plus H2O, and a bulk-silicate composition that gives K, U, Th (ppm), or a stellar Th/Eu scaling. Outputs: (a) XUV/X-ray optical depth of the column, which sets how much stellar high-energy light reaches the mesosphere or the ground; (b) proton stopping energies and ranges, which set the SEP threshold; (c) GCR dose versus X and dipole cutoff; (d) SEP dose per event; (e) ground gamma, radon and radiolysis rates; (f) ion-pair production profile for NOx/HOx chemistry (part 04); (g) a dose in mSv/yr to compare with the limits table.

Key findings that change the current code (details in sections 2 to 4):

1. The dipole cutoff scales as M/r^2, not M/r^3 (COMPUTED); the code uses (R_E/R_p)^3.
2. The GCR dose-versus-column curve is 1.5 to 5 times below the proton-primary CORSIKA results (READ) and the aviation measurement (READ) between X = 100 and 1033 g/cm2. The flat magnetic factor 0.6 is wrong for thin atmospheres: the CORSIKA ratio for an Earth-strength field is 0.31 at X = 100 and 0.52 at X = 1036 (READ).
3. The range law R = 0.0022 E^1.77 g/cm2 is within 13% of PSTAR from 10 to 300 MeV but overstates the range 1.2 times at 1 GeV and 5 times at 10 GeV (COMPUTED vs DATA). It is adequate for SEP thresholds, poor for GCR secondaries.
4. A hydrogen-helium envelope absorbs EUV about as well per gram as N2 or CO2 below 36 nm (COMPUTED). Per molecule it is 3 to 6 times more transparent.
5. The code's 0.48 mSv/yr ground term equals the UNSCEAR 2008 occupancy-weighted terrestrial external dose (READ) and its 0.5 mSv/yr radon term is a full-time outdoor value, while UNSCEAR gives 0.095 outdoors and 1.0 indoors (READ).
6. Composition changes the dose through lambda_I and X0; [R0187] simulated proton primaries in C/N/O air only and found a few percent effect (READ). H2/He columns are unstudied (NOT REACHED).

## 2. Tables

### 2.1 Current code constants against the evidence

| Code constant (radiation.py) | Value | Evidence found | Verdict | Source | Status |
|---|---|---|---|---|---|
| GCR_POINTS (X, mSv/yr) | (0, 500), (16.4, 237), (100, 110) | Moon 1369 uSv/day = 500 mSv/yr (solar-minimum upper bound); Mars 22.5 g/cm2 about 0.7 mSv/day = 255 mSv/yr (range 128 to 420); X = 100 CORSIKA 171 (M=1) to 553 (M=0) | Moon and Mars points sound; X = 100 point is 2.6 (with field) to 5 (no field) times below proton-primary CORSIKA | [R4884] PMC7518862; [R1730] 1507.03473; [R0187] 1307.4704 | READ |
| GCR_ATTENUATION_G_CM2 | 164 | CORSIKA e-folding: 124 (X 200 to 500), 165 (500 to 1036) for M=1 | Within the spread above 500, short below | [R0187] Table XV, my fit | COMPUTED (READ) |
| magnetic_factor, floor 0.6 | 0.6 at Rc >= 15 GV | Earth-strength field ratio to no-field: 0.31 (X=100), 0.37 (X=200), 0.51 (X=500), 0.61 (X=700), 0.52 (X=1036) | Factor depends on X; 0.6 is right only near sea level | [R0187] | READ |
| EARTH_CUTOFF_GV | 14.9 | Vertical dipole cutoff c (mu0/4 pi) M / (4 r^2) = 14.34 GV for M = 7.768e22, r = 6371 km; 14.9 implies M = 8.1e22 | 4% high, acceptable | formula form RECALLED, constants evaluated by me | COMPUTED (RECALLED) |
| cutoff scaling | (R_E/R_p)^3 | Stormer: Rc proportional to M/r^2 | Wrong exponent | F1 | COMPUTED (RECALLED) |
| RANGE_COEFFICIENT, EXPONENT | 0.0022, 1.77 | PSTAR air: ratio 0.92 (10 MeV), 0.87 (100), 0.92 (300), 1.22 (1 GeV), 2.16 (3 GeV), 5.2 (10 GeV) | Good below 300 MeV | NIST PSTAR, my comparison | COMPUTED (DATA) |
| SEP_R0_GV | 0.2 | Soft event value; hard 1956-type events are harder | Needs spread | no source | NOT REACHED |
| SEP_FREE_MSV_YR | 20 | No event-level fluence mapping in code | Replace by fluence law F4 | [R0183] 1606.07027 | READ |
| GROUND_MSV_YR | 0.48 | UNSCEAR 2008 terrestrial external 0.48 mSv/yr (0.07 outdoors, 0.41 indoors) | Equals indoor-weighted world mean | UNSCEAR 2008 Annex B Table 12 | READ |
| RADON_MSV_YR | 0.5 | UNSCEAR: radon inhalation 1.15 total, indoor 1.0, outdoor 0.095 | 0.5 equals a full-time outdoor value (COMPUTED 0.47) | UNSCEAR 2008 | READ |
| RADON_MIN_PRESSURE_PA | 1000 | No source found | Needs source | none | NOT REACHED |

### 2.2 Photon absorption: XUV and EUV (10 to 100 nm) photoionisation cross sections

Mb = 1e-18 cm2. Mean over a window of plus or minus 3% about each wavelength (smooths H2 and N2 line structure). Origin per cell: L = Leiden measured or evaluated file, H = Henke atomic sum sigma = 2 r0 lambda f2, V = [R4542] fit. Henke fails beyond about 70 nm (no outer-shell physics), so H entries at 70 nm or above are not given. Henke and Leiden agree within about 10 to 30% below 70 nm; Henke and Verner agree within 3 to 10% for He.

| Gas | 5 nm | 10 nm | 20 nm | 30 nm | 50 nm | 70 nm | 90 nm | Source | Status |
|---|---|---|---|---|---|---|---|---|---|
| H2 | 0.0022 H | 0.0199 H | 0.255 L | 0.978 L | 4.4 L | 10.1 L | 13.9 L | Leiden cross_sections/H2 [R1873]; Henke f2 (xraydb 4.5.8) | DATA |
| He | 0.0288 V | 0.22 V | 1.27 V | 2.97 V | 4.43 V | 0 | 0 | [R4542] astro-ph/9601009; edge 50.4 nm | COMPUTED (READ) |
| N2 | 0.234 H | 1.27 H | 6.2 L | 11.5 L | 24.5 L | 25.5 L | 30 L | Leiden N2 (range 11.6 to 100 nm); Henke | DATA |
| O2 | 0.575 L | 2.43 L | 9.37 L | 16.9 L | 24.4 L | 25.3 L | 9.4 L | Leiden O2 (4.2 to 240 nm) | DATA |
| CO2 | 0.478 H | 2.57 H | 9.7 H | 19 H | 30.4 L | 35.2 L | 37.9 L | Leiden CO2 (34.6 to 201 nm); Henke below | DATA |
| H2O | 0.213 H | 1.08 L | 4.81 L | 9.59 L | 19.1 L | 22.2 L | 21.1 L | Leiden H2O (6.2 to 207 nm) | DATA |
| CH4 | 0.0609 H | 0.321 L | 1.91 L | 5.92 L | 22.6 L | 39.3 L | 48.4 L | Leiden CH4 (5.7 to 154 nm) | DATA |
| NH3 | 0.12 H | 0.581 L | 3.19 L | 8.1 L | 21.8 L | 31.6 L | 20.9 L | Leiden NH3 (8 to 226 nm) | DATA |
| Ar | 2.36 H | 1.05 H | 1.17 H | 3.08 H | 34.9 H | 36.2 H | 0 | Henke only (K edge 3.2 keV, L edge near 250 eV, valence edge 78.7 nm); outer-shell data NOT REACHED | COMPUTED |
| SO2 | 2.12 L | 2.88 L | 9.97 L | 16.5 L | 40.8 L | 61.7 L | 46.1 L | Leiden SO2 (1.7 to 417 nm) | DATA |
| CO | 0.267 H | 1.38 L | 6.25 L | 11.8 L | 21 L | 25.6 L | 48.8 L | Leiden CO (6.2 to 167 nm) | DATA |
| O3 | 0.632 H | 3.39 H | 12.5 H | 23.9 H | 37.2 H | 34.3 L | 13.9 L | Leiden O3 (52.6 to 1100 nm); Henke below | DATA |

Not tabulated: H2S (Leiden file exists, 41.6 to 238 nm). No cross-section source reached for HCN, H2CO, HCl, HF, CH3Cl, Ne, Kr, Xe, SiO, Na, K, TiO, VO, NOx (NOT REACHED); the Henke atomic sum covers their 5 to 60 nm part.

Ionisation edges and thresholds ([R4542] fits, COMPUTED with the printed parameters; He I parameters E0 = 13.61 eV, s0 = 949.2 Mb, ya = 1.469, P = 3.188, yw = 2.039, y0 = 0.4434, y1 = 2.136):

| Edge | Wavelength (nm) | Energy (eV) | Threshold sigma (Mb) | tau = 1 column of that absorber (cm-2) | tau = 1 mass (g/cm2) | Source | Status |
|---|---|---|---|---|---|---|---|
| H I Lyman continuum | 91.2 | 13.60 | 6.35 | 1.6e17 | 2.6e-7 | [R4542] | COMPUTED (READ) |
| He I | 50.4 | 24.59 | 7.43 | 1.35e17 | 9.0e-7 | [R4542] | COMPUTED (READ) |
| He II | 22.8 | 54.42 | 1.59 | 6.3e17 | 4.2e-6 | [R4542] (hydrogenic scaling 6.3 Mb / 4) | COMPUTED (READ) |

### 2.3 Tau = 1 columns for gas mixtures in the XUV and EUV

Column N(tau = 1) in molecules/cm2, then g/cm2. H/He is a solar-composition envelope (x_He = 0.145 by number, from [R0169] log eps_He = 10.93, READ in 0909.0948). Earth-like air is 78% N2, 21% O2, 1% Ar. Venus-like CO2+N2, steam and Titan-like cases were computed and behave like their main component (not shown). Origin: sigma table in 2.2.

| Mixture | 10 nm | 22.5 nm | 30 nm | 50 nm | 90 nm | Source | Status |
|---|---|---|---|---|---|---|---|
| H2/He envelope | 2.0e19 / 7.8e-5 | 1.8e18 / 6.7e-6 | 7.9e17 / 3.0e-6 | 2.3e17 / 8.7e-7 | 8.4e16 / 3.2e-7 | table 2.2 weighted by mole fraction | COMPUTED (DATA) |
| N2 | 7.9e17 / 3.7e-5 | 1.2e17 / 5.6e-6 | 8.7e16 / 4.0e-6 | 4.1e16 / 1.9e-6 | 3.3e16 / 1.5e-6 | same | COMPUTED (DATA) |
| CO2 | 3.9e17 / 2.8e-5 | 8.4e16 / 6.1e-6 | 5.3e16 / 3.9e-6 | 3.3e16 / 2.4e-6 | 2.6e16 / 1.9e-6 | same | COMPUTED (DATA) |
| Earth-like air | 6.6e17 / 3.2e-5 | 1.1e17 / 5.4e-6 | 8.0e16 / 3.8e-6 | 4.1e16 / 2.0e-6 | 3.9e16 / 1.9e-6 | same | COMPUTED (DATA) |

Band means (molecules/cm2 / g/cm2) for the H/He, N2 and CO2 comparison:

| Mixture | 0.5 to 10 nm | 10 to 36 nm | 36 to 50.4 nm | 50.4 to 91.2 nm | Source | Status |
|---|---|---|---|---|---|---|
| H2/He envelope | 7.7e19 / 2.9e-4 | 1.3e18 / 5.0e-6 | 2.9e17 / 1.1e-6 | 1.1e17 / 4.3e-7 | band mean of 2.2 | COMPUTED (DATA) |
| N2 | 1.9e18 / 8.9e-5 | 1.2e17 / 5.8e-6 | 4.6e16 / 2.2e-6 | 3.9e16 / 1.8e-6 | same | COMPUTED (DATA) |
| CO2 | 9.6e17 / 7.1e-5 | 7.7e16 / 5.6e-6 | 3.5e16 / 2.6e-6 | 3.2e16 / 2.4e-6 | same | COMPUTED (DATA) |

Reading: every case reaches tau = 1 at 1e-6 to 1e-4 g/cm2, which is about 1e-9 to 1e-7 of an Earth atmosphere. EUV never reaches the ground of any body with a real atmosphere. It is absorbed in the thermosphere, so the relevant output for the generator is the energy deposition altitude (part 04, escape), not a surface dose. [R4127] (astro-ph/0308311): atmospheres thicker than about 100 g/cm2 block nearly all X-rays but transmit and reprocess gamma rays; thick atmospheres turn absorbed energy into UV and visible re-emission (READ).

### 2.4 X-ray and gamma mass attenuation coefficients mu/rho (cm2/g), 0.1 keV to 10 MeV

Mixture rule mu/rho = sum of mass fraction times elemental mu/rho. Elements from NIST XCOM (total with coherent scattering, DATA, 1 keV to 20 MeV); below 1 keV from the Henke photoabsorption sum (xraydb 4.5.8). Molecular binding and near-edge structure are ignored. O3 equals O2 by construction (row omitted). Values within 100 eV to 300 eV of the C, N and O K edges (284, 410, 543 eV) are edge-sensitive.

| Material | 0.1 keV | 0.3 keV | 1 keV | 3 keV | 10 keV | 30 keV | 100 keV | 0.3 MeV | 1 MeV | 3 MeV | 10 MeV | Status |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| H2 | 1.15e4 | 358 | 7.22 | 0.561 | 0.385 | 0.357 | 0.294 | 0.211 | 0.126 | 0.0692 | 0.0325 | COMPUTED (DATA) |
| He | 6.5e4 | 2.61e3 | 60.8 | 2.01 | 0.248 | 0.184 | 0.149 | 0.106 | 0.0636 | 0.035 | 0.017 | COMPUTED (DATA) |
| N2 | 4.35e4 | 3.06e3 | 3.31e3 | 146 | 3.88 | 0.307 | 0.153 | 0.107 | 0.0636 | 0.0358 | 0.0202 | COMPUTED (DATA) |
| O2 | 6.79e4 | 4.95e3 | 4.59e3 | 217 | 5.95 | 0.378 | 0.155 | 0.107 | 0.0637 | 0.036 | 0.0209 | COMPUTED (DATA) |
| CO2 | 5.6e4 | 1.61e4 | 3.94e3 | 183 | 4.98 | 0.345 | 0.154 | 0.107 | 0.0637 | 0.0359 | 0.0205 | COMPUTED (DATA) |
| H2O | 6.16e4 | 4.43e3 | 4.08e3 | 193 | 5.33 | 0.376 | 0.171 | 0.119 | 0.0707 | 0.0397 | 0.0222 | COMPUTED (DATA) |
| CH4 | 2.12e4 | 3.44e4 | 1.66e3 | 67.8 | 1.87 | 0.282 | 0.187 | 0.133 | 0.0794 | 0.0441 | 0.0228 | COMPUTED (DATA) |
| NH3 | 3.79e4 | 2.58e3 | 2.72e3 | 120 | 3.26 | 0.316 | 0.178 | 0.125 | 0.0748 | 0.0417 | 0.0224 | COMPUTED (DATA) |
| Ar | 1.99e4 | 4.67e4 | 3.18e3 | 170 | 63.2 | 2.7 | 0.204 | 0.0995 | 0.0576 | 0.0338 | 0.0245 | COMPUTED (DATA) |
| SO2 | 4.12e4 | 2.31e4 | 3.51e3 | 779 | 28.1 | 1.25 | 0.179 | 0.108 | 0.0637 | 0.0366 | 0.0234 | COMPUTED (DATA) |
| CO | 4.92e4 | 2.25e4 | 3.57e3 | 163 | 4.42 | 0.326 | 0.154 | 0.107 | 0.0637 | 0.0358 | 0.0203 | COMPUTED (DATA) |
| SiO2 (regolith) | 5.11e4 | 1.77e4 | 3.18e3 | 573 | 19 | 0.872 | 0.168 | 0.108 | 0.0637 | 0.0363 | 0.0226 | COMPUTED (DATA) |
| Fe | 5.78e4 | 1.4e4 | 9.09e3 | 558 | 171 | 8.18 | 0.372 | 0.11 | 0.06 | 0.0362 | 0.0299 | COMPUTED (DATA) |

Henke photoabsorption alone against XCOM total (percent, at 1, 3, 10, 30 keV): N -0.1, -1.4, -11.7, -69; O +0.6, +1.6, -6.8, -58; C +0.3, +0.1, -13, -80; Ar +5.1, +3.2, +0.8, -8.1; H -6, -70, -99, -100 (COMPUTED). Henke is safe to about 3 keV for C, N, O and Ar; use XCOM above.

Vertical tau = 1 mass columns 1/(mu/rho) in g/cm2 (COMPUTED from the table): N2 3.0e-4 (1 keV), 0.26 (10 keV), 6.5 (100 keV), 15.7 (1 MeV), 49 (10 MeV); CO2 2.5e-4, 0.20, 6.5, 15.7, 49; H2O 2.5e-4, 0.19, 5.9, 14.1, 45; H2 0.139, 2.6, 3.4, 7.9, 31; SiO2 3.2e-4, 0.053, 5.9, 15.7, 44. A 1 bar atmosphere (1033 g/cm2) is 66 optical depths thick at 1 MeV, so ground gamma and X-ray dose is secondary production plus ground-sourced gamma, never primary stellar photons. At 100 keV to 1 MeV, H2 absorbs and scatters about 2 times more per gram than N2 or CO2 (mu/rho 0.294 against 0.153 at 100 keV, 0.126 against 0.064 at 1 MeV), because Compton scattering scales with Z/A. Below 10 keV, H2 and He are 100 to 1000 times more transparent per gram.

### 2.5 Proton stopping power and range (NIST PSTAR, ICRU 49)

PSTAR electronic plus nuclear stopping, CSDA range, 10 MeV to 10 GeV (20 energies), pulled for 34 materials. Nuclear interactions are not in CSDA, so the practical stopping length for the primary beam is about lambda_I (section 2.6). I values (eV, ICRU 49, DATA): H2 19.2, He 41.8, N2 82.0, O2 95.0, Ne 137, Ar 188, Kr 352, Xe 482, CO2 85.0, CH4 41.7, air 85.7, SiO2 139.2, liquid water 75.0, Al2O3 145.2, Fe 286, graphite 78, polyethylene 57.4. PDG gives I = 79.7 eV for water and 136.4 eV for standard rock against PSTAR 75.0 and (SiO2) 139.2: a spread of 6% in I changes the stopping power by about 1%.

CSDA range, g/cm2 (DATA, physics.nist.gov/PhysRefData/Star):

| Material | 10 MeV | 30 MeV | 100 MeV | 300 MeV | 1 GeV | 3 GeV | 10 GeV | Source | Status |
|---|---|---|---|---|---|---|---|---|---|
| H2 | 0.05346 | 0.4035 | 3.633 | 24.69 | 158.7 | 637.7 | 2254 | PSTAR matno 001 | DATA |
| He | 0.1211 | 0.8945 | 7.929 | 53.36 | 340 | 1356 | 4757 | PSTAR 002 | DATA |
| N2 | 0.1393 | 0.9982 | 8.68 | 57.78 | 365 | 1446 | 5030 | PSTAR 007 | DATA |
| O2 | 0.1439 | 1.024 | 8.859 | 58.81 | 370.8 | 1466 | 5091 | PSTAR 008 | DATA |
| CO2 | 0.1404 | 1.004 | 8.722 | 58.02 | 366.3 | 1450 | 5044 | PSTAR 134 | DATA |
| CH4 | 0.09747 | 0.7178 | 6.356 | 42.76 | 272.5 | 1087 | 3812 | PSTAR 197 | DATA |
| Air | 0.1408 | 1.007 | 8.744 | 58.16 | 367.2 | 1454 | 5054 | PSTAR 104 | DATA |
| Liquid water (ice) | 0.123 | 0.8853 | 7.718 | 51.45 | 325.4 | 1304 | 4700 | PSTAR 276 | DATA |
| SiO2 (regolith) | 0.1573 | 1.099 | 9.39 | 61.87 | 388.1 | 1545 | 5493 | PSTAR 245 | DATA |
| Fe | 0.2064 | 1.377 | 11.37 | 73.71 | 459.6 | 1818 | 6389 | PSTAR 026 | DATA |

Electronic stopping powers at the same energies are in the scratch data (DATA) and follow from F2 within 1% above 30 MeV. PDG minimum-ionisation values (MeV cm2/g): H2 4.103, air 1.815, CO2 1.819, rock 1.688 (DATA, pdg.lbl.gov AtomicNuclearProperties).

Bethe-Bloch check (COMPUTED, formula F2 with PSTAR I values, no shell, Barkas, Bloch or density-effect terms): ratio PSTAR / simple Bethe at 10, 100, 1000, 10000 MeV: air 0.992, 0.999, 1.000, 1.000; H2 1.000 at all four; CO2 0.992, 0.999, 1.000, 1.000; Ar 0.980, 0.998, 1.000, 1.000; liquid water 0.993, 0.999, 0.999, 0.922; SiO2 0.988, 0.999, 0.996, 0.931; Fe 0.970, 0.994, 0.980, 0.920. Range ratio PSTAR / integrated Bethe (starting from PSTAR range at 10 MeV): within 0.4% for gases at all energies, within 6% for condensed matter at 10 GeV (the missing density-effect term). So for gases the plain Bethe formula is enough above 10 MeV; for ice, water, rock and iron above 3 GeV the density correction must be added (shell and Barkas corrections matter below 10 MeV, outside this table).

Power-law fits R = a E^b (g/cm2, MeV), COMPUTED from PSTAR, maximum error against the table:

| Material | a, b (10 to 300 MeV) | max error | a, b (300 MeV to 10 GeV) | max error | Status |
|---|---|---|---|---|---|
| Air (N2, O2, CO2 within 2%) | 0.002371, 1.779 | 4.0% | 0.05512, 1.259 | 24.9% | COMPUTED (DATA) |
| H2 | 0.0008476, 1.811 | 5.1% | 0.02164, 1.273 | 25.0% | COMPUTED (DATA) |
| He | 0.001971, 1.797 | 4.7% | 0.04855, 1.267 | 24.9% | COMPUTED (DATA) |
| CH4 | 0.001589, 1.796 | 4.6% | 0.0389, 1.267 | 25.0% | COMPUTED (DATA) |
| Water | 0.002058, 1.783 | 4.1% | 0.04459, 1.273 | 23.6% | COMPUTED (DATA) |
| SiO2 | 0.002721, 1.765 | 3.6% | 0.05609, 1.266 | 23.7% | COMPUTED (DATA) |
| Fe | 0.003767, 1.737 | 2.5% | 0.06946, 1.259 | 23.9% | COMPUTED (DATA) |

The code's pair (0.0022, 1.77) is the 10 to 300 MeV fit for N2 or CO2 to within 7%; it is 2.1 times too long for H2 across that range, because H2 has a larger electron-to-mass ratio (range per gram is 2.6 times shorter than N2 at 10 MeV). Ratio of the code's range to PSTAR at 10, 30, 100, 300, 1000, 3000, 10000 MeV: air 0.92, 0.90, 0.87, 0.92, 1.22, 2.16, 5.23; H2 2.42, 2.24, 2.10, 2.16, 2.83, 4.92, 11.7; water 1.05, 1.02, 0.99, 1.04, 1.38, 2.41, 5.63; SiO2 0.82, 0.82, 0.81, 0.86, 1.16, 2.03, 4.82 (COMPUTED).

Proton kinetic energy (MeV) needed to reach a vertical column X (CSDA, COMPUTED by inverting PSTAR; rigidity in GV below):

| X (g/cm2) | Air | CO2 | H2 | He | Water | SiO2 | Code's proton_cutoff (as T, MeV) | Status |
|---|---|---|---|---|---|---|---|---|
| 1 | 30 | 30 | 49 | 32 | 32 | 28 | 32 | COMPUTED (DATA) |
| 3 | 55 | 55 | 90 | 58 | 59 | 53 | 59 | COMPUTED (DATA) |
| 10 | 108 | 108 | 177 | 114 | 116 | 104 | 117 | COMPUTED (DATA) |
| 22.5 | 171 | 172 | 284 | 181 | 184 | 165 | 184 | COMPUTED (DATA) |
| 50 | 274 | 275 | 463 | 289 | 295 | 264 | 289 | COMPUTED (DATA) |
| 100 | 420 | 421 | 727 | 443 | 453 | 405 | 428 | COMPUTED (DATA) |
| 300 | 868 | 870 | 1610 | 917 | 945 | 835 | 796 | COMPUTED (DATA) |
| 1033 | 2241 | 2245 | 4656 | 2381 | 2466 | 2133 | 1601 | COMPUTED (DATA) |

Proton kinetic energy to rigidity: 10, 30, 100, 300, 1000, 3000, 10000 MeV = 0.137, 0.239, 0.445, 0.808, 1.696, 3.825, 10.90 GV (COMPUTED). Mars check: PSTAR gives 171 MeV for 22.5 g/cm2; [R1728] (1705.06763) finds 140 to 160 MeV with angular spread and secondaries (READ). The last column is the code's own formula (X / 0.0022)^(1/1.77) evaluated at the same columns (COMPUTED): it tracks air within 8% up to 300 g/cm2 and is 29% low in energy at 1033 g/cm2.

Alpha particles: ASTAR tables for the same 34 materials (12 energies, 10 MeV to 1 GeV) are in the scratch data (DATA). A He-4 ion at four times a proton's kinetic energy has the same velocity and, since mass 4 divided by charge squared 4 equals 1, about the same range as that proton (Bethe scaling, shell and Barkas terms omitted). The helium-ion shielding length in Mars air is 41.2 g/cm2 ([R1729], 1712.06987, READ).

### 2.6 Hadronic, electromagnetic and muon properties of gases and solids (PDG)

DATA from pdg.lbl.gov AtomicNuclearProperties pages (static HTML). Values for helium and argon gas pages were not found (page-name guesses returned 404; NOT REACHED).

| Material | X0 (g/cm2) | lambda_I (g/cm2) | lambda_T collision (g/cm2) | Critical energy e- (MeV) | Muon critical energy (GeV) | Source | Status |
|---|---|---|---|---|---|---|---|
| H2 gas | 63.04 | 52.0 | 42.8 | 344.8 | not fetched | PDG hydrogen_gas | DATA |
| CH4 | 46.47 | 73.8 | 54.0 | 146.9 | not fetched | PDG methane | DATA |
| NH3 | 40.87 | 79.5 | 56.8 | 121.7 | not fetched | PDG ammonia | DATA |
| Water or ice | 36.08 | 83.3 | 58.5 | 78.6 | 1031 | PDG water_ice | DATA |
| N2 | 37.99 | 89.7 | 61.1 | 91.74 | 1154 | PDG nitrogen_gas | DATA |
| Dry air | 36.62 | 90.1 | 61.3 | 87.92 | 1115 | PDG air_dry_1_atm | DATA |
| CO2 | 36.19 | 88.9 | 60.7 | 86.17 | not fetched | PDG carbon_dioxide_gas | DATA |
| O2 | 34.24 | 90.2 | 61.3 | 81.45 | not fetched | PDG oxygen_gas | DATA |
| Fused quartz (SiO2) | 27.05 | 97.8 | 65.2 | 50.58 | not fetched | PDG silicon_dioxide | DATA |
| Standard rock | 26.54 | 101.3 | 66.8 | 49.13 | 693 | PDG standard_rock | DATA |
| CaCO3 | 24.03 | 99.4 | 66.1 | 44.20 | not fetched | PDG calcium_carbonate | DATA |
| Sulfur | 19.50 | 112.4 | 72.2 | 35.85 | not fetched | PDG sulfur_S | DATA |
| Iron | 13.84 | 132.1 | 81.7 | 21.68 | 347 | PDG iron_Fe | DATA |
| Polyethylene | 44.77 | 78.5 | 56.1 | 101.8 | not fetched | PDG polyethylene | DATA |

Muon CSDA range, g/cm2 (PDG muon tables, DATA; interpolated): air 508 (1 GeV), 4204 (10 GeV), 3.41e4 (100 GeV), 2.31e5 (1 TeV); water ice 472, 4268, 3.63e4, 2.43e5; standard rock 553, 4920, 4.08e4, 2.45e5; iron 640, 5532, 4.31e4, 2.06e5. Rock and ice stop a 1 GeV muon in 2.1 m and 5.1 m respectively with densities 2.65 and 0.92 (COMPUTED). Mean muon ionisation loss near minimum is about 2 MeV cm2/g (PDG a_mu).

Light-gas consequence (design inference F8, not a measurement): a nucleon interacts after 52 g/cm2 of H2 but 90 g/cm2 of air, so hadronic cascades develop at smaller X in H2 or CH4, while electromagnetic cascades need 1.7 times (H2) or 1.27 times (CH4) more g/cm2 (X0 ratio to air); the Pfotzer maximum, an electromagnetic and muon feature, should sit deeper in g/cm2 for light gases.

### 2.7 Dose versus column: simulations, measurements and the code's curve

[R1466] (1307.4704) Table XV, as tabulated by [R1687] (1603.06500) Table 2: total biological dose (mSv/yr) from muons, electrons, photons and neutrons in CORSIKA air showers. Primaries are protons only (20 million per run, READ); hadronic secondaries are modelled for primaries above 8 GeV, so low-column values are lower bounds and helium and heavier nuclei are absent. Planet columns use Earth-like N2/O2 air. Magnetic moment in Earth units; the code's curve is evaluated with its own function (COMPUTED).

| X (g/cm2) | Atri M = 0 | Atri M = 0.15 | Atri M = 0.5 | Atri M = 1 | Code, no field | Code, 14.9 GV field | Measured anchor | Source | Status |
|---|---|---|---|---|---|---|---|---|---|
| 0 | not simulated | not simulated | not simulated | not simulated | 500 | 300 | Moon, Chang'E-4 LND: 1369 uSv/day = 500 mSv/yr, solar-minimum bound | [R4884, R4888] PMC7518862 | READ |
| 22.5 | not simulated | not simulated | not simulated | not simulated | 208 | 125 | Mars Gale: 0.7 mSv/day mean = 255 mSv/yr; 0.35 (solar maximum) to 1.15 (solar minimum) mSv/day; 220 uGy/day absorbed | [R1730] 1507.03473; [R1729] 1712.06987 | READ |
| 100 | 553 | 510 | 257 | 171 | 110 | 66 | none reached | [R0187] | READ |
| 122.7 (Pfotzer, Earth) | 444 (log interpolation) | | | 145 (log interpolation) | 96 | 58 | Pfotzer maximum 15.0 +/- 0.32 km, 122.7 +/- 6.3 g/cm2 | [R3919] 1707.00275 | READ |
| 200 | 263 | 260 | 141 | 97 | 60 | 36 | none reached | [R0187] | READ |
| 265 (10 km) | 114 (log interpolation) | | | 46 (log interpolation) | 40 | 24 | 5.46 uSv/h effective dose = 47.9 mSv/yr (aviation) | [R3826] 2112.07911 | READ |
| 500 | 17.2 | 18.2 | 12.0 | 8.7 | 9.6 | 5.8 | none reached | [R0187] | READ |
| 700 | 3.02 | 3.38 | 2.53 | 1.84 | 2.8 | 1.7 | none reached | [R0187] | READ |
| 1033 to 1036 | 0.65 | 0.46 | 0.42 | 0.34 | 0.37 | 0.22 | sea-level cosmic 0.27 ionising + 0.05 neutron = 0.32 mSv/yr (UNSCEAR 2008 Table 5); UNSCEAR 2000: 31 nSv/h = 270 uSv/yr ionising plus 48 uSv/yr neutrons | UNSCEAR 2000, 2008 | READ |

Reading: with proton primaries the CORSIKA curve reproduces the aviation point (46 against 47.9) and sea level (0.34 against 0.32). The code's curve with its 0.6 field factor is 2.6 times low at X = 100, 2.7 at X = 200, 1.9 at X = 265 and 1.5 at X = 500 and at sea level. Mars is 255 mSv/yr measured against 125 to 208 in the code. The muon share at sea level is 0.13 to 0.18 of 0.34 mSv/yr (READ). Atri states that composition changes the result by only a few percent for C/N/O atmospheres and flags pure H2 as unstudied (READ). Spread of the Pfotzer maximum: 122.7 +/- 6.3 g/cm2 [R3919], about 100 g/cm2 [R1729], 15 to 26 km (an altitude). The Mars absorbed dose peaks near 6 g/cm2 of atmosphere (protons near 12 g/cm2) and falls toward the surface, with a pressure anticorrelation of 0.28 uGy/day per Pa; the Gale column is 18.9 to 25.7 g/cm2 ([R1729], READ). The solar-cycle spread of the Mars dose equivalent is a factor 3.3. ISS: 731 uSv/day total, 523 from GCR ([R4888], READ).

Magnetic-moment factor, ratio to M = 0 at the same column. The values (READ); ratios COMPUTED from them. M = moment in Earth units.

| M | 0 | 0.1 | 0.15 | 0.25 | 0.5 | 0.75 | 1 | 2 | 3 | 6 | 10 | Status |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| B1036 (mSv/yr) | 0.65 | 0.48 | 0.46 | 0.44 | 0.42 | 0.39 | 0.34 | 0.28 | 0.23 | 0.19 | 0.1 | READ |
| B100 (mSv/yr) | 553 | 527 | 510 | 405 | 257 | 216 | 172 | 53 | 15 | 5.7 | 2.3 | READ |
| Ratio at X = 1036 | 1 | 0.74 | 0.71 | 0.68 | 0.65 | 0.60 | 0.52 | 0.43 | 0.35 | 0.29 | 0.15 | COMPUTED (READ) |
| Ratio at X = 100 | 1 | 0.95 | 0.92 | 0.73 | 0.47 | 0.39 | 0.31 | 0.096 | 0.027 | 0.010 | 0.004 | COMPUTED (READ) |

The READ ratios themselves come from numerical trajectory tracing ([R0187], 1307.4704): 28 million protons, 14 energies from 64 MeV to 524 GeV, four dipole strengths, leapfrog integration through a magnetospheric model including the tail, counting arrivals at a shell 100 km above the surface (READ). They are not Stormer cutoffs.

Independent check (COMPUTED with formula F1, the vertical Stormer cutoff Rc = 14.34 GV (M/M_E)(R_E/r)^2 cos^4(latitude), averaged over the sphere with area weight, and a force-field-modulated proton spectrum whose local interstellar shape is RECALLED, Usoskin-type). Fraction of the unshielded primary protons above the cutoff, for modulation potential phi = 0.3 / 0.6 / 1.0 GV:

| M | 0.15 | 0.5 | 1 | 2 | 10 | Status |
|---|---|---|---|---|---|---|
| Number flux fraction | 0.67 / 0.76 / 0.83 | 0.34 / 0.40 / 0.48 | 0.22 / 0.27 / 0.32 | 0.15 / 0.18 / 0.21 | 0.063 / 0.074 / 0.087 | COMPUTED (RECALLED spectrum) |
| Energy-weighted fraction | 0.84 / 0.91 / 0.95 | 0.57 / 0.67 / 0.75 | 0.43 / 0.52 / 0.60 | 0.31 / 0.38 / 0.45 | 0.14 / 0.175 / 0.21 | COMPUTED (RECALLED spectrum) |
| Griessmeier ratio at X = 100 | 0.92 | 0.47 | 0.31 | 0.096 | 0.004 | COMPUTED (READ) |
| Griessmeier ratio at X = 1036 | 0.71 | 0.65 | 0.52 | 0.43 | 0.15 | COMPUTED (READ) |

Pattern: the X = 1036 ratio follows the energy-weighted fraction (0.52 against 0.52 for M = 1, 0.65 against 0.67 for M = 0.5, 0.15 against 0.175 for M = 10), as expected for dose from high-energy primaries. The X = 100 ratio follows the number-flux fraction within about 0.1 for M up to 1 (0.31 against 0.27, 0.47 against 0.40) but the model is too transparent for M = 2 (0.18 against 0.096) and for M = 0.15 (0.76 against 0.92). Spread to report: within 0.1 to 0.2 for M up to 1, within a factor 2 above. Use the READ ratios (log-interpolated in X between 100 and 1036 and in M); use the model only beyond M = 10. Not valid for non-dipolar or compressed fields, or for X below 100 (no READ ratio exists there).

### 2.8 Stellar energetic particles and flares

| Quantity | Value | Source | Status |
|---|---|---|---|
| Belov relation, flux of protons above 10 MeV | Ip = (4.8 +/- 1.3)e7 Ix^(1.14 +/- 0.14) pfu, Ix = peak 1 to 8 Angstrom X-ray flux in W/m2 | [R4000] (1006.0022) Eq. 1 | 2ND-HAND (read in Segura and Tilley) |
| Check: GOES X13 flare (Ix = 1.3e-3 W/m2) | 2.5e4 pfu; Tilley quotes 24,611 pfu | my evaluation; [R4365] (1711.08484) | COMPUTED (READ) |
| Oct 1989 event, observed peak | about 4e4 pfu above 10 MeV, duration about 27 h | [R4365] | READ |
| Carrington event fluence above 30 MeV | about 1e10 cm-2 ([R3789], Segura); 1.1e10 (Tilley); about 1e10 for 1e32 erg [R0183] | [R0183, R4000, R4365] 1606.07027 | READ |
| Carrington total flare energy | 1e31.9 erg (Tilley, calibrated to fluence); 1e32 (Atri); 5e32 (cited by Tilley) | as listed | READ (spread factor 6) |
| Fluence law (Atri) | F30 = about 1e9 cm-2 for 1e31 erg and 1e10 for 1e32 erg at 1 AU (two points, a linear law in flare energy, COMPUTED); AD 775 class about 1e33 erg | [R0183, R0184] | READ |
| AD Leo great flare, scaled | X-ray 9 W/m2 at 0.16 AU gives 5.9e8 pfu and fluence 1.5e12 cm-2 (200 times Carrington, upper limit) | [R4000] | READ |
| Tilley scaling | 1e30.5 erg flare: 3.25e5 pfu, fluence 1.95e8 cm-2; 1e34 erg: 6.2e8 pfu, 9.6e12 cm-2 | [R4365] | READ |
| Flare frequency (GJ 1243, Kepler, U band) | log10 nu(>E) = -1.01 log10 E + 31.65 (cumulative per day, E in erg); amplitude log10(dF/F) = 0.7607 log10 E - 25.855; duration log10 t = 0.395 log10 E - 9.269 s | [R1851] via [R4365] | 2ND-HAND |
| One 1e34 erg flare per | about 489 days (that star) | [R4365] | READ |
| Proton event given flare | 100% for GOES X2 or stronger (Tilley) | [R4365] | 2ND-HAND |
| Probability a CME hits the planet | P = (Delta + delta) sin[(Delta + delta)/2] / (2 pi sin Theta); 0.083 (restrictive) or 0.25 (permissive) for Delta = 5 pi/183 or pi/2, Theta = pi/4, delta = 1e-8 | [R4365] | 2ND-HAND |
| NOx yield per ion pair | 1.25; 1.3 to 1.6 [R3852]; HOx 2 per pair | [R4288] 1511.04920; [R0190, R4365] | 2ND-HAND |
| Ion pair production rate in the middle stratosphere | 3 to 30 cm-3 s-1 [R4525] | [R4288] | 2ND-HAND |
| Extinction-level SEP doses | 5 to 10 Sv at 700 g/cm2 or less for a hard 1956-type event at Proxima b distance; dose range for M = 1: 2.57e3 to 1.34e-8 Sv across columns 30 to 1000 g/cm2; M = 0.05: 1.46e4 to 5.77e-8 Sv | [R0183] | READ |
| Weak composition dependence of SEP dose | C/N/O only | [R0183] | READ |
| Event spectra | Band function for Feb 1956, Aug 1972, Oct 1989; Weibull [R4789] for Oct 1989; neutron-monitor fit for Sep 1989. Parameters are in figures only | [R0183, R1728] 1705.06763 | READ (parameters NOT REACHED) |
| Share of Mars surface dose from primary protons above 1 GeV | under 4% | [R1727, R1728] | READ |
| Minimum proton energy for a surface-detectable cascade, Earth | about 450 MeV | Shea and Smart | 2ND-HAND |

SEP survival thresholds in the atmosphere come from the PSTAR inversion table in section 2.5: only protons above 171 MeV (22.5 g/cm2), 420 MeV (100 g/cm2) or 2.2 GeV (1033 g/cm2) reach those depths as primaries in an air-like column; 285, 730 and 4660 MeV in H2 for the same columns (COMPUTED). The secondary cascade extends the dose slightly below these thresholds ([R1727] finds 140 to 160 MeV at 22.5 g/cm2).

### 2.9 Natural radionuclides: decay data, heating and abundances

| Nuclide | Half-life (yr) | Q or heat-retaining Q_h (MeV) | Heating per kg element (uW/kg) | Natural abundance and notes | Source | Status |
|---|---|---|---|---|---|---|
| U-238 | 4.4683e9 | Q 51.694; Q_h 47.6 | U element 98.293 (with U-235) | 99.2742% of U; 8 alphas to Pb-206 | [R2987] (1912.04655) | READ |
| U-235 | 7.0348e8 | Q 46.397; Q_h 44.4 | in the 98.293 | 0.7204% of U; 7 alphas | [R2987] | READ |
| Th-232 | 1.41e10 | Q 42.646; Q_h 40.4 | 26.180 | 100% of Th; 6 alphas | [R2987] | READ |
| K-40 | 1.262e9 (geochronology value); NUBASE2020 quotes 1.248e9 | Q 1.331; Q_h 0.676 | K element 0.003387 | 0.01167% of K (117 ppm); 89.44% beta-minus to Ca-40, 10.56% electron capture to Ar-40 (McDonough); 89.6% beta in Wikipedia citing NUBASE | [R2987]; Wikipedia Potassium-40 (NUBASE2020) | READ (spread 1.1% in half-life, 0.2 percentage point in branching) |
| Rb-87 | 4.961e10 | beta only | Rb element 0.01136 | 27.83% of Rb | [R2987] | READ |
| Al-26 | 7.17e5 (0.717 Myr) | 1809 keV gamma | overall 1.777e4 nW/kg-element = 17.8 uW per kg Al | initial 26Al/27Al = 5.25e-5 | [R2987]; Wikipedia Aluminium-26 for the half-life | READ |
| Fe-60 | 2.62e6 | | 0.3635 nW/kg-element | initial 60Fe/56Fe = 1.01e-8 | [R2987] | READ |

Heating formula [R2987]: h [uW/m3] = rho [kg/m3] (3.387e-3 K + 0.01139 Rb + 0.04595 Sm + 26.18 Th + 98.29 U), with abundances in kg of element per kg of rock (so 1 ppm U gives 98.3 pW/kg of rock). First-principles check (COMPUTED from activities times Q_h): U 98.2 uW/kg element (94.1 from U-238 plus 4.05 from U-235), Th 26.2, K 0.00336, which reproduces the READ values to under 1%.

Calibrated bulk compositions (calibration points, not presets). U and Th in ppm, K in ppm or wt%.

| Reservoir | U | Th | K | Th/U | K/U | Source | Status |
|---|---|---|---|---|---|---|---|
| Bulk silicate Earth (BSE), McDonough and Sun / Arevalo | 0.020 | 0.0754 to 0.077 | 280 ppm | 3.77 (mass) | 1.4e4 | [R2987] | READ |
| Upper continental crust | 2.7 | 10.5 | 2.32 wt% | 4.0 | 8900 | [R4726] (1907.12184) | READ |
| Middle continental crust | 0.82 to 0.94 | 3.6 to 4.6 | 1.4 to 1.7 wt% | about 4.4 to 4.9 | | [R4726] | READ |
| Lower continental crust | 0.15 to 0.19 | 0.81 to 1.14 | 0.7 to 0.9 wt% | about 5.4 to 6.0 | | [R4726] | READ |
| Bulk continental crust (Wipperfurth) | 1.3 | 5.8 | 1.8 wt% | 4.5 | 11800 | [R4726] | READ |
| Bulk continental crust (Rudnick and Gao) | 1.3 | 5.6 | 15000 ppm | 4.3 | 11500 | [R2988] (2209.13746) | READ |
| Bulk continental crust (Hacker range) | 1.09 to 1.33 | 4.2 to 5.3 | 11900 to 18800 ppm | | | [R2988] | READ |
| Oceanic crust (MORB) | | | | | | [R4687] | NOT REACHED |

Reservoir power and stellar-composition anchors:

| Item | Value | Source | Status |
|---|---|---|---|
| BSE radiogenic power | 19.9 TW (U-238 38.4%, U-235 1.6%, Th 40.2%, K-40 19.3%) | [R2987] | READ |
| BSE model spread | 10 TW (low), 20 +/- 5 TW (medium), 30 to 38 TW (high) | [R2988] | READ |
| BSE power recomputed | 19.95 TW from 20 ppb U, 77 ppb Th, 280 ppm K, 0.6 ppm Rb and mass 4.042e24 kg | my sum of the heating law | COMPUTED |
| Chondrite refractory enrichment | 1 to 2.2 times CI | [R2988] | READ |
| Solar photosphere / meteorite log eps | K 5.03 / 5.08; Rb 2.52 / 2.36; Th 0.02 / 0.06; U none / -0.54; Eu 0.52 / 0.51; Nd 1.42 / 1.45 | [R0169] (0909.0948) Table 1 | READ |
| CI chondrite derived from the meteorite eps with Si 10.7 wt% | K about 550 ppm, Th about 31 ppb, U about 8 ppb | my conversion | COMPUTED (RECALLED Si) |

Stellar abundance scaling (Th and U follow r-process Eu, not s-process Nd):

| Quantity | Value | Source | Status |
|---|---|---|---|
| log eps_Th in solar twins and analogues, range | 59% to 251% of solar | [R4469] (1505.00280) | READ |
| Thin-disc solar twins: [Th/H] at ZAMS | mean +0.080 dex, scatter 0.058 dex, range -0.117 to +0.257 (53 stars) | [R0479] (1810.10413) | READ |
| Th follows Eu, not Nd (Nd about 60% s-process) | qualitative | [R0479] | READ |
| Th and U 50% condensation temperatures | 1604 K and 1647 K | [R2726] via [R4469] | 2ND-HAND |
| Model Th for example planets | 198.3, 79 and 46.6 ppb; assumed Th/U = 4, K/U = 1e4 | [R4469] | READ |
| Back-extrapolated heating at t = 0 (BSE-like) | 51.6 TW | [R4469] | READ |
| Th/U growth with cosmic age; K scaling with [Fe/H] | | Frebel / Hill cosmochronometers; K trend | NOT REACHED |

Heating-versus-age law (COMPUTED from the READ shares above, heating of a fixed rock mass, relative to its value at 4.5 Gyr, radiogenic K, U, Th, Rb only): age 0.01, 0.1, 0.5, 1, 2, 3, 4.5, 6, 8, 10, 13 Gyr gives 4.89, 4.65, 3.77, 2.97, 1.98, 1.43, 1.00, 0.771, 0.595, 0.484, 0.374. Al-26 adds 4.2e-7 W/kg at formation for BSE rock with Al at 2.35 wt% (RECALLED), 8.5e4 times the 4.9e-12 W/kg of the long-lived set at 4.5 Gyr, halving every 0.717 Myr (COMPUTED). McDonough reports early heating of 1e3 to 1e4 times the later value for the planetesimal-forming epochs (READ).

### 2.10 Ground gamma dose, radon and noble-gas outgassing

Specific activities from half-lives and atomic masses (COMPUTED): U-238 12.35 Bq/mg U (natural U), U-235 0.569 Bq/mg U, Th-232 4.043 Bq/mg, K-40 0.0310 Bq/mg K (0.0317 with the NUBASE half-life), Rb-87 0.868 Bq/mg Rb. Therefore 1 ppm U = 12.35 Bq/kg (U-238 chain), 1 ppm Th = 4.04 Bq/kg, 1 wt% K = 314 Bq/kg, 1 ppm Rb = 0.87 Bq/kg.

| Item | Value | Source | Status |
|---|---|---|---|
| Gamma dose rate in air 1 m above ground per Bq/kg | K-40 0.0417, U series 0.462, Th series 0.604 nGy/h | UNSCEAR 2000 Annex B Table 6; UNSCEAR 2008 | READ |
| Absorbed to effective dose | 0.7 Sv/Gy; occupancy 0.8 indoors, 0.2 outdoors | UNSCEAR 2000, 2008 | READ |
| Per unit abundance | 5.70 nGy/h per ppm U, 2.44 per ppm Th, 13.1 per wt% K | the two rows above times specific activity | COMPUTED |
| Bulk continental crust example (1.3 U, 5.6 Th, 1.5 wt% K) | 40.7 nGy/h; 0.25 mSv/yr effective if outdoors full time | COMPUTED | COMPUTED |
| Soil activity and dose | 2000 medians K 400, U 35, Th 30 Bq/kg give 51 nGy/h; population-weighted 420, 33, 45 give 60; 2008: K 400 +/- 24, U-238 37 +/- 4, Th-232 33 +/- 3 Bq/kg, outdoor mean 54 nGy/h, indoor/outdoor 1.4 (indoor 84, outdoor 59 in 2000) | UNSCEAR 2000, 2008 | READ |
| Annual effective dose, world mean | cosmic 0.39; terrestrial external 0.48 (0.07 outdoors, 0.41 indoors); inhalation 1.26 (Rn-222 1.15, Rn-220 0.1, U/Th series 0.006); ingestion 0.29 (K-40 0.17, U/Th 0.12); total 2.4 mSv/yr | UNSCEAR 2008 Table 12 | READ |
| Principal gamma lines | K-40 1460.82 keV (10.66% per decay); Bi-214 609.3 keV (reference 100), 1764.5 (33.6), 1120.3 (32.8), 297.8 (24.0), 579.1 (16.0), 1238.1 (12.8); Ac-228 911.2 (25.8), 969.0 (15.8), 338.3 (11.3); Tl-208 2614 keV (99.754% per Tl-208 decay, from an ENSDF comment); Ra-226 186.2 keV | NNDC ENSDF (nndc.bnl.gov) | READ |
| Th-232 chain branching to Tl-208 | 35.94% | | RECALLED |
| Rn-222 half-life; emanation coefficient range | 3.824 d; 0.05 to 0.7 | UNSCEAR 2008 | READ |
| Pore concentration | C_Rn = C_Ra f rho_s C^-1 (1 - C)(m[K_T - 1] + 1)^-1; example 78 kBq/m3 (warm, moist) against 21 kBq/m3 (cold, dry) | UNSCEAR 2008 Eq. 1 | READ |
| Radon flux and concentration | worldwide mean flux 0.016 Bq m-2 s-1; outdoor 10 Bq/m3 (range 1 to over 100); indoor 40 Bq/m3; equilibrium factor 0.6 outdoors, 0.4 indoors; 9 nSv per Bq h m-3 of equilibrium-equivalent concentration; doses 1.0 mSv/yr indoors and 0.095 outdoors | UNSCEAR 2008 | READ |
| Radon model checks | J = C_Ra f rho_b sqrt(lambda D_e) with 40 Bq/kg, f = 0.2, D_e = 2e-6 m2/s gives 0.026 Bq m-2 s-1 (0.016 observed); box model C = J tau / h with tau = 5.5 d, h = 1 km gives 7.6 Bq/m3 (10 observed); full-time outdoor dose 9e-9 * 6 * 8766 Sv = 0.47 mSv/yr | my evaluation | COMPUTED |
| Ar-40 from K-40 (BSE, 280 ppm K, mass 4.042e24 kg) | 2.04e42 K-40 atoms at 4.5 Gyr; 2.34e42 Ar-40 atoms produced = 1.55e17 kg over 4.5 Gyr (EC branch 0.1056) | my evaluation | COMPUTED |
| Atmospheric Ar | 0.934% by volume, 1.288% by mass; atmosphere 5.15e18 kg; so 6.6e16 kg Ar, which is 43% of the Ar-40 produced | Wikipedia Argon, Atmosphere of Earth (Lide handbook) | 2ND-HAND (fraction COMPUTED) |
| He-4 from the chains (BSE, U 20 ppb, Th 75.4 ppb) | 3.7e42 atoms over 4.5 Gyr (6.1e18 mol, 2.4e16 kg); production rate 1.6e25 atoms/s = 8.2e8 mol/yr; per kg rock with 1 ppm U and 4 ppm Th, 1.1e-5 mol per Gyr | my evaluation (8, 7, 6 alphas) | COMPUTED |
| Radiogenic dose rate in rock | BSE: 4.9e-12 Gy/s = 0.16 mGy/yr; bulk continental crust: 3.25e-10 W/kg = 10 mGy/yr (upper bound, gamma partly escapes) | COMPUTED from the heating law | COMPUTED |

### 2.11 Radiolysis, neutrons, muons at depth, Europa and trapped belts

| Item | Value | Source | Status |
|---|---|---|---|
| Subsurface GCR energy deposition, peak depth | Mars about 0.6 m, Europa about 1 m, Enceladus about 2 m: 1e6 to 1e7 eV g-1 s-1; deep mine 1e5 to 1e6 | [R0188] (2207.14675); Atri | READ (mine value 2ND-HAND) |
| Same in absorbed dose | 1 Gy/yr = 1.98e8 eV g-1 s-1, so 1e6 to 1e7 is 5 to 51 mGy/yr and the mine range is 0.5 to 5 mGy/yr | conversion | COMPUTED |
| Peak depth in solid cover | 0.6 m of Mars regolith, 1 m of ice, 2 m of porous ice all correspond to roughly 90 to 120 g/cm2 if the densities are about 1.5 to 2, 0.92 and 0.5 g/cm3 (densities are assumptions) | my conversion of the READ depths | COMPUTED (assumed densities) |
| Comet-ice dose peak depth | 5 m for 300 GeV protons in a low-density nucleus | [R1693] (2012.05772) | READ |
| Shielding material ranking for GCR | hydrogen-rich (water, polyethylene) best; regolith intermediate; aluminium helpful only combined with low-Z; total dose 0.213 mGy/day on Mars surface, 0.240 on ISS | Gakis and Atri (2205.13786) | READ |
| Europa surface charged-particle dose | Figure 1 axes span about 1e-1 to 1e7 Gy/yr over depths 1e-6 to 1 m; proton dose dominates at all depths in the amorphisation study; oxygen and sulfur ion doses disagree by over three orders between studies; helium negligible | [R3038] (2505.09490), using [R0879] | READ (curve values not extracted) |
| Europa electrons | below 25 MeV hit the trailing hemisphere preferentially, above 25 MeV the leading one; trailing electron energy confined to the upper 0.1 cm or less, leading deposition at tens of cm | [R3038] (citing Cooper, Paranicas, Nordheim); [R4827] (2604.01927) | READ |
| Europa "5.4 Sv/day" | traced through Wikipedia to a course note; not a primary source; as a conversion 5.4 Sv/day = 1.97e3 Sv/yr, a depth-dependent number, not a surface constant | Wikipedia Europa; [R0879, R3308] not reached | NOT REACHED |
| Electron and bremsstrahlung dominance below about 1 m on trailing side, protons on leading side | search summary only | web search | SNIPPET |
| Radiolytic G-value definition | number of molecules formed per 100 eV absorbed; products H2, O2, H2O2 | [R2177] (1804.10589) | READ |
| G(H2) liquid water, low-LET | 0.45 molecules/100 eV, and 2.4 with a hydrogen-atom donor (formic acid) | Wikipedia Radiolysis citing [R3931] (OSTI 7347831) | 2ND-HAND |
| SI form of G | 1 molecule/100 eV = 0.1036 umol/J, so G(H2) = 0.45 gives 0.0466 umol/J | conversion | COMPUTED |
| H2 per kg of ice at 5 mGy/yr over 4.5 Gyr | 2.3e-10 mol/kg/yr, about 1.0 mol/kg cumulative | COMPUTED | COMPUTED |
| O2 yield from irradiated ice | steady-state yield roughly exponential in temperature, apparent activation energy 0.06 eV; precursor trapping builds up with fluence; the net O2 yield per ion depends on energy density and temperature; a numerical G(O2) was not given | [R2177, R4325] | READ (numbers NOT REACHED) |
| Neutron share at sea level | 0.05 of 0.32 mSv/yr (about 16%); dose from neutrons grows with altitude (0.12 of 0.46) | UNSCEAR 2008 Table 5 | READ |
| Muon vertical intensity, deep rock | I(h) = exp(A1 + A2 h) + exp(A3 + A4 h) + I_nu, A1 = -11.22, A2 = -0.00262, A3 = -14.10, A4 = -0.001213, h in hg/cm2 (1 hg/cm2 = 1 m of water equivalent), I_nu = 2.17e-13 cm-2 s-1 sr-1; valid for at least about 1 km w.e. | [R0933] in [R0581] (hep-ph/9803488) | 2ND-HAND (fit READ in Bugaev) |
| Muon flux and dose below 1 km w.e. | 1 km w.e.: I = 1.2e-6 cm-2 s-1 sr-1, dose about 4e-5 mSv/yr; 3 km w.e.: 2.5e-8 and 8e-7 mSv/yr; 10 km w.e.: 4.3e-12 and 1.4e-10 mSv/yr (omnidirectional flux taken as pi I, mean loss 2 MeV cm2/g); below 1 km w.e. the formula is out of range (100 m w.e. would give 3.5e-4 mSv/yr, an extrapolation) | COMPUTED from the fit | COMPUTED |
| Muon energy loss and critical energy | a_mu about 2.0 MeV cm2/g; b_mu 3.5e-6 cm2/g; critical energies in table 2.6 | [R0581]; PDG | 2ND-HAND / DATA |
| Trapped belts, Earth | inner belt protons from 100 keV to over 400 MeV, outer belt electrons 0.1 to 10 MeV, inner edge near 200 km in the South Atlantic Anomaly; relativistic electrons above 5 MeV excluded inside a sharp inner edge of the outer belt | Wikipedia Van Allen radiation belt (ECSS standard) | 2ND-HAND |
| Belt fluxes, Jovian and other planets | AP8/AE8, AE9/AP9, [R2954] | NOT REACHED | NOT REACHED |

### 2.12 Energy per ion pair W and the ionisation profile

DATA: xraydb 4.5.8 table ionization_potentials (from Knoll, Radiation Detection and Measurement, Table 5-1; fast electrons). Ion pairs per kg per Gy = 1/(W e), COMPUTED.

| Gas | W (eV per ion pair) | Ion pairs per kg per Gy | Source | Status |
|---|---|---|---|---|
| H2 | 36.5 | 1.71e17 | xraydb (Knoll) | DATA |
| He | 41.3 | 1.51e17 | xraydb (Knoll) | DATA |
| N2 | 34.8 | 1.79e17 | xraydb (Knoll) | DATA |
| O2 | 30.8 | 2.03e17 | xraydb (Knoll) | DATA |
| CO2 | 33.0 | 1.89e17 | xraydb (Knoll) | DATA |
| Ar | 26.4 | 2.36e17 | xraydb (Knoll) | DATA |
| CH4 | 27.3 | 2.29e17 | xraydb (Knoll) | DATA |
| Air | 33.8 | 1.85e17 | xraydb (Knoll) | DATA |
| Ne, Kr, Xe | 35.4, 24.4, 22.1 | | xraydb (Knoll) | DATA |
| Air, project value used by [R0190] | 35 | | [R0190] (0804.3207) | READ |
| NH3, H2O vapour, SO2, CO, HCN, CH3Cl | not available | | | NOT REACHED |
| Alpha and proton W in each gas (about 3 to 10% above electron W) | | ICRU 31 | | NOT REACHED |

Reading: W varies by a factor 1.6 across the gases; the ionisation rate per unit absorbed dose is within 10% of air for H2, N2, O2 and CO2, and differs by +24% (CH4), +28% (Ar) and -18% (He). Applying the Atri M = 1 value at X = 100 (171 mSv/yr, treated as an upper bound of 171 mGy/yr), air density 0.145 kg/m3 at that column gives q = D rho / (W e) of about 145 ion pairs cm-3 s-1 (COMPUTED), against 3 to 30 cm-3 s-1 for the middle stratosphere ([R4525], 2ND-HAND). The 5 to 50 times gap shows that the dose equivalent in mSv, which includes neutron and quality weighting, must be reduced to absorbed dose and that the stratospheric range refers to higher altitude (lower density); the factor needs a calibration with an absorbed-dose profile (Boss list item 7).

### 2.13 Biological dose limits and tolerances

| Quantity | Value | Source | Status |
|---|---|---|---|
| Occupational total effective dose, annual limit (NRC) | 0.05 Sv (5 rem); organ 0.5 Sv; lens of the eye 0.15 Sv; skin 0.5 Sv | 10 CFR 20.1201 via law.cornell.edu | READ |
| Public annual limit (NRC) | 1 mSv (0.1 rem) | 10 CFR 20.1301 | READ |
| ICRP occupational | 20 mSv/yr averaged over 5 years, no single year above 50 mSv (some sources add 100 mSv per 5 years); public 1 mSv/yr | Wikipedia Radiation protection (cites ICRPedia), Wikipedia Sievert (cites ICRP 103) | 2ND-HAND |
| ICRP 103 retained these limits; detriment about 5% per Sv; dose and dose-rate effectiveness factor 2 for solid cancers | qualitative text and the 5% value | ICRP 103 free extract (icrp.org) | READ (numeric limits not in the extract) |
| Acute exposure above 8 Gy | almost always lethal even with care | Wikipedia Acute radiation syndrome | 2ND-HAND |
| Lethal doses by taxon | humans 5 to 10 Sv; small mammals 6 to 10; goldfish 100; bats 150; insects 100 to 1000; viruses 1e4; Deinococcus radiodurans 1e4 to 1e5 Sv | [R0183] (1606.07027) | READ |
| LD50/30 for humans without care | 4 to 5 Gy | | RECALLED |
| Project text: Deinococcus above 10,000 to 15,000 Gy; tardigrade (Ramazzottius) 4000 to 5000 Gy | | project documents Naturally Occuring Ionizing Radiation.md, Environmental.Extremes.and.Protection.Research.md | 2ND-HAND |
| World mean natural background | 2.4 mSv/yr | UNSCEAR 2008 | READ |
| Career limit for astronauts, deep-space dose rates, D10 of organisms | | | RECALLED |

## 3. Formulas and parametrisations with validity ranges

F1. Vertical dipole cutoff rigidity of a body (GV): Rc(lat) = 14.34 (M / 7.768e22 A m2) (R_E / r)^2 cos^4(lat), with r the body radius for the surface value. The coefficient is c mu0 M / (4 pi r^2) / 4, evaluated for Earth (COMPUTED); the form is textbook (RECALLED). The scaling is M / r^2. Validity: a dipole field with a standoff distance of several body radii; for a compressed magnetosphere the effective cutoff is lower. Replaces the code's (R_E/R_p)^3.

F2. Bethe stopping power (charge z, mass M): -dE/dx = K z^2 (Z/A) (1/beta^2) [0.5 ln(2 m_e c^2 beta^2 gamma^2 T_max / I^2) - beta^2], K = 0.307075 MeV cm2/mol, T_max = 2 m_e c^2 beta^2 gamma^2 / (1 + 2 gamma m_e/M + (m_e/M)^2). For a compound, Z/A = sum of w_i Z_i/A_i and ln I = sum of w_i (Z_i/A_i) ln I_i / (Z/A). Validity (COMPUTED against PSTAR): protons in gases within 1% above 30 MeV and 2% at 10 MeV (Ar); in ice, rock and iron within 2% to 1 GeV and 8% at 10 GeV (density effect omitted); below 10 MeV shell, Barkas and Bloch terms are needed (flagged, not implemented). For alphas use z = 2, M = 3727 MeV.

F3. GCR dose under a column X with field M: D = D0(X) f_B(M, X) f_sol, where D0(X) is the no-field curve through the READ points (X = 100: 553, 200: 263, 500: 17.2, 700: 3.02, 1036: 0.65 mSv/yr) interpolated linearly in log-log, f_B is the READ ratio table in section 2.7, and f_sol is the solar-cycle factor (Mars 0.5 to 1.64 of the mean, a total spread of 3.3). D0 is proton-primary only and so is a lower bound by the heavier-primary contribution (NOT REACHED). Below X = 100 use the Moon and Mars anchors and treat the shape as unknown.

F4. SEP event: F30 = 1e9 cm-2 (E_flare / 1e31 erg) at 1 AU, scaled by (1 AU / d)^2 (COMPUTED from two READ points), or the X-ray route Ip = 4.8e7 Ix^1.14 pfu. Event rate from log10 nu = -1.01 log10 E + 31.65 per day (GJ 1243, READ) times P_CME (2.8). Spectrum J(>R) = J0 exp(-R/R0) with R0 spread unknown (the code's 0.2 GV is a soft event; Band and Weibull parameters NOT REACHED). Validity: M dwarf use extrapolates two to four decades in flare energy.

F5. Column threshold for a vertical proton: T_min(X) from inverting the PSTAR range (table in 2.5; for X from 1 to 1033 g/cm2 in air T_min = 30, 55, 108, 171, 274, 420, 868, 2241 MeV at X = 1, 3, 10, 22.5, 50, 100, 300, 1033). Light gases need more energy per g/cm2: multiply by 1.65 for H2, 1.06 for He, 1.07 for water, 0.96 for rock (ratios at 22.5 g/cm2, COMPUTED). Transmitted SEP share = fraction of the SEP spectrum above T_min (primaries only; secondaries add a few tens of percent at the edge, [R1727]).

F6. Photon optical depth: tau(lambda) = sum_i x_i sigma_i(lambda) N, with N the column in molecules/cm2 = X N_A / mu_bar, mu_bar the mean molar mass in g/mol (N = 6.02e23 X / mu_bar). The tau = 1 mass column is mu_bar / (N_A sum x_i sigma_i). Validity: EUV 10 to 100 nm with table 2.2 (line cores of H2, N2, O2 are much stronger than the window means); X-rays to 3 keV with Henke, above with F7.

F7. X-ray and gamma: mu/rho = sum of w_i (mu/rho)_i from the XCOM element tables; dust (iron oxides, perchlorates, carbonates, sulfides, silicates) use their element list; all bound-molecule effects ignored. Energy absorption coefficients mu_en/rho are in the same data file for dose.

F8. Composition shift of cascade depth (design inference, not a measurement): the electromagnetic part of the curve and the Pfotzer maximum scale with X0 and the hadronic part with lambda_I. Take X_EM = X (36.6 / X0_gas) and X_had = X (90.1 / lambda_I_gas) and evaluate the air curve at those columns; for H2, X_EM = 0.58 X, X_had = 1.73 X. For CO2, N2, O2 the factors are within 4% of 1. Validity: C/N/O atmospheres are covered by the READ statement of a few percent. H2, He, CH4 and NH3 columns are unverified; a CORSIKA or Geant4 run is item 6 of the Boss list.

F9. Ground gamma dose in air at 1 m: D [nGy/h] = 5.70 U + 2.44 Th + 13.08 K, U and Th in ppm, K in wt% (COMPUTED from READ coefficients); effective annual dose = D x 8766 h x 0.7 Sv/Gy x occupancy (0.2 outdoors, 0.8 indoors, UNSCEAR). For a dense atmosphere multiply by exp(-(mu/rho) rho_gas h) with h = 100 cm and (mu/rho) from table 2.4 at about 1 MeV (design inference, equals 0.66 for CO2 at 65 kg/m3); valid for rho_gas below about 0.03 g/cm3 where the Earth coefficients apply to first order.

F10. Heating and age: h = rho (3.387e-3 K + 0.01139 Rb + 0.04595 Sm + 26.18 Th + 98.29 U) uW/m3 with abundances as mass fractions in kg/kg; evolve each element by exp(lambda t) using the half-lives in 2.9. The age law in 2.9 gives heating relative to 4.5 Gyr from 4.9 at 0.01 Gyr to 0.37 at 13 Gyr for BSE shares. Add Al-26 and Fe-60 for the first 20 Myr only.

F11. Radon: J = C_Ra f rho_b sqrt(lambda D_e) (f = emanation 0.05 to 0.7, D_e effective diffusion coefficient, lambda = ln 2 / 3.824 d); C = J tau / h_mix with tau = 5.52 d; E [mSv/yr] = 9e-6 F C T_h with F = 0.6 outdoors, 0.4 indoors, T_h occupied hours. Validity: Earth-like soil and h_mix of order 1 km; D_e versus pressure is NOT REACHED, so treat it as an input with at least a factor-2 uncertainty.

F12. Noble-gas outgassing: n(Ar-40) = n(K-40 at 4.5 Gyr) (exp(lambda_K t) - 1) x 0.1056; n(He-4) = 8 n(U-238)(exp(lambda_238 t) - 1) + 7 n(U-235)(exp(lambda_235 t) - 1) + 6 n(Th-232)(exp(lambda_232 t) - 1). The fraction released to the atmosphere is a calibration: Earth 0.43 for Ar-40 (COMPUTED from 2.10). Helium escapes thermally, so the atmospheric He inventory needs part 04.

F13. Ionisation rate profile: q(z) = rho(z) D(z) / (W e), with D in Gy/s from the column curve D(X) (absorbed, not dose equivalent) and W from section 2.12; or q(z) = (dE/dz) / W for a single beam, with dE/dz from PSTAR in 2.5. NOx production is 1.25 per ion pair, HOx 2 per pair. Validity: Chapman-style, no recombination; calibrate (Boss item 7).

F14. Radiolysis: n(H2) = G E_dep / 100 eV, with the SI conversion 0.1036 umol/J per unit G, and E_dep from the dose rates in 2.11 (Gy = J/kg).

F15. Muon dose underground: Crouch fit (2.11), valid from about 1 km w.e. down; shallower, interpolate toward the surface muon dose of 0.13 to 0.18 mSv/yr (sea-level READ).

## 4. Surface dose recipe and code status

Surface dose rate D_s [mSv/yr] = GCR + SEP + gamma + radon (+ UV index, separate quantity). Each part, whether it exists in `src/planetgen/physics/radiation.py` (READ in this session) and what is new.

| Part | Formula | Exists? | What to change | Constants and labels |
|---|---|---|---|---|
| Column | X = P0 / g | exists (`column_g_cm2`) | none | trivial |
| GCR | D0(X) f_B(M, X) f_sol [x f_comp] | exists (`gcr_dose_msv_yr`, `magnetic_factor`, `helio_multiplier`) | replace curve and field factor with the READ tables; add F8 composition scaling; add f_sol with 0.5 to 1.64 spread | READ (Atri, Griessmeier, Guo, Zhang); F8 design inference |
| Cutoff | F1 | exists but wrong exponent (`cutoff_rigidity_gv`) | use M / r^2 | COMPUTED, RECALLED form |
| SEP | event fluence F4 x spectrum above T_min (F5) | exists in simplified form (`sep_dose_msv_yr`, `free_sep_msv_yr`, `proton_cutoff_gv`) | range law is fine to 300 MeV; use the PSTAR table above and a fluence-based free dose | READ (Atri, Segura, Tilley); spectra NOT REACHED |
| Ground gamma | F9 with K, U, Th | exists as a constant 0.48 times abundance times heat ratio (`ground_dose_msv_yr`) | derive from K, U, Th; star Th scaling; dense-gas factor | READ (UNSCEAR), COMPUTED |
| Radon | F11 | exists as constant 0.5 if P > 1 kPa and land (`generate_dose`) | exhalation, mixing height, occupancy | READ (UNSCEAR), D_e NOT REACHED |
| UV index | not in this part | exists (`uv_index`) | see part 02 | |
| Ice, rock or ocean cover | F15 and [R0188] deposition | new | depth-dose for subsurface habitats; GCR peak near 90 to 120 g/cm2 of cover | READ [R0188], assumed densities |
| Subsurface radiolysis | F14, section 2.11 | new | H2 and O2 budgets for chemotrophs | G(H2) 2ND-HAND, G(O2) NOT REACHED |
| Radiogenic noble gases | F12 | new | Ar-40 and He-4 inventories | COMPUTED |
| Ionisation profile | F13 | new | feeds NOx, HOx and ozone-loss chemistry (part 04) | DATA (W), 2ND-HAND (NOx yield) |
| XUV and X-ray column | F6, F7 | new | sets the altitude of energy deposition and the XUV that reaches the lower atmosphere | COMPUTED (DATA) |
| Trapped belts | | new | for moons inside giant-planet magnetospheres | NOT REACHED |

Worked calibration (COMPUTED; Earth as a calibration point only): X = 1033, M = 1, f_sol = 1 gives GCR 0.34 mSv/yr (READ 0.32 measured). Ground: continental crust (1.3, 5.6, 1.5 wt%) gives 40.7 nGy/h, 0.25 mSv/yr outdoors full time, 0.07 with 0.2 occupancy (READ 0.07). Radon: 9e-6 * 0.6 * 10 * 1760 = 0.095 mSv/yr (READ 0.095). Total outdoor with occupancy 0.2 is about 0.34 + 0.07 + 0.095 = 0.51 mSv/yr, the code's total outdoor value on Earth is 0.22 + 0.48 + 0.5 = 1.2 mSv/yr because it carries the indoor-weighted ground term and the full-time radon term.

## 5. Boss must research

1. Band-function parameters, [R4789] Weibull parameters, and a flare-energy to SEP-spectrum map for hard and soft events (Youngblood and Airapetian flare-proton links). Replaces the single R0 = 0.2 GV.
2. [R1843] RAD surface dose (Science; domain refused). The 0.64 mSv/day figure is SNIPPET; [R1730] (READ) is used instead.
3. [R0879, R3308] (Europa dose versus depth, origin of 5.4 Sv/day); tables; belt flux models AP8/AE8, AE9/AP9, [R1118, R2954].
4. [R4687] (oceanic crust U, Th, K); measured K/U and Th/U for the Moon, Mars, Venus, Mercury and chondrite classes (K volatile depletion).
5. Stellar Th/U/Eu and K trends: Frebel and Hill cosmochronometry; K versus [Fe/H] and [alpha/Fe] for K-40 scaling (Unterborn and Botelho cover Th only).
6. CORSIKA or Geant4 runs for H2, He, CH4, NH3 and CO2 columns from 0.01 to 100 bar with helium and heavier primaries, to test F8; the helium and HZE share of the dose.
7. Absorbed-dose (Gy) and ionisation-rate profiles in air to calibrate F13 (Bazilevskaya; [R4525]); W values for alphas and protons and for NH3, H2O vapour, SO2, CO, HCN, CH3Cl.
8. Photoabsorption cross sections 10 to 100 nm for Ar, Ne, Kr, Xe, HCN, H2CO, HCl, HF, CH3Cl, SiO, TiO, VO, Na, K, Mg, Fe (outer-shell data; Henke gives the inner-shell part only).
9. ICRP 103 numeric limits (not in the free extract), ICRP 132 and 137, NCRP 153; primary sources for human LD50/30 and organism D10 values.
10. Radon transport (effective diffusion, advective release) in thin and dense atmospheres and porous regolith; origin of the code's 1 kPa threshold.
11. Secondary and albedo neutron dose at depth in regolith and ice (Gakis and Atri, Gronoff partly mined).
12. PDG helium and argon pages (X0, lambda_I); the unidentified "Ribbon fluence" reference from the task notes; origin of the code's helio-multiplier values.

## 6. Recalled / to verify

| Value | Where used | What to check |
|---|---|---|
| Vertical dipole cutoff form 14.9 cos^4(lat) / r^2 GV and the 59.6 GV coefficient | F1 | a textbook |
| Local interstellar proton spectrum J ~ P^-2.78 / (1 + 0.4866 P^-2.51) | section 2.7 check | no source |
| Bulk silicate Earth Al 2.35 wt% (Al2O3 4.45 wt%) | Al-26 heating | McDonough and Sun 1995 |
| CI chondrite Si 10.7 wt% | CI K, Th, U | Lodders tables |
| Th-232 to Tl-208 branching 35.94% | gamma line yield | NNDC |
| Ar first ionisation energy 15.76 eV (78.7 nm) | Ar cross-section edge | NIST ASD |
| Human LD50/30 of 4 to 5 Gy; astronaut career limits; deep-space cruise dose rates; D10 values | section 2.13 | ICRP, NASA-STD-3001 |
| Concorde and polar-flight dose rates near the Pfotzer maximum | comparison to Atri at X = 100 | EURADOS reports |
| Knoll W values apply to electrons; alpha values are somewhat higher | section 2.12 | ICRU 31 |
| Bazilevskaya ionisation maximum near the Pfotzer layer | section 2.12 | no source |
| Early radiogenic heating "1e3 to 1e4 times" statement | section 2.9 | [R2987] section on short-lived nuclides |

## 7. Cross-links to neighbouring parts

- Part 01 (gas table, thermodynamics): molar masses and Z/A for F2 and F6; gas density for F9 and F13.
- Part 02 (UV and visible optics): stellar spectrum for table 2.2, the 200 to 300 nm UV index, ozone cross sections that overlap the Leiden O3 file.
- Part 03 (infrared, thermal): radiogenic heating F10 feeds interior heat flux.
- Part 04 (photochemistry, ozone, escape): XUV deposition altitude, ion-pair NOx and HOx (F13), SEP ozone loss (Segura, Tilley), escape of H2 and He.
- Part 05 (mixing ratios): mole fractions for tables 2.3 and 2.4, Ar-40 and radon content, noble gases.
- Part 06 (condensates, hazes, dust): regolith and dust composition for F7 and F9; hazes do not shield ionising radiation at their column masses.
- Repository: `docs/design/habitability/activity-magnetism-radiation-hydrosphere.md` sections 2 to 4, 6, 7; `src/planetgen/physics/radiation.py`.

## 8. Domains that worked in this part

Worked: export.arxiv.org (PDF), arxiv.org (PDF, search), physics.nist.gov (PSTAR, ASTAR, XCOM), pdg.lbl.gov, www.unscear.org, nndc.bnl.gov, home.strw.leidenuniv.nl, europepmc.org, api.crossref.org, en.wikipedia.org (raw), pypi.org and files.pythonhosted.org, law.cornell.edu, www.icrp.org (page and free-extract PDF). Reachable, unused: ntrs.nasa.gov landing page, osti.gov, zenodo.org. Refused or unusable: api.github.com, un-ilibrary.org, ADS, inis.iaea.org, www.iaea.org, journals.sagepub.com, nature.com, www.ecfr.gov, a nasa.gov STD link (404).
