# Part 03: Infrared absorption, greenhouse efficiency and thermal structure from composition

Research lane for the atmosphere-composition work (v9.0). Read-only research; nothing in the repo was edited. All Solar System bodies appear only as TEST points and are labelled so. Status labels follow `_worker-rules.md` (READ, 2ND-HAND, DATA, COMPUTED, SNIPPET, RECALLED, NOT REACHED).

Reachability record for this lane (what actually worked): export.arxiv.org and arxiv.org PDFs (abs and pdf, but arxiv.org/search is a bot page and the export API answers 429); hitran.org (CIA files, MT_CKD page); raw.githubusercontent.com (MT_CKD_H2O netCDF, refractiveindex.info database); pypi.org and files.pythonhosted.org (fair, exo_k, petitRADTRANS sdist); webbook.nist.gov (Shomate cp, Shimanouchi vibrational frequencies); ntrs.nasa.gov API (abstracts); nature.com (abstract only); ipcc.ch (AR6 Ch.7 PDF); api.crossref.org; dspace.library.uvic.ca (author copies of Byrne and Goldblatt 2014 and Goldblatt 2013); courses.seas.harvard.edu (Etminan 2016 PDF); faculty.washington.edu/dcatling (Robinson and Catling 2012); users.physics.ox.ac.uk (Wordsworth and Pierrehumbert 2013). Refused or failing: dash.harvard.edu (429), agupubs, science.org, doi.org. No opacity data tables ship inside the exo_k wheel or the petitRADTRANS sdist (they download them), and the PyPI package named "clima" is an unrelated command-line library, not a climate model.

---

## 1. Scope and how the numbers feed the generator

Today (`physics/planets.py` 696-716, and the trace in `reports/composition-density-link-2026-10-17.md` item 4) the surface temperature is `T_s = [(1-A) S (1+G) / (4 sigma)]^(1/4)` with `G = min(cap, base_ratio * class_multiplier)`, where `base_ratio` is the drawn molar density over the CO2 value. The mix, the pressure and the gravity never enter G. This part supplies what is needed to replace G by a function of composition, pressure, gravity and absorbed flux:

| Input (all available after composition is drawn) | Output | Where in this file |
|---|---|---|
| mole fractions of the radiatively active gases, total pressure P_s, gravity g, mean molar mass | equivalent grey infrared optical depth `tau0` at the surface | sections 2, 3 (needs line data or a precomputed grid; CIA part is coded here) |
| absorbed flux F = (1-A) S f / 4-type budget, `tau0`, adiabatic index gamma, lapse scaling alpha, exponent n | surface temperature T_s, radiative-convective boundary, tropopause pressure, downwelling IR | section 5 (closed enough to code) |
| composition, P_s, g, stellar spectrum | Bond albedo from Rayleigh scattering plus surface plus cloud | section 7 |
| stored T_s, P_s, S | consistent (A, tau0) pairs, and the implied "greenhouse factor" G = (T_s/T_e)^4 - 1 for the existing field | section 6 |
| T_eq, P_s, tau0 | day-to-night heat redistribution factor | section 8 |

The existing multiplier is exactly `G + 1 = (T_s/T_e)^4`, so every result below can be reported in the repo's own variable. COMPUTED from the test inputs in section 6: Earth G = 0.625, Venus G = 99.6 (the trace report quotes about 101 for Venus's ratio), Titan G = 0.69.

---

## 2. Greenhouse gases: band centres, strength class, forcing

### 2.1 Band centres and strength (per gas)

Band positions are the Shimanouchi 1972 gas-phase fundamentals as tabulated at NIST (DATA, accessed at webbook.nist.gov, mask 800, public domain US government data). Where NIST returned no polyatomic table (CO2, O3, N2O, SO2, OCS, NO2, HCl, HF, CO) the textbook values are RECALLED and repeated in the to-verify list. Strength class column: `T` = in the Byrne and Goldblatt 2014b (arXiv:1409.1880) HITRAN-2012 saturated blackbody absorption ranking (READ: "H2O ranks 11th strongest, CO2 and CH4 rank 16th and 27th" of the ~40 HITRAN gases tested, so per molecule many minor gases beat CO2 and CH4; their abundance, not strength, is what limits them). `F` = forcing statements from the same paper (READ): "radiative forcings of up to a few to 10 W m-2 from concentrations of 0.1-1 ppmv for many gases".

| Gas | Fundamentals, cm-1 (um) | Strength class and notes | Source | Status |
|---|---|---|---|---|
| H2O | nu2 bend 1594.59 (6.27); nu1 3656.65 (2.73); nu3 3755.79 (2.66); pure rotation band below about 600 cm-1 (beyond 17 um); window 800-1250 cm-1 held by the continuum (section 3.2) | very strong; rank 11; dominant on warm worlds; condensable, so its column is set by T and RH | NIST/Shimanouchi 1972; rank from Byrne and Goldblatt 2014b | DATA; rank READ; rotation edge RECALLED |
| CO2 | nu2 667 (15.0); nu3 2349 (4.26); 2nu2/nu1 Fermi dyad about 1388 and 1285 (7.2, 7.8); hot bands about 960 and 1063 (10.4, 9.4) | strong in band centre, saturated; rank 16; far wings plus CIA matter at tens of bar | recalled | RECALLED; rank READ |
| CH4 | nu4 1306.2 (7.66); nu2 1533.6 (6.52); nu1 2917 (3.43); nu3 3018.9 (3.31) | rank 27 in the saturated metric; strong near-IR absorption also absorbs sunlight, so forcing saturates: peak 10.3, 9, 8.3 W m-2 for 0.5, 1, 2 bar N2 (READ in B&G 2014b) | NIST; B&G 2014b | DATA; READ |
| NH3 | nu2 umbrella 932.5/968.3 (10.3, 10.7); nu4 1626 (6.15); nu1 3336; nu3 3444 | strong; window-filling at 10 um; broadening-sensitive: 55-60 percent of the forcing difference between 0.5 and 2 bar N2 at 5 ppm is pressure broadening (READ) | NIST; B&G 2014b | DATA; READ |
| H2S | nu2 1182.7 (8.45); nu1 2614.6 (3.82); nu3 2626 (3.81) | strong band but narrow; listed among the weaker per-molecule gases in the B&G ordering | NIST | DATA |
| HCN | nu2 bend 711.98 (14.04); nu3 CN 2096.85 (4.77); nu1 CH 3311.47 (3.02) | strong, overlaps the 15 um CO2 band; in the "forcings of order 10 W m-2 at under 1 ppmv" group | NIST; B&G 2014b | DATA; READ |
| C2H6 | CH3 rock nu12 821.6 (12.2); CH3 deform 1379-1469 (6.8-7.3); CH stretch 2896-2985 (3.4); torsion 289 | moderate; fills the 12 um window; important for Titan-like and methane-photochemistry worlds | NIST | DATA |
| C2H2 | nu5 CH bend 730.3 (13.7); nu2 CC 1973.8 (5.07); nu3 3289 (3.04) | strong (rank 6 by the saturated metric); photochemical product | NIST; B&G 2014b | DATA; READ |
| H2CO | nu2 CO 1746.1 (5.73); nu3 CH2 scis 1500.1; nu5 1249.1; nu6 1167.3; nu1/nu4 2782.5/2843.1 | HITRAN lacks the 1000-1550 cm-1 band; B&G note it would be optically thick above 5e-6 (READ) | NIST; B&G 2014b | DATA; READ |
| CH3Cl | nu3 CCl 732.1 (13.7); nu6 rock 1017.3 (9.8); nu2 1354.9 (7.4); nu5 1452.1 | in the strong group of the saturated ranking | NIST | DATA |
| C2H4 | nu7 wag 949.3 (10.5); nu10 rock 826.0 (12.1); nu12 1443.5; nu11 2988.66 | strong at 10.5 um window; HITRAN underestimates (B&G apply a factor 0.1 on abundance) | NIST; B&G 2014b | DATA; READ |
| PH3 | nu2 992.1 (10.1); nu4 1118.3 (8.9); nu1 2322.9; nu3 2327.7 | window absorber; broadening-sensitive | NIST | DATA |
| CH3OH | nu8 CO str 1033 (9.68); nu6 OH bend 1345; CH3 modes 1455-1477, 2844-3000; OH str 3681 | rank 1 in the saturated metric (strongest per molecule) | NIST; B&G 2014b | DATA; READ |
| O3 | nu3 1042 (9.6); nu1 1110; nu2 705 (14.2) | strong, but it is also a strong UV/visible absorber, so it heats the stratosphere (affects tropopause, section 5) | recalled | RECALLED |
| N2O | nu2 589 (17); nu1 1285 (7.8); nu3 2224 (4.5) | moderate; overlaps CO2 | recalled | RECALLED |
| SO2 | nu2 518 (19.3); nu1 1151 (8.7); nu3 1362 (7.3) | strong (rank 8); also a strong UV absorber | recalled; rank B&G 2014b | RECALLED; rank READ |
| CO | fundamental about 2143 (4.67) | poor absorber, still optically thin at 1 bar for a pure-CO atmosphere (READ in Koll 2019, arXiv:1907.13145, his Fig. 4 text) | Koll 2019 | READ; band RECALLED |
| OCS, NO2, HCl, HF, HNO3 | OCS about 520, 859, 2062; NO2 about 750, 1318, 1618; HCl about 2886; HF about 3962 | all listed by B&G among HITRAN gases that reach 10 W m-2 below 1e-5 for the strongest of them | recalled | RECALLED |
| N2, O2, H2, He, Ar, Ne, Kr, Xe | no permanent dipole: no fundamental absorption; only collision-induced bands (section 3.1) and, for O2, weak magnetic-dipole bands | transparent in lines; important only as CIA partners and as broadening gas | physics | READ in Wordsworth and Pierrehumbert 2013 text for H2 and N2 |
| SiO | fundamental about 1242 (8.05) | hot-body species (lava worlds); no reachable line-strength data in this lane | recalled | RECALLED |
| TiO, VO, Na, K, Fe, Mg | electronic bands and resonance lines in the visible and near-IR; they act as shortwave absorbers (large gamma in the Guillot sense) | thermal inversions in hot Jupiters when present (Guillot 2010: strong TiO/VO absorption gives gamma large and a temperature inversion) | Guillot 2010 | READ |

Species in the project text that have no standard band entry here: perchlorates, carbonates, sulfuric acid and mineral dust (they act as aerosols; see section 7 and the aerosol lane), and noble gases Ne, Kr, Xe (transparent; no fundamentals).

Mid-IR forcing regime (READ in Koll 2019 text, citing Pierrehumbert 2010): weak-line limit tau proportional to p; strong-line (Lorentz wing) limit proportional to sqrt(p); band saturation gives logarithmic dependence on pressure and on abundance; CIA gives p squared. This is the physical reason for the linear, square-root and logarithmic forms below.

### 2.2 Forcing parametrisations (coefficients, validity)

C is CO2 in ppm, M is CH4 in ppb, N is N2O in ppb unless stated. F in W m-2, relative to the baseline subscript 0.

| Name | Expression | Validity | Source | Status |
|---|---|---|---|---|
| Myhre 1998 / IPCC 1990 CO2 | F = 5.35 ln(C/C0) | 180-1000 ppm; B&G 2014 find it underestimates above 1000 ppm | Byrne and Goldblatt 2014a (GRL 41, 152) Table 1, citing IPCC 1990 and Myhre et al. 1998 | 2ND-HAND |
| CO2 | F = 4.841 ln(C/C0) + 0.0906 (sqrt C - sqrt C0) | to about 1e4 ppm within 10 percent | same Table 1 | 2ND-HAND |
| WMO 1999 CO2 | F = 3.35 [g(C) - g(C0)], g(C) = ln(1 + 1.2 C + 0.005 C^2 + 1.4e-6 C^3) | very close to line-by-line to 50000 ppm | same Table 1 | 2ND-HAND |
| IPCC 1990 CH4 | F = 0.036 (sqrt M - sqrt M0) - [f(M,N0) - f(M0,N0)] | fails above 20 ppm | same Table 1 | 2ND-HAND |
| IPCC 1990 N2O | F = 0.12 (sqrt N - sqrt N0) - [f(M0,N) - f(M0,N0)] | fails above 10 ppm | same Table 1 | 2ND-HAND |
| overlap f(M,N) | f = 0.47 ln[1 + 2.01e-5 (MN)^0.75 + 5.31e-15 M (MN)^1.52] | all above | same Table 1 | 2ND-HAND |
| Byrne and Goldblatt 2014a CO2 | F = 5.32 ln(C/C0) + 0.39, concentrations as mole fractions, C0 = 278e-6 | 200-10000 ppm | their Table 2 (authoritative version with the erratum brackets) | READ |
| B&G 2014a CH4 | F = 1173 (sqrt M - sqrt M0) - 71636 (sqrt M - sqrt M0)^2 for 0.1-2.5 ppm; F = 0.824 + 5 ln(M/M1)/4 + 5 [ln(M/M1)]^2/1 .. see note | 0.1-100 ppm; M0 = 715e-9, M1 = 2.5e-6 | their Table 2 (the extracted text of the high-range line is garbled: it reads "0.824 + (5/4) ln(M/M1) + (1/5)... [ln(M/M1)]^2"; verify before coding) | READ for 0.1-2.5 ppm line; high-range line to verify |
| B&G 2014a N2O | F = 3899 (sqrt N - sqrt N0) + 38256 (sqrt N - sqrt N0)^2 for 0.1-2.5 ppm; F = 4.182 + 3 ln(N/N1) + 0.5469 [ln(N/N1)]^2 for 2.5-100 ppm | N0 = 270e-9, N1 = 2.5e-6 | their Table 2 | READ |
| B&G 2014a overlaps | N2O-CO2: dF = -16.16 exp(-0.036 (ln(C-C0) - 0.0024)^2 - 0.05 (ln(N-N0) + 6.5)^2); N2O-CH4: dF = -24 exp(-0.02 (ln(M-M0) - 0.01)^2 - 0.044 (ln(N-N0) + 7.73)^2) | overlaps below 0.1 W m-2 up to RCP8.5-2250 levels | their Table 2 | READ |
| Etminan et al. 2016 CO2 | F = [a1 (C-C0)^2 + b1 abs(C-C0) + c1 Nbar + 5.36] ln(C/C0); a1 = -2.4e-7, b1 = 7.2e-4, c1 = -2.1e-4 (units W m-2 ppm-2, ppm-1, ppb-1) | 180-2000 ppm; fit error 3.6 percent (0.15 W m-2) | Etminan 2016 (GRL 43, 12614) Table 1; coefficients also found in `fair` 2.2.5 `forcing/ghg.py` (Apache-2.0) | READ; DATA |
| Etminan 2016 N2O | F = [a2 Cbar + b2 Nbar + c2 Mbar + 0.117] (sqrt N - sqrt N0); a2 = -8.0e-6, b2 = 4.2e-6, c2 = -4.9e-6 | 200-525 ppb | same | READ; DATA |
| Etminan 2016 CH4 | F = [a3 Mbar + b3 Nbar + 0.043] (sqrt M - sqrt M0); a3 = -1.3e-6, b3 = -8.2e-6 | 340-3500 ppb | same | READ; DATA |
| AR6 statement | the Etminan fit is the AR6 standard; Etminan CH4 includes shortwave forcing and a water-continuum update that raised CH4 forcing by 25 percent over Myhre 1998 | n/a | IPCC AR6 WGI Ch.7 text (ipcc.ch PDF) | READ |
| Byrne and Goldblatt 2014b: CO2 | at 1e-2 mole fraction CO2 on 2, 1, 0.5 bar N2: 35, 26, 15 W m-2; shortwave forcing -2, -4, -18 W m-2 at 0.01, 0.1, 1 (Rayleigh scattering of the larger atmosphere plus near-IR absorption) | LBL SMART, HITRAN 2012, Earth-like T-p, N2-H2O background | arXiv:1409.1880 | READ |
| generic gas, many | F up to about 10 W m-2 at 0.1-1 ppmv; the ranking of concentration needed for 10 W m-2 is: CH3OH, HCOOH, HNO3, CH3Cl, COF2, HCN, H2O2, PH3, CH3Br, C2H4, SO2, OCS, NH3, HOCl, O3, N2O ... (panel order of their Fig. 11, extracted from a two-column layout, so the order beyond the first few is uncertain) | 1e-8 to 1e-5, 0.5-2 bar N2 | arXiv:1409.1880 | READ (order partly uncertain) |

Use: the forms are only valid as perturbations around a reference state (Earth-like profile, H2O vapour present, N2 background). They cannot be applied to a CO2-dominated, H2-dominated or steam atmosphere; for those, section 3 and section 5 apply.

### 2.3 Regime-specific results to calibrate against (all from a full radiative-convective or LBL calculation)

| Result | Value | Source | Status |
|---|---|---|---|
| H2-N2 warming of an Archean-like Earth: 2-3 times present N2, H2 mixing ratio 0.1 | surface above 0 C at 75 percent of present flux with CO2 only 2-25 times present; H2-N2 forcing above 24 W m-2 at 280 K or less in the high scenario; at 1 PAL N2 only a few W m-2, 2-3 K | Wordsworth and Pierrehumbert 2013, Science 339, 64 | READ |
| H2-N2 CIA in the window | N2-H2 peaks in the 750-1200 cm-1 window; H2-H2 broad but weak because it scales with the square of H2; N2-N2 far from the blackbody peak at terrestrial T | same | READ |
| Early Mars | 1.3-4 bar CO2 plus water with 5-20 percent H2 raises mean T above freezing | no source | READ (abstract only) |
| Volcanic hydrogen habitable zone | outer edge at 50 percent H2: 2.4 AU for the Sun, Seff 0.1723; at 10 percent: 2.02 AU; inner edge moves only 0.1-4 percent | Ramirez and Kaltenegger 2017, arXiv:1702.08618, Tables 1 and 2 | READ; COMPUTED distances |
| Hot desert world inner edge | as close as 0.38 AU for a solar-like star with about 1 percent relative humidity and raised albedo; 1-10 bar needed to limit water loss | Zsom et al. 2013, arXiv:1304.3714 | READ |
| Water loss with CO2-rich air | efficient loss around G stars only for CO2 0.1-1 bar near the runaway limit; around M stars less than about 270 W m-2 absorbed gives under one ocean lost | Wordsworth and Pierrehumbert 2013 (arXiv:1306.3266) | READ |

---

## 3. Continuum and collision-induced absorption

### 3.1 CIA pairs (HITRAN CIA files, DATA)

Source: hitran.org/cia, files under `/data/CIA/main/`, listed with their year tag. Coefficients are binary absorption coefficients k(nu, T) in cm^5 molecule^-2 (HITRAN convention; sanity check: the H2-N2 grey optical depth computed below for the Wordsworth and Pierrehumbert 2013 scenario comes out at order unity, consistent with their tens of W m-2 forcing). The absorption coefficient in cm-1 is `k * n1 * n2` with n in cm-3. The files carry licence terms of HITRAN (see hitran.org/citepolicy; cite Karman et al. 2019 for the CIA compilation, RECALLED reference). Every row COMPUTED means: I parsed the file, took the block nearest the stated temperature for each separate band block (blocks for the same pair come at different temperatures per band), and formed the Planck-weighted mean `kP(T) = integral k B_nu dnu / integral B_nu dnu` over the full Planck function (zero outside the file's range).

| Pair (file) | T grid in file (K) | Spectral coverage (cm-1) | kP at selected T, cm^5 molecule^-2 | Notes | Status |
|---|---|---|---|---|---|
| N2-N2 (N2-N2_2021) | 70-400 (rotational, 0-450/550/650 cm-1), 200-362 (fundamental 1850-3000 cm-1), 200-330 (overtone 4300-5000) | 0-5000 | 1.7e-45 (100 K), 6.8e-46 (150 K), 3.9e-46 (200 K), 2.6e-46 (250 K), 1.9e-46 (300 K), 1.5e-46 (400 K) | strong at Titan temperature, small at 288 K; the far-IR rotational part dominates | DATA; COMPUTED |
| N2-H2 (N2-H2_2024) | 40-400 (0-1886 cm-1), 296 (3900-5100) | 0.02-5100 | 7.6e-45 (100 K), 9.7e-45 (150 and 200 K), 8.8e-45 (250 K), 7.6e-45 (300 K), 5.9e-45 (400 K) | about 40 times stronger per pair than N2-N2 near 300 K; peaks in the window | DATA; COMPUTED |
| H2-H2 (H2-H2_2011) | 200-3000 | 20-10000 | 3.0e-45 (200 K), 3.1e-45 (250 K), 3.0e-45 (300 K), 2.7e-45 (400 K), 2.4e-45 (600 K), 2.3e-45 (800 K), 2.8e-45 (1000 K), 3.5e-45 (1500 K), 3.9e-45 (2000 K), 2.6e-45 (3000 K) | file does not go below 200 K, so Uranus/Neptune/Saturn upper levels are outside it; alternates at hitran.org/data/CIA/alternate (H2-H2_eq_2018, _norm_2018, ortho-para fractions _fp0.3 to _fp1) exist but were not parsed | DATA; COMPUTED |
| H2-He (H2-He_2011) | 200-9900 | 20-20000 | 1.1e-45 (200 K), 1.1e-45 (250 K), 1.2e-45 (300 K), 1.2e-45 (400 K), 1.1e-45 (600 K), 1.3e-45 (800 K), 1.9e-45 (1000 K), 3.0e-45 (1500 K), 3.6e-45 (2000 K), 3.2e-45 (3000 K) | alternate H2-He_exp_2024 exists; 147 MB file | DATA; COMPUTED |
| CO2-CO2 (CO2-CO2_2024) | 192-800 (several band blocks) | 1-4500 | 4.1e-45 (200 K), 3.0e-45 (250 K), 3.0e-45 (300 K), 3.5e-45 (400 K), 3.4e-45 (600 K), 8e-47 (800 K, only the 1-750 cm-1 block exists there) | covers CIA only; far-wing "continuum" of the 15 um band is not in this file | DATA; COMPUTED |
| CO2-H2 (CO2-H2_2024) | 100-600 | 5-1500 | 6.7e-44 (100 K), 3.8e-44 (200 K), 2.7e-44 (300 K), 2.0e-44 (400 K), 1.2e-44 (600 K) | about 10 times CO2-CO2 per pair; the basis of the early-Mars warming | DATA; COMPUTED |
| CO2-CH4 (CO2-CH4_2024) | 100-600 | 0-720 | 1.2e-43 (100 K), 1.8e-44 (200 K), 7.2e-45 (300 K), 3.5e-45 (400 K), 1.5e-45 (600 K) | steep fall with T | DATA; COMPUTED |
| N2-CH4 (N2-CH4_2024) | 70-400 | 0-800 | 1.4e-44 (100 K), 6.8e-45 (150 K), 4.2e-45 (200 K), 3.0e-45 (250 K), 2.3e-45 (300 K), 1.6e-45 (400 K) | Titan-relevant | DATA; COMPUTED |
| H2-CH4 (H2-CH4_eq_2011, _norm_2011) | 40-400 | 0.02-1946 | 1.5e-44 (111 K), 1.9e-44 (144 K), 1.9e-44 (186 K), 1.7e-44 (240 K), 1.5e-44 (310 K), 1.1e-44 (400 K) | giant-planet and Titan | DATA; COMPUTED |
| O2-O2, O2-N2 (2024 files) | 193-356 | 1150-33670 | O2-O2: 8e-47 (250 K), 1.7e-46 (300 K), 4e-46 (350 K); O2-N2: 8e-47 (250 K), 1.7e-46 (300 K), 3.8e-46 (360 K) | mostly near-IR and visible (about 1.27 and 1.06 um); heats air, irrelevant for the greenhouse | DATA; COMPUTED |
| CO2-H2O (CO2-H2O_2024), N2-H2O (N2-H2O_2018) | 296; 250-350 | 5700-6300; 1930-2830 | negligible in the thermal window (about 1e-55 and 1e-47 weighted) | near-IR pairs only; the CO2-H2O and N2-H2O far-IR continuum is not given here | DATA; COMPUTED |
| CH4-He, CO2-He, N2-He, CO2-Ar, N2-Ar, N2-air, O2-air, H2-H, He-H | see files | see files | not parsed except CH4-He: 1.2e-45 (100 K), 1.0e-45 (150 K), 1.1e-45 (200 K), 8.6e-46 (300 K) | present in HITRAN; Ar, He, H partners are weak | DATA (listing); CH4-He COMPUTED |
| CO2-N2, CO2-O2 far-IR, H2O-H2O, NH3, HCN, C2H6 pairs | none in the HITRAN CIA list | n/a | n/a | the H2O-H2O continuum is separate (3.2); CO2-N2 CIA is missing | NOT REACHED |

Grey conversion for coding (COMPUTED derivation, assumes an isothermal hydrostatic layer, ideal gas, mole fractions x1 and x2 of the pair):

    tau_CIA,thin = kP(T) * x1 * x2 * n0^2 * H / 2
    n0 = P_s / (k_B T)  [cm-3],   H = k_B T / (m_bar g)  [cm]

For optically thick bands the Planck-thin value overestimates. Koll 2019 (his Eq. 13, READ) defines an equivalent grey thickness `tau_LW = -ln[ integral exp(-tau_nu) B_nu dnu / integral B_nu dnu ]`; I applied it to CIA alone, using the file blocks at the stated T:

| Case (isothermal column) | tau_eq (Koll definition), CIA only | Planck-thin | Role | Status |
|---|---|---|---|---|
| 3.3 bar, 90 percent N2 + 10 percent H2, 273 K, g = 9.81 (Wordsworth and Pierrehumbert 2013 scenario, TEST) | 1.26 | 2.6 | explains their large forcing | COMPUTED |
| 1.1 bar, same mix | 0.25 | 0.29 | few W m-2 forcing at 1 PAL N2 | COMPUTED |
| 1 bar N2 alone, 288 K | 0.040 | 0.057 | negligible | COMPUTED |
| Titan-like 1.47 bar N2, 94 K, g = 1.35 (TEST), N2-N2 only | 1.24 | 26 | Titan's N2 CIA is optically thick in its rotational band | COMPUTED |
| Titan-like N2-CH4 pair, 5 percent CH4 (TEST) | 2.5 | 11 | comparable to N2-N2 | COMPUTED |
| 3 bar CO2, 250 K, g = 3.71 (early-Mars-like, TEST) | 0.28 | 14.6 | CIA alone is a modest window closer | COMPUTED |
| 3 bar, 90 percent CO2 + 10 percent H2, 250 K, g = 3.71 (scenario) | 2.4 | 16.8 | an H2 admixture multiplies the CIA effect by about 9 | COMPUTED |
| Mars-like 6.4 mbar CO2, 220 K (TEST) | 8.5e-5 | 8.5e-5 | CIA irrelevant when thin | COMPUTED |
| Venus-like 92 bar CO2, 735 K (TEST, CIA file stops at 800 K and below 750 cm-1) | 0.056 | 65 | in the file's range the CIA is far from the Planck peak; lines and window continuum dominate | COMPUTED |
| 10 bar, 85 percent H2 (rest He), 300 K, g = 12 (hycean-like) | 6.5 (H2-H2) plus 3.6 (H2-He part) | 569 and 39 | a dense H2 envelope is optically thick in CIA alone | COMPUTED |
| 1 bar H2/He, 165 K, g = 24.8 (giant-like; nearest file block used, outside the 25 percent temperature tolerance in spirit) | 1.9 | 4.6 | | COMPUTED |

The tau_eq values are lower bounds on the real grey thickness because they omit lines, the water continuum and the non-isothermal profile, and they add (not exactly) to the line part.

### 3.2 Water vapour continuum (MT_CKD 4.3)

Data: `data/absco-ref_wv-mt-ckd.nc` from github.com/AER-RC/MT_CKD_H2O (read via raw.githubusercontent.com; netCDF classic, 2003 points, 10 cm-1 steps from -20 to 20000; reference 1013 mbar and 296 K; licence file in the repo, AER; general reference Mlawer et al. 2012, doi 10.1098/rsta.2011.0295; version note says changes to the foreign continuum above 4000 cm-1 in 2023). The file stores the self and foreign continuum coefficients C(nu) in cm^2 molecule^-1 cm-1; the scaling I used is `sigma(nu) = nu tanh(h c nu / 2 k T) * [C_self (n_h2o/n_ref) (296/T)^texp + C_for (n_for/n_ref)]` in cm^2 per H2O molecule, which follows the text of the repo README (density scaling, self temperature exponent file `self_texp`) but whose exact radiation-term convention I did not verify in the user guide (to verify).

| nu (cm-1) | wavelength (um) | C_self (296 K) | C_foreign | self T exponent | sigma_self at n_h2o = n_ref (cm2) | sigma_foreign at n_for = n_ref (cm2) | Status |
|---|---|---|---|---|---|---|---|
| 200 | 50 | 3.03e-22 | 2.20e-23 | 2.76 | 2.7e-20 | 2.0e-21 | DATA; COMPUTED |
| 400 | 25 | 2.22e-23 | 2.34e-25 | 1.63 | 6.7e-21 | 7.0e-23 | DATA; COMPUTED |
| 600 | 16.7 | 2.64e-24 | 1.27e-26 | 3.06 | 1.4e-21 | 6.8e-24 | DATA; COMPUTED |
| 800 | 12.5 | 4.62e-25 | 1.33e-27 | 5.06 | 3.5e-22 | 1.0e-24 | DATA; COMPUTED |
| 1000 | 10 | 1.33e-25 | 2.45e-28 | 5.64 | 1.3e-22 | 2.4e-25 | DATA; COMPUTED |
| 1200 | 8.3 | 7.0e-26 | 4.2e-28 | 6.11 | 8.4e-23 | 5.0e-25 | DATA; COMPUTED |
| 2500 | 4.0 | 3.4e-27 | 7.3e-31 | 6.21 | 8.5e-24 | 1.8e-27 | DATA; COMPUTED |
| 3600 | 2.8 | 3.9e-25 | 3.2e-26 | 2.42 | 1.4e-21 | 1.2e-22 | DATA; COMPUTED |

Reading: the self continuum dominates over the foreign continuum by factors of 100-1000 in the window, and its temperature exponent rises to about 5-6 in the window, so the window closes steeply with warming: this is the physical origin of the water runaway (Goldblatt, Robinson, Zahnle and Crisp 2013: at about 280-310 K the surface starts to radiate through the window only weakly; above 310 K the emitting level stays between 250 and 300 K, READ in the arXiv-less author copy, Nature Geosci. 6, 661). Kopparapu et al. 2013 (READ) found that omitting the continuum moves the runaway limit from 0.97 AU to 0.76 AU.

---

## 4. Specific heat, lapse rates, latent heat

Heat capacities from the NIST Shomate fits (DATA; NIST-JANAF; valid ranges differ per species and the CO2, H2, CH4, NH3, SO2, H2S, CO, N2O, HCN fits start at 298 K, N2 and O2 at 100 K, H2O at 500 K). Computed from the Shomate polynomial; cp/R and gamma at 300 K are COMPUTED from them. Molar masses from IUPAC (RECALLED values to 4 digits, trivial).

| Gas | M (g/mol) | cp (J/mol/K) at 300 K | at 500 K | at 1000 K | cp/R at 300 K | gamma (cp/cv) at 300 K | cp (J/kg/K) at 300 K | Status |
|---|---|---|---|---|---|---|---|---|
| N2 | 28.013 | 29.12 | 29.58 | 32.69 | 3.50 | 1.400 | 1040 | DATA; COMPUTED |
| O2 | 31.998 | 29.39 | 31.08 | 34.86 | 3.54 | 1.394 | 919 | DATA; COMPUTED |
| CO2 | 44.010 | 37.22 | 44.61 | 54.30 | 4.48 | 1.288 | 846 | DATA; COMPUTED |
| H2 (equilibrium-type fit, 298 K and up) | 2.016 | 28.85 | 29.26 | 30.20 | 3.47 | 1.405 | 14311 | DATA; COMPUTED |
| H2O | 18.015 | not in range | 35.22 | 41.27 | n/a | about 1.33 steam at 500 K | n/a | DATA; COMPUTED |
| CH4 | 16.043 | 35.71 | 46.35 | 71.79 | 4.30 | 1.303 | 2226 | DATA; COMPUTED |
| NH3 | 17.031 | 35.70 | 42.03 | 56.50 | 4.29 | 1.304 | 2096 | DATA; COMPUTED |
| SO2 | 64.064 | 39.94 | 46.56 | 54.47 | 4.80 | 1.263 | 623 | DATA; COMPUTED |
| H2S | 34.081 | 34.22 | 37.20 | 45.75 | 4.12 | 1.321 | 1004 | DATA; COMPUTED |
| CO | 28.010 | 29.15 | 29.82 | 33.18 | 3.51 | 1.399 | 1041 | DATA; COMPUTED |
| N2O | 44.013 | 38.70 | 45.82 | 54.87 | 4.65 | 1.274 | 879 | DATA; COMPUTED |
| HCN | 27.025 | 35.93 | 41.72 | 50.23 | 4.32 | 1.301 | 1329 | DATA; COMPUTED |
| O3 | 47.998 | 39.33 | 47.22 | 55.00 | 4.73 | 1.268 | 819 | DATA; COMPUTED |
| He, Ne, Ar, Kr, Xe | monatomic | 5R/2 = 20.786 | same | same | 2.5 | 1.667 | He: 5193 | COMPUTED (kinetic theory) |

H2 caveat: below about 200 K the ortho-para ratio and rotational freeze-out change cp; Shomate here is not valid there; cold-giant lapse rates need ortho/para-resolved cp (NOT REACHED in this lane). Real-gas effects (cp rising with pressure) are significant at tens of bar for CO2 and are not in these ideal-gas values.

Dry adiabatic lapse rate `Gamma_d = g / cp_mass`, and the exponent kappa = R/cp. COMPUTED test points (the gravity and composition are the TEST bodies' own):

| TEST case | cp_mass (J/kg/K) | Gamma_d (K/km) | R/cp | Status |
|---|---|---|---|---|
| Earth-like air, 300 K, g 9.81 | 1011 | 9.70 | 0.285 | COMPUTED |
| CO2 at 300 K, g 3.71 (Mars-like) | 846 | 4.39 | 0.223 | COMPUTED |
| 96.5 percent CO2 + N2 at 735 K, g 8.87 (Venus-like) | 1141 | 7.77 | 0.168 | COMPUTED |

Moist adiabat (standard textbook form, Pierrehumbert 2010 Chapter 2, NOT REACHED in that source, so derived and checked): `Gamma_m = g (1 + L r / (R_d T)) / (cp + L^2 r epsilon / (R_d T^2))`, r = epsilon e_s/(p - e_s), epsilon = M_v/M_d. COMPUTED test: T = 288 K, p = 1 atm, L = 2.5e6 J/kg, e_s = 1.70 kPa gives 4.7 K/km (observed tropical-troposphere mean about 5 K/km, RECALLED). Robinson and Catling 2012 replace this by `alpha * Gamma_d` with alpha 0.6-0.9 (READ): Earth 0.6, Venus 0.8, Titan 0.77, Jupiter 0.85, Saturn 0.94, Uranus 0.83, Neptune 0.87 (Robinson and Catling 2014 Table 1, READ). Their text warns this fails for very moist and one-component condensing atmospheres, where T follows Clausius-Clapeyron, T proportional to 1/ln(p/p0).

Latent heats and saturation curves for H2O, CO2, NH3, CH4 are in the condensates lane; the generator already has vapour pressures.

---

## 5. Thermal structure from an analytic atmosphere

### 5.1 Guillot 2010 (semi-grey, irradiated), READ and DATA

arXiv:1006.4702 (NOT 1006.4780; the id given in the brief is a mathematics paper). Variables: `tau = kappa_th * m`, `m = P/g` (column mass), `gamma = kappa_v / kappa_th`, `T_irr = (R_star/D)^(1/2) T_star`, redistribution factor f = 1 (substellar), 1/2 (dayside), 1/4 (global).

Eq. 29 (isotropic irradiation, mu* = 1/sqrt 3, f_K = 1/3, f_H = 1/2):

    T^4 = (3 T_int^4 / 4) (2/3 + tau)
        + (3 T_irr^4 / 4) f [ 2/3 + 1/(gamma sqrt 3) + (gamma/sqrt 3 - 1/(gamma sqrt 3)) exp(-gamma tau sqrt 3) ]

Also READ: collimated version Eq. 27 with mu*; skin temperature `T(tau=0)^4 = T_int^4/2 + T_irr^4 mu* (1 + gamma/(2 mu*))/2` (Eq. 28); `P = tau g / kappa_th` for constant g; global-average solution Eq. 49 (needs the exponential integral E2). The same Eq. 29 is implemented in petitRADTRANS 3.4.0 `physics.py` function `temperature_profile_function_guillot` (DATA; MIT-licensed code, read as text only; note it sets `T_irr = T_eq * sqrt 2`, so f = 1/4 reproduces T_eq).

COMPUTED check with the opacities of Guillot's HD 209458 b example (kappa_th = 1e-2 cm2/g, kappa_v = 4e-3, T_int = 100 K, T_eq = 1469 K, g = 1000 cm/s2, f = 1/4): T = 1473 K at 66.7 mbar (tau = 2/3; Guillot's text puts tau = 2/3 at 65 mbar), 1647 K at 1 bar, 1648 K at 10 bar, so the deep isotherm is about 1650 K. Guillot quotes for other irradiation levels `kappa_v = 6e-3 sqrt(T_irr/2000 K)` cm2/g with kappa_th = 1e-2 (READ).

Limits (READ): for gamma = 1, T^4 = 3 T_int^4 (tau + 2/3)/4 + 0.93 f T_irr^4; high gamma (strong visible absorbers like TiO and VO) can give a negative gradient (inversion); low gamma means deep heating and an isothermal radiative zone.

Validity: radiative-equilibrium only; Guillot's own comparison shows that at pressures above about 10 bar real CIA and electron opacity make kappa proportional to P, so constant kappa underestimates the deep T (READ).

### 5.2 Parmentier and Guillot 2014 (two thermal and two visible opacities), partly READ

arXiv:1311.6597 (Paper I; Paper II, the analytic-versus-numerical comparison, is arXiv:1311.6322 and was not downloaded). The general solution (Eq. 98) is `T^4 = (3 T_int^4/4)[tau + A + B exp(-tau/tau_lim)] + sum_i (3 beta_i T_irr^4 mu*/4)[C_i + D_i exp(-tau/tau_lim) + E_i exp(-gamma_i tau/mu*)]` with coefficients defined in their Eqs. 79-81 (READ only that the form exists; the coefficient expressions were not transcribed: NOT REACHED in this lane). Exact non-irradiated skin temperature for a grey atmosphere: `T_skin^4 = (sqrt 3 / 4) T_int^4` (READ, Eq. 17). The non-grey parameter beta controls the "blanketing" effect: for beta above 0.5 the temperature structure is affected down to tau about 1 and deeper (abstract, READ).

### 5.3 Robinson and Catling 2012 and 2014 (radiative-convective, analytic), READ

arXiv:1209.1833 (ApJ 757, 104) and arXiv:1312.6859 (Nature Geosci. 7, 12). Model: grey thermal, two-stream with diffusivity D = 1.66; `tau(p) = tau0 (p/p0)^n`; convective profile `T = T0 (p/p0)^(alpha (gamma-1)/gamma)`; radiative region above. Using the simplest version (no shortwave attenuation, no internal flux), Eq. 24-25: in the radiative region `sigma T^4 = F_net (1 + D tau)/2` and `F_up = F_net (2 + D tau)/2`; the surface temperature for a purely radiative atmosphere with a blackbody ground is therefore

    T_s^4 = T_e^4 (1 + D tau0 / 2)           [ D = 1.66 ]

(the Eddington-style `1 + 3 tau/4` with D = 3/2; Guillot's gamma to 0 limit with the D = sqrt 3 convention gives the same form). For a convective troposphere the radiative-convective boundary tau_rc is fixed by continuity of temperature and upwelling flux (their Eqs. 27 and 30; I re-derived the continuity condition from their Eqs. 12, 24, 25 and 27 and solved it numerically, see section 5.4). The 2014 paper adds: `p_tp` (tropopause, temperature minimum) obeys `p_tp proportional to tau0^(-1/n)`, about 0.1 bar for worlds with strong-enough pressure-dependent opacity and an absorbing stratosphere. READ Table 1 of the 2014 paper (alpha, absorbed fluxes F_strato and F_tropo, k_strato = tau_sw,s/tau_IR, k_tropo, tau_rc, tau0, tau_tp, model and observed p_tp), TEST points only:

| TEST world | p0 (bar) | T0 (K) | alpha | F_strato + F_tropo (W m-2) | k_strato | k_tropo | tau_rc | tau0 | p_tp model / observed (bar) | Status |
|---|---|---|---|---|---|---|---|---|---|---|
| Earth | 1 | 288 | 0.6 | 7 + 233 | 90 | 0.16 | 0.15 | 1.9 | 0.16 / 0.16 | READ |
| Jupiter | 1 | 166 | 0.85 | 1.3 + 7.0 (+ F_i 5.4) | 90 | 0.06 | 0.34 | 6.3 | 0.10 / 0.14 | READ |
| Saturn | 1 | 135 | 0.94 | 0.41 + 2.04 (+ 2.01) | 180 | 0.03 | 0.44 | 9.2 | 0.066 / 0.08 | READ |
| Titan | 1.5 | 94 | 0.77 | 1.49 + 1.12 | 120 | 0.20 | 4.4 | 5.6 (2.5 at 1 bar) | 0.18 / 0.2 | READ |
| Uranus | 1 | 76 | 0.83 | 0.24 + 0.41 | 220 | 0.08 | 0.62 | 8.7 | 0.070 / 0.11 | READ |
| Neptune | 1 | 72 | 0.87 | 0.09 + 0.18 (+ 0.43) | 580 | 0.20 | 0.41 | 3.0 | 0.075 / 0.1 | READ |

The same paper (2012, Table 1, READ) gives the Venus run: p0 = 92 bar, T0 = 730 K, gamma = 1.3, alpha = 0.8, absorbed 160 W m-2 (Bond albedo 0.76) no stratospheric inversion; solved tau0 = 400 and tau_rc = 1 for n = 1 (p_rc about 0.2 bar), tau0 = 2e5 and tau_rc = 0.1 for n = 2 (p_rc 0.07 bar), vs observed 0.1-0.3 bar. Agreement with SMART line-by-line fluxes is about 10 percent (READ). Titan with tau0 = 3, k1 = 160, n = 4/3 (2ND-HAND).

Physical content for the generator: the dependence `tau proportional to p^n` with n between 1 (Doppler-broadened, abundance independent) and 2 (pressure-broadened, CIA), larger n if the absorber mixing ratio itself rises with depth (water on Earth).

### 5.4 Surface temperature from tau0: solver, table, test points (COMPUTED)

Implementation (about 25 lines): unknowns tau_rc and T0 for given tau0, F, gamma, alpha, n, with m = 4 beta / n, beta = alpha (gamma-1)/gamma, D = 1.66:

1. Continuity of upwelling flux: `[ exp(-D(tau0 - tau_rc)) + D * integral_{tau_rc}^{tau0} (t/tau0)^m exp(-D(t - tau_rc)) dt ] / (tau_rc/tau0)^m = (2 + D tau_rc)/(1 + D tau_rc)`; root-find tau_rc.
2. Continuity of temperature: `sigma T0^4 (tau_rc/tau0)^m = F (1 + D tau_rc)/2`.
3. If no root exists (tau0 small) the atmosphere is purely radiative and `T0^4 = (F/sigma)(1 + D tau0/2)`.
4. For tau0 above about 1e4 the root becomes independent of tau0; scale `T0 proportional to tau0^(beta/n)`.

COMPUTED table of the ratio T_s / T_e (T_e is the airless grey temperature, so the existing generator multiplier is `(T_s/T_e)^4 - 1`), with D tau_rc in the last rows:

| tau0 | 0.01 | 0.1 | 0.3 | 1 | 2 | 5 | 10 | 30 | 100 | 300 | 1000 | 3000 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| air-like gamma 1.4, alpha 0.6, n 2 | 1.00 | 1.01 | 1.03 | 1.08 | 1.14 | 1.23 | 1.30 | 1.43 | 1.59 | 1.75 | 1.94 | 2.13 |
| diatomic gamma 1.4, alpha 0.8, n 2 | 1.00 | 1.01 | 1.03 | 1.10 | 1.17 | 1.30 | 1.41 | 1.59 | 1.83 | 2.07 | 2.38 | 2.70 |
| CO2 gamma 1.3, alpha 0.8, n 1 | 1.00 | 1.01 | 1.04 | 1.13 | 1.23 | 1.44 | 1.63 | 2.00 | 2.50 | 3.06 | 3.83 | 4.69 |
| CO2 gamma 1.3, alpha 0.8, n 2 | 1.00 | 1.01 | 1.03 | 1.09 | 1.15 | 1.25 | 1.33 | 1.47 | 1.64 | 1.82 | 2.03 | 2.25 |
| H2/He gamma 1.4, alpha 0.85, n 2 | 1.00 | 1.01 | 1.03 | 1.10 | 1.18 | 1.32 | 1.43 | 1.64 | 1.89 | 2.16 | 2.50 | 2.86 |

Asymptotic D tau_rc (tau0 above 5): 0.15 (air, alpha 0.6, n 2), 0.36 (diatomic alpha 0.8, n 2), 2.08 (CO2 n 1), 0.19 (CO2 n 2), 0.44 (H2/He). Hence for tau0 above about 5, `T_s/T_e = ((1 + D tau_rc_inf)/2)^(1/4) (tau0/tau_rc_inf)^(beta/n)`. Rough shortcut `D tau_rc = 1` (a statement in Robinson and Catling 2012) is not accurate for these parameter sets and gives Earth's T_s/T_e low by about 3 percent: use the table or the solver.

TEST points (inputs labelled READ from the 2012/2014 papers; outputs COMPUTED here; single-channel model, n = 2 unless stated):

| TEST body | absorbed F (W m-2), T0 (K) | tau0 solved | paper value | tau_rc solved | G = (T_s/T_e)^4 - 1 | Comment | Status |
|---|---|---|---|---|---|---|---|
| Earth | 240, 288 (T_e 255.1) | 1.82 at 1 bar | 1.9 (READ) | 0.086 | 0.625 | 4 percent from the paper value even though the paper uses two shortwave channels | COMPUTED |
| Venus, n = 1 | 160, 730 (T_e 230.5) | 360 | 400 (READ) | 1.25 | 99.6 | tau_rc 1.25 vs 1 in paper; matches the repo's ratio about 101 | COMPUTED |
| Venus, n = 2 | 160, 730 | extrapolated 2e5 gives 763 K | 2e5 (READ) | 0.11 | same | 4 percent high at the paper's tau0 | COMPUTED |
| Titan | 2.61, 94 (T_e 82.4) | 1.58 at 1.5 bar | 5.6 (READ) | 0.175 | 0.69 | factor 3.5 below the paper: a haze world with a hot stratosphere needs the two-channel model and n varying with p (McKay: n 4/3) | COMPUTED |
| Mars-like (inputs RECALLED: S = 586 W m-2, Bond albedo 0.25, T_s 210-220 K) | 110, 215 (T_e 209.8) | 0.12 (radiative regime) for 215 K, 0.25 for 220 K | n/a | none | 0.10 to 0.21 | mean annual T_s is within a few K of T_e, so tau0 is tiny and the greenhouse is small | COMPUTED with RECALLED inputs |

Inverse use: given T_s, P_s, S, g and an albedo budget (section 7), solve `tau0` by bisection (monotone). The same `tau0` then fixes the sky radiance (5.5) and the tropopause (`p_tp ~ tau0^(-1/n)`).

**Built (GEN.812, `physics/infrared_opacity.py`; the module and its cloud depth `CLOUD_LONGWAVE_TAU0` were removed in GEN.850, which put the radiative-convective model's fitted table in their place, `composition/greenhouse-table.md`, with the cloud as IR12).** The module holds the surface's vapour pressure and thins it with height (Manabe and Wetherald's humidity, the troposphere's `T_s (p/p_s)^beta`), so an Earth twin whose vapour stands at 0.77 of saturation at its own temperature holds the Earth's 25 kg/m2 of water and not three times that. With the surface vapour pressure of 1.3 kPa the module gives `tau0` 1.70 and `G` 0.625 for the Earth (targets 1.8 and 0.625). It is a clear-sky greenhouse (outgoing 266 W/m2 at 288 K, the measured clear-sky value), so the all-sky albedo needs the cloud longwave effect (section 7.3, +26 to +30 W/m2): one grey cloud depth, `CLOUD_LONGWAVE_TAU0` 0.35, scaled with the surface humidity (GEN.725 rows GH1 to GH4), puts the Earth twin at its own fixed point at 288 K. The runaway starts between 255 and 285 W/m2 absorbed (Goldblatt 2013: 282 for a saturated air; Kopparapu 2013: 262 at 1.1 times the present flux). Pressure broadening, checked against Byrne and Goldblatt 2014b (1 percent CO2 on 0.5, 1 and 2 bar of N2: 15, 26 and 35 W/m2): 16, 20 and 27 W/m2, a little weak, not strong.

### 5.5 Downwelling infrared and effective sky temperature (COMPUTED from the same model)

Downwelling flux at the ground: in a purely radiative atmosphere `F_dn(tau0) = F D tau0/2` (Eq. 26 of R&C 2012, READ); with a convective layer `F_dn(tau0) = F_dn(tau_rc) exp(-D(tau0 - tau_rc)) + D sigma T0^4 integral_{tau_rc}^{tau0} (t/tau0)^m exp(-D(tau0 - t)) dt` (their Eq. 14 with the boundary at tau_rc). Effective sky temperature `T_sky = (F_dn/sigma)^(1/4)`, emissivity `eps = F_dn/(sigma T_s^4)`.

| TEST body | F_dn (W m-2) | eps_sky | T_sky (K) | observed comparison | Status |
|---|---|---|---|---|---|
| Earth, tau0 = 1.82 | 326 | 0.83 | 275 | about 340 W m-2 global mean back radiation (RECALLED, AR6 energy budget figure; not located in the PDF text) | COMPUTED |
| Earth, tau0 = 1.9 | 334 | 0.84 | 277 | same | COMPUTED |
| Venus, n = 1 | 16087 | 0.999 | 729.8 | black-body-like at the surface | COMPUTED |
| Titan | 3.4 | 0.77 | 88 | no check | COMPUTED |
| Mars-like (215 K, 220 K) | 11, 23 | 0.09, 0.17 | 100, 120 | no check | COMPUTED with RECALLED inputs |

Empirical clear-sky formulae for Earth air only (not for other compositions): `F_dn = 5.31e-13 T_air^6` W m-2, and `eps = 1 - (1 + w) exp(-sqrt(1.2 + 3 w))`, `w = 46.5 e/T` with e in hPa (both RECALLED, to verify). Outside Earth air use 5.5.

---

## 6. Rosseland and Planck means, and the inverse problem

### 6.1 Freedman et al. 2014 Rosseland fit (H2/He-dominated, solar to 50x metallicity)

arXiv:1409.0026 (ApJS 214, 25). `kappa_gas = kappa_lowP + kappa_highP`, cm2/g, T in K, P in dyne/cm2, `met` = [M/H]. Equations as extracted from the PDF (the layout of Eq. 4 is garbled in extraction; I reconstructed it with the form that the paper says it modifies, so mark: reconstructed):

    log10 kappa_lowP  = c1 atan(log10 T - c2) - c3/(log10 P + c4) * exp[(log10 T - c5)^2] + c6 met + c7
    log10 kappa_highP = c8 + c9 log10 T + c10 (log10 T)^2 + log10 P (c11 + c12 log10 T)
                        + c13 met [1/2 + (1/pi) atan((log10 T - 2.5)/0.2)]

Table 2 coefficients (READ): c1 10.602, c2 2.882, c3 6.09e-15, c4 2.954, c5 -2.526, c6 0.843, c7 -5.490 (all T); for T below 800 K: c8 -14.051, c9 3.055, c10 0.024, c11 1.877, c12 -0.445, c13 0.8321; for T above 800 K: c8 82.241, c9 -55.456, c10 8.754, c11 0.7048, c12 -0.0414, c13 0.8321. The paper states the fit reproduces the tables within a factor 2 for 90 percent of points at solar metallicity (80 percent at [M/H] = +1.7), typical error about 50 percent. COMPUTED sanity values (solar, cm2/g): 1 bar: 1.2e-2 (100 K), 1.7e-2 (200 K), 2.1e-2 (300 K), 2.7e-2 (500 K), 1.4e-2 (1000 K), 1.0e-2 (1500 K), 3.1e-2 (2000 K); 100 bar: 1.1 (100 K), 0.73 (300 K), 0.20 (1000 K), 0.16 (2000 K); the value rises linearly with P at 100 K between 1 and 100 bar, as expected for H2 CIA. Status: READ coefficients, COMPUTED values, equation form reconstructed.

Validity: giant planets and brown dwarfs, 75-4000 K (their abstract and search summary), 1 ubar to 300 bar, chemical-equilibrium cloud-free gas; the paper cautions Rosseland versus Planck choice and weighting temperature (their section 4). Not applicable to CO2-, N2- or H2O-dominated air: NOT REACHED for those; see to-research list. Tables are at CDS (J/ApJS/214/25).

### 6.2 What a single grey tau0 means, and how to build it from composition

Three routes, in order of defensibility:

1. Equivalent grey thickness of the real spectrum (Koll 2019 Eq. 13): `tau_LW = -ln[ integral e^(-tau_nu) B_nu(T_s) dnu / integral B_nu dnu ]`. READ definition; needs per-gas spectra; correct for thin atmospheres (surface-dominated emission). Koll's Fig. 4 (READ text) shows for pure gases at an LHS 3844 b-like T_s: H2O optically thick above about 1e-2 bar; CO2 and O3 above about 1 bar (tau_LW of order unity at 1 bar for CO2); CO optically thin even at 1 bar; at high pressure CIA (p squared) takes over; band-edge closure gives log dependence on p.
2. Sum of CIA (computed here, section 3.1) plus line contributions from published k-tables, calibrated by a precomputed radiative-convective grid (offline, using exo_k or petitRADTRANS with ExoMol/HITRAN tables: see to-research list).
3. Calibrate `tau0(P)` against TEST bodies. This is allowed only as a test, not as a preset: Earth about 1.9 at 1 bar, Venus about 400 at 92 bar (n 1), Titan about 2.5 at 1 bar, Jupiter 6.3, Saturn 9.2 (all at 1 bar for the giants).

### 6.3 Inverse problem recipe (given stored T_s, P_s, S, g, mix)

1. Compute T_e for a trial albedo A from section 7: `T_e = [(1-A) S f/(4 sigma)]^(1/4)` with f from section 8 if the planet is tidally locked.
2. Choose gamma and alpha from the dominant gas (section 4 table; alpha from R&C 0.6-0.9), and n from the regime (1 thin/Doppler, 2 pressure-broadened or CIA).
3. Solve section 5.4 for `tau0` so that T0 = T_s at p0 = P_s. Report `G = (T_s/T_e)^4 - 1` for the existing field.
4. Non-uniqueness: one equation, two unknowns (A, tau0). The pair family is `T_s = F(A, tau0)`; fix A from composition (Rayleigh plus surface plus cloud) and solve tau0, or fix tau0 from composition (6.2) and solve A. The residual between the two routes is the diagnostic that the draw is inconsistent.
5. A stored T_s above the runaway limit (section 9) with liquid water must be rejected or turned into a steam atmosphere.

---

## 7. Albedo: Rayleigh, surface, clouds, haze

### 7.1 Rayleigh cross-sections (DATA plus RECALLED depolarisation)

Computed from the refractive-index dispersion fits of the refractiveindex.info database (public domain, CC0; files `main/N2/nk/Peck-0C.yml`, `main/CO2/nk/Bideau-Mehu.yml`, `main/H2/nk/Peck.yml`, `main/He/nk/Borzsonyi.yml`, `main/Ar/nk/Bideau-Mehu.yml`, all at 273.15 K and 1 atm) with `sigma = 24 pi^3 (n^2-1)^2 / (lambda^4 N0^2 (n^2+2)^2) * F_K`, N0 = 2.6868e19 cm-3. The King depolarisation factors F_K are RECALLED: N2 1.034, CO2 1.1364, H2 1.0, He 1.0, Ar 1.0.

| Gas | (n-1) at 0.55 um, STP | sigma (cm2) at 0.4 / 0.55 / 0.7 / 1.0 um | Status |
|---|---|---|---|
| N2 | 2.99e-4 | 1.7e-26 / 4.6e-27 / 1.7e-27 / 4.1e-28 | DATA; COMPUTED; F_K RECALLED |
| CO2 | 4.50e-4 | 4.3e-26 / 1.2e-26 / 4.3e-27 / 1.0e-27 | DATA; COMPUTED; F_K RECALLED |
| H2 | 1.39e-4 | 3.6e-27 / 9.7e-28 / 3.6e-28 / 8.6e-29 | DATA; COMPUTED |
| He | 3.44e-5 | 2.2e-28 / 5.9e-29 / 2.3e-29 / 5.4e-30 | DATA; COMPUTED |
| Ar | 2.82e-4 | 1.5e-26 / 4.0e-27 / 1.5e-27 / 3.6e-28 | DATA; COMPUTED |

Column Rayleigh optical depth `tau_R = sigma P N_A / (m_bar g)`. COMPUTED at 0.55 um: Earth-like 1.013 bar air 0.099 (matches the standard value of about 0.1, so the chain is sound); 92 bar CO2 with g 8.87: 16.4; 6.4 mbar CO2 with g 3.71: 0.0027; 1.47 bar N2 with g 1.35: 1.08; 1 bar H2/He (86/14) with g 24.8: 0.089. Water vapour Rayleigh constants: Kopparapu 2013 adopt `sigma_R,H2O(lambda) = 4.577e-... * (6 + 3D)/(6 - 7D) ...` (garbled in extraction; READ that it exists; the Goldblatt 2013 supplementary table 1 with new constants was not reached). Neighbouring lane: sky colour.

### 7.2 Bond albedo of a conservative Rayleigh layer over a surface

Closure (COMPUTED, simple): layer reflectance `R = (3/4) tau/(1 + (3/4) tau)` (Eddington, g = 0), transmittance `1 - R` (no absorption), surface albedo `A_s`, `A = R + (1-R)^2 A_s / (1 - R A_s)`, solar-weighted over 0.3-2.5 um with a blackbody. Accuracy about 20 percent against exact multiple-scattering is expected, not verified here. Absorption by CO2, H2O, CH4 near-IR bands and UV absorbers is ignored (it lowers the albedo, mostly for cool stars).

| Atmosphere | P (bar) | Sun-like (5772 K), A_s = 0 | A_s = 0.3 | 3000 K blackbody, A_s = 0 | A_s = 0.3 | Status |
|---|---|---|---|---|---|---|
| N2, g 9.81 | 0.1 | 0.009 | 0.304 | 0.001 | 0.301 | COMPUTED |
| N2 | 1 | 0.070 | 0.337 | 0.011 | 0.306 | COMPUTED |
| N2 | 10 | 0.297 | 0.475 | 0.079 | 0.342 | COMPUTED |
| N2 | 92 | 0.621 | 0.699 | 0.303 | 0.478 | COMPUTED |
| CO2, g 8.87 | 1 | 0.108 | 0.358 | 0.019 | 0.309 | COMPUTED |
| CO2 | 10 | 0.377 | 0.528 | 0.116 | 0.364 | COMPUTED |
| CO2 | 30 | 0.540 | 0.641 | 0.227 | 0.430 | COMPUTED |
| CO2 | 92 | 0.695 | 0.754 | 0.385 | 0.531 | COMPUTED |
| H2/He 86/14, g 24.8 | 1 | 0.064 | 0.333 | 0.010 | 0.305 | COMPUTED |
| H2/He | 10 | 0.279 | 0.464 | 0.071 | 0.338 | COMPUTED |

Consistency checks: the Earth-like clear-sky Rayleigh albedo over a dark surface of 0.07 agrees with the usual value of about 0.06-0.07 (RECALLED); Kopparapu 2013 (READ) states that Rayleigh scattering by condensing and dense CO2 raises albedo at the outer habitable zone so that Seff goes through a minimum, and that for a 2600 K star Rayleigh never matters (as in the 3000 K column here). Byrne and Goldblatt 2014b (READ) quote a shortwave forcing of -18 W m-2 for a 1-bar-CO2 (1.0 mole fraction) addition, which is the same effect.

### 7.3 Clouds, haze, anti-greenhouse

| Item | Value | Source | Status |
|---|---|---|---|
| Titan haze anti-greenhouse | haze absorbs sunlight but is transparent in the thermal infrared; removing it would raise the surface temperature by over 20 K; N2, CH4 and H2 pressure-induced opacity provide the greenhouse | McKay, Pollack and Courtin 1991, Science 253, 1118 (abstract at ntrs.nasa.gov, id 19910067426) | READ (abstract) |
| Titan split 9 K cooling and 21 K warming | as reported in a search summary of the same paper | secondary summary | SNIPPET (do not rely) |
| Titan in the analytic model | F_strato 1.49 W m-2 absorbed high up, F_tropo 1.12 below, k_strato 120, k_tropo 0.20, tau_rc 4.4 (radiative zone reaches to 0.18 bar tropopause and deeper, because of the weak shortwave penetration) | Robinson and Catling 2014 Table 1 | READ |
| One-channel estimate of haze effect | feeding only F_tropo = 1.12 into the one-channel model with the same tau0 gives T_s = 76 K, with half of the stratospheric heat returned as downward IR 86 K, with the full 2.61 W m-2 94 K: a one-channel model overstates the cooling; use the two-channel formulas (their Eqs. 18-21) | COMPUTED here | COMPUTED |
| Earth cloud radiative effect | short-wave about -45 to -50 and long-wave about +26 to +30 W m-2, net about -20 W m-2 | recalled from CERES/AR6 energy-budget figures; not located in the AR6 PDF text | RECALLED |
| Venus | Bond albedo 0.76; absorbed about 160 W m-2 | Robinson and Catling 2012 text | 2ND-HAND |
| Cloud albedo of a conservative cloud, optical depth tau_c, asymmetry g | use the same two-stream R with `(1-g) tau_c` in place of tau; typical water cloud g about 0.85 (RECALLED) | derivation | COMPUTED form; g RECALLED |
| Water cloud and clear-sky coverage in 1D HZ models | cloud-free 1D limits: the Earth-like case falls right on the runaway limit in the Kopparapu 2014 model, "reflecting 1D models' inability to account for relative humidity and clouds, which move the inner HZ to higher stellar fluxes" | Kopparapu et al. 2014 text | READ |

Aerosol optical constants (tholins, soot, sulfate, silicate, iron, KCl, ZnS, NH4SH and so on) belong to the aerosol lane; here only their role (shortwave absorber means high gamma and an inversion, Guillot 2010) is used.

---

## 8. Day-night contrast and heat redistribution

Koll 2019 (arXiv:1907.13145, READ). Let `X = tau_LW^(1/3) (p_s/1 bar)^(2/3) (T_eq/600 K)^(-4/3)`. The redistribution factor (the f of the energy budget, with `T_day,eff^4 = f (1-A) S/sigma`, f = 2/3 no redistribution, 1/4 uniform):

    f = 2/3 - (5/12) * X / (k + X)                      (his Eq. 10)

Limits READ: f goes to 2/3 as tau_LW and p_s go to zero and to 1/4 as they go to infinity. `k` collects the planet properties (his definition with the drag coefficient, the dissipation efficiency and the heat-transport parameters); fitted to idealised GCM runs: k = 1.2 for TRAPPIST-1b, 1.9 for GJ 1132 b, 2.3 for LHS 3844 b (READ). Nightside flux `F_n` (his Eq. 8) and observed dayside flux Eq. 9 (the 8/3 and 5/3 factors) are in the paper; the scaling is for the surface-pressure range about 0.01 to 10 bar, with `tau_LW` from Koll's Eq. 13. Conclusions READ: redistribution becomes efficient above about 1 bar for high mean-molecular-weight air; much less efficient for H2, and strongly gas dependent (steam and CO atmospheres behave differently); typical k varies only by a factor 2 across planets.

Wordsworth 2015 (arXiv:1412.5575, READ): in an idealised grey GCM the boundary layer, not the large-scale circulation, controls the day-night temperature difference; theoretical expressions for the temperature difference and surface wind reproduce the GCM; empirical collapse formula `log10(p_crit/1 Pa) = c1 F_s^3 + c2 F_s^2 + c3 F_s + c4` with F_s the stellar flux in W m-2 (as understood from his section 3; units to verify): for 1 M_E c = (-0.84e-10, 0.73e-6, -0.0022, 6.01), for 10 M_E c = (-0.81e-10, 0.63e-6, -0.0017, 6.34); critical CO2 collapse pressure at Earth's stellar flux about 0.14 bar, about five times previous grey-gas estimates; results are an upper limit because clouds, trace gases, ocean transport and tidal heating make atmospheres more stable. COMPUTED evaluation of the 1 M_E formula at F_s = 1000 W m-2: log10 p = -0.84e-10*1e9 + 0.73e-6*1e6 - 2.2 + 6.01 = -0.084 + 0.73 - 2.2 + 6.01 = 4.456, i.e. 2.9e4 Pa = 0.29 bar (at 1000 W m-2; the 0.14 bar quoted in the abstract refers to Earth's flux, 1360 W m-2; the formula should return about that near 1360 W m-2, which I did not evaluate: to verify units).

Column radiative time scale (COMPUTED from the definition, simple): `t_rad = cp p_s / (g 4 sigma T_e^3)`, the e-folding time of the whole column heat capacity `cp p/g` against linearised emission. Earth-like (cp 1005, p 1.013e5 Pa, g 9.81, T_e 255 K): 2.8e6 s, about 32 days. Mars-like (cp 846, p 636 Pa, g 3.71, T_e 210 K, inputs RECALLED): 7e4 s, about 0.8 sol, so a thin atmosphere cannot buffer the diurnal cycle: the ground sets the day-night contrast through its thermal inertia `I = sqrt(k rho c)`. A table of regolith, rock, ice and ocean-mixed-layer thermal inertia is NOT REACHED in this lane (Mars TES maps by Putzig and Mellon, lunar regolith, and so on); the formula, the column heat capacity and Koll's f are the part that can be coded now.

Pierrehumbert and Hammond 2019 (Annu. Rev. Fluid Mech. 51, 275; the weak-temperature-gradient and the substellar-to-nightside scalings) was not reachable (no arXiv copy found; the related arXiv:1810.11231 and arXiv:2009.00358 of Hammond and Pierrehumbert exist but were not read). Tidal heating link (brief): the surface heat flux from tides enters as an internal flux F_i in the Robinson and Catling formulas (their Eq. 17, READ) and as a floor under T_s; the formula `E_dot = (21/2)(k2/Q) G M_star^2 R^5 n e^2/a^6` (Peale et al. 1979) is RECALLED and belongs to the interior lane.

---

## 9. Runaway greenhouse, moist greenhouse, habitable-zone limits

### 9.1 Radiation limits (READ, mostly arXiv)

| Quantity | Value | Source | Status |
|---|---|---|---|
| Nakajima/Simpson limit (moist troposphere, outgoing thermal flux limit, pure steam, g about 10 m s-2) | 282 W m-2 (spectrally resolved LBL), vs 310 W m-2 (Kasting 1988) and 290 W m-2 (Nakajima 1992 grey) | Goldblatt et al. 2013 (author copy at dspace.library.uvic.ca; Goldblatt and Watson 2012 for the older values, arXiv:1201.1593) | READ |
| Kopparapu 2013 steam-atmosphere limit | 291 W m-2, matching Fig. 4.37 of Pierrehumbert 2010 for a pure-water-vapour planet at 10 m s-2 | arXiv:1301.6674 footnote | READ |
| Absorbed solar flux of a runaway pure-steam, cloud-free atmosphere at present solar constant | 294 W m-2 (vs 222 W m-2 earlier) so a cloudless steam atmosphere could be stable at Earth's flux: bistability | Goldblatt et al. 2013 | READ |
| Komabayashi-Ingersoll limit (moist stratosphere) | about 385 W m-2 in the example; occurs only if tropopause T rises to 255 K or more; effectively impossible to reach on Earth | Goldblatt and Watson 2012, section 2 | READ |
| Moist greenhouse onset | stratospheric H2O mixing ratio about 3e-3 reached at surface temperature 340 K in the Kopparapu 2013 model; water-loss time scale then approaches the age of the Earth | arXiv:1301.6674 | READ |
| Runaway limit without H2O continuum absorption | moves from 0.97 AU (with continuum) to 0.76 AU | arXiv:1301.6674 | READ |
| Kasting 1993 reference limits for the Sun | water loss 0.95 AU, runaway 0.84 AU, maximum greenhouse 1.67 AU, recent Venus 0.75 AU, early Mars 1.77 AU | as quoted in Kopparapu 2013 | 2ND-HAND |
| Hot desert world minimum distance | 0.38 AU for a solar-like star (1 percent relative humidity, higher albedo); with 1-10 bar | Zsom et al. 2013 | READ |
| Grey/analytic summary of the physics | the OLR of a saturated atmosphere approaches `sigma T_(tau=1)^4` with T_(tau=1) set by the Clausius-Clapeyron curve; with background non-condensing gas the limit is overshot (lapse rate increases as the water fraction falls) | Goldblatt and Watson 2012 | READ |

### 9.2 Habitable-zone flux coefficients (code-ready)

`Seff = Seff_sun + a T* + b T*^2 + c T*^3 + d T*^4`, `T* = Teff - 5780 K`, valid 2600 K to 7200 K, cloud-free N2-CO2-H2O 1D model; distance `d = sqrt(L/Lsun / Seff)` AU.

Kopparapu et al. 2013 (arXiv:1301.6674, Table 3, READ), 1 Earth mass, 1 bar N2 background:

| Limit | Seff_sun | a | b | c | d |
|---|---|---|---|---|---|
| Recent Venus | 1.7753 | 1.4316e-4 | 2.9875e-9 | -7.5702e-12 | -1.1635e-15 |
| Runaway greenhouse | 1.0512 | 1.3242e-4 | 1.5418e-8 | -7.9895e-12 | -1.8328e-15 |
| Moist greenhouse | 1.0140 | 8.1774e-5 | 1.7063e-9 | -4.3241e-12 | -6.6462e-16 |
| Maximum greenhouse | 0.3438 | 5.8942e-5 | 1.6558e-9 | -3.0045e-12 | -5.2983e-16 |
| Early Mars | 0.3179 | 5.4513e-5 | 1.5313e-9 | -2.7786e-12 | -4.8997e-16 |

COMPUTED check for the Sun: 0.751, 0.975, 0.993, 1.706, 1.774 AU (the paper's Table: 0.75, 0.97, 0.99, 1.70, 1.77, READ). Values at Teff 4000 K: Seff 1.561, 0.891, 0.892, 0.256, 0.237; at 3000 K: 1.494, 0.864, 0.853, 0.226, 0.209 (COMPUTED).

Kopparapu et al. 2014 (arXiv:1404.5292 v2, Table 1, READ; the abstract of a published erratum exists, check it, RECALLED): the runaway limit depends on planetary mass because the N2 background pressure is scaled with radius and so the H2O column changes.

| Mass | Seff_sun (runaway) | a | b | c | d |
|---|---|---|---|---|---|
| 0.1 M_E | 0.99 | 1.209e-4 | 1.404e-8 | -7.418e-12 | -1.713e-15 |
| 1 M_E | 1.107 | 1.332e-4 | 1.58e-8 | -8.308e-12 | -1.931e-15 |
| 5 M_E | 1.188 | 1.433e-4 | 1.707e-8 | -8.968e-12 | -2.084e-15 |

Maximum greenhouse (all masses) 0.356, 6.171e-5, 1.698e-9, -3.198e-12, -5.575e-16; recent Venus 1.776, 2.136e-4, 2.533e-8, -1.332e-11, -3.097e-15; early Mars 0.32, 5.547e-5, 1.526e-9, -2.874e-12, -5.011e-16. The paper reports the Earth then sits right on the 1D runaway limit; the inner HZ flux changes by about 10 percent (0.1 M_E) and about 7 percent (5 M_E) relative to Earth.

Ramirez and Kaltenegger 2017 (arXiv:1702.08618, Tables 1 and 2, READ): H2-enriched (1 to 50 percent) limits: outer-edge Seff_sun 0.3129, 0.2671, 0.2452, 0.2271, 0.2009, 0.1723 and inner-edge 0.9975, 0.9946, 0.9884, 0.9775, 0.9655, 0.9377 for 1, 5, 10, 20, 30, 50 percent H2, with the same polynomial in T*; coefficient columns are in the paper (50 percent H2 outer: 1.4029e-5, -3.2318e-10, -2.8689e-13, 4.2926e-17; inner: 1.0123e-4, 7.9833e-9, -2.2204e-12, 1.2098e-16). COMPUTED distances for the Sun: outer edge 2.41 AU (50 percent H2), 2.02 AU (10 percent); inner 1.03 AU (50 percent).

Use in the generator: these give flux thresholds for runaway, moist-greenhouse and maximum-greenhouse limits as functions of stellar Teff, so a stored habitable-zone class and T_s can be cross-checked; they do not give surface temperatures.

---

## 10. Which results can be code, which need full radiative transfer

Can be turned into code now (all inputs available after the composition draw, no radiative transfer needed):

1. T_e and f redistribution (Koll Eq. 10; Guillot f).
2. Radiative-convective T_s from tau0 (section 5.4, solver or table), inverse by bisection (6.3), and the repo's G as a derived readout.
3. Sky radiance and tropopause pressure from tau0 (5.5, `p_tp ~ tau0^(-1/n)`).
4. CIA part of tau0 (section 3.1 formula with a small table of kP(T) per pair; for fidelity keep the HITRAN blocks in a data file and compute tau_eq by the Koll definition).
5. CO2, CH4, N2O forcing as perturbations (Etminan, Byrne and Goldblatt forms) for Earth-like profiles only, and the ordering of strong minor gases (section 2).
6. Rayleigh albedo with the two-stream closure (7.2) and layer adding for clouds and surface.
7. Habitable-zone flux thresholds (9.2) and the runaway and moist limits as numbers.
8. Gas cp, gamma, dry and moist lapse rates, column heat capacity and radiative time (section 4 and 8).
9. Guillot 2010 profile for hot, irradiated atmospheres (5.1) with T_int, kappa_th and gamma as inputs.

Need full radiative transfer (or precomputed grids from it):

1. The line contribution to tau0 for CO2, H2O, CH4, NH3, SO2, O3 and mixtures, including overlap, pressure broadening and the water window; the generator should read a lookup grid tau0(P_s, T_s, x_CO2, x_H2O, x_CH4, x_H2, ...) built offline from correlated-k tables.
2. Real shortwave heating (stratospheric inversions from O3, SO2, haze, TiO/VO): the two-channel analytic form gives the trend only.
3. The runaway and bistability thresholds as functions of composition (the 282 W m-2 Simpson-Nakajima limit holds for pure steam at 10 m s-2; mixtures and other gravities need runs).
4. Non-grey effects (beta above 0.5), super-Lorentzian wings of CO2, the H2O continuum at tens to hundreds of bar, supercritical steam.
5. Clouds as an interactive component and their longwave effect, and 3D effects (relative humidity, wind-driven dry zones).

---

## 11. "Boss must research" (not reachable or only partly reached)

2. Original papers behind 2ND-HAND rows: Myhre et al. 1998 (GRL 25, 2715; agupubs blocked), Kasting et al. 1993 (Icarus 101, 108), Nakajima et al. 1992 (J. Atmos. Sci. 49, 2256), Wordsworth et al. 2010 (CO2 CIA, Icarus), Paynter and Ramaswamy 2011 (BPS water continuum), Goldblatt et al. 2013 supplement (water-vapour Rayleigh constants, Table 1 of the supplement), Pierrehumbert and Hammond 2019 review.
3. Parmentier and Guillot 2014 coefficient formulas (their Eqs. 79-81) and Parmentier et al. 2015 (arXiv:1311.6322) Paper II with the fitted two-band opacity tables for T_irr; needed to give hot Jupiter and lava-world temperature profiles with realistic non-grey behaviour.
4. Rosseland and Planck means for CO2-dominated, N2-O2-H2O (Earth-like), and steam atmospheres, as functions of T and P (Freedman 2014 covers H2/He with metals only). Build offline from HITRAN/HITEMP/ExoMol line lists and CIA; compare against Venus and Earth profiles.
5. The line and continuum spectra needed to compute tau0 for CO2 and H2O: HITRAN line parameters (hitran.org line-by-line is reachable but requires a login), or ExoMol/HITEMP cross-section tables for hot cases; MT_CKD user guide (github.com AER-RC/MT_CKD_H2O docs) to verify the radiation-term convention used in 3.2.
6. CIA gaps: CO2-N2, CO2-O2 and N2-N2 above 400 K, H2-H2 below 200 K (needed for ice giants; alternates at hitran.org/data/CIA/alternate not parsed), H2-He alternates, NH3 and H2O self-CIA beyond MT_CKD, hot H2-H, H-He; CIA of C2H6 and C2H2 with N2 for Titan-like worlds.
7. Surface thermal inertia and heat capacity table for regolith, rock, ice and ocean mixed layers, with sources (Mars TES, Moon Diviner, laboratory).
8. Observational cloud radiative effect for Earth (CERES EBAF) with numbers, and cloud albedo versus optical depth and asymmetry parameter for H2O, NH3, H2SO4, CO2, KCl, ZnS, silicate and iron clouds (aerosol lane).
9. Two-channel analytic model with a hot stratosphere for non-Titan haze worlds: Titan is the only calibrated case; exoplanet haze worlds have no calibration. Also the 2012 Robinson and Catling Titan row (the table text was garbled in extraction, only the 2014 version was read in full).
10. Ortho-para-resolved H2 cp at 50-200 K and real-gas cp for dense CO2, N2, H2O (steam) above 10 bar.
11. Water-world and steam-atmosphere radiative structure above 400 K at tens to hundreds of bar with supercritical water; the Goldblatt 2013 curves only go to 2000 K surface temperature and are figures, not tables.
12. Rocky exoplanet tropopause: Robinson and Catling 2014 has only Solar System calibration; no data for CO2-dominated or O3-free rocky exoplanets.

---

## 12. Recalled or to verify

| Item | Value as recalled | Why it matters | Check against |
|---|---|---|---|
| CO2, O3, N2O, SO2, OCS, NO2, HCl, HF, CO band centres | CO2 667, 2349, 1388/1285, hot bands 960 and 1063; O3 1042, 1110, 705; N2O 589, 1285, 2224; SO2 518, 1151, 1362; OCS 520, 859, 2062; NO2 750, 1318, 1618; HCl 2886; HF 3962; CO 2143; SiO 1242 | band list table | NIST webbook other masks, HITRAN molecule pages |
| Pure rotation band of H2O below about 600 cm-1 | stated qualitatively | window definition | HITRAN or MT_CKD plot |
| King depolarisation factors | N2 1.034, CO2 1.1364, H2 1.0, He 1.0, Ar 1.0 | Rayleigh cross-sections | Sneep and Ubachs 2005 (J. Quant. Spectrosc. Radiat. Transfer 92, 293) |
| Mars-like test inputs | S = 586 W m-2, Bond albedo 0.25, mean T_s 210-220 K, P 636 Pa, g 3.71 | Mars test row, radiative time | any Mars fact sheet (PDS is now reachable) |
| Earth moist lapse rate observed about 5 K/km; Earth back radiation about 340 W m-2 | test targets | section 4 and 5.5 checks | AR6 Ch.7 energy budget figure |
| Swinbank and Prata clear-sky emissivity formulas | as quoted | sky radiance for Earth air | original papers |
| Earth cloud radiative effect about -45/+28 W m-2 | test only | section 7.3 | CERES EBAF documentation |
| Water cloud asymmetry parameter 0.85 | cloud albedo example | section 7.3 | Mie calculation |
| Peale tidal heating formula | as quoted | tidal link | Peale, Cassen and Reynolds 1979 |
| HITRAN CIA unit cm^5 molecule^-2 | convention | section 3.1 | hitran.org documentation page for CIA (Karman et al. 2019) |
| MT_CKD radiation-term convention `nu tanh(h c nu/2kT)` and density scaling by n/n_ref | used in 3.2 | sigma values | MT_CKD user guide PDF |
| Existence and content of the Kopparapu 2014 erratum | not checked | whether Table 1 numbers changed | ApJL 787 L29 and its erratum |
| Wordsworth 2015 collapse-formula units (F_s) and the 0.14 bar quote | formula evaluated at 1000 W m-2 only | section 8 | his Fig. 11 and Table 2 |
| B&G 2014a CH4 high-range equation | garbled in extraction | forcing form above 2.5 ppm | paper, Table 2 |
| Freedman Eq. 4 layout | reconstructed (Valencia form) | opacity values | paper Eq. 4 |

---

## 13. Cross-links to neighbouring parts

- Composition and mixing ratios (gas lists, mean molar mass): the bulk-composition and atmosphere-mix parts; mean molar mass `m_bar` here replaces `atm_molar_density` in `scale_height_m` and in `tau0` (`tau = kappa P/g`, `H = k T/(m_bar g)`), see the trace report items 1 to 4.
- Condensates, saturation vapour pressures and latent heats: the condensate/ice part (moist adiabat, cloud base, water loss inventory).
- Aerosols, hazes, tholins, silicate and iron clouds and their optical constants (single-scattering albedo, asymmetry, mass absorption coefficients): the aerosol part; this file only supplies the albedo layering and the anti-greenhouse sense.
- Atmospheric colour and Rayleigh scattering of the sky: the sky-colour part (shares the Rayleigh cross-section table, section 7.1).
- Escape, water loss and XUV: the escape part (moist greenhouse boundary, section 9, hands over to energy-limited and diffusion-limited escape).
- Radiation environment at the surface (UV, ozone, atmospheric and surface radiation): the radiation part; O3 and SO2 shortwave absorption affect the stratospheric temperature here.
- Interior heat flux and tidal heating: the interior part (the internal flux F_i in section 5.3).
- Tidal locking and rotation: the orbit and rotation part (f in section 8, collapse pressure).
