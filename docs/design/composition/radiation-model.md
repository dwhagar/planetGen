# Radiation model: transmission and greenhouse by band, shielding by column, and the ground's own dose

Status: GEN.459 (the radiogenic module, `physics/radiogenic.py`) and the ground dose of section 6 are built (2026-10-23, epoch 67, schema v136); the rest waits on its items. Approved by Boss as recommended (2026-10-17, 9:40 to 9:55 PM EDT). Item GEN.436 (owner Blocker 1, built after TEST.186; GEN.437 and GEN.472 wait on it). Evidence date 2026-10-17. This file is docs/design/composition/radiation-model.md and changes no code. Revised after Boss's papers report (docs/design/research/boss-papers-2026-10-17.md, read from origin/main, sections 2.4, 2.5, 2.8, 3 and 4): where a paper contradicts a second-hand value the paper wins and both values are stated.

Short names. master = docs/design/composition/atmosphere-composition.md. Part NN = docs/design/composition/atmosphere-composition-parts/NN-*.md (08A, 08B and 10B are the lane checks of Boss's two research notes). Body part E, F, G2 = docs/design/composition/body-composition-parts/. CM = composition-model.md. GFP = gap-filling-plan.md. SC = atmosphere-sources-check.md. Code = src/planetgen/physics/.

Verification states (CM 4b, Decided 3:02 AM EDT): **verified** (read in the source, data, or computed by a lane), **partly verified** (second-hand or a PARTLY verdict, with its limit named), **open** (no source; a seeded bracket meanwhile), **stand-in** (a measured Solar System shape used until a derived law replaces it). Every number here cites a source or comes from the code in section 12; no value is fixed per class or copied from a Solar System body as a preset.

Conflict rules (Boss, 12:39 PM EDT, Decided): when in doubt verified wins; when two verified sources are in doubt the newest wins; when two sources are both in question, 'most permission take wins'. The third rule means the widest value or range of the two sources: Boss confirmed this at 1:12 PM EDT (Decided), as composition-model.md section 4a item 12 records it. Each place this note settles a disagreement carries 'Rule applied: verified / newest / most permissive', and section 15 lists them all.

## 1. Decisions

Decided (Boss, with the time as the coordinator relayed it, all EDT on 2026-10-17):

- **R1 Precedence.** Decided 3:02 AM EDT (CM 4a item 10): verified and partly verified values are usable unless contradicted; the verified source wins; contradicted and unsupported values do not ship. Decided 12:39 PM EDT: the three conflict rules above (verified, then newest, then most permissive).
- **R2 Downloaded data.** Decided between 11:19 AM and 11:23 AM EDT (GFP 6.5): line lists (HITRAN, ExoMol) for the correlated-k grids (AD3) are downloaded at build time, no derived table is stored in the repository, tests run offline from fixtures or a cache outside it. So no infrared line grid is in the repository, and the line part of the infrared opacity is open (R10).
- **R3 Spectra and ICRP 103.** Decided between 11:19 AM and 11:23 AM EDT (AD5, BC13): particle spectra come from the reachable NTRS reports (Band rows) and wR, wT from the secondary sources (08B section 4). **BC13 is closed from the free IAEA GSR Part 3 except the neutron weighting function**; ICRP 103 is not needed. Decided (Boss, contradiction sweep card, 7:29 AM EDT 2026-10-21). Boss's two ICRP uploads (a 35-page free extract and a strict subset of it) hold no wR, wT or dose limit; those tables are in ICRP 103 section 4 (printed pp. 61 to 78) and section 5.10 (from p. 98), which the extract does not reach, and IAEA GSR Part 3 gives the same values. The extract does give the ICRP 60 Q(L) function and the definitions (section 5). Hassler 2014 has been read (section 5).
- **R4 Targets, never presets.** Decided 1:43 AM EDT (CM 4a item 7): measured Solar System shapes enter only as flagged stand-ins with a replacing item; Earth, Mars, Venus, the Moon and Titan are test targets (section 13). The Moon and Mars dose anchors below 100 g/cm2 are stand-ins (R13).
- **R5 Seeded brackets.** Decided 1:43 AM EDT and 3:36 AM EDT (CM 4a items 6 and 11, stands as Decided unless Boss overrules): unsourced quantities are seeded draws in labelled brackets; haze yield is a nuisance parameter with the Brande scatter as prior (D5); life-linked gases are flagged design brackets, default off (D4), so this note adds no O2, only reads the stored one.

Default taken (drafted here, listed in section 16 for Boss):

- **R6 Band edges.** Ionising below 100 nm (plus X-rays and gamma); UV-C 100 to 280 nm; UV-B 280 to 315; UV-A 315 to 400; visible 400 to 700; near infrared 700 nm to 4 micron; thermal infrared beyond 4 micron. The UV names are partly verified (part 02 2.6, second-hand); the 700 nm edge is a convention. Evidence for 4 micron (computed, section 12): 0.98 percent of a 5772 K blackbody lies beyond it and 0.14 percent of a 288 K blackbody lies below it; at 3000 K 5.5 percent of starlight lies beyond it, which the model counts as thermal (stated limit).
- **R7 Transmission.** Beer-Lambert per wavelength from stored partial pressures, Rayleigh by the verified table or the generic polarizability formula with the per-gas error of Sneep and Ubachs (2.5 to 4 percent for N2 and Ar, 7 for CO2, 10 to 18 for CO, CH4, N2O; replaces the generic 5 to 9 percent of part 02; Rule applied: verified), Chapman slant path (never sec chi beyond 75 degrees), Eddington two-stream for thick air (verified, within 1.3 percent for mu0 of 0.5 and 1); O2 Schumann-Runge bands by k-distribution, not band means.
- **R12 Stored dose quantity.** Effective dose, ICRP 103 weights, for a reference person in the open air full time, all derived. The world-mean occupancy of UNSCEAR (0.2 outdoors, indoor gamma 1.4 times outdoor, indoor radon 40 Bq/m3) is rebuilt only in the Earth test (section 11), because a wild planet has no buildings.
- **R19 Ordering.** GEN.453 corrects the constants first; GEN.443 and GEN.444 replace the proxies together with a generator epoch bump and a migration of the dose columns; a test fails if a proxy returns (GEN.444).

Open (seeded bracket stated):

- **R8 Ozone.** Median from the Kozakis 2022 Sun-row fit (partly verified: 1-bar N2 runs, 0.01 to 150 percent PAL; it reproduces the paper's own values), host dependence by interpolating the five hosts in T_eff, PAL taken from p_O2 over 21.2 kPa (design inference). Spread open: the lane's plus or minus 30 percent does not reach Earth's 279 to 313 DU (fit 192 DU); seeded range times 0.7 to times 1.65 of the median (134 to 316 DU for Earth air; Rule applied: most permissive, the range widened to cover the observed 279 to 313 DU). Below 0.01 percent PAL: 0 to the clamp value (2.3 DU).
- **R9 Aerosols and cloud.** Cloud reflectance from the Monte Carlo table (part 06 3.10; partly verified: scalar, spherical, black floor). Column optical depth tau = 3 Q m_c / (4 rho_p r_eff); haze tau from the flux balance Mdot / v_f (engineering, not validated). Depths are drawn by GEN.414; bracket for cloud tau_c is the table's 0.5 to 100. Mars dust omega0: UV band 0.67 to 0.74, solar mean 0.89 to 0.94, red 0.93 to 0.975 (partly verified, 08B 1.1; Wolff 2009 gives figures only, no tables of omega0, n, k, so the ranges stay brackets); Titan haze extinction law confirmed (08B 2.2).
- **R10 Greenhouse opacity.** The grey tau0 of the Robinson and Catling closure is verified (my solver gives Earth 1.82 against the paper's 1.9, Venus 360 against 400). tau0 from composition is CIA (computed from HITRAN files) plus the water continuum (MT_CKD, radiation-term convention unverified) plus lines (correlated-k grids, R2). The line part has no seeded range in the research. Interim: the stored GEN.412 multiplier stays authoritative; the check reports tau0 implied and requires it above the CIA floor.
- **R11 Closure below tau0 of about 0.6.** The continuity root and the purely radiative closure give tau0 implied differing by a factor 2 (thin-air world: 0.25 against 0.12) for the same T_s; G is the same. Report both, bracket 0.12 to 0.25, until a lane settles it (Rule applied: most permissive).
- **R13 GCR curve.** D = D0(X_eff) f_B(M, X) f_sol. Verified tables for X of 100 to 1036 g/cm2 (Atri 2013); below 100 stand-in between the Moon anchor (500 mSv/yr) and the Mars point (232 mSv/yr, GCR only, Hassler 2014 Table 3, measured at about 21 g/cm2, taken at 22.6 g/cm2 from Guo's 840 Pa; verified) until a shower run (AD5) replaces it. **Built (GEN.746):** below 100 g/cm2 the curve is the primaries thinned by nuclear collisions (28.3 g/cm2, which the fit returns and which matches the CNO collision length) plus a Gaisser-Hillas-style build-up, fitted to the Moon, the Mars dose and its slope (-0.19 over 6.8 g/cm2, Guo 2018) and the shower peak of the first Atri row (553 mSv/yr at 100 g/cm2); stand-ins CR1 to CR4 in the gap-filling plan 6.22. Guo's 255 mSv/yr and the bracket 128 to 420 mSv/yr (0.35 to 1.15 mSv/day over the solar cycle) are kept as the stand-in's range (Rule applied: verified for the centre, most permissive for the range). The code's Mars point (16.4 g/cm2, 237 mSv/yr) matches the value but not the column; above 1036 hadronic cascade gone, muons (Crouch fit). f_sol 0.5 to 1.64 of the mean (Mars, spread 3.3). Hydrogen, helium, methane and ammonia columns are open (F8 is a design inference; bracket X_EM 0.58 X and X_had 1.73 X for H2).
- **R14 Stellar particles.** Integral rigidity Band form with parameters (verified; the published form wins over Boss's discontinuous one). R0 drawn in 0.019 to 0.364 GV (34 events, READ). Fluence law F30 = 1e9 cm-2 (E_flare / 1e31 erg) / d^2 (two READ points). Annual rate open (the flare frequency law is one star, GJ 1243, second-hand); fluence to dose open, interim primaries only (section 5). **Built (GEN.747):** the yearly dose is the sum of 47 measured events (Vainio 2017 Tables A.1 and A.2, Band fits of 1997 to 2015) at the blood-forming organs, 265 mSv/yr free in space at 1 AU for the Sun, less what an air column and the field keep out (`physics/solar_particles.py`); the NOAA 255-event list is the check; stand-ins SE1 to SE5 in 6.22.
- **R15 Ground abundances.** K, U, Th from GEN.459: Th from the star (offset sigma 0.22 dex, truncated plus or minus 0.5, or a catalogue value), U = Th times 0.571 (atoms; Dauphas 2005 via Wang 2020, partly verified) with an extra N(0, 0.05) dex placeholder, K by a depletion d_K log-uniform 0.03 to 1 on the CI K/U of 67407 (placeholder tied to the large K differences of Solar System bodies). The `a_rad` draw goes.
- **R16 Radon.** J = C_Ra f rho_b sqrt(lambda D_e), box concentration C = J tau / h with h of 1 km scaled by the scale height (design inference). D_e against pressure is not reached: seeded factor 1 to P0/P on D_e, so J between 1 and sqrt(P0/P) times the 1-bar value. The 1 kPa gate has no source and goes.
- **R17 Stellar UV.** The blackbody 200 to 300 nm proxy overpredicts the Sun by 3.9 times in UV-C, 1.6 in UV-B, 1.27 in UV-A (section 12). No source for UV against T_eff and activity was read; seeded UV-C share between 0.1 and 1 of the blackbody value (the Sun's measured ratio 0.26 is a target, not a preset).
- **R18 Ionising band and belts.** No surface dose from XUV or X-rays; the output is the deposition column for escape (part 04). W(air) is 33.97 eV (ICRU 90). Outer-shell cross sections for Ar, Kr, Xe, HCN, H2CO, HCl, HF, CH3Cl, SiO, TiO, VO are open. Trapped belts: AP9/AE9 exist for Earth only. **Built (GEN.749):** `physics/radiation_belts.py` gives the belts of any dipole from its moment, the wind-set magnetopause and the spin (a wind-limited outer edge at 0.75 of the standoff, a spin-limited one at 9 synchronous distances, an inner proton and an outer electron belt) and the dose rate behind 100 mil of aluminium at an orbit, scaled from Earth's by the moment to the 0.8 and checked on Earth and Jupiter; stand-ins are in the gap-filling plan 6.15 (RB1 to RB8). No belt dose is stored on the body; it is asked for an orbit.

## 2. The bands

| Band | What the stored mixture does to it | Source and state |
|---|---|---|
| Ionising (below 100 nm, X-ray, gamma) | Photoionisation reaches tau = 1 at 1e-6 to 1e-4 g/cm2 for every mixture, so none reaches the ground; X-ray tau = 1 at 15.7 g/cm2 at 1 MeV, 0.26 at 10 keV | part 07 2.2 to 2.4, verified (DATA, COMPUTED); gaps in R18 |
| UV-C, UV-B, UV-A | O2 (Schumann-Runge, Herzberg), O3 (Hartley), SO2, H2S, S8 absorb; CO2 and H2O stop below 200 and 220 nm; ppb of SO2 or O3 block 200 to 290 nm in a 1-bar column; Rayleigh tau at 300 nm is 12.6 times that at 550 nm | part 02 2.3, 2.7; part 04 2.5; verified (DATA) |
| Visible | Rayleigh, O3 Chappuis (peak 5.1e-21 cm2), cloud, haze, dust, ground albedo | part 02 2.1, 3.5; verified |
| Near infrared | stellar light absorbed by H2O, CO2, CH4 bands (shortwave heating); shortwave forcing of CO2 on N2 is -2, -4, -18 W/m2 at 0.01, 0.1, 1 | part 03 2.2 (B&G 2014b), verified |
| Thermal infrared | the planet's emission against the grey tau0 of CIA, lines and continuum | part 03 sections 3 and 5; verified/open (R10) |

## 3. Transmission by band from the stored mixture

Column and optical depths: X = 0.1 P / g in g/cm2; mean molar mass mu = sum x_i M_i; N_i = x_i X N_A / mu in cm-2. Rayleigh tau_R(lambda) = N sum x_i sigma_i, with sigma_i from part 02 table 2.1 (verified for 13 gases) or the generic 1.3057e-28 alpha^2 F_K lambda^-4 (1 + B/lambda^2)^2 (mean ratio 0.988, sd 0.046; 1.9 percent low on Earth air). Absorption tau_abs = sum N_j sigma_j(lambda, T) from part 02 2.3 and 2.4 (O3 column from GEN.413).

Direct beam F cos(chi) exp(-m tau), m the Chapman function of R/H. Thick air and the diffuse field use the Eddington transmission T(tau, mu0) = [(2/3 + mu0) + (2/3 - mu0) exp(-tau/mu0)] / (4/3 + tau), R = 1 - T, spherical albedo s = integral of 2 mu (1 - T), and with ground albedo A the reflectance R + T A (1 - s)/(1 - A s) and ground flux T/(1 - A s) (part 02 3.5). Band transmission is the star-weighted mean of T(lambda) over the band edges of R6.

Aerosols and cloud add tau_aer and (1 - g) tau_c to the two-stream, or use the Monte Carlo reflectance table (g 0.85, omega 1, mu0 0.5: 0.462 at tau 5, 0.870 at tau 50; omega 0.99 saturates near 0.60). Dust tau is a drawn clear-to-storm column; haze and cloud masses come from GEN.414.

The star's side: the Sun at 1 au gives 6.85, 17.44 and 86.00 W/m2 in UV-C (200 to 280), UV-B and UV-A, an erythemally weighted 666 and a Setlow-weighted 8.66 W/m2 (ATLAS-1, part 02 2.7). Other stars need a UV source (R17).

## 4. Greenhouse from the infrared opacity, and the check of the GEN.412 multiplier

The repo variable is G = (T_s/T_e)^4 - 1 with T_s^4 = (1 - A) S (1 + G) / (4 sigma). The closure of Robinson and Catling 2012 (part 03 5.3, 5.4, verified) has tau(p) = tau0 (p/p0)^n, diffusivity D = 1.66, lapse scaling alpha, beta = alpha (gamma - 1)/gamma, and T0^4 = T_e^4 (1 + D tau0 / 2) when purely radiative; otherwise tau_rc and T0 follow from the continuity of flux and temperature. The tropopause obeys p_tp proportional to tau0^(-1/n).

tau0 from composition: CIA by the Planck-weighted pair coefficients of part 03 3.1 (computable), water self-continuum from MT_CKD (exponent 5 to 6 in the window, the origin of the runaway), lines from the correlated-k grids of R2. The Koll equivalent grey thickness of CIA alone is a lower bound. Forcing laws (Byrne and Goldblatt 2014a, Etminan 2016) apply only as perturbations of an Earth-like state; Myhre 5.35 ln(C/C0) fails above 1 bar of CO2. Myhre 1998 Table 3 (read): CO2 5.35 and N2O 0.12 are Myhre's own values; IPCC 1990 gave 6.3 and 0.14, so the 'IPCC 1990' labels in part 03 should read Myhre 1998. Myhre's CO2 doubling is 3.71 W/m2; its Table 2 total (2.247) does not equal its rows (2.2537). Runaway limits: 282 (Simpson-Nakajima), 291 (Kopparapu) and 310 W/m2 (Kasting); a stored T_s above the limit with liquid water is rejected or made a steam atmosphere (part 03 6.3). A haze that absorbs sunlight and passes infrared cools the surface (Titan, over 20 K); that needs the two-channel form (Titan tau0 comes out 3.5 times low in one channel).

The check of a stored multiplier, in order: (1) from stored A, S, T_s compute T_e and G; (2) invert the closure for tau0 implied (monotone, bisection); (3) compute tau0 from composition; (4) hard test: tau0 implied is not below the CIA floor; soft test: ratio to the composition value within a factor set when the grids exist (open); (5) when the grids exist T_s is derived and G becomes the readout. The doubling test: a CO2 forcing of 5.32 ln 2 = 3.69 W/m2 (the printed `+ 0.39` term of part 03 2.2 is probably a squared-log term lost in extraction and is not used) needs tau0 raised by 5.0 percent and warms Earth's state by 1.12 K at fixed absorbed flux (section 11). Myhre's 5.35 ln 2 = 3.71 W/m2 gives 5.1 percent and 1.13 K; the two laws differ by 0.5 percent, so the test takes 3.7 W/m2 (Rule applied: newest, Byrne and Goldblatt 2014a, drives the computed shift; Myhre 1998 is the check).

## 5. Shielding by column and its composition

**Proton cutoff.** The kinetic energy that crosses X g/cm2 vertically comes from inverting the PSTAR range (part 07 2.5, verified): air 2241 MeV at 1033 g/cm2, 420 at 100, 171 at 22.5; CO2 165 MeV at 21 g/cm2 (247 MeV at 60 degrees, 10B 2.1); hydrogen needs 1.65 times, helium 1.06, water 1.07, rock 0.96 (ratios at 22.5 g/cm2). Mixtures use the Bragg sum of mass stopping powers. Slant paths and the nuclear loss exp(-X / lambda_I) (lambda_I 90.1 for air, 88.9 CO2, 52.0 H2, PDG) apply to primaries. Above about 100 g/cm2 the dose is cascade secondaries, not primaries.

**Galactic cosmic rays.** Dipole cutoff Rc = 14.34 GV (M/M_E)(R_E/r)^2 cos^4(lat) (scaling M/r^2, replaces (R_E/r)^3; partly verified: textbook form, constants computed). D0(X) is the Atri table at M = 0 (553, 263, 17.2, 3.02, 0.65 mSv/yr at 100, 200, 500, 700, 1036 g/cm2), f_B the Griessmeier ratios (M = 1: 0.31 at X = 100, 0.52 at 1036; valid X of 100 or more and M up to 10), proton primaries only, so a lower bound. Composition enters through X_eff = sum X_k (90.1 / lambda_k) (body part E 2.9d; ground gram = 0.89 air, ice 1.08). Neutron share: 2 to 7 percent of absorbed dose; secondary neutrons dominate only the dose equivalent (Mars about 50 percent at the surface, over 80 percent below 50 cm; absorbed share 2 to 9 percent; 08B 4.3, 10B 2.3).

**Stellar particles.** J(>R) = J0 R^-g1 exp(-R/R0) for R up to R1 = (g2 - g1) R0, and J0 R1^-g1 exp(-R1/R0) (R/R1)^-g2 above (continuous). The event of 2006-12-13 has J0 1.33e8, g1 1.045, g2 5.80, R0 0.1765 GV; the 1976 event has 7.62e6, 1.710, 6.33, 0.2054 (10B 1.2a). Transmitted share = the fluence above max(R_atm, Rc) over the >30 MeV fluence. Photons and gamma: mass attenuation by sum w_i (mu/rho)_i (part 07 2.4).

**Dose quantities.** Absorbed (Gy), equivalent (Sv, wR) and effective (Sv, wT). ICRP 103 wR (secondary source, R3): photons, electrons, muons 1; protons and charged pions 2; alpha particles, fission fragments and heavy ions 20; neutrons 2.5 + 18.2 exp(-[ln E]^2 / 6) below 1 MeV, 5.0 + 17.0 exp(-[ln 2E]^2 / 6) from 1 to 50, 2.5 + 3.25 exp(-[ln 0.04E]^2 / 6) above (peak 20.7 at 1 MeV, floor 2.5; Boss's 5 to 20 is corrected). wT (sum 1.00) belongs to tissues, not particles: marrow, colon, lung, stomach, breast, remainder 0.12; gonads 0.08; thyroid, oesophagus, liver, bladder 0.04; bone surface, skin, brain, salivary glands 0.01. The Moon anchor of 500 mSv/yr against 80 to 99 mGy/yr absorbed implies Q near 5.5. Mars, Hassler 2014 Tables 1 and 3 (read, GCR only, about 300 sols, column about 21 g/cm2): 0.21 +- 0.04 mGy/day, mean Q 3.05 +- 0.26, 0.64 +- 0.12 mSv/day, 76 mGy/yr and 232 mSv/yr; MSL cruise 0.48 mGy/day, Q 3.82, 1.84 mSv/day. Hassler uses the ICRP 60 function Q(L) (Q = 1 below L = 10 keV/micron, 0.32 L - 2.2 from 10 to 100, 300/sqrt(L) above 100; also in the ICRP extract), not ICRP 103 weighting factors, so the Mars and Moon anchors are dose equivalents, not effective doses, and the stored quantity (R12) differs from them by a definition. The paper's SEP sentence (about 25 percent of the GCR dose equivalent) holds for absorbed dose only (24 percent; 8 percent for dose equivalent); one SEP event is 0.025 mGy in Table 2 and about 50 microGy in the text (unexplained); the mission total is 982 mSv from the tables against about 1.01 Sv in the text. W(air) = 33.97 eV (ICRU 90) for ionisation rates.

**Muons.** Crouch fit I(h) = exp(-11.22 - 0.00262 h) + exp(-14.10 - 0.001213 h) + 2.17e-13 (h in hg/cm2 = m water equivalent), valid from about 1 km w.e.; dose = pi I x 2 MeV cm2/g. Below that it is an extrapolation (Titan's 10850 g/cm2 gives 3.4e-7 Gy/yr against 3.2e-6, a factor 9 low).

## 6. The ground: radiogenic dose and radon

Gamma dose in air at 1 m: D = 5.70 U + 2.44 Th + 13.1 K nGy/h (U and Th in ppm, K in wt percent; from UNSCEAR 0.462, 0.604, 0.0417 nGy/h per Bq/kg and 12.35, 4.04 Bq/kg per ppm, 314 Bq/kg per wt percent; verified, part 07 2.10). Effective annual dose = D x 8766 h x 0.7 Sv/Gy x occupancy. Age: scale each isotope by exp(-lambda (t - 4.567 Gyr)) with half-lives K-40 1.262e9, U-238 4.4683e9, Th-232 1.41e10 yr (part 07 2.9). Under dense gas multiply by exp(-0.0637 rho h) (rho in g/cm3, h = 100 cm); the Earth coefficients are valid for rho below about 0.03 g/cm3 (F9).

Radon: J = 12.35 U f rho_b sqrt(lambda D_e) (f 0.2, rho_b 1600 kg/m3 recovered from the 0.026 of part 07 2.10, D_e 2e-6 m2/s), C = J tau / h with tau 5.52 d, dose 9 nSv per Bq h m-3 times equilibrium factor 0.6 outdoors. Heat uses the same abundances: h = rho (3.387e-3 K + 26.18 Th + 98.29 U) microW/m3 (kg/kg), so the dose and the radiogenic flux of GEN.459 cannot drift apart (primitive mantle with the 2003 edition values below: 21.0 TW for 4.042e24 kg; the older 20 ppb U, 77 ppb Th, 280 ppm K of part 07 2.9 gave 19.95 TW with Rb; Rule applied: verified). Reference abundances: Palme and O'Neill **2003 edition** Table 4, K 260 ppm, Th 83.4 ppb, U 21.8 ppb (read twice); the 'Palme and O'Neill 2014' Th 84.9 and U 22.9 ppb of our files (1.8 and 5.0 percent higher) cannot be confirmed from a paper and are not used (Rule applied: verified; newest cannot apply because the 2014 edition is not in hand). The primitive-mantle oxides sum to 98.41 weight percent, not 100 (1.34 minor oxides, under 0.15 trace), so a normalisation to 100 must say so. Crust: Rudnick and Gao 2014 Table 11 (the 2003 text), upper crust Th 10.5, U 2.7 ppm, K2O 2.80 weight percent (K 23,244 ppm); total crust 5.6, 1.3 and 15,026 ppm; every crust value in our files matches it, and Palme and O'Neill's own crust column (Th 3.5, U 0.91 ppm, K 9100 ppm) disagrees (Rule applied: verified, the dedicated crustal review). A rocky surface is crust, not mantle: the primitive mantle gives 0.67 nGy/h, the upper crust 71.5 and the UNSCEAR soil 59.9, so worlds B, C, D carry the primitive mantle as the low bracket and B also the upper-crust value. Ground type and band albedo come from GEN.441 (lab reflectance, not Bond albedo; the conversion is open). An airless ice moon has almost no K, U, Th and so almost no ground dose.

## 7. Surface ultraviolet from transmission and reflectance

E(lambda) at the ground = direct + diffuse + ground-reflected, from section 3. UVI = 40 m2/W x integral S_er E dlambda, with the CIE erythemal S_er = 1 for 250 to 298 nm, 10^(0.094 (298 - lambda)) to 328 nm, 10^(0.015 (140 - lambda)) to 400 nm (partly verified, the NOAA text has a typo). DNA dose = integral S_Setlow (lambda/300) E dlambda (energy basis; the file says unity at 300 nm but stores 0.0333, reported not resolved). The direct UVI falls as DU^-1.08 for 100 to 500 DU. A planet with O2 but no ozone has lethal UV-C (Setlow dose about 280 times the 300 DU case). Ground reflectance enters through 1/(1 - A s), an upper bound because the return path is not ozone-filtered.

## 8. Derived against stored

| Quantity | Stored | Derived (pure function, no new column) |
|---|---|---|
| Pressure, gravity, partial pressures (GEN.411) | yes | X, mean molar mass, mole fractions |
| Albedo, greenhouse multiplier, T_s (GEN.412) | yes | T_e, G, tau0 implied, tau0 from composition, the check |
| Ozone column and p_O3 (GEN.413) | yes | band absorption, UVI, DNA dose |
| Cloud, dust, haze depths (GEN.414) | yes | reflectance, band transmission |
| Rayleigh tau by band (GEN.417) | no | from the table |
| Band fluxes at the ground (GEN.443) | no | from the star, section 3 |
| K, U, Th at formation, age (GEN.459) | yes | present activities, gamma, radon, heat |
| Ground type, grain, band albedos (GEN.441) | yes | ground reflectance |
| Magnetic moment, star T_eff, L, flare rate | yes | Rc, share, event fluence |
| `dose_gcr_msv_yr`, `dose_sep_msv_yr`, `dose_ground_msv_yr`, `surface_dose_msv_yr`, `uv_surface_index`, `ozone_loss_flag` | columns kept | values become derived (epoch bump, migration); `dose_helio_mult` is stored 1.0 and has no source |

## 9. Proxies replaced (radiation.py, habitability.py, planets.py as read)

| Proxy today | Replaced by | Item |
|---|---|---|
| `ozone_ratio` = sqrt(p_O2 / 21.2 kPa), capped 1.2, 0 under 0.02 kPa | Kozakis column by O2 and host (R8); the fit's own floor is 0.01 percent PAL, 2.1 Pa | GEN.413, used by GEN.443 |
| `uv_index` = `uv_ratio` x flux x ozone^-1.6 (a ratio to the Sun, not the CIE index) | UVI and DNA dose of section 7 | GEN.445. **Built (epoch 74):** `uv_surface_index` = `body_radiation.dna_dose_index(body)`, the Setlow dose of the ground spectrum over Earth's; the proxy and its tests are gone |
| a loss event kept the stored ozone column as it was (`ozone_loss_flag` only said so) | (none: a helper `radiation.apply_air_loss` scaled the column by the share of the air kept and worked the UV index again) | GEN.668. **Removed by TEST.524:** nothing called it (no generation step removes air after the ozone is drawn: escape and the water loss act on the air before it, and the final-age refresh, GEN.839, works the dose again), so no generated value changed |
| `_band_fraction`, `uv_ratio` (blackbody 200 to 300 nm, floor 0.1) | band shares of R6 and the stellar UV of R17 | GEN.443. **Built:** removed from `radiation.py` with GEN.445 |
| `a_rad` log-uniform 0.5 to 2 times 0.48 mSv/yr times heat ratio, plus radon 0.5 mSv/yr behind a 1 kPa gate | gamma and radon from K, U, Th (section 6) | GEN.444 with GEN.459 |
| `GCR_POINTS` (500, 237, 110), `GCR_ATTENUATION_G_CM2` 164, `magnetic_factor` (cut 0.4, flat) | Atri and Griessmeier tables (R13) | GEN.453 then GEN.444. **Built (GEN.453):** `GCR_SHOWER_POINTS`, `GCR_ANCHORS`, `magnetic_factor(cutoff, column)` |
| `cutoff_rigidity_gv` (R_E/R)^3 | M / r^2 | GEN.453. **Built:** 14.34 GV for Earth |
| `proton_cutoff_gv` range law `R = 0.0022 E^1.77` | NIST PSTAR range tables | GEN.453. **Built:** `physics/proton_range.py` (mixtures by Bragg's rule for GEN.444) |
| `helio_multiplier` 1 + 1.5 exp(-X/100) c (no source) | dropped or sourced | GEN.453, GEN.256 |
| `sep_dose_msv_yr`, `free_sep_msv_yr` (20 mSv/yr, R0 0.2 GV, range law 0.0022 E^1.77) | Band spectrum, PSTAR inversion, fluence law (R14) | GEN.444. **Built (GEN.747):** `solar_particles.FREE_MSV_YR` and `dose_share` |
| habitability.py second GCR curve (0.45 exp(-X/22) + 0.05 exp(-X/200) + 0.002 Sv/yr) | reads the derived dose | GEN.445. **Built (epoch 78):** `habitability.bare_cosmic_ray_dose_sv_yr` reads `radiation.no_field_dose_msv_yr`; the second curve is gone |
| the Radiation row said a dose and nothing of what made it | the terms and the composition behind each | GEN.445. **Built:** `physics/habitability_explain.radiation_terms(body)` gives cosmic rays, stellar particles, the ground and the ultraviolet with the air column and gases, the cutoff, the stored potassium, thorium and uranium and the ozone, and which of the dose and the ultraviolet sets the rating; the planet's details carry them as a Radiation table (`describe/build._radiation_block`) and `/help/habitability` says what sets each (`habitability_ref.radiation_text`); tests move a gas, the ozone, the air column, the field and the ground's potassium and check the terms and the rating follow |
| planets.py per-class `greenhouse_multiplier_range` times the molar-mass ratio to 0.04345, ceiling 500 | mixture-derived tau0 and the check (R10) | GEN.423 with GEN.443 |

Descriptions of the dose, the UV and the transmission read only these derived values (GEN.446). `magnetism.py` (dipole share, standoff, XUV, wind pressure) stays and feeds the field and the ionising band.

## 10. Boss-supplied radiation research and how the lanes checked it

| Claim (Boss's note) | Check | State |
|---|---|---|
| ICRP 103 wR: photons and electrons 1, protons 2, alpha and heavy 20, neutrons 5 to 20 (1:22 AM EDT) | values confirmed; neutrons are a continuous function, 2.5 to 20.7; wT belongs to tissues, not particles (body lane F item 7, 1:28 AM EDT; 08B 4, about 2:00 AM EDT) | partly verified (secondary source); BC13 closed from IAEA GSR Part 3 except the neutron weighting function (R3) |
| Secondary neutrons dominate subsurface dose | holds for dose equivalent only: Mars about 50 percent at the surface, over 80 percent below 50 cm; absorbed 2 to 9 percent; MSL RAD 2.2 to 7 percent (08B 4.3, 10B 2.3) | partly verified |
| AP9, AE9, SPM (Ginet 2013) for trapped belts | the model exists (Crossref) and is for Earth; other bodies not reached (08A) | open |
| W values and EUV outer-shell cross sections | W(air) 33.97 eV (ICRU 90) replaces 33.85; Kr, Xe, HF, CH3Cl, TiO, VO cross sections not reached (10B 5) | partly verified / open |
| Band function for solar particles (2:11 AM EDT) | discontinuous as written (CONTRADICTED); the published integral rigidity form with parameters passes: F30/F200 = 28.0 against 26.1 plus or minus 6.8 (10B 1.2a) | verified |
| Mars 0.21 mGy/day, 0.64 mSv/day; 150 MeV cutoff at 20 to 22 g/cm2 | doses read in Hassler 2014 (0.21 +- 0.04 mGy/day, Q 3.05, 0.64 +- 0.12 mSv/day, 76 mGy/yr, 232 mSv/yr, column about 21 g/cm2, ICRP 60 Q(L)); solar-cycle factor 3.3; cutoff is 165 MeV vertical, 247 at 60 degrees (10B 2) | doses verified; cutoff partly verified |
| Europa 1e3 to 1e4 Gy/yr and 5.4 Sv/day | not a surface value; true at 2 to 4 cm depth in ice; 5.4 Sv/day has no primary source (10B 3) | contradicted |
| Titan column 1000 g/cm2, dose 1.5 mGy/yr | column is 10850 g/cm2; dose about 3e-6 Gy/yr (10B 4) | contradicted |
| K/U and Th/U by reservoir, Venus K/U | Venus K/U is 6 to 8 thousand; Solar System values are targets only (body lane section 15, G2) | partly verified |

Counts for the second note: 148 sub-claims, 39 confirmed, 57 partly, 30 contradicted, 15 unsupported, 7 not reached (master 3c).

## 11. Worked example

Four test draws, each sized to a target (Earth, Mars, Venus, the Moon) and never a preset. Class mixtures are the templates of `atmosphere.py` (M, K, N); the Sun at 1 au is S = 1361.17 W/m2 (computed from 3.828e26 W). Inputs: **A** Earth-like: 101.325 kPa, g 9.80665, 1.0 au, Bond albedo 0.30, T_s 288 K, dipole 1 Earth, soil of the UNSCEAR means (33, 45, 420 Bq/kg). **B** thin air: class K mixture at 0.8 kPa, g 3.71, 1.524 au, albedo 0.25 (recalled, part 03), T_s 215 K (recalled), no field, primitive-mantle ground (Palme and O'Neill 2003 edition: U 0.0218, Th 0.0834 ppm, K 260 ppm). **C** thick CO2: class N mixture at 9.2 MPa, g 8.87, 0.723 au, cloud tau_c 30 with omega 0.999 giving albedo 0.752 (table interpolation), T_s 737 K, same ground. **D** airless: 1.0 au, albedo 0.12, same ground, no field.

Steps per world: X, mu, tau_R; T_e from A and S; G and tau0 implied; the stored class-mid multiplier of `tuning.py` (molar mass over 0.04345 times the middle of the multiplier range) and the tau0 it implies; ozone, UV; cutoff and event; GCR; gamma and radon; then the proxies as `radiation.py` computes them. Piece checks pass first (section 12 output): dry air tau_R 0.0970 against Bucholtz Table 4 0.09721 (-0.2 percent; part 02 computed 0.0969); Kozakis 191.8, 199.4, 108.0 DU; Earth tau0 1.819, tau_rc 0.0863, G 0.6254; Venus tau0 359.9, G 99.6; PSTAR 2237 MeV (table 2241) and 165; Band ratio 28.0; muon 3.8e-5 at 1 km w.e. (4e-5); soil gamma 59.9 nGy/h (60); GCR at M = 1 0.341 mSv/yr at sea level (0.32 measured) and 45.2 at 265 g/cm2 (47.9 aviation).

| Output | A Earth-like | B thin air | C thick CO2 | D airless |
|---|---|---|---|---|
| X (g/cm2); mean molar mass | 1033; 28.87 | 21.56; 43.45 | 1.037e5; 43.45 | 0 |
| tau_R at 550 and 300 nm | 0.0970; 1.23 | 0.0034; 0.043 | 16.4; 207 | 0 |
| T_e (K) | 254.6 | 209.8 | 231.0 | 269.6 |
| G implied | 0.638 | 0.103 | 102.6 | none |
| tau0 implied (tau_rc) | 1.86 (0.087) | 0.25 root, 0.12 radiative | 374 (1.25) | none |
| CIA floor tau_eq (part 03 3.1) | 0.040 | 1.3e-4 (scaled from 8.5e-5 at 6.4 mbar) | 0.056 (file range limited) | none |
| Stored class-mid G, T_s it gives | 1.168, 308.9 K | 0.034, 211.6 K | 277, 943 K | none |
| tau0 that stored G implies | 4.3 | 0.083 | 1426 | none |
| Ozone (DU) | 192 (134 to 316) | 0 to 2.3 | 0 | 0 |
| UVI overhead Sun | direct 12.0, global 23.9 (300 DU: 7.27, 13.9) | 271 | 0 (SO2 tau 9e4) | 666 |
| Setlow dose (W/m2) | 0.0105 | 3.44 | 0 | 8.66 |
| Proton cutoff (MeV) | 2237 | 167 | 2.1e5 to 2.3e5 | none |
| Event 1e32 erg, primaries only (mSv) | 5e-7 | 448 | 0 | 3.9e4 |
| GCR (mSv/yr) | 0.341 | 232 stand-in (Hassler; range 128 to 420) | 3.5e-5 (muons) | 500 stand-in |
| Gamma, full time (mSv/yr) | 0.368 | 0.0041 (upper crust 0.44) | 0.0027 | 0.0041 |
| Radon (mSv/yr) | 0.488 | 0.045 (x1 to x11 on J; upper crust 5.6) | 0.0004 | none |
| Total derived (mSv/yr) | 1.20 | 232 | 0.0032 | 500 |
| GCR as `radiation.py` gives | 0.223 | 211.6 | 3.7e-4 | 500 |
| Ground term today (`a_rad` 0.5 to 2) | 0.49 to 1.96 | 0.24 to 0.96 | 0.49 to 1.96 | 0.24 to 0.96 |
| UV index today (ratio to Earth) | 1.0 | 431 | 1913 | 1000 |

Readings. A: the closure reproduces Earth's 288 K at tau0 1.86; plus or minus 0.1 in tau0 moves T_s by plus or minus 1.2 K. The world-mean rebuild is GCR 0.34 + external 0.485 + radon 1.108 + ingestion 0.29 (READ) = 2.22 mSv/yr against UNSCEAR 2.4 (7 percent low; Rn-220 0.1 is omitted). The Kozakis median makes the overhead UVI 12.0 direct, 1.7 times what 300 DU gives; the ozone spread therefore sets the UV. B: no ozone, so the Setlow dose is 356 times A's and UV-C reaches the ground. The stored class multiplier is 3 times below the G the inverse needs for a 215 K state, and its tau0 is below the recalled 0.12 to 0.25. C: no UV and no primaries reach the ground; the class multiplier of 277 is 2.7 times the G that Venus-like inputs imply, because the class tuning also absorbs the orbital distance (tuning.py comment on class N); `uv_index` gives the highest value to the world with zero UV. The gamma attenuation factor 0.66 is outside the 0.03 g/cm3 validity of F9 (gas is 0.065). D: no shielding, the Moon anchor of 500 is a stand-in; a CO2-doubling check on A (tau0 +5.0 percent, +1.12 K) is the sensitivity target for the line grids.

## 12. Code

Run as `python3 -I`; all constants are those cited in the comments.

```python
import math, numpy as np
from scipy.integrate import quad, trapezoid
from scipy.optimize import brentq
from scipy.interpolate import CubicSpline
NA, SIG, D = 6.02214076e23, 5.670374419e-8, 1.66
S0 = 3.828e26/(4*math.pi*149597870700.0**2)                       # W m-2 at 1 au
MM = dict(n2=28.014,o2=31.998,ar=39.95,h2o=18.015,co2=44.009,ch4=16.043,so2=64.058,h2s=34.076,co=28.010)   # part 01
SR = dict(n2=4.635e-27,o2=4.072e-27,ar=3.991e-27,co2=1.164e-26,ch4=9.905e-27,h2o=3.38e-27,so2=2.41e-26,h2s=1.98e-26,co=5.734e-27)  # part 02 2.1, cm2 at 550 nm
BSE = (0.0218, 0.0834, 0.026)                  # primitive mantle, U ppm, Th ppm, K wt%: Palme and O'Neill, 2003 edition, Table 4 (21.8 ppb, 83.4 ppb, 260 ppm)
UCC = (2.7, 10.5, 2.80*.8301)                  # upper crust: Rudnick and Gao Table 11 (the 2003 text), U ppm, Th ppm, K2O 2.80 wt% x 0.8301
SOIL = (33/12.35, 45/4.043, 420/314.)          # UNSCEAR soil Bq/kg over Bq/kg per ppm (part 07 2.10)
STAND = dict(B=232., D=500.)                   # GCR stand-ins, mSv/yr (R13): Hassler 2014 Table 3 at 22.6 g/cm2 for B, the Moon for D
W = dict(  # d au; alb Bond; Ts K (test input); M dipole in Earth units; gr ground; cls = class-mid (molar mass kg/mol, multiplier) of tuning.py; rc = (gamma, alpha, n)
 A=dict(mix=dict(n2=.7703,o2=.21,ar=9.3e-3,h2o=.01,co2=4e-4,ch4=1.9e-6),P=101325.,g=9.80665,d=1.,alb=.30,Ts=288.,M=1.,gr=SOIL,cls=(.0290,1.75),rc=(1.4,.6,2)),
 B=dict(mix=dict(co2=.95,n2=.028,ar=.019,o2=1.5e-3,co=7e-4,h2o=3e-4),P=800.,g=3.71,d=1.524,alb=.25,Ts=215.,M=0.,gr=BSE,cls=(.04265,.035),rc=(1.3,.8,2)),
 C=dict(mix=dict(co2=.965,n2=.035,so2=1.5e-4,h2s=1e-4,ar=7e-5,co=1.7e-5,h2o=2e-5),P=9.2e6,g=8.87,d=.723,alb=.752,Ts=737.,M=0.,gr=BSE,cls=(.0434,277.5),rc=(1.3,.8,1)),
 D=dict(mix={},P=0.,g=1.62,d=1.,alb=.12,Ts=None,M=0.,gr=BSE,cls=None,rc=None))
def column(w):
    if w['P'] <= 0: return dict(X=0., mu=0., N=0., tau=0., x={})
    s = sum(w['mix'].values()); x = {k: v/s for k, v in w['mix'].items()}
    mu = sum(x[k]*MM[k] for k in x); X = w['P']/w['g']*0.1; N = X*NA/mu                  # g/mol, g/cm2, cm-2
    return dict(X=X, mu=mu, N=N, tau=sum(x[k]*SR[k] for k in x)*N, x=x)
ray = lambda nm: (550/nm)**4*((1+.0075/(nm/1e3)**2)/(1+.0075/.55**2))**2          # tau_R(nm)/tau_R(550), part 02 3.1
KS = [6.203e16,7.716e17,3.147e18,4.157e18,8.412e17,4.640e18]                     # Kozakis 2022, Sun row
ozone_du = lambda po2: float(np.polyval(KS, math.log10(0.21*min(max(po2/21.2e3, 1e-4), 1.5)))/2.687e16)
DUg = np.log([100,200,300,500.]); UD, UG, DN = (np.log(a) for a in ([23.4,11.5,7.27,4.11],[48.9,22.8,13.9,7.4],[.0339,9.7e-3,4.4e-3,1.6e-3]))
uv = lambda du, t: float(np.exp(np.interp(math.log(du), DUg, t)))               # part 02 2.7, 100 to 500 DU only
Tedd = lambda tau, mu: ((2/3+mu)+(2/3-mu)*math.exp(-tau/mu))/(4/3+tau)
sph = lambda tau: quad(lambda m: 2*m*(1-Tedd(tau, m)), 0, 1)[0]
def bb(T, l1, l2):                                                               # blackbody share of sigma T^4 between l1 and l2 (m)
    f = lambda l: 2*6.62607015e-34*299792458.**2/l**5/(math.exp(min(700, 6.62607015e-34*299792458./(l*1.380649e-23*T)))-1)
    return quad(f, l1, l2, limit=200)[0]*math.pi/(SIG*T**4)
def rc(tau0, F, gm, al, n):                    # Robinson and Catling 2012 (part 03 5.4): T0 and tau_rc
    m = 4*al*(gm-1)/gm/n
    def f(t):
        I = quad(lambda s: (s/tau0)**m*math.exp(-D*(s-t)), t, tau0)[0]
        return (math.exp(-D*(tau0-t))+D*I)/(t/tau0)**m-(2+D*t)/(1+D*t)
    xs = np.logspace(-6, math.log10(tau0), 200); v = [f(x) for x in xs]
    for a, b, fa, fb in zip(xs, xs[1:], v, v[1:]):
        if fa*fb < 0:
            t = brentq(f, a, b); return ((F*(1+D*t)/2)/SIG/(t/tau0)**m)**.25, t
    return ((F/SIG)*(1+D*tau0/2))**.25, None
tau0_for = lambda Ts, F, *p: brentq(lambda t: rc(t, F, *p)[0]-Ts, 1e-6, 1e6, xtol=1e-9)
PS = dict(air=[.1408,1.007,8.744,58.16,367.2,1454,5054], co2=[.1404,1.004,8.722,58.02,366.3,1450,5044])  # PSTAR CSDA range g/cm2, part 07 2.5
EN = np.array([10,30,100,300,1000,3000,10000.])
fR = lambda g: CubicSpline(np.log(EN), np.log(PS[g])); iR = lambda g: CubicSpline(np.log(PS[g]), np.log(EN))
Tmin = lambda X, g: float(math.exp(iR(g)(math.log(X))))                          # proton MeV that crosses X g/cm2
rig = lambda T: math.sqrt(T*T+2*T*938.272)/1000                                  # GV
J0, G1, G2, R0 = 1.33e8, 1.045, 5.80, 0.1765; R1 = (G2-G1)*R0                    # GLE 2006-12-13
J = lambda R: J0*R**-G1*math.exp(-R/R0) if R <= R1 else J0*R1**-G1*math.exp(-R1/R0)*(R/R1)**-G2
def share(Rt, Rc0):                            # area-weighted share of the >30 MeV fluence above max(Rt, Rc0 cos^4 lat)
    lat = np.linspace(0, math.pi/2, 721); w = np.cos(lat)
    return trapezoid(np.array([J(max(Rt, Rc0*math.cos(l)**4)) for l in lat])*w, lat)/trapezoid(w, lat)/J(rig(30))
def Sw(T):                                     # Bethe stopping power in water, MeV cm2/g (part 07 F2, I = 75 eV)
    me = .511; g = 1+T/938.272; b2 = 1-1/g**2; tm = 2*me*b2*g*g/(1+2*g*me/938.272+(me/938.272)**2)
    return .307075*.55509/b2*(.5*math.log(2*me*1e6*b2*g*g*tm*1e6/75.**2)-b2)
def sep_event(X, gas, F30, lam=88.9):          # primaries only: fluence (cm-2) and dose (mSv, wR 2) in 1 g/cm2 of water
    T = np.logspace(math.log10(30), 4.3, 600); dN = -np.diff([J(rig(t)) for t in T])/J(rig(30))*F30
    Tm = np.sqrt(T[:-1]*T[1:]); surv = math.exp(-X/lam)
    if X > 0:
        r = np.exp(fR(gas)(np.log(Tm)))-X; rr = np.maximum(r, PS[gas][0])
        res = np.where(r > PS[gas][0], np.exp(iR(gas)(np.log(rr))), 10*(np.maximum(r, 0)/PS[gas][0])**(1/1.79))
    else: res = Tm
    dep = np.array([min(e, Sw(e)) if e > 0 else 0. for e in res])
    return float(np.sum(dN[res > 0])*surv), float(np.sum(dN*surv*dep)*1.602e-10*2e3)
XD = np.array([100,200,500,700,1036.]); D0 = np.array([553,263,17.2,3.02,.65])    # Atri 2013, M = 0 (mSv/yr)
MG = np.array([0,.1,.15,.25,.5,.75,1,2,3,6,10.])                                  # ratios to M = 0
Q100 = np.array([1,.95,.92,.73,.47,.39,.31,.096,.027,.010,.004]); Q1036 = np.array([1,.74,.71,.68,.65,.60,.52,.43,.35,.29,.15])
def gcr(X, M):                                 # mSv/yr, X from 100 to 1036 g/cm2
    t = (math.log(X)-math.log(100))/(math.log(1036)-math.log(100)); a, b = np.interp(M, MG, Q100), np.interp(M, MG, Q1036)
    return float(math.exp(np.interp(math.log(X), np.log(XD), np.log(D0)))*(a+(b-a)*t))
def crouch(X):                                 # muons, mSv/yr (wR 1); the fit holds from 1 km w.e. = 1e5 g/cm2
    h = X/100; I = math.exp(-11.22-.00262*h)+math.exp(-14.10-.001213*h)+2.17e-13
    return math.pi*I*2.0*1.602e-10*3.156e7*1e3
gam = lambda U, Th, K: 5.70*U+2.44*Th+13.1*K                                     # nGy/h in air at 1 m, part 07 2.10
eff = lambda g, occ: g*8766e-6*0.7*occ                                            # mSv/yr
LAM = math.log(2)/(3.824*86400)
heat = lambda U, Th, K, M=4.042e24: M*(9.829e-5*U*1e-6+2.618e-5*Th*1e-6+3.387e-9*K*1e-2)/1e12   # TW, part 07 2.9 law
def radon(U, fD=1.):                           # flux Bq m-2 s-1, box concentration Bq m-3 (h 1 km), outdoor full-time mSv/yr
    j = 12.35*U*.2*1600.*math.sqrt(LAM*2e-6*fD); c = j*5.52*86400/1000.
    return j, c, 9e-6*.6*c*8766
def today(X, cut, po2, d):                     # radiation.py as read: gcr*magnetic_factor, ozone_ratio, uv_index (Sun)
    pts = ((0., 500.), (16.4, 237.), (100., 110.))
    if X >= 100: g = max(110*math.exp(-(X-100)/164.), 1e-3*110*math.exp(-933/164.))
    else:
        (x0, d0), (x1, d1) = (pts[0], pts[1]) if X <= 16.4 else (pts[1], pts[2])
        g = math.exp(math.log(d0)+(math.log1p(X)-math.log1p(x0))/(math.log1p(x1)-math.log1p(x0))*(math.log(d1)-math.log(d0)))
    g *= 1 if cut <= 1 else 1-.4*min(1, math.log10(cut)/math.log10(15))
    o = 0. if po2 < 20. else min(math.sqrt(po2/21.2e3), 1.2)                     # 0.02 kPa threshold
    return g, (1e3 if o <= 0 else min(1e3, o**-1.6))/d**2
# ---- pieces against their sources
dry = column(dict(mix=dict(n2=.78084,o2=.20946,ar=.00934,co2=4e-4),P=101325.,g=9.80665))
print("Rayleigh dry air %.4f (Bucholtz Table 4 0.09721; part 02 0.0969); Kozakis %.1f %.1f %.1f DU (191.8 199.4 108.0)" % ((dry['tau'],)+tuple(ozone_du(21.2e3*p) for p in (1,.1,.01))))
t0 = tau0_for(288., 240., 1.4, .6, 2)
print("R&C Earth tau0 %.3f tau_rc %.4f G %.4f; Venus tau0 %.1f G %.1f" % (t0, rc(t0, 240., 1.4, .6, 2)[1], (288/(240/SIG)**.25)**4-1, tau0_for(730., 160., 1.3, .8, 1), (730/(160/SIG)**.25)**4-1))
print("PSTAR air 1033 %.0f MeV (2241), CO2 21 %.0f (165); Band F30/F200 %.1f (28.0), share above 165 MeV %.3f (0.057); muon 1 km w.e. %.1e (4e-5); soil gamma %.1f nGy/h (60); GCR M=1 sea level %.3f (0.32), 265 g/cm2 %.1f (47.9)" % (Tmin(1033.23,'air'), Tmin(21,'co2'), J(rig(30))/J(rig(200)), J(rig(165))/J(rig(30)), crouch(1e5), gam(*SOIL), gcr(1033.23,1.), gcr(265.,1.)))
print("heat TW: primitive mantle 2003 edition %.1f; U 20 ppb, Th 75.4 ppb, K 280 ppm %.1f (part 07 19.95 with Rb); gamma nGy/h: primitive mantle %.2f, upper crust %.1f, soil %.1f" % (heat(*BSE), heat(.020, .0754, .028), gam(*BSE), gam(*UCC), gam(*SOIL)))
print("Sun blackbody W/m2 UV-C %.1f UV-B %.1f UV-A %.1f (ATLAS-1: 6.85 17.44 86.00); beyond 4 um: 5772 K %.4f, 3000 K %.4f; 288 K below 4 um %.4f" % (*(S0*bb(5772, a*1e-9, b*1e-9) for a, b in ((200,280),(280,315),(315,400))), bb(5772,4e-6,1e-3), bb(3000,4e-6,1e-3), bb(288,1e-9,4e-6)))
# ---- the four worlds
for k, w in W.items():
    c = column(w); F = (1-w['alb'])*S0/w['d']**2/4; Te = (F/SIG)**.25; po2 = w['mix'].get('o2', 0)/sum(w['mix'].values() or [1])*w['P']
    print(f"== {k}: X {c['X']:.4g} g/cm2, mu {c['mu']:.2f}, tau_R550 {c['tau']:.4g}, tau_R300 {c['tau']*ray(300):.4g}, F {F:.1f} W/m2, T_e {Te:.1f} K, pO2 {po2:.4g} Pa")
    if w['Ts']:
        G = (w['Ts']/Te)**4-1; t0 = tau0_for(w['Ts'], F, *w['rc']); Gs = w['cls'][0]/.04345*w['cls'][1]; Ts_s = Te*(1+Gs)**.25
        print(f"   greenhouse: G implied {G:.4g}, tau0 {t0:.4g}, tau_rc {rc(t0, F, *w['rc'])[1]:.3g}, radiative closure 2G/D {2*G/D:.4g}; stored class-mid G {Gs:.4g} -> T_s {Ts_s:.1f} K, tau0 it implies {tau0_for(Ts_s, F, *w['rc']):.4g}")
    X = c['X']; gas = 'air' if k == 'A' else 'co2'; U, Th, K = w['gr']; g0 = gam(U, Th, K); att = 1.; fD = None; gc = STAND.get(k); dr = 0.
    if k == 'A':
        du = ozone_du(po2); t3 = c['tau']*ray(300); gc = gcr(X, w['M']); fD = 1.
        print(f"   ozone {du:.1f} DU (x0.7 to x1.65: {0.7*du:.0f} to {1.65*du:.0f}); UVI direct {uv(du, UD):.1f}, crude global {uv(du, UG):.1f}; at 300 DU {uv(300, UD):.2f}, {uv(300, UG):.1f}; Setlow direct {uv(du, DN):.4f} W/m2; s(300 nm) {sph(t3):.3f}, ground factor A 0.05 {1/(1-.05*sph(t3)):.3f}, A 0.8 {1/(1-.8*sph(t3)):.3f}")
    elif k == 'B': fD = 101325/w['P']; print(f"   ozone 0 to {ozone_du(po2):.1f} DU (clamped, upper bound); UVI {666/w['d']**2*math.exp(-c['tau']*ray(280)):.0f}; Setlow {8.66/w['d']**2*math.exp(-c['tau']*ray(260)):.2f} W/m2")
    elif k == 'C': fD = 101325/w['P']; gc = crouch(X); rho = w['P']*c['mu']/(8.314*w['Ts'])/1e6; att = math.exp(-.06369*rho*100); print(f"   UVI 0: SO2 optical depth at 250 to 290 nm {c['x']['so2']*c['N']*4.3e-19:.3g}, tau_R300 {c['tau']*ray(300):.0f}; gas density {rho:.4f} g/cm3")
    else: print(f"   UVI {666/w['d']**2:.0f}; Setlow {8.66/w['d']**2:.2f} W/m2 (unfiltered)")
    if k == 'C': print(f"   shielding: cutoff {10+(X-5044)*2.06/1000:.0f} to {10+(X-5044)*2.27/1000:.0f} GeV, X/lambda_I {X/88.9:.0f}; muons (Crouch, edge of validity) {gc:.2g} mSv/yr")
    else:
        Tc = Tmin(X, gas) if X > 0 else 30.; phi, dose = sep_event(X, gas, 1e10/w['d']**2)
        print(f"   shielding: proton cutoff {Tc:.0f} MeV ({rig(Tc):.2f} GV), share of >30 MeV fluence {share(rig(Tc), 14.34*w['M']):.3g}; 1e32 erg event fluence {phi:.3g} cm-2, primaries-only dose {dose:.3g} mSv")
    gm = eff(g0, 1)*att; print(f"   ground: gamma {g0:.2f} nGy/h x{att:.2f} = {gm:.4f} mSv/yr full time", end="")
    if fD:
        j, cc, dr = radon(U, fD); print(f"; radon J {j:.3g} Bq/m2/s, C {cc:.3g} Bq/m3, {dr:.3g} mSv/yr (D_e factor {fD:.3g}, J range 1 to {math.sqrt(fD):.1f} times the 1-bar value)")
    else: print("; no atmosphere, no radon")
    if k == 'B': print(f"   ground bracket, upper crust in place of primitive mantle: gamma {eff(gam(*UCC), 1):.3f} mSv/yr, radon {radon(UCC[0], fD)[2]:.2f} mSv/yr")
    print(f"   total derived {gc + gm + dr:.4g} mSv/yr" + (" (GCR is a stand-in)" if k in STAND else ""))
    gt, ut = today(X, 14.9*w['M'], po2, w['d']); rn = .5 if w['P'] >= 1000 else 0.
    print(f"   today: GCR {gt:.3g} mSv/yr, ground {.5*(.48+rn):.2f} to {2*(.48+rn):.2f} mSv/yr, UV index {ut:.4g}")
    if k == 'A': print(f"   GCR {gc:.3f}; open-air full-time total {gc+gm+dr:.3f} mSv/yr; world-mean rebuild: GCR + external {eff(.2*g0+.8*1.4*g0, 1):.3f} + radon {40*.4*9e-6*.8*8766+dr*.2:.3f} + ingestion 0.29 = {gc+eff(.2*g0+.8*1.4*g0, 1)+40*.4*9e-6*.8*8766+dr*.2+.29:.2f} mSv/yr (UNSCEAR 2.4)")
```

## 13. Test targets (Solar System bodies are targets, never presets)

| Target | Value and source | Proposed tolerance |
|---|---|---|
| Earth air Rayleigh tau at 550 nm | 0.09721 (Bucholtz 1995 Table 4, read; replaces the computed 0.0969, which is 0.3 percent lower) | 1 percent |
| Earth greenhouse | T_s 288 K at tau0 1.82 to 1.9 (R&C, READ); G = 0.625 | T_s within 3 K |
| Ozone column | 279 to 300 DU (Cooke 279, Segura 313, Kozakis 192) | inside the seeded 0.7 to 1.65 range |
| CO2 doubling | 3.7 W/m2 (Myhre 1998 3.71; Byrne and Goldblatt 3.69); tau0 +5.0 to +5.1 percent, +1.12 to +1.13 K at fixed absorbed flux (computed) | 3 percent on the forcing, 20 percent on the shift |
| Venus-like | tau0 about 400, G about 101 (n = 1, R&C); solver 360 | 15 percent in tau0 |
| Surface UV | 300 DU overhead crude global UVI 13.9, direct 7.27 (part 02 2.7) | plus or minus 50 percent (crude) |
| GCR | sea level 0.32 mSv/yr; 265 g/cm2 aviation 47.9 mSv/yr (Roy 2021); Moon 500, Mars 232 (Hassler, GCR only, 21 g/cm2; solar-cycle range 128 to 420) | 10 percent; anchors per R13 |
| Mars cutoff | 165 MeV at 21 g/cm2 vertical, 247 at 60 degrees; absorbed 0.21 mGy/day (76 mGy/yr), 0.64 mSv/day (232 mSv/yr), Q 3.05 (Hassler 2014) | 3 percent on the cutoff, 15 percent on the doses |
| Ground | 60 nGy/h soil; 0.016 Bq/m2/s flux; 10 Bq/m3 outdoor radon; total 2.4 mSv/yr world mean | 10 percent on the total |
| Heat | primitive mantle 21.0 TW (2003 edition values, computed); 19.95 TW from the older values of part 07 2.9 | 10 percent |
| Upper crust gamma | 71.5 nGy/h from Rudnick and Gao Table 11 (U 2.7, Th 10.5 ppm, K 2.32 wt percent), against soil 60 | 10 percent |
| Titan | column 10850 g/cm2, GCR dose about 3e-6 Gy/yr, tau0 about 2.5 at 1 bar | order of magnitude |
| Airless ice moon | almost no radiogenic dose | below 1 percent of the GCR term |
| Runaway | 282, 291, 310 W/m2 | reject or steam above |

## 14. Contradictions and gaps found between sources

1. SC states an Earth air column of 2.148e26 cm-2; the same tau_R of 0.0969 follows from 2.148e25 (part 04 uses 2.22e25 for N2). SC has the exponent typo.
2. Ozone: the Kozakis median (192 DU) sits below Earth's 279 to 313 DU; the lane's plus or minus 30 percent does not span it (R8).
3. Part 03 5.4: the table rows below tau0 of about 0.6 follow the continuity root, while its Mars rows use the radiative formula (R11).
4. Part 03 2.2 prints the Byrne and Goldblatt CO2 form with a constant 0.39; check before coding.
5. W(air): part 07 2.12 lists 33.8 (Knoll) and 35 (Atri); ICRU 90 gives 33.97 (10B 5) and wins.
6. The master says the SEP spectra are not reached; the Band form and two Atwell events are read (10B 1.2a). Open: the Weibull parameters and 32 of 34 events.
7. The Moon (500) and Mars (232, Hassler) anchors are dose equivalents with the ICRP 60 Q(L), not ICRP 103 effective doses; the shape below 100 g/cm2 is unknown, and the code's falling curve disagrees with CORSIKA, which rises to the Pfotzer maximum.
8. UNSCEAR soil values are READ in part 07 and RECALLED in body part E; READ wins. Part E's 0.5 mSv/yr radon equals a full-time outdoor value, not a world mean.
9. `uv_index` is a ratio to Earth with an exponent -1.6 on a square-root ozone proxy (net p_O2^-0.8), not the CIE index and not DU^-1.08.
10. Papers against our files (boss-papers report, section 3): the Myhre labels (item in section 4), the Earth tau_R 0.0969 against 0.09721, the generic Rayleigh error, the Mars column of the code point (16.4 against 21 g/cm2), Palme and O'Neill 2014 against 2003 values, Wolff without tables. Printing faults to avoid when coding: Sneep and Ubachs Eq. 13 (CO2) is 1000 off as printed, Eq. 17 (CO) gives half its own table, Eq. 23 (O2) lacks 10^12, Table 1 prints the N2 delta as 0.974 for 0.0974; Bates Table 1 F'_K 1.047 at 0.5 micron against 1.049; Myhre Table 2 total; the Hassler 25 percent SEP sentence (absorbed dose only). The working forms of SC already avoid most Rayleigh faults.
11. BC13 (wR, wT, dose limits) is closed from the free IAEA GSR Part 3 except the neutron weighting function; ICRP 103 is not needed (Decided (Boss, contradiction sweep card, 7:29 AM EDT 2026-10-21)). The uploaded ICRP extract has no tables. The values here stay partly verified until read from the IAEA text.
12. Open with no seeded range: the line part of tau0, belts of other bodies, Xapsos spectra, annual SEP rate, radon D_e against pressure, X_EM for light gases.

## 15. Source conflicts and the rule that decided each

Rules (Boss, 12:39 PM EDT): V = verified wins; N = newest wins between two verified sources in doubt; M = most permissive (the widest value or range) when both are in question. M is 'most permission take wins', which Boss confirmed at 1:12 PM EDT means the widest range (Decided; composition-model.md section 4a item 12).

| Conflict | Sources | Settled as | Rule |
|---|---|---|---|
| Earth Rayleigh tau at 550 nm | computed 0.0969 (SC, part 02) against Bucholtz Table 4 0.09721 | 0.09721 is the test target | V |
| Rayleigh error of the generic formula | part 02 'about 5 to 9 percent' against Sneep and Ubachs per gas (2.5 to 18 percent) | per-gas errors | V |
| Mars GCR dose | Guo 255 mSv/yr (range 128 to 420) against Hassler 232 mSv/yr at 21 g/cm2 | centre 232, range 128 to 420 kept | V for the centre, M for the range |
| Mars column of the code point | 16.4 g/cm2 against Hassler 21 g/cm2 | 21 g/cm2 | V |
| Myhre coefficients | part 03 'IPCC 1990' labels against Myhre Table 3 (IPCC 1990: 6.3, 0.14) | labels read Myhre 1998 | V |
| CO2 doubling | Byrne and Goldblatt 3.69 against Myhre 3.71 | 3.7, B&G drives the shift, Myhre is the check | N |
| Primitive-mantle Th, U | 'Palme and O'Neill 2014' 84.9, 22.9 ppb (not confirmed) against 2003 edition 83.4, 21.8 ppb (read) | 2003 edition, labelled | V (N cannot apply: 2014 unread) |
| Crust Th, U, K | Palme and O'Neill crust column against Rudnick and Gao Table 11 | Rudnick and Gao | V |
| Primitive-mantle oxide total | 100 against 98.41 | 98.41 | V |
| Soil activities and radon | body part E RECALLED against part 07 READ UNSCEAR | part 07 | V |
| Gamma coefficients | 5.72 (part E) against 5.70 (part 07 2.10) | 5.70 | V |
| W(air) | 33.8 (Knoll), 35 (Atri), 33.85 (Boss) against ICRU 90 33.97 | 33.97 eV | N |
| Neutron wR | Boss '5 to 20' against the ICRP 103 function (secondary source) | continuous function, 2.5 to 20.7 | V (partly) |
| SEP spectrum form | Boss's Band form (discontinuous) against the integral rigidity form | Atwell form | V |
| Mars dose (SNIPPET 0.7 mSv/day against 10B 0.64) | 10B and Hassler | 0.64 mSv/day | V |
| Titan column and dose | Boss 1000 g/cm2 and 1.5 mGy/yr against 10850 g/cm2 and 3e-6 Gy/yr | the latter | V |
| Europa surface dose | Boss 1e3 to 1e4 Gy/yr against a 2 to 4 cm depth value | not a surface value | V |
| Ozone spread | Kozakis 192 DU with plus or minus 30 percent against observed 279 to 313 DU | times 0.7 to 1.65 | M |
| Small-tau0 closure | continuity root 0.25 against radiative 0.12 for the thin-air world | report both, bracket 0.12 to 0.25 | M |
| Radon against pressure | 1 kPa gate (no source) against D_e bracket | gate dropped, J factor 1 to sqrt(P0/P) | M |
| Mars dust omega0 | 0.92 to 0.95 (Boss) against solar mean 0.89 to 0.94, UV 0.67 to 0.74 | ranges stated per band, Wolff has no tables | M |
| Moon and Mars anchors below 100 g/cm2 | dose-equivalent anchors with different Q against each other | both kept as stand-ins with the range between | M |

## 16. Approved by Boss

**Decided (Boss, 2026-10-17, 9:40 to 9:55 PM EDT, cards answered in the project chat, relayed by the coordinator): rows 1 to 3 by tap, rows 4 to 12 and 14 by his typed "I approve all your defaults" at 9:54 PM EDT; every ruling is approved as recommended.** Row 13 needed none. Row 6 approves fetching ICRP 103 sections 4 and 5.10; the pages are still Boss's to supply (the secondary values are used meanwhile). Each item was a decision taken as default or left open, with the recommendation in the second sentence:

1. **R6 band edges** (100, 280, 315, 400, 700 nm, 4 micron). Approve; 4 micron leaks under 1 percent for Sun-like and cool-Earth blackbodies.
2. **R12 stored dose quantity**, open-air full-time effective dose, with the UNSCEAR occupancy rebuilt only as the Earth test. Approve; the rebuild gives 2.22 against 2.4 mSv/yr.
3. **R8 ozone spread**, times 0.7 to 1.65 of the Kozakis median. Approve; the symmetric 30 percent cannot reach Earth.
4. **R13 stand-ins** for the GCR dose below 100 g/cm2 (Moon 500, Mars 232 from Hassler, range 128 to 420) until a shower run. Approve as flagged stand-ins; ask the Boss list item 6 (CORSIKA or Geant4 runs) to replace them.
5. **R10 and R11 greenhouse**: keep the stored GEN.412 multiplier authoritative until correlated-k grids are built at build time (R2), report tau0 implied and the CIA floor. Approve; the per-class multiplier ranges then become labelled stand-ins and are retired when the grids arrive.
6. **R3 BC13 weights and limits**: read wR, wT and the dose limits from the free IAEA GSR Part 3; only the neutron weighting function stays open, and ICRP 103 is not needed (Decided (Boss, contradiction sweep card, 7:29 AM EDT 2026-10-21)). The secondary values are used meanwhile. Hassler 2014 is read.
7. **R14 SEP**: replace the 20 mSv/yr constant by an event dose with the Band spectrum; the annual rate stays open. Approve; do not store an annual SEP dose until the flare-frequency law is sourced for non-solar stars. (Built by GEN.747: the Sun's yearly dose from 47 measured events, other stars scaled by flare rate; stored as before.)
8. **R16 radon**: drop the 1 kPa gate and use the D_e bracket. Approve.
9. **R17 stellar UV**: seeded UV-C share 0.1 to 1 of blackbody. Approve as an open bracket and ask for a UV-against-T_eff source.
10. **R18 belts**: omit trapped-belt dose and flag it. Approve. (Superseded by GEN.749: a belt model for any dipole, asked per orbit, not stored.)
11. **R15 K depletion** (log-uniform 0.03 to 1) and the U/Th placeholder dex. Approve as placeholders, replaced by GEN.459 sourcing.
12. **R19 ordering**: GEN.453 first, then GEN.443 and GEN.444 together with the epoch bump. Approve.
13. **Conflict rule M**: Decided by Boss at 1:12 PM EDT: 'most permission take wins' means the widest range of the two sources (section 15 lists the six rows it decided or shaped). Nothing to approve.
14. **Ground reference**: use the Palme and O'Neill 2003 edition values (labelled) for the primitive mantle and Rudnick and Gao for the crust, and hold the 2014 values until the edition arrives. Approve.
