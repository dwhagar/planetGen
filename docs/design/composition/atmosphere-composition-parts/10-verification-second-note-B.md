# 10 Verification of the second boss-added note, part B (claims 1 to 9: particles, doses, ionisation, kinetics)

Scope: read-only verification of nine claims in the owner's second research note. Labels follow the worker rules (READ, 2ND-HAND, DATA, COMPUTED, SNIPPET, RECALLED, NOT REACHED). Verdicts: CONFIRMED, PARTLY, CONTRADICTED, UNSUPPORTED (searched, nothing found either way), NOT REACHED.

## 0. Summary table (details and sources in sections 1 to 9)

| # | Claim | Verdict | Evidence in one line | Label | Corrected value |
|---|---|---|---|---|---|
| 1a | Band form J(>E) as written | CONTRADICTED | high branch has an extra x^-(g2-g1); discontinuous at break unless g2-g1 = 1; slope -(2g2-g1) | COMPUTED, 2ND-HAND | integral in rigidity, J0 R^-g2 [(g2-g1)R0]^(g2-g1) e^(g1-g2) above R1 = (g2-g1)R0 |
| 1b | Weibull/exp-in-rigidity underestimate above 100 MeV | PARTLY | exp in rigidity under-estimates in 12 of 12 GLEs; Weibull varies, agrees only Feb 1956; Xapsos is Weibull, not rigidity-exponential | READ, 2ND-HAND | event dependent |
| 1c | NM data support Band | CONFIRMED | F80/F200 6.02 vs 6.15 (PAMELA, GLE 71) | READ | |
| 1d | Original ICRC paper | NOT REACHED | form read in four later papers, two by the authors | | |
| 2a | Column 20 to 22 g/cm2 | CONFIRMED | Gale 18.9 to 25.7 | READ | |
| 2b | Cutoff 150 MeV/nuc | PARTLY | [R1727, R1728] says about 150 (140 to 160); PSTAR CO2: 165 MeV at 21 g/cm2 | READ, COMPUTED | 165 vertical, 247 at 60 deg |
| 2c | 0.21 mGy/d, 0.64 mSv/d, 76 and 232 per yr | CONFIRMED | +/-0.04 and +/-0.12; 76.7 and 233.6 per yr | 2ND-HAND | solar-cycle factor 3.3 |
| 2d | Neutrons 15 to 20 percent of dose | CONTRADICTED | RAD 2.2 to 7 percent of dose, 3.8 to 10 percent of dose equivalent; model DE 33 to 50 percent | READ | see 2.3 |
| 3a | Trailing-hemisphere electrons | CONFIRMED | below 25 MeV trailing | READ | |
| 3b | 1 to 10 MeV electrons reach mm | CONTRADICTED | CSDA in ice 0.48 to 5.4 cm | COMPUTED | mm is 0.1 to 0.5 MeV |
| 3c | 1e3 to 1e4 Gy/yr | CONTRADICTED | literature compilation 1e7 at surface, 3e5 at 1 mm | READ (figure) | 1e3 to 1e4 at 2 to 4 cm |
| 3d | X-rays decimetres | CONFIRMED | 6 to 15 cm (0.1 to 1 MeV), yield 0.4 to 4 percent | COMPUTED | |
| 3e | 5.4 Sv/day from Nordheim | CONTRADICTED | traced to Zubrin book table (forum quote); equals 1.97e3 Gy/yr for wR = 1 | SNIPPET | no primary source |
| 3f | [R0879, R3308, R3452] tables | NOT REACHED | only via [R3038] figure | | |
| 4a | Titan column 1000 g/cm2 | CONTRADICTED | 1.467 bar / 1.352 = 10,850 g/cm2 | COMPUTED | 10.85 kg/cm2 |
| 4b | 1 to 3 GeV | CONTRADICTED | Earth column 2.25 GeV; Titan about 22 GeV ionisation, 121 interaction lengths, secondaries from over 100 to 350 GeV primaries | COMPUTED, READ | |
| 4c | 1.5 mGy/yr, 3 to 5 mSv/yr | CONTRADICTED | 1e-11 rad/s = 3.2e-6 Gy/yr | READ | factor 475 lower |
| 4d | Muons and neutrons dominate | PARTLY | muons and photons; neutrons do not reach the surface | READ | |
| 5a | W 33.85 eV | PARTLY | ICRU 31 value; ICRU 90 33.97 +/- 0.12 | 2ND-HAND | 33.97 |
| 5b | W above 1 keV constant, rises below 100 eV | PARTLY | W(air) 35.5 at 1 keV, 42.5 at 100 eV, 125 at 20.7 eV | READ | |
| 5c | Photoionisation sources | PARTLY | Leiden has HCN, H2CO, HCl, SiO, 15 metals; no noble gas, HF, CH3Cl, TiO, VO | READ | |
| 5d | Kr, Xe, HF, CH3Cl, TiO, VO cross sections | NOT REACHED | | | |
| 6a | N = 0.75 - 1.27 log10(10F) | CONTRADICTED | N = -0.238 for F = 0.6 | COMPUTED | N = 1.032 (IUPAC) or 1 (JPL) |
| 6b | "Standard JPL" N formula | PARTLY | that is IUPAC/Troe; JPL fixes N = 1, Fc = 0.6 | READ | |
| 6c | F about 0.6 | CONFIRMED | JPL fixed 0.6 | READ | |
| 7a | 19-5 updated ClO + O, + NO, dimer | CONTRADICTED | evaluated 10-6, 10-6, 10-6 (bimolecular); dimer 15-10 | READ | |
| 8a | New PH3 VUV cross sections | UNSUPPORTED | all work uses [R0762] | READ | |
| 8b | PH2 + H + M -> PH3 + M | PARTLY | in [R0233] as reverse of thermolysis; two-body in Sousa-Silva | READ | |
| 9 | Harman, Zahnle, Hunten updated by new cross sections | UNSUPPORTED | not meaningful as written; H2O cross sections revised by [R3662] | READ | |

Counts (29 sub-claims): CONFIRMED 6, PARTLY 8, CONTRADICTED 10, UNSUPPORTED 2, NOT REACHED 3 (1d, 3f, 5d). The JPL 19-5 licence question (section 7) was answered: not public domain.

## 1. Claim 1: Band function for SEP spectra

### 1.1 Published form


| Source | What it shows | Label |
|---|---|---|
| [R2355] slides, NTRS 20100005687 ("Radiation dose assessments of SPEs with spectral representation at high energies") | Band function in RIGIDITY, 4 parameters (J0, gamma1, gamma2, R0). Integral: Phi(>R) = J0 R^-g1 exp(-R/R0) for R <= (g2-g1) R0; Phi(>R) = J0 R^-g2 [(g2-g1) R0]^(g2-g1) exp(g1-g2) for R >= (g2-g1) R0. | 2ND-HAND (text extracted from the slide PDF; formula typography garbled but parseable) |
| [R4478], "Why is solar cycle 24 an inefficient producer of high-energy particle events?" (A&A 604, A47; [R4478]) eq. 4; Tylka and Dietrich are co-authors | the integral rigidity spectra of proton fluences in GLEs are double power laws: J(>R) = J0 R^-g1 exp(-R/R0) for R < R1 = (g2-g1) R0; J(>R) = J0 R1^-g1 exp(-R1/R0) (R/R1)^-g2 for R >= R1. "The Band function and its derivative are continuous." | 2ND-HAND (READ in arXiv PDF; by the original authors) |
| [R0562] eq. 2 | Same shape in kinetic energy, break (gb-ga)E0 (a different parametrisation, different fitted values) | 2ND-HAND (READ) |
| [R0173] | Band is "a double power law in rigidity with a smooth exponential junction"; Tylka and Dietrich fitted about 59 GLEs | 2ND-HAND (READ) |

Corrected statement: the fits are INTEGRAL event-integrated fluence spectra in RIGIDITY R (GV): J0 in cm^-2, gamma1 low-rigidity index, gamma2 high-rigidity index, R0 the exponential rigidity, break R1 = (gamma2 - gamma1) R0.

Normalisation identity (COMPUTED): with C = J0 R0^-g1 and x = R/R0, the correct high branch is C (g2-g1)^(g2-g1) e^(g1-g2) x^-g2; the note's high branch is this times x^-(g2-g1), hence discontinuous at the break and too steep above it.

### 1.2a Published Band parameters (fills the part 07 gap) and a validation (READ, COMPUTED)

[R3803], ICES (NTRS 20150009484), Table 1 (read from the rendered page), Band fit of the event-integrated integral fluence in rigidity, J(>R) with R in GV, J0 in protons/cm2, R0 in GV:

| Event | J0 (cm-2) | gamma1 | gamma2 | R0 (GV) | break R1 = (g2-g1) R0 |
|---|---|---|---|---|---|
| GLE 1976-04-30 | 7.620e6 | 1.710 | 6.33 | 0.2054 | 0.949 GV (396 MeV) |
| GLE 2006-12-13 | 1.330e8 | 1.045 | 5.80 | 0.1765 | 0.839 GV (321 MeV) |
| sub-GLEs (14 events 1974 to 2006) | 2.4e4 to 2.0e8 | -0.19 to 4.08 | 4.35 to 10.20 | 0.084 to 0.364 | |

Validation (COMPUTED): with the published form and the GLE 2006-12-13 row, F(>30 MeV) = 1.53e8 cm-2, F(>200 MeV) = 5.46e6, ratio F30/F200 = 28.0; [R0173: Table 1] ([R0173], READ) lists 26.1 +/- 6.8 for GLE 70 (13 Dec 2006). The 7 percent agreement confirms that the parameters are integral, rigidity-based and that the published form above is the one in use. With the note's formula the same event is lower than the published curve by 5.6e-4 at 330 MeV, 1.7e-4 at 500 MeV, 2.1e-5 at 1 GeV and 2.0e-6 at 2 GeV (COMPUTED), so the note's version would remove essentially every GeV proton. Proton rigidity R = sqrt(T^2 + 2 T m_p) / 1 GV with m_p = 938.272 MeV (COMPUTED): 30, 100, 200, 500, 1000 MeV = 0.239, 0.445, 0.644, 1.09, 1.70 GV.

Generator reading (COMPUTED): for the 2006 GLE the share of above-30-MeV fluence passing a Mars-like column (primaries above 165 MeV) is 5.7 percent. Part 07's single exponential with R0 = 0.2 GV normalised at 30 MeV is 2.2, 3.7 and 7.2 times the Band curve at 100, 200 and 500 MeV for this event, and Band R0 spans 0.019 to 0.364 GV over 34 events (READ).

### 1.2 The note's version fails continuity (COMPUTED)

The note's low branch C (E/E0)^-g1 exp(-E/E0) is correct; its high branch carries an extra (E0/E)^(g2-g1). Computed at the break E_b = (g2-g1)E0:

| g1, g2 | note: high/low ratio at break | published: ratio at break | note's asymptotic slope | published slope |
|---|---|---|---|---|
| 1.0, 4.0 | 0.0370 (= (g2-g1)^(g1-g2)) | 1.0000 | -(2 g2 - g1) = -7 | -4 |
| 1.0, 2.0 | 1.000 (only continuous when g2-g1 = 1) | 1.0000 | -3 | -2 |

Verdict claim 1a (formula): CONTRADICTED as written (not the published form; discontinuous at the break except for g2 - g1 = 1; wrong high-energy slope). The low-energy branch and the break definition (g2-g1)E0 are correct. The published form is continuous in value and in log-slope (low-branch slope -g1 - R/R0 equals -g2 at R1; COMPUTED).

### 1.3 Weibull and exponential-in-rigidity versus Band above 100 MeV

| Statement | Evidence | Label |
|---|---|---|
| Exponential in rigidity grossly under-estimates both Band and Weibull in all 12 GLEs studied | no source; Weibull varies considerably against Band, reasonable agreement only for 23 Feb 1956 | READ (NTRS PDF) |
| Exponential/Weibull extrapolation in proton energy underestimates the proton spectrum above 100 MeV | [R2355] slides concluding remarks: "exponential/Weibull extrapolation in proton energy underestimates the actual proton spectrum above 100 MeV"; but 23 Feb 1956 exponential OVER-estimates fluence, Weibull and Band agree; Aug 1972 Band LOWER at 60-200 MeV because satellite data there (IMP8 channel) read high; Sept/Oct 1989 effective dose similar for all three | 2ND-HAND (NTRS 20100005687) |
| "Weibull/exponential-in-rigidity [R4789]" | Mislabelled: [R4789] is the Weibull (stretched exponential in kinetic energy, Phi = Phi0 exp(-k E^alpha)); the exponential in rigidity is [R1454] / Malitson and Webber. Two different fits | READ |
| GLE neutron-monitor data support Band | Band is fitted with satellite data below about 100 MeV and world-wide neutron monitor data above about 430 MeV (1 GV); fit of about 59 to 66 GLEs; [R0173] tests the Band parametrisation against PAMELA direct data for GLE 71: F80/F200 ratio 6.15 +/- 0.9 (PAMELA) vs 6.02 +/- 1.38 (Band); agreement below 1 GeV, PAMELA errors grow above 1 GeV from GCR subtraction; PAMELA was not used in the fit | READ [R0173] |
| Limit of Band | Reduces to a single power law at high energy and misses rollovers seen by PAMELA; [R0562] multiplies Band by exp(-E/Er) | 2ND-HAND [R0562] |

Verdict claim 1b: PARTLY (underestimate above 100 MeV is true for exponential-in-rigidity in all 12 tested events and for Weibull in most, not all; Xapsos is Weibull, not rigidity-exponential). Verdict claim 1c (NM data support Band): CONFIRMED with the single-event PAMELA test.

Reconciliation with part 07 section 2.8: its Band and Weibull parameter gap is partly filled by the rows in 1.2a (two events plus ranges); the full 34-event table is in NTRS 20150009484 and was not transcribed.

## 2. Claim 2: Mars RAD numbers

### 2.1 Proton energy needed to cross the Gale column (COMPUTED from NIST PSTAR, material 134 carbon dioxide, CSDA range)

PSTAR table pulled directly (physics.nist.gov, ap_table-t.pl, matno 134, 132 energy rows 1 keV to 10 GeV) and inverted by log-log interpolation:

| Column (g/cm2 CO2) | Path | Proton energy to reach end of CSDA range (MeV) | Label |
|---|---|---|---|
| 17.8 | vertical | 150 | COMPUTED (DATA) |
| 18.9 (Gale low, [R1729]) | vertical | 155 | COMPUTED (DATA) |
| 21 | vertical | 165 | COMPUTED (DATA) |
| 22.5 (part 07) | vertical | 172 | COMPUTED (DATA) |
| 25.7 (Gale high, [R1729]) | vertical | 185 | COMPUTED (DATA) |
| 21 | zenith 45 deg (path 29.7) | 202 | COMPUTED (DATA) |
| 21 | zenith 60 deg (path 42.0) | 247 | COMPUTED (DATA) |

So the note's 150 MeV is the CSDA energy for about 17.8 g/cm2, 10% low in energy for 21 g/cm2 (165 MeV) and 38% low at 60 degrees (247 MeV).

### 2.2 What the papers state

| Item | Source statement | Label |
|---|---|---|
| Cutoff | [R1728, R1729]: at about 20 g/cm2 "Primary protons with energies less than about 150 MeV ... lack the range to reach the surface"; for He-4 "cutoff energy for incoming particles at about 150 MeV/nuc" (21 g/cm2); the atmosphere stops protons "(energies <= 140 to about 160 MeV depending on elevation)"; protons below the cutoff still feed the surface through secondary electrons and neutrons. Matches part 07 ("140 to 160 MeV") | READ |
| Column | [R1727, R1728] uses about 22.5 g/cm2 at Gale mean pressure (footnote: 5 to 6 g/cm2 more than at mean Mars elevation); [R1729]: Gale 18.9 to 25.7 g/cm2 (from part 07) | READ |
| Absorbed dose | Reitz (DLR) ICES-2016-437 (elib.dlr.de/106090): "absorbed dose measured in the E detector accounts to 0.21 +/- 0.04 mGy/d"; Hassler abstract in the same text; Q = 3.1 +/- 0.3; "dose equivalent of 0.64 +/- 0.12 mSv/d", 7 Aug 2012 to 1 Jun 2013; cruise 1.84 +/- 0.3 mSv/day. Also 2208.00892 (Crewed missions to Mars) quotes 0.210 +/- 0.040 mGy/day citing [R1843] | 2ND-HAND (science.org refused; the Crossref abstract has no numbers) |
| [R1730] | 0.7 mSv/day mean surface dose equivalent, 220 uGy/day absorbed (E detector fit 267 uGy/day, B detector 224 at reference), range 0.35 to 1.15 mSv/day over solar cycle | READ [R1730] |

Annual conversion (COMPUTED): 0.21 x 365 = 76.7 mGy/yr; 0.64 x 365 = 233.6 mSv/yr. The note's 76 and 232 are rounded consistently (about 0.5 percent low). CONFIRMED (numbers); the uncertainty is +/- 19 percent (0.04 of 0.21) and +/- 19 percent (0.12 of 0.64), and the solar-cycle spread is a factor 3.3 (0.35 to 1.15 mSv/day, [R1730]), so quoting 232 mSv/yr as one value is a solar-maximum anchor, not a constant.

### 2.3 Neutron share: "15 to 20 percent of total surface dose"

| Quantity | Value | Source | Label |
|---|---|---|---|
| Neutrons 8 to 740 MeV, first inversion | 14 +/- 4 uGy/d and 61 +/- 15 uSv/d: 7 percent of surface dose rate and 10 percent of dose equivalent | [R2419] (Crossref abstract, [R2419]) | READ (abstract) |
| Neutrons 7 to 740 MeV, improved inversion | 5.1 +/- 1.0 uGy/d and 23.6 +/- 4.1 uSv/d: 2.2 percent of dose, 3.8 percent of dose equivalent; extrapolated to 1 to 1000 MeV: 6.76 uGy/d and 31.24 uSv/d, about 3 percent and 5 percent | Koehler et al. [R1728] | READ |
| Model neutron share of absorbed dose in water at surface | 4, 8, 9 percent (82, 753, 1200 Pa) | [R4890] via part 08B section 4.3 | READ (in 08B) |
| Model neutron share of effective dose | about 50 percent at surface, above 80 percent below 50 cm | [R4890] via 08B | READ (in 08B) |
| GEANT4 calculation of dose equivalent | "more than one third in the dose equivalent calculation is due to the neutron component" with calculated neutron Q greater than 10 versus measured Q about 4 | Reitz ICES-2016-437 | 2ND-HAND |

Verdict claim 2d (15 to 20 percent of total surface dose): CONTRADICTED for absorbed dose (2 to 9 percent) and PARTLY for dose equivalent (RAD-measured 4 to 10 percent; model dose equivalent 33 to 50 percent depending on weighting). The statement is only defensible if it is stated as a model dose-equivalent fraction with ICRP 103 weighting, not as a fraction of dose. This agrees with 08B section 4.3 and adds the measured RAD numbers (2.2 to 7 percent of dose; 3.8 to 10 percent of dose equivalent).

Verdicts claim 2: column 20 to 22 g/cm2 CONFIRMED (Gale 18.9 to 25.7, mean about 22.5); 150 MeV cutoff PARTLY (right as Guo's quoted figure and as a rule of thumb, but PSTAR gives 165 MeV at 21 g/cm2 vertically and 247 MeV at 60 degrees; Guo's own range is 140 to 160); doses CONFIRMED (2ND-HAND); neutron share CONTRADICTED/PARTLY as above.


## 3. Claim 3: Europa

Primary papers: [R3308] is Nature Astronomy 2, 673 to 679 (published 23 July 2018; the note's year is right, the volume is not given). Only its abstract was reached (nature.com, preview); the dose-rate tables and figures are behind the paywall. [R0879], and [R3452] are on refused hosts; their numbers were read second-hand through [R3038] ([R3038: Figure 1], which combines [R0879]) and [R4827].

| Statement in the note | Evidence | Label | Verdict and corrected value |
|---|---|---|---|
| Trailing hemisphere is the electron target | [R3426] GRL abstract (Crossref): 10 keV to 10 MeV electrons "precipitate primarily on the trailing hemisphere". [R3038]: electrons below 25 MeV hit the trailing hemisphere preferentially, above 25 MeV the leading one (citing [R0879, R3308]); ion difference between hemispheres is only a few times. [R4827]: trailing deposition confined to the upper 0.1 cm or less; leading deposition reaches tens of cm | READ (abstract, arXiv) | CONFIRMED |
| Electrons of 1 to 10 MeV penetrate "mm depths" | ESTAR liquid water CSDA range (physics.nist.gov, matno 276): 1, 2, 5, 10, 25 MeV = 0.437, 0.979, 2.55, 4.98, 11.3 g/cm2; in ice of 0.917 g/cm3 = 0.48, 1.07, 2.8, 5.4, 12 cm. Mm depths belong to about 0.1 to 0.5 MeV electrons (0.02 to 0.19 cm) | COMPUTED (DATA) | CONTRADICTED for 1 to 10 MeV: mm to about 5 cm. Nordheim abstract: processing significant to about 1 cm at mid to high latitudes and 10 to 20 cm in the "radiation lenses" centred on the leading and trailing hemispheres (READ abstract) |
| Dose "about 1e3 to 1e4 Gy per year" | [R3038: Figure 1] (compilation of three studies, maximum of the interpolated values; READ figure, values eyeballed, uncertainty a factor of 3): electron dose rate about 1e7 Gy/yr down to 10 micrometres, 4e6 at 0.1 mm, 3e5 at 1 mm, 4e4 at 1 cm, 7e1 at 10 cm, 2 at 1 m; protons about 1e5 at 1 mm, 1e2 at 1 cm, 1 at 10 cm. The 1e3 to 1e4 Gy/yr band is reached at roughly 2 to 4 cm depth for electrons | READ (figure, approximate) | CONTRADICTED as a surface or mm-depth value (3 to 4 orders too low); correct only as a rate at about 2 to 4 cm in ice. Oxygen and sulfur ion doses differ by over 3 orders between the studies |
| Bremsstrahlung X-rays penetrate decimetres | XCOM water (physics.nist.gov, ComTab/water): mu/rho 0.171, 0.0969, 0.0707, 0.0222 cm2/g at 0.1, 0.5, 1, 10 MeV; attenuation length in ice 6.4, 11, 15, 49 cm. ESTAR radiation yield 0.36 percent at 1 MeV, 4.1 percent at 10 MeV, 10 percent at 25 MeV. [R4827] adds that radiative loss is small below the electron critical energy of water (78 MeV) | COMPUTED (DATA); READ (Yoffe) | CONFIRMED for range (decimetres), but the dose it adds is small (under 4 percent of the electron energy at 10 MeV); |
| "5.4 Sv/day" (part 07 flag) | Not in the Nordheim abstract and [R3308] is later than the figure's circulation. A 2002 forum post quotes Zubrin, "Entering Space" p. 167, table: Io 3600 rem/day, Europa 540 rem/day, Ganymede 8, Callisto 0.01; 540 rem = 5.4 Sv. A Centauri Dreams post (2021) repeats it and a reader comment says it is per day | SNIPPET (forum and blog text); not a primary source | CONTRADICTED as a Nordheim figure; origin is a popular book table (RECALLED title and page from the forum quote), no primary source found. 5.4 Sv/day = 1.97e3 Sv/yr (1.97e3 Gy/yr if wR = 1), 3 to 4 orders below the surface electron dose rate, so not a surface constant; it matches the dose at 2 to 4 cm depth |

Generator reading: a Europa-like surface dose is a depth profile (values above), with trailing deposition within 0.1 cm and leading deposition to tens of cm; ion (O, S) doses unresolved.

## 4. Claim 4: Titan

### 4.1 Column mass (COMPUTED)

X = P / g = 1.467e5 Pa / 1.352 m/s2 = 1.085e5 kg/m2 = 10,850 g/cm2 (10.85 kg/cm2). This is 10.9 times the note's 1000 g/cm2 and 10.5 times Earth's 1033 g/cm2 (101325 Pa / 9.80665 = 1033 g/cm2). P = 1.467 bar and g = 1.352 m/s2 are the note's inputs (the Huygens and Cassini values are RECALLED, not re-read here).

| Quantity | Value | Label |
|---|---|---|
| PSTAR N2 CSDA range reached by 10 GeV protons | 5030 g/cm2 (PSTAR upper limit 10 GeV) | DATA |
| Proton kinetic energy for 1033 g/cm2 N2 (Earth) | 2.25 GeV (agrees with part 07: 2241 MeV) | COMPUTED (DATA) |
| Proton energy to span 10,850 g/cm2 of N2 by ionisation alone (CSDA, extrapolated past 10 GeV with stopping power 2.06 to 2.27 MeV cm2/g) | about 22 to 23 GeV | COMPUTED, extrapolated |
| Nuclear interaction length N2 | 89.7 g/cm2 (PDG, part 07 table 2.6); 10,850 / 89.7 = 121 interaction lengths, survival exp(-121) = 1e-53 | DATA, COMPUTED |
| Earth for comparison: 1033 / 90.1 | 11.5 interaction lengths | COMPUTED |
| Primary energy giving significant surface energy deposit (Sittler GEANT3 and MGEANT/GCALOR model) | greater than about 350 GeV; mean 334 GeV; secondaries at the surface are mostly muons and photons, and "secondary neutrons, pions and protons do not reach the surface" | READ (NTRS) |

Verdict on "E_proton >= 1 to 3 GeV for the 1000 g/cm2 column": CONFIRMED for Earth's column (2.25 GeV) but CONTRADICTED when applied to Titan, whose column is 10.9 times larger: ionisation alone needs about 22 GeV and nuclear attack makes primaries irrelevant; the surface sees only secondaries from primaries above about 100 to 350 GeV.

### 4.2 Surface dose

| Source | Value | Label |
|---|---|---|
| no source | total dose rate about 1e-11 rad/s at the surface (tholin, density 1 g/cm3); 142 krad accumulated over 450 Myr (30 krad to 1350 krad over 100 Myr to 4.5 Gyr); surface energy flux 8.01e-6 erg/cm2/s; model lower limit, factor of 2 allowance for heavy ions and quality factor | READ |
| Converted | 1e-11 rad/s = 3.2e-6 Gy/yr = 0.003 mGy/yr; 142 krad / 450 Myr = 3.2e-6 Gy/yr (the two agree) | COMPUTED |
| Note | 1.5 mGy/yr and 3 to 5 mSv/yr | claim |
| Ratio | note / Sittler = 475 | COMPUTED |
| Earth sea level | 0.32 mSv/yr (part 07, READ), so the note's 1.5 mGy/yr would exceed Earth under 10.5 times less air | part 07 |
| Order-of-magnitude check | muons at about 1e4 g/cm2 with an assumed intensity of 1e-4 cm-2 s-1 sr-1 (RECALLED) and 2 MeV cm2/g give about 3e-6 Gy/yr, the Sittler scale | COMPUTED from a RECALLED input |

Verdicts: 1.5 mGy/yr and 3 to 5 mSv/yr CONTRADICTED (about 3e-6 Gy/yr in the only Titan model reached, GCR only); "muons and secondary neutrons dominate" PARTLY (muons and photons yes, neutrons do not reach the surface). [R1694], 2011 (A&A, refused host) and Dartnell Titan dose papers: NOT REACHED. Sittler: GCR deposition peaks near 100 km altitude; muons reach about 100 m into the ice.

## 5. Claim 5: W value, photoionisation sources

### 5.1 W for air (electrons)

| Item | Value | Source | Label |
|---|---|---|---|
| ICRU 31 (1979) dry air, electrons | 33.85 eV ("a mean value of 33.85 eV was recommended"; also cited as "34 eV" when comparing with Monte Carlo) | Fourkal, Nahum, Ma, Paskalev, [R1431] | 2ND-HAND (READ in that paper) |
| ICRU 90 (2016) | W_air/e = 33.97 +/- 0.12 eV (0.35 percent), "asymptotic value", electrons; protons 34.44 +/- 0.14 eV (previous 34.23); C ions 34.71 +/- 0.52 eV; ICRU 73 recommends 34.2 +/- 0.1 for protons | Seltzer CIRMS 2016 slides (ICRU 90 summary) and [R0592] table (33.97 eV +/- 0.12 eV) | 2ND-HAND (READ) |
| Spread | 33.85 to 33.97 = 0.35 percent; ICRU 90 supersedes | COMPUTED | |
| Low-energy electrons, air | W(1000 eV) 35.54; W(500) 36.81; W(100.7) 42.45; W(50.7) 51.32; W(20.7) 125.6; W(15.7) 640; threshold near 12.07 eV (O2). Ratio to the asymptote 33.9: 1.05 at 1 keV, 1.09 at 500 eV, 1.25 at 100 eV, 1.51 at 51 eV | [R0858], Radiat. Res. 84, 189, Table I and II (measured 5 to about 500 eV, joined to 1 keV) | READ (Helmholtz Munich PDF) |

Verdict: 33.85 eV CONFIRMED as the ICRU 31 value but superseded by 33.97 eV (ICRU 90; PARTLY). "Electrons above 1 keV W about 33.85": PARTLY: W(air) is 35.5 eV at 1 keV (5 percent above the asymptote), reaching the asymptote only above several keV (ICRU 90 states an asymptotic value). "W rises below 100 eV": CONFIRMED in direction but it rises gradually already below about 5 keV and is 25 percent high at 100 eV. For the generator: use the asymptotic W for a fast primary and its whole cascade (the asymptotic W already integrates the secondary electrons); use W(E) only for monoenergetic low-energy electrons (photoelectrons, auroral electrons below about 5 keV).

### 5.2 W for the other gases (Mars, Titan relevance); fills part 07 table 2.12 gaps

| Gas | High-energy W (eV) | W at 100 eV | W at 1 keV | Threshold (eV) | Relevance | Label |
|---|---|---|---|---|---|---|
| N2 (Titan main) | 34.8 | 43.3 | 36.7 | 15.57 | ion pairs per Gy 1.79e17 (part 07) | READ |
| CO2 (Mars, Venus main) | 33.0 | 40.3 (at 100.7) | 37.7 at 400 eV, no 1 keV row | 13.77 | 1.89e17 | READ |
| CH4 (Titan 5 percent) | 27.3 | 32.3 | 30.9 | 12.55 | 13 percent above asymptote at 1 keV | READ |
| H2O vapour | 29.6 | 35.2 | about 30.4 at 900 eV | 12.61 | fills the part 07 NOT REACHED entry | READ |
| H2 | 36.4 | not tabulated (41.2 at 200 eV) | 38.6 at 500 eV | 15.42 | | READ |
| Ar, Kr, Xe | 26.4, 24.3, 22.0 | 32.1, 28.2, 25.4 | Xe 23.8 | 15.76, 14.00, 12.13 | Ne not in Combecher | READ |

Not reached: NH3, SO2, CO, HCN, CH3Cl (no W source); alpha W in air (ICRU 31, RECALLED value not used). ICRU 90 proton W in air (34.44 eV) replaces part 07's guess of "3 to 10 percent above electron W": it is +1.4 percent for air.

### 5.3 Outer-shell photoionisation cross sections: which species have a reachable source

Leiden / [R1873] [R1873: Table 1] (READ, 102 species; threshold wavelengths and Ly-alpha cross sections listed; data and rates at the Leiden database, host not tested):

| Species | In Leiden Table 1? | Other reachable source | Verdict on source |
|---|---|---|---|
| Ar, Ne, Kr, Xe | No (no noble gases in the table) | [R4542] fit (all subshells, Z up to 30, so Ne and Ar only; used in part 07 as "V"); Henke f2 fails beyond 70 nm; Kr and Xe outer-shell data (Samson and Stolte; Marr and West) are RECALLED, not reached | Ne, Ar reachable via Verner; Kr, Xe NOT REACHED |
| HCN | Yes (H + CN dissociation, with ionisation entry) | part 02: HCN 0.1 to 187.9 nm in the photochem data | reachable |
| H2CO | Yes (H2 + CO) | part 02: 6.3 to 375.9 nm | reachable |
| HCl | Yes (also HCl+) | part 02: 4.5 to 229.9 nm | reachable |
| SiO | Yes | Leiden only | reachable |
| HF | No | Mainz atlas file HF.spc listed in part 02 (not read here) | NOT REACHED here |
| CH3Cl | No | JPL 19-5 photochemical section (host refused in this pass; the JPL 19-5 PDF itself was reachable at science.jpl.nasa.gov, CH3Cl cross-section table not extracted) | NOT REACHED |
| TiO, VO | No | no VUV photoionisation data; line lists only (ExoMol, RECALLED) | NOT REACHED |
| Metals | Li, Na, Mg, Al, Si, K, Ca, Ti, Cr, Mn, Fe, Co, Ni, Zn, Rb listed | NIST ASD line data (part 02) | reachable for these elements |


## 6. Claim 6: Troe falloff form

### 6.1 N for Fc = 0.6 (COMPUTED)

| Form | N | Remark |
|---|---|---|
| IUPAC/Troe: N = 0.75 - 1.27 log10(Fc) | 1.0317 | [R0182] (ACP 4, 1461, READ, p. 1469): log F = log Fc / (1 + [log(k0/kinf)/N]^2), N = [0.75 - 1.27 log Fc] |
| JPL 15-10 and 19-5 | fixed 1 (and Fc fixed at 0.6) | READ (see 6.2) |
| Note: N = 0.75 - 1.27 log10(10 F), F = 0.6 | -0.238 | wrong (log10(6) = 0.778, so 0.75 - 0.988); the factor 10 has no basis in the Troe or IUPAC forms |

Effect on F_eff = Fc^(1 / (1 + (x/N)^2)), x = log10(k0[M]/kinf), at x = 0.5, 1, 2, 3 (COMPUTED): N = 1 gives 0.6645, 0.7746, 0.9029, 0.9502; N = 1.0317 gives 0.6612, 0.7684, 0.8982, 0.9474; the note's N gives 0.9098, 0.9729, 0.9929, 0.9968. Worked test, JPL 15-10 OH + NO2 + M (k0 = 1.8e-30, n = 3.0, kinf = 2.8e-11, m = 0; READ, OCR), 298 K, [M] = 2.46e19 cm-3: JPL form 1.049e-11 cm3/s, IUPAC N 1.048e-11 (0.1 percent), the note's N 1.270e-11 (+21 percent), Lindemann 1.715e-11 (+63 percent). Verdicts: note's N CONTRADICTED; "standard JPL N formula" PARTLY (that is IUPAC/Troe; JPL fixes N = 1 and Fc = 0.6, within 1 percent in k for Fc = 0.6); F about 0.6 CONFIRMED (JPL fixed; [R0182] also lists 0.35 to 0.85).

### 6.2 JPL form as printed (READ)

JPL Publication 15-10 (Evaluation 18, NTRS 20160006735, rendered page 2-1): k0(T) = k0^300 (T/300)^-n, kinf(T) = kinf^300 (T/300)^-m; k_f([M], T) = {k0(T)[M] / (1 + k0(T)[M] / kinf(T))} 0.6^{ {1 + [log10(k0(T)[M] / kinf(T))]^2}^-1 }. Text: "The fixed value 0.6 that appears in this formula fits the data for all listed reactions adequately, although in principle this quantity may be different for each reaction". JPL Publication 19-5 (Evaluation 19, May 2020; science.jpl.nasa.gov/documents/1487/NASA-JPL_Evaluation_19-5.pdf, 1610 pages, reachable in this pass), eq. 2.3, text extracted: k_f = {kinf(T) k0(T)[M] / (kinf(T) + k0(T)[M])} 0.6^{ {1 + [log10(k0(T)[M]/kinf(T))]^2}^-1 }, same fixed 0.6; the temperature form changed to k0(T) = k0^298 (298/T)^n and kinf = kinf^298 (298/T)^m (reference 298 K, not 300 K), a coefficient-reading trap for the generator. 19-5 adds Table 2-2 for chemical-activation (CA) reactions (OH + CO etc.), which have a second pressure-dependent product channel. This corrects part 04's statement that the Lindemann form is the "upper bound": JPL itself is Lindemann times the 0.6-broadening factor, so the shipped Lindemann entries (no Troe keys) overestimate k near the centre by up to 1/0.6 = 1.67 (COMPUTED).

## 7. Claim 7: JPL 19-5 changes and licence

| Item | Evidence | Label | Verdict |
|---|---|---|---|
| "updated ClO + O, ClO + NO and ClO + ClO dimer kinetics" | JPL 19-5 Section 1.7 "Noteworthy changes": entire NOx section re-evaluated; FOx, ClOx, BrOx and IOx reviewed "with particular emphasis on the reactions of OH and Cl with hydrocarbons and halocarbons"; sulfur, metal (Na + CO, NaOH + H), Criegee, CA reactions moved to Table 2-2. ClO + O, ClO + NO and ClO + ClO are not named | READ (19-5 text) | not supported as stated |
| Per-entry evaluation tags in 19-5 | O + ClO (F1): "Table 10-6, Note 10-6, Evaluated 10-6"; ClO + NO (F130): Table 82-57, Note 10-6, Evaluated 10-6; ClO + ClO bimolecular (F137): Table 94-26, Note 97-4, Evaluated 10-6 (Nickolaisen); ClO + ClO + M termolecular (F10): Table 15-10, Note 15-10, Evaluated 15-10 | READ | CONTRADICTED: the three bimolecular entries were last evaluated in JPL 10-6; the ClO dimer termolecular entry was updated in JPL 15-10, not 19-5 |
| Values in 19-5 (298 K) | O + ClO: A 2.8e-11, E/R -85, k 3.7e-11, f 1.05; ClO + NO: A 6.4e-12, E/R -290, k 1.7e-11, f 1.15; ClO + ClO bimolecular: k1 (Cl2 + O2) A 1.0e-12, E/R 1590, k 4.8e-15; k2 (ClOO + Cl) A 3.0e-11, E/R 2450, k 8.0e-15; k3 (OClO + Cl) A 3.5e-13, E/R 1370, k 3.5e-15; branching 0.29 : 0.50 : 0.21 at 298 K; ClO + ClO + M: k0 1.9e-32, n 3.6, kinf 3.7e-12, m 1.6, k 3.2e-13, f 1.15 (same as 15-10) | READ | these close the ClOx gap in the part 04 mechanism: data are reachable and listed here |
| Authors | Burkholder, Sander and 12 co-authors; May 2020; NTRS 20210006316 is metadata only | READ | CONFIRMED (08B had this as RECALLED) |
| Copyright | JPL 19-5 title page: "Copyright 2019, California Institute of Technology. U.S. Government sponsorship acknowledged. All rights reserved." JPL 15-10: "Copyright 2015, California Institute of Technology. U.S. Government sponsorship acknowledged. All rights reserved." NTRS metadata for 19-5: "Items in TRS are protected by copyright, but are furnished with U.S. government purpose use rights"; for 15-10 the NTRS determination is PUBLIC_USE_PERMITTED | READ | NOT public domain; no open licence stated. Rate constants are facts, but the compilation and notes are copyrighted; cite and, for embedded tables, seek permission or regenerate from primary papers |

## 8. Claim 8: PH3

| Item | Evidence | Label | Verdict |
|---|---|---|---|
| "Updated VUV cross sections beyond [R4163] show shorter PH3 lifetimes" | [R4163] take PH3 UV cross sections at 295 K from [R0762] (via [R3661] and the Mainz atlas), quantum yield 1 below 230 nm, branching unity. [R0233] use the same [R0762] 295 K set and the same yield. Mainz atlas PH3 files ([R0762] at 295 K 120 to 230 nm and 155 K 127 to 200 nm) were seen only in a search listing. arXiv searches for PH3 photodissociation, VUV, cross section, lifetime returned no new cross-section paper; a 2024 Frontiers paper found is about PH+ and PO and is not relevant | READ (two papers); SNIPPET (Mainz listing) | UNSUPPORTED: no measurement newer than [R0762] was found; the "shorter lifetime" results ([R0233]: below 1 s above 78 to 98 km; transport-limited about 400 to 700 yr) are modelling with old cross sections. [R0233] states direct photolysis is not the dominant PH3 sink below 100 km on Venus |
| "needing continuous production" | PH3 destruction by OH, H, O, Cl and photolysis is the standing argument of [R0233, R4163] (surface or biological source on UV-shielded or CO2 worlds); lifetimes of 5 to 28 h with OH (part 04) | READ | CONFIRMED in general (not tied to new cross sections) |
| PH2 + H + M -> PH3 + M | [R0233] SI: this is "included as the reverse reaction for thermolytic decay of PH3" (PH3 + M -> PH2 + H + M), and they neglect other PH3 recombination. [R4163] give the two-body form PH2 + H -> PH3, k = 3.7e-10 exp(-340 K / T) cm3/s = 1.1e-10 at 288 K, important for H2-rich worlds. COMPUTED: 3.7e-10 exp(-340/288) = 1.14e-10 | READ; COMPUTED | PARTLY: both forms exist (three-body reverse in Bains, two-body in Sousa-Silva); neither is "new" |

## 9. Claim 9: "[R1814, R2050, R4860] are updated by incorporating revised cross sections"

| Paper | What it is | Cross-section content | Label |
|---|---|---|---|
| [R1814] | One-dimensional photochemical study of abiotic O2 from CO2 photolysis around F, G, K and M stars: O2 should not accumulate on F and G planets; K and M planets may build O2 because of low NUV; O3 detectable in UV for a broader set. Text states that the H2O cross section is extrapolated from 196 nm out to 240 nm (all models, including Hu et al.) | its own input H2O and CO2 sets; later revised | READ (IOP abstract and text) |
| [R4860] (ambiguous) | Two candidates: Zahnle K., 2006, "Earth's Earliest Atmosphere", Elements 2, 217 (review); or [R4860], Geobiology 4, 271 (1-D photochemistry of O2, CH4 and S; S-MIF loss by methane collapse; abstract names S8, not haze) | no cross-section update | READ (Crossref abstract of the Geobiology paper) |
| [R2050], J. Atmos. Sci. 30, 1481 | analytic escape theory of light gases | none | READ (Crossref title); cited by [R1814] for diffusion-limited H escape |

Judgement: as written the claim is not meaningful: published papers are not updated, [R2050] has no cross-section input, and "[R4860]" is ambiguous with neither candidate a cross-section paper. What is real: the H2O near-UV cross sections that [R1814] extrapolated have since been measured. [R3662] (ApJ 896, 148; Crossref abstract, READ): first NUV H2O absorption measurement at 292 K, "orders of magnitude more than previously assumed"; enhanced OH recombines CO and O2, suppressing both "by orders of magnitude relative to past predictions" and "eliminating the low-outgassing false-positive scenario for O2 around solar-type stars". Consequence for part 04: abiotic O2 and CO numbers from Hu 2012, [R1511, R1814] (2ND-HAND) pre-date these cross sections and overestimate O2 for anoxic CO2 atmospheres containing water. Verdict: UNSUPPORTED as stated; the underlying fact is CONFIRMED for [R1814] only.

## 10. What this changes in the generator items

Radiation items (part 07 F4, F5 and Mars, Titan, Europa rows):
- SEP source: use the published integral rigidity Band form (section 1.1), not the note's; event library from the rows (GLE 2006: J0 1.33e8, g1 1.045, g2 5.80, R0 0.1765 GV). Keep R0 a drawn parameter (0.019 to 0.364 GV).
- Atmospheric cutoff from the PSTAR CSDA inversion of X = P/g in the local gas, not a fixed 150 MeV: 165 MeV at 21 g/cm2 CO2 vertical, 247 MeV at 60 degrees. Dose below cutoff comes from secondaries.
- Thick columns: Titan X = 10,850 g/cm2 (not 1000). Ionisation alone needs about 22 GeV; the surface sees muon and photon secondaries only. The one calibration point (labelled as such) is a single value at 3.2e-6 Gy/yr. The note's 1.5 mGy/yr is not to be used.
- Neutron share: RAD 2.2 to 7 percent of dose and 3.8 to 10 percent of dose equivalent; 33 to 50 percent only as a model effective-dose share.
- W: asymptotic 33.97 +/- 0.12 eV for air (ICRU 90); 33.85 is superseded; W(E) table (35.5 eV at 1 keV, 42.5 at 100 eV) only for electrons below a few keV; per-gas W in 5.2.
- Europa-like: dose is a depth profile (1e7 Gy/yr at the surface to about 1e3 to 1e4 at 2 to 4 cm, electrons 1 to 10 MeV reaching 0.5 to 5 cm); drop "5.4 Sv/day".

Photochemistry: code Troe as k0[M]/(1 + k0[M]/kinf) times 0.6^(1/(1 + log10(k0[M]/kinf)^2)), reference 298 K in JPL 19-5 (300 K in 15-10); use the IUPAC N formula only for non-0.6 Fc, never the note's. ClOx bimolecular and dimer rate data in section 7 close the part 04 gap; JPL text is copyrighted, so regenerate or seek permission. Abiotic O2 and CO values for CO2 worlds with water are upper bounds [R3662]. PH3 cross sections remain [R0762] and the lifetime table of part 04 stands; PH2 + H has both a two-body and a three-body form.

Part 09 row d (retention regimes d1 to d5): add a derived column X = P/g so d3 to d5 pick the cutoff and dose regime above; no weight change. Row p (host UV regime p1 to p6): no weight change; add a flag that abiotic-O2 outcomes for p2 and p3 (K and M) are not reliable until [R3662] cross sections are in the chemistry.

## 11. NOT REACHED

- Weibull parameters.
- [R0879, R3308, R3452] dose tables (paywall or refused); [R1694], 2011; Dartnell Titan dose.
- Photoionisation of Kr, Xe, HF, CH3Cl, TiO, VO; W for NH3, SO2, CO, HCN, CH3Cl, Ne.

## 12. Recalled / to verify

- Zubrin, Entering Space, p. 167 as the origin of 5.4 Sv/day (title and page from a forum quote).
- Titan P = 1.467 bar and g = 1.352 m/s2 (the note's inputs) and a vertical muon intensity of about 1e-4 cm-2 s-1 sr-1 at 1e4 g/cm2.
- Kr and Xe outer-shell data (Samson and Stolte; Marr and West); TiO and VO line lists in ExoMol.
