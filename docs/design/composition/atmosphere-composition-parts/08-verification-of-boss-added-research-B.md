# 08 Verification of boss-added research note, part B (claims 1 to 6)

The owner-supplied note (machine-generated, stray "[cite: n]" markers) was treated as unverified. Labels follow _worker-rules.md: READ (original text), 2ND-HAND (read in a reachable paper citing it), DATA (data file), COMPUTED (re-derived or digitised by me), SNIPPET, RECALLED (listed in section 9), NOT REACHED. Verdicts: CONFIRMED, PARTLY, CONTRADICTED, UNSUPPORTED (no source found anywhere reachable), NOT REACHED. No repo file was edited.

---

## 1. Mars dust (claim 1)

### 1.1 omega0 about 0.92 to 0.95 at 670 nm. Verdict: PARTLY

[R4743] (JGR 114, E00D04, [R4743]) full text is paywalled (AGU/Wiley refused, Unpaywall closed); only its abstract was read. Values reached:

| Item | Value | Source | Status |
|---|---|---|---|
| [R4743] abstract | literature "dark" solar-band omega0 0.89 to 0.90, "bright" 0.92 to 0.94; new solar-band omega0 0.94 | Crossref abstract | READ (abstract) |
| Solar-weighted omega0 | [R3326]: 0.89; Dust value 0.94 (also [R0809, R4380] "nearly equal to 0.94"); dust 0.86 (g 0.79); GFDL 0.92 (g 0.75) | [R2800] (HAL hal-01120102), sections 2 and 3 | 2ND-HAND |
| Channel means 0.1-0.5 um / 0.5-5 um | Ockert-Bell 0.665 / 0.928; Wolff (gamma dist., reff 1.5 um, neff 0.3, cylinders) 0.738 / 0.955 | [R2800: Table 2] | 2ND-HAND |
| [R4743] near 650 nm | 0.975 ("average of results ... surface-orbit observations"); alternatives tested 0.940 [R4395], 0.930 [R3567] | Chen-Chen et al., [R0765, R0766] | 2ND-HAND |
| Near-IR 1.0 to 2.6 um | 0.97 (confirmed) | Vincendon et al., [R4555] | 2ND-HAND |

0.92 to 0.95 is the solar-band mean of the "bright" literature, not the 670 nm value (second-hand about 0.975), and it hides the UV drop (0.1 to 0.5 um band mean 0.67 to 0.74). Supported generator ranges: red/NIR 0.93 to 0.975, solar mean 0.89 to 0.94, UV band 0.67 to 0.74. The spectral omega0 table of [R4743] is NOT REACHED.

### 1.2 k greater than 0.01 below 350 nm. Verdict: CONFIRMED (figure reading); Fe3+ attribution NOT REACHED

[R2343] (JGR Planets, [R2343]; open copy atmos.washington.edu/~sgw/PAPERS/2021_Mars_dusty_ice.pdf) plot [R4743] imaginary index in Figure 1. I digitised the log-axis curve at 400 dpi (frame calibrated to 12 decades; uncertainty about 15 to 20 percent). COMPUTED from a 2ND-HAND figure:

| lambda (um) | 0.30 | 0.35 | 0.40 | 0.45 | 0.50 | 0.55 | 0.60 | 0.70 | 0.80 | 1.0 | 1.5 | 2.0 | 2.5 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| k | 0.0129 | 0.0114 | 0.0093 | 0.0074 | 0.0054 | 0.0036 | 0.0017 | 0.00095 | 0.0011 | 0.0018 | 0.0027 | 0.0040 | 0.0051 |

k crosses 0.01 near 0.38 um; the curve starts near 0.27 um. Minimum about 1e-3 at 0.7 to 0.8 um. The "Fe3+" cause appears in no text read: [R3326] (abstract, READ) "does not assume a particular mineralogical composition" and adopts the index of bright surface units. [R4395] abstract (READ): imaginary index rises steeply from 670 nm to shorter wavelengths. Lab contrast (DATA, ApJS 268:47, Table 3): analog k is only 0.0001 to 0.0013 over 200 to 2000 nm (JSC Mars-1: 0.00130 at 300 nm, 0.00047 at 650 nm); they state it is "much smaller" than airborne-dust literature values. Airborne k therefore spans about one decade depending on grain-shape assumptions.

### 1.3 Row "Volcanic SO2 aerosol: [R2645, R4395], 0.2-0.4 um". Verdict: CONTRADICTED (mis-citation)

[R4395]: "Properties of dust in the Martian atmosphere from the Imager on Mars Pathfinder" (mean radius 1.6 +/- 0.15 um). [R2645]: "Atmospheric Imaging Results from the Mars Exploration Rovers" (visible optical depth 0.9; dust radius 1.47 +/- 0.21 um at Gusev, 1.52 +/- 0.18 um at Meridiani; visible to 9 um optical depth ratio 2.0 +/- 0.2). Abstracts READ via Crossref/OpenAlex. Neither concerns SO2 or sulphate aerosol.

---

## 2. Titan haze (claim 2)

### 2.1 "Equivalent-sphere Mie underpredicts forward scattering by up to 40 percent". Verdict: UNSUPPORTED

Searched [R4392, R4394] team slides, [R0230, R0768, R2008, R4004] ([R0768], earlier mis-attributed here to Gao and Stam; corrected in part 10C section 0), snippets of [R2584, R4673]. No "40 percent" anywhere. Established: DISR data need fractal aggregates, 4000 monomers of 0.04 um radius ([R2008], READ) or N = 3000 of r = 0.05 um (Tomasko team slides, READ); the slides compare equal-volume sphere radius 0.72 um with equal-projected-area radius 2.03 um. Cheng et al. (Pluto, READ): aggregates give "a too small backscatter lobe by factor about 3" (different body and quantity). Treat the 40 percent as invented.

### 2.2 Absolute extinction law. Verdict: CONFIRMED (earlier exponents correct)

Source: Tomasko and the DISR team, CHARM workshop slides 2008-02-26 (atmos.nmsu.edu; READ; the PSS 56, 669 paper is NOT REACHED). Cumulative optical depth, lambda in nm:

- above 80 km: 1.012e7 lambda^-2.339; 30 to 80 km: 2.029e4 lambda^-1.409; below 30 km: 6.270e2 lambda^-0.9706.
- vertical shape above 80 km: exp(-(h-80)/65 km); linear below (Buenzli, [R0296] A.4, 2ND-HAND).

Unit check (COMPUTED): at 1500 nm the laws give 0.38, 0.68, 0.52 while [R0230] (Table 4, READ) fit 0.31, 0.34, 0.16 using stated scale factors 0.8, 0.5, 0.3; the ratios 0.82, 0.50, 0.31 match, so nm is correct.

| lambda (nm) | 355 | 430 | 500 | 550 | 675 | 830 | 934 | 1078 | 1288 | 1583 |
|---|---|---|---|---|---|---|---|---|---|---|
| tau above 80 km | 10.97 | 7.01 | 4.92 | 3.94 | 2.44 | 1.51 | 1.14 | 0.82 | 0.54 | 0.33 |
| tau total | 18.2 | 12.7 | 9.6 | 8.1 | 5.7 | 4.0 | 3.3 | 2.6 | 2.0 | 1.45 |

The slide plot reaches 18 to 20 at 355 nm (consistent). Single-scattering albedo is plotted (0.8 to 1.0) but not tabulated in text ([R4394: Table 2] NOT REACHED). [R1144] (Icarus 270, 355, "...radiance measurements made inside and outside the atmosphere"; NTRS 20170002309; abstract READ, no PDF): optical-depth wavelength dependence "similar ... from 500 to 700 nm, but smaller at shorter wavelengths and larger toward longer wavelengths"; scale height 65 km below 140 km, 45 km above; extra layer near 55 km; thin layers at 21, 11, 7 km ([R2263], real title "Rain and dewdrops on Titan based on in situ imaging"). New coefficients, [R2584, R3663] full text NOT REACHED ([R2584] abstract SNIPPET: charge 15 e/um, imaginary index larger than lab tholins).

---

## 3. Comet volatiles (claim 3)

Read: Harrington [R3535] ([R3535]; HP22); Biver, Dello Russo, Opitom, Comets III ([R0410: Table 3], updates [R0435], Dello [R1063, R3835]); [R4705] ([R4705], Dello [R1063] ranges); [R0859] ([R0859], SOHO/SWAN, 61 comets); [R0048] (NTRS 20120008332; text only, Table 1 OCR damaged). [R3206] (abstract only) and [R0436] NOT REACHED.

| Species | Claim | Found (percent of H2O) | Status | Verdict |
|---|---|---|---|---|
| CO | 0.15 to 30 | within 2.5 au 0.3 to 26, median 3.3 +/- 1.3 (32 comets); mean 6.4 (HP22), 5.2 +/- 1.3 (24 comets, Dello [R1063]); IR 0.3 to 26, radio < 1.23 to 35; OCC 0.4 to 26.2 (mean 6.1), JFC 0.3 to 4.3 (mean 1.6) | READ | PARTLY |
| CO2 | 2 to 30 | 2 to 30, median 12 +/- 2 (23 comets); IR 4 to 30; AKARI median 17 (quoted); CO+CO2 median 18 +/- 4 | READ | CONFIRMED |
| CH4 | 0.2 to 1.5 | IR 0.15 to 2.7; OCC 0.15 to 1.57 (mean 0.88); JFC < 0.25 to 0.54 (mean 0.31) | READ | PARTLY |
| NH3 | 0.2 to 1.0 | IR 0.1 to 3.6 (3.63 in C/2012 S1); radio 0.18 to 0.60; OCC mean 0.91, JFC mean 0.59, all 0.80 | READ | CONTRADICTED at top |
| Q(H2O) | 1e26 to 1e31 /s | SWAN peak values 5.11e26 (C/2012 E2) to 1.41e31 (Hale-Bopp) over 141 table entries; Garradd at 2 au 8e28 (IR) to 2.9e29 (Herschel, SWAN) | READ; COMPUTED min/max | PARTLY (low end is a SWAN detection limit) |

Also: CO2/CO varies by more than two decades and the JFC range encloses the LPC range ([R0048], READ); Q(CO2) 0.06 to 740 and Q(CO) 0.2 to 3000, units 1e26 /s, 0.71 to 6.18 au, 25 comets (HP22); Q(CO2)/Q(H2O) goes as r^2.2.

### 3.1 Differences against 05-mixing-ratios-by-class.md section 9 (not edited)

| Row in 05 | 05 says | Reached sources | Note |
|---|---|---|---|
| CO2 | median about 12, 25 comets | median 12 +/- 2 from 23 comets within 2.5 au; range 2 to 30 | add range |
| CO | 0.5 to 23, median about 3 | 0.3 to 26; median 3.3; Hale-Bopp 26.2 [R4705] vs 23.0 [R4704] | widen; source spread |
| NH3 | 0.1 to 0.73 | 0.1 to 3.6 (IR), 0.18 to 0.60 (radio) | upper bound too low |
| CH4, C2H6, H2CO | not extracted | CH4 0.15 to 2.7; C2H6 IR 0.1 to 2.7 (OCC mean 0.63, JFC 0.34); H2CO IR < 0.02 to 1.1, radio 0.13 to 1.4 | fill in |
| CH3OH | 0.2 to 2.6 | radio 0.7 to 6.1; IR < 0.13 to 4.3 | widen |
| HCN | 0.09 to 0.29 | radio 0.05 to 0.25; IR 0.03 to 0.5 | widen |
| H2S | 1.10 (67P) | radio 0.09 to 1.5 | cross-comet range missing |
| 67P | [R3835] snapshots | mission-integrated: CO 0.3 to 3, CO2 7.0, CH4 0.4, C2H6 0.8, NH3 0.4, H2S 2.0, O2 2.0 | CO2 4.7 vs 7.0 is a spread, not an error |
| 9.2 flux | 1.54e22 /m2/s; 4.8e28 /s for 1 km | COMPUTED: 1.543e22 and 4.85e28 (flux times pi R^2). Hale-Bopp 1.41e31 /s needs about 900 km2 sunlit area at the 1 au flux | consistent; L = 51 kJ/mol still RECALLED |

---

## 4. ICRP 103 weighting factors and subsurface neutron dose (claim 4)

### 4.1 wR. Verdict: CONFIRMED (numbers 2ND-HAND)

ICRP free extract (icrp.org, READ), summary (h): neutron and proton values changed, "values for neutrons given as a continuous function of neutron energy", charged pions added, "photons, electrons, muons, and alpha particles are unchanged". The ICRP 103 table is not reachable; values from Wikipedia "Equivalent dose" quoting it: photons, electrons, muons 1; protons and charged pions 2; alpha, fission fragments, heavy ions 20; neutrons (E in MeV): 2.5 + 18.2 exp(-[ln E]^2/6) for E < 1; 5.0 + 17.0 exp(-[ln 2E]^2/6) for 1 to 50; 2.5 + 3.25 exp(-[ln 0.04E]^2/6) for E > 50. COMPUTED: continuous at 1 MeV (20.69) and 50 MeV (5.50); 0.01 MeV 3.03; 0.1 MeV 10.0; 0.5 MeV 19.3; 2 MeV 17.3; 5 MeV 12.0; 10 MeV 8.8; 20 MeV 6.8; 100 MeV 4.86; 1 GeV 2.84; floor 2.5. Gakis and Atri's "Q > 20" for neutrons is an ICRP 60 style statement; ICRP 103 peaks at 20.7 near 1 MeV.

### 4.2 wT. Verdict: CONFIRMED (2ND-HAND)

[R3386] (Radiat. Environ. Biophys. 64, 637, PMC12701026, Table 1) and Wikipedia agree: bone marrow, colon, lung, stomach 0.12; breast 0.12; remainder 0.12 (adrenals, ET region, gall bladder, heart, kidneys, lymph nodes, muscle, oral mucosa, pancreas, prostate/uterus-cervix, small intestine, spleen, thymus); gonads 0.08; thyroid, oesophagus, liver, bladder 0.04; bone surface, skin, brain, salivary glands 0.01. COMPUTED sum 1.00. ICRP 60 values (gonads 0.20, breast 0.05) must not be mixed in.

### 4.3 "Subsurface dose dominated by secondary neutrons rather than primary protons". Verdict: PARTLY

| Source | Finding | Status |
|---|---|---|
| [R4890], JGR Planets 127, e2021JE007157 (open) | neutrons about 50% of effective dose at the Mars surface and above 80% below about 50 cm; neutron flux and effective dose peak near 30 cm depth (60 to 70 g/cm2 with atmosphere), neutrons "primarily induced by GCR protons"; but neutrons only 4, 8, 9% of absorbed dose in water at the surface (82, 753, 1200 Pa); 100 mSv/yr needs 1 to 1.6 m regolith | READ |
| Koehler et al., [R1728] (MSL RAD) | 7 to 740 MeV neutrons 5.1 uGy/d and 23.6 uSv/d vs totals 233 and 610: 2.2% of dose, 3.8% of dose equivalent; 1 to 1000 MeV extrapolation about 3% and 5% | READ |
| [R4888], Sci. Adv. 6, eaaz1334 (Chang'E-4 LND, PMC7518862) | Moon surface: neutral particles 3.1 +/- 0.5 of 13.2 +/- 1 uGy/h in Si (23 +/- 8%, includes gamma); charged dose equivalent 57.1 uSv/h; neutron contribution "expected to increase" under regolith; Apollo 17 neutron flux rises to about 150 g/cm2 | READ |
| [R0996], Biogeosciences Discuss. 4, 455 | with depth "muons and pions, and the electrons ... become increasingly dominant"; neutrons ionise only via recoil protons | READ |
| Gakis and Atri, [R1493] | water-rich regolith attenuates neutrons below 10 MeV; model weak on neutrons vs MSL RAD; surface shielding, not subsurface dose | READ |
| no source | H in water-rich regolith cuts neutrons below 10 MeV | 2ND-HAND (via [R4890]) |

Corrected statement: in the first metre of regolith, dose equivalent and effective dose are neutron-dominated (Mars: about 50% surface, above 80% below 50 cm); absorbed dose at the surface is charged-particle dominated (neutrons 2 to 9%); at metres to tens of metres muons, pions and cascade electrons dominate. Not supported without these qualifiers.

---

## 5. Dataset identities and licences (claim 5)

| Item | Identity | Terms found | Status | Verdict |
|---|---|---|---|---|
| "LAMDA / [R1873]" | [R1873] is the Leiden photodissociation and photoionisation database (102 species), home.strw.leidenuniv.nl/~ewine/photo/. LAMDA (home.strw.leidenuniv.nl/~moldata/) is the collisional line database, A&A 432, 369 | Photo site, exact text: "When using material from this site, please cite one of the overview papers. In cases of specific molecules, please cite the original source of the cross sections as appropriate." LAMDA: "If you use the data files in your work please refer to the publication by Schoier ... 2005, A&A 432, 369-379 ... Please mention the date when you accessed the database". No formal licence on either | READ | CONTRADICTED (conflated) |
| [R2040] | citation correct: Planet. Space Sci. 106, 11-45, [R2040]; Crossref funder NASA award NNX08AF03G; article Elsevier licences | data page phidrates.space.swri.edu connection reset; "public domain NASA legacy" not shown by anything reached | NOT REACHED (terms) | PARTLY |
| JPL 19-5 | JPL Publication 19-5, Evaluation No. 19 | jpldataeval.jpl.nasa.gov refused (502); not in NTRS or Crossref | NOT REACHED | NOT REACHED |
| Mainz atlas [R2303] | ESSD 5, 365; 8203 cross-section files, 1402 species, 154 quantum-yield files (site count) | live site: "please cite our publication: Keller-Rudek ... (2013)", no licence. Zenodo record 6951 ([R2302], frozen 2013-06-24, 126 MB): cc-by-4.0, open access. ESSD paper CC BY 3.0. Files derive from third-party papers (rights not individually cleared, my caution) | READ; DATA (Zenodo metadata) | terms found |
| Bremen / Serdyuchenko O3 | iup.uni-bremen.de .../o3spectra2011: 11 temperatures 193 to 293 K, 213 to 1100 nm, zips of 2012-07-25 and 2013-07-22; papers AMT 7, 625 and 609 (2014) | no licence or terms text on the data page (only a "Data Privacy" link); papers CC BY 3.0 (Crossref), not stated for the files | READ | no licence stated |

Rule for planetGen: do not vendor these tables without a licence line; read user-supplied copies or ask the data owners. Extends the licence rows of 02.

---

## 6. NH3 chain (claim 6)

### 6.1 "[R0847]". Verdict: PARTLY

No standalone [R0847] NH3 paper exists in Crossref. C2018 is the preliminary room-temperature CoYuTe, named in [R2098] (Icarus 321, 572, "Analysis of gaseous ammonia (NH3) absorption in the visible spectrum of Jupiter - Update", [R2097]; READ). Published list: [R4323], "ExoMol molecular line lists XXXV. A rotation-vibration line list for hot ammonia" ([R0847]; READ): 16.9 billion transitions, to 20 000 cm-1 (beyond 0.5 um), up to 1500 K (visible part only 300 K per [R2098]), improving BYTe (ends at 0.833 um). Authors caution: visible line positions off by several to tens of cm-1; absolute visible intensities not independently validated.

### 6.2 "Resolve anomalies near the 500 mbar condensation level". Verdict: UNSUPPORTED as stated

- [R2098] (READ): VLT/MUSE 0.48 to 0.93 um, NEMESIS correlated-k with CH4 from [R2265]; chi2/n 1.93 to 1.60 (Equatorial Zone) and 2.06 to 1.76 (NEB) switching from [R2777] to C2018; the 0.75 um feature is an artefact. No 500 mbar anomaly mentioned.
- [R2100] ([R2100], MNRAS; READ): CoYuTe below 0.977 um, HITRAN2022 above, CH4 [R2265] below 0.93 um and TheoReTS above; upper cloud at about 0.55 bar "coinciding with the ammonia condensation level", lower cloud at 1 to 2 bar at the start of the NH3 drop.
- [R0411] ([R0411]; READ): reflecting layer 570 +/- 30 mbar from 5.32 um NH3 lines (NH3 cloud); NH3 200 +/- 50 ppm at 0.7 to 5 bar (GRS).
- Not reached: Li 2017 Juno [R2097] first paper full text.

### 6.3 Visible NH3 absorption: what is established

| Fact | Value | Source | Status |
|---|---|---|---|
| Bands | 0.552, 0.648, 0.765, 0.79, 0.825, 0.93 um | [R2098] | READ |
| Line list (6475 A band) | 33 strongest lines, 36 m, 1 atm, 294 K, about 17% of band absorption, no lower-state energies | [R2098] | 2ND-HAND |
| [R2777] | 6475 and 5520 A cross sections, room T, self-broadened; resolution about 10 A not 2 A; strengths inconsistent with longer-wavelength data | [R2098] | 2ND-HAND |
| Numeric band strengths 500 to 800 nm | not in any reachable text (ExoMol file contents needed) | | NOT REACHED |
| Widths used | H2 to J = 9 then 0.0788 cm-1/atm; He 0.0282; T exponents 0.59 (H2), 0.37 (He) | [R2098] | READ |
| "[R4186], Jupiter's NH3 absorption" | not found; real: [R4186], "The source of 3-um absorption in Jupiter's clouds: Reanalysis of ISO observations using new NH3 absorption models" (infrared) | Crossref | READ (title) |

Reachable data: [R2258] cold methane absorption coefficients, 300 to 1000 nm, 1750 rows, at https://atmos.nmsu.edu/planetary_datasets/infraredsaturn.txt (DATA; terms not stated; PDS Atmospheres node). Sample (vacuum nm, 1/km-am): 500 0.0014; 540 0.0718; 600 0.022; 620 0.683; 650 0.0379; 700 0.244; 726.8 4.10; 790 0.550; 840 0.602; 860 2.51; 890 24.1; 900 12.7; 990 12.3.

---

## 7. Final table: claim, verdict, corrected number

| # | Claim | Verdict | Corrected or supported value |
|---|---|---|---|
| 1a | omega0 0.92 to 0.95 at 670 nm | PARTLY | solar mean 0.89 to 0.94; about 0.975 near 650 nm (2ND-HAND); UV band 0.67 to 0.74 |
| 1b | k > 0.01 below 350 nm | CONFIRMED | 0.013 (0.30 um), 0.011 (0.35), 0.0093 (0.40); min about 1e-3 at 0.7 to 0.8 um (digitised) |
| 1c | Fe3+ attribution | NOT REACHED | no source states it |
| 1d | SO2 row cites [R2645, R4395] | CONTRADICTED | both Mars dust imaging papers |
| 2a | Mie underpredicts forward scattering 40% | UNSUPPORTED | no source |
| 2b | haze extinction law (our gap) | CONFIRMED | 1.012e7 L^-2.339, 2.029e4 L^-1.409, 6.27e2 L^-0.9706 (nm); total 9.6 at 500 nm |
| 3a | CO 0.15 to 30% | PARTLY | 0.3 to 26, median 3.3 |
| 3b | CO2 2 to 30% | CONFIRMED | 2 to 30, median 12 |
| 3c | CH4 0.2 to 1.5% | PARTLY | 0.15 to 2.7 |
| 3d | NH3 0.2 to 1.0% | CONTRADICTED | 0.1 to 3.6; mean 0.8 |
| 3e | Q(H2O) 1e26 to 1e31 | PARTLY | 5.1e26 to 1.41e31 (SWAN sample) |
| 4a | wR values and neutron function | CONFIRMED | formula in 4.1 |
| 4b | wT values | CONFIRMED | sum 1.00 |
| 4c | subsurface neutron dominance | PARTLY | dose equivalent only, 50% surface to above 80% below 50 cm (Mars); absorbed dose 2 to 9% |
| 5a | "LAMDA / [R1873]" | CONTRADICTED | [R1873] = Leiden photo database; LAMDA = a molecular database; citation request only |
| 5b | [R2040] public domain | PARTLY | citation right; terms NOT REACHED |
| 5c | JPL 19-5 terms | NOT REACHED | |
| 5d | Mainz terms | terms found | cite [R2302, R2303]; Zenodo CC BY 4.0 |
| 5e | Bremen terms | none stated | |
| 6a | [R0847] | PARTLY | [R4323] |
| 6b | 500 mbar anomaly resolved | UNSUPPORTED | chi2/n 1.93 to 1.60 and 2.06 to 1.76; NH3 cloud 0.55 to 0.57 bar |
| 6c | K and T 2009/2010/2011 titles | PARTLY | only 2010 is the CH4 paper |

Counts (23 rows): CONFIRMED 5; PARTLY 8; CONTRADICTED 3 (1d, 3d, 5a); UNSUPPORTED 2 (2a, 6b); NOT REACHED 2 (1c, 5c); terms-found or none-stated 3 (5d, 5e plus the 5b terms).

---

## 8. NOT REACHED (Boss must research)

1. [R4743] tables of omega0(lambda) and (n, k)(lambda); [R0809, R3326, R3567] tables.
2. [R4394] (PSS 56, 669) Tables 2 and 3 and phase functions; [R1144] coefficients; [R2584, R3663] full text. Any source for a Mie versus aggregate forward-scattering percentage.
3. [R0048, R0436, R3206] Table 1, Dello [R1063] original table; Q(H2O) of small comets.
4. ICRP 103 body tables and ICRP 116; Chang'E-4 neutron-only dose equivalent.
5. Terms: Huebner and Mukherjee data page, JPL 19-5, Bremen files, Leiden photo database (citation request only), NMSU Karkoschka files.
6. NH3 numeric band strengths 500 to 800 nm ([R2777], CoYuTe contents); [R2097] first paper; Li 2017; [R2265] supplementary table.

## 9. Recalled / to verify

JPL Publication 19-5 title and authors; L = 51 kJ/mol (from 05); Fe3+ as the UV absorber; that the ICRP 103 body table matches the Wikipedia reproduction (extract confirms narrative only).

## 10. Cross-links

02 (optical constants, dataset licences), 05 section 9 (differences in 3.1), 06 (dust, hazes, sulphate; the misplaced Lemmon/Tomasko row is Mars dust), 07 (radiation weighting, shielding, subsurface neutrons).
