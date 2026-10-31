# NH3 photodissociation rate at 1 AU: why 5.6e3 s and 6.7e3 s differ, and what to use

Scope: one number for the generator, J(NH3, 1 AU) for a quiet and an active Sun, with its spread. Labels as in _worker-rules.md. All files were read as untrusted data (own directories, python3 -I). Rate J in 1/s, lifetime tau = 1/J, scaled to distance r (AU) by r^-2.

## 1. Answer in four lines

1. Neither value is wrong. They are 17 percent apart (1.54e-4 vs 1.8e-4 per s), inside the cross-section and solar-flux uncertainty.
2. The gap is NOT the 200 nm cut. The Leiden "Solar" NH3 rate (4.08e-9 x 37700 = 1.538e-4 per s) uses the Leiden solar file, which runs 67 to 420 nm, with the NH3 cross section to 226 nm. I reproduced 1.5378e-4 to four digits only when integrating to 226 nm; the part below 200 nm is 7.35e-5, and the 200 to 226 nm part (8.0e-5) is half of the total. A 200 nm cut would give tau about 1.4e4 s, a factor 2 too long. The statement "covers only 91.2 to 200 nm" applies to the interstellar and normalisation fields, not to the NH3 solar rate.
3. The gap is solar state plus input data: Leiden uses a quiet-Sun spectrum (UARS SOLSTICE Woods 1996 above) and the cross section, which is my highest cross section (+11 percent vs Cheng 2006). Huebner et al. 1992 (H92) and Huebner and Mukherjee 2015 (HM15) inputs were NOT REACHED.
4. Recomputed from modern measured inputs: quiet 1.28 to 1.54e-4, active 1.35 to 1.70e-4 per s. H92/HM15 1.8e-4 sits at the upper edge. Recommendation in section 5.

## 2. Literature values (all reached values)

| Source | Value for NH3 at 1 AU | Detail | Status |
|---|---|---|---|
| Drozdovskaya et al. 2023, arXiv:2307.11486, footnote 9 | 1.8e-4 per s (5.6e3 s) | "average rate for a quiet and an active Sun from Huebner et al. (1992) and Huebner and Mukherjee (2015)", for the NH2(X) + H channel only | READ (2ND-HAND for H92/HM15) |
| same footnote | 1.5e-4 per s (6.7e3 s) | "value obtained in Heays et al. (2017) for the Solar radiation field at 1 au" | READ |
| Biver et al. 2012, arXiv:1201.4318 | 1.8e-4 per s | "for the quiet Sun at rh = 1 AU (Huebner et al. 1992)" | 2ND-HAND |
| Kawakita and Mumma 2011, ApJ 727, 91 (IOPscience full text) | NH3 1.8e-4, NH2 2.15e-6 per s | cites Huebner et al. 1992, no activity state given | 2ND-HAND |
| arXiv:0904.0592 (Tempel 1, Deep Impact) | NH2 rate 2.2e-6 (quiet) to 3.4e-6 (active) per s from H92; "literature" NH3 rate 6.6e-5 to 2.1e-4 per s at 1.51 AU | NH3 rate source not named; times r^2 gives 1.5e-4 to 4.8e-4 at 1 AU | 2ND-HAND |
| Hatchell et al. 2005, astro-ph/0505439 | 5400 to 5800 s (1.72 to 1.85e-4) | source not named; "probably somewhat longer than average" at solar minimum | READ |
| Biver et al. 2007 (Odin), astro-ph/0610779 | 6700 s (1.49e-4) | no source named | READ |
| no source abstract | 7e3 s | value used in a vectorial-model fit of NH2 profiles, not an independent measurement | READ (abstract) |
| Heays et al. 2017, arXiv:1701.04459, Table 22 "Solar" | NH2 + H 4.0e-9 (99 percent), NH + 2H/H2 5.6e-11, per s in the 37700-scaled field | x 37700 gives 1.51e-4 and 2.1e-6 per s | READ, COMPUTED |
| Leiden photo_rates.html, "Solar" column | dissociation 4.08e-9, ionisation 9.57e-12 | x 37700: 1.538e-4 (tau 6.5e3 s) and 3.6e-7 per s | READ, COMPUTED |
| sbpy 0.6.0, activity/gas/data/__init__.py | no NH3 entry | holds H2O, OH, HCN, CH3OH, H2CO, CO, CO2, CN only | DATA |
| Huebner and Mukherjee 2015, Crovisier 1994 | NH3 value not found | NTRS record 23332861698053 is metadata only (downloadsAvailable false); phidrates.space.swri.edu unreachable | NOT REACHED |

Reading the Huebner values: Biver 2012 and Kawakita and Mumma pair 1.8e-4 with an NH2 rate of 2.15e-6, which is the quiet end of the H92 NH2 range (2.2e-6 quiet), so 1.8e-4 is probably the H92 quiet-Sun NH3 rate, not an average (INFERENCE from 2ND-HAND). Drozdovskaya calls it an average. If 1.8e-4 is the quiet value, the H92 active value is higher still. This cannot be settled without the H92 or HM15 table.

## 3. My computation (COMPUTED)

Formula: J = integral of sigma_diss(lambda) F(lambda) d lambda, F in photons per cm2 per s per nm = E x lambda / (h c), E from the irradiance file at 1 AU. Dissociation quantum yield 1 at all wavelengths from 100 to 230 nm: Heays 2017 section 4.3.46 assumes 100 percent dissociation below the ionisation limit (READ), and the photochem_clima_data NH3.h5 channel yields (cited Huebner 2015) sum to 1.000 at every point (DATA). Ionisation above 123.1 nm is zero; below it, the Leiden dissociation part is used.

Cross sections (cm2, 0.1 nm bins): 100 to 123.5 nm Leiden photodissociation (photochem_clima_data NH3.h5 from raw.githubusercontent.com/Nicholaswogan/photochem_clima_data, Heays-derived, GPL-3 package); 123.5 to 140 nm Wu et al. 2007 (Mainz file); 140 to 230 nm one of: Cheng 2006 (primary, 0.02 nm, 298 K), or the Leiden file (per Heays 2017 Table 12). All Mainz files from uv-vis-spectral-atlas-mainz.org/uvvis_data/cross_sections/ (DATA). The A-band integral over 170 to 230 nm (1e-17 cm2 nm) is: Suto-Lee 22.1, Chen 23.9, Watanabe 24.7, Limao-Vieira 24.9, Cheng 25.1, Leiden 27.0. Peak value about 2.0e-17 cm2 at 193.8 nm in all.

Solar irradiance at 1 AU (DATA, sbuv.gsfc.nasa.gov/solar, files by M. DeLand; units mW/m2/nm for the h5 files, checked against the W/m2/nm WHI file, ratio 1.01 to 1.06): SORCE SOLSTICE v16 1 nm (120 to 300 nm), UARS SOLSTICE v18 and SUSIM v22 1 nm, daily medians of FillFlag 0 days; WHI 2008 reference spectra of Woods et al. 2009 (0.1 nm, TIMED SEE below 116 nm); ATLAS-1 Thuillier 2004. Below 120 nm the WHI quiet spectrum is used, scaled by the Lyman-alpha ratio to the file (1.0, 1.35, 1.48). Rebinning check: WHI at 1 nm gives 1.329e-4 vs 1.321e-4 at 0.1 nm.

Result with the Cheng 2006 cross section (units 1e-4 per s):

| Solar spectrum (Lyman-alpha mW/m2/nm) | Total | Up to 200 nm | 200 to 217 nm | 217 to 230 nm | tau (s) |
|---|---|---|---|---|---|
| Quiet: SORCE 2008-09 minimum (6.03) | 1.389 | 0.697 | 0.673 | 0.019 | 7200 |
| Quiet: SUSIM 1996 minimum (6.56) | 1.373 | 0.690 | 0.662 | 0.021 | 7280 |
| Quiet: WHI 2008 Apr 10-16, F10.7 68.9 | 1.324 | 0.672 | 0.632 | 0.019 | 7550 |
| Moderate: ATLAS-1 1992 | 1.509 | 0.767 | 0.720 | 0.022 | 6630 |
| Active: SORCE 2014 maximum (8.12) | 1.460 | 0.746 | 0.694 | 0.020 | 6850 |
| Active: SUSIM 2001-02 maximum (9.07) | 1.468 | 0.753 | 0.693 | 0.022 | 6810 |
| Active: UARS SOLSTICE 2000-01 maximum (8.91) | 1.529 | 0.769 | 0.738 | 0.022 | 6540 |
| Leiden file x 37700, Leiden cross section | 1.538 | 0.735 | 0.802 (200 to 230) | | 6500 |

Cross-section spread (total, 1e-4 per s; SORCE minimum | UARS SOLSTICE maximum): Cheng 1.389 | 1.529; Chen 1.334 | 1.470; Suto-Lee 1.280 | 1.410; Watanabe 1.374 | 1.511; Limao-Vieira 1.381 | 1.520; Leiden 1.544 | 1.698.

Where the rate comes from (Cheng, quiet SORCE | active UARS; 1e-6 per s): 100 to 123.5 nm 4.3 | 6.3; 123.5 to 170 nm 1.0 | 1.3; 170 to 180 3.2 | 3.7; 180 to 190 18.9 | 19.7; 190 to 200 42.2 | 45.8; 200 to 210 50.8 | 56.0; 210 to 217 16.5 | 17.8; 217 to 230 1.9 | 2.2. So 95 percent of J is from 170 to 217 nm, and Lyman-alpha supplies only 2 to 3 percent. Solar photon flux 170 to 217 nm at 1 AU: 4.46e13 (quiet), 4.89e13 (active) per cm2 per s; the active-to-quiet ratio is 1.04 to 1.15 per band there but 1.3 to 1.5 below 170 nm.

Other effects: gas temperature 175 K instead of 295 K (truncated at 218 nm) lowers J by 4 to 6 percent (1.278e-4 vs 1.334e-4). NH3 photoionisation adds 3.6e-7 (0.25 percent). Water check of the method (Leiden H2O file, 100 to 230 nm only, SORCE minimum): 1.07e-5 per s vs HM15 quiet 1.21e-5 (tau 8.29e4 s from sbpy), consistent given the missing EUV. For H2O the active-to-quiet ratio is 1.33 (mine) vs 1.68 (HM15); the HM15 "active Sun" is more extreme than any cycle-maximum spectrum used here, which hints that H92 and HM15 also differ from modern spectra in the UV.

## 4. How much of the gap is what

| Step (Cheng unless stated) | J (1e-4 per s) | Change |
|---|---|---|
| Leiden published (quiet, Liang cross section, to 226 nm) | 1.54 | start |
| Same spectrum cut at 200 nm | 0.74 | factor 0.48 (not what the Leiden number does) |
| Swap Liang/Leiden for Cheng cross section, SORCE quiet | 1.39 | -10 percent |
| Quiet to cycle-maximum, Cheng | 1.46 to 1.53 | +5 to +10 percent |
| Quiet to maximum, Leiden cross section | 1.54 to 1.70 | +5 to +10 percent |
| Remaining to H92/HM15 1.8 | 0.1 to 0.5 | NOT explained: needs H92 flux and cross section |

Wavelength coverage explains none of the 6.7e3 vs 5.6e3 gap. Solar state explains at most 5 to 10 percent. The cross-section choice explains up to 10 to 15 percent. The rest (5 to 25 percent) is the unknown H92/HM15 inputs (older irradiance scale, older cross sections, an "active" Sun that is stronger than modern maxima).

## 5. Recommendation for the generator

| Quantity | J at 1 AU (per s) | tau (s) | Range |
|---|---|---|---|
| Quiet Sun | 1.4e-4 | 7.1e3 | 1.3e-4 to 1.55e-4 (6.5e3 to 7.7e3 s), cross-section spread |
| Mid or default | 1.5e-4 | 6.7e3 | 1.3e-4 to 1.8e-4 (5.6e3 to 7.7e3 s) |
| Active Sun | 1.55e-4 | 6.5e3 | 1.35e-4 to 1.7e-4 (5.9e3 to 7.4e3 s); 1.8e-4 (5.6e3 s) as literature upper bound |

Reasoning: the Leiden value is a transparent, reproducible integral of measured quantities and my independent recomputation reproduces it, so 1.5e-4 is the anchor. The H92/HM15 value cannot be reproduced and one 2ND-HAND chain labels it quiet, another an average, so it is kept as the upper edge, not the centre. NH3 is a near-UV (170 to 217 nm) absorber, so solar activity changes J by only 5 to 10 percent; draw the quiet/active choice as a +/-5 percent term and the cross-section/flux scale as a +/-12 percent term. The old note's 2.5e4 s is still far too long (factor 3.5 to 4.5).

Code recipe: J(r, star) = J_sun_1AU x [F_star(170 to 217 nm at planet) / 4.5e13 photons per cm2 per s] x 1.0 + about 5 percent from shorter wavelengths scaled by the star's FUV (Lyman-alpha and 100 to 170 nm); for the Sun J_sun_1AU = 1.4e-4 to 1.55e-4. Validity: optically thin, T 100 to 300 K; the cross-section peak 2e-17 cm2 gives tau = 1 at a column near 5e16 cm-2, so self-shielding matters in thick NH3 columns. Products: NH2 + H about 99 percent in the solar field (Heays 2017 Table 22).

## 6. NOT REACHED / Boss must research

- Huebner, Keady and Lyon 1992 (Ap&SS 195, 1) NH3 quiet and active rates and their solar and cross-section inputs. Decides whether 1.8e-4 is quiet or average.
- Huebner and Mukherjee 2015 (Planet. Space Sci. 106, 11) NH3 channel table (Elsevier refused; NTRS metadata only; phidrates.space.swri.edu unreachable).
- Crovisier 1994 (JGR 99, 3777) NH3 value (Wiley/AGU refused; the 2013 survey arXiv:1308.6282 quotes it for other species only).
- The NH3 cross-section file itself (Mainz holds only the 15NH3 file) and the Leiden NH3 file as such (I used the same data via photochem_clima_data).
- Source of Hatchell 5400 to 5800 s, Odin 6700 s, and the 6.6e-5 to 2.1e-4 range in arXiv:0904.0592.
- Temperature-dependent NH3 cross sections below 175 K and above 295 K for the 200 to 217 nm tail (the largest part of J).

## 7. Cross-links

Part 04 (photochemistry lifetimes and escape) and the earlier note C section 3.1 (NH3 row: verdict CONTRADICTED stays, but its explanation "91.2 to 200 nm only" should be corrected as above); Part 02 (NH3 cross-section tables, Leiden file).
