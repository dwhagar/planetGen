# Part 06: condensates, clouds, hazes, dust and minerals

Research note for the planetGen atmosphere and body composition effort (v9.0 headline). Read-only worker output; nothing in the repository was changed. Solar System bodies appear only as labelled calibration or test points, never as presets.

Verification labels: READ (original paper text), 2ND-HAND (read in a paper or code that cites the original), DATA (data file in a package or database), COMPUTED (re-derived here with Python from stated inputs), SNIPPET, RECALLED, NOT REACHED. Where two sources disagree both are shown. Approximate table-row counts by primary label: DATA 85, READ 66, COMPUTED 55, 2ND-HAND 28, RECALLED 15 (section 5), SNIPPET 1, NOT REACHED 3 (plus the 11 items of section 4).

## 0. Sources reached and refused (for the next worker)

Reached and used: export.arxiv.org (PDF text), raw.githubusercontent.com (LX-MIE data), pypi.org and files.pythonhosted.org (virga-exo 2.0.2, miepython, chemicals 1.5.2, thermo 0.6.1), zenodo.org (virga optical data, record 15886530, CC-BY-4.0), webbook.nist.gov (phase-change pages), webmineral.com, the HITRAN refractive-index mirror (HITRAN-RI ascii files; per-file references exist, overall licence not read), refractiveindex.info (YAML, abbreviated RII), atmos.washington.edu, NTRS (abstract records only), WebSearch. Refused or unusable: doi.org, science.org, sciencedirect, agupubs/wiley, aanda, oup, cambridge, acs, researchgate, github.com web pages, usgs, mindat, lpi.usra.edu; the arXiv API rate-limited. NIST WebBook pages for SiO, TiO, H2SO4, Fe, CaCO3, Cr2O3, S8, N2O4 and HNO3 gave no phase-change table.

## 1. Scope and how the numbers feed the generator

Input to this part (from parts 01 to 05 and the orbit/star model): a pressure-temperature profile T(P), gravity g, metallicity [Fe/H] (or [M/H]), partial pressures or mixing ratios of the condensable gases, stellar UV flux, surface type, wind and surface-grain properties.
Output from this part, per body: (1) which condensates exist, their base pressure (T-P profile crossing the condensation curve) and the gas removed by condensation (feeds back into part 05 mixing ratios); (2) cloud and haze optical depth, particle radius, density, complex refractive index at 450, 550, 650 nm, single scattering albedo and Angstrom exponent (feeds atmosphere colour GEN.414 and radiation parts); (3) dust lifting threshold, airborne dust radius and column optical depth; (4) the surface or crust mineral list and vapour species for lava worlds.
Method order for a generator: (a) compute T_cond(P, [Fe/H]) for each species from 2.1 to 2.5 using the rules in 3.1 to 3.4; (b) sort by T_cond at the profile, find base pressures, deplete the limiting element in order (relative masses at the end of 2.1); (c) set particle size from the Ackerman and Marley sedimentation rule (3.5, 3.6) with a drawn f_sed; (d) take optical constants from 2.8 and compute optical depth from column mass and radius; (e) add photochemical haze using the production rule in 3.7; (f) add dust from 3.8 and 3.9; (g) convert optical depth to cloud albedo with 3.10.

Coverage of the compounds named by the project text (where each condensed phase sits; label of the main data):

| Compound | Condensed form or role | Where | Status |
|---|---|---|---|
| H2, He | no condensate at planetary temperatures (H2 triple point 13.95 K, 0.0721 bar; He T_c 5.2 K) | 2.5 | DATA |
| N2, O2, Ar, Ne, Kr, Xe, CO, CH4, ethane, C2H2, C2H4, HCN, SO2, H2S | ices | 2.4 | 2ND-HAND |
| O3, radon, HCl, HF, CH3Cl, formaldehyde, NO2, OCS | liquid Antoine sets and anchors (polymers, hydrates NOT REACHED) | 2.5 | DATA |
| CO2 | ice curve and optical constants | 2.4, 2.8 | 2ND-HAND, READ |
| H2O | liquid and ice; ice-salt eutectics NOT REACHED | 2.4, 2.5 | 2ND-HAND, DATA |
| NH3, NH4SH, NH4Cl | ice and NH4SH fit; NH4Cl enthalpy only | 2.1, 2.4, 2.5 | READ, 2ND-HAND |
| sulfides, H2SO4 | MnS, Na2S, ZnS (FeS suppressed); H2SO4-H2O vapour curve NOT REACHED | 2.1, 2.3, 2.8 | READ, NOT REACHED |
| Na, K | Na2S, KCl, metal vapour sets | 2.1, 2.5 | READ, DATA |
| mineral dust, carbonates | densities and optics; decomposition equilibria NOT REACHED | 2.8, 2.12 | DATA |
| perchlorates, NOx hydrates, HNO3 | NOT REACHED (HNO3 optics files exist, unread) | section 4 | NOT REACHED |
| tholins, soot, haze | production, size, optics | 2.8, 2.9 | DATA, READ |
| SiO, Mg, Fe, TiO, VO | lava-world vapour | 2.11 | READ for SiO, Mg, Fe; NOT REACHED for TiO and VO |

## 2. Tables

### 2.1 Condensation temperature fits (T_cond in K from 10^4/T_cond; P_T total pressure in bar; solar-scaled by [Fe/H])

All fits are for condensation from an undepleted gas. Each is valid only where that phase is the first (highest-temperature) member of its group, and only for solar-like element ratios scaled together (Wakeford 2017 and Visscher 2010 state this). Fits come from the cited papers; the numerical T_cond values are tabulated in 2.2.

| Condensate | 10^4 / T_cond = | Source | Status |
|---|---|---|---|
| Fe | 5.44 - 0.48 log P_T - 0.48 [Fe/H] | Visscher 2010, arXiv:1001.3639 | READ |
| Mg2SiO4 (forsterite) | 5.89 - 0.37 log P_T - 0.73 [Fe/H] | Visscher 2010 | READ |
| MgSiO3 (enstatite) | 6.26 - 0.35 log P_T - 0.70 [Fe/H] | Visscher 2010 | READ |
| Cr | 6.576 - 0.486 log P_T - 0.486 [Fe/H] | Morley 2012, arXiv:1206.4313 | READ |
| MnS | 7.447 - 0.42 log P_T - 0.84 [Fe/H] (Visscher 2006: 7.45 - 0.42 (log P_T + [Mn/H] + [S/H])) | Morley 2012; Visscher 2006, astro-ph/0511136 | READ |
| Na2S | 10.045 - 0.72 log P_T - 1.08 [Fe/H] (Visscher 2006: 10.05 - 0.72 (log P_T + [Na/H] + 0.5 [S/H])) | Morley 2012; Visscher 2006 | READ |
| ZnS | 12.527 - 0.63 log P_T - 1.26 [Fe/H] (Visscher 2006: 12.52 - 0.63 (log P_T + [Zn/H] + [S/H])) | Morley 2012; Visscher 2006 | READ |
| NH4SH | 48.91 - 4.15 log P_T - 4.15 [M/H] | Visscher 2006 eq. 32 | READ |
| H2S ice | 86.49 - 8.54 log P_T - 8.54 [S/H] | Visscher 2006 eq. 34 | READ |
| NH4H2PO4 | 29.99 - 0.20 (11 log P_T + 15 [M/H]) | Visscher 2006 eq. 52 | READ |
| CaAl4O7 (grossite) | 4.990 - 0.2394 log P_T + 1.398e-3 (log P_T)^2 - 0.595 [Fe/H] | Wakeford 2017, arXiv:1610.03325 eq. 2 | READ |
| CaAl12O19 (hibonite) | 4.997 - 0.2262 log P_T + 1.888e-3 (log P_T)^2 - 0.585 [Fe/H] | Wakeford 2017 eq. 3 | READ |
| Al2O3 (corundum) | 5.014 - 0.2179 log P_T + 2.264e-3 (log P_T)^2 - 0.580 [Fe/H] | Wakeford 2017 eq. 4 | READ |
| Ca3Ti2O7 | 5.114 - 0.284 log P_T - 0.568 [Fe/H] | Wakeford 2017 eq. 5 | READ |
| Ca4Ti3O10 | 5.113 - 0.281 log P_T - 0.563 [Fe/H] | Wakeford 2017 eq. 6 | READ |
| CaTiO3 (perovskite) | 5.125 - 0.277 log P_T - 0.554 [Fe/H] | Wakeford 2017 eq. 7 | READ |

Further READ anchor points: Visscher 2010 gives quartz T_cond about 1550 K but SiO2 cannot form in the presence of enstatite (which condenses at about 1600 K in the same comparison); MgO condenses only above 10^3.5 bar; forsterite melts at 2163 K and enstatite at 1830 K; the Fe/Fe(OH)2 boundary is log P_T = 11.94 - 12088/T - 2 [Fe/H]. Visscher 2006 states that FeS does not form in substellar atmospheres because Fe metal clouds have removed the iron at depth (Fe/S about 1.9 in the gas), so sulfur condenses as MnS, Na2S, ZnS and then NH4SH. Wakeford 2017 Table 1 gives initial condensation temperatures at 0.1 bar, solar, after Lodders 2002: hibonite 1914, Ca4Ti3O10 1854, grossite 1848, Ca3Ti2O7 1818, gehlenite 1796, perovskite 1707, Fe 1682, spinel 1663, akermanite 1626, Ti3O5 and Ti2O3 1610, forsterite 1592, diopside 1560, anorthite 1513, enstatite 1508, Cr 1416, Ti4O7 1290, MnS 1271, Na2S 925 K (READ; the paper marks grossite, Ca3Ti2O7, perovskite and several others as first appearances in the full sequence after removal of higher-temperature phases at 0.1 bar, not as initial condensation temperatures, and states that corundum does not appear in the sequence at 0.1 bar, which explains why the fit curves differ from those entries).

Relative cloud masses (limiting-element mass if all the element condenses; Na2S = 1000 in Wakeford; Morley 2012 uses Na2S = 1.00 with a different normalisation): Wakeford 2017 Table 1: hibonite 2091, grossite 2442, corundum 1915, Ca4Ti3O10 170, Ca3Ti2O7 180, perovskite 149, gehlenite 3678, Fe 21032, spinel 2673, akermanite 3657, forsterite 32180, diopside 5809, anorthite 5227, enstatite 44585, Cr 303, MnS 360 (READ). Morley 2012 Table 3 (Na2S = 1.00): Cr 0.30, MnS 0.36, Na2S 1.00, ZnS 0.05, KCl 0.12, Fe 20.85 (READ). Corundum mass scales with metallicity: 1915, 3822, 6057, 19155 for solar, [M/H] 0.3, 0.5, 1.0 (READ).

### 2.2 Computed condensation temperatures (K) at [Fe/H] = 0 (and +0.5 where shown)

All entries COMPUTED from the 2.1 fits (KCl from the Gao/Morley expression with gas pressure 2.55e-7 P_T) at total pressure P_T (bar). Valid only where the phase is the first of its group; the two Mg silicates are both listed although enstatite only forms below the forsterite curve at pressures where forsterite is already gone.

| Condensate | [Fe/H] | 1e-4 | 1e-3 | 1e-2 | 0.1 | 1 | 10 | 100 |
|---|---|---|---|---|---|---|---|---|
| CaAl4O7 | 0 | 1675 | 1748 | 1827 | 1912 | 2004 | 2104 | 2214 |
| CaAl12O19 | 0 | 1686 | 1757 | 1833 | 1914 | 2001 | 2095 | 2197 |
| Al2O3 | 0 | 1689 | 1758 | 1832 | 1911 | 1994 | 2084 | 2180 |
| Al2O3 | +0.5 | 1776 | 1853 | 1935 | 2023 | 2117 | 2218 | 2327 |
| Ca3Ti2O7 | 0 | 1600 | 1676 | 1760 | 1853 | 1955 | 2070 | 2200 |
| Ca4Ti3O10 | 0 | 1603 | 1679 | 1762 | 1854 | 1956 | 2070 | 2197 |
| CaTiO3 | 0 | 1604 | 1679 | 1761 | 1851 | 1951 | 2063 | 2188 |
| Fe | 0 | 1359 | 1453 | 1562 | 1689 | 1838 | 2016 | 2232 |
| Fe | +0.5 | 1404 | 1506 | 1623 | 1761 | 1923 | 2119 | 2358 |
| Mg2SiO4 | 0 | 1357 | 1429 | 1508 | 1597 | 1698 | 1812 | 1942 |
| MgSiO3 | 0 | 1305 | 1368 | 1437 | 1513 | 1597 | 1692 | 1799 |
| Cr | 0 | 1174 | 1245 | 1325 | 1416 | 1521 | 1642 | 1784 |
| MnS | 0 | 1096 | 1149 | 1207 | 1271 | 1343 | 1423 | 1514 |
| Na2S | 0 | 774 | 819 | 871 | 929 | 996 | 1072 | 1162 |
| Na2S | +0.5 | 807 | 857 | 914 | 978 | 1052 | 1138 | 1240 |
| ZnS | 0 | 665 | 694 | 725 | 760 | 798 | 841 | 888 |
| KCl (see 2.3) | 0 | 625 | 662 | 702 | 749 | 801 | 862 | 933 |

Checks: the 0.1 bar solar values of Fe (1689), hibonite (1914), Ca4Ti3O10 (1854), Cr (1416), MnS (1271) and Na2S (929) agree with Wakeford Table 1 (1682, 1914, 1854, 1416, 1271, 925) within 7 K. For MnS, ZnS and MgSiO3 the ratio of the [Fe/H] coefficient to the log P_T coefficient is 2, so +0.5 dex of metallicity equals a tenfold pressure increase (visible in the table: MnS +0.5 at 0.1 bar equals solar at 1 bar, 1343 K).

### 2.3 Effective saturation expressions and microphysics inputs used by Gao 2020 and virga

Gao 2020 (arXiv:2005.11939, supplementary tables 2, 3, 4) and virga-exo 2.0.2 (`virga/pvaps.py`, `gas_properties.py`) give log10 P_sat in bar versus T. These expressions for the sulfides, Cr, Mg2SiO4, Fe, TiO2 and Al2O3 are effective equilibrium expressions for the net reaction with the metal gas as the controlling species (the metallicity term carries the H2S or oxygen scaling), not pure vapour pressures of the compound. KCl illustrates this: at 1400 K the Morley/Gao expression gives 0.30 bar, while the NIST Antoine sets for KCl vapour give 0.09 bar (COMPUTED from the two NIST sets in 2.5); extrapolating the NIST sets to the dilute gas pressure 2.55e-7 bar (outside their fit range) gives about 793 K against 801 K from the Gao/Morley expression (COMPUTED, one data point, extrapolated).

| Species | log10 P_sat (bar) | Surface energy (erg/cm2) | Density (g/cm3) | Latent heat (1e10 erg/g) | Contact angle on TiO2 at 800 K (deg) | Status |
|---|---|---|---|---|---|---|
| KCl | 7.611 - 11382/T | 160.4 - 0.07 (T - 273.15) | 1.988 | 2.923 | none (acts as nucleus for ZnS) | READ |
| ZnS | 12.812 - 15873/T - [Fe/H] | 860 | 4.04 | 3.118 | 81.7 (on KCl) | READ |
| Na2S | 8.550 - 13889/T - 0.5 [Fe/H] | 1033 | 1.856 | 11.57 | 61 | READ |
| MnS | 11.532 - 23810/T - [Fe/H] | 2326 | 4.0 | 8.297 | 77.6 | READ |
| Cr | 7.490 - 20592/T | 1642 - 0.2 (T - 2133.15) | 7.15 | 7.582 | 74.8 | READ |
| Mg2SiO4 | 4.88 (printed) or 14.88 (virga) - 32488/T - 1.4 [Fe/H] - 0.2 log P | 436 | 3.21 | 25.59 | 0.1 (imposed minimum) | READ and DATA, discrepancy below |
| Fe | 7.23 - 20995/T (Gao); 7.09 - 20995/T (virga) | 1862 - 0.39 (T - 1803.15) | 7.87 | 7.197 | 77.2 | READ and DATA |
| TiO2 | 9.5489 - 32457/T (Gao); 13.95 - 38266/T - [Fe/H] (virga) | 535.124 - 0.04396 T | 4.25 | 7.78 | none (homogeneous nucleus, about 2000 K) | READ and DATA |
| Al2O3 | 17.7 - 45893/T - 1.66 [Fe/H] (Gao); 15.24 - 41481/T - 1.66 [Fe/H] (virga, diamondback) | 690 | 3.99 | 32.56 | 43.6 | READ and DATA |

The Mg2SiO4 constant: with a Mg gas partial pressure of 3.55e-5 P (Mg/H from A(Mg) = 7.55, RECALLED; see section 5) the 14.88 form gives T_cond = 1681 K at 1 bar against 1698 K from the Visscher fit (COMPUTED), while 4.88 gives 3482 K, so the Gao Table 2 constant is a misprint of 14.88 (COMPUTED consistency check; confirm with the authors if used).
Other virga constants (DATA): MgSiO3 13.43 - 28665/T - [Fe/H]; CaTiO3 30.24 - 72160/T - log10 P - 2 [Fe/H]; CaAl12O19 16.46 - 44021/T - 0.083 log10 P - 1.67 [Fe/H]; SiO2 13.168 - 28265/T (the code comment says "MADE UP", do not use); NH3 6.900 - 1588/T below 195.4 K and 5.201 - 1248/T above (check: 0.0594 bar at 195.4 K against NIST triple point 0.0606 bar, COMPUTED); H2O; CH4 by Clausius-Clapeyron from its triple point (90.68 K, 0.11719 bar).
Lower boundary mixing ratios used by Gao 2020 (ppm as printed, mass or mole basis not stated in the extracted text; solar, multiplied by 10 for 10x solar): ZnS 0.076, Na2S 3.34, MnS 0.541, Cr 0.887, Mg2SiO4 59.36, Fe 57.8, Al2O3 4.937 (READ). Morley 2012 Table 3 gas partial pressures below the cloud, per unit P_T x metallicity: Cr 8.87e-7, Mn 6.32e-7, Na 3.97e-6, Zn 8.45e-8, KCl 2.55e-7, Fe 5.78e-5 (READ).

### 2.4 Ice sublimation temperatures versus partial pressure (K)

Fray and Schmitt 2009 polynomials ln P(bar) = A0 + sum A_i / T^i, as quoted in Lisse et al. arXiv:2009.02277 SOM Table 2 (2ND-HAND, all rows; the inversions are COMPUTED). Units were established as bar by checking each fit at its triple point (COMPUTED, all within 2 percent except N2-2: 0.1217 bar against 0.125). T at which the ice sublimation pressure equals the species partial pressure (columns, bar); ">Ttr" means the partial pressure exceeds the triple-point pressure, so a liquid forms instead (use 2.5). Values for the Fray and Schmitt fits below their stated lower limit are extrapolations (the lower limits of the fits were not in the quoted table: NOT REACHED).

| Ice (fit) | 1e-9 | 1e-7 | 1e-5 | 1e-3 | 1e-2 | 0.1 | 1 | T_tr (K) | P at T_tr (bar) |
|---|---|---|---|---|---|---|---|---|---|
| Ne | 8.4 | 9.9 | 12.2 | 15.6 | 18.2 | 21.7 | >Ttr | 24.5 | 0.3955 |
| N2 (beta, N2-2) | 26.3 | 30.6 | 36.7 | 46.0 | 52.8 | 62.2 | >Ttr | 63.15 | 0.1217 |
| CO (CO-2) | 29.1 | 33.6 | 40.0 | 49.5 | 56.5 | 65.9 | >Ttr | 68.1 | 0.1538 |
| Ar | 30.4 | 35.5 | 42.6 | 53.4 | 61.2 | 71.7 | >Ttr | 83.8 | 0.6882 |
| O2 (O2-2) | 31.6 | 36.4 | 42.9 | 52.2 | >Ttr | >Ttr | >Ttr | 54.4 | 0.0023 |
| CH4 | 37.4 | 43.7 | 52.7 | 66.3 | 76.1 | 89.5 | >Ttr | 90.69 | 0.1188 |
| Kr | 42.3 | 49.3 | 59.1 | 73.8 | 84.4 | 98.7 | >Ttr | 115.8 | 0.7352 |
| Xe | 58.9 | 68.6 | 82.1 | 102.4 | 117.0 | 136.6 | >Ttr | 161.4 | 0.8188 |
| C2H6 | 67.8 | 77.3 | 89.9 | >Ttr | >Ttr | >Ttr | >Ttr | 90.35 | 1.14e-5 |
| C2H4 | 64.6 | 73.7 | 85.9 | 103.1 | >Ttr | >Ttr | >Ttr | 103.9 | 0.0012 |
| C2H2 | 74.3 | 85.9 | 101.8 | 124.9 | 140.8 | 161.5 | 189.3 | 192.4 | 1.245 |
| H2S | 80.3 | 93.0 | 110.5 | 136.1 | 153.9 | 177.1 | >Ttr | 187.7 | 0.2365 |
| CO2 | 86.1 | 98.1 | 114.0 | 136.3 | 151.2 | 170.1 | 194.6 | 194.7 | 1.009 |
| NH3 | 101.5 | 115.7 | 134.6 | 161.0 | 178.6 | >Ttr | >Ttr | 195.4 | 0.0611 |
| SO2 | 111.9 | 126.9 | 146.7 | 174.3 | 192.7 | >Ttr | >Ttr | 197.6 | 0.0171 |
| HCN | 125.2 | 142.5 | 165.8 | 198.9 | 221.4 | 250.1 | >Ttr | 259.0 | 0.1828 |
| CH3OH | 141.5 | 160.0 | >Ttr | >Ttr | >Ttr | >Ttr | >Ttr | 175.6 | 2.3e-6 |
| C6H6 | 148.7 | 169.2 | 196.2 | 233.5 | 257.9 | >Ttr | >Ttr | 278.7 | 0.0513 |
| H2O ice (Feistel and Wagner 2007) | 161.1 | 183.3 | 212.6 | 252.8 | >Ttr | >Ttr | >Ttr | 273.16 | 0.006117 |

Fray and Schmitt coefficients (A0, A1, ...), as used (2ND-HAND from arXiv:2009.02277; two typos in the quoted table were corrected: SO2 "3.5.08" read as 3.508e3 and Kr "-s1.210e8" read as -1.210e8; both confirmed only by the triple-point check):
CH4: 10.51, -1110, -4341, 1.035e5, -7.910e5. CO-1: 10.43, -721.3, -1.074e4, 2.341e5, -2.392e6, 9.478e6; CO-2: 10.25, -748.2, -5843, 3.939e4. CO2-1: 14.76, -2571, -7.781e4, 4.325e6, -1.207e8, 1.350e9; CO2-2: 18.61, -4154, 1.041e5 (use it only from 195.3 K up; below that it is 4 K too cold at 100 Pa and 17 K too cold at 1 mPa, so use CO2-1: part 01 section 7.4, DOC.74). C2H6: 15.11, -2207, -2.411e4, 7.744e5, -1.161e7, 6.763e7. C2H4: 15.40, -2206, -1.216e4, 2.843e5, -2.203e6. C2H2: 13.40, -2536. HCN: 13.93, -3624, -1.325e5, 6.314e6, -1.128e8. N2-1 (alpha, below 35.6 K): 12.40, -807.4, -3926, 6.297e4, -4.633e5, 1.325e6; N2-2 (beta): 8.514, -458.4, -1.987e4, 4.800e5, -4.524e6. NH3: 15.96, -3537, -3.310e4, 1.742e6, -2.995e7. H2S: 12.98, -2707. SO2: 15.60, -3508, -9.401e4, 4.152e6, -6.946e7. Ne: 9.886, -269.9, 128.3, -162.4. Ar: 10.69, -893.2, -3567, 6.574e4, -4.280e5. Kr: 10.77, -1223, -8903, 2.635e5, -4.260e6, 3.575e7, -1.210e8. Xe: 10.98, -1737, -1.332e4, 4.349e5, -7.027e6, 4.447e7. O2-1: 15.41, -1148, 334.9, 60.21; O2-2: 13.35, -1012, -2971, 2.926e4. CH3OH: 19.18, -5648. C6H6: 17.35, -5663.
Water ice (Feistel and Wagner 2007 as quoted): ln(P/P_tr) = 1.5 ln(T/T_tr) + (1 - T_tr/T) h(theta), theta = T/273.16, P_tr = 611.657 Pa, h = sum e_i theta^i with e0..e6 = 20.9969665107897, 3.72437478271362, -13.9205483215524, 29.6988765013566, -40.1972392635944, 29.7880481050215, -9.13050963547721. Computed checks: 611.7 Pa at 273.16 K, 76.03 Pa at 250 K, 0.1626 Pa at 200 K.

### 2.5 Liquid-branch and trace-species anchors (NIST WebBook and chemicals/thermo)

Liquid saturation temperature T_sat (K) at total or partial pressure P (bar), from correlations in the `chemicals` 1.5.2 / `thermo` 0.6.1 packages (all rows DATA; method named per row; "-" means below the triple-point pressure, above the critical point, or outside the correlation). Triple and critical points from the same packages.

| Species | Method | T_tr (K) | P_tr (bar) | T_c (K) | 0.01 | 0.1 | 1 | 10 | 100 |
|---|---|---|---|---|---|---|---|---|---|
| H2O | WAGNER_MCGARRY | 273.16 | 0.006117 | 647.1 | 280.1 | 319.0 | 372.8 | 453.0 | 584.1 |
| NH3 | WAGNER_POLING | 195.49 | 0.06053 | 405.6 | - | 202.0 | 239.6 | 298.1 | 398.4 |
| CH4 | WAGNER_MCGARRY | 90.69 | 0.117 | 190.6 | - | - | 111.5 | 149.1 | - |
| CO2 | WAGNER_MCGARRY | 216.59 | 5.18 | 304.1 | - | - | - | 233.0 | - |
| N2 | WAGNER_MCGARRY | 63.15 | 0.1252 | 126.2 | - | - | 77.2 | 103.7 | - |
| Ar | WAGNER_MCGARRY | 83.81 | 0.6889 | 150.7 | - | - | 87.2 | 116.6 | - |
| CO | WAGNER_MCGARRY | 68.16 | 0.1554 | 132.9 | - | - | 81.5 | 109.0 | - |
| H2S | DIPPR_PERRY_8E | 187.70 | 0.2326 | 373.1 | - | - | 212.7 | 272.1 | - |
| SO2 | DIPPR_PERRY_8E | 197.70 | 0.01666 | 430.6 | - | 221.2 | 262.9 | 328.3 | - |
| HCN | DIPPR_PERRY_8E | 259.86 | 0.1872 | 457.0 | - | - | 298.5 | 375.9 | - |
| C2H6 | WAGNER_MCGARRY | 90.37 | 1.14e-5 | 305.3 | 127.8 | 150.3 | 184.3 | 241.1 | - |
| O2 | WAGNER_MCGARRY | 54.36 | 0.001463 | 154.6 | 61.3 | 72.7 | 90.1 | 119.6 | - |
| HCl | DIPPR_PERRY_8E | 159.07 | 0.1383 | 324.7 | - | - | 187.8 | 241.1 | - |
| H2SO4 (pure acid only; not for Venus-type 75 wt percent droplets) | AMBROSE_WALTON | 277.31 | 4.1e-9 | 924.0 | 459.8 | 522.1 | 610.7 | 748.6 | - |

Cross-check: the NH3 value 239.6 K at 1 bar coincides with the switch point between the two NIST (Stull 1947) Antoine sets for NH3, 239.6 K, which is its normal boiling point (COMPUTED vs DATA).

NIST WebBook phase-change anchors and Antoine sets for trace and volcanic species (DATA, NIST Chemistry WebBook SRD 69, phase-change pages; log10 P(bar) = A - B/(T + C); Antoine sets describe the liquid, so using them below the triple or melting point or outside the fit range is an extrapolation):

| Species | Normal boiling / triple / critical | Antoine set(s): T range K; A; B; C; reference | Status |
|---|---|---|---|
| NaCl | melting point not in the parsed page | 860 to 1240; 3.56682; 5200.904; -317.409; 1138 to 1738; 5.07184; 8388.497; -82.638 (Stull 1947) | DATA |
| KCl | boils about 1680 K (from the fit: 1 bar at 1679 K) | 1170 to 1466; 4.61668; 6910.833; -176.083; 1094 to 1680; 4.78236; 7440.691; -122.709 (Stull 1947) | DATA |
| Na | T_boil 1156 K, T_fus 370.96 K | 924 to 1118; 2.46077; 1873.728; -416.372 | DATA |
| K | T_boil 1037 K, T_triple 336.35 K | 679.4 to 1033; 4.45718; 4691.58; 24.195 | DATA |
| Rn | NOT in the page text (no T_boil parsed) | 129 to 211.8; 3.1908; 558.293; -36.803 (Stull 1947); gives 1 bar at 212 K | DATA |
| O3 | T_boil 161.3 K, T_fus 80.7 K, T_c 261.15 K | 92.8 to 162; 4.23637; 712.487; 6.982 (Stull 1947) | DATA |
| HCl | T_fus 161.15 K; P_tr 0.138 bar; T_c 324.68 K | 122.3 to 188.3; 3.60765; 535.172; -39.847; 188.3 to 309.4; 4.57389; 868.358; 1.754 (Stull 1947) | DATA |
| HF | T_boil 292.7 K, T_fus 190 K | 198.5 to 292.9; 4.16129; 1142.985; -17.993 (Stull 1947) | DATA |
| CH3Cl | T_boil 247 K, T_triple 175.43 K, P_tr 0.0087 bar, T_c 416 K | 183 to 249.4; 4.15454; 916.223; -28.466 | DATA |
| H2CO | T_boil 253.85 K, T_triple 155.1 K | 163.76 to 250.86; 4.28176; 959.43; -29.758 | DATA |
| OCS | T_triple 134.3 K, T_c 378.8 K | 161.8 to 223.84; 4.0436; 808.49; -22.72 | DATA |
| NO2 | T_boil 295.08 K | 217.6 to 294; 3.35248; 540.635; -131.93 (Stull 1947) | DATA |
| NH4Cl | sublimation enthalpy 177 kJ/mol (no Antoine set) | none | DATA (anchor only) |
| H2, He | H2 T_tr 13.95 K, P_tr 0.0721 bar, T_c 33.18 K; He T_c 5.2 K | H2: 21.01 to 32.27; 3.54314; 99.395; 7.726 | DATA |

NIST pages for SiO, TiO, VO, H2SO4, Fe, CaCO3, Cr2O3, S8, N2O4 and HNO3 gave no usable phase table (NOT REACHED for these through the WebBook; JANAF condensed-phase tables were not fetched).

### 2.6 Condensate sequence by temperature window (what to put in a profile)

All rows COMPUTED from the READ fits in 2.1 and the 2ND-HAND ice curves in 2.4 (solar metallicity unless stated; pressure range 0.1 to 1 bar for the hot species, 0.1 to 10 bar for NH4SH and H2S). Overlaps in the hot rows are real: a profile crosses several curves, and the deepest condensate removes its limiting element first (3.3).

| Window of local T (K) | Condensate | T_cond range (K) | Limiting element, relative mass note |
|---|---|---|---|
| 1850 to 2000 | hibonite, grossite, corundum, Ca-Ti oxides, perovskite | 1850 to 2004 | Al and Ti; Al-group about 2000 relative to Na2S = 1000, Ti-group about 150 to 180 |
| 1690 to 1840 | Fe metal | 1689 to 1838 | Fe, relative 21032 (Na2S = 1000) |
| 1600 to 1700 | forsterite | 1597 to 1698 | Mg, Si; relative 32180 |
| 1510 to 1600 | enstatite | 1513 to 1597 | relative 44585 |
| 1420 to 1520 | Cr | 1416 to 1521 | Cr, relative 303 |
| 1270 to 1340 | MnS | 1271 to 1343 | Mn, relative 360 |
| 930 to 1000 | Na2S | 929 to 996 | Na, relative 1000 |
| 750 to 800 | KCl and ZnS | 749 to 801 and 760 to 798 | K (relative 0.12 in Morley units) and Zn (0.05) |
| 190 to 225 | NH4SH | 188.5 to 223.4 (0.1 to 10 bar, [M/H] = 0); 204.5 to 246.2 at [M/H] = +1 | S or N, whichever is lower |
| 105 to 130 | H2S ice (only if NH3 is absent or less abundant than H2S) | 105.2 to 128.3 at [S/H] = 0 | S |
| 160 to 255 | water ice (partial pressure 1e-9 to 1e-3 bar) | 161 to 253 | O |
| 115 to 195 | CO2 ice (partial pressure 1e-5 to 1 bar) | 114 to 194.6 | C |
| 135 to 180 | NH3 ice (partial pressure 1e-5 to 1e-2 bar) | 134.6 to 178.6 | N |
| 52 to 90 | CH4 ice (partial pressure 1e-5 to 0.1 bar) | 52.7 to 89.5 | C |
| 37 to 62 | N2 ice (partial pressure 1e-7 to 0.1 bar), CO ice 34 to 66 K | 30.6 to 62.2 | N |

NH4H2PO4 (Visscher 2006 eq. 52, READ) would condense at 311 K (0.1 bar), 333 K (1 bar), 360 K (10 bar) at [M/H] = 0 (COMPUTED) if phosphorus were in equilibrium; Visscher states this equilibrium is not attained in the giant planets of the Solar System, so treat it as an upper-bound flag, not a default cloud.
Between about 300 and 750 K the reached curves cover only liquid water (2.5) and, at the top of the range, the salts and sulfides; the vapour curves for H2SO4-H2O, NH4Cl, ammonium salts, NaCl, perchlorates and carbonates are NOT REACHED (section 4). For T_eq below about 950 K photochemical organic haze dominates the aerosol opacity over silicate clouds in the Gao 2020 grid (READ), and KCl is hidden below the haze unless the haze is absent.
Sudarsky class sequence, Cahoy colour table and Venus and Earth cloud optical depths (20 to 30 and about 10) are quoted in the project note docs/design/composition/atmosphere-colour.md section 4 and its source check; they were not re-read here (2ND-HAND at best, not re-verified).

### 2.7 Cloud structure parameters (sedimentation efficiency and mixing)

| Item | Value | Source | Status |
|---|---|---|---|
| Ackerman and Marley baseline rain factor f_rain (f_sed in the notation of Gao 2018) | 3 | A and M 2001, astro-ph/0103423 table 1 | READ |
| Baseline K_min 1e5 cm2/s, mixing-length floor Lambda 0.1, sigma_g 2; mixing length L = H max(Lambda, Gamma / Gamma_adiabatic) | as stated | A and M 2001 table 1, eq. 6 | READ |
| Free-convection eddy diffusivity | K = (H/3) (L/H)^(4/3) (R F / (mu rho_a c_p))^(1/3), coefficient 1/3 | A and M 2001 eq. 5 | READ |
| Calibration point only: giant-planet cloud base | H = L = 20 km, K = 2e8 cm2/s, w* = K/L = 1.1 m/s | A and M 2001 | READ |
| f_sed range used in the literature | 0.01 (super-Earths, Morley 2015) to nearly 10 (brown dwarfs, Saumon and Marley 2008) | Gao, Marley and Ackerman 2018, arXiv:1802.06241 | 2ND-HAND |
| f_sed from microphysics, KCl, homogeneous nucleation | 0.125, 0.093, 0.025 for K_zz = 1e6, 1e7, 1e8 cm2/s at log g = 5.25; 0.050 and 0.036 at log g 4.25 and 3.25 (K_zz 1e8) | Gao et al. 2018 table 1 | READ |
| f_sed vs surface energy of KCl (K_zz 1e7) | 0.036 (0.25x), 0.059 (0.5x), 0.093 (1x), 0.150 (2x), 0.353 (3x), above 10 (4x) | Gao et al. 2018 table 1 | READ |
| f_sed with heterogeneous nuclei (K_zz 1e8) | above 10 (1 nm nuclei at 10 cm-2 s-1), 7.7 (100), 0.94 (1000), 0.94 (10 nm nuclei, 1000), 1.31 (0.1 nm, 1000) | Gao et al. 2018 table 1 | READ |
| Model cloud albedo | cloudy 400 to 800 K sub-Neptunes: 0.05 to 0.20; hazy: 0.00 to 0.03 | Morley 2015, arXiv:1511.01492 abstract | READ |
| Observed cloudiness trends with T_eq (transmission) | cloudiness rises for T_eq above 2300 K and for 1100 to 1600 K, falls for 1600 to 2300 K and below 1100 K | Gao 2020 | READ |
| Silicate versus haze crossover | silicates dominate aerosol opacity for T_eq above about 950 K | Gao 2020 | READ |

Nightside brightness temperature of hot Jupiters stays near 1100 K while dayside rises with T_eq, tied by one interpretation to the condensation temperature of a cloud (Gao, Wakeford, Moran and Parmentier 2021, arXiv:2102.03480): 2ND-HAND, interpretive, not a generator input.

### 2.8 Optical constants, density and Mie scattering of condensates, hazes and dust

Method (COMPUTED): m = n - ik interpolated at 450, 550 and 650 nm (k floor 1e-9) from the named file; Lorenz-Mie for homogeneous spheres (miepython 3.3.0, checked against an independent Bohren and Huffman code to 6 digits); each result averaged over radius +-10 percent (41 points) to damp Mie ripple. Angstrom exponent = -d ln Q_ext / d ln lambda between 450 and 650 nm (positive: bluer extinction; near 0: grey; negative: redder). SSA at 550 nm. One radius per case, no size distributions, mixing rules or aggregates; radii 0.1, 1, 5 micron shown (0.05, 0.5, 20 micron in the scratch results). Optical data are DATA files: Kitzmann and Heng 2018 LX-MIE compilation (raw GitHub exoclime/LX-MIE; arXiv:1710.04946 table 1 READ for provenance), HITRAN-RI ascii files, refractiveindex.info YAML (catalog paths main/H2O/nk/Hale.yml, main/H2O/nk/Warren-2008.yml, organic/CH4 - methane/nk/Martonchik-solid-90K.yml, main/Fe2O3/nk/Querry-o.yml, main/Fe3O4/nk/Querry.yml, main/Al2O3/nk/Malitson-o.yml, main/CaCO3/nk/Ghosh-o.yml, main/NaCl/nk/Querry.yml, main/KCl/nk/Querry.yml, main/ZnS/nk/Querry.yml, main/Fe/nk/Querry.yml, main/SiO2/nk/Malitson.yml), virga refrind files (Zenodo 15886530).

| Material | rho | n, k at 550 nm | SSA at r 0.1/1/5 um | Angstrom at r 0.1/1/5 um | Optical-constant source | Status |
|---|---|---|---|---|---|---|
| water, liquid | 0.997 R | 1.333, 2.0e-9 | 1/1/1 | 3.22/0.15/-0.04 | Hale and Querry 1973, RII | DATA |
| water ice Ih | 0.93 V | 1.311, 2.3e-9 | 1/1/1 | 3.25/-0.55/-0.05 | Warren and Brandt 2008, HITRAN-RI | DATA |
| NH3 ice | 0.84 V | 1.436, 2.9e-5 | 1/0.999/0.996 | 3.41/0.53/-0.03 | virga file | DATA |
| CH4 ice (90 K) | 0.49 V | 1.323, 1.0e-9 | 1/1/1 | 3.19/-0.19/-0.00 | RII | DATA |
| CO2 ice | 1.56 R | 1.413, 8.2e-7 | 1/1/1 | 3.27/0.90/-0.07 | no source | READ |
| H2SO4 75 wt% | 1.67 T | 1.431, 1.0e-9 | 1/1/1 | 3.17/0.65/-0.02 | Palmer and Williams 1975, HITRAN-RI | DATA |
| KCl | 1.99 G | 1.493, 1.0e-9 | 1/1/1 | 3.32/-0.54/-0.05 | Palik/Querry, LX-MIE | DATA |
| ZnS | 4.04 G | 2.389, 3.9e-6 | 1/1/0.999 | 3.49/-0.06/-0.05 | Palik/Addamiano, LX-MIE | DATA |
| Na2S | 1.86 G | 1.827, 0.0146 | 0.938/0.748/0.576 | 4.29/-0.05/-0.04 | LX-MIE | DATA |
| MnS | 4 G | 2.869, 1.0e-3 | 0.998/0.968/0.899 | 0.02/-0.12/-0.04 | LX-MIE | DATA |
| FeS (troilite) | 4.6 W | 1.415, 1.42 | 0.399/0.610/0.653 | -0.10/-0.14/-0.06 | Pollack 1994, LX-MIE | DATA |
| Cr | 7.15 G | 3.116, 4.42 | 0.711/0.780/0.799 | -0.02/0.04/0.02 | Lynch and Hunter, LX-MIE | DATA |
| Fe | 7.87 G | 2.081, 2.69 | 0.595/0.700/0.733 | -0.13/-0.06/-0.06 | Lynch and Hunter, LX-MIE | DATA |
| forsterite (amorph.) | 3.21 G | 1.618, 1.5e-4 | 0.999/0.995/0.983 | 3.29/-0.11/-0.03 | sol-gel, LX-MIE | DATA |
| enstatite (amorph.) | 3.2 W | 1.553, 1.4e-4 | 0.999/0.996/0.984 | 3.26/-0.59/-0.05 | sol-gel, LX-MIE | DATA |
| enstatite glass | 3.2 W | 1.576, 2.2e-5 | 1/0.999/0.997 | 3.34/-0.43/-0.02 | LX-MIE | DATA |
| fayalite | 4.39 W | 1.851, 9.7e-4 | 0.996/0.971/0.904 | 3.77/-0.22/-0.04 | LX-MIE | DATA |
| Al2O3 amorphous | 3.99 G | 1.580, 0.0373 | 0.773/0.582/0.541 | 2.80/-0.29/-0.04 | LX-MIE | DATA |
| Al2O3 sapphire | 3.99 G | 1.770, 1.0e-9 | 1/1/1 | 3.66/0.11/-0.05 | RII (k set 1e-9) | DATA |
| TiO2 anatase | 4.25 G | 2.406, 1.3e-3 | 0.996/0.963/0.881 | 3.30/0.01/-0.04 | LX-MIE | DATA |
| CaTiO3 | 4 W | 2.332, 1.8e-3 | 0.994/0.950/0.843 | 3.74/-0.03/-0.04 | LX-MIE | DATA |
| spinel | 3.6 W | 1.719, 1.0e-9 | 1/1/1 | 3.55/0.11/-0.03 | LX-MIE | DATA |
| periclase | 3.8 W | 1.740, 6.6e-8 | 1/1/1 | 3.62/0.16/-0.04 | Roessler and Huffman, LX-MIE | DATA |
| SiO2 | 2.65 W | 1.546, 6.4e-6 | 1/1/0.999 | 3.29/-0.60/-0.07 | LX-MIE | DATA |
| SiC | 3.2 R | 2.524, 9.1e-4 | 0.997/0.971/0.908 | 2.66/-0.13/-0.04 | Laor and Draine 1993, LX-MIE | DATA |
| NaCl | 2.17 W | 1.547, 1.0e-9 | 1/1/1 | 3.36/-0.57/-0.07 | Eldridge and Palik, LX-MIE | DATA |
| hematite (Triaud) | 5.3 W | 3.025, 0.441 | 0.603/0.594/0.630 | -0.50/-0.06/-0.02 | Triaud (Jena), LX-MIE | DATA |
| hematite (Querry) | 5.3 W | 3.318, 0.498 | 0.563/0.610/0.644 | -1.59/-0.06/-0.02 | RII | DATA |
| magnetite | 5.15 W | 2.343, 0.105 | 0.760/0.555/0.591 | 2.68/-0.08/-0.04 | RII | DATA |
| calcite | 2.71 W | 1.661, 1.0e-9 | 1/1/1 | 3.47/-0.05/-0.05 | RII (n only, k set 1e-9) | DATA |
| Titan tholin (Khare) | 1 A | 1.702, 0.0256 | 0.873/0.668/0.552 | 3.80/0.07/-0.04 | Khare 1984, HITRAN-RI | DATA |
| Titan analogue (Ramirez) | 1 A | 1.568, 0.002 | 0.984/0.944/0.830 | 3.50/-0.44/-0.04 | HITRAN-RI | DATA |
| Pluto tholin | 1 A | 1.694, 0.00907 | 0.950/0.827/0.611 | 4.05/-0.03/-0.04 | HITRAN-RI | DATA |
| flame soot | 1.8 R | 1.773, 0.626 | 0.353/0.535/0.578 | 0.89/-0.09/-0.04 | Chang 1990, HITRAN-RI | DATA |
| amorphous carbon | 1.8 R | 2.334, 0.979 | 0.451/0.580/0.619 | 0.27/-0.09/-0.04 | Stagg 1993, HITRAN-RI | DATA |
| graphite | 2.2 W | 2.567, 1.19 | 0.476/0.599/0.637 | 0.03/-0.09/-0.04 | LX-MIE | DATA |
| Saharan dust, Burkina Faso | 2.6 R | 1.530, 0.00673 | 0.942/0.864/0.646 | 3.31/-0.65/-0.06 | Wagner 2012, HITRAN-RI | DATA |
| Saharan dust, Cairo | 2.6 R | 1.530, 0.00351 | 0.969/0.920/0.749 | 3.20/-0.69/-0.05 | Wagner 2012, HITRAN-RI | DATA |
| sea salt (Shettle) | 2.2 R | 1.500, 1.0e-8 | 1/1/1 | 3.27/-0.57/-0.04 | Shettle, AFCRL 1987, HITRAN-RI | DATA |
| ammonium sulfate (Shettle) | 1.77 R | 1.530, 1.0e-7 | 1/1/1 | 3.38/-0.63/-0.04 | Shettle, AFCRL 1987, HITRAN-RI | DATA |
| carbonaceous (Shettle) | 1.8 R | 1.750, 0.44 | 0.370/0.522/0.566 | 1.40/-0.09/-0.04 | Shettle, AFCRL 1987, HITRAN-RI | DATA |
| volcanic dust (Shettle) | 2.6 R | 1.500, 0.008 | 0.924/0.854/0.620 | 3.01/-0.56/-0.03 | Shettle, AFCRL 1987, HITRAN-RI | DATA |
| meteoric dust (Shettle) | 3 R | 1.513, 2.6e-4 | 0.997/0.993/0.972 | 3.17/-0.73/-0.03 | Shettle, AFCRL 1987, HITRAN-RI | DATA |
| hematite (Shettle) | 5.3 R | 3.102, 0.0925 | 0.857/0.600/0.631 | -1.06/-0.10/-0.03 | Shettle, AFCRL 1987, HITRAN-RI | DATA |
| sand, 10% hematite (Shettle) | 2.8 R | 1.664, 0.00474 | 0.971/0.888/0.702 | 3.53/-0.01/-0.05 | Shettle, AFCRL 1987, HITRAN-RI | DATA |
| dust-like (Shettle) | 2.6 R | 1.530, 0.008 | 0.931/0.844/0.621 | 3.05/-0.64/-0.05 | Shettle, AFCRL 1987, HITRAN-RI | DATA |
| water-soluble (Shettle) | 1.7 R | 1.530, 0.006 | 0.948/0.876/0.663 | 3.05/-0.67/-0.05 | Shettle, AFCRL 1987, HITRAN-RI | DATA |

Key: rho in g/cm3 with code G = Gao 2020 supplementary table 3, W = webmineral, V = virga gas_properties.py, T = Timmermans via HITRAN-RI (20 C), A = assumed value used by the haze papers, R = RECALLED (section 5). Status refers to the optical constants (DATA, or READ for CO2 ice); every Mie number is COMPUTED. RII = refractiveindex.info YAML, LX-MIE = Kitzmann and Heng compilation.

Colour effect from the k values at 450, 550 and 650 nm (COMPUTED from the same files): hematite absorbs blue strongly, k falling from 450 to 650 nm from 0.99 to 0.24 (Triaud), 1.0 to 0.057 (Querry), 0.31 to 0.0052 (Shettle), so hematite-bearing dust is reddish whatever the dataset; tholins also absorb blue, k from 0.057 to 0.012 (Khare), 0.034 to 0.0024 (Pluto), 0.006 to 0.0009 (Ramirez), giving yellow to brown haze; Saharan dust k falls from 0.016 to 0.0032 (Burkina Faso) and 0.0065 to 0.0024 (Cairo); MnS k falls from 0.039 to 0.00028; Na2S is nearly flat (0.02 to 0.013); magnetite k rises toward the red (0.032 to 0.14); Fe, Cr, graphite, amorphous carbon and soot are dark and grey (SSA 0.35 to 0.8, Angstrom near 0 at 1 micron). Non-absorbing white or grey condensates (SSA above 0.999 at 1 micron): water, ice, CH4, CO2, H2SO4, KCl, NaCl, ZnS, sapphire, spinel, periclase, quartz, calcite, sea salt, ammonium sulfate.

CO2 ice, (READ from an OCR of the PDF at atmos.washington.edu; the visible k values come from Egan and Spagnolo, which Warren flags as uncertain, possibly impurity-limited): lambda (um): n, k = 0.363: 1.432, 6.78e-7; 0.40: 1.426, 6.83e-7; 0.452: 1.420, 7.04e-7; 0.490: 1.417, 7.37e-7; 0.553: 1.413, 8.27e-7; 0.616: 1.411, 9.92e-7; 0.679: 1.409, 1.17e-6; 0.742: 1.407, 1.36e-6. Density 1.56 g/cm3 is RECALLED (section 5).

Where sources disagree (generator should draw within these spreads):
- Hematite k(550 nm): Triaud via LX-MIE 0.44, Querry ordinary ray 0.50, Shettle 0.0925. Hematite Angstrom exponent for 0.1 micron grains: -0.50, -1.59, -1.06 respectively (reddish extinction in all three; magnitude is uncertain).
- Alumina k(550): amorphous/porous Begemann 0.037 against sapphire about 0 (Malitson, k floor 1e-9), a spread of many orders of magnitude; choose by sample type (porous or amorphous condensate versus crystalline).
- TiO2 k(550): LX-MIE anatase 1.3e-3 against virga refrind file about 0.41. CaTiO3: LX-MIE k 1.8e-3; the virga CaTiO3 refrind file has flat placeholder rows.
- virga refrind files have placeholder rows (CaAl12O19 zeros, MnS k = 1e-30 at 650 nm); prefer LX-MIE where it exists.
- Troilite (FeS) from LX-MIE reads n = 1.415, k = 1.4 at 550 nm; the n looks low for a metallic sulfide, so treat the FeS row as unverified (section 5).
- CH4 ice has no measured visible k in the Martonchik 90 K file; the table shows the floor value 1e-9, not a measurement. NH3 ice and CO2 ice k are laboratory upper limits to values of order 1e-5 to 1e-6.
- Wagner 2012 Saharan dust has n fixed at 1.53; Shettle volcanic, dust-like and water-soluble classes have k 0.006 to 0.008 at 550 nm.

### 2.9 Photochemical and thermochemical haze: production, size, density, colour

Laboratory plasma production rates in mg/h (Hörst et al. 2018, arXiv:1801.06512, figure 1 as printed; 72 h runs; "film" is a visible deposit too thin to weigh; standard Titan gas, 5 percent CH4 in N2, gives 7.4 mg/h in the same setup; the authors state that lab rates are comparable between cases but cannot be converted to planetary rates). READ.

| Metallicity of the starting gas | 300 K | 400 K | 600 K |
|---|---|---|---|
| 100x (about 70 % H2, 6 to 8 % H2O, 3 to 4.5 % CH4) | film | 0.25 | 0.04 |
| 1000x (300 K: 66 % H2O, 6.6 % CH4, 6.5 % N2, 4.9 % CO2; 400 K: 55.5 % H2O, 10.5 % CH4, 9.8 % CO2, 6.4 % N2; 600 K: 42 % H2, 20 % CO2, 16 % H2O, 5 % N2, 1.9 % CO, 1.7 % CH4) | 10.43 | 10.00 | 0.15 |
| 10000x (about 66 % CO2, 12 to 13 % N2, 6 to 16 % H2O; no CH4) | film | film | 0.31 |

The paper's text says the H2-dominated runs give 0.04 and 0.25 mg/h while the figure lists 0.15 for 1000x at 600 K; both are as printed. Quoted conclusions: the two highest rates are the H2O-dominated 1000x runs at 300 and 400 K; H2 lowers the rate; CO can raise it; the rate is not a function of carbon abundance alone; the 10000x run at 600 K makes haze with no CH4. He et al. 2018 (arXiv:1803.01706, READ): film colours dark gold (1000x at 300 and 400 K, like Titan standard), olive green (100x, 600 K), chocolate brown (100x, 400 K), light yellow (100x at 300 K, 1000x and 10000x at 600 K), none (10000x at 300 and 400 K); particle diameters 20 to 180 nm (100x: 20 to 55 nm; others 25 to 180 nm; Titan standard 25 to 80 nm, mean 46.7 nm), Rayleigh regime in the visible. He et al. 2018b (arXiv:1805.10488, READ): UV-driven rates are far lower, at most 0.060 mg/h (1000x, 300 K) and falling with temperature except for 10000x; density 1.38 g/cm3 assumed. He et al. 2020 (arXiv:2008.09700, READ): 800 K, 100x and 1000x: diameters 20 to 140 nm, 2.4e-5 to 9.7e-5 mg cm-3 h-1 (plasma about three times UV).

Model-side production and size rules:

| Quantity | Value or rule | Source | Status |
|---|---|---|---|
| Titan reference monomer mass production | 1e-14 g cm-2 s-1 (0.5e-14 to 2e-14) | Kawashima and Ikoma 2018, arXiv:1712.02808 | 2ND-HAND |
| Scaling with stellar Lyman-alpha | Mdot = beta (I_Lya / I_Lya,Titan) Mdot_Titan; I_Lya,Titan = 6.2e9 photons cm-2 s-1 (Trainer 2006); beta a free efficiency | K and I 2018 eqs. 33, 35 | 2ND-HAND |
| Vertical distribution | proportional to (f_HCN + f_C2H2) N(z); paper II uses photodissociation rates of CH4, HCN, C2H2 | K and I 2018 eq. 34; arXiv:1904.09986 | READ |
| Monomer radius; density | 1 nm; 1.0 g/cm3 (K and I II, Adams 2019, Gao 2020); 1.38 (He 2018b); 1.4 quoted for He 2017 | as named | READ, 2ND-HAND |
| Diffusion-limited CH4 supply | P(z) = N K_zz d f_CH4 / dz with d f / dz about f_CH4 / H_a, formed at 1 microbar; always below the Lyman-alpha-limited rate | Adams et al. 2019, arXiv:1902.05231 eq. 7; Gao 2020 | READ |
| Adams 2019 grid | 1e8 to 1e11 CH4 molecules cm-2 s-1 (2.55e-15 to 2.55e-12 g cm-2 s-1); K_zz 1e8 to 1e10 cm2/s; monomers 1 and 10 nm; CH4 at 1 microbar 30 ppm in equilibrium, diffusion-limited 1.3e11 molecules cm-2 s-1 (3.5e-12 g cm-2 s-1, COMPUTED conversion) | Adams et al. 2019 | READ |
| Aggregates | fractal dimension 1.5 (dimers) rising to 2.4 above about 2000 monomers; opacity flatter (grey) than spheres | Adams et al. 2019 | READ |
| Gao 2020 haze | soot-like, 1 g/cm3, minimum radius 50 nm | Gao 2020 | READ |
| Calibration only: Pluto | column mass flux about 1.2e-14 g cm-2 s-1 (equals column CH4 destruction); bimodal haze in later work | Gao 2017, arXiv:1610.01679; Fan 2022, arXiv:2201.04392 | READ; abstract only |

Haze generator rule (engineering, COMPUTED/RECALLED mixture, not validated): use the diffusion-limited CH4 supply as an upper bound and multiply by a yield factor drawn between about 0.01 and 1 (lab rates span a factor of 250); near zero for carbon-free gas; lower for H2-rich gas and above about 800 to 1000 K unless CO-rich; monomers 1 to 50 nm that coagulate to the sizes in the microphysics papers (coagulation not coded here).

### 2.10 Dust: size classes, lifting threshold, settling, Mars-type optical constants

Dust size classes and emission (Kok 2012, arXiv:1201.4353, READ): long-term suspension below about 20 micron diameter, short-term suspension 20 to 70 micron, saltation 70 to 500 micron. Emitted dust has three size modes with median diameters 1.5, 6.7 and 14.2 micron; the threshold is lowest for grains of roughly 75 to 100 micron on a thick-air planet.

Threshold friction velocity (formulas in 3.8). COMPUTED values for the four reference conditions of Kok 2012 table 2.3, used only as calibration points (Kok lists the fluid-threshold u* as about 0.2 m/s for Earth, 1.5 for Mars, 0.02 for Venus, 0.04 for Titan; Tokano 2010 infers 0.036 from Titan dune orientation):

| Calibration condition (Kok table 2.3) | Grain density (kg/m3) | Air density (kg/m3) | Viscosity (Pa s) | g (m/s2) | Shao and Lu minimum u* (m/s) at D | Shao and Lu with gamma 1.65e-4 to 5e-4 N/m | Iversen and White minimum u* | Status |
|---|---|---|---|---|---|---|---|---|
| quartz on dense thin air (Earth-like) | 2650 | 1.2 | 1.8e-5 | 9.81 | 0.239 at 107 um | 0.206 to 0.272 | 0.207 at 75 um | READ formula, COMPUTED |
| basalt on thin CO2 (Mars-like) | 3000 | 0.02 | 1.2e-5 | 3.71 | 1.50 at 165 um | 1.29 to 1.70 | 1.44 at 123 um | READ formula, COMPUTED |
| basalt on very dense CO2 (Venus-like) | 3000 | 66 | 3.2e-5 | 8.87 | 0.0323 at 107 um | 0.0278 to 0.0367 | 0.0247 at 73 um | READ formula, COMPUTED |
| tholin or ice grains on cold N2 (Titan-like) | 1000 | 5.1 | 6.3e-6 | 1.35 | 0.0554 at 472 um | 0.0477 to 0.0629 | 0.0348 at 240 um | READ formula, COMPUTED |

The computed minima agree with the Kok values to within a factor of about 1.6 (the Venus-like and Titan-like cases are the loosest); the two formulas differ from each other by 15 percent (Earth-like) to 60 percent (Titan-like), so the generator should draw the threshold inside the span of the two formulas and the surface-energy range, and the cohesion term matters most for small grains and weak gravity.

Settling and residence time of spherical particles (COMPUTED). Stokes fall speed v = (2/9)(rho_p - rho_a) g r^2 C / eta with the Cunningham factor of 3.6 (Kawashima and Ikoma 2018 eq. 21, READ; others use C = 1 + 1.26 Kn, which differs by 0.03 percent at Kn = 0.1, 5 percent at Kn = 1 and 21 percent at Kn = 10, COMPUTED); gas viscosity from the `chemicals`/`thermo` Lucas correlation (within 3 percent of DIPPR for N2, CO2 and H2 at the test temperatures). Particle density 2000 kg/m3 (scales linearly). Cells show v (m/s) / tau, the time to fall one scale height H = R T / (M g) without mixing. Particle Reynolds number stayed below 0.2 in every cell.

| Conditions (illustrative, not presets) | eta (uPa s) | lambda (um) | H (km) | r = 0.1 um | r = 1 um | r = 10 um |
|---|---|---|---|---|---|---|
| N2, 1 bar, 288 K, g 9.8 | 17.4 | 0.064 | 8.7 | 4.6e-6 / 2.2e4 d | 2.7e-4 / 373 d | 0.025 / 4.0 d |
| N2, 0.2 bar, 217 K, g 9.8 | 13.8 | 0.22 | 6.6 | 1.4e-5 / 5.6e3 d | 4.0e-4 / 189 d | 0.032 / 2.4 d |
| CO2, 600 Pa, 210 K, g 3.7 | 10.7 | 4.5 | 10.7 | 1.1e-4 / 1.1e3 d | 1.2e-3 / 101 d | 0.024 / 5.1 d |
| N2, 1.45 bar, 94 K, g 1.35 | 6.4 | 0.0092 | 20.7 | 1.1e-6 / 2.3e5 d | 9.5e-5 / 2.5e3 d | 9.4e-3 / 25 d |
| CO2, 1 bar, 340 K, g 8.9 | 17.1 | 0.055 | 7.2 | 4.0e-6 / 2.1e4 d | 2.5e-4 / 339 d | 0.023 / 3.6 d |
| H2, 1 bar, 165 K, g 25 | 6.1 | 0.063 | 27.2 | 3.4e-5 / 9.4e3 d | 2.0e-3 / 160 d | 0.18 / 1.7 d |
| H2, 1 mbar, 1500 K, g 10 | 27.9 | 871 | 619 | 0.023 / 311 d | 0.23 / 31 d | 2.3 / 3.1 d |
| N2, 1 Pa, 70 K, g 0.6 | 4.7 | 857 | 34.6 | 8.0e-3 / 50 d | 0.080 / 5.0 d | 0.80 / 0.5 d |

Mars-type airborne dust (optical properties):

| Quantity | Value | Source | Status |
|---|---|---|---|
| Effective radius, airborne dust | 1.85 +- 0.3 micron (low dust, effective variance 0.5) and 1.52 +- 0.3 micron (dust storm) | Pollack et al. 1995, NTRS record 19950051173 (abstract) | READ (abstract only) |
| Effective radius | 1 to 2 micron, Angstrom exponent near 0 | Lemmon et al. 2015, arXiv:1403.4234 | READ |
| Single scattering albedo at 650 nm | 0.975 (Wolff 2009), asymmetry parameter 0.673 +- 0.081, sizes 1.2 to 1.8 micron | quoted in Mars sky-brightness paper arXiv:1905.01074 | 2ND-HAND |
| Solar-band single scattering albedo | 0.86 to 0.88 (Pollack-type) against 0.93 to 0.94 (TES-type); the CRISM retrieval depends on the assumed size | Wolff abstract 9125.pdf (NTRS) | READ (abstract only) |
| Single scattering albedo, 443 and 671 nm | 0.84 and 0.92 | Tomasko 1999 | SNIPPET, do not rely |
| Mars analogue refractive indices | palagonite 1.52 + 0.0005i, basalt 1.52 + 0.001i, JSC-1A 1.65 + 0.003i; laboratory sample single scattering albedos 0.960, 0.892, 0.708 are for samples with effective radius 4.5, 6.9, 15.85 micron (not airborne sizes) | arXiv:1905.01074 | READ |
| Retrieved imaginary index, Mars simulants (20 to 40 micron fractions) | about 7e-4 (450 nm) and 6e-4 (550 nm) for JSC Mars-1; about 3e-4 for MGS-1 | no source | READ |
| k needed for omega0 = 0.975 at 650 nm, n = 1.5 | 8.9e-4 (r_eff 1.5 micron, variance 0.4); 7.4e-4 (1.85 micron, 0.5); 8.9e-4 (1.52 micron, 0.5); 1.1e-3 (1.2 micron, 0.4); 5.5e-4 (2.5 micron, 0.4); asymmetry parameter 0.71 to 0.76 | COMPUTED with miepython 3.3.0 and a gamma size distribution (Hansen and Travis form) | COMPUTED |
| Check of the same code on analogue indices at 650 nm, r_eff 1.5 micron, variance 0.4 | omega0 0.986 (palagonite), 0.972 (basalt), 0.922 (JSC-1A) | COMPUTED | COMPUTED |

Dust emission flux scalings (Kok 2012, READ; coefficients not transcribed): Kawamura F_d = C_K (rho_a / g) u*^3 (1 - u*_it^2 / u*^2)(1 + u*_it / u*); Gillette and Passi F_d = C_GP u*^4 (1 - u*_it / u*); the vertical dust flux scales with saltation flux, F_d = epsilon Q, epsilon from about 1e-5 to 1e-2 m-1 (the dust emission efficiency depends on soil type).

### 2.11 Lava worlds and hot rocky planets: vapour species and cloud mineralogy

Miguel et al. 2011 (arXiv:1110.2426, READ): magma outgassing with the MAGMA code, assuming volatile-free melt (H, C, N, S lost), Bond albedo A = 0.01, surface temperature T_p = T_eff,* ((1-A) R_*^2 / (4 a^2))^(1/4) (their eq. 1). Magma compositions (wt percent; Table 1, row labels not in the extracted text, rows identified by composition: row 1 komatiite, row 2 bulk silicate Earth; BSE row equals Schaefer and Fegley 2009 Table 1):

| Oxide | Komatiite (Miguel 2011 row 1) | Bulk silicate Earth |
|---|---|---|
| SiO2 | 47.10 | 45.97 |
| MgO | 29.60 | 36.66 |
| Al2O3 | 4.04 | 4.77 |
| TiO2 | 0.24 | 0.18 |
| Fe2O3 | 12.80 | 0.0 |
| FeO | 0.0 | 8.24 |
| CaO | 5.44 | 3.78 |
| Na2O | 0.46 | 0.35 |
| K2O | 0.09 | 0.04 |

Atmosphere types by surface temperature (komatiite; dominant gas order from most to least abundant), READ:

| Type | T_p range (K) in the text | Dominant gases |
|---|---|---|
| I | below 2033 | Na, O2, O, Fe |
| II | 2033 to 2588 | Na, O2, O, SiO, Fe, Mg |
| III | 2588 to 2890 | Na, O2, O, SiO, Mg, Fe (O2, O and SiO column densities nearly equal) |
| IV | 2890 to 3168 | SiO more abundant than O2; silicate and Na dominated |
| V | above 3168 | SiO, then O, O2, Na |

The caption of their partial-pressure figure uses transition temperatures 2096, 2460, 2735, 2974 and 3170 K for the same transitions, which differ from the text by up to 128 K, and shift slightly for the BSE composition (values garbled in the extraction). Columns are nearly independent of planet mass between 1 and 10 Earth masses, so surface temperature is the main control.
Schaefer and Fegley 2009 (arXiv:0906.1204, READ): for BSE, continental crust, oceanic crust and Moon-like melts, atmospheres are Na, O2, O and SiO in decreasing order; total pressure falls from about 1e-1 bar at 3000 K to about 1e-7 bar at 1500 K, and at 2200 K ranges from 5.5e-4 to 1e-3 bar across the four compositions. In BSE melt Na is completely lost from the system after only 2 percent vaporisation, after which Al, Ca, TiO2 and AlO gases dominate; continental crust keeps Na longer (monatomic Na to about 11 percent vaporisation). In a volatile-bearing case SO2 and NaCl are the major S and Cl species. A simple condensation calculation for a BSE atmosphere at 2200 K (g = 36.2 m/s2, dry adiabat) gives cloud layers: Mg as liquid enstatite below 2125 K (about 2 km, over 90 percent condensed below 4 km); Al as corundum and spinel (90 percent condensed below about 6 km); Ca as wollastonite (7 km); Si as MgSiO3, CaSiO3 and SiO2 (8 km); Fe as liquid FeO (10 km); Ti as CaTiO3 and MgTiO3 (15 km); Na and O as Na2O (28 km); K remains and condenses as metal near 50 km. Most condensed Mg, Al, Ca, Si and Fe are expected to rain back into the magma; Na and K may form extended clouds. VO and TiO gas abundances for such atmospheres: NOT REACHED (Wakeford 2017 notes Ti condensation removes TiO from the gas in giant planets).
LavAtmos (arXiv:2210.10463) and Schaefer, Lodders and Fegley 2012 (Vaporization of the Earth) were downloaded or identified but their partial-pressure tables were not extracted: NOT REACHED.

### 2.12 Mineral inventory densities (webmineral.com data pages, DATA)

Use for dust and crust grain density draws (g/cm3; ranges as printed; wuestite not on the pages fetched).

| Group | Mineral (formula) | Density |
|---|---|---|
| Silicates | forsterite Mg2SiO4 3.21 to 3.33; fayalite Fe2SiO4 4.39; olivine 3.27 to 3.37; enstatite MgSiO3 3.1 to 3.3; quartz SiO2 2.6 to 2.65; cristobalite SiO2 2.27; kaolinite 2.6 | 2.27 to 4.39 |
| Oxides | periclase MgO 3.67 to 3.9; corundum Al2O3 4 to 4.1; rutile TiO2 4.25; perovskite CaTiO3 4.0; spinel MgAl2O4 3.57 to 3.72; hibonite CaAl12O19 3.84; grossite CaAl4O7 2.88; hematite Fe2O3 5.3; magnetite Fe3O4 5.1 to 5.2; goethite FeOOH 3.3 to 4.3; eskolaite Cr2O3 5.18; chromite 4.5 to 5.09 | 2.88 to 5.3 |
| Carbonates | calcite CaCO3 2.71; dolomite 2.8 to 2.9; magnesite MgCO3 3.0; siderite FeCO3 3.96 | 2.71 to 3.96 |
| Sulfates | gypsum 2.3; bassanite 2.7; anhydrite 2.96 to 2.98; kieserite 2.57; epsomite 1.67 to 1.68; thenardite 2.67 to 2.7; mirabilite 1.49; jarosite 2.9 to 3.3 | 1.49 to 3.3 |
| Sulfides | troilite FeS 4.58 to 4.65; pyrite FeS2 5.0; sphalerite ZnS 3.9 to 4.2; alabandite MnS 3.95 to 4.04 | 3.9 to 5.0 |
| Salts | halite NaCl 2.17; sylvite KCl 1.99; natron 1.42 to 1.47; nahcolite 2.21 | 1.42 to 2.21 |
| Carbon | graphite 2.09 to 2.23 | 2.1 to 2.2 |
| Ices and droplets | water ice 0.93 (virga, 213 K); CH4 ice 0.49 and NH3 ice 0.84 (virga gas_properties.py); H2SO4 75 wt percent 1.67 (Timmermans via HITRAN-RI, 20 C) | 0.49 to 1.67 |

All entries DATA; no stability or decomposition data were extracted.


## 3. Formulas the generator can code (with validity)

3.1 Condensation from a fit (2.1): T_cond = 10^4 / (a - b log10 P_T + d (log10 P_T)^2 - c [Fe/H]), P_T in bar. Valid for scaled-solar element ratios and where the phase is the first of its group; the fit pressure ranges are not printed in the extracted texts (tabulated for 1e-4 to 100 bar, beyond that is extrapolation). For KCl, ZnS and the other effective expressions of 2.3 solve log10 P_sat(T) = log10 (x_i P_T) with x_i from part 05 (Morley 2012 table 3 lists the K, Zn, Na, Mn, Cr, Fe gas pressures). Changing C/O, Mg/Si or Fe/Mg breaks the fits; a Gibbs minimisation (FastChem, GGchem, LavAtmos) is outside this note.
3.2 Ices and liquids: T_cond solves P_sub(T) = x_i P (2.4); if x_i P exceeds the triple-point pressure use the liquid curve (2.5). For gases absent from 2.4 and 2.5 use ln P = ln P_tr - (L / R)(1/T - 1/T_tr) with NIST triple-point values (standard thermodynamics, RECALLED form; poor far from the triple point).
3.3 Cloud base and depletion: for each condensate sorted by T_cond, the base is the first crossing of T(P) with T_cond(P) going up from depth. Condensate mass mixing ratio at base is the limiting-element mass mixing ratio times M_condensate / (n M_element), minus what hotter condensates already removed (Visscher 2010, Wakeford 2017). Above the base the condensing gas follows saturation, x_i(P) = P_sat(T(P)) / P, which is the correction returned to part 05. FeS is suppressed once Fe metal has settled (Visscher 2006).
3.4 Sedimentation (Ackerman and Marley 2001, READ): q_t(z) = q_below exp(-f_sed (q_c / q_t) z / L) from -K dq_t/dz - f_sed w* q_c = 0, with w* = K / L, K from 2.7 bounded below by K_min. Draw K_zz between about 1e5 and 1e10 cm2/s (Gao 2018: 1e6 to 1e8; Adams 2019: 1e8 to 1e10) and f_sed between 0.01 and 10: about 0.03 to 0.15 for homogeneous KCl nucleation at K_zz 1e6 to 1e8, above 1 for high surface-energy materials or few nuclei (2.7).
3.5 Particle size (A and M eqs. 9, 11, 13, 14, READ): log-normal dn/dr = N / (r sqrt(2 pi) ln sigma_g) exp(-ln^2(r / r_g) / (2 ln^2 sigma_g)); fall speed v_f = w* (r / r_w)^alpha with v_f(r_w) = w*; r_g = r_w f_sed^(1/alpha) exp(-((alpha + 6) / 2) ln^2 sigma_g); N = (3 epsilon rho_a q_c / (4 pi rho_p r_g^3)) exp(-(9/2) ln^2 sigma_g), epsilon the condensate-to-air molecular-weight ratio. Area-weighted radius r_eff = r_g exp((5/2) ln^2 sigma_g) (COMPUTED algebra; eq. 17 was garbled in the extraction). alpha = 2 in Stokes flow without slip and about 1 in slip flow (COMPUTED from 3.6); the drag-law fits for Reynolds numbers 1 to 1000 (A and M appendix B) were not transcribed. Layer optical depth for condensate column mass m_c (kg/m2): tau = 3 Q_ext m_c / (4 rho_p r_eff) (geometric optics; Q_ext about 2 for r above about 1 micron; small particles from 2.8).
3.6 Settling (2.10): v = (2/9)(rho_p - rho_a) g r^2 C / eta, C = 1 + 1.257 Kn + 0.400 Kn exp(-1.10 / Kn) (READ), Kn = lambda / r, lambda = eta / (0.499 rho_a sqrt(8 R T / (pi M))) (kinetic theory, RECALLED coefficient; gives 0.064 micron for N2 at 1 bar and 288 K). A particle stays aloft where the mixing time H^2 / K_zz is shorter than the fall time H / v.
3.7 Haze: Mdot = beta (I_Lya / I_Lya,Titan) Mdot_Titan, or the diffusion-limited supply N K_zz f_CH4 / H_a as an upper bound (2.9). Steady-state haze mass column = Mdot / v_f, hence tau about 3 Q_ext Mdot / (4 rho_p r_eff v_f) (flux balance, COMPUTED/engineering, not validated against a microphysical model). Optical constants: Khare, Ramirez or Pluto tholin for cold organic haze, soot for carbon-rich hot haze (2.8).
3.8 Dust lifting threshold (Kok 2012, READ): Shao and Lu u*_ft = A_N sqrt((rho_p - rho_a) g D / rho_a + gamma / (rho_a D)), A_N = 0.111, gamma = 1.65e-4 to 5.0e-4 N/m (2.9e-4 from Kok and Renno 2006). Minimum at D_opt = sqrt(gamma / (rho_p g)): u*_min = A_N sqrt(2 sqrt(gamma rho_p g) / rho_a) for rho_p much larger than rho_a (COMPUTED, reproduces the numerical minima in 2.10). Bagnold: u*_ft = A_ft sqrt((rho_p - rho_a) g D / rho_a), A_ft about 0.10 (fluid) and 0.082 (impact). (Kok eq. 2.6): with R*_t = rho_a u*_ft D / mu and K = (1 + 0.006 / (rho_p g D^2.5))^(1/2): A_ft = 0.2 K / (1 + 2.5 R*_t)^(1/2) for 0.03 <= R*_t <= 0.3; 0.129 K / (1.928 R*_t^0.092 - 1)^(1/2) for 0.3 < R*_t < 10; 0.12 K [1 - 0.0858 exp(-0.0671 (R*_t - 10))] for R*_t >= 10; solve implicitly (constant 0.006 taken in cgs units, units not stated in the extracted text). Friction velocity from wind: u* = kappa U / ln(z / z0), kappa = 0.4, z0 about D/30 to 1e-3 m (both RECALLED). Impact-to-fluid threshold ratio (Kok table 2.3): about 0.8 for thick-ish air, about 0.1 for thin air, above 1 for dense air.
3.9 Dust emission (2.10): F_d = epsilon Q or the Kawamura or Gillette-Passi power laws; emitted dust has three size modes with median diameters 1.5, 6.7 and 14.2 micron; below about 20 micron grains stay aloft long term. Airborne column optical depth is a drawn parameter between clear and storm values; effective radii 1.5 to 1.85 micron; single scattering albedo from the inversion in 2.10.
3.10 Cloud albedo versus optical depth (Monte Carlo, COMPUTED): plane-parallel slab, Henyey-Greenstein phase function, single scattering albedo omega, black lower boundary, 20000 to 60000 photons per entry (noise about 0.002 to 0.004); checked against first-order single scattering (isotropic R = tau / (2 mu0) at tau = 0.01: 0.0100 and 0.0051 against 0.0100 and 0.0050) and R + T = 1 for omega = 1. Plane albedo R at solar-zenith cosine mu0, g = 0.85:

| g | omega | mu0 | tau 0.5 | 1 | 2 | 5 | 10 | 20 | 50 | 100 |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.85 | 1 | 0.2 | 0.306 | 0.419 | 0.520 | 0.634 | 0.734 | 0.825 | 0.911 | 0.952 |
| 0.85 | 1 | 0.5 | 0.088 | 0.163 | 0.281 | 0.462 | 0.602 | 0.740 | 0.870 | 0.930 |
| 0.85 | 1 | 1.0 | 0.020 | 0.041 | 0.092 | 0.238 | 0.426 | 0.616 | 0.811 | 0.903 |
| 0.85 | 0.999 | 0.2 | 0.306 | 0.420 | 0.511 | 0.631 | 0.727 | 0.808 | 0.876 | 0.898 |
| 0.85 | 0.999 | 0.5 | 0.087 | 0.165 | 0.280 | 0.455 | 0.594 | 0.717 | 0.823 | 0.852 |
| 0.85 | 0.999 | 1.0 | 0.021 | 0.041 | 0.092 | 0.235 | 0.416 | 0.592 | 0.747 | 0.782 |
| 0.85 | 0.99 | 0.2 | 0.301 | 0.403 | 0.490 | 0.589 | 0.653 | 0.703 | 0.714 | 0.708 |
| 0.85 | 0.99 | 0.5 | 0.084 | 0.159 | 0.263 | 0.416 | 0.512 | 0.583 | 0.604 | 0.593 |
| 0.85 | 0.99 | 1.0 | 0.019 | 0.041 | 0.085 | 0.211 | 0.347 | 0.438 | 0.474 | 0.471 |

For omega = 1, mu0 = 0.5: g = 0 gives R = 0.501, 0.656, 0.820, 0.900, 0.943 at tau 1, 2, 5, 10, 20; g = 0.7 gives 0.271, 0.411, 0.605, 0.741, 0.842, 0.929 at tau 1, 2, 5, 10, 20, 50. Fit for the conservative case when (1 - g) tau is at least about 1: R(mu0) = 1 - (mu0 + 0.71) / ((1 - g) tau + 1.42), with errors against the table for tau of 5 or more of at most 0.02 (mu0 = 0.5), 0.05 (mu0 = 0.2), 0.03 (mu0 = 1); it fails below that (R = 0 at small tau for mu0 = 1). Absorbing clouds saturate (omega = 0.99 near 0.60 at mu0 = 0.5); use the table. Kitzmann et al. 2013 (arXiv:1306.3940, READ abstract) find that two-stream schemes overestimate reflected and transmitted radiation for CO2 ice clouds relative to a discrete-ordinate solver, so prefer this table to a two-stream formula. 

## 4. Boss must research (not reached, or reached only second-hand)

1. Original Fray and Schmitt 2009 (Planet. Space Sci. 57, 2053): lower validity limits of every sublimation fit and the coefficients themselves (read through a SOM table of arXiv:2009.02277 with two typos); H2 and He sublimation. Needed for every cold ice.
2. Saturation or solubility data not found: H2SO4-H2O vapour pressure for 60 to 98 wt percent (Ayers 1980, Zeleznik 1991, Tabazadeh 1997; only the statement that vapour is below 1e-6 Pa at stratospheric temperatures was read, in Loftus 2019 citing Kulmala and Laaksonen 1990); NH4Cl and NH4HSO4; HNO3 and HCl hydrates; solid NaCl and KCl vapour (JANAF); Cr2O3, Fe2O3, Fe3O4 against the metal; carbonate and sulfate decomposition (CaCO3, MgCO3, CaSO4; NIST-JANAF); the radon sublimation branch.
3. Perchlorate brines: Ca(ClO4)2 and Mg(ClO4)2 eutectic temperatures, deliquescence humidity, water activity (Chevrier, Hanley, Gough and related papers). Needed for cold dry rocky surfaces.
4. Optical constants: numeric Mars dust spectrum and imaginary index (Wolff 2009, Ockert-Bell 1997, Lemmon 2004, Tomasko 1999); Titan haze absolute spectrum (Tomasko 2008, Doose 2016); visible k of CH4 and NH3 ice; H2SO4 at Venus temperatures; hydrated sea salt; perchlorates; carbonate k; the FeS row; the HITRAN-RI H2SO4, HNO3 and organic-acid files (unread).
5. Haze yield: lab mg/h to planetary mass flux (Hörst 2018 says impossible), production rates, K_zz(P) per class, coagulation outcome tables.
6. Cloud microphysics: A and M appendix B drag law; Powell 2018 and 2019 (silicate and titanium clouds) and Mbarek and Kempton 2016 (species by regime) were downloaded but their tables not extracted; haze grid.
7. Dust: Kok 2011 fragmentation size distribution (cited only), coefficients C_K and C_GP, wet or crusted soil thresholds, cohesion and electrostatics at low gravity, dust-devil lofting.
8. Lava worlds: partial pressures of Na, K, SiO, Mg, Fe, TiO, VO against temperature and composition (Schaefer, Lodders and Fegley 2012; LavAtmos arXiv:2210.10463), sodium and potassium escape.
9. Validated cloud radiative transfer with polarisation and non-spherical particles (the Monte Carlo here is scalar and spherical).
10. Equilibrium condensation for non-solar C/O, Mg/Si, Fe/Mg (FastChem or GGchem); the fits in 2.1 scale metallicity only.
11. Widening the proxy for doi.org, science.org, agupubs, aanda, sciencedirect, usgs and mindat would give the journal versions of the above.

## 5. Recalled / to verify

Each item here is RECALLED or an assumption of this note; none is presented as verified.

| Item | Value or statement | Used in | Where to verify |
|---|---|---|---|
| CO2 ice density | 1.56 g/cm3 | 2.8 optical table | Fray and Schmitt or a CO2 phase-diagram source |
| Liquid water density | 0.997 g/cm3 | 2.8 | any water property table (CRC) |
| Soot and amorphous carbon density | 1.8 g/cm3 | 2.8 | Maricq 2004, Rissler 2013 (cited in Adams 2019) |
| SiC density | about 3.2 g/cm3 | 2.8 | webmineral moissanite page |
| Saharan dust, volcanic, dust-like, sand densities | 2.6, 2.6, 2.6, 2.8 g/cm3; sea salt 2.2; ammonium sulfate 1.77; carbonaceous 1.8; meteoric 3.0; water-soluble 1.7 | 2.8 | no source |
| A(Mg) solar abundance | 7.55 (Mg/H 3.55e-5) | check of the Mg2SiO4 constant in 2.3 | Lodders 2010 (arXiv:1010.2746, downloaded, not extracted) |
| Clausius-Clapeyron fallback | ln P = ln P_tr - (L/R)(1/T - 1/T_tr) | 3.2 | standard thermodynamics text |
| Kinetic mean free path | lambda = eta / (0.499 rho_a v_mean) | 2.10, 3.6 | Jennings 1988 or Seinfeld and Pandis |
| Gamma size distribution | n(r) proportional to r^((1 - 3 b) / b) exp(-r / (a b)) with a = r_eff, b = v_eff | Mars inversion in 2.10 | Hansen and Travis 1974 |
| von Karman constant and roughness length | kappa = 0.4; z0 about D/30 | 3.8 | any boundary-layer text |
| Iversen-White cgs units | the constant 0.006 in cgs | 3.8, 2.10 | no source |
| FeS optical row | n 1.415, k 1.4 at 550 nm from LX-MIE (data, but physically doubtful) | 2.8 | Pollack 1994 |
| Mg2SiO4 constant in the Gao 2020 saturation expression | 14.88 is inferred (printed as 4.88; 14.88 matches the Visscher fit within 17 K at 1 bar) | 2.3 | Gao 2020 authors |
| CO2 ice table | n and k values read from OCR text of the PDF | 2.8 | the scanned table |
| Hematite colour direction | negative Angstrom exponent in all three hematite sets | 2.8 | spectra; Triaud |

## 6. Cross-links to neighbouring parts

- Part 05 (mixing ratios by planet class, file 05-mixing-ratios-by-class.md): supplies x_i for every condensable gas; this part returns the depleted profiles x_i(P) above cloud base (3.3) and the sulfur, sodium, potassium and carbon sinks (MnS, Na2S, KCl, NH4SH, haze).
- Gas table and thermodynamics part: triple points, critical points, enthalpies and the Antoine anchors in 2.5; molecular weights for epsilon in 3.5.
- Gas optics (Rayleigh, UV and visible scattering) part: cloud and haze tau, SSA and Angstrom exponent from 2.8 and 3.5 to 3.7 add to the Rayleigh optical depth of the mixture for atmosphere colour (GEN.414, GEN.417 to GEN.419); the project note atmosphere-colour.md says dust is grey, haze goes as lambda^-1 to lambda^-2.3, and cloud depth follows Venus (20 to 30) and Earth (about 10): compare with the Angstrom values computed here (small haze particles 3 to 4, 1 micron droplets near 0, dust 3 at 0.1 micron falling to negative at 1 micron).
- Infrared greenhouse and thermal part: cloud and haze infrared opacity, scattering greenhouse of CO2 ice (Kitzmann 2013), cloud-top temperatures; the cloud base rule here uses T(P) from that part.
- Photochemistry, ozone and escape part: precursor photolysis rates (CH4, HCN, C2H2), Lyman-alpha flux for the haze rule, H2SO4 formation from SO2 (Loftus 2019), UV shielding by haze and dust; escape of Na, K and lava-world species.
- Radiation shielding and ground sources part: UV and particle shielding by aerosol columns; radon (2.5) and mineral grain hosts (2.12) for ground radiation sources.
- Bulk composition chain (GEN.425 to GEN.430, docs/design/composition/bulk-composition.md): the mineral list and densities in 2.12 and the lava-world compositions in 2.11 feed the crust and dust composition; dust grain density draws use 2.12.
