# 01 Gas table: thermodynamics, phase data, vapour pressure curves, real-gas and mixture formulas

Part 01 of the v9.0 atmosphere and composition notes (research worker output; the repository was not edited).

## 0. Provenance and which domains worked

Worked: export.arxiv.org PDFs (arXiv 2009.02277 Lisse et al. 2021, 2309.05078 Grundy et al. 2024, astro-ph/0511136 Visscher et al. 2006, 2410.11138 Lodders and Fegley 2024, 1312.3750 Gao et al. 2014; 1104.2131, 2103.15978 and 2209.03513 were fetched but supply no number here); pypi.org and files.pythonhosted.org wheels (chemicals 1.5.2, thermo 0.6.1, periodictable 2.1.0, cantera data files nasa_gas.yaml and nasa_condensed.yaml); webbook.nist.gov (species pages and the fluid server: saturation tables, 1e-4 MPa isotherms, Z points); exomol.com partition-function files (*.pf); ntrs.nasa.gov (Blakley et al. 2024 (PSS 244, 105863, arXiv:2403.09917; NASA NTRS record dated 2025), equilibrium vapour pressures of NH3 and O2 ices). Refused or not reached: ScienceDirect (Fray and Schmitt 2009), NIST SRD 156, journals. Downloaded files were read as data only; Python reading them ran under -I.

Cp tags: J JANAF 1998 (chemicals), N NASA-7 (cantera), S Shomate, f NIST fluid at 1e-4 MPa, x NASA-7 extrapolated below its 200 K range (unreliable), C computed 5R/2, P Poling polynomial. Phase-point tags: h heos_constants (EOS values), e NIST EOS table, n NIST WebBook median, y Yaws, c CRC; "sub" marks a sublimation point at 101325 Pa. Row counts by label are in section 15.

## 1. Scope and generator feed

Input: species, T, total pressure P, partial pressures. Output: molar mass, saturation pressure or condensation temperature, ideal and real density, Cp and gamma, scale height, Henry solubility, viscosity and flags (condenses, absorbs IR, photolysed). The species are the 72 gases of T1 plus the extras of section 5. Section 7 replaces the two-point rows of `CONDENSATION_K` in `src/planetgen/physics/atmosphere.py` (lines 77 to 84; `vapour_pressure_kpa`, line 137, interpolates ln P between the temperature at 100 kPa and at 0.001 kPa). Solar System bodies appear only as labelled test points (section 12).

Flags (defined here, not sourced): IR column: monatomic none, homonuclear diatomics CIA only (collision-induced absorption), all others IR (greenhouse candidate; strengths are in the infrared part). Volatility class from the normal boiling or sublimation point: UV below 30 K, SV 30 to 150 K, V 150 to 350 K, SmV 350 to 1000 K, Ref above 1000 K; "condensible" on a body when the class is SV or V and x_i P exceeds the saturation curve. NUV = Y when the thermodynamic photolysis threshold is at or beyond 300 nm. This is a necessary condition only (threshold = lowest bond-energy channel from NASA enthalpies, lambda_nm = 119626.565 / dH in kJ/mol); H2, CH4 and H2O have long thresholds but absorb only in the vacuum UV. True fragility needs cross sections (part 03).

## 2. T1: identity, mass, flags, photolysis threshold

M (g/mol) is COMPUTED from periodictable 2.1.0 element masses (isotope masses for 3He, HD, HDO). CAS from NIST WebBook (DATA). Threshold: nm and products of the lowest channel (cantera nasa_gas.yaml enthalpies), COMPUTED.

| Species | CAS | M | IR | Vol | Threshold nm (channel) | NUV | Status |
|---|---|---|---|---|---|---|---|
| H2 | 1333-74-0 | 2.0160 | CIA | UV | 274 (H+H) | - | COMPUTED |
| HD | 13983-20-5 | 3.0219 | IR(weak) | ? | 272 (H+D) | - | COMPUTED |
| He | 7440-59-7 | 4.0026 | - | UV | - | - | COMPUTED |
| 3He | 14762-55-1 | 3.0160 | - | UV | - | - | COMPUTED |
| Ne | 7440-01-9 | 20.1797 | - | UV | - | - | COMPUTED |
| Ar | 7440-37-1 | 39.9500 | - | SV | - | - | COMPUTED |
| Kr | 7439-90-9 | 83.7980 | - | SV | - | - | COMPUTED |
| Xe | 7440-63-3 | 131.2930 | - | V | - | - | COMPUTED |
| Rn | 10043-92-2 | 222.0176 | - | V | - | - | COMPUTED |
| N2 | 7727-37-9 | 28.0140 | CIA | SV | 127 (N+N) | - | COMPUTED |
| O2 | 7782-44-7 | 31.9980 | CIA | SV | 240 (O+O) | - | COMPUTED |
| O3 | 10028-15-6 | 47.9970 | IR | V | 1114 (O2+O) | Y | COMPUTED |
| H2O | 7732-18-5 | 18.0150 | IR | SmV | 240 (H+OH) | - | COMPUTED |
| HDO | 14940-63-7 | 19.0209 | IR | ? | - | - | COMPUTED |
| H2O2 | 7722-84-1 | 34.0140 | IR | SmV | 556 (OH+OH) | Y | COMPUTED |
| CO | 630-08-0 | 28.0100 | IR | SV | 111 (C+O) | - | COMPUTED |
| CO2 | 124-38-9 | 44.0090 | IR | V | 225 (CO+O) | - | COMPUTED |
| CH4 | 74-82-8 | 16.0430 | IR | SV | 272 (CH3+H) | - | COMPUTED |
| C2H2 | 74-86-2 | 26.0380 | IR | V | 215 (C2H+H) | - | COMPUTED |
| C2H4 | 74-85-1 | 28.0540 | IR | V | 257 (C2H3+H) | - | COMPUTED |
| C2H6 | 74-84-0 | 30.0700 | IR | V | 317 (CH3+CH3) | Y | COMPUTED |
| C3H8 | 74-98-6 | 44.0970 | IR | V | 323 (CH3+C2H5) | Y | COMPUTED |
| C4H10 | 106-97-8 | 58.1240 | IR | V | - | - | COMPUTED |
| C6H6 | 71-43-2 | 78.1140 | IR | SmV | 253 (C6H5+H) | - | COMPUTED |
| HCN | 74-90-8 | 27.0260 | IR | V | 229 (H+CN) | - | COMPUTED |
| CH2O | 50-00-0 | 30.0260 | IR | V | 325 (H+HCO) | Y | COMPUTED |
| CH3OH | 67-56-1 | 32.0420 | IR | V | 309 (CH3+OH) | Y | COMPUTED |
| NH3 | 7664-41-7 | 17.0310 | IR | V | 264 (NH2+H) | - | COMPUTED |
| N2H4 | 302-01-2 | 32.0460 | IR | SmV | 423 (NH2+NH2) | Y | COMPUTED |
| N2O | 10024-97-2 | 44.0130 | IR | V | 714 (N2+O) | Y | COMPUTED |
| NO | 10102-43-9 | 30.0060 | IR | SV | 190 (N+O) | - | COMPUTED |
| NO2 | 10102-44-0 | 46.0050 | IR | V | 391 (NO+O) | Y | COMPUTED |
| N2O4 | 10544-72-6 | 92.0100 | IR | V | - | - | COMPUTED |
| HNO3 | 7697-37-2 | 63.0120 | IR | SmV | 577 (OH+NO2) | Y | COMPUTED |
| HNCO | 75-13-8 | 43.0250 | IR | ? | 241 (H+NCO) | - | COMPUTED |
| C3O2 | 504-64-3 | 68.0310 | IR | V | 436 (C2O+CO) | Y | COMPUTED |
| CN | 2074-87-5 | 26.0180 | IR | ? | 160 (C+N) | - | COMPUTED |
| C2N2 | 460-19-5 | 52.0360 | IR | V | 209 (CN+CN) | - | COMPUTED |
| CH3CN | 75-05-8 | 41.0530 | IR | SmV | 229 (CH3+CN) | - | COMPUTED |
| HC3N | 1070-71-9 | 51.0480 | IR | V | - | - | COMPUTED |
| C4H2 | 460-12-8 | 50.0600 | IR | V | - | - | COMPUTED |
| PH3 | 7803-51-2 | 33.9978 | IR | V | - | - | COMPUTED |
| H2S | 7783-06-4 | 34.0760 | IR | V | 317 (SH+H) | Y | COMPUTED |
| SO2 | 7446-09-5 | 64.0580 | IR | V | 217 (SO+O) | - | COMPUTED |
| SO3 | 7446-11-9 | 80.0570 | IR | V | 344 (SO2+O) | Y | COMPUTED |
| H2SO4 | 7664-93-9 | 98.0720 | IR | SmV | - | - | COMPUTED |
| OCS | 463-58-1 | 60.0700 | IR | V | 392 (CO+S) | Y | COMPUTED |
| CS2 | 75-15-0 | 76.1310 | IR | V | 272 (CS+S) | - | COMPUTED |
| S2 | 23550-45-0 | 64.1200 | CIA | ? | 281 (S+S) | - | COMPUTED |
| S8 | 10544-50-0 | 256.4800 | IR | ? | - | - | COMPUTED |
| SO | 13827-32-2 | 48.0590 | IR | ? | 230 (S+O) | - | COMPUTED |
| CH3SH | 74-93-1 | 48.1030 | IR | V | - | - | COMPUTED |
| (CH3)2S | 75-18-3 | 62.1300 | IR | V | - | - | COMPUTED |
| HCl | 7647-01-0 | 36.4580 | IR | V | 277 (H+Cl) | - | COMPUTED |
| HF | 7664-39-3 | 20.0064 | IR | V | 210 (H+F) | - | COMPUTED |
| Cl2 | 7782-50-5 | 70.9000 | CIA | V | 493 (Cl+Cl) | Y | COMPUTED |
| CH3Cl | 74-87-3 | 50.4850 | IR | V | 342 (CH3+Cl) | Y | COMPUTED |
| HBr | 10035-10-6 | 80.9120 | IR | V | 327 (H+Br) | Y | COMPUTED |
| HClO4 | 7601-90-3 | 100.4540 | IR | SmV | - | - | COMPUTED |
| CF4 | 75-73-0 | 88.0046 | IR | SV | 221 (CF3+F) | - | COMPUTED |
| SF6 | 2551-62-4 | 146.0504 | IR | V | 306 (SF5+F) | Y | COMPUTED |
| SiH4 | 7803-62-5 | 32.1170 | IR | V | 304 (SiH3+H) | Y | COMPUTED |
| SiO | 10097-28-6 | 44.0840 | IR | ? | 150 (Si+O) | - | COMPUTED |
| Na | 7440-23-5 | 22.9898 | - | Ref | - | - | COMPUTED |
| K | 7440-09-7 | 39.0983 | - | Ref | - | - | COMPUTED |
| Ca | 7440-70-2 | 40.0780 | - | Ref | - | - | COMPUTED |
| Mg | 7439-95-4 | 24.3050 | - | Ref | - | - | COMPUTED |
| Fe | 7439-89-6 | 55.8450 | - | Ref | - | - | COMPUTED |
| TiO | 12137-20-1 | 63.8660 | IR | ? | 179 (Ti+O) | - | COMPUTED |
| NaCl | 7647-14-5 | 58.4398 | IR | Ref | 292 (Na+Cl) | - | COMPUTED |
| KCl | 7447-40-7 | 74.5483 | IR | Ref | 281 (K+Cl) | - | COMPUTED |
| NH4Cl | 12125-02-9 | 53.4890 | IR | SmV | - | - | COMPUTED |

## 3. T2: phase points and critical constants

Tt (K), Pt (Pa), Tb (K), Tc (K) [min-max over sources / n sources, outliers beyond 15% dropped], Pc (MPa) likewise, acentric factor omega. Where sources disagree, draw within the spread. Warnings: single-source metal and salt Tc, Pc are estimates (Mg Pc 280 MPa, Fe Pc 1015 MPa and Tc 9340 K are implausible as measurements); NH4Cl omega 3.92 is meaningless (it dissociates); HClO4 and H2SO4 values are Yaws estimates for compounds that decompose on boiling; 3He Tb is Yaws; HD Tt is a NIST median; Rn Tt is NOT REACHED.

| Species | Tt K | Pt Pa | Tb K | Tc K [range/n] | Pc MPa [range/n] | omega | Status |
|---|---|---|---|---|---|---|---|
| H2 | 13.957h | 7358e | 20.369h | 33.15 [33.15-33.20/5] | 1.296 [1.296-1.313/5] | -0.219 | DATA |
| HD | 16.600n | -- | -- | 36.00 [36.00-36.00/1] | 1.479 [1.479-1.479/1] | - | DATA |
| He | 2.1768h | 5039e | 4.2238h | 5.20 [5.19-5.30/5] | 0.228 [0.227-0.232/5] | -0.383 | DATA |
| 3He | -- | -- | 3.2000y | 3.31 [3.31-3.31/2] | 0.114 [0.114-0.117/2] | -0.471 | DATA |
| Ne | 24.556h | 43355e | 27.100h | 44.40 [44.40-44.40/5] | 2.662 [2.653-2.756/5] | -0.036 | DATA |
| Ar | 83.806h | 68891e | 87.302h | 150.69 [150.65-150.86/8] | 4.863 [4.855-4.898/7-1out] | -0.002 | DATA |
| Kr | 115.78h | 73534e | 119.73h | 209.48 [209.35-209.48/6] | 5.525 [5.502-5.525/6] | -0.001 | DATA |
| Xe | 161.41h | 81772e | 165.05h | 289.73 [289.70-289.74/6] | 5.842 [5.836-5.842/5] | 0.004 | DATA |
| Rn | -- | -- | 211.45y | 377.00 [377.00-377.50/3] | 6.282 [6.280-6.282/3] | - | DATA |
| N2 | 63.151h | 12520e | 77.355h | 126.19 [126.19-128.45/9] | 3.396 [3.394-3.460/7-2out] | 0.037 | DATA |
| O2 | 54.361h | 146.3e | 90.188h | 154.58 [154.58-155.15/8] | 5.043 [4.985-5.046/9] | 0.022 | DATA |
| O3 | -- | -- | 161.30n | 261.10 [261.00-261.15/6] | 5.573 [5.532-5.573/4] | 0.212 | DATA |
| H2O | 273.16h | 611.7e | 373.12h | 647.10 [647.00-647.30/6] | 22.064 [22.046-22.064/10] | 0.344 | DATA |
| HDO | -- | -- | -- | - [-] | - [-] | - | NOT REACHED |
| H2O2 | -- | -- | 423.35y | 728.00 [728.00-730.15/3] | 22.000 [21.700-22.000/3] | 0.358 | DATA |
| CO | 68.160h | 15537e | 81.638h | 132.86 [132.86-134.45/6] | 3.494 [3.475-3.521/9] | 0.050 | DATA |
| CO2 | 216.59h | 5.180e5e | 194.67y sub | 304.13 [304.10-304.35/10] | 7.377 [7.340-7.383/8] | 0.225 | DATA |
| CH4 | 90.694h | 11696e | 111.67h | 190.56 [190.56-190.60/8] | 4.599 [4.599-4.610/8] | 0.011 | DATA |
| C2H2 | 191.75h | 1.235e5h | 189.30n sub | 308.30 [308.30-309.70/11] | 5.988 [5.988-6.240/8] | 0.191 | DATA |
| C2H4 | 103.99h | 122.0e | 169.38h | 282.35 [282.34-282.50/8] | 5.042 [5.032-5.060/8] | 0.087 | DATA |
| C2H6 | 90.368h | 1.142e | 184.57h | 305.32 [305.30-305.42/8] | 4.872 [4.872-4.900/8] | 0.100 | DATA |
| C3H8 | 85.525h | 1.724e-4e | 231.04h | 369.89 [369.82-369.95/8] | 4.251 [4.246-4.251/8] | 0.153 | DATA |
| C4H10 | 134.90h | 0.6657e | 272.66h | 425.12 [425.00-425.20/8] | 3.796 [3.790-3.800/8] | 0.201 | DATA |
| C6H6 | 278.67h | 4785e | 353.22h | 562.02 [562.00-562.16/8] | 4.907 [4.890-4.907/8] | 0.211 | DATA |
| HCN | 259.86n | -- | 299.00n | 456.80 [456.65-457.00/5] | 5.390 [5.390-5.400/5] | 0.410 | DATA |
| CH2O | 155.10n | -- | 253.85n | 408.00 [408.00-408.00/2] | 6.586 [6.586-6.586/2] | 0.282 | DATA |
| CH3OH | 175.61h | 0.1864e | 337.63h | 513.38 [512.50-513.38/7] | 8.216 [8.010-8.216/7] | 0.561 | DATA |
| NH3 | 195.49h | 6053e | 239.83h | 405.56 [405.40-405.65/7] | 11.363 [11.277-11.363/6] | 0.256 | DATA |
| N2H4 | 274.62n | -- | 386.65n | 653.00 [653.00-653.15/4] | 14.692 [14.692-14.700/4] | 0.314 | DATA |
| N2O | 182.33h | 87837e | 184.68h | 309.52 [309.49-309.65/9] | 7.245 [7.235-7.270/8] | 0.162 | DATA |
| NO | -- | -- | 121.38y | 180.00 [180.00-180.15/4] | 6.485 [6.480-6.485/4] | 0.583 | DATA |
| NO2 | -- | -- | 295.08n | 431.15 [431.00-431.15/2] | 10.133 [10.133-10.133/2] | 0.851 | DATA |
| N2O4 | -- | -- | 294.31y | 431.00 [431.00-431.15/3] | 10.133 [10.133-10.133/3] | 0.853 | DATA |
| HNO3 | -- | -- | 356.15y | 520.00 [520.00-520.00/2] | 6.890 [6.890-6.890/2] | 0.714 | DATA |
| HNCO | -- | -- | -- | - [-] | - [-] | - | NOT REACHED |
| C3O2 | -- | -- | 279.95y | 431.80 [431.80-431.80/1] | 5.119 [5.119-5.119/1] | - | DATA |
| CN | -- | -- | -- | - [-] | - [-] | - | NOT REACHED |
| C2N2 | 245.32n | 75100n | 252.00n | 400.00 [397.00-400.15/6] | 5.978 [5.924-6.303/6] | 0.312 | DATA |
| CH3CN | 229.32n | -- | 354.80n | 545.00 [545.00-548.00/6] | 4.870 [4.833-4.880/6] | 0.338 | DATA |
| HC3N | -- | -- | 315.70n | - [-] | - [-] | - | NOT REACHED |
| C4H2 | -- | -- | 284.00n | 478.02 [478.02-478.02/1] | 5.863 [5.863-5.863/1] | 0.100 | DATA |
| PH3 | 139.41n | -- | 185.41y | 324.80 [324.75-324.80/3] | 6.535 [6.535-6.540/3] | 0.045 | DATA |
| H2S | 187.70h | 23259e | 212.85h | 373.10 [372.80-373.53/7] | 9.000 [8.937-9.000/7] | 0.100 | DATA |
| SO2 | 197.70h | 1666e | 263.14h | 430.64 [430.34-430.80/6] | 7.887 [7.883-7.887/5] | 0.256 | DATA |
| SO3 | -- | -- | 317.90y | 491.00 [490.85-491.00/4] | 8.207 [8.207-8.210/4] | 0.424 | DATA |
| H2SO4 | -- | -- | 610.00y | 924.00 [924.00-924.00/1] | 6.400 [6.400-6.400/1] | 0.494 | DATA |
| OCS | 134.30h | 64.43e | 222.99h | 378.77 [375.00-378.80/6] | 6.370 [5.877-6.370/6] | 0.097 | DATA |
| CS2 | 161.35n | -- | 319.20n | 552.00 [552.00-552.49/5] | 7.903 [7.329-7.903/4] | 0.111 | DATA |
| S2 | -- | -- | -- | - [-] | - [-] | - | NOT REACHED |
| S8 | -- | -- | -- | - [-] | - [-] | - | NOT REACHED |
| SO | -- | -- | -- | - [-] | - [-] | - | NOT REACHED |
| CH3SH | 150.15n | -- | 279.10n | 469.90 [469.90-470.00/7] | 7.230 [7.230-7.240/6] | 0.158 | DATA |
| (CH3)2S | 174.85n | -- | 311.00n | 503.00 [503.00-503.04/8] | 5.530 [5.400-5.532/6] | 0.194 | DATA |
| HCl | 159.07h | 13828h | 188.17h | 324.68 [324.60-324.68/6] | 8.313 [8.256-8.313/6] | 0.132 | DATA |
| HF | -- | -- | 292.70n | 461.00 [461.00-461.15/3] | 6.485 [6.480-6.485/3] | 0.382 | DATA |
| Cl2 | 172.17h | 1381h | 239.20h | 416.87 [416.87-417.15/7] | 7.642 [7.642-7.991/7] | 0.069 | DATA |
| CH3Cl | 175.51h | 870.0n | 249.17h | 416.30 [416.00-416.30/7] | 6.690 [6.677-7.400/10] | 0.153 | DATA |
| HBr | -- | -- | 206.45y | 363.20 [363.15-363.20/4] | 8.552 [8.552-8.552/4] | 0.069 | DATA |
| HClO4 | -- | -- | 385.00y | 631.00 [631.00-631.00/1] | 3.860 [3.860-3.860/1] | 0.050 | DATA |
| CF4 | 120.00h | 11294e | 145.10h | 227.51 [227.50-227.54/8] | 3.750 [3.680-3.750/7] | 0.181 | DATA |
| SF6 | 223.56h | 2.314e5e | 209.30n sub | 318.72 [318.69-318.76/8-2out] | 3.755 [3.755-3.760/6-1out] | 0.218 | DATA |
| SiH4 | -- | -- | 161.00y | 269.70 [269.00-269.70/3] | 4.840 [4.840-4.855/3] | 0.094 | DATA |
| SiO | -- | -- | -- | - [-] | - [-] | - | NOT REACHED |
| Na | 370.98n | -- | 1156.0n | 2573.00 [2573.00-2573.00/1] | 35.464 [35.464-35.464/1] | -0.105 | DATA |
| K | 336.35n | -- | 1037.0n | 2223.00 [2223.00-2223.00/1] | 16.000 [16.000-16.000/1] | -0.203 | DATA |
| Ca | 1112.9n | -- | 1757.0n | 3267.00 [3267.00-3267.00/1] | 15.300 [15.300-15.300/1] | 0.117 | DATA |
| Mg | 922.00n | -- | 1363.0n | 2241.04 [2241.04-2241.04/1] | 280.240 [280.240-280.240/1] | - | DATA |
| Fe | -- | -- | 3133.3y | 9340.00 [9340.00-9340.00/1] | 1015.000 [1015.000-1015.000/1] | -0.011 | DATA |
| TiO | -- | -- | -- | - [-] | - [-] | - | NOT REACHED |
| NaCl | -- | -- | 1738.2y | 3400.00 [3400.00-3400.00/1] | 35.500 [35.500-35.500/1] | 0.189 | DATA |
| KCl | -- | -- | 1689.0y | 3470.00 [3470.00-3470.00/1] | 18.000 [18.000-18.000/1] | -0.121 | DATA |
| NH4Cl | -- | -- | 612.00y | 882.00 [882.00-882.00/1] | 1.640 [1.640-1.640/1] | 3.920 | DATA |

Sources: chemicals 1.5.2 (Misc/heos_constants.tsv, Critical Properties, Triple Properties, Vapor Pressure), NIST fluid saturation tables, NIST WebBook phase-change pages.

## 4. T3: heat capacity, gamma, enthalpies

Ideal-gas Cp (J/mol/K) at 100, 300, 1000 K; gamma = Cp/(Cp - R) at 300 K (R = 8.31446261815 J/mol/K); spread = percent spread among other reachable sources (JANAF, NASA-7, Shomate) at 300 K / 1000 K; dHvap (kJ/mol) at Tb (e NIST EOS table, c CRC); dHsub = NIST WebBook range. Cautions: x values are extrapolations; C6H6 at 100 K is withheld (unphysical extrapolation); C2H6 and C3H8 at 100 K are x values because the NIST 1e-4 MPa isotherm returns liquid there; H2 at 100 K depends on the ortho-para state (NIST 22.57, NASA-7 extrapolation 24.67, JANAF not read), so treat H2 Cp below 250 K as uncertain by about 10%. Monatomic gases are 5R/2 = 20.786 and gamma = 5/3 (COMPUTED).

| Species | Cp 100 | Cp 300 | Cp 1000 | gamma 300 | spread % | dHvap | dHsub | Status |
|---|---|---|---|---|---|---|---|---|
| H2 | 22.57f | 28.85N | 30.16N | 1.405 | 0.0/0.1 | 0.90e | - | DATA |
| HD | 29.29J | 29.20J | 30.71J | 1.398 | 0.0/0.0 | - | - | DATA |
| He | 20.786C | 20.786C | 20.786C | 1.667 | 0.0/0.0 | 0.08e | - | COMPUTED |
| 3He | 20.786C | 20.786C | 20.786C | 1.667 | -/- | - | - | COMPUTED |
| Ne | 20.786C | 20.786C | 20.786C | 1.667 | 0.0/0.0 | 1.73e | - | COMPUTED |
| Ar | 20.786C | 20.786C | 20.786C | 1.667 | 0.0/0.0 | 6.44e | - | COMPUTED |
| Kr | 20.786C | 20.786C | 20.786C | 1.667 | 0.0/0.0 | 8.97e | - | COMPUTED |
| Xe | 20.786C | 20.786C | 20.786C | 1.667 | 0.0/0.0 | 12.55e | - | COMPUTED |
| Rn | 20.786C | 20.786C | 20.786C | 1.667 | -/- | - | - | COMPUTED |
| N2 | 29.10f | 29.13N | 32.68N | 1.400 | 0.0/0.0 | 5.58e | - | DATA |
| O2 | 29.11f | 29.39N | 34.88N | 1.395 | 0.0/0.1 | 6.82e | - | DATA |
| O3 | 33.29J | 39.33J | 55.02J | 1.268 | 0.3/1.9 | - | - | DATA |
| H2O | 33.30J | 33.60J | 41.27J | 1.329 | 0.0/0.1 | 40.65e | - | DATA |
| HDO | 33.30J | 33.80J | 43.28J | 1.326 | 0.0/0.0 | - | - | DATA |
| H2O2 | 33.53J | 43.22J | 62.84J | 1.238 | 1.8/0.8 | - | - | DATA |
| CO | 29.10f | 29.14J | 33.18J | 1.399 | 0.0/0.1 | 6.01e | 7.6-8.1 | DATA |
| CO2 | 29.21J | 37.22J | 54.31J | 1.288 | 0.0/0.0 | - | 25.2-27.2 | DATA |
| CH4 | 33.28f | 35.71J | 71.80J | 1.304 | 0.1/2.5 | 8.19e | 0.1-10.0 | DATA |
| C2H2 | 29.35J | 44.23J | 68.28J | 1.232 | 0.2/0.5 | - | 21.8-25.2 | DATA |
| C2H4 | 33.27J | 43.06J | 93.90J | 1.239 | 0.1/0.3 | 13.53e | 15.0-18.4 | DATA |
| C2H6 | 35.52x | 52.72N | 122.69N | 1.187 | -/- | 14.72e | 20.5-22.6 | DATA (x) |
| C3H8 | 41.58x | 73.95N | 174.84N | 1.127 | -/- | 18.77e | 28.5 | DATA (x) |
| C4H10 | 58.29x | 99.11N | 227.75N | 1.092 | -/- | 22.42e | 35.9 | DATA (x) |
| C6H6 | - | 82.48N | 211.32N | 1.112 | -/- | 30.75e | 38.0-53.9 | DATA |
| HCN | 29.17J | 35.93J | 50.23J | 1.301 | 0.0/0.0 | - | 35.6-37.6 | DATA |
| CH2O | 33.26J | 35.46J | 61.99J | 1.306 | 0.1/0.0 | - | - | DATA |
| CH3OH | 39.71x | 44.15N | 89.66N | 1.232 | -/- | 35.28e | - | DATA (x) |
| NH3 | 33.28J | 35.70J | 56.49J | 1.304 | 0.1/0.4 | 23.33e | 31.2 | DATA |
| N2H4 | 34.37J | 51.02J | 96.36J | 1.195 | 5.0/4.3 | 41.80c | - | DATA |
| N2O | 29.35J | 38.70J | 54.87J | 1.274 | 0.0/0.7 | 16.47e | 23.6-25.1 | DATA |
| NO | 32.30J | 29.84J | 33.99J | 1.386 | 0.1/0.0 | 13.83c | - | DATA |
| NO2 | 33.28J | 37.03J | 52.17J | 1.290 | 0.6/1.6 | - | - | DATA |
| N2O4 | 46.26J | 77.49J | 119.21J | 1.120 | 2.5/0.5 | 38.12c | - | DATA |
| HNO3 | 33.92J | 53.53J | 90.42J | 1.184 | 1.5/1.6 | - | - | DATA |
| HNCO | 33.47J | 45.16J | 67.52J | 1.226 | 0.0/0.1 | - | - | DATA |
| C3O2 | 45.31J | 67.16J | 97.64J | 1.141 | 0.6/0.1 | - | - | DATA |
| CN | 29.10J | 29.16J | 33.46J | 1.399 | 0.0/0.1 | - | - | DATA |
| C2N2 | 36.40J | 56.86J | 76.30J | 1.171 | 0.6/0.7 | 23.33c | 32.4-34.4 | DATA |
| CH3CN | 36.77x | 52.39N | 98.42N | 1.189 | -/- | 29.75c | - | DATA (x) |
| HC3N | - | - | - | - | -/- | - | - | NOT REACHED |
| C4H2 | 32.91x | 73.92N | 111.25N | 1.127 | -/- | - | 36.2-36.3 | DATA (x) |
| PH3 | 33.26J | 37.18J | 64.30J | 1.288 | 0.1/0.0 | 14.60c | 17.2 | DATA |
| H2S | 33.26J | 34.21J | 45.79J | 1.321 | 0.0/0.1 | 18.62e | 22.5-25.4 | DATA |
| SO2 | 33.53J | 39.95J | 54.48J | 1.263 | 0.0/0.0 | 24.96e | - | DATA |
| SO3 | 34.08J | 50.80J | 75.97J | 1.196 | 0.1/0.0 | 40.69c | - | DATA |
| H2SO4 | 42.10J | 84.05J | 132.69J | 1.110 | 0.0/0.1 | - | - | DATA |
| OCS | 29.61J | 41.59J | 56.99J | 1.250 | 0.1/0.0 | 18.56e | - | DATA |
| CS2 | 31.00J | 45.75J | 59.28J | 1.222 | 0.0/0.0 | 26.74c | - | DATA |
| S2 | 29.37J | 32.52J | 37.28J | 1.343 | 0.2/0.3 | - | - | DATA |
| S8 | 87.16J | 156.31J | 180.06J | 1.056 | 0.2/0.2 | - | - | DATA |
| SO | 29.11J | 30.20J | 36.05J | 1.380 | 0.0/0.1 | - | - | DATA |
| CH3SH | 37.24P | 50.51P | 93.97P | 1.197 | -/- | - | - | DATA |
| (CH3)2S | - | 74.38P | 146.87P | 1.126 | -/- | 27.00c | - | DATA |
| HCl | 29.12J | 29.14J | 31.63J | 1.399 | 0.0/0.0 | 16.15c | 19.6-19.7 | DATA |
| HF | 29.13J | 29.14J | 30.17J | 1.399 | 0.0/0.0 | - | - | DATA |
| Cl2 | 28.16x | 33.98N | 37.43N | 1.324 | 0.0/0.0 | 20.41c | - | DATA (x) |
| CH3Cl | 33.28J | 40.86J | 78.84J | 1.255 | 0.0/0.2 | 21.40c | 31.6 | DATA |
| HBr | 29.11J | 29.14J | 32.32J | 1.399 | 0.0/0.0 | - | - | DATA |
| HClO4 | - | - | - | - | -/- | - | - | NOT REACHED |
| CF4 | 34.74J | 61.29J | 98.80J | 1.157 | 0.0/0.0 | 11.83e | 14.0-17.0 | DATA |
| SF6 | 38.72J | 97.40J | 149.27J | 1.093 | 0.0/0.0 | - | 23.2-23.3 | DATA |
| SiH4 | 33.27J | 42.98J | 84.51J | 1.240 | 0.0/0.0 | 12.10c | - | DATA |
| SiO | 29.11J | 29.92J | 35.72J | 1.385 | 0.0/0.0 | - | - | DATA |
| Na | 20.79J | 20.79J | 20.79J | 1.667 | 0.0/0.0 | - | - | DATA |
| K | 20.79J | 20.79J | 20.79J | 1.667 | 0.0/0.0 | - | - | DATA |
| Ca | 20.79J | 20.79J | 20.79J | 1.667 | 0.0/0.0 | - | - | DATA |
| Mg | 20.79J | 20.79J | 20.79J | 1.667 | 0.0/0.0 | - | - | DATA |
| Fe | 21.39J | 25.68J | 22.49J | 1.479 | 0.0/0.2 | - | - | DATA |
| TiO | 33.88J | 32.50J | 37.74J | 1.344 | 0.0/0.0 | - | - | DATA |
| NaCl | 30.39J | 35.81J | 37.92J | 1.302 | 0.1/0.0 | - | - | DATA |
| KCl | 31.70J | 36.51J | 38.06J | 1.295 | 0.0/0.0 | - | - | DATA |
| NH4Cl | - | - | - | - | -/- | - | 168.6-176.6 | NOT REACHED |

## 5. Extras: refractory and hot-atmosphere species

Cp from ExoMol partition functions Q(T) (exomol.com *.pf), COMPUTED: Cp = 5R/2 + R (2 T dlnQ/dT + T^2 d2lnQ/dT2), finite differences of 5 K.

The twenty-five rows from H to HCOOH (added with GEN.454 and GEN.499) are the atoms, radicals and oxides that a hot or lava-world gas holds in equilibrium (`physics/gas_equilibrium.py`), so that a mixture drawn there is in the table. Their Cp is the NASA Glenn 9-coefficient fit of `physics/species_thermo.py` (McBride, Zehe and Gordon 2002; the fits start at 200 or 300 K, so Cp 100 is not given), COMPUTED; their molar mass is the element table's.

| Species | M g/mol | Cp 100 | Cp 300 | Cp 1000 | Phase or vapour data | Status |
|---|---|---|---|---|---|---|
| FeH | 56.853 | 33.19 | 36.33 | 43.90 | NOT REACHED (formation Gibbs data) | COMPUTED |
| CrH | 53.004 | 27.00 | 29.70 | 35.01 | NOT REACHED | COMPUTED |
| AlH | 27.990 | 29.13 | 29.38 | 35.06 | NOT REACHED | COMPUTED |
| VO | 66.941 | 29.14 | 30.73 | 36.49 | NOT REACHED | COMPUTED |
| CaOH | 57.085 | 31.93 | 46.85 | 53.59 | NOT REACHED | COMPUTED |
| TiO | 63.866 | 33.89 | 31.88 | 36.98 | JANAF Cp in T3 agrees within 4% at 300 K; vapour curve in T6 | COMPUTED |
| NH4SH | 51.107 | - | - | - | condensation only, section 7.4 | NOT REACHED (Cp) |
| H | 1.008 | - | 20.79 | 20.79 | no phase data; free energy in species_thermo.py | COMPUTED |
| O | 15.999 | - | 21.90 | 20.92 | no phase data; free energy in species_thermo.py | COMPUTED |
| OH | 17.007 | - | 29.88 | 30.68 | no phase data; free energy in species_thermo.py | COMPUTED |
| C | 12.011 | - | 20.84 | 20.79 | no phase data; free energy in species_thermo.py | COMPUTED |
| N | 14.007 | - | 20.79 | 20.79 | no phase data; free energy in species_thermo.py | COMPUTED |
| S | 32.060 | - | 23.67 | 21.49 | no phase data; free energy in species_thermo.py | COMPUTED |
| Si | 28.085 | - | 22.23 | 20.97 | no phase data; free energy in species_thermo.py | COMPUTED |
| Ti | 47.867 | - | 24.40 | 21.32 | no phase data; free energy in species_thermo.py | COMPUTED |
| V | 50.942 | - | 25.98 | 25.24 | no phase data; free energy in species_thermo.py | COMPUTED |
| Cr | 51.996 | - | 20.79 | 20.84 | no phase data; free energy in species_thermo.py | COMPUTED |
| P2 | 61.948 | - | 32.07 | 36.83 | no phase data; free energy in species_thermo.py | COMPUTED |
| FeO | 71.844 | - | 31.44 | 36.69 | no phase data; free energy in species_thermo.py | COMPUTED |
| MgO | 40.304 | - | 32.15 | 54.74 | no phase data; free energy in species_thermo.py | COMPUTED |
| SiO2 | 60.083 | - | 44.14 | 58.56 | no phase data; free energy in species_thermo.py | COMPUTED |
| TiO2 | 79.865 | - | 44.16 | 55.71 | no phase data; free energy in species_thermo.py | COMPUTED |
| NaOH | 39.997 | - | 48.03 | 54.54 | no phase data; free energy in species_thermo.py | COMPUTED |
| KOH | 56.105 | - | 49.24 | 54.71 | no phase data; free energy in species_thermo.py | COMPUTED |
| Na2 | 45.980 | - | 37.58 | 39.39 | no phase data; free energy in species_thermo.py | COMPUTED |
| NaO | 38.989 | - | 35.11 | 39.72 | no phase data; free energy in species_thermo.py | COMPUTED |
| MgOH | 41.312 | - | 46.56 | 54.28 | no phase data; free energy in species_thermo.py | COMPUTED |
| Mg(OH)2 | 58.319 | - | 80.80 | 96.12 | no phase data; free energy in species_thermo.py | COMPUTED |
| Fe(OH)2 | 89.859 | - | 71.74 | 98.26 | no phase data; free energy in species_thermo.py | COMPUTED |
| CaO | 56.077 | - | 32.49 | 37.25 | no phase data; free energy in species_thermo.py | COMPUTED |
| SiS | 60.145 | - | 32.31 | 36.91 | no phase data; free energy in species_thermo.py | COMPUTED |
| HCOOH | 46.025 | - | 45.84 | 89.40 | no phase data; free energy in species_thermo.py | COMPUTED |
| Al | 26.982 | - | 21.38 | 20.84 | no phase data; free energy in species_thermo.py | COMPUTED |
| K2 | 78.196 | - | 37.99 | 37.05 | no phase data; free energy in species_thermo.py | COMPUTED |
| KO | 55.097 | - | 35.38 | 38.09 | no phase data; free energy in species_thermo.py | COMPUTED |
| K2O | 94.195 | - | 54.22 | 57.76 | no phase data; free energy in species_thermo.py | COMPUTED |
| Na2O | 61.979 | - | 56.82 | 61.71 | no phase data; free energy in species_thermo.py | COMPUTED |
| AlO | 42.981 | - | 30.91 | 37.14 | no phase data; free energy in species_thermo.py | COMPUTED |
| Al2O | 69.963 | - | 52.04 | 60.65 | no phase data; free energy in species_thermo.py | COMPUTED |
| AlO2 | 58.980 | - | 51.73 | 60.70 | no phase data; free energy in species_thermo.py | COMPUTED |
| Al2O2 | 85.962 | - | 68.23 | 83.80 | no phase data; free energy in species_thermo.py | COMPUTED |

## 6. Saturation vapour pressure correlations

### 6.1 T4: best published correlation per species, tested against Tb and Tt

All candidates in chemicals (NIST Antoine, Poling Antoine, Landolt-Bornstein Antoine, Perry DIPPR101, Alcock and Wagner forms) were evaluated at Tb and Tt; the table shows the best by error at Tb (solid fits penalised). Forms: A log10(P/bar) = a - b/(T + c); a log10(P/Pa) = a - b/(T + c); L ln(P/Pa) = a - b/(T + c); P ln(P/Pa) = a + b/T + c ln T + d T^e; suffix S = sublimation. Error columns: pressure error at Tb against 101325 Pa (percent, with the equivalent temperature shift) and at Tt against Pt. Several fits (Xe, K, Rn, NO2, N2O4, C4H2) sit at or beyond their published range at Tb and must not be used there; N2O4 and NO2 interconvert and the tabulated curve belongs to the mixture.

| Species | Form | Coefficients | Range K | Err at Tb % (K) | Err at Tt % | Source | Status |
|---|---|---|---|---|---|---|---|
| H2 | A | 3.5431, 99.395, 7.7260 | 21.01-32.27 | -0.09 (+0.00 K) | - | NIST van Itterbeek, Verbeke | DATA |
| He | P | 11.533, -8.9900, 0.67240, 0.27430, 1.0000 | 1.760-5.200 | +0.59 (-0.01 K) | -0.18 | Perry | DATA |
| 3He | a | 6.3975, 5.5940, 0.69000 | 1.120-4.410 | -10.10 (+0.13 K) | - | Poling | DATA |
| Ne | P | 29.755, -271.06, -2.6081, 5.2700e-4, 2.0000 | 24.56-44.40 | +0.82 (-0.03 K) | +0.88 | Perry | DATA |
| Ar | P | 42.127, -1093.1, -4.1425, 5.7254e-5, 2.0000 | 83.78-150.9 | +0.20 (-0.02 K) | +0.05 | Perry | DATA |
| Kr | a | 8.7556, 416.38, -8.7000 | 111.3-129.2 | -0.02 (+0.00 K) | +0.10 | Poling | DATA |
| Xe | A | 2.8388, 326.60, -49.796 | 161.4-162.6 | -0.13 (+0.02 K) | -0.01 | NIST Chen, Lim, et al., 197 | DATA |
| Rn | a | 9.6204, 884.41, -18.150 | 158.0-222.9 | +9.49 (-1.65 K) | - | Poling | DATA |
| N2 | L | 20.024, 611.76, -5.3660 | 63.50-78.00 | -0.00 (+0.00 K) | +0.21 | Landolt | DATA |
| O2 | L | 20.284, 732.89, -6.5060 | 72.00-93.00 | -0.00 (+0.00 K) | - | Landolt | DATA |
| O3 | A | 4.2364, 712.49, 6.9820 | 92.80-162.0 | -0.74 (+0.13 K) | - | NIST Stull, 1947 | DATA |
| H2O | A | 5.0835, 1663.1, -45.622 | 344.0-373.0 | -0.09 (+0.02 K) | - | NIST Bridgeman and Aldrich, | DATA |
| H2O2 | L | 23.242, 4344.4, -52.550 | 330.0-446.0 | +0.00 (-0.00 K) | - | Landolt | DATA |
| CO | P | 45.698, -1076.6, -4.8814, 7.5673e-5, 2.0000 | 68.15-132.9 | -0.06 (+0.01 K) | -0.52 | Perry | DATA |
| CO2 | A | 6.8123, 1301.7, -3.4940 | 154.3-195.9 | -0.51 (+0.06 K) | - | NIST | DATA |
| CH4 | P | 39.205, -1324.4, -3.4366, 3.1019e-5, 2.0000 | 90.69-190.6 | +0.01 (-0.00 K) | -0.02 | Perry | DATA |
| C2H2 | LS | 25.522, 2645.6, -0.31000 | 148.0-190.0 | -0.26 (+0.03 K) | -2.13 | Landolt sub | DATA |
| C2H4 | A | 3.8726, 584.15, -18.307 | 149.4-188.6 | +0.05 (-0.01 K) | - | NIST Michels and Wassenaar, | DATA |
| C2H6 | a | 8.9541, 663.72, -16.469 | 133.8-198.2 | -0.01 (+0.00 K) | - | Poling | DATA |
| C3H8 | A | 4.0116, 834.26, -22.763 | 166.0-231.4 | +0.06 (-0.01 K) | - | NIST | DATA |
| C4H10 | L | 20.568, 2154.7, -34.361 | 196.0-298.0 | +0.01 (-0.00 K) | - | Landolt | DATA |
| C6H6 | L | 20.689, 2726.8, -55.623 | 280.0-374.0 | +0.04 (-0.01 K) | -0.89 | Landolt | DATA |
| HCN | L | 21.886, 2886.9, -20.307 | 261.0-300.0 | +0.17 (-0.05 K) | - | Landolt | DATA |
| CH2O | P | 101.51, -4917.2, -13.765, 0.022031, 1.0000 | 181.2-408.0 | -0.73 (+0.17 K) | - | Perry | DATA |
| CH3OH | a | 10.203, 1580.1, -33.650 | 262.6-356.0 | -0.20 (+0.05 K) | - | Poling | DATA |
| NH3 | P | 90.483, -4669.7, -11.607, 0.017194, 1.0000 | 195.4-405.6 | +0.43 (-0.09 K) | +1.60 | Perry | DATA |
| N2H4 | a | 9.9268, 1679.1, -45.450 | 298.9-408.4 | +0.00 (-0.00 K) | - | Poling | DATA |
| N2O | L | 21.020, 1506.5, -25.990 | 185.0-197.0 | +0.08 (-0.01 K) | +0.07 | Landolt | DATA |
| NO | a | 10.868, 682.94, -4.8800 | 106.9-127.6 | +0.01 (-0.00 K) | - | Poling | DATA |
| NO2 | A | 3.3525, 540.63, -131.93 | 217.6-294.0 | +7.90 (-1.61 K) | - | NIST Stull, 1947 | DATA |
| N2O4 | L | 21.897, 2730.2, -38.970 | 264.0-320.0 | -27.48 (+7.91 K) | - | Landolt | DATA |
| HNO3 | L | 22.190, 3237.4, -52.150 | 274.0-376.0 | +1.41 (-0.40 K) | - | Landolt | DATA |
| C2N2 | P | 81.565, -4808.9, -9.3748, 1.3901e-5, 2.0000 | 245.2-400.1 | +0.15 (-0.03 K) | -1.29 | Perry | DATA |
| CH3CN | L | 21.923, 3436.3, -24.208 | 280.0-354.0 | +0.26 (-0.08 K) | - | Landolt | DATA |
| HC3N | L | 22.254, 3384.8, 0 | 279.0-315.0 | +0.68 (-0.20 K) | - | Landolt | DATA |
| C4H2 | L | 19.203, 1628.5, -71.198 | 237.0-283.0 | +2.45 (-0.67 K) | - | Landolt | DATA |
| PH3 | a | 8.8405, 645.51, -17.080 | 137.4-199.5 | -0.01 (+0.00 K) | - | Poling | DATA |
| H2S | L | 21.250, 1858.0, -21.761 | 185.0-228.0 | +0.09 (-0.02 K) | -0.14 | Landolt | DATA |
| SO2 | a | 9.4072, 999.90, -35.960 | 199.7-279.5 | +0.02 (-0.00 K) | +0.77 | Poling | DATA |
| SO3 | a | 11.176, 1735.3, -36.650 | 284.5-332.0 | +0.01 (-0.00 K) | - | Poling | DATA |
| OCS | L | 20.773, 1845.7, -23.427 | 161.0-284.0 | -0.15 (+0.03 K) | - | Landolt | DATA |
| CS2 | L | 20.843, 2672.5, -32.498 | 255.0-354.0 | -0.49 (+0.15 K) | - | Landolt | DATA |
| CH3SH | A | 4.1920, 1031.4, -32.720 | 221.9-279.1 | -0.01 (+0.00 K) | - | NIST Russell, Osborne, et a | DATA |
| (CH3)2S | P | 84.390, -5740.6, -9.6454, 1.0073e-5, 2.0000 | 174.9-503.0 | +1.72 (-0.49 K) | - | Perry | DATA |
| HCl | L | 21.402, 1717.3, -14.268 | 160.0-201.0 | +0.14 (-0.03 K) | +0.84 | Landolt | DATA |
| HF | L | 22.579, 3397.7, 14.730 | 212.0-313.0 | +0.09 (-0.02 K) | - | Landolt | DATA |
| Cl2 | P | 71.334, -3855.0, -8.5171, 0.012378, 1.0000 | 172.1-417.1 | +0.17 (-0.04 K) | -0.61 | Perry | DATA |
| CH3Cl | A | 4.1545, 916.22, -28.466 | 183.0-249.4 | -0.57 (+0.13 K) | - | NIST Beersmans and Jungers, | DATA |
| HBr | a | 8.4124, 540.82, -47.710 | 185.1-221.5 | -0.06 (+0.01 K) | - | Poling | DATA |
| CF4 | a | 8.9589, 510.60, -15.950 | 106.2-155.5 | -0.03 (+0.00 K) | -0.26 | Poling | DATA |
| SF6 | a | 10.541, 1096.5, -11.150 | 162.0-220.4 | +0.35 (-0.05 K) | +3.32 | Poling | DATA |
| K | A | 4.4572, 4691.6, 24.195 | 679.4-1033 | +7.26 (-7.25 K) | - | NIST Fiock and Rodebush, 19 | DATA |
| NH4Cl | a | 11.481, 3703.7, -41.150 | 494.0-640.5 | -2.98 (+1.16 K) | - | Poling | DATA |

No reachable correlation: HDO, HNCO, CN, S2, S8, SO, SiO, TiO, HClO4, SO3 above 332 K. Na, Ca, Mg, Fe, NaCl, KCl, S2, TiO and H2SO4 use the NASA route (T6).

### 6.2 T5: Wagner fits to NIST EOS saturation tables (liquid branch)

P = Pc exp((a x + b x^1.5 + c x^2.5 + d x^5)/Tr), x = 1 - T/Tc, Tr = T/Tc, P in Pa. Tc is fixed from heos_constants; Pc, a to d are least-squares fits to the NIST saturation table between Tt and Tc; max err is the largest relative deviation from the table; dHvap(Tb) in kJ/mol from the same table. Valid from Tt to Tc only.

| Species | Range K | Tc | Pc Pa | a | b | c | d | max err % | dHvap Tb | Status |
|---|---|---|---|---|---|---|---|---|---|---|
| H2 | 13.957-33.15 | 33.145 | 1296303.5 | -4.89948 | 1.05654 | 0.74762 | 0.04196 | 0.015 | 0.90 | COMPUTED |
| He | 2.177-5.20 | 5.195 | 228365.6 | -4.11659 | 1.15844 | 1.01625 | 0.29612 | 0.064 | 0.08 | COMPUTED |
| Ne | 24.556-44.40 | 44.400 | 2661458.2 | -5.78536 | 1.14227 | -0.07076 | -1.30581 | 0.007 | 1.73 | COMPUTED |
| Ar | 83.806-150.69 | 150.687 | 4863032.1 | -5.92648 | 1.21423 | -0.52935 | -1.54635 | 0.004 | 6.44 | COMPUTED |
| Kr | 115.775-209.48 | 209.480 | 5525731.7 | -5.98453 | 1.34328 | -0.66158 | -1.31519 | 0.008 | 8.97 | COMPUTED |
| Xe | 161.405-289.73 | 289.733 | 5841283.1 | -5.99720 | 1.34294 | -0.72630 | -1.33548 | 0.012 | 12.55 | COMPUTED |
| N2 | 63.151-126.19 | 126.192 | 3395710.6 | -6.12290 | 1.25875 | -0.75882 | -1.79283 | 0.007 | 5.58 | COMPUTED |
| O2 | 54.361-154.58 | 154.581 | 5044902.0 | -6.07788 | 1.31005 | -0.72890 | -1.49664 | 0.110 | 6.82 | COMPUTED |
| H2O | 273.160-647.10 | 647.096 | 22071748.8 | -7.87987 | 1.92877 | -2.32932 | -2.05762 | 0.101 | 40.65 | COMPUTED |
| CO | 68.160-132.86 | 132.860 | 3497225.4 | -6.18485 | 1.29617 | -0.91665 | -2.03164 | 0.028 | 6.01 | COMPUTED |
| CO2 | 216.592-304.13 | 304.128 | 7376475.3 | -7.01631 | 1.48903 | -2.15456 | -3.18597 | 0.011 | - | COMPUTED |
| CH4 | 90.694-190.56 | 190.564 | 4599091.9 | -6.02354 | 1.26790 | -0.57083 | -1.36956 | 0.010 | 8.19 | COMPUTED |
| C2H4 | 103.986-282.35 | 282.350 | 5044941.7 | -6.42322 | 1.47150 | -1.24789 | -1.99668 | 0.138 | 13.53 | COMPUTED |
| C2H6 | 90.368-305.32 | 305.322 | 4881269.7 | -6.58230 | 1.69617 | -1.48179 | -1.56810 | 0.212 | 14.72 | COMPUTED |
| C3H8 | 85.525-369.89 | 369.890 | 4261023.3 | -6.88163 | 1.86363 | -1.92801 | -1.86160 | 0.231 | 18.77 | COMPUTED |
| C4H10 | 134.895-425.13 | 425.125 | 3798336.6 | -7.11781 | 1.88057 | -2.10643 | -2.20494 | 0.090 | 22.42 | COMPUTED |
| CH3OH | 175.610-513.38 | 513.380 | 8194162.4 | -8.67443 | 1.35531 | -2.74480 | -0.58626 | 0.262 | 35.28 | COMPUTED |
| NH3 | 195.490-405.56 | 405.560 | 11367055.2 | -7.28265 | 1.57703 | -1.88641 | -2.28486 | 0.072 | 23.33 | COMPUTED |
| N2O | 182.330-309.52 | 309.520 | 7240558.2 | -6.77061 | 1.52561 | -1.74908 | -2.49144 | 0.057 | 16.47 | COMPUTED |
| H2S | 187.700-373.10 | 373.100 | 8993628.0 | -6.48162 | 1.46237 | -1.29817 | -2.06320 | 0.057 | 18.62 | COMPUTED |
| SO2 | 197.700-430.64 | 430.640 | 7884667.6 | -7.25016 | 1.62324 | -2.21507 | -2.85747 | 0.024 | 24.96 | COMPUTED |
| OCS | 134.300-378.77 | 378.770 | 6359936.9 | -6.46670 | 1.44102 | -1.21326 | -2.18411 | 0.137 | 18.56 | COMPUTED |
| C6H6 | 278.674-562.02 | 562.020 | 4905466.6 | -7.11487 | 1.80003 | -2.11804 | -3.45642 | 0.037 | 30.75 | COMPUTED |
| CF4 | 120.000-226.93 | 227.510 | 3772442.6 | -6.99615 | 1.77979 | -1.87108 | -2.12372 | 0.051 | 11.83 | COMPUTED |
| SF6 | 223.555-318.72 | 318.723 | 3754644.4 | -7.07990 | 1.59103 | -1.93730 | -2.38244 | 0.009 | - | COMPUTED |

### 6.3 T6: NASA-polynomial equilibrium vapour pressure (metals, salts, TiO, sulfur, H2SO4)

P = exp(-dG/RT) x 1 bar for the condensed phase and gas from NASA-7 polynomials (cantera nasa_gas.yaml, nasa_condensed.yaml), COMPUTED from DATA. Columns: T (K) at 1 mPa, 1 Pa, 100 Pa, 10 kPa, 100 kPa; Tb (Yaws); fit ln(P/Pa) = a + b/T + c ln T with a, b, c, its range in K and its maximum fit error in percent. Check: the 100 kPa temperature agrees with Yaws Tb within 0.1 to 1.2% for Na, K, Ca, Mg and Fe and 3.0 to 3.5% for NaCl and KCl (100 versus 101.325 kPa explains up to 0.4%). The ideal-vapour assumption fails at high pressure: H2O is within 1% to 348 K, -4% at 423 K, -12% at 498 K and -25% at 572 K against the NIST table. S2 is the S2 partial pressure only; S8 was computed to 100 Pa only.

| Species | 1 mPa | 1 Pa | 100 Pa | 10 kPa | 100 kPa | Tb | Fit a, b, c; range K; max err % | Status |
|---|---|---|---|---|---|---|---|---|
| H2O | - | 213 | 253 | 319 | 373 | - | 94.6794, -8514.7, -10.1913; 200-394; 7.8 | COMPUTED |
| Na | 422 | 552 | 697 | 951 | 1168 | 1156 | 31.6630, -13041.8, -1.2735; 306-1253; 10.2 | COMPUTED |
| K | 359 | 472 | 601 | 834 | 1039 | 1037 | 32.0734, -10944.2, -1.4459; 261-1122; 10.5 | COMPUTED |
| Ca | 674 | 868 | 1077 | 1453 | 1778 | 1757 | 46.6015, -23075.2, -2.9532; 499-1899; 7.9 | COMPUTED |
| Mg | 548 | 699 | 859 | 1135 | 1364 | 1363 | 45.1718, -19114.1, -2.7187; 413-1442; 9.5 | COMPUTED |
| Fe | 1383 | 1729 | 2094 | 2683 | 3139 | 3133 | 60.7921, -53555.1, -4.0033; 1189-3311; 4.8 | COMPUTED |
| NaCl | 752 | 935 | 1129 | 1489 | 1791 | 1738 | 91.0374, -32883.4, -8.1684; 571-1902; 16.3 | COMPUTED |
| KCl | 729 | 907 | 1096 | 1447 | 1748 | 1689 | 91.7153, -31937.9, -8.2971; 560-1867; 16.0 | COMPUTED |
| TiO | 1774 | 2188 | 2629 | 3320 | 3841 | - | 88.0351, -78033.7, -6.8179; 1372-4034; 15.2 | COMPUTED |
| H2SO4 | - | 332 | 400 | 513 | 609 | - | 81.4810, -11911.9, -7.8600; 300-647; 1.9 | COMPUTED |
| S2 | 406 | 504 | 605 | 763 | 882 | - | 65.2805, -16872.9, -5.1118; 312-925; 8.8 | COMPUTED |
| S8 | 307 | 374 | 454 | - | - | - | - | COMPUTED |

## 7. Replacement for the two-point `CONDENSATION_K` rows

### 7.1 Recipe

Psat(T) is composite: Wagner liquid (T5) for Tt <= T <= Tc; below Tt the solid curve (T9). The pieces join within the "Err at Tt" column of T9 (O2, CO, Ar, H2O, H2 within 0.15%; CO2 +0.5%; CH4 +1.6%; H2S +1.7%; NH3 +1.8%; N2 -2.8%; SO2 +3.8%, accepted as the published spread). Condensation temperature at pressure p is the root of Psat(T) = p (bisection; Psat is monotonic); for a gas of mole fraction x use p = x P. H2 solid: Clausius-Clapeyron from the triple point (13.957 K, 7358 Pa) with dHsub = dHvap(Tt) + dHfus (CRC via chemicals), COMPUTED and not checked against a measurement.

### 7.2 T9: Fray and Schmitt 2009 sublimation polynomials

ln(P/bar) = A0 + A1/T + ... + A6/T^6, transcribed from SOM Table 2 of Lisse et al. 2021 (arXiv:2009.02277), 2ND-HAND. That paper's caption does not state the unit; bar was established here by reproducing the triple pressures (Err at Tt, COMPUTED). The original Fray and Schmitt paper was not reached, so per-polynomial validity ranges are NOT REACHED. Not used: CH3OH-1 (+23%), C6H6-1 (+7%), N2O-1 (+5.7%), N2O-2 (-100% at Tt).

| Name | A0 | A1 | A2 | A3 | A4 | A5 | A6 | Err at Tt % | Use | Status |
|---|---|---|---|---|---|---|---|---|---|---|
| O2-1 | 15.41 | -1148 | 334.9 | 60.21 | 0 | 0 | 0 | +154.22 | alpha, T<23.14 K | 2ND-HAND |
| O2-2 | 13.35 | -1012 | -2971 | 2.926e+04 | 0 | 0 | 0 | +54.95 | beta, 23.14-43.78 K | 2ND-HAND |
| O2-3 | 10.18 | -803.5 | -7080 | 7.553e+04 | 0 | 0 | 0 | +0.09 | gamma, 43.78 K to Tt | 2ND-HAND |
| O3-1 | 17.46 | -2352 | 0 | 0 | 0 | 0 | 0 | - | | 2ND-HAND |
| CO-1 | 10.43 | -721.3 | -1.074e+04 | 2.341e+05 | -2.392e+06 | 9.478e+06 | 0 | +3.32 | T<60.7 K | 2ND-HAND |
| CO-2 | 10.25 | -748.2 | -5843 | 3.939e+04 | 0 | 0 | 0 | +0.15 | 60.7 K to Tt | 2ND-HAND |
| CO2-1 | 14.76 | -2571 | -7.781e+04 | 4.325e+06 | -1.207e+08 | 1.35e+09 | 0 | -3.92 | alternative (Tt error -3.9%) | 2ND-HAND |
| CO2-2 | 18.61 | -4154 | 1.041e+05 | 0 | 0 | 0 | 0 | +0.54 | recommended only from 195.3 K up (Tt error +0.5%); in the cold it is 4 K too cold at 100 Pa and 17 K too cold at 1 mPa against this table's own sublimation heat (25.2 to 27.2 kJ/mol), so use CO2-1 below 195.3 K (GEN.452 lane, #2593) | 2ND-HAND |
| CH3OH-1 | 19.18 | -5648 | 0 | 0 | 0 | 0 | 0 | +23.47 | | 2ND-HAND |
| CH3OH-2 | 17.06 | -5314 | 0 | 0 | 0 | 0 | 0 | -0.72 | | 2ND-HAND |
| CH4-1 | 10.51 | -1110 | -4341 | 1.035e+05 | -7.91e+05 | 0 | 0 | +1.64 | all T<Tt | 2ND-HAND |
| C2H2-1 | 13.4 | -2536 | 0 | 0 | 0 | 0 | 0 | - | | 2ND-HAND |
| C2H4-1 | 15.4 | -2206 | -1.216e+04 | 2.843e+05 | -2.203e+06 | 0 | 0 | +0.43 | all T<Tt | 2ND-HAND |
| C2H6-1 | 15.11 | -2207 | -2.411e+04 | 7.744e+05 | -1.161e+07 | 6.763e+07 | 0 | +0.18 | all T<Tt | 2ND-HAND |
| C6H6-1 | 17.35 | -5663 | 0 | 0 | 0 | 0 | 0 | +7.09 | | 2ND-HAND |
| HCN-1 | 13.93 | -3624 | -1.325e+05 | 6.314e+06 | -1.128e+08 | 0 | 0 | - | all T<Tt | 2ND-HAND |
| HC3N-1 | 13.01 | -4426 | 0 | 0 | 0 | 0 | 0 | - | | 2ND-HAND |
| C2N2-1 | 16.53 | -4109 | 0 | 0 | 0 | 0 | 0 | - | | 2ND-HAND |
| C4N2-1 | 19.09 | -6036 | 0 | 0 | 0 | 0 | 0 | - | | 2ND-HAND |
| N2-1 | 12.4 | -807.4 | -3926 | 6.297e+04 | -4.633e+05 | 1.325e+06 | 0 | +153.55 | alpha T<36.08 K | 2ND-HAND |
| N2-2 | 8.514 | -458.4 | -1.987e+04 | 4.8e+05 | -4.524e+06 | 0 | 0 | -2.75 | beta, 36.08 K to Tt | 2ND-HAND |
| NH3-1 | 15.96 | -3537 | -3.31e+04 | 1.742e+06 | -2.995e+07 | 0 | 0 | +1.83 | all T<Tt | 2ND-HAND |
| N2O-1 | 16.22 | -2971 | 0 | 0 | 0 | 0 | 0 | +5.66 | | 2ND-HAND |
| N2O-2 | 6.566 | -1271 | -6.684e+05 | 4.496e+07 | -1.097e+09 | 0 | 0 | -100.00 | | 2ND-HAND |
| H2S-1 | 12.98 | -2707 | 0 | 0 | 0 | 0 | 0 | +1.67 | all T<Tt | 2ND-HAND |
| H2S-2 | 8.933 | -726 | -3.504e+05 | 2.724e+07 | -8.582e+08 | 0 | 0 | +0.57 | low T alternative | 2ND-HAND |
| SO2-1 | 15.6 | -3508 | -9.401e+04 | 4.152e+06 | -6.946e+07 | 0 | 0 | +3.79 | all T<Tt | 2ND-HAND |
| Ne-1 | 9.886 | -269.9 | 128.3 | -162.4 | 0 | 0 | 0 | -6.54 | | 2ND-HAND |
| Ne-2 | 10.61 | -308.6 | 986 | -9069 | 3.514e+04 | 0 | 0 | -0.21 | | 2ND-HAND |
| Ar-1 | 10.69 | -893.2 | -3567 | 6.574e+04 | -4.28e+05 | 0 | 0 | -0.02 | all T<Tt | 2ND-HAND |
| Kr-1 | 10.77 | -1223 | -8903 | 2.635e+05 | -4.26e+06 | 3.575e+07 | -1.21e+08 | -0.26 | | 2ND-HAND |
| Xe-1 | 10.98 | -1737 | -1.332e+04 | 4.349e+05 | -7.027e+06 | 4.447e+07 | 0 | +0.17 | | 2ND-HAND |

### 7.3 T10: laboratory curves versus F and S at low temperature

ln(p/bar) = A - B/T + C ln T + D2 T + D3 T^2 + D4 T^3 from Grundy et al. 2024 (arXiv:2309.05078 Table 2: alpha-N2, beta-N2, alpha-CO, beta-CO, CH4 I) and Blakley et al. 2024 (PSS 244, 105863, arXiv:2403.09917; NASA NTRS record dated 2025) (NTRS 20250005811 Table 1: NH3, beta-O2), READ. Script check against the Grundy reference pressures at T0: alpha-N2 at 27 K 0.79e-9 bar (listed 0.76), alpha-CO at 31 K 1.66e-9 (1.64), CH4 I at 37 K 0.48e-9 (0.43), so the transcribed coefficients reproduce Grundy's Table 1 reference pressures within 11%. F and S exceeds the lab curves by up to 1.4 to 2.3 (alpha-N2, 25 to 33 K), 1.9 (CH4, 30 K), 4 (beta-O2, 24 K) and 6 (alpha-CO, 31 K); they agree near 60 K (CO) and 70 K (CH4), and within 12% at 40 K (O2). Prefer the lab curves for bodies colder than 40 K.

| Curve | Ratio lab / F and S | Status |
|---|---|---|
| alpha-N2 lab/F&S (Grundy 2024) | 25 K: 0.72; 27 K: 0.61; 30 K: 0.51; 33 K: 0.43 | READ |
| alpha-CO lab/F&S (Grundy 2024) | 31 K: 0.17; 40 K: 0.38; 50 K: 0.70; 60 K: 1.02 | READ |
| CH4 I lab/F&S (Grundy 2024) | 30 K: 0.52; 37 K: 0.65; 50 K: 0.83; 70 K: 1.01 | READ |
| beta-O2 lab/F&S (Blakley 2024) | 24 K: 0.24; 30 K: 0.46; 40 K: 0.88 | READ |
| NH3 lab/F&S (Blakley 2024) | 100 K: 0.68; 150 K: 0.92; 195.49 K: 1.05 | READ |

### 7.4 New entries: NH3, H2S, NH4SH, HCN, C2H6, H2SO4

NH3, H2S, C2H6: composite recipe (T5, T9). HCN: Perry DIPPR101 liquid ln P = 36.75 - 3927.1/T - 2.1245 ln T + 3.8948e-17 T^6 (259.83 to 456.65 K; DATA chemicals) above Tt = 259.86 K, F and S HCN-1 below (3.3% step at Tt); Tb from the fit 298.81 K against NIST 298.79 K and Landolt 298.95 K (COMPUTED). H2SO4: T6 NASA fit.

NH4SH condenses from NH3 + H2S; the condition needs both partial pressures. READ in Lodders and Fegley 2024 (2410.11138): log10 Kvap = 14.89 - 4732/T, Kvap = p(NH3) p(H2S) over the solid; 14.82 - 4705/T (condensation temperatures within 0.5 K). Cloud forms when log Kvap(T) <= log(x_NH3 x_H2S) + 2 log P_T (bar). Closed forms, valid for solar N/S only (about 5): Visscher et al. 2006 (astro-ph/0511136 eq. 32) 10^4/Tcond = 48.91 - 4.15 (log P_T + 0.5 [N/H] + 0.5 [S/H]); Lodders and Fegley 2024 10^4/Tcond = 49.18 - 4.225 (same bracket). COMPUTED with Kelley and the combined log(x_NH3 x_H2S) = -8.378 of that paper: Tcond = 151.3, 161.7, 173.5, 187.3, 203.4, 222.5, 245.6 K at log P_T = -4 to 2 (printed table 151, 162, 174, 187, 203, 222, 246). The 2006 fit is up to 1.3 K higher. At [M/H] = +0.5 and 1 bar the 2006 fit gives 213.5 K; each +1 in [M/H] adds about 17 K at 1 bar. Pure-solid dissociation pressure 2 sqrt(Kvap): 0.121 bar at 273.15 K, 0.646 bar at 298.15 K (COMPUTED; Kelley's pressure unit not stated, effect below 1%).

H2SO4 (2ND-HAND, Gao et al. 2014, arXiv:1312.3750 eqs. 4a to 4c, read from the page image): ln p(atm) = ln p0 + 10156 [-1/T + 1/T0 + 0.38/(Tc - T0) (1 + ln(T0/T) - T0/T)] - H/(R T), T0 = 340 K, Tc = 905 K, ln p0 = -10156/T0 + 16.259 (Ayers 1980 as modified by Kulmala and Laaksonen 1990); H = 4.184e7 [23624.8 - 1.14208e8/(4798.69 + (W - 105.315)^2)] erg/mol, W the acid weight percent from Tabazadeh 1997 (not reached). With H = 0 (pure acid) at 300, 340, 400, 500, 600 K this gives 0.0022, 0.124, 10.1, 1124, 20716 Pa (COMPUTED); Ayers alone gives 0.0023, 0.124, 11, 1760, 52000 Pa; the NASA H2SO4(L) route gives 0.048, 1.87, 99.6, 6750, 84800 Pa (20 times higher at 300 K). Solution factor exp(-H/RT) at 340 K for W = 98, 90, 80, 70, 60 %: 0.88, 0.25, 0.020, 9.1e-4, 3.4e-5 (COMPUTED from the READ formula). A droplet of 75 to 85 percent acid holds the acid at about 2e-3 to 8e-2 of the pure-acid pressure at 300 to 340 K; W(T, pH2O) is a Boss item.

| Entry | T 100 kPa | 10 kPa | 1 kPa | 100 Pa | 10 Pa | 1 Pa | 1 mPa | Max dev two-point K | Recipe | Status |
|---|---|---|---|---|---|---|---|---|---|---|
| H2 | 20.3 | 14.5 | 11.4 | 9.4 | 8.0 | 7.0 | 5.0 | 0.2 | Wagner liquid above Tt 13.957 K; solid by Clausius-Clapeyron from the triple point | COMPUTED |
| N2 | 77.2 | 62.2 | 52.8 | 46.0 | 40.8 | 36.7 | 28.7 | 1.1 | Wagner liquid; F&S N2-2 (beta) 36.08 K to Tt; N2-1 (alpha) below | COMPUTED |
| CO | 81.5 | 65.9 | 56.6 | 49.6 | 44.3 | 40.0 | 31.0 | 1.6 | Wagner liquid; F&S CO-2 above 60.7 K, CO-1 below | COMPUTED |
| Ar | 87.2 | 71.7 | 61.2 | 53.4 | 47.4 | 42.6 | 32.7 | 0.4 | Wagner liquid; F&S Ar-1 | COMPUTED |
| O2 | 90.1 | 72.7 | 61.3 | 53.2 | 47.4 | 42.9 | 33.8 | 1.2 | Wagner liquid; F&S O2-3 (43.78 K to Tt), O2-2 (23.14-43.78 K), O2-1 below | COMPUTED |
| CH4 | 111.5 | 89.5 | 76.1 | 66.3 | 58.7 | 52.7 | 40.3 | 1.7 | Wagner liquid; F&S CH4-1 | COMPUTED |
| CO2 | 194.4 | 169.2 | 148.8 | 131.8 | 117.5 | 105.0 | 74.4 | 7.9 | Wagner liquid above 216.592 K; F&S CO2-2 solid below (CO2-2 is wrong in the cold, use CO2-1 below 195.3 K: see the CO2-2 row) | COMPUTED |
| SO2 | 262.8 | 221.0 | 192.7 | 174.3 | 159.3 | 146.7 | 118.9 | 6.9 | Wagner liquid; F&S SO2-1 | COMPUTED |
| H2O | 372.8 | 318.9 | 280.1 | 252.8 | 231.0 | 212.6 | 171.5 | 6.3 | Wagner liquid above 273.16 K; Feistel-Wagner ice formula below | COMPUTED |
| H2S | 212.6 | 177.1 | 153.9 | 136.1 | 122.0 | 110.5 | 86.2 | 2.3 | Wagner liquid; F&S H2S-1 | COMPUTED |
| NH3 | 239.6 | 202.0 | 178.6 | 161.0 | 146.6 | 134.6 | 108.1 | 5.3 | Wagner liquid; F&S NH3-1 | COMPUTED |
| C2H6 | 184.3 | 150.4 | 127.9 | 111.8 | 99.5 | 89.9 | 72.2 | 3.4 | Wagner liquid; F&S C2H6-1 | COMPUTED |
| HCN | 298.5 | 250.1 | 221.4 | 198.9 | 180.8 | 165.8 | 133.2 | 7.2 | Perry DIPPR101 liquid above 259.86 K; F&S HCN-1 solid below (+3.3% step at Tt) | COMPUTED |
| H2SO4 | 608.4 | 512.9 | 448.0 | 399.9 | 362.5 | 332.3 | - | 8.8 | NASA H2SO4(L) to gas fit, 300-647 K only | COMPUTED |

Columns give the checked condensation temperature (K) at each pressure; the H2SO4 row is the NASA fit (the pure-acid K and L curve is 42 K warmer at 100 Pa, 442 K versus 400 K, and 71 K warmer at 100 kPa, 679 K versus 608 K; it is extrapolated far above the temperatures of the Venus aerosol layer). "Max dev two-point" is the largest temperature error between 1 mPa and 100 kPa of a two-point Clausius-Clapeyron interpolation anchored at this table's 100 kPa and 1 Pa temperatures (the form `vapour_pressure_kpa` uses): 0.2 to 2.3 K for H2, Ar, N2, O2, CO, CH4, H2S; 3.4 K for C2H6; 5.3 to 8.8 K for NH3, H2O, SO2, HCN, CO2 and H2SO4. Corrected two-point anchors are enough for the first group and not for the second, which should use the composite curve.

### 7.5 T8: old two-point rows against the checked curves

Each cell: checked / project temperature (K) at that pressure; last column checked minus project at 1 Pa (the project's 0.001 kPa anchor). H2S in the project is a copy of the SO2 row: 20 to 50 K warm. CO2 project row is 0.6 to 7 K warm, worst at low pressure. N2, CO, Ar, O2, CH4 are 0.7 to 5 K cold, worst at low pressure. SO2 is 2 to 16 K cold. H2O is 3 to 6 K warm at 100 Pa to 10 kPa (its pressure low by 18 to 31% mid-range). NH3, C2H6 have no project row.

| Species | 100 kPa | 10 kPa | 1 kPa | 100 Pa | 1 Pa | 1 mPa | dT at 1 Pa K | Status |
|---|---|---|---|---|---|---|---|---|
| H2 | 20.32 / 20 | 14.51 / 13.64 | 11.40 / 10.34 | 9.42 / 8.33 | 6.98 / 6.00 | 5.0 / 4.2 | +1.0 | COMPUTED |
| N2 | 77.24 / 77 | 62.18 / 60.79 | 52.84 / 50.22 | 46.03 / 42.78 | 36.70 / 33.00 | 28.7 / 24.6 | +3.7 | COMPUTED |
| CO | 81.52 / 81 | 65.94 / 64.80 | 56.55 / 54.00 | 49.64 / 46.29 | 39.98 / 36.00 | 31.0 / 27.0 | +4.0 | COMPUTED |
| Ar | 87.18 / 87 | 71.68 / 69.16 | 61.18 / 57.40 | 53.40 / 49.05 | 42.62 / 38.00 | 32.7 / 28.4 | +4.6 | COMPUTED |
| O2 | 90.06 / 90 | 72.68 / 72.00 | 61.28 / 60.00 | 53.24 / 51.43 | 42.85 / 40.00 | 33.8 / 30.0 | +2.9 | COMPUTED |
| CH4 | 111.51 / 112 | 89.50 / 88.42 | 76.15 / 73.04 | 66.29 / 62.22 | 52.70 / 48.00 | 40.3 / 35.7 | +4.7 | COMPUTED |
| CO2 | 194.45 / 195 | 169.22 / 169.83 | 148.80 / 150.41 | 131.85 / 134.98 | 104.98 / 112.00 | 74.4 / 89.2 | -7.0 | COMPUTED |
| SO2 | 262.84 / 263 | 220.98 / 218.89 | 192.71 / 187.45 | 174.33 / 163.91 | 146.71 / 131.00 | 118.9 / 100.7 | +15.7 | COMPUTED |
| H2O | 372.77 / 373 | 318.95 / 323.82 | 280.12 / 286.09 | 252.81 / 256.24 | 212.58 / 212.00 | 171.5 / 168.4 | +0.6 | COMPUTED |
| H2S | 212.60 / 263 | 177.13 / 218.89 | 153.94 / 187.45 | 136.11 / 163.91 | 110.52 / 131.00 | 86.2 / 100.7 | -20.5 | COMPUTED |
| NH3 | 239.57 / none | 201.99 / - | 178.60 / - | 161.04 / - | 134.64 / - | 108.1 / - | - | COMPUTED |
| C2H6 | 184.32 / none | 150.44 / - | 127.92 / - | 111.76 / - | 89.93 / - | 72.2 / - | - | COMPUTED |

## 8. Henry constants

kH(T) = kH(298.15 K) exp(d (1/T - 1/298.15)), kH in mol/(kg bar), d = -dln kH/d(1/T) in K. Columns: NIST WebBook (compiling Sander 2015) median with [min-max, n], its median d, and the independent thermo 0.6.1 fit (Sander_henry_const.json, Sander_henry_T_dep.json; converted from mole-fraction constants assuming dilute solution). For dilute physical solubility between 273 and 323 K. CO2, SO2, NH3, HCl, H2S and HCN dissolve more with pH or dissociation (not modelled); the HCl and H2O2 spread shows this.

| Species | kH NIST median [range, n] | d K | kH thermo | d thermo K | Status |
|---|---|---|---|---|---|
| CO2 | 0.034 [0.031-0.045, n=15] | 2400 | 0.0375 | 2357 | DATA |
| SO2 | 1.2 [1.2-1.5, n=14] | 3100 | 2.26 | 2853 | DATA |
| NH3 | 58 [10-78, n=15] | 4100 | 52.5 | 3256 | DATA |
| H2S | 0.1 [0.001-0.1, n=9] | 2150 | 0.0952 | 1976 | DATA |
| O2 | 0.0013 [0.0012-0.0013, n=6] | 1700 | 0.00128 | 1602 | DATA |
| N2 | 0.000625 [0.0006-0.00065, n=2] | 1300 | 0.000628 | 1427 | DATA |
| CH4 | 0.0014 [0.00097-0.0092, n=10] | 1750 | 0.00159 | 1731 | DATA |
| HCl | 19.5 [1.1-2.5e+03, n=6] | 2000 | 6.05 | 1521 | DATA |
| H2 | 0.00078 [0.00078-0.00078, n=4] | 500 | 0.000782 | 560 | DATA |
| He | 0.000375 [0.00037-0.00038, n=2] | 226 | 0.00039 | 135 | DATA |
| Ne | 0.00045 [0.00045-0.00045, n=2] | 490 | 0.000457 | 488 | DATA |
| Ar | 0.0014 [0.0014-0.0014, n=2] | 1300 | 0.0014 | 1621 | DATA |
| Kr | 0.00245 [0.0024-0.0025, n=2] | 1700 | 0.00199 | 1840 | DATA |
| Xe | 0.0043 [0.0043-0.0043, n=2] | 2050 | 0.0045 | 2386 | DATA |
| CO | 0.00095 [0.00086-0.0074, n=5] | 1300 | 0.00101 | 1503 | DATA |
| NO | 0.0019 [0.0014-0.0019, n=5] | 1500 | 0.00166 | 1446 | DATA |
| N2O | 0.025 [0.024-0.026, n=7] | 2650 | 0.0219 | 2683 | DATA |
| O3 | 0.0115 [0.0088-0.013, n=10] | 2400 | 0.00155 | 2413 | DATA |
| C2H6 | 0.002 [0.0018-0.011, n=6] | 2350 | 0.00171 | 2438 | DATA |
| C2H2 | 0.041 [0.039-0.042, n=3] | 1800 | 0.0584 | 1713 | DATA |
| HCN | 9.3 [7.5-12, n=3] | 5000 | 6.56 | 4623 | DATA |
| H2O2 | 8.95e+04 [6.9e+04-1.4e+05, n=8] | 7300 | 9.02e+04 | 6973 | DATA |
| HNO3 | 2.1e+05 [8.9e+04-2.1e+05, n=3] | 8700 | 2.22e+05 | 8750 | DATA |
| Cl2 | 0.088 [0.061-0.095, n=6] | 2550 | 0.0937 | 2360 | DATA |
| CH3Cl | 0.11 [0.029-0.13, n=11] | 2900 | 0.0806 | 2560 | DATA |
| (CH3)2S | 0.55 [0.08-0.7, n=11] | 3500 | 0.587 | 3611 | DATA |
| OCS | 0.021 [0.015-0.033, n=9] | 3300 | 0.00324 | 2787 | DATA |
| CS2 | 0.055 [0.044-0.075, n=5] | 3400 | 0.0352 | 3351 | DATA |
| NO2 | 0.024 [0.007-0.041, n=7] | 2150 | 0.0123 | 2377 | DATA |
| PH3 | 0.0081 [0.0081-0.0081, n=1] | 2000 | 0.00616 | 2871 | DATA |
| Rn | 0.0093 [0.0093-0.0093, n=1] | 2600 | 0.00893 | 2983 | DATA |
| CH3OH | 220 [1.4e+02-230, n=7] | 5400 | 47.3 | 5275 | DATA |
| C6H6 | 0.18 [0.12-0.22, n=23] | 4050 | 0.18 | 3303 | DATA |
| CF4 | 0.0002 [0.00019-0.00021, n=4] | 1650 | 0.000221 | 971 | DATA |
| SF6 | 0.00024 [0.00024-0.00024, n=1] | 2400 | 0.00021 | 2745 | DATA |

## 9. Viscosity

Chapman-Enskog: mu = 2.6693 sqrt(M T)/(sigma^2 Omega22) in uPa s (M g/mol, T K, sigma Angstrom; the part first printed 26.693, which is micropoise: GEN.410 corrected it to 2.6693 and the table below reproduces), Omega22 = 1.16145 T*^-0.14874 + 0.52487 exp(-0.7732 T*) + 2.16178 exp(-2.43787 T*), T* = kT/eps (Neufeld constants RECALLED, validated against NIST fluid values at 300 K). sigma and eps/k from chemicals (Poling table). H2O is 10% off (polar): use NIST. Dilute gas only. Diffusion coefficients are not tabulated (Boss item 14).

| Species | sigma A | eps/k K | mu 300 K | NIST 300 K | diff % | Status |
|---|---|---|---|---|---|---|
| H2 | 2.827 | 59.7 | 8.89 | 8.94 | -0.6 | COMPUTED |
| He | 2.551 | 10.22 | 20.23 | 19.93 | +1.5 | COMPUTED |
| Ne | 2.82 | 32.8 | 31.24 | 31.84 | -1.9 | COMPUTED |
| Ar | 3.542 | 93.3 | 22.82 | 22.72 | +0.4 | COMPUTED |
| Kr | 3.655 | 178.9 | 25.24 | 25.45 | -0.8 | COMPUTED |
| Xe | 4.047 | 231.0 | 23.09 | 23.15 | -0.2 | COMPUTED |
| N2 | 3.798 | 71.4 | 17.70 | 17.88 | -1.0 | COMPUTED |
| O2 | 3.467 | 106.7 | 20.57 | 20.63 | -0.3 | COMPUTED |
| H2O | 2.641 | 809.1 | 10.76 | 9.77 | +10.1 | COMPUTED |
| H2O2 | 4.196 | 289.3 | 9.80 | - | - | COMPUTED |
| CO | 3.69 | 91.7 | 17.68 | 17.84 | -0.9 | COMPUTED |
| CO2 | 3.941 | 195.2 | 15.18 | 14.99 | +1.3 | COMPUTED |
| CH4 | 3.758 | 148.6 | 11.19 | 11.12 | +0.6 | COMPUTED |
| C2H2 | 4.033 | 231.8 | 10.34 | - | - | COMPUTED |
| C2H4 | 4.163 | 224.7 | 10.22 | 10.36 | -1.4 | COMPUTED |
| C2H6 | 4.443 | 215.7 | 9.46 | 9.28 | +2.0 | COMPUTED |
| C3H8 | 5.118 | 237.1 | 8.27 | 8.18 | +1.1 | COMPUTED |
| C4H10 | 4.687 | 531.4 | 7.46 | 7.45 | +0.1 | COMPUTED |
| HCN | 3.63 | 569.1 | 8.20 | - | - | COMPUTED |
| CH3OH | 3.626 | 481.8 | 9.74 | 9.72 | +0.2 | COMPUTED |
| NH3 | 2.9 | 558.3 | 10.30 | 10.18 | +1.1 | COMPUTED |
| N2O | 3.828 | 232.4 | 14.90 | 14.90 | +0.0 | COMPUTED |
| NO | 3.492 | 116.7 | 19.14 | - | - | COMPUTED |
| PH3 | 3.981 | 251.5 | 11.67 | - | - | COMPUTED |
| H2S | 3.623 | 301.1 | 12.89 | 12.45 | +3.5 | COMPUTED |
| SO2 | 4.112 | 335.4 | 12.97 | 12.97 | -0.0 | COMPUTED |
| OCS | 4.13 | 336.0 | 12.43 | 12.45 | -0.2 | COMPUTED |
| CS2 | 4.483 | 467.0 | 9.98 | - | - | COMPUTED |
| HCl | 3.339 | 344.7 | 14.62 | - | - | COMPUTED |
| HF | 3.148 | 330.0 | 12.47 | - | - | COMPUTED |
| Cl2 | 4.217 | 316.0 | 13.38 | - | - | COMPUTED |
| CH3Cl | 4.182 | 350.0 | 10.88 | - | - | COMPUTED |
| SiH4 | 4.084 | 207.6 | 11.77 | - | - | COMPUTED |
| C2N2 | 4.361 | 348.6 | 10.18 | - | - | COMPUTED |
| C6H6 | 5.349 | 412.3 | 7.58 | 7.62 | -0.6 | COMPUTED |
| CF4 | 4.662 | 134.0 | 17.62 | 17.29 | +1.9 | COMPUTED |
| SF6 | 5.128 | 222.1 | 15.45 | 15.29 | +1.0 | COMPUTED |
| HBr | 3.353 | 449.0 | 18.77 | - | - | COMPUTED |

## 10. Ideal versus real gas

Z = P/(rho_molar R T). van der Waals: a = 27 R^2 Tc^2/(64 Pc), b = R Tc/(8 Pc). Peng-Robinson: a = 0.45724 R^2 Tc^2 alpha/Pc, b = 0.07780 R Tc/Pc, alpha = (1 + kappa (1 - sqrt(Tr)))^2, kappa = 0.37464 + 1.54226 omega - 0.26992 omega^2 (RECALLED constants, validated by the table). Largest real root. omega from the NIST EOS saturation table (T2). Z NIST = fluid EOS (DATA). At 90 bar CO2 and 1000 bar H2, Peng-Robinson is within 0.3 to 1.1% (CO2 at 735 and 400 K), 7% (CO2 at 300 K, liquid-like) and 0.3 to 3.9% (H2); van der Waals is off by 2.6 to 48%. CO2 at 735 K and 90 bar has Z = 1.0053, so the ideal density 64.81 kg/m3 is 0.5% high. H2 at 300 K and 1000 bar has Z = 1.64: the ideal gas overstates molar density by 64%. N2 there has Z = 1.97 and PR is 10% low. Mixtures need mixing rules and interaction parameters (Boss item 13).

| Species | T K | P bar | Z NIST | Z vdW (err) | Z PR (err) | omega | Status |
|---|---|---|---|---|---|---|---|
| CO2 | 735 | 90.0 | 1.0053 | 0.9789 (-2.6%) | 1.0081 (+0.3%) | 0.225 | COMPUTED |
| CO2 | 735 | 92.1 | 1.0055 | 0.9785 (-2.7%) | 1.0083 (+0.3%) | 0.225 | COMPUTED |
| CO2 | 300 | 90.0 | 0.2033 | 0.3000 (+47.6%) | 0.2176 (+7.0%) | 0.225 | COMPUTED |
| CO2 | 400 | 90.0 | 0.8366 | 0.7969 (-4.7%) | 0.8273 (-1.1%) | 0.225 | COMPUTED |
| H2 | 300 | 1000.0 | 1.6415 | 1.9727 (+20.2%) | 1.5771 (-3.9%) | -0.219 | COMPUTED |
| H2 | 1000 | 1000.0 | 1.1865 | 1.2988 (+9.5%) | 1.1824 (-0.3%) | -0.219 | COMPUTED |
| N2 | 300 | 1000.0 | 1.9673 | 2.2446 (+14.1%) | 1.7633 (-10.4%) | 0.037 | COMPUTED |
| H2O | 700 | 220.0 | 0.6793 | 0.6720 (-1.1%) | 0.6694 (-1.5%) | 0.344 | COMPUTED |
| CH4 | 300 | 90.0 | 0.8663 | 0.8288 (-4.3%) | 0.8445 (-2.5%) | 0.011 | COMPUTED |

## 11. Formulas for the generator

1. M_bar = sum x_i M_i (g/mol), x_i = p_i / sum p_j.
2. Molar density n = P/(Z R T) (mol/m3); mass density rho = n M_bar/1000 (kg/m3, M_bar in g/mol). Z = 1 when P is below about 1 bar and T above 1.5 Tc of the major species; otherwise PR Z (section 10) with Kay's rule for mixture critical constants (RECALLED, untested).
3. n_i = x_i n.
4. Scale height H = R T/(M_bar g) (M_bar in kg/mol); P(z) = P0 exp(-integral dz/H); for tall columns g(z) = g0 (Rb/(Rb + z))^2 and local T(z). Well-mixed below the homopause; above it each species has its own H_i.
5. Cp = sum x_i Cp_i, gamma = Cp/(Cp - R), dry adiabat Gamma = g/cp with cp = Cp/M_bar (J/kg/K). Dry air at 300 K from T3: Cp = 29.10, gamma = 1.400, Gamma = 9.76 K/km (COMPUTED).
6. Saturation: composite (section 7); condensation level is the lowest height where x_i P(z) >= Psat_i(T(z)). Local slope dlnP/dT = dH/(R T^2) (accurate to a few percent over 30 K).
7. Henry: dissolved mol per kg water = kH(T) x_i P(bar).
8. Photolysis threshold lambda_nm = 119626.565/dH(kJ/mol).

## 12. Test-case table

All rows are TEST POINTS (known outcome). Planet inputs are calibration, not presets.

| Test point | Input | Reference | Computed here | Status |
|---|---|---|---|---|
| TEST POINT: Earth dry air M | N2 0.7812, O2 0.2096, Ar 0.0092; or nine-species N2 0.78084, O2 0.209476, Ar 0.00934, CO2 314e-6, Ne, He, CH4, Kr, H2 (RECALLED fractions) | 28.96 g/mol (task) | 28.9589 (three gases), 28.9648 (nine species) | COMPUTED |
| TEST POINT: Earth sea-level density | 101325 Pa, 288.15 K | 1.225 kg/m3 (RECALLED) | 1.2247 kg/m3 | COMPUTED |
| TEST POINT: H2O at 273.16 K | T = 273.16 K | 611.657 Pa (task) | NIST EOS table 611.6544 Pa (-0.0004%); Feistel-Wagner ice formula 611.6577 (equal by construction); Wagner fit T5 611.04 (-0.10%); NASA route 612.24 (+0.10%) | DATA / COMPUTED |
| TEST POINT: ice below the triple point | 273.15, 250, 200 K | none | 611.15, 76.03, 0.1626 Pa | COMPUTED |
| TEST POINT: CO2 triple point | 216.58 K | 517.95 kPa (task) | EOS 518.0 kPa at 216.592 K; Wagner fit 517.99 kPa (216.592 K), 517.72 kPa (216.58 K, -0.05%); F and S CO2-2 520.76 kPa (+0.54%); CO2-1 497.7 kPa (-3.9%); NIST Antoine extrapolated 505.3 kPa (-2.4%) | DATA / COMPUTED |
| TEST POINT: CO2 sublimation at 101325 Pa | solve Psat = 101325 | 194.686 K (RECALLED) | 194.607 (CO2-2), 194.747 (CO2-1), 194.733 K (NIST Antoine) | COMPUTED |
| TEST POINT: boiling points | N2, O2, Ar, CH4, CO, H2O | heos Tb 77.355, 90.188, 87.302, 111.67, 81.638, 373.12 K | composite at 100 kPa: 77.24, 90.06, 87.18, 111.51, 81.52, 372.77 K (0.12 to 0.35 K lower, as 1.3% pressure implies) | COMPUTED |
| TEST POINT: scale heights | listed inputs | Earth 8.4, Mars 11, Venus 16, Titan 21 km (RECALLED) | Earth 288.15 K, 9.80665: 8.44 km; Mars 210 K, 3.71, CO2 0.9532 N2 0.027 Ar 0.016 O2 0.0013 CO 0.0008 (M 43.48): 10.82; Venus 737 K, 8.87, CO2 0.965 N2 0.035 (M 43.45): 15.90; Titan 94 K, 1.352, N2 0.95 CH4 0.05 (M 27.42): 21.09 km | COMPUTED, inputs RECALLED |
| TEST POINT: ideal vs real CO2 | 735 K, 90 bar | none | Z 1.0053; ideal 64.81, real 64.47 kg/m3 | COMPUTED / DATA |
| TEST POINT: conversion | 79.1 kPa N2, 21.2 kPa O2, 0.9 kPa Ar, 288.15 K | none | P 101.2 kPa; x 0.7816, 0.2095, 0.0089; n 42.2 mol/m3 | COMPUTED |

## 13. Boss must research

1. Fray and Schmitt 2009 (Planet. Space Sci. 57, 2053): validity per polynomial, enthalpies, data. ScienceDirect refused; T9 rests on it.
2. Solid H2 vapour pressure (Woolley); HD, HDO, D2O, 3He, 4He vapour pressures.
3. Jager and Span 2012 and Span and Wagner 1996 (CO2 sublimation and EOS; 0.5% spread at Tt).
4. H2SO4-H2O: Tabazadeh 1997 W(T, pH2O), Zeleznik 1991, Marti 1997, and the Ayers 1980 original. Venus cloud acid pressure depends on W by orders of magnitude.
5. NH4SH constant, units and range; Lodders 2003 abundances.
6. Clathrate hydrate dissociation curves (NIST SRD 156 subscription only; Sloan and Koh; Choukroun et al. 2010).
7. Perchlorates and carbonates: HClO4 has no NIST data or Cp; the salts need decomposition and vapour data.
8. Ortho-para equilibrium H2 Cp, 20 to 300 K (JANAF H2); Lobo-Ferreira data ranges (openknowledge.nau.edu/id/eprint/6246).
9. Gibbs formation data for FeH, CrH, CaOH, VO, AlH (FastChem fits).
10. Ideal-gas Cp for HC3N, HClO4, NH4Cl.
11. Rn triple point and a better vapour curve (Poling Antoine is +9.5% at Tb).
12. Vapour pressure for SiO, SiH4, SO, HNCO, C3O2, CN, S2 (total sulfur vapour), SO3 above 332 K; real critical constants for metals.
13. Binary interaction parameters and tested mixing rules (H2-He, H2-CH4, CO2-N2, CO2-H2O).
14. Photolysis channels, cross sections and yields for NUV-flagged species; diffusion coefficients (parts 03 and neighbours).

## 14. Recalled / to verify

1. Neufeld Omega22 constants (matched NIST viscosity within 2% for non-polar species).
2. Peng-Robinson constants 0.45724, 0.07780, 0.37464, 1.54226, 0.26992 (matched NIST Z within a few percent).
3. h c N_A = 119626.565 kJ nm/mol.
4. Mars, Venus, Titan T, g, compositions; Earth fractions; 28.9644 g/mol; 1.225 kg/m3; scale-height references.
5. The task figures 611.657 Pa, 517.95 kPa and 194.686 K (task statement; not read in a source here).
6. A search snippet showed Lewis as 14.83/4715; the READ value is 14.82/4705.
7. Flag rules (IR, volatility, NUV) are definitions of this note, not sourced.
8. Kay's rule for mixtures.

## 15. Row counts by status

Table rows by label (a row carries its dominant label): DATA 206; DATA with extrapolated Cp (x) 7; COMPUTED 195; 2ND-HAND 32 (T9); READ 5 (T10); NOT REACHED 13; RECALLED 0 and SNIPPET 0 (RECALLED items are listed only in section 14, one SNIPPET noted there). Formula-level READ or 2ND-HAND items outside tables: NH4SH constants (READ), H2SO4 Ayers, Kulmala-Laaksonen and Giauque terms (2ND-HAND).

## 16. Cross-links

Feeds: 02 and 05 composition (05-mixing-ratios-by-class.md uses T1, T2 and the condensation test), infrared greenhouse (T3, IR flag), photochemistry (T1 threshold, Boss 14), gas optics (CIA flags, densities), condensates, clouds, hazes (T5, T6, T9, NH4SH and H2SO4 entries), radiation shielding (column mass, formulas 1 to 4), volcanic outgassing and escape (Henry constants, scale height, critical constants).
