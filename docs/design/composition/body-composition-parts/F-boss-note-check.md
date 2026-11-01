# Part F: check of Boss's additional research note (2026-10-17, 1:22 AM EDT)

Boss posted a 14-item research note (Tier 1 to Tier 5) at 1:22 AM EDT. This part checks the items that touch body composition against primary or open sources, relabels them with the same labels as the other parts (READ, 2ND-HAND, COMPUTED, SNIPPET, RECALLED, NOT REACHED), and says where the rest belongs (mostly the atmosphere lane). The note itself is a summary with few numbers and carries `[cite: n]` tags that point to no readable source; it is treated as a pointer list, not as a source.

Files: `parts/F-files/fray_schmitt_table2_transcribed_and_checked.py` (the transcribed polynomial table and the checks) and the two cropped images of the table it was read from.

## 1. Verdicts

| # | Item in the note | Verdict | Label | Evidence |
|---|---|---|---|---|
| 1 | Water-ice sublimation, [R1345] form with e0 to e6, P_triple 6.116577e-3 bar, T_triple 273.16 K | **Correct as written** | COMPUTED | Evaluated against IAPWS R14-08 (`iapws` 1.5 package): ratio 0.9997 to 1.0006 from 60 K to 273.16 K (611.657 Pa at the triple point, 0.1626 Pa at 200 K, 6.0956e-6 Pa at 150 K, 1.0857e-14 Pa at 100 K). Against [R3217] it is within 0.3 percent down to 80 K and 1 percent at 60 K (Murphy-Koop is itself an extrapolation there). The coefficients are the same as in the NASA/Icarus supplement (below) to the digits given |
| 2 | [R1445] polynomial `ln(P_sat[bar]) = A0 + sum A_i/T^i` for N2, CO, CH4, CO2, SO2, H2S, NH3 | **Form and units confirmed; coefficients checked only through a secondary table; primary paper not reachable** | COMPUTED from a 2ND-HAND table | Fray and Schmitt (Planet. Space Sci. 57, 2053) is on ScienceDirect (refused) and the HAL/ISTEX copy is behind an institutional login. The coefficients below were read from Table 2 of the supplementary material of the Arrokoth/Pluto ices paper (NTRS 20210019502, image, because the text layer garbles exponents) |
| 3 | "Secondary preprints carry sign inversions and truncated precision" | **Partly supported** | COMPUTED | No sign inversion found in the transcribed Table 2: every polynomial gives the right triple-point pressure (section 2). The text layer of that supplement does contain typesetting errors (SO2 A1 printed "-3.5.08x10+3", C2H4 A1 "-2.206 10+3", Kr A6 "-s1.210x10+8"); read from the image these are -3.508e3, -2.206e3, -1.210e8. What the check does show: the two-parameter fits (the older, more familiar set) are off by large factors in the cold range: CO is 6 times too high at 40 K and 30 times at 30 K against Fray and Schmitt, and H2S is 2.5 times too low at its triple point against NIST (section 2) |
| 4 | Validity ranges in the note (N2 20 to 63 K, CO 20 to 68 K, CH4 30 to 90 K, CO2 70 to 195 K, SO2 100 to 200 K, H2S 80 to 180 K, NH3 100 to 195 K) | **Not verified** | NOT REACHED | The supplement gives the polynomials without ranges. Where a species has two polynomials the crossing points tell which is which (N2-1 crosses N2-2 at 36.1 K, the alpha-beta transition at 35.6 K, so N2-1 is the alpha phase) but CO-1 and CO-2 cross twice (39 K and 61 K) and CO2/H2S pairs do not cross, so the selection rule needs the primary paper |
| 5 | Magma ocean: Fe3+/(Fe3+ + Fe2+) = 0.037 +/- 0.005 at 2173 K and 1 bar, fO2 = IW + 0.5, CO-rich gas, CO2-N2 atmosphere after cooling, Earth and Venus alike at the start | **Confirmed** | READ | [R4154], Sci. Adv. 6, eabd1387, read on PMC (PMC7688334): the abstract and Fig. 1 text give exactly these figures, f_H2O/f_H2 = 1.05 and f_CO2/f_CO = 0.21 at 2173 K and IW + 0.5 |
| 6 | "Water retains over 99 percent in the magma ocean; carbon has extremely low solubility and goes to the gas" | **Overstated for carbon** | READ | The paper computes that about 99 percent of the bulk silicate Earth's water and **about 35 percent of its CO2** are initially dissolved at 2173 K, and N2 stays almost entirely in the gas. The numbers use the BSE inventories of H, C and N (Earth's), so they are Earth's case, not a law for other planets |
| 7 | ICRP 103 radiation and tissue weighting factors (photons and electrons 1, protons 2, alpha and heavy ions 20, neutrons a continuous function; bone marrow, lung, colon 0.12, gonads 0.08) | **Values right, table mapping wrong, neutron range off** | RECALLED | In ICRP 103 the tissue weights w_T belong to tissues, not to particle types (the note's table pairs each particle with one tissue). The neutron weighting runs from 2.5 to about 20.7 (peak near 1 MeV), not 5 to 20. ICRP 103 itself was not read here |
| 8 | Comet volatile ranges: CO/H2O 0.15 to 30 percent, CO2/H2O 2 to 30, CH4 0.2 to 1.5, NH3 0.2 to 1.0 | **Consistent with what part D read** | READ (D 2.6) | Harrington [R3535]: CO/H2O median 3.3, range 0.3 to 26; CO2/H2O median 12, range 2 to 30. 67P [R3835]: CH4 0.34, NH3 0.67 percent, both inside the ranges. The note's lower CO bound 0.15 percent is not in what was read. [R3206] not reached |
| 9 | [R4352]: M_Z against M_p with 0.3 to 0.5 dex scatter | **Exponent matches; scatter does not** | READ (C 2.6) | M_Z = 57.9 (M/MJ)^0.61 Earth masses with scatter factor 1.82 (0.26 dex). The wider scatter in the note is not in the paper as read; keep 0.26 dex |
| 10 | Mars dust single scattering albedo 0.92 to 0.95 near 670 nm, [R2645, R4395] under "Volcanic SO2 aerosol" | **Misattributed row** | RECALLED | [R2645, R4395] are Mars dust papers (the atmosphere sources check lists them for Mars dust, gap 1). They give no SO2 aerosol optics. Atmosphere lane |
| 11 | Tiers 2 and 3 other items (giant visible absorption, Titan haze fractal aggregates, Parmentier-Guillot profile, Myhre forcing, King factors, trapped-particle models, lava-planet vaporisation, brines, clathrates, photochemistry, licensing) | **Not checked here** | n/a | Atmosphere lane or later items (section 3). The Rayleigh item (King factor F = (6 + 3 rho)/(6 - 7 rho), N2 F = 1.034 + 3.17e-4 / lambda^2) is a standard form; the atmosphere sources check already covers Bates, Bucholtz and Sneep and Ubachs |

## 2. The vapour-pressure table, transcribed and checked

Source of the numbers: Table 2 of the supplementary material of the Icarus 2019 paper on Arrokoth and Pluto ices (NTRS 20210019502, page 2, read from the image). It reproduces [R1445]. Unit: P in bar, T in K, `ln P = A0 + A1/T + A2/T^2 + ... + A6/T^6`.

| Polynomial | A0 | A1 | A2 | A3 | A4 | A5 | Check at the NIST triple point (T, P) | Ratio table/NIST | Label |
|---|---|---|---|---|---|---|---|---|---|
| N2-1 (alpha) | 12.40 | -8.074e2 | -3.926e3 | 6.297e4 | -4.633e5 | 1.325e6 | below 35.6 K only | n/a (crosses N2-2 at 36.1 K) | COMPUTED |
| N2-2 (beta) | 8.514 | -4.584e2 | -1.987e4 | 4.800e5 | -4.524e6 | 0 | 63.14 K, 0.1252 bar | 0.970 | COMPUTED vs READ |
| CO-1 / CO-2 | 10.43 / 10.25 | -7.213e2 / -7.482e2 | -1.074e4 / -5.843e3 | 2.341e5 / 3.939e4 | -2.392e6 / 0 | 9.478e6 / 0 | 68.1 K: 0.160 and 0.155 bar (NIST triple pressure for CO not extracted) | n/a | COMPUTED |
| CH4-1 | 10.51 | -1.110e3 | -4.341e3 | 1.035e5 | -7.910e5 | 0 | 90.67 K, 0.1169 bar | 1.013 | COMPUTED vs READ |
| CO2-1 / CO2-2 | 14.76 / 18.61 | -2.571e3 / -4.154e3 | -7.781e4 / 1.041e5 | 4.325e6 / 0 | -1.207e8 / 0 | 1.350e9 / 0 | 216.58 K, 5.185 bar | 0.959 / 1.004 | COMPUTED vs READ |
| NH3-1 | 15.96 | -3.537e3 | -3.310e4 | 1.742e6 | -2.995e7 | 0 | 195.4 K, 0.0606 bar | 1.008 | COMPUTED vs READ |
| H2S-1 / H2S-2 | 12.98 / 8.933 | -2.707e3 / -7.260e2 | 0 / -3.504e5 | 0 / 2.724e7 | 0 / -8.582e8 | 0 | 187.66 K, 0.232 bar | 1.016 / 1.005 | COMPUTED vs READ |
| SO2-1 | 15.60 | -3.508e3 | -9.401e4 | 4.152e6 | -6.946e7 | 0 | 197.64 K, 0.0167 bar | 1.029 | COMPUTED vs READ |

NIST triple points were read from the NIST WebBook phase-change pages (N2 63.14 K and 0.1252 bar; CH4 90.67 K and 0.1169 bar; CO2 216.58 K and 5.185 bar; H2S 187.66 K and 0.232 bar; SO2 197.64 K and 0.0167 bar; NH3 0.0606 bar with T about 195.4 K, the page lists a fusion temperature of 194.95 K). The agreement (at most 3 percent, one polynomial 4 percent) shows the units and the exponent layout are right; it does **not** verify the cold end (below 50 K), which is the end that matters for Pluto-type bodies. For that end the table agrees with the older fits for CH4 (within 12 percent from 40 K to 90 K) and for CO2 above 130 K, and disagrees for CO and NH3, where the note's claim that the Prialnik fits fail in the cold is what the check finds.

## 3. What this changes in the design

New items, grouped for the TODO sweep (they attach to existing clusters in section 10 of the main draft):

| Title | Scope | Attaches to | Depends on | Size |
|---|---|---|---|---|
| F1 Vapour-pressure function | One `p_sat(species, T)` for the ices (Fray and Schmitt polynomials, water by Feistel and Wagner), feeding the volatile-ice table (C), comet and KBO volatile retention (D), freeze-out in the atmosphere chain (GEN.422) and the hydrosphere | Cluster 2 (materials table) and the atmosphere chain | Boss to supply or the lane to read the primary paper; until then the transcribed table above is the working copy | small |
| F1t Test | Triple-point anchors from NIST (3 percent) and the IAPWS check for water (0.1 percent) | Cluster 11 | F1 | small |
| F2 Magma-ocean redox state | Store a mantle oxidation state per rocky planet (Fe3+/total Fe or delta-IW), tied to B's mantle FeO; relations f_H2O/f_H2 and f_CO2/f_CO from fO2 and T (Sossi Fig. 1); H, C, N inventories and solubility laws for the first outgassed gas. The value for Earth (IW + 0.5) is a calibration target only; the draw range is a parametrisation (Sossi's Fig. 1 explores IW - 4.5 to + 3.1) | Cluster 7 and the atmosphere chain (GEN.410 onward) | B stored mantle oxides; atmosphere lane for the species table | medium |
| F3 Lava-planet vaporisation (later) | Silicate vapour equilibria over molten surfaces (SiO2, Na2O, FeO to gas) for ultra-hot rocky planets; [R3945] and LavAtmos named in the note, neither read here | B / atmosphere | F2 | medium |
| F4 Clathrates and brines (later) | Clathrate dissociation pressures and eutectic brines for icy moons and Titan-type crusts ([R0779] named, neither read) | Cluster 8 (hydrosphere) with D ocean moons | C, D tables | medium |
| F5 Licence column | Every data file loaded into the repo states its source and licence (NIST, PubChem, USGS public domain; published compilations with citation); the materials-table loader refuses rows without one | Cluster 2 | none | small |

The remaining items in the note (giant visible opacities, aerosols, infrared optical depth, photochemistry, Rayleigh, trapped particles, brine acid vapour pressures) belong to the atmosphere lane and are not duplicated here; the note's Tier 2 matches gaps 1 to 3 and 5 of the atmosphere sources check, which it names but gives no numbers for.

## 4. Gaps still open (Boss must research)

1. [R1445] itself (Planet. Space Sci. 57, 2053; ScienceDirect refused; HAL copy behind a login): the validity range of each polynomial and the selection rule between the two sets for N2, CO, CO2 and H2S; and an independent check of the cold end.
2. [R4154] for planets other than Earth: how fO2 varies with planet mass, composition and water ([R1491] and 2021 named in the note; [R1953]); no distribution exists in anything read.
3. ICRP Publication 103 (the neutron weighting function and the tissue table), for the surface-dose conversion of part E.
