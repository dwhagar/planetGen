# Part A: Stars and stellar remnants (composition research, 2026-10-17)

Lane A (STARS) of the v9.0 body-composition research. Read-only research: nothing in the repo, `docs/` or `TODO.md` was touched. Scratch scripts and downloaded catalogues are outside the repo; every COMPUTED row states its inputs and method. Labels are the brief's: READ, 2ND-HAND, COMPUTED, SNIPPET, RECALLED, NOT REACHED.

Abbreviations for sources (arXiv id in brackets): A21 Asplund, Amarsi, Grevesse 2021 [2105.01661]; A09 Asplund, Grevesse, Sauval, Scott 2009 [0909.0948]; L19 Lodders 2019 [1912.00844]; L09 Lodders, Palme, Gail 2009 [0901.1149]; BF16 Brewer and Fischer 2016 [1608.06286] with catalogue Brewer et al. 2016 [1606.07929]; B14 Bensby, Feltzing, Oey 2014 [1309.2631] (CDS table c3); N14 Nissen et al. 2014, A&A 568 A25 (CDS J/A+A/568/A25); BB16 Battistini and Bensby 2016 [1511.00966] (CDS J/A+A/586/A49); H15 Hayden et al. 2015 [1503.02110]; Hk14 Hinkel et al. 2014 [1405.6719]; T15 Thiabaud et al. 2015 [1507.01343].

## 1 Scope and interfaces

**Covers.** The composition of every kind of star and stellar remnant, and the statistics needed to generate it from quantities the project already stores (mass, age, [Fe/H], galactic position, population, spectral type, luminosity class): (1) solar abundance sets and how to turn a [Fe/H] into X, Y, Z and per-element mass fractions; (2) how element ratios vary with [Fe/H] and population in the Galaxy, and the resulting Fe/Si, Mg/Si, C/O distributions; (3) interior and surface composition by type and phase; (4) white dwarfs, neutron stars, black holes, brown dwarfs.

**The solar set is a reference scale, not a preset.** Every star's composition is `solar-set abundance x 10^([X/H])` with [X/H] = [Fe/H] + [X/Fe] drawn from the Galactic relations below. The Sun's own numbers appear only as the zero point of that scale and as a test target (section 4). No Solar System body is a default.

**Takes from other parts / the generator.** Mass, age, population label (young, intermediate, old, bulge, halo, cluster, satellite, nuclear), galactic radius R and height z, [Fe/H] (`metallicity_fe_h`, schema v95), spectral type and Yerkes class.

**Hands over.**
- To the planet lanes (rocky, volatile-rich, giants): the host's birth mass fractions X, Y, Z and the element-ratio vector (`[Mg/Fe]`, `[Si/Fe]`, `[O/Fe]`, `[C/Fe]` -> Fe/Si, Mg/Si, C/O, plus O/Si, which the bulk-composition chain GEN.425 to GEN.430 asks for). Number ratios and mass ratios are both given (conversion factors 1.988, 0.8654, 0.7507 for Fe/Si, Mg/Si, C/O, section 2.3).
- To the ground-radiation/radiogenic lane: birth abundances of Th, U, K (and Eu as the Th tracer) per unit Si or per unit mass of rock, with decay laws.
- To the atmosphere lane (interface only, nothing designed here): a star's photospheric composition is the stored X, Y, Z plus the [X/Fe] vector; spectral-class chemistry of stellar atmospheres belongs to that lane.
- To the compact-remnant design notes: remnant composition class (WD core type, NS crust sequence), radius laws, mass distributions.

**Current code (checked in `/home/user/planetGen/src/planetgen`).** Only [Fe/H] is stored per star (`generation/star.py` `metallicity_fe_h`, drawn by `physics/metallicity.py:sample_fe_h`). Mean is population offset + disc gradient, normal scatter, clamp to `METALLICITY_RANGE_FE_H = (-3.5, 0.6)` (`tuning.py` lines 661-700): young (0.00, sigma 0.10), intermediate (0.00, 0.15), old (-0.20, 0.25), bulge (-0.05, 0.40, no gradient), nuclear cluster (+0.20, 0.40), nuclear disc (+0.20, 0.30), halo (-1.6 inside 25 kpc, -2.2 outside, sigma 0.5); gradient 0.06 dex/kpc zero at the Sun's distance, marked "recalled, to verify" in the source. [Fe/H] is used only for planet occurrence (`solids_factor` = 10^[Fe/H], giant factor 10^(2[Fe/H])) and disk surface density (`formation.py`). There is no X, Y, Z, no alpha enhancement, no element ratios, no age-metallicity link. White dwarfs already use the Kalirai 2008 IFMR (`tuning.py` 600-610, verified below); NS mass 1.1-2.2 Msun and radius 10-13 km, BH 5-20 Msun spin 0-0.998 are plain ranges (`tuning.py` 2023-2091).

## 2 Findings

Counting rule for section 7: one label per table row; a row that carries two labels (inputs READ, derived quantity COMPUTED) is counted once under each.

### 2.1 Solar abundance sets and what each is good for

| Set | Quantity | Value | Unit | Source | Label |
|---|---|---|---|---|---|
| A21 photosphere | X, Y, Z | 0.7438, 0.2423, 0.0139 | mass fraction | A21 Table 2 log eps summed with atomic weights; reproduces the paper's X, Y, Z | READ+COMPUTED |
| A21 proto-solar | X, Y, Z | 0.7121, 0.2725, 0.0154 | mass fraction | A21 log eps + 0.070 dex (He) and +0.064 dex (all metals) diffusion correction from A21 sect. 5, re-summed | COMPUTED |
| A09 photosphere | X, Y, Z | 0.7374, 0.2492, 0.0134 | mass fraction | A09 Table 1 log eps summed (He 10.93 is helioseismic) | COMPUTED |
| A09 diffusion correction | He, metals | +0.05 / +0.04 | dex | A09 text, 0.04-0.05 dex | READ |
| L09 present-day | Z/X, X, Y, Z, A(He) | 0.0191, 0.7390, 0.2469, 0.0141, 10.925 | mass fraction | L09 Table 7 (X from Basu and Antia helioseismology) | READ |
| L09 proto-solar | Z0/X0, X0, Y0, Z0, A(He)0 | 0.0215, 0.7112, 0.2735, 0.0153, 10.986 | mass fraction | L09 Table 9; offsets A(He)0 = A(He)+0.061, A(X)0 = A(X)+0.053 dex | READ |
| L19 present-day | X, Y, Z | 0.7389, 0.2463, 0.0148 | mass fraction | L19 Tables 7-8 log eps summed; reproduces L19 Table 7 | READ+COMPUTED |
| GS98 (Grevesse and Sauval 1998) | X, Y, Z, Z/X | 0.7347, 0.2483, 0.0169, 0.0231 | mass fraction | L09 Table 7 (primary Space Sci. Rev. 85, 161 not opened) | 2ND-HAND |
| GS98 proto-solar | X0, Y0, Z0 | 0.7086, 0.2750, 0.0163 | mass fraction | L09 Table 9 | 2ND-HAND |
| Lodders 2003 | log eps table (Fe/Si mass 1.69 quoted in the GEN.425 sources check) | not opened: the IOP page returned a bot-check page and was not routed around; L09 is its documented successor and gives 1.692 for the same ratio | | Lodders 2003, ApJ 591, 1220 | NOT REACHED |
| Solar initial helium | Y_ini | 0.278 +- 0.006 (standard models); 0.273 +- 0.006 (extra mixing) | mass fraction | Serenelli and Basu 2010 [1006.0244] abstract | READ |
| Solar surface helium (helioseismic) | Y_surf | 0.2485 +- 0.0035 | mass fraction | same abstract | READ |

Which to use: A21 (log eps per element, 83 elements, Table 2) as the scale; L09 and L19 as the cross-check and as the CI-chondrite link. Photosphere values describe the Sun now; the birth composition of a star is the "proto" row.

Differences that matter for planets: L09/L19 O is 0.04 dex higher than A21 and C is 0.07 lower in L09; Ne ranges 7.93-8.15 and Ar 6.38-6.50 across sets (the noble gases are not measured in the photosphere; they come from corona/SEP/B-star data), so Ne and Ar mass fractions are the least certain entries of Z.

### 2.2 Per-element table (log eps on H = 12, and mass fractions)

Mass fractions are COMPUTED from the A21 log eps (Table 2) times standard atomic weights (`periodictable`), proto-solar rows add the A21 sect. 5 offsets. The log eps columns are READ from the cited tables (A09 Table 1, L19 Tables 7-8, L09 Table 4). The CI column is L09 Table 4 for the Orgueil CI chondrite on the same H = 12 scale (use for ratios only; volatile elements are depleted). Th and U in the photosphere are upper limits in L09; A21 and L19 take them from meteorites. Row label: eps READ, mass fractions COMPUTED.

| El | A21 log eps | A21 mass fraction (photosphere) | A21 mass fraction (proto-solar) | A09 log eps | L19 log eps | L09 log eps | CI (L09) | Label |
|---|---|---|---|---|---|---|---|---|
| H | 12.000 | 7.4382e-01 | 7.1207e-01 | 12.00 | 12.000 | 12.000 | n/a | READ+COMPUTED |
| He | 10.914 | 2.4230e-01 | 2.7252e-01 | 10.93 | 10.924 | 10.925 | n/a | READ+COMPUTED |
| C | 8.460 | 2.5561e-03 | 2.8356e-03 | 8.43 | 8.470 | 8.390 | 7.410 | READ+COMPUTED |
| N | 7.830 | 6.9880e-04 | 7.7519e-04 | 7.83 | 7.850 | 7.860 | 6.280 | READ+COMPUTED |
| O | 8.690 | 5.7823e-03 | 6.4144e-03 | 8.69 | 8.730 | 8.730 | 8.420 | READ+COMPUTED |
| Ne | 8.060 | 1.7097e-03 | 1.8966e-03 | 7.93 | 8.150 | 8.050 | n/a | READ+COMPUTED |
| Na | 6.220 | 2.8154e-05 | 3.1232e-05 | 6.24 | 6.270 | 6.300 | 6.290 | READ+COMPUTED |
| Mg | 7.550 | 6.3636e-04 | 7.0592e-04 | 7.60 | 7.520 | 7.540 | 7.550 | READ+COMPUTED |
| Al | 6.430 | 5.3589e-05 | 5.9447e-05 | 6.45 | 6.420 | 6.470 | 6.450 | READ+COMPUTED |
| Si | 7.510 | 6.7063e-04 | 7.4394e-04 | 7.51 | 7.510 | 7.520 | 7.530 | READ+COMPUTED |
| S | 7.120 | 3.1187e-04 | 3.4596e-04 | 7.12 | 7.150 | 7.140 | 7.170 | READ+COMPUTED |
| Ar | 6.380 | 7.0717e-05 | 7.8448e-05 | 6.40 | 6.500 | 6.500 | n/a | READ+COMPUTED |
| K | 5.070 | 3.3897e-06 | 3.7603e-06 | 5.03 | 5.070 | 5.120 | 5.100 | READ+COMPUTED |
| Ca | 6.300 | 5.9008e-05 | 6.5459e-05 | 6.34 | 6.270 | 6.330 | 6.310 | READ+COMPUTED |
| Fe | 7.460 | 1.1885e-03 | 1.3184e-03 | 7.50 | 7.450 | 7.450 | 7.470 | READ+COMPUTED |
| Ni | 6.200 | 6.8643e-05 | 7.6147e-05 | 6.22 | 6.200 | 6.230 | 6.220 | READ+COMPUTED |
| Eu | 0.520 | 3.7132e-10 | 4.1191e-10 | 0.52 | 0.510 | 0.520 | 0.530 | READ+COMPUTED |
| Th | 0.030 | 1.8347e-10 | 2.0353e-10 | 0.02 | 0.040 | <0.08 (limit) | 0.080 | READ+COMPUTED |
| U | -0.540 | 5.0656e-11 | 5.6194e-11 | -0.54 | -0.540 | <-0.47 (limit) | -0.520 | READ+COMPUTED |

Shares of the metal mass Z (A21 photosphere):

| Quantity | Value | Unit | Source | Label |
|---|---|---|---|---|
| Share of Z: O, C, Ne, Fe | 0.416, 0.184, 0.123, 0.086 | fraction of Z | A21 mass fractions / Z | COMPUTED |
| Share of Z: N, Si, Mg, S | 0.050, 0.048, 0.046, 0.022 | fraction of Z | same | COMPUTED |

### 2.3 Element ratios of the reference sets (what planets inherit)

Number ratio = 10^(eps_a - eps_b); mass ratio = number ratio x (A_a/A_b); factors Fe/Si 1.9884, Mg/Si 0.8654, C/O 0.7507 (atomic weights 55.845, 28.0855, 24.305, 12.011, 15.999). T15 quotes mass ratios, so compare mass columns with its envelope.

| Set | Fe/Si number | Fe/Si mass | Mg/Si number | Mg/Si mass | C/O number | C/O mass | Label |
|---|---|---|---|---|---|---|---|
| A21 | 0.891 | 1.772 | 1.096 | 0.949 | 0.589 | 0.442 | COMPUTED |
| A09 | 0.977 | 1.943 | 1.230 | 1.065 | 0.550 | 0.413 | COMPUTED |
| L19 | 0.871 | 1.732 | 1.023 | 0.886 | 0.550 | 0.413 | COMPUTED |
| L09 | 0.851 | 1.692 | 1.047 | 0.906 | 0.457 | 0.343 | COMPUTED |
| CI chondrite (L09 Table 4) | 0.871 | 1.732 | 1.047 | 0.906 | 0.098 | 0.073 | COMPUTED |

The spread among sets for the Sun itself is 1.69-1.94 in Fe/Si (mass), 0.89-1.07 in Mg/Si, 0.34-0.44 in C/O. This is the floor on how well any "scaled solar" statement can be made.

### 2.4 Helium, and X, Y, Z from [Fe/H]

| Quantity | Value | Unit | Source | Label |
|---|---|---|---|---|
| Primordial Y_P (Planck 2018, PArthENoPE) | 0.24672 (+0.00011 / -0.00012 from omega_b; +-0.00061 incl. theory, 95%) | mass fraction | Planck 2018 VI [1807.06209] eq. 72a | READ |
| Primordial Y_P (PRIMAT) | 0.24714 (+-0.00049, 95%) | mass fraction | same, eq. 72b | READ |
| Delta Y / Delta Z, solar calibration, A21 proto, Y_P = 0.24672 | (0.2725 - 0.24672) / 0.0154 = 1.68 | dimensionless | A21 proto values above + Planck | COMPUTED |
| Same, L09 proto (0.2735, 0.0153) | 1.75 | dimensionless | L09 Table 9 + Planck | COMPUTED |
| Same, Y_ini = 0.278 | 2.03 | dimensionless | Serenelli and Basu + A21 Z0 | COMPUTED |
| Observed Delta Y / Delta Z, K dwarfs | 2.1 +- 0.9 | dimensionless | [astro-ph/0703766] abstract | READ |
| no source | 2.1 +- 0.4 | dimensionless | quoted | 2ND-HAND |
| Closed form | Y = Y_P + g Z; with r = Z/X: X = (1 - Y_P) / (1 + (1+g) r), Z = r X, Y = 1 - X - Z | | algebra; g = 1.68 reproduces the A21 proto X, Y, Z at r = 0.02163 | COMPUTED |
| r from [Fe/H] and [X/Fe] | r = r0 * sum_i s_i 10^([Fe/H] + [X_i/Fe]) with s_i the A21 proto mass fraction of metal i divided by the proto Z0, r0 = 0.02163 | | algebra on 2.2 | COMPUTED |

Validity: [Fe/H] here is a birth value. A main-sequence star's photosphere loses 0.04-0.07 dex of He and metals to diffusion over its life (A09, A21, L09 offsets), so the stored [Fe/H] should be read as initial.

### 2.5 What a given [Fe/H] means for Z, with real alpha enhancement

Method (COMPUTED): the fitted mean [X/Fe] relations of section 2.7 for C, O, Mg, Si, Ca (thin-disk sequence for [Fe/H] > -0.8, high-alpha sequence otherwise); Ne follows O (Ne column) or Mg (second column), S, Ar, Ti follow Ca, everything else [X/Fe] = 0; X, Y, Z from the closed form in 2.4.

| [Fe/H] | Sequence | X | Y | Z | Z / (Z0 10^[Fe/H]) with Ne=O | Same with Ne=Mg | Label |
|---|---|---|---|---|---|---|---|
| +0.2 | thin | 0.6970 | 0.2820 | 0.02105 | 0.862 | 0.884 | COMPUTED |
| 0.0 | thin | 0.7129 | 0.2720 | 0.01511 | 0.981 | 0.988 | COMPUTED |
| 0.0 | high-alpha | 0.6998 | 0.2802 | 0.02000 | 1.298 | 1.274 | COMPUTED |
| -0.5 | thin | 0.7353 | 0.2580 | 0.00672 | 1.380 | 1.386 | COMPUTED |
| -0.5 | high-alpha | 0.7267 | 0.2634 | 0.00995 | 2.043 | 1.962 | COMPUTED |
| -1.0 | high-alpha | 0.7410 | 0.2544 | 0.00459 | 2.982 | 2.760 | COMPUTED |
| -2.0 | high-alpha | 0.7520 | 0.2475 | 0.00047 | 3.027 | 2.798 | COMPUTED |
| -3.0 | high-alpha | 0.7532 | 0.2468 | 0.00005 | 3.031 | n/a | COMPUTED |

Generic comparison with a single [alpha/Fe] (alpha = O, Ne, Mg, Si, S, Ar, Ca, Ti raised together):

| Case | Z / (Z0 10^[Fe/H]) | Source | Label |
|---|---|---|---|
| [Fe/H] -0.5 with [alpha/Fe] +0.3; -1.0 and -2.0 with +0.4 | 1.674, 2.029, 2.034 | A21 mass fractions re-summed | COMPUTED |
| Factor f = 0.638 x 10^[alpha/Fe] + 0.362, same three cases | 1.635, 1.965, 1.965 | formula from memory | RECALLED |

Reading: the project's rule "solids scale as 10^[Fe/H]" underestimates the metal mass of metal-poor high-alpha stars by 1.7 to 3 times, mostly because O is +0.6 dex there. The 3.0 depends on the NLTE oxygen trend (N14) and on Ne following O; the generic alpha model gives 2.0. Keep both numbers in the work item that decides the Ne and O choice.

### 2.6 Galactic abundance patterns: metallicity distributions

| Quantity | Value | Unit | Source | Label |
|---|---|---|---|---|
| Midplane (abs(z) < 0.5 kpc) mean [Fe/H], R = 3-5, 5-7, 7-9, 9-11, 11-13, 13-15 kpc | +0.08, +0.11, +0.01, -0.11, -0.21, -0.31 | dex | H15 Table 1 | READ |
| Midplane sigma[Fe/H] same bins | 0.24, 0.22, 0.20, 0.19, 0.18, 0.18 | dex | H15 Table 1 | READ |
| Midplane skewness same bins | -1.68, -1.26, -0.53, -0.02, +0.17, +0.47 | dimensionless | H15 Table 1 | READ |
| 0.5 < abs(z) < 1 kpc mean, same bins | -0.19, -0.12, -0.10, -0.15, -0.23, -0.33 | dex | H15 Table 1 | READ |
| 1 < abs(z) < 2 kpc mean, same bins | -0.42, -0.36, -0.31, -0.29, -0.29, -0.39 | dex | H15 Table 1 | READ |
| Midplane radial gradient, linear fit of the bin means to bin centres 6-14 kpc | -0.053 | dex/kpc | H15 Table 1 means | COMPUTED |
| Same fit including the 3-5 kpc bin | -0.043 | dex/kpc | same | COMPUTED |
| Gradient above the plane | 0.5-1 kpc: -0.015; 1-2 kpc: +0.005 | dex/kpc | same | COMPUTED |
| Thin-disk [Fe/H] vs age, local sample | -0.027 dex/Gyr, scatter 0.265 dex | | B14 table c3, thin stars (TD/D < 0.1, Teff > 5400 K) linear fit | COMPUTED |
| Mean ages | thin 5.7 Gyr, thick 10.9 Gyr | Gyr | B14 table c3 | COMPUTED |
| Fraction of stars kinematically thick (TD/D>1) by age 0-2, 2-4, 4-6, 6-8, 8-10, 10-12, 12-15 Gyr | 0.08, 0.11, 0.17, 0.35, 0.33, 0.54, 0.68 | fraction | B14; sample over-represents the thick disk, shape only | COMPUTED |
| Chemically high-alpha fraction (mean of [Mg,Si,Ca,Ti/Fe] > 0.1) vs [Fe/H] bins -1..-0.7, -0.7..-0.5, -0.5..-0.35, -0.35..-0.2, -0.2..-0.1, -0.1..0, 0..0.1, 0.1..0.2 | 1.00, 0.84, 0.66, 0.41, 0.18, 0.05, 0.01, 0.00 | fraction | B14 | COMPUTED |
| Logistic fit to the above: P(high-alpha) = 1 / (1 + exp(([Fe/H] + 0.350) / 0.132)) | f0 = -0.350, w = 0.132 | dex | B14 fit | COMPUTED |
| Bulge [Fe/H] distribution (90 microlensed dwarfs and subgiants) | peaks at -1.09, -0.63, -0.20, +0.12, +0.41 | dex | Bensby et al. 2017 [1702.02971] abstract | READ |
| Bulge ages | at [Fe/H] > 0 more than 35% younger than 8 Gyr; at [Fe/H] <~ -0.5 most stars 10 Gyr or older; star-formation episodes 3, 6, 8, 11 Gyr ago | | same | READ |
| Bulge alpha knee | tentatively at slightly higher [Fe/H] than in the local thick disk | | same | READ |
| Very metal-poor stars carbon-rich ([C/Fe] > +0.7) | 20-30% for -4 < [Fe/H] < -2; most stars below -4 | fraction | Frebel and Norris 2015 [1501.06921] sect. on C-rich stars | READ |
| CEMP stars below [Fe/H] -3 | almost all CEMP-no | | same | READ |

### 2.7 [X/Fe] versus [Fe/H] by population (drawable models)

All fits here are COMPUTED by me from CDS catalogue rows (B14 table c3 for Mg, Si, Ca, Ti, Ni, Cr, Na, Al and Bensby's own O; N14 tables 4-5 for C and NLTE O/C; BB16 table 4 for Eu), not read off figures. Forms: thin = `c + b min([Fe/H] - k, 0)` ("plateau above k, rising with slope -b below it"); high-alpha = `P + s max([Fe/H] - k, 0)` ("plateau P below k, falling with slope s above"); sd is the residual scatter after the fit, dex. Population assignment: thin = Pop D (N14) or TD/D < 0.1 (B14); high-alpha = thick-disk plus high-alpha halo (N14 pops C and A) or B14 stars with TD/D > 1 and not halo.

| Element | Thin disk mean | sd | High-alpha mean | sd | n thin / high | Source | Label |
|---|---|---|---|---|---|---|---|
| [Mg/Fe] | 0.006 - 0.343 min([Fe/H] - 0.03, 0) | 0.066 | 0.357 - 0.373 max([Fe/H] + 0.73, 0) | 0.073 | 256 / 127 | B14 | COMPUTED |
| [Si/Fe] | 0.017 - 0.166 min([Fe/H] + 0.07, 0) | 0.042 | 0.253 - 0.199 max([Fe/H] + 0.99, 0) | 0.073 | 256 / 126 | B14 | COMPUTED |
| [Ca/Fe] | -0.041 - 0.163 min([Fe/H] - 0.42, 0) | 0.038 | 0.282 - 0.251 max([Fe/H] + 0.99, 0) | 0.046 | 256 / 128 | B14 | COMPUTED |
| [Ti/Fe] | 0.015 - 0.181 min([Fe/H] - 0.01, 0) | 0.048 | 0.258 - 0.382 max([Fe/H] + 0.54, 0) | 0.065 | 256 / 128 | B14 | COMPUTED |
| [O/Fe] (B14 own O) | -0.497 - 0.516 min([Fe/H] - 0.95, 0) | 0.089 | 0.582 - 0.663 max([Fe/H] + 0.78, 0) | 0.113 | 217 / 102 | B14 | COMPUTED |
| [O/Fe] (N14 NLTE O I triplet) | -0.011 - 0.373 [Fe/H] (linear) | 0.052 | 0.606 - 0.495 max([Fe/H] + 0.91, 0) | 0.065 | 57 / 50 | N14 | COMPUTED |
| [C/Fe] (N14 NLTE) | -0.030 - 0.163 [Fe/H] (linear) | 0.044 | 0.191 - 0.521 max([Fe/H] + 0.10, 0) | 0.050 | 57 / 58 | N14 | COMPUTED |
| [Na/Fe] | 0.018 - 0.176 min([Fe/H] + 0.25, 0) | 0.056 | 0.042 - 0.063 max([Fe/H] + 0.39, 0) | 0.084 | 256 / 122 | B14 | COMPUTED |
| [Al/Fe] | 0.032 - 0.029 min([Fe/H] - 0.06, 0) | 0.089 | 0.225 - 0.298 max([Fe/H] + 0.48, 0) | 0.105 | 253 / 107 | B14 | COMPUTED |
| [Ni/Fe] | 0.029 + 0.074 min([Fe/H] - 0.51, 0) | 0.037 | -0.033 + 0.039 max([Fe/H] + 1.22, 0) | 0.047 | 256 / 127 | B14 | COMPUTED |
| [Cr/Fe] | 0.008 + 0.034 min([Fe/H] - 0.33, 0) | 0.028 | 0.005 - 0.006 max([Fe/H] + 0.78, 0) | 0.033 | 256 / 126 | B14 | COMPUTED |
| [Eu/Fe] (r-process, Th tracer) | 0.060 - 0.404 [Fe/H] (linear) | 0.077 | 0.446 - 0.290 max([Fe/H] + 0.93, 0) | 0.084 | 207 / 124 | BB16 table 4 matched on HIP to B14 membership | COMPUTED |
| Low-alpha halo [C/Fe], [O/Fe] | mean -0.043, +0.400 | 0.097, 0.113 | | | n = 33, 27 | N14 pop B | COMPUTED |

Notes. (a) The hinge fits are valid inside the sampled [Fe/H] ranges (thin -0.9 to +0.5, high-alpha -1.6 to +0.1); below -1.6 hold the plateau, above +0.5 hold the thin value at +0.5 (extrapolation, flagged). (b) C and O from N14 are used in preference to B14's O because N14 treats NLTE; the two O sets differ by up to 0.1 dex on the plateau. (c) The low-alpha halo has C and O only; Mg and Si for it are NOT fitted (gap, section 5). (d) The O plateau of +0.5 to +0.6 on the high-alpha sequence is the single largest lever on Z and on C/O at low [Fe/H].

Membership: P(high-alpha) logistic in [Fe/H] from section 2.6; the sample over-represents the thick disk, so the absolute local fraction must be set by a stored target and only the shape taken from the fit.

Residual correlations (matter because Mg/Si is a ratio of two correlated terms):

| Pair | Value | Source | Label |
|---|---|---|---|
| Thin disk, residual corr (O,Mg,Si): O-Mg 0.55, O-Si 0.78, Mg-Si 0.77 | sd 0.089, 0.064, 0.043; Mg-Si ratio residual 0.041 dex | B14 | COMPUTED |
| High-alpha, residual corr: O-Mg 0.66, O-Si 0.62, Mg-Si 0.71 | sd 0.106, 0.069, 0.075; Mg-Si ratio residual 0.055 dex | B14 | COMPUTED |
| BF16 (844 dwarfs), residual about a quadratic in [Fe/H]: C-O 0.56, C-Mg 0.69, C-Si 0.71, O-Mg 0.36, O-Si 0.37, Mg-Si 0.77 | sd 0.052, 0.074, 0.049, 0.046 (C, O, Mg, Si) | BF16 [Fe/H]-referenced table 9 | COMPUTED |

### 2.8 Resulting Fe/Si, Mg/Si, C/O distributions

BF16 values are number-ratio abundances [X/H]; I rebuilt log eps with A21 as reference (log eps_A21 + [X/H]) and converted to mass ratios. Selection: log g >= 3.5, 4800 < Teff < 6100 K, vsini <= 20 km/s, S/N > 100, finite C, O, Mg, Si, Fe: n = 844, [Fe/H] 5th/50th/95th percentile -0.36 / +0.08 / +0.32.

| Quantity | Median [P2.5, P97.5] (mass ratio) | Source | Label |
|---|---|---|---|
| Observed BF16 Fe/Si | 1.81 [1.17, 2.39] | BF16 catalogue, A21 scale | COMPUTED |
| Observed BF16 Mg/Si | 0.93 [0.79, 1.12] | same | COMPUTED |
| Observed BF16 C/O | 0.385 [0.24, 0.47] | same | COMPUTED |
| Thiabaud 2015 envelope (mass ratios) | Fe/Si 1.34-2.81, Mg/Si 0.66-1.85, C/O 0.12-0.75 | T15 [1507.01343] as quoted in the GEN.425 sources check | 2ND-HAND |
| Hinkel 2014 catalogue-to-catalog spread of the same star | mean 0.14 dex, median 0.11 dex over all elements | Hk14 sect. 3.2 | READ |
| Model draw, solar-neighbourhood-like [Fe/H] ~ N(-0.10, 0.22), high-alpha share 0.21, Fe/Si | 1.63 [1.09, 2.02] | model of 2.7 with the BF16 residual correlations of 2.7, n = 2e5 | COMPUTED |
| Same, Mg/Si | 1.02 [0.80, 1.36] | same | COMPUTED |
| Same, C/O | 0.40 [0.26, 0.54] | same | COMPUTED |
| Share of that draw inside the whole Thiabaud box | 0.88 (Fe/Si 0.88, Mg/Si 1.00, C/O 1.00 inside) | same | COMPUTED |
| By [Fe/H] bin, median Fe/Si, Mg/Si, C/O | -1.0..-0.6: 1.18, 1.31, 0.23; -0.6..-0.3: 1.40, 1.17, 0.32; -0.3..0: 1.63, 1.03, 0.40; 0..0.3: 1.70, 0.93, 0.45; 0.3..0.6: 1.70, 0.93, 0.51 | same | COMPUTED |
| Thin branch only, [Fe/H] ~ N(0.0, 0.2) (a BF16-like sample), Fe/Si, Mg/Si, C/O | 1.68 [1.36, 2.04]; 0.97 [0.78, 1.25]; 0.42 [0.32, 0.56] | same | COMPUTED |
| Project populations at the Sun's R (young, intermediate: thin only; old, bulge: logistic; halo: high-alpha branch). Fe/Si median [P2.5, P97.5] | young 1.69 [1.33, 2.06]; intermediate 1.68 [1.28, 2.05]; old 1.57 [0.98, 2.00]; bulge 1.61 [0.96, 2.02]; halo 1.00 [0.72, 1.44] | same, with project [Fe/H] means and sigmas | COMPUTED |
| Same, Mg/Si | young 0.96 [0.78, 1.21]; intermediate 0.97 [0.78, 1.24]; old 1.06 [0.81, 1.45]; bulge 1.01 [0.79, 1.44]; halo 1.21 [0.97, 1.52] | same | COMPUTED |
| Same, C/O | young 0.42 [0.34, 0.53]; intermediate 0.42 [0.32, 0.54]; old 0.38 [0.21, 0.53]; bulge 0.41 [0.18, 0.61]; halo 0.17 [0.13, 0.26] | same | COMPUTED |
| Share inside the Thiabaud box | young 0.97, intermediate 0.96, old 0.79, bulge 0.82, halo 0.05 | same | COMPUTED |

Reading. (0) Against BF16 the thin-branch model medians are 8% low in Fe/Si, 5% high in Mg/Si and 10% high in C/O, i.e. 0.03 to 0.04 dex, smaller than the 0.11 dex catalogue-to-catalogue offsets of Hk14: the model was fitted to B14 and N14 on their own zero points, BF16 is a third catalogue. (1) The homogeneous BF16 spread (Mg and Si residual sd about 0.05 dex each, correlated 0.77, so the Mg/Si ratio scatter is narrower still) is much tighter than the literature-compilation range in T15 because catalogue-to-catalogue offsets (Hk14: 0.11-0.14 dex) are larger than the real star-to-star scatter; the T15 envelope is therefore a compilation envelope, not a physical one. (2) The full draw (21% high-alpha, more metal-poor stars than BF16) has Fe/Si 1.63 vs 1.81 observed, Mg/Si 1.02 vs 0.93, C/O 0.40 vs 0.385. (3) Metal-poor high-alpha stars have lower Fe/Si (core mass fraction falls), higher Mg/Si and lower C/O: halo rocky bodies would be Si- and Mg-rich with less iron, and carbon planets are not expected at low [Fe/H] on this model (C/O stays below 0.3), except CEMP stars (2.6). (4) The Thiabaud box excludes the halo, as expected: it was built from the solar neighbourhood.

### 2.9 Radiogenic elements in stars

| Quantity | Value | Unit | Source | Label |
|---|---|---|---|---|
| Th follows Eu: [Th/Eu] at ZAMS vs age | +0.014 mean, sd 0.045; slope +0.0042 +-0.0016 dex/Gyr | dex | Botelho et al. 2019 [1810.10413] Tab. 4 and text | READ |
| [Th/Fe] at ZAMS vs age | +0.086 - 0.0009 t (+-0.0015); sd 0.047 | dex | same Tab. 4 | READ |
| [Th/H] at ZAMS vs age | +0.138 - 0.0117 t; mean +0.080, sd 0.058, range -0.117 to +0.257 | dex | same | READ |
| Th decay correction used | half-life 14.05 Gyr, t_sun = 4.56 Gyr | | same | READ |
| Th abundance spread in solar twins and analogues | 59% to 251% of solar (log eps Th) | | Unterborn et al. 2015 [1505.00280] abstract | READ |
| Photospheric Th, U in the Sun | Th 0.03, U -0.54 (log eps) | | A21 (meteoritic for U) | READ |
| Present atomic U/Th, Sun | 0.269 | | A21 eps | COMPUTED |
| Birth atomic U/Th, Sun (4.567 Gyr back; 238U 4.468 Gyr, 235U 0.7038 Gyr, 232Th 14.05 Gyr; 235U fraction 0.7204%) | 0.572 (log -0.242) | | A21 eps with decay factors 1.253 (Th), 2.031 (238U), 89.8 (235U) | COMPUTED |
| Birth atomic Th/Eu, Sun | 0.405 | | A21 eps, Th decay 1.253 | COMPUTED |
| Half-lives of 238U, 235U, 40K | 4.468, 0.7038, 1.248 | Gyr | standard values, not checked this session (232Th is READ above) | RECALLED |
| 40K isotopic fraction | 0.0117% | | not checked | RECALLED |
| [Eu/Fe] shape in the disc | Nd, Sm and Eu follow alpha-like trends | | BB16 abstract | READ |
| Remaining fraction of Th, 238U, 235U, 40K after tau = 1, 4.567, 8, 12 Gyr | Th 0.952/0.798/0.674/0.553; 238U 0.856/0.492/0.289/0.155; 235U 0.374/0.011/0.0004/~0; 40K 0.574/0.079/0.012/0.001 | | decay law, half-lives above | COMPUTED |

### 2.10 Stellar interiors and surfaces by type and phase

| Class / phase | Fact | Value | Source | Label |
|---|---|---|---|---|
| Main sequence, interior | Fully convective below | 0.35 Msun (minimum mass for a radiative core, all [M/H] -2 to 0) | Chabrier and Baraffe 1997 [astro-ph/9704118] | READ |
| Main sequence, mixed star | Interior is homogeneous in fully convective stars; surface = birth composition | | same | READ |
| MS, Sun as a model of the process | Central H fraction X_c: 0.708 (ZAMS), 0.340 (4.57 Gyr), 0.145 (6.5 Gyr), 0.008 (8 Gyr) | mass fraction | [astro-ph/0010346] sect. 3.4 | READ |
| MS, fit to the three post-ZAMS points | X_c / X_0 = 1 - (t / t_x)^1.15 with t_x = 8.0 Gyr for 1 Msun | | three points above | COMPUTED |
| MS, convective core above about 1.2 Msun; CNO overtakes pp near 1.3 Msun | mass thresholds | | not checked | RECALLED |
| Red giant branch, first dredge-up (Z = 0.02, initial Y = 0.28) | Y 0.2800 -> 0.3001 (0.85 Msun), 0.3025 (1.0), 0.2957 (1.5), 0.2888 (2.0), 0.2915 (2.5) | mass fraction | Boothroyd and Sackmann 1999 [astro-ph/9512121] Table 2 | READ |
| Same, C/Z and N/Z | C/Z 0.1733 -> 0.1530 (1.0 Msun), 0.1273 (1.5), 0.1167 (2.0), 0.1143 (2.5); N/Z 0.0531 -> 0.0771, 0.1073, 0.1197, 0.1324 | | same | READ |
| Same, 12C/13C | 90 -> 29.8 (1.0), 23.5 (1.5), 21.7 (2.0), 21.1 (2.5); O/Z unchanged to 2 Msun | number ratio | same | READ |
| Same, Z = 0.001 (Y 0.2400) | Y -> 0.2598 (1.0 Msun), 0.2606 (1.5), 0.2545 (2.0); C/Z 0.0765 -> 0.0654, 0.0518, 0.0449 | | same, Table 4 | READ |
| Extra mixing on the RGB (cool bottom processing) | lowers 12C/13C to 11-15 for 1-2 Msun at solar Z; strength rises roughly as M^-2 and Z^-1 | | same, abstract and Table 2 | READ |
| Horizontal branch, red clump | core He burning, core C/O grows from the 3-alpha reaction; surface Y after RGB; colour set by [Fe/H] and mass loss (already in `stellar_evolution.py`) | | project tables checked in `stellar-evolution-tables-check-2026-10-17.md`; core C/O values | RECALLED |
| AGB, carbon stars | Third dredge-up makes C/O > 1 above a minimum initial mass: 1.5-1.6 Msun (solar Z), 1.5 Msun at Z=0.02 and 1.4 Msun at Z=0.01 (Cristallo 2011) | Msun | Karakas and Lugaro 2016 [1604.02178] sect. 2 | 2ND-HAND |
| AGB, hot bottom burning | Base of envelope hot enough for proton capture in intermediate-mass AGB stars (grid 1-7.5/8 Msun, Z = 0.007, 0.014, 0.03); raises He and N; produces Li; threshold mass not extracted | | same, abstract | READ |
| AGB, s-process | Slow neutron capture enriches Sr-Ba-Pb peaks in the envelope, yields tabulated C to Bi for 1-8 Msun | | same; yield tables are not extracted here | READ |
| Wolf-Rayet WN (H-free early type) | X_H ~ 0, X_N ~ 1% by mass, X_C ~ 0.05% | mass fraction | Crowther 2007 [astro-ph/0610356] sect. 3.5 | READ |
| WN late type | X_H ~ 15 +- 10% | | same | READ |
| WN/C transition | X_C ~ 5%, X_N ~ 1% | | same | READ |
| WC | C/He = 0.04-0.7 by number (10-60% C by mass); refined 0.1-0.5 (20-55%) | | same | READ |
| WO | rare, oxygen-rich, the most advanced He-burning stage; strong O v-vi lines; mass fractions not extracted, so the code's carbon 30-55% and oxygen 10-40% are a stand-in (BC24) | | same, sect. 2.1 | READ |
| WR progenitor and WR masses | typical masses 10-25 Msun, up to 80 for H-rich WN; progenitors above ~25 (H-free WN), ~40 (WC), ~75 Msun (H-rich WN) | | same abstract | READ |
| Population III | initial mass 10-1000 Msun in a sample of 100 first stars | Msun | Hirano et al. 2014 [1308.4456] abstract | READ |
| Pop III composition | pure BBN: X = 1 - Y_P = 0.7533, Y = 0.2467, Z = 0 | | Planck eq. 72a | COMPUTED |
| Extremely metal-poor stars | C-rich fraction 20-30% for [Fe/H] -4 to -2 | | Frebel and Norris 2015 | READ |
| Core He fraction vs age | Y_c = 1 - X_c - Z_c; no table beyond the Sun | | derived from X_c above | COMPUTED |

### 2.11 White dwarfs

| Quantity | Value | Unit | Source | Label |
|---|---|---|---|---|
| IFMR currently in the project | M_f = (0.109 +- 0.007) M_i + 0.394 +- 0.025 | Msun | Kalirai et al. 2008 [0706.3894] line 895; matches `tuning.py` | READ |
| IFMR, PARSEC-based, piecewise | M_f = 0.0873 M_i + 0.476 (0.87-2.80); 0.181 M_i + 0.210 (2.80-3.65); 0.0835 M_i + 0.565 (3.65-8.20) | Msun | Cummings et al. 2018 [1809.01673] eqs. 1-3 | READ |
| IFMR, MIST-based, piecewise | 0.080 M_i + 0.489 (0.83-2.85); 0.187 M_i + 0.184 (2.85-3.60); 0.107 M_i + 0.471 (3.60-7.20) | Msun | same, eqs. 4-6 | READ |
| Difference between Kalirai and Cummings PARSEC at M_i = 1, 2, 3, 4, 6, 8 | M_f 0.503/0.563, 0.612/0.651, 0.721/0.753, 0.830/0.899, 1.048/1.066, 1.266/1.233 | Msun | both relations | COMPUTED |
| DA mass distribution | extremely narrow peak at 0.59 Msun, shoulder at 0.7-0.9, near absence of M > 1 among cool DAs | Msun | Kilic et al. 2020 [2006.00323] abstract | READ |
| Spectral types | about 80% of white dwarfs are DA | | Xu, Rogers, Blouin 2024 [2404.15425] sect. 2 | READ |
| Polluted fraction | 25-50% of white dwarfs show heavy elements | | same | READ |
| Core types | He core for very low mass; C/O core for typical; O/Ne core above about 1.05 Msun | Msun | Althaus et al. 2010 [1007.2659] text on core composition | READ |
| Ultra-massive WD | O/Ne cores from the super-AGB phase; crystallisation affects the majority below 40 000 K | | Camisassa et al. 2019 [1807.03894] abstract | READ |
| C/O split in the core | lower oxygen abundance is expected in more massive progenitors (Althaus); the numeric C:O split by mass was not extracted | | Althaus 2010 text | RECALLED |
| Mass-radius, zero temperature, mu_e = 2 | M = 0.223, 0.510, 0.897, 1.096, 1.273, 1.395 Msun at central x = 0.5, 1, 2, 3, 5, 10; R = 14049, 9707, 6477, 4998, 3513, 2071 km | | relativistic degenerate electron gas (Chandrasekhar EOS), CODATA 2018 constants, hydrostatic equilibrium integrated with SciPy; no Coulomb or temperature corrections | COMPUTED |
| Same, interpolated | R(0.5) = 0.01411, R(0.6) = 0.01271, R(1.0) = 0.00821, R(1.2) = 0.00599, R(1.3) = 0.00466 Rsun; limiting mass 5.83/mu_e^2 = 1.457 Msun | Rsun | same | COMPUTED |
| Sirius B test | dynamical mass 1.018 +- 0.011 Msun (Bond 2017), radius 0.803 +- 0.011 x 10^-2 Rsun, gravitational redshift mass 1.017 +- 0.025 | | Joyce et al. 2018 [1809.01240] | READ |
| Model for Sirius B | R(1.018) = 0.00802 Rsun (within 0.1% of the measurement; real stars have thin H/He layers and Coulomb corrections, so a 3-5% tolerance is the honest test limit) | Rsun | my integration | COMPUTED |
| Polluted WD geochemistry, six WDs | accreted rocks have oxygen fugacities like Earth, Mars and typical asteroids | | Doyle et al. 2019 [1910.12989] abstract | READ |
| Polluted WD geochemistry, 16 WDs with Al, Ca, Si, Mg, Fe, O all detected | most oxidising; about 25% consistent with more reduced parent bodies | | Doyle et al. 2020 as reported in Xu et al. 2024 sect. 4 | 2ND-HAND |
| Mass of heavy elements in a polluted WD convection zone (bulk-Earth composition) | 1e20 to 1e25 g; parent-body masses 1e23 to 7e25 g | g | Xu et al. 2024 (citing Harrison 2021) | 2ND-HAND |

### 2.12 Neutron stars

| Quantity | Value | Unit | Source | Label |
|---|---|---|---|---|
| Precisely measured masses | about 35 stars from 1.17 to 2.0 Msun; companion of J0453+1559 1.174 +- 0.004 | Msun | Ozel and Freire 2016 [1603.02698] sect. 2 | READ |
| Mass populations | double NS 1.33 +- 0.09; recycled 1.54 +- 0.23; slow pulsars 1.49 +- 0.19 (Gaussians); two-peak recycled fit 1.388 +- 0.058 and 1.814 +- 0.152 | Msun | same | READ |
| Heaviest | 2.35 +- 0.17 (J0952-0607); 2.14 +0.10 -0.09 (J0740+6620 first); 2.08 +- 0.07 (J0740 refined); 2.01 +- 0.04 (J0348+0432) | Msun | Romani 2022 [2207.05124]; Cromartie 2020 [1904.06759]; Fonseca 2021 [2104.00880]; Antoniadis 2013 [1304.6875] (as quoted in Cromartie) | READ |
| Maximum nonrotating mass | 2.01 <= M_TOV <= 2.16 +0.17 -0.15 (GW170817 argument); causality bound 4.09 Msun | Msun | Rezzolla et al. 2018 [1711.00314]; Lattimer and Prakash 2016 [1512.07820] eq. 10 | READ |
| Radii | R(1.4) = 12.45 +- 0.65 km and R(2.08) = 12.35 +- 0.75 km (full range of the 1-sigma intervals over three EOS frameworks); J0740+6620 R_eq = 13.7 +2.6 -1.5 km (68%) | km | Miller et al. 2021 [2105.06979] abstract | READ |
| J0030+0451 NICER | M = 1.34 +0.15 -0.16 Msun, R_eq = 12.71 +1.14 -1.19 km | | Riley et al. 2019 [1912.05702] abstract | READ |
| Minimum mass | no sharp minimum; forecast: with 50 merging binary NS the lower mass cutoff can be measured to about 0.2 Msun and the upper to 0.1 Msun | Msun | Chatziioannou and Farr 2020 [2005.00482] abstract | READ |
| Outer crust ground-state sequence (nucleus, density at which it stops, g/cm3) | 56Fe 8.02e6; 62Ni 2.71e8; 64Ni 1.33e9; 66Ni 1.50e9; 86Kr 3.09e9; 84Se 1.06e10; 82Ge 2.79e10; 80Zn 6.07e10; 82Zn 8.46e10; 128Pd 9.67e10; 126Ru 1.47e11; 124Mo 2.11e11; 122Zr 2.89e11; 120Sr 3.97e11; 118Kr 4.27e11 | | Chamel and Haensel 2008 [0812.3955] Table 2 (Ruster et al. data + BSk8) | READ |
| Neutron drip | about 4e11 g/cm3 (SLy), 6e11 (FPS) | | same | READ |
| Inner crust clusters (Negele and Vautherin) | Z = 40 (Zr) at 4.67e11-2.66e12; Z = 50 (Sn) 6.24e12-3.41e13; Zr 7.94e13; Ge (Z = 32) 1.32e14 g/cm3, cluster counts N 140 to 1460 | | same Table 3 | READ |
| Nuclear saturation density | 2.7-2.8e14 g/cm3, 0.16 fm^-3 | | Miller et al. 2021 [2105.06979] intro | READ |
| Crust thickness about 1 km, mass about 0.01 Msun | | | typical values, not checked | RECALLED |
| Core | n, p, e, mu in beta equilibrium; hyperons or quarks unconstrained; EOS family fixed by M-R | | Lattimer and Prakash 2016 | 2ND-HAND |

### 2.13 Black holes (no composition)

| Quantity | Value | Unit | Source | Label |
|---|---|---|---|---|
| Stationary black holes are described by mass, angular momentum and charge only (Kerr-Newman family) | no hair | | Chrusciel, Costa, Heusler 2012 [1205.6112] | READ |
| Schwarzschild radius r_s = 2GM/c^2 | 8.86 (3 Msun), 14.77 (5), 29.53 (10), 59.07 (20) | km | CODATA 2018 constants | COMPUTED |
| Hawking temperature T = hbar c^3 / (8 pi G M k) | 2.06e-8 (3 Msun), 1.23e-8 (5), 6.17e-9 (10), 3.09e-9 (20) | K | same | COMPUTED |
| Kerr horizon r_+ = M (1 + sqrt(1 - a^2)) and prograde ISCO | a* = 0, 0.5, 0.9, 0.998: r_+ = 2.000, 1.866, 1.436, 1.063; ISCO = 6.000, 4.233, 2.321, 1.237 (GM/c^2); radiative efficiency 0.057, 0.082, 0.156, 0.321 | | Bardeen formulae | COMPUTED |
| Spin cap 0.998 | in `BLACK_HOLE_SPIN_RANGE` | | Thorne 1974 | RECALLED |
| Astrophysical charge | effectively zero (store 0) | | standard argument | RECALLED |

### 2.14 Brown dwarfs and the lowest-mass stars

| Quantity | Value | Unit | Source | Label |
|---|---|---|---|---|
| Deuterium burning mass | about 13 MJ; 50% of D burns at 13.0 +- 0.8 MJ; full range 11.0 MJ (3x solar Z, 10% burned) to 16.3 MJ (Z = 0, 90% burned) | MJ | Spiegel, Burrows, Milsom 2011 [1008.5150] abstract | READ |
| For He 0.25-0.30, D/H 1-2e-5, 0.5-2 solar Z | narrower range, 50% criterion | | same | READ |
| Lithium burning mass | about 0.06 Msun (63 MJ) | Msun | Burrows et al. 2001 [astro-ph/0103383] sect. IV | READ |
| Hydrogen burning limit, solar metallicity | edge mass near 0.075 Msun; BHAC15 models start at 0.07 Msun as the H-burning limit | Msun | Burrows 2001; Baraffe et al. 2015 [1503.04107] | READ |
| Fully convective | all brown dwarfs and stars below 0.35 Msun | | Chabrier and Baraffe 1997 | READ |
| Radius (Rsun) at log t = 7 (10 Myr) | 0.157 (0.010 Msun), 0.245 (0.020), 0.284 (0.030), 0.292 (0.040), 0.315 (0.050), 0.337 (0.060), 0.359 (0.070), 0.383 (0.080), 0.416 (0.100) | Rsun | BHAC15 track file, interpolated | COMPUTED |
| Radius (Rsun) at 100 Myr | 0.125 (0.020), 0.124 (0.030), 0.130 (0.040), 0.139 (0.050), 0.148 (0.060), 0.158 (0.070), 0.167 (0.080), 0.186 (0.100) | Rsun | same | COMPUTED |
| Radius (Rsun) at 1 Gyr | 0.0960 (0.050), 0.0940 (0.060), 0.0944 (0.070), 0.0971 (0.075), 0.102 (0.080), 0.113 (0.090), 0.124 (0.100), 0.174 (0.150), 0.291 (0.300), 0.460 (0.500) | Rsun | same | COMPUTED |
| Radius at 10 Gyr | 0.0900 (0.070), 0.0910 (0.075), 0.0990 (0.080), 0.175 (0.150), 0.304 (0.300) | Rsun | same | COMPUTED |
| Teff at 1 Gyr | 1383 K (0.050 Msun), 1586 (0.060), 2051 (0.070), 2401 (0.080), 2810 (0.100), 3130 (0.150) | K | same | COMPUTED |
| Surface lithium, log(Li/Li0) at 1 Gyr | 0 to -0.03 (0.050 Msun), -4.68 (0.060), -10.9 (0.070), -10.7 (0.080); at 10 Gyr -11 (0.070) | | same | COMPUTED |
| Track coverage | tracks for M < 0.07 Msun end before 10 Gyr in the structure file (0.010 Msun ends at 41 Myr, 0.020 at 123 Myr, 0.050 at 1.5 Gyr) | | same | COMPUTED |
| Luminosity at 10 Gyr, log L/Lsun | -4.30 (0.070), -3.89 (0.075), -3.58 (0.080), -2.58 (0.150) | | same | COMPUTED |

## 3 How the generator should compute or draw it

### 3.1 Two routes

**Route 1 (keep drawn bulk figures, fill a composition that fits).** Keep `sample_fe_h` as it is. After [Fe/H] and population are known, draw the extra ratios `[Mg/Fe], [Si/Fe], [O/Fe], [C/Fe]` (and optionally `[Eu/Fe]`) from section 2.7, then compute X, Y, Z (2.4) and the planet-facing ratios (2.3). Nothing already stored moves; planets, orbits and habitability keep their current inputs.

**Route 2 (draw composition first, derive the rest).** Draw the star's birth gas composition from a Galactic chemical evolution model (age, R, z -> mean [Fe/H] and [alpha/Fe] from the enrichment history), then derive [Fe/H] from the composition. This is physically cleaner (age-metallicity-alpha become linked) but needs a chemical evolution model the project does not have, and it would change every existing star's [Fe/H].

**Recommendation: Route 1, with two cheap links from Route 2.** The measured distributions (2.6, 2.7) already encode the age-alpha-metallicity link in the form of the P(high-alpha | [Fe/H], age) membership, so a seeded draw from them is "published distribution, not made up". Link A: make the high-alpha probability depend on population and age (young and intermediate populations thin only; old, bulge, halo follow the logistic; the age table in 2.6 gives the shape). Link B: optionally set the mean [Fe/H] of the thin disk from age with the -0.027 dex/Gyr slope (scatter 0.265 dex), only if the project wants an age-metallicity relation; at present ages and metallicities are drawn independently, so metal-poor young stars can occur.

### 3.2 Draw order (per star)

1. Inputs: mass, age, population, galactic R and z, [Fe/H] (all stored).
2. Branch: `branch = thin` for young and intermediate populations. For old and bulge: `high` with probability `p = 1 / (1 + exp(([Fe/H] + 0.350) / 0.132))` times a population scale `f_pop` (stored target; start 1.0 for old, since the old population's own mean of -0.2 already encodes the thick/thin mix); halo: `high` with probability 1 (add a low-alpha halo branch when Mg, Si are fitted). Cluster and satellite stars: use the branch of the component they belong to (they already carry peaks).
3. Means: `mu_Mg, mu_Si, mu_O, mu_C (, mu_Eu)` from the 2.7 table for the branch at the clamped [Fe/H]. Residuals: draw four standard normals, correlate with Cholesky of the BF16 correlation matrix in 2.7 (C, O, Mg, Si), scale by the branch sd (thin: C 0.044, O 0.052, Mg 0.066, Si 0.042; high: 0.050, 0.065, 0.073, 0.073). Store results.
4. Other elements: Ca, Ti, S, Ar follow Ca (2.7), Ne follows Mg (the lower-Z option; switch to O when the project decides), N [N/Fe] = 0 (not fitted), Na, Al, Ni, Cr from 2.7 if ever needed, everything else [X/Fe] = 0.
5. Compute r = Z/X, X, Y, Z (2.4). Check X + Y + Z = 1.
6. Ratios: log eps_i = log eps_sun,i + [Fe/H] + [X_i/Fe]; Fe/Si, Mg/Si, C/O and O/Si as number and mass ratios.
7. Radiogenic birth abundances: Eu from [Eu/Fe] (2.7), Th = Eu x (Th/Eu)_birth,sun x 10^[Th/Eu] with [Th/Eu] ~ N(0.014, 0.045); U/Th from the Sun's birth ratio 0.572 (scatter not measured); K from [K/Fe] = 0 plus the alpha-like trend (not fitted, gap). A system of age tau then has decayed fractions from the table in 2.9.
8. CEMP branch: for [Fe/H] < -2, with probability 0.25 (-4 < [Fe/H] < -2) rising to 1 below -4, add [C/Fe] ~ U(0.7, 2.0) (the C-rich definition is [C/Fe] > +0.7; the upper bound is a parametrisation). Metal-poor low-mass stars only.

### 3.3 Quantities: stored versus new

| Quantity | Status | Unit | Range | Per |
|---|---|---|---|---|
| mass, age, [Fe/H], population, R, z, spectral type, class | already stored | | | star |
| `abundance_branch` (thin / high / halo_low / cemp) | new, small enum | | 4 values | star |
| `mg_fe`, `si_fe`, `o_fe`, `c_fe` | new, float | dex | -0.3 to +0.7 (C up to +2.5 for CEMP) | star |
| `eu_fe` (only if the radiogenic lane needs Th/U) | new, float | dex | -0.2 to +0.6 | star |
| X, Y, Z, Fe/Si, Mg/Si, C/O, O/Si, Th, U, K birth abundances | derived on demand, not stored (cache allowed; invariant below) | mass fraction or ratio | | star |
| core H fraction X_c(t) | derived on demand | mass fraction | 0 to X | star |
| WD core type, NS crust zone, remnant radius | derived on demand from stored mass | | | remnant |

### 3.4 Stellar phases: what to compute

- **Main sequence.** Surface = birth composition (X, Y, Z). Interior: fully convective and homogeneous below 0.35 Msun; above that the core burns H. Central hydrogen: `X_c = X (1 - (t / t_x)^1.15)` with `t_x` the hydrogen-exhaustion time (use the project's `main_sequence_lifetime_gy` scaled so that 1 Msun gives 8 Gyr, or keep 10 Gyr and re-fit the exponent; this is a decision item because the project's 10 Gyr and the Sun's 8 Gyr differ by 25%). Core helium `Y_c = 1 - X_c - Z`. Valid for roughly solar-mass stars; no mass-dependence is read.
- **Red giant branch.** After first dredge-up the surface Y rises by +0.01 to +0.023 (0.85 to 1.5 Msun), C/Z falls by 12-34%, N/Z rises by 45-145%, 12C/13C falls from 90 to 21-30 (table 2.10). Parametrise by interpolating the table in mass and scaling with Z (the Z = 0.001 table gives the low-Z behaviour). Store nothing: apply as an additive shift to Y, C, N when the evolved-star class is III.
- **Horizontal branch / red clump.** Surface equals post-RGB composition; core is He-burning; values of core C/O are not read (RECALLED, flag).
- **AGB and carbon stars.** If M exceeds the carbon-star minimum (1.4-1.6 Msun at solar Z, lower at low Z) and the star is on the thermally pulsing AGB, set C/O > 1 (carbon star). Hot bottom burning in intermediate-mass AGB stars raises He and N and produces Li. Amounts from Karakas and Lugaro 2016 tables are not extracted here (work item).
- **Wolf-Rayet.** Use section 2.10 mass fractions: early WN H-free, X_N ~ 1%; late WN X_H ~ 15%; WN/C X_C ~ 5%; WC C 10-60% by mass; draw within the quoted ranges.
- **Population III.** X = 1 - Y_P, Y = Y_P, Z = 0 exactly; mass 10-1000 Msun; no planets.
- **White dwarfs.** Final mass from the IFMR (keep Kalirai, which the project already uses and which is READ; Cummings PARSEC is 0.03 to 0.07 Msun higher below 4 Msun, see 2.11). Core type from mass: He for very low mass remnants (threshold not read; about 0.45 Msun is RECALLED), CO for the typical case, ONe above about 1.05 Msun. Atmosphere type DA with probability 0.8, else DB/DC/DZ (atmosphere detail is another lane). Radius: interpolate the 11-point zero-temperature Chandrasekhar table in 2.11 in log mass (or fit it once and store the fit); it reproduces Sirius B to 0.1%. Pollution: 25-50% of old WDs carry accreted rock; the planet lanes supply the accreted composition.
- **Neutron stars.** Mass from the Gaussian mixture in 2.12 (recycled 1.54 +- 0.23, slow 1.49 +- 0.19), floored at the measured 1.17 and capped near the maximum mass (2.0 to 2.4 Msun). Radius: the NICER and EOS full range is 11.6 to 13.1 km at 1.4 and 2.08 Msun, nearly mass independent, so draw radius from that band rather than the project's 10-13 km (the 10-11.6 km part has no support in the cited constraints). Layers by density using the Chamel-Haensel sequence (outer crust 56Fe to 118Kr, neutron drip near 4e11 g/cm3, inner crust Zr/Sn clusters in a neutron gas, uniform core).
- **Black holes.** Store mass, spin, charge (0). r_s, Hawking temperature, horizon radius, ISCO and efficiency from the formulas in 2.13. The stellar BH mass range is the project's 5-20 Msun; nothing here constrains it.
- **Brown dwarfs.** Mass classes by burning limit: below ~13 MJ no D burning, 13 MJ to 0.06 Msun D burning only, 0.06-0.075 Msun lithium burned, above 0.075 Msun hydrogen burning (solar Z; the limit rises with lower Z, value not read). Radius from the BHAC15 table at the object's age (interpolate in log mass, log age); fully convective; composition is the birth composition with Li depleted past the Li-burning mass.

### 3.5 Consistency invariants a test could assert

1. `abs(X + Y + Z - 1) < 1e-9` for every star.
2. At `[Fe/H] = 0` and all `[X/Fe] = 0`, the derived (X, Y, Z) equal the reference set's proto-solar values to 1e-6 (self-consistency of the scale, not a preset).
3. `Y >= Y_P` for every star and `dY/dZ` between stars equals `g` (1.68) within 1e-9 for the birth composition.
4. Recomputing Fe/Si, Mg/Si, C/O from the stored `[X/Fe]` equals the cached ratios.
5. Over 1e5 draws of the thin branch, the sample sd of the residuals equals the table sd within 5%, and the Mg-Si correlation within 0.03 of 0.77.
6. WD radius within 5% of the Chandrasekhar integration at the same mass; NS mass in [1.1, 2.4] Msun; BH `r_s = 2GM/c^2` to 1e-6.
7. Pop III star has Z = 0 and Y = Y_P.

### 3.6 What would change in the code

- `physics/metallicity.py`: add a `sample_abundance_ratios` function taking [Fe/H], population, age and the RNG; keep `sample_fe_h` unchanged. The disc gradient 0.06 dex/kpc is 0.01 above the 0.053 dex/kpc COMPUTED from H15 (R > 6 kpc); update the "recalled" comment with the source and change the constant to 0.053 once the decision is made.
- `tuning.py`: new tables holding the hinge parameters, sds and Cholesky matrix; solar-set constants (A21 log eps, atomic weights, Y_P, g); remnant tables.
- `generation/star.py`: new fields (3.3); the `solids_factor` rule `10^[Fe/H]` can stay as the disk-solids law but should be re-based on `Z/Z0` once Z is available (see 2.5: the two differ by up to a factor 3 at low [Fe/H]).
- Database: one migration adding five floats and one enum to stars and bright_stars.

## 4 Calibration and test targets

Solar System or measured bodies used only to test the law, never as defaults.

| Target | Expected | Tolerance | Source | Label |
|---|---|---|---|---|
| Reference scale reproduces itself | proto-solar (X, Y, Z) = (0.7121, 0.2725, 0.0154) at [Fe/H] = 0, [X/Fe] = 0 | 1e-6 | A21 | COMPUTED |
| Local FGK dwarfs, ratio medians | Fe/Si 1.81, Mg/Si 0.93, C/O 0.385 (mass) | within 0.2, 0.1, 0.05 | BF16 | COMPUTED |
| Local FGK dwarfs, 95% spread | Fe/Si 1.17-2.39, Mg/Si 0.79-1.12, C/O 0.24-0.47 | factor 1.3 on width | BF16 | COMPUTED |
| Thin-disk Mg/Si residual | 0.041 dex | 20% | B14 | COMPUTED |
| Midplane gradient | -0.053 dex/kpc for R 6-14 kpc | 0.01 | H15 | COMPUTED |
| Midplane sigma[Fe/H] | 0.18-0.24 dex | 0.03 | H15 | READ |
| Sun's central H fraction | 0.340 at 4.57 Gyr, 1 Msun | 0.03 | no source | READ |
| First dredge-up | 1 Msun, Z = 0.02: Y +0.0225, C/Z x0.88, N/Z x1.45, 12C/13C 29.8 | 10% | Boothroyd and Sackmann | READ |
| Carbon star minimum mass | 1.4-1.6 Msun at Z = 0.01-0.02 | 0.2 Msun | Karakas and Lugaro | 2ND-HAND |
| Sirius B | M = 1.018 Msun gives R = 0.0080 Rsun | 5% | Joyce 2018 | READ |
| WD mass peak | 0.59 Msun for cool DA | 0.05 | Kilic 2020 | READ |
| Heaviest NS | 2.35 +- 0.17 Msun must be reachable by the draw tail | | Romani 2022 | READ |
| J0030 radius | 12.7 km at 1.34 Msun lies inside the model band | 1 km | Riley 2019 | READ |
| BH r_s | 29.5 km at 10 Msun | 1e-3 | COMPUTED | COMPUTED |
| Brown dwarf radius | 0.0944 Rsun at 0.07 Msun, 1 Gyr | 5% | BHAC15 | COMPUTED |
| Th/Eu constancy | [Th/Eu] flat in age, 0.014 +- 0.045 | | Botelho 2019 | READ |

## 5 Gaps

### 5.1 Boss must research (still unreachable)

1. Lodders 2003, ApJ 591, 1220, Table 1 (log eps and CI): the IOP page returned a bot-check page and was not routed around. It gives the exact L03 set the GEN.425 report quotes (Fe/Si 1.69 mass). L09 Table 4 here supersedes it, so this is low priority. Label NOT REACHED.
2. Grevesse and Sauval 1998, Space Sci. Rev. 85, 161: only its X, Y, Z via L09 (2ND-HAND). Needed only if the project wants GS98 as a selectable set. NOT REACHED.

### 5.2 To verify (RECALLED items, ranked by impact)

1. The O plateau and Ne choice that drive Z at low [Fe/H] (2.5): check [O/Fe] in the halo/thick disk from forbidden-line (6300 A) and NLTE triplet work (N14 supplies both columns; I used the NLTE triplet) and a measured [Ne/Fe] trend. Impact: factor 1.4 in Z at [Fe/H] <= -1.
2. Salaris et al. 1993 alpha-enhancement factor (RECALLED), used only as a comparison.
3. Mass thresholds: convective core above ~1.2 Msun; CNO over pp near 1.3 Msun. Solar Fusion II [1004.2318] is downloaded but I did not extract a crossover value.
4. Half-lives of 238U, 235U, 40K and the isotopic fractions 235U 0.72%, 40K 0.0117% (232Th 14.05 Gyr is READ in Botelho).
5. Horizontal-branch core C/O and CO white-dwarf core C/O split (Althaus 2010 text; table values not extracted).
6. Black hole spin cap 0.998 (Thorne 1974), zero charge, and NS crust thickness (1 km) and mass (0.01 Msun).
7. K abundance trend with [Fe/H] (needed for 40K heating): not fitted.

### 5.3 Reachable but not yet read or extracted (agent time-box)

- Karakas and Lugaro 2016 yield and surface-abundance tables (C/O and [s/Fe] vs thermal pulse) [1604.02178]; Lagarde et al. 2012 grid [1204.5193] for thermohaline/rotation effects.
- Low-alpha halo Mg, Si (Nissen and Schuster 2010, not located by arXiv id); Eu in the halo.
- Hydrogen-burning limit vs metallicity (Chabrier and Baraffe 1997 has it; not extracted).
- NS mass-radius relation functions; Lattimer and Prakash [1512.07820] tables of EOS families.
- Mass-dependent X_c(t) for stars other than the Sun (needs a track grid with central abundances; BHAC15 files have Mrad and Rrad but no X_c).
- Hayden et al. [1503.02110] alpha distribution functions (Table 2) for the alpha side of the gradient.

## 6 Proposed work items

| Title | Scope | Depends on | Size |
|---|---|---|---|
| A-R1 Decide the reference scale | Choose A21 proto-solar log eps as the scale; record the A21/L09/L19 spread as the uncertainty on ratios | this note | small |
| A-R2 Decide Ne and O treatment | Settle Ne following O or Mg, and the O plateau (NLTE triplet vs forbidden line) | A-R1 | small |
| A-T1 Solar-set table module | `tuning.py` constants: 83-element log eps, atomic weights, Y_P, dY/dZ, CI-chondrite scale | A-R1 | small |
| A-T2 [X/Fe] model table | Hinge parameters, sds, Cholesky matrix from 2.7, with provenance comments | A-R2 | small |
| A-T3 Remnant tables | WD IFMR choice, NS crust sequence, NS mass mixture, BD burning limits and BHAC15 radius grid | none | medium |
| A-D1 Branch and ratio draw | `sample_abundance_ratios` (branch, means, correlated residuals, CEMP) seeded through `draw.py` | A-T1, A-T2 | medium |
| A-D2 X, Y, Z derivation | Closed form 2.4 from [Fe/H] and ratios | A-T1, A-D1 | small |
| A-D3 Age link (optional) | Thin-disk [Fe/H]-age relation and age-gated membership | A-D1 | medium |
| A-D4 Evolved-star surface shifts | FDU table interpolation, carbon-star flag, WR compositions, Pop III | A-D2 | medium |
| A-D5 Radiogenic birth abundances | Eu, Th, U, K with decay law for the radiation lane | A-D1, A-R3 | medium |
| A-D6 Remnant structure | WD mass-radius, core type, NS layers, BH quantities, BD radius | A-T3 | medium |
| A-S1 Storage | Schema migration: branch enum and four to five floats on stars and bright_stars | A-D1 | small |
| A-S2 Rebase solids law | Replace `10^[Fe/H]` by `Z/Z0` in disk surface density where wanted | A-D2 | small |
| A-V1 Tests | Invariants 3.5 and targets section 4 as unit tests | A-D1 to A-D6 | medium |
| A-V2 Distribution test | 1e5 draws compared with BF16 quantiles and Thiabaud box shares | A-D1 | small |
| A-X1 Docs | Design section from this part; correct the "recalled" comments in `tuning.py` | all | small |
| A-R3 Extract remaining tables | Items in 5.3 | none | medium |

## 7 Verification summary

Counts of table rows by label (a row with two labels counts once under each; counted from the tables in section 2, with the section 4 targets counted in their own column).

| Label | Section 2 rows | Section 4 targets | Meaning for this part |
|---|---|---|---|
| READ | 92 | 8 | paper table, equation or abstract read in the arXiv PDF or CDS table (includes 19 element rows and 2 solar-set rows counted also under COMPUTED) |
| 2ND-HAND | 9 | 1 | GS98 and A05 sets through L09, Jimenez via Casagrande, Thiabaud envelope via the GEN.425 report, Doyle 2020 via Xu 2024, WD pollution masses, NS core, carbon-star minimum mass |
| COMPUTED | 103 | 7 | re-derived with Python from READ inputs (mass fractions, ratios, fits to CDS catalogues, simulations, Chandrasekhar integration, BH formulas, BHAC15 interpolation) |
| SNIPPET | 0 | 0 | no search-result-only numbers |
| RECALLED | 9 | 0 | Salaris factor, half-lives, 40K fraction, MS mass thresholds, HB core C/O, WD core C:O split, BH spin cap, BH charge, NS crust thickness (all listed in 5.2) |
| NOT REACHED | 1 row (Lodders 2003) plus GS98 primary named in 5.1 | 0 | IOP page gave a bot-check page; GS98 primary not opened |

Totals: 193 labelled rows in section 2 and 16 in section 4. Grounded labels (READ + 2ND-HAND + COMPUTED) are 204 label uses against 10 for RECALLED + NOT REACHED in section 2. Every RECALLED item is either a comparison number (Salaris), a standard constant (half-lives), or a threshold the generator can leave as a stored tunable until checked.

Cross-check with the GEN.425 sources check (`/mnt/project-files/reports/bulk-composition-sources-check-2026-10-17.md`): agree that Thiabaud's ranges are mass ratios and that Lodders 2003 gives Fe/Si 1.69 (L09 1.692); differ in that this part shows the Thiabaud box is a compilation envelope (Hk14 spread 0.11-0.14 dex) wider than the homogeneous-catalogue spread of BF16, and that metal-poor stars fall outside its Fe/Si lower edge.
