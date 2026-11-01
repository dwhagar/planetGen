# 10D Verification of Boss's second research note, part D (Rayleigh table, Parmentier and Guillot equation)

Checked by the research lane itself, between about 2:15 and 2:35 AM EDT as the lane clock reports. The note arrived at 2:11 AM EDT. Parts 10A (ices, redox, lava, CH4/NH3), 10B (energetic particles, kinetics) and 10C (aerosols, comets, isotopes, licences) are separate files. Labels and verdict words as in 08A. Scratch code: `/tmp/claude-0/pk/ray.py` (Rayleigh) and the Parmentier source `Parmentier2013b.tex` ([R3439] LaTeX source).

## 1. Rayleigh section of the note

Method: every dispersion formula in the note was evaluated at 532 nm with sigma = 24 pi^3 (n^2-1)^2 / (lambda^4 N^2 (n^2+2)^2) F_K and compared with the formulas already verified in the first report (section 2 of `reports/atmosphere-sources-check-2026-10-17.md`, 2ND-HAND from colour-science, eradiate, exo_k, [R0623, R2384]) and with the measured values of He et al. 2021 (ACP 21, 14927, abstract read through Crossref; its Table 2 values were read in the first report).

| # | Claim in the note | Verdict | Evidence | Label |
|---|---|---|---|---|
| D1 | sigma(nu) = 24 pi^3 nu^4 / N_s^2 ((n^2-1)/(n^2+2))^2 F_k; N_s = 2.54743e19 cm^-3 at 15 C and 101.325 kPa | CONFIRMED (and fixes the first note's 2.68678e19) | p/(k T) at 288.15 K gives 2.54692e19; the note's value is 0.02 percent higher (it is the figure printed in Sneep and Ubachs). The first note used the 273.15 K value for every gas | COMPUTED |
| D2 | F_k = (6 + 3 rho)/(6 - 7 rho) | CONFIRMED | Same law as 08A section 4.1 | READ |
| D3 | "Generic constant depolarization factor introduces +-5 to +-9 percent" | PARTLY | The 5 to 9 percent size is reproduced by the reference-density bias (5.5 percent, 08A 4.2) and by the generic polarizability formula (-9 to +6 percent, part 02 section 3.1), so the number is real but its stated cause is one of two | COMPUTED |
| D4 | N2: (n-1)e8 = 6498.2 + 3072900/(144 - lambda^-2); F = 1.034 + 3.17e-4 lambda^-2; sigma(532 nm) = 5.10e-27 | PARTLY | The dispersion is the verified one (3074330 verified against 3072900: 0.05 percent, n-1 = 2.8375e-4 against 2.8385e-4). The King factor is the verified one. The cross section computed from them is **5.31e-27**, so the printed 5.10e-27 is 4 percent low | COMPUTED |
| D5 | CO2: (n-1)e8 = 57918.17 + 167909/(57.362 - lambda^-2); F about 1.15; sigma = 1.24e-26 | **CONTRADICTED** (formula) | The formula gives n-1 = 6.104e-4 at 532 nm, 43 percent above the verified 4.278e-4 (Bideau-Mehu form, part 02), and sigma = 2.7e-26, twice the measured value. The F value 1.15 is right (1.1364 + 2.53e-3/lambda^2 = 1.145). The printed sigma 1.24e-26 is 7 percent below He et al. 2021's measured 13.3e-27, the same snippet value that the first report closed in favour of 13.3e-27 | COMPUTED, READ |
| D6 | O2: (n-1)e8 = 20563 + 2480990/(132.274 - lambda^-2); F = 1.096 + 1.38e-3 lambda^-2; sigma = 4.41e-27 | **CONTRADICTED** (formula) | The verified form is 20564.8 + 248090/(40.9 - lambda^-2) at 273.15 K (above 546 nm 21351.1 + 218567/(40.9 - lambda^-2)). The note's denominator constant is not that of any form in the sources read: n-1 = 3.98e-4 against the true 2.72e-4 (46 percent too high), sigma 1.1e-26 (2.4 times too high). F lacks the lambda^-4 term (1.448e-4). Printed sigma 4.41e-27 is 3 percent below He 2021's measured 4.55e-27 (n-based 4.67e-27) | COMPUTED, READ |
| D7 | Ar: (n-1)e8 = 6432.8 + 2860600/(144 - lambda^-2); F = 1; sigma = 4.68e-27 | PARTLY | Dispersion correct (Peck and Fisher form, n-1 = 2.6798e-4, identical to the verified one). The formula gives sigma = 4.57e-27; the printed 4.68e-27 is 2.4 percent above it | COMPUTED |
| D8 | CH4: (n-1)e8 = 46662 + 4.02e8/(111 - lambda^-2); F = 1; sigma = 1.15e-26 (He et al. 2021) | **CONTRADICTED** (formula), CONFIRMED (cross section) | The printed dispersion gives n-1 = 3.79e-2 (90 times too large; it is a garbled version of the older coded form n-1 = 4.6662e-4 + 4.02e-6/lambda^2, which part 02 already shows to be 30 percent high). He 2021's own dispersion (4869.8 + 4.1023e14/(1.133e10 - nu^2) at 288.15 K) gives 1.136e-26 at 532 nm, and the measured 11.3e-27; the printed 1.15e-26 is 1.8 percent above it | COMPUTED, READ |
| D9 | F_k(CH4) = 1.000 "spherical symmetry" | PARTLY | CH4 is a spherical top with zero intrinsic anisotropy, so F_K = 1 is the physical limit, and He 2021 uses it; a measured depolarization of order 1e-3 would change sigma by under 0.2 percent | READ |
| D10 | SO2: rho about 0.019, F = 1.05, sigma(532) = 3.12e-26, source [R0443] / [R4138] | **CONTRADICTED** (arithmetic) and NOT REACHED (data) | F(rho = 0.019) = 1.032, not 1.05 (F = 1.05 needs rho = 0.030). Part 02 DFT gives alpha = 3.88 A^3, rho = 0.039 and sigma(550) = 2.41e-26 (about 2.6e-26 at 532), 17 percent below the note's 3.12e-26. [R0443] was not reached. The note gives no dispersion formula ("via Bogaard") | COMPUTED |
| D11 | H2S: rho about 0.0015, F = 1.004, sigma(532) = 1.85e-26 | **CONTRADICTED** (arithmetic) and NOT REACHED (data) | F(rho = 0.0015) = 1.0025, not 1.004. Part 02 DFT gives sigma(550) = 1.98e-26 (about 2.2e-26 at 532), so 1.85e-26 is 16 percent low and gives H2S a smaller cross section than the polarizability suggests. [R0443] not reached | COMPUTED |
| D12 | "CRDS and BBCES validate the refractive-index-based cross sections within +-0.4 percent to +-2.8 percent over 307 to 725 nm" | PARTLY | He 2021 abstract (READ) is BBCES only (not CRDS) and 307 to 725 nm is right; the numbers 0.4 to 2.8 percent are the *differences* from n-based values (CO2 0.4 +- 1.2, N2O -0.6 +- 1.1, SF6 0.9 +- 1.4, O2 2.8 +- 1.2, CH4 0.9 +- 2.2 percent), not uncertainties | READ |
| D13 | The measurements "refute historical claims of King factors above 1.20 ... anomalous cross sections resulted from non-physical spherical grain extrapolations" | UNSUPPORTED | The He 2021 abstract and Table 2 say nothing about King factors above 1.20 or grain extrapolation. For simple gases F_K above 1.20 would need rho above 0.105, which no measurement of N2, O2, CO2, Ar, CH4 shows; that part of the statement is consistent with 08A 4.4 | COMPUTED |

Net for the generator (GEN.410, GEN.417): no coefficient from the second note replaces a verified one. Use the first-report formulas and the part 02 table. The note's N_s is right and removes the objection raised in 08A 4.2.

## 2. Parmentier and Guillot equation (note section 3)

The note now writes T^4(tau) = (3 T_int^4/4)(2/3 + tau) + (3 T_irr^4/4) mu_* [2/3 + 1/(gamma sqrt 3) + (gamma/sqrt 3 - 1/(gamma sqrt 3)) exp(-gamma tau sqrt 3)], gamma = kappa_v/kappa_th, and cites [R3439] (2015).

| # | Claim | Verdict | Evidence | Label |
|---|---|---|---|---|
| D14 | The equation is the double-grey profile | PARTLY | It is the semi-grey (double-grey) solution, equal to equation 30 of [R3439] (LaTeX line 505: `(3 T_irr^4/4) mu_* [2/3 + mu_*/gamma_v + (gamma_v/(3 mu_*) - mu_*/gamma_v) exp(-gamma_v tau/mu_*)]`) only when mu_* = 1/sqrt 3 inside the bracket. That is the [R1715] form with f replaced by mu_*. Written with a free mu_* in the prefactor and sqrt 3 inside the bracket the note's formula is self-inconsistent unless mu_* = 0.577; use eq. 30 with mu_* throughout | READ |
| D15 | It is the [R3439] (2015) model | PARTLY | The 2014 paper's own new result is the picket-fence model (eq. 98, coefficients A to E, read and coded in 08A section 3.1); the 2015 paper (Parmentier et al., [R3440]) supplies the opacity-based gamma and tau_lim tables. The semi-grey equation is the limit of eq. 98 | READ |
| D16 | "In steam atmospheres a two-band distribution with a line-to-continuum opacity ratio beta; beta < 0.1 atomic line cores, 0.1 < beta < 0.9 molecular bands" | **CONTRADICTED** | beta is defined in the paper (LaTeX line 324) as the equivalent bandwidth beta = delta nu/Delta nu, the fractional spectral width of band 1, not an opacity ratio (the ratio is kappa_1/kappa_2). The paper's examples are beta = 0.01 and 0.5. No division into "atomic" below 0.1 and "molecular" 0.1 to 0.9 exists in the text read | READ |
| D17 | T_eff relation | CONFIRMED (not in the note) | The paper's effective temperature is T_eff,mu^4 = T_int^4 + mu_* T_irr^4 (LaTeX line 240) | READ |

Net for the generator (GEN.412, GEN.423, giants): code eq. 30 with mu_* for the double-grey case and eq. 98 for the picket-fence case; define beta as the band fraction.

## 3. Counts for parts D

Tally of the 16 claims D1 to D16 (D17 is an extra check, not a claim of the note): CONFIRMED 2 (D1, D2), PARTLY 7 (D3, D4, D7, D9, D12, D14, D15), CONTRADICTED 6 (D5, D6, D8, D10, D11, D16), UNSUPPORTED 1 (D13).

## 4. Not reached (Boss must research)

1. [R0443] (Proc. R. Soc. A 363, 155): dispersion and depolarization of SO2, H2S and others.
2. [R1715] full text: only quoted through [R3439].

## 5. Recalled / to verify

Polarizability of SO2 (3.88 A^3) and H2S (3.63 to 3.68 A^3) come from the part 02 DFT calculation, which is COMPUTED and unverified against measurement. No other recalled values.
