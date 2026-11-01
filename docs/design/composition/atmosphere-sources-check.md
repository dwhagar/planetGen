# Atmosphere colour chain: sources check (GEN.408)

**Status (2026-10-17):** the result of GEN.408 (done), copied from the read-only thread's report so that every lane has it in the repository. It feeds [atmosphere-colour.md](atmosphere-colour.md) and GEN.409 to GEN.419; the lane that builds each item reads its section first. Where this report says a source was not reached, the input is open and is Boss's to research (section 7). Times in it are as the lane's clock reported.

Date: 2026-10-17 (task started 11:35 PM EDT on 10-16 as the session clock reported; lane clocks drift). Read-only research, no repo change. Chain: goal issue #1944, plan note `docs/design/composition/atmosphere-colour.md`, base `docs/design/composition/sky-colour-analysis.md`, traces `sky-color-analysis-2026-10-17.md` and `composition-density-link-2026-10-17.md`.

## 0. How far the check could go (read this first)

The session's proxy lets through only a short list of hosts. **Reachable:** `export.arxiv.org` (abstracts **and full PDFs**: `curl -o x.pdf https://export.arxiv.org/pdf/<id>`), `ar5iv.labs.arxiv.org` (some papers only; it redirects others to arxiv.org), `raw.githubusercontent.com` (source files of open-source code), `pypi.org` and `files.pythonhosted.org` (package sources, including data files). **Refused (403):** `arxiv.org` itself, doi.org, ADS, NASA NTRS, Wikipedia, readthedocs, github.com, Optica/OSA, Science, ScienceDirect, Copernicus (AMT/ACP), PMC, NIST, ResearchGate, Zenodo, Semantic Scholar, Crossref, the MPI Mainz UV/VIS atlas and the Bremen ozone page. Web search still returns snippets.

So **no journal-only paper could be read** (Sneep and Ubachs, Bates, Peck and Huang, Penndorf, Bucholtz, Kasting and Donahue, [R2645], [R4392, R4394], [R2258]). Where a number is given for one of them it was read in a reachable arXiv paper that quotes it, or in the source of an open-source code that implements it with the citation. Those are marked **2nd-hand**. Nothing below is quoted from memory. A snippet from a search result is marked **snippet** and is not verified.

Labels used:

| Label | Meaning |
|---|---|
| READ | The paper's own text was read (arXiv PDF or ar5iv) |
| 2ND-HAND | The number or formula was read in a reachable paper or code that cites the original; the original was not reached |
| COMPUTED | Re-derived here with Python from READ or 2ND-HAND inputs (inputs stated) |
| SNIPPET | Seen only in a search result summary; not verified |
| NOT REACHED | No reachable page gave it |

Tip for later lanes: the arXiv PDFs and package sources are the way to read a paper here; `export.arxiv.org` is the host to use.

## 1. Verdict and counts

19 sources were on the list (the 13 named in the brief plus 6 recalled ones in the sky analysis report).

| Result | Count | Sources |
|---|---|---|
| READ and verified | 3 | Betremieux and Kaltenegger [R0375, R4250], [R3014, R3016] |
| Verified 2nd-hand (value or formula read in citing paper or code) | 10 | [R0282, R0575, R3480, R4138] (test value also COMPUTED), [R2279, R2645], [R4392, R4394], [R2258, R2274] |
| Original not reached, no substitute for what the brief needed | 4 | [R3470] (H2 dispersion; a 2nd-hand cross-section law replaces it), [R0599, R4023] (data and licence) |
| Citation verified only, content not needed | 2 | [R0560, R1940] |

**Corrections to what the notes said (9):**

1. Betremieux and Kaltenegger [R0375] is a 2015 paper (the notes say 2014). The Table 1 values are right.
2. [R4138] eq. 13 (CO2 refractive index) has a typo: last numerator 0.1218145e-6, not e-4 ([R2384] says so). Negligible in the visible.
3. The widely coded Sneep and Ubachs formula for **CO** gives a cross section about half the measured one; do not use it (section 3).
4. [R0599] covers the **Chappuis band only, 410-760 nm** (snippet), not 220-700 nm. The Hartley and Huggins bands need another source.
5. [R4023] has **eleven temperatures between about 193 and 293 K** (snippet), not "11 K". It is the Part 2 (temperature dependence) paper.
6. Sudarsky class bounds: Class III is T_eff above about 350 K (not 250 K), Class IV is 900-1500 K (not 900-1400 K).
7. The prior claims for H2 and He relative to N2 are slightly high: H2 is 0.20 (not 0.22), He is 0.013 (not 0.014). CO2 at 2.5 is confirmed.
8. [R3014, R3016] is the **colour of the planet seen from space** (reflected disk), not the sky colour from the ground. It supports the method only in that Rayleigh scattering plus the star's spectrum decide a clear world's look.
9. Venus total cloud optical depth is about 20-30 at 0.63 micron, not 30-40.

**One finding that changes the design:** ozone column does not follow a simple power of O2. The Chapman steady state gives [O3] proportional to [O2] times a square root locally, but the *column* flattens or peaks with O2 for hot stars. Use the [R2466] fit (section 4), which was checked, instead of a power law.

## 2. Rayleigh scattering

### 2.1 The law (verified 2nd-hand, four independent readings)

    sigma(lambda) = 24 pi^3 (n^2 - 1)^2 / ( lambda^4 N^2 (n^2 + 2)^2 ) * F_K

with F_K = (6 + 3 d) / (6 - 7 d), d the depolarization ratio ([R2384] eq. 6-7, [R0623] eq. 12, which attributes it to Sneep and Ubachs). N is the number density at which (n - 1) is defined.

**[R4138]**: JQSRT 92, 293 (volume and page 2nd-hand from [R2384]; the [R4138] is from the same, a search snippet gave a different one, so confirm the DOI before printing it). Original NOT REACHED. Verdict: formula and citation verified 2nd-hand; CO2 typo corrected.

**The reference-density hazard (found while checking the codes).** Each coded refractivity formula is defined at a particular temperature, so N must match:

- N2, Ar, CO2, CO, CH4: coefficients correspond to 288.15 K, N = 2.546899e19 cm^-3.
- O2: coefficients correspond to 273.15 K, N = 2.68678e19 cm^-3.
- Evidence: a mole-fraction mix built this way reproduces Peck-Reeder air at 550 nm to 0.05 percent; using the wrong N for O2 is 1.1 percent off. Betremieux and Kaltenegger's 1-micron table values equal the 288.15 K formula times 288.15/273.15 for N2, Ar and CO2 and are unscaled for O2. TauREx and exo_k mix conventions and are off by about 8 percent in sigma; do not copy them blindly.
- Our table must store, per gas, the reference temperature with the formula.

### 2.2 Per-gas formulas (n - 1 as a function of wavenumber s in cm^-1; (n-1)e8 = A + B/(C - s^2) unless noted)

All 2ND-HAND (read in colour-science, eradiate, exo_k, TauREx, petitRADTRANS, pyratbay sources and in [R0623] Appendix D, [R2384]); originals NOT REACHED. Sigma at 550 nm with the N above, from COMPUTED values:

| Gas | Formula and King factor | sigma(550 nm), cm2 | ratio to N2 | Status |
|---|---|---|---|---|
| N2 | (n-1)e8 = 6498.2 + 3.0743305e14/(1.44e10 - s^2); F = 1.034 + 3.17e-4/lam^2 (lam in um); valid above about 468 nm per Betremieux and Kaltenegger | 4.635e-27 | 1.000 | 2ND-HAND [R0282] |
| O2 | lam > 546 nm: 21351.1 + 2.18567e13/(4.09e9 - s^2); shorter: 20564.8 + 2.480899e13/(4.09e9 - s^2); F = 1.096 + 1.385e-3/lam^2 + 1.448e-4/lam^4; N at 273.15 K | 4.072e-27 | 0.879 | 2ND-HAND [R0282] |
| Ar | 6432.135 + 2.8606021e14/(1.44e10 - s^2); F = 1 | 3.991e-27 | 0.861 | 2ND-HAND (the B term reconstructed from the Betremieux consistency check, not read cleanly) |
| CO2 | (n-1)e8 = 1.1427e3 * [5799.25/(128908.9^2 - s^2) + 120.05/(89223.8^2 - s^2) + 5.3334/(75037.5^2 - s^2) + 4.3244/(67837.7^2 - s^2) + 1.218145e-6/(2418.136^2 - s^2)] (with the corrected last term); F = 1.1364 + 2.53e-3/lam^2 | 1.164e-26 | 2.51 | 2ND-HAND (Sneep and Ubachs); I re-evaluated it: n-1 = 4.272e-4, sigma 1.164e-26 |
| CH4 | n - 1 = 4.6662e-4 + 4.02e-6/lam^2 (lam in um); F = 1 | 1.28e-26 | 2.76 | 2ND-HAND; [R3661] says the Sneep and Ubachs formalism overestimates CH4 at 532.5 nm by 15 percent and [R2446] scales it down by 15 percent |
| CO | **do not use the coded Sneep and Ubachs form (22851 + 0.456e12/(71427^2 - s^2), F = 1.016); it gives 2.98e-27, half the measured value (exo_k's own comment says the same).** Use (n-1) = 32.7e-5 (1 + 8.1e-3/lam^2) at 273.15 K (Allen/Cox form in TauREx): 5.73e-27 at 550 nm, about 6.6e-27 at 532 nm | 5.7e-27 | 1.24 | 2ND-HAND, flagged; closer to the snippet value 6.19e-27 at 532 nm |
| H2O vapour | n - 1 = 4.92303e-2/(238.0185 - 1/lam^2) + 1.42723e-3/(57.3262 - 1/lam^2) for lam > 0.23 um; depolarization 2.99e-4 ([R3219] via Ranjan), F = 1.0010. Do not use d = 0.17 (that is liquid water) | n-1 = 2.362e-4 at 550 nm | about 0.5 | 2ND-HAND; reference temperature of the formula unknown |
| NH3 | n - 1 = 37.0e-5 (1 + 12.0e-3/lam^2) at 273.15 K (TauREx, citing Allen/Cox); no King factor found, F = 1 as coded | 7.4e-27 | 1.6 | 2ND-HAND, one code only, weakest entry |
| H2 | cross section directly [R0980]: sigma = 8.14e-13/L^4 + 1.28e-6/L^6 + 1.61/L^8 cm2, L in Angstrom (same in Caldas and petitRADTRANS; TauREx writes 8.14e-13 L^-4 (1 + 1.572e6 L^-2 + 1.981e12 L^-4), which agrees) | 9.38e-28 | 0.202 | 2ND-HAND |
| He | sigma = 5.484e-14 L^-4 (1 + 2.44e5 L^-2) cm2, L in Angstrom (pyratbay and TauREx agree); consistent with n-1 = 3.48e-5 at 1 um | 6.04e-29 | 0.0130 | 2ND-HAND |
| H2S, SO2 | NOT FOUND in refractivity form. [R3661] gives the route: sigma = 1.306e20 * KCF * alpha^2 / lam^4 (sigma in cm2, lam in um, alpha in cm3) with alpha and delta from [R0443]; those values were NOT REACHED | | | NOT REACHED |

**[R3470]** (H2 refractive index and dispersion, JOSA 67, 1550): NOT REACHED and the dispersion constants were not found anywhere reachable. The Dalgarno and Williams cross section replaces it. The implied H2 refractivity from it (1.369e-4 at 273.15 K) matches Betremieux and Kaltenegger's 1.37e-4, so the substitute is consistent with the one source we did read.

### 2.3 Betremieux and Kaltenegger [R0375]: VERIFIED (READ twice)

Table 1 caption "STP refractivities of various gases and planetary atmospheres at 1 um": H2 1.37e-4, He 3.48e-5, N2 2.95e-4, O2 2.68e-4, Ar 2.79e-4, CO2 4.43e-4 (also Jupiter 1.23e-4 and Earth 2.90e-4). Sources it states: N2, O2, CO2, Ar from its 2013 Table 2 ([R0282, R0374] read); He; H2 derived from the cross section. This is a 2015 paper by its arXiv number.

### 2.4 Relative cross sections against the prior claim (COMPUTED)

| Gas | Prior claim | Result | Verdict |
|---|---|---|---|
| CO2 | about 2.5 | 2.51 (2.48 using the 1-micron refractivities) | confirmed |
| H2 | about 0.22 | 0.202 (Dalgarno and Williams, 550 nm); 0.209 from 1-micron refractivities times King | corrected: use 0.20 |
| He | about 0.014 | 0.0130 (Kurucz/Dalgarno, 550 nm); 0.0135 from 1-micron | corrected: use 0.013 |

Ratios hardly change at 532 nm (CO2/N2 2.513, O2/N2 0.879, Ar/N2 0.861).

### 2.5 [R3480] (JOSA 47, 176), [R0282, R0575]

- **[R3480]**: NOT REACHED. 2ND-HAND (colour-science source): (n-1)e8 = 6432.8 + 2949810/(146 - s^2) + 25540/(41 - s^2) with s = 1/lam in um^-1, and a constant air King factor 1.0608. Its n is Edlen's 1953 form. At 550 nm Penndorf, (8342.13 + 2406030/(130 - s^2) + 15997/(38.9 - s^2)) and [R3472] (8060.51 + 2480990/(132.274 - s^2) + 17455.7/(39.32957 - s^2)) all give n - 1 = 2.778e-4. **Recommendation: do not cite Penndorf as the standard; use the Peck-Reeder air and King factor only for the Earth test, and the per-gas table above for the generator.**
- **[R0282]** (Planet. Space Sci. 32, 785): NOT REACHED. King factors above are 2ND-HAND from five independent codes that agree. Bates's own analytic cross-section fit was not found.
- **[R0575]** (Appl. Opt. 34, 2765): citation verified from reference lists; NOT REACHED. Its fit form is beta_s(lambda) = A lambda^(-B - C lambda - D/lambda), tau = beta_s (Ts/Ps) * integral of P/T dz with Ts = 288.15 K, Ps = 1013.25 mbar (READ in [R2441] appendix A and [R4333] eq. 9), but the A, B, C, D numbers (its Table 3) and its tabulated tau(550 nm) were NOT REACHED. The snippet "4.5 percent lower than Penndorf" is unverified; the coded formula difference accounts for only 1.1 percent.
- **Earth test value, COMPUTED independently by two of us (the worker and the lane):** Peck-Reeder air, F_air = 1.0488 (400 ppm CO2, 550 nm), N = 2.546899e19 cm^-3 gives sigma = 4.510e-27 cm2. Air column P/(m g) with P = 101325 Pa, m = 28.964 u, g = 9.80665 gives 2.148e25 cm^-2, so **tau_R(550 nm) = 0.0969**. With Penndorf's F = 1.0608 it is 0.0980. colour-science's own doctest is 0.0936 at 555 nm (different g). **The test target "about 0.097" stands as a computed value, not a read one.** Use tolerance 2 percent.

## 3. Directly usable cross-check table for tests

| Test | Value | Status |
|---|---|---|
| tau_R(550) of 1 bar Earth air | 0.0969 (0.0980 with Penndorf's F) | COMPUTED |
| tau_R at 350, 443, 550, 670, 865 nm (1 bar) | 0.62, 0.23, 0.10, 0.043, 0.015 | READ ([R1697: Table 1]); ratio 350/865 = 41 against (865/350)^4 = 37, consistent with an exponent about -4 |
| sigma_N2(550 nm) | 4.635e-27 cm2 | COMPUTED from the formula above |
| sigma_CO2 / sigma_N2 | 2.51 | COMPUTED |
| sigma_H2 / sigma_N2 | 0.20 | COMPUTED |
| doubling the pressure doubles tau | exact in the model | by construction |

## 4. Ozone

### 4.1 Chapman scaling and units: VERIFIED (READ)

- Reactions ([R0875], eq. 3.1; [R0739], Phil. Mag. XXXV): O2 + hv -> O + O; O + O2 + M -> O3 + M; O3 + hv -> O2 + O; O3 + O -> 2 O2.
- Steady state (Grenfell et al., [R1680] eq. 2): [O3] = [O2] * ( j_O2 * k(O+O2+M) * [M] / ( k(O3+O) * j_O3 ) )^0.5. It is a local density; it is not "ozone goes as the square root of O2". Chapman over-estimates ozone because catalytic loss is left out (same paper).
- 1 DU = 2.687e16 molecules cm^-2 (Cooke; 2.687e20 m^-2). Modern Earth 279 DU in Cooke's pre-industrial run. "About 300 DU" is a snippet-level statement; 279 to 300 DU is the test range.
- Hartley band becomes optically thick at about 1e17 molecules cm^-2 against about 1e19 for Earth ([R3785], READ), and about 0.1 percent PAL O2 gives a shielding layer there.

### 4.2 Column against O2 and the host star: the best checked relation

[R2466] ([R2466], READ; one-dimensional Atmos code, G0V, Sun, K2V, K5V, M5V, O2 from 0.01 to 150 percent PAL), Table 3: O3 column (cm^-2) = a x^5 + b x^4 + c x^3 + d x^2 + e x + f with x = log10(O2 mixing ratio, absolute, = 0.21 times PAL). Valid only for 0.01 to 150 percent PAL; clamp outside.

| Host | a | b | c | d | e | f |
|---|---|---|---|---|---|---|
| G0V | 4.582e16 | 6.021e17 | 2.481e18 | 2.618e18 | -1.126e18 | 5.742e18 |
| Sun | 6.203e16 | 7.716e17 | 3.147e18 | 4.157e18 | 8.412e17 | 4.640e18 |
| K2V | 3.417e16 | 3.846e17 | 1.283e18 | 7.391e17 | -8.266e17 | 3.734e18 |
| K5V | 6.826e15 | 2.379e16 | -3.424e17 | -1.848e18 | -9.440e17 | 4.898e18 |
| M5V | -2.144e16 | -3.381e17 | -1.933e18 | -4.458e18 | -1.941e18 | 4.542e18 |

I re-evaluated the Sun row: 192 DU at 100 percent PAL, 199 DU at 10 percent, 108 DU at 1 percent (the worker also got G0V 263, 277, 169 DU). It reproduces the paper's own quoted Sun values (4.91e18 at 150 percent against 4.92e18 in the text; 4.82e18 at 5 percent against 4.80e18). Kozakis finds hot stars peak at 25-55 percent PAL and cooler stars rise monotonically.

Cross-checks of the shape (all 2ND-HAND or READ):

| Source | Result |
|---|---|
| [R0875] (3-D WACCM6, READ) | 279, 169, 66 and 18 DU at 100, 10, 1 and 0.1 percent PAL |
| Cooke (2ND-HAND) | 266 DU at 10 percent PAL, 124 DU at 1 percent PAL |
| [R2279] via Cooke (2ND-HAND) | about 100 DU at 0.5 percent PAL; their profiles are an upper limit; Cooke's lower columns by factors 1.76 (10 percent) and 2.24 (1 percent) imply about 297 and 148 DU (that implication is the worker's arithmetic, not quoted) |
| [R3999] ([R3999], READ, Table 2, at Earth O2) | Sun 8.4e18 cm^-2 (313 DU); AD Leo 4.4e18 (164 DU); GJ 643 1.3e19 (484 DU); quiescent 3100 K 1.2e18 (45 DU), 3400 K 2.4e18 (89 DU), 3650 K 3.2e18 (119 DU). Quiescent M stars under 40 percent of Earth's |
| [R1680] (astro-ph/0610460, READ) | K2V 209-236 DU; Sun 296-327 DU |
| [R3845, R3846] ([R3846], READ, Table 1, cm^-2 inferred) | modern Earth F0V 1.15e19, F7V 8.62e18, Sun 5.29e18, G8V 5.56e18, K2V 3.62e18, K7V 3.46e18, M1V 4.03e18, M3V 3.59e18, M8V 2.47e18; 2.0 Ga (1 percent PAL) Sun 6.01e17, F0V 5.13e18; prebiotic Sun 1.84e15 |

The spread between models (Cooke 3-D lower than the 1-D codes by a factor about 2 at low O2; the Sun at Earth O2 from 192 to 313 DU) is real. The generator should draw within it rather than pretend precision: **GEN.413 should use the Kozakis fit for the median, host-star dependence by interpolating between the five hosts in T_eff, and a seeded spread of about +-30 percent** (a design choice for the lane, not a source value).

**[R2279]** (JGR 85, 3255-3263): original NOT REACHED; citation and numbers 2ND-HAND from Cooke. Use only as "1-D result; roughly 1 percent PAL is the usual shielding threshold; ozone rises with O2". Do not attribute the plateau between 0.1 and 1 percent PAL to it.

### 4.3 Cross sections

**[R0599]** (GRL 21, 581-584, [R0599] per snippet): original and the MPI Mainz files NOT REACHED. Snippet: Chappuis band, **410-760 nm**, at 220, 240, 260, 280 K relative to room temperature; negligible temperature effect at the peak (550-660 nm), about 40 percent drop at 420 nm from 298 to 220 K; scale from [R0110]. It does not contain the Hartley band.

**Stand-in actually read:** the Malicet/Brion 295 K ozone cross section (header "1995Malicet_1998Brion_O3_295.txt", 195-830 nm in 0.01 nm steps) inside the TUV-x package data (`musica` 0.17.1 sdist on PyPI; the package is Apache-2.0, the licence of the data files is NOT found). Computed from it:

| Feature | Value |
|---|---|
| Hartley peak | 1.149e-17 cm2 at 255.4 nm |
| Huggins (310-360 nm) | maximum 1.0e-19 cm2 at the 310.2 nm edge, falling steeply |
| Chappuis peak | 5.23e-21 cm2 at 602.7 nm |
| 450, 500, 550, 600, 650 nm | 1.9e-22, 1.2e-21, 3.4e-21, 5.1e-21, 2.5e-21 cm2 |

Band positions agree with [R0374] ([R0374]: Hartley 255 nm, Huggins 310-360 nm, Chappuis 602 nm) and [R3785] (255 nm Hartley, 400-700 nm Chappuis).

**[R4023]** (AMT 7, 625-636, [R4023]; title per snippet "High spectral resolution ozone absorption cross-sections - Part 2: Temperature dependence"; Serdyuchenko, Gorshelev, Weber, Chehade, Burrows): 213-1100 nm, resolution 0.02-0.24 nm, **eleven temperatures from about 193 to 293 K** (snippet, not the abstract; confirm). Data host is the IUP Bremen "Molecular Spectroscopy" page (footnote of [R4808]); blocked here; **licence not found**. Recommendation: ship a small fitted table (a few dozen wavelengths over Hartley, Huggins, Chappuis) and get the licence from Bremen before using the data.

### 4.4 Does ozone change the sky colour? COMPUTED (indicative)

Inputs: Malicet/Brion 295 K cross sections above; 300 DU = 8.061e18 cm^-2; Rayleigh as in section 2; 5772 K blackbody; single scattering below the ozone layer; ozone slant factor from a thin shell at 22 km (an assumption); Wyman et al. CIE fit from memory. Vertical optical depth for 300 DU:

| lambda (nm) | tau_O3 | tau_Rayleigh |
|---|---|---|
| 450 | 0.0015 | 0.221 |
| 500 | 0.0097 | 0.143 |
| 550 | 0.0271 | 0.097 |
| 600 | 0.0416 | 0.068 |
| 650 | 0.0202 | 0.049 |
| 700 | 0.0070 | 0.036 |

Hartley and Huggins vertical tau for 300 DU: 91 at 255 nm, 3.2 at 300 nm, 0.26 at 320 nm, 0.016 at 340 nm. Zenith sky transmission at 600 nm: 0.92 (sun at zenith), 0.88 (60 degrees), 0.77 (80 degrees), 0.58 (90 degrees); chromaticity shift (-0.004, -0.006) to (-0.024, -0.037), luminosity factor 0.95 to 0.73.

Conclusion: Chappuis absorption removes orange and red and shifts the sky slightly bluer/purer. By day it is minor (5-8 percent of luminance, shift 0.005-0.01 in xy); it matters at twilight and at low sun. [R2047] (JOSA 43, 113; two thirds of the zenith blue at 90 degrees solar zenith is Chappuis absorption) is SNIPPET only; NOT REACHED. Adams and Dobson: nothing found. **Design rule: apply ozone as Beer-Lambert absorption along the sun path (Chappuis in version 1, Hartley/Huggins only for the UV index); it is not a visible-colour driver at noon.**

## 5. Dust, haze, cloud and giants

### 5.1 Mars dust: [R2645]: original NOT REACHED

READ (Lemmon et al. on five Mars years of MER data, [R2644], the Icarus paper): effective radius 1 micron in northern summer, 2 micron at the onset of a dust lifting event; size distributions around 1.5-1.65 micron effective radius, variance 0.2-0.5 for Viking, Pathfinder and MER (citing [R2645, R3567, R4395]); opacity ratios correspond to an **Angstrom exponent near zero** (large particles); such aerosols have about 7-11 percent less opacity at 440 nm than at 880 nm; rover storm optical depths 1-1.5 (local), 1.24, and a peak of 1.98 at Opportunity. SNIPPET only: [R2645] abstract tau 0.9, radii 1.47 +- 0.21 and 1.52 +- 0.18 micron, decline 0.6-0.7 percent per sol (the decline is confirmed 2ND-HAND by 1403.4234), scale height 11.56 km.
**Not reached:** the single scattering albedo spectrum (only hobbyist numbers 0.690 at 450, 0.797 at 551, 0.922 at 615 nm turned up; they are SNIPPET and must not be cited). "Dust absorbs blue" is plausible (the imaginary index rises towards short wavelengths, snippet citing [R4395]) but not READ.
**Use:** dust extinction tau(lambda) = tau_ref (lambda/lambda_ref)^(-alpha) with alpha about 0 (grey), a 7-11 percent tilt between 880 and 440 nm; tau 0.5-1 typical, 1-1.5 for local storms, to about 5 for global storms (the last from the earlier analysis, not checked here); radius 1-2 micron. Colour comes from the single scattering albedo, which is **Boss's to research**.

### 5.2 Titan haze: [R4394] (Planet. Space Sci. 56, 669): original NOT REACHED

2ND-HAND: aerosols are fractal aggregates of about 4000 monomers of radius 0.04 micron ([R2008] review, [R2008], READ); optical depth above 80 km grows with decreasing height with a scale height of 65 km; wavelength fits from Figure 47 of [R4392, R4394] read in Buenzli et al. ([R0296], appendix A.4): tau80 = 1.012e7 lambda^-2.339, tau(80-30 km) = 2.029e4 lambda^-1.409, tau(0-30 km) = 6.270e2 lambda^-0.9706 (lambda in nm presumably). **Do not use these as absolute values:** at 500 nm they give 4.9 + 3.2 + ... above 8, against a snippet DISR extinction of about 4.5 at 531 nm. Normalisation or units are unclear and [R1144, R4394] were not reached. Verified for use: the exponents (haze extinction roughly lambda^-1 to lambda^-2.3, steeper aloft) and "orange" colour (aerosols give Titan its orange colour, [R2008]).

### 5.3 Giant planets: [R2258]: original NOT REACHED

2ND-HAND ([R0622: Table 3], READ): CH4 absorption features near 0.46, 0.48, 0.54, 0.62, 0.78, 0.79, 0.84, 0.86, 0.89, 0.91, 0.99 micron; Na 0.59; K 0.40 and 0.77; H2O 0.94. [R4250] (READ): methane bands shortward of 1 micron, strongest near 0.73, 0.86, 0.89 micron. **Band strengths and absorption coefficients were not obtained:** Boss must supply them ([R2258] tables, [R2263] and 2010) or ask for a file.

READ [R2101]: Uranus and Neptune are blue because of methane (about 4 percent by mole, citing [R2263], 2011), which absorbs red and infrared, plus Rayleigh scattering of a low-aerosol atmosphere; with no methane and no aerosol they would be pearly white; aerosols make them paler (Uranus is paler than Neptune from a thicker micron-sized Aerosol-2 layer). Gas giants are less bright in the blue; ice giants are brightest in the blue (Cahoy, citing [R2258] data).

### 5.4 Cloud classes: [R4250] (ApJ 538, 885, astro-ph/9910504): READ, VERIFIED with corrections

| Class | T_eff (2000 paper) | Cloud | Look |
|---|---|---|---|
| I | below about 150 K | NH3 above a deeper water cloud | Bond albedo about 0.4-0.65; reflective in the visible, lower in the IR |
| II | about 250 K | water cloud only | highest visible albedo, Bond up to nearly 0.9 |
| III ("clear") | above about 350 K (**2003 paper: T_eq 350-800 K**) | none; Na and K absorb | spherical albedo under 0.6 in UV/blue, under 0.1 in the red; dark |
| IV | 900-1500 K | silicate deck at 5-10 bar, deep iron cloud | lowest albedo; geometric albedo about 0.03 at 0.48 micron (tau Boo b) |
| V | above about 1500 K and/or very low gravity | silicate layer high up | geometric albedo 0.39 at 0.48 micron; Bond about 0.51-0.57 |

Transitions at 150, 250, 350 K are approximate, not sharp. Whole sample Bond albedo about 0.02 to 0.8. "Less condensation, larger particle sizes and wider size distributions result in lower albedos." [R4249] (ApJ 588, 1121, astro-ph/0210216, READ): classes defined by T_eq; Class I profile meets the NH3 curve near 1 bar and the water curve near 10 bar; Class II water clouds at about 1-2 AU. Jupiter's ammonia deck about 0.3-0.7 bar.

[R0622] (READ), solar-analog star, geometric albedo B, V, R, I and Bond albedo for 1x-solar Jupiter-like planets: 0.8 AU (clear) 0.556, 0.322, 0.156, 0.047, Bond 0.408; 2 AU (water cloud) 0.733, 0.742, 0.727, 0.673, Bond 0.914; 5 AU (NH3 over water) 0.609, 0.567, 0.531, 0.453, Bond 0.675; 10 AU 0.450, 0.386, 0.358, 0.313, Bond 0.489. A 10x-Neptune-metallicity planet at 5 AU: 0.434, 0.326, 0.238, 0.145, Bond 0.407. Mie scattering extinction is essentially constant for radii over 1 micron (so clouds are white/grey). Water cloud optical thickness can reach about 10,000. Water-cloud base about 500 mbar (2 AU), NH3 near 50 mbar and water near 1 bar (5 AU).

 (astro-ph/9810073, abstract READ): Rayleigh scatter dominates in the blue; molecular absorption removes red and infrared; cloud-free atmospheres above about 400 K are dark beyond 0.6 micron.

**Use for GEN.419 (giants):** pick the class from T_eff (or T_eq) by the table, take the B, V, R, I colour from the Cahoy table nearest in class, interpolating in orbital distance and metallicity; do not fix a colour per class letter without that table (the project rule). Cahoy Table 6 has all 16 models.

### 5.5 Condensation: partly verified

READ ([R2731], astro-ph/0511136, eq. 32 and 34, P in bar, solar composition): NH4SH: 10^4/T_cond = 48.91 - 4.15 log10(P) - 4.15 [M/H] (1 bar solar: 204.5 K, COMPUTED by the worker); H2S: 10^4/T_cond = 86.49 - 8.54 log10(P) - 8.54 [S/H] (1 bar solar: 115.6 K); H2S condenses only if NH3 is absent or less abundant than H2S, but NH3 is about 4.4 times H2S in protosolar gas, so NH4SH takes the sulfur first. READ [R2101]: methane condenses at about 1.5 bar on both ice giants; H2S above about 3 bar; NH4SH near 40 bar. **NOT REACHED:** the NH3, H2O and CH4 condensation curves from [R2726]/2010, [R2725, R4653]. Boss to research, or use the Clausius-Clapeyron forms the repo already has for H2O (`vapour_pressure_kpa`) and take NH3 and CH4 from NIST/Lodders.

### 5.6 Clouds in general

- Venus: [R0157]: optical depth at 0.63 micron: upper haze 0.2-1, upper cloud 6-8, middle 8-10, lower 6-12, lower haze 0.1-0.2; the three main layers sum to about 20-30. **"30-40" is not verified.** Droplets 75-85 percent H2SO4; modes 0.49, 1.04, 3.85 micron; Bond albedo 0.76 typical, range 0.7-0.82, 0.81 in models ([R1569], READ). [R3014, R3016]: clouds reflect over 70 percent and only about 3 percent reaches the surface.
- Earth clouds ([R1697: Table 1], READ): model cloud optical thickness bins 0, 5, 10, 20 (MODIS bins), droplet radius 10-17.5 micron, tau 9.775 to 10.002 from 350 to 865 nm (grey). The ISCCP global mean about 3.8 is SNIPPET.
- Mie extinction of large particles is grey (Cahoy, Lemmon READ). Titan's small monomers are steep. **A reachable Angstrom-exponent review [R3982] was NOT REACHED**; the range "0-2, under 1 for dust and sea salt" is SNIPPET.

### 5.7 [R3014]: READ (Sect. 4.6), VERIFIED with a nuance

"Worlds with clear sky paths, i.e. without global cloud or haze cover, appear blue despite the redness of Proxima Centauri" (small blue enhancement in the star's spectrum, steep Wien tail, blue photons Rayleigh-scattered as on Earth); for 10-bar O2 cases the blue enhancement overwhelms the redness; modern Earth "appears lavender"; "the haze-covered Archean Earth and the Venus-like cases are redder than modern Earth and are surprisingly lilac-colored since they do not exhibit a significant spectral impact from Rayleigh scattering"; Proxima b "might appear as a pale blue, or pale purple, dot". This is the **planet's disk colour from space** (top-of-atmosphere spectrum with the CIE matching functions), not the colour of the sky from the ground. **The sky analysis note's wording "a direct precedent for the method" overstates it; it is a precedent for the star-plus-Rayleigh reasoning only.**

## 6. The other recalled items

| Item | Result |
|---|---|
| [R2274] (Appl. Opt. 28, 4735-4738) | 2ND-HAND from pvlib source (raw.githubusercontent.com): m = 1/(cos z + 0.50572 (6.07995 + (90 - z))^-1.6364), z in degrees, apparent zenith. COMPUTED: m = 0.9997 at 0, 1.994 at 60, 10.31 at 85, 37.9 at 90 degrees. It is an Earth fit (ISO 1972 standard atmosphere, refraction included), so it is **only a regression target for Earth**; other worlds need the Chapman function with the planet's R/H. |
| Chapman function ([R0739], Proc. Phys. Soc. 43, 483) | existence and use confirmed by a snippet and by [R4517] (Earth, Planets and Space; best analytic approximations named: Green and Barnum); the formula itself NOT REACHED. The lane must take the exact integral, which is easy to evaluate numerically for an exponential atmosphere. |
| A textbook | textbook; the ground-bounce factor 1/(1 - A_g s) was not looked up. NOT REACHED. Any radiative transfer text gives it; Boss may prefer a two-stream pass instead. |
| [R0560] | Computer Graphics Forum 27(4), 1079-1086 (search metadata, not EGSR); [R1940]: CGF 39(4), 13-22, EGSR 2020 (search metadata). Citations verified; contents not needed for the chain. Hillaire's reachable copy at diglib7.eg.org is refused here. |
| [R4788] | already validated numerically in `docs/design/orbits/sky-view.md`; not rechecked. |
| [R3472] (JOSA 62, 958) | appeared as a reference in a data README (search); formula read in colour-science (section 2.5). |

## 7. What Boss must research himself (what could not be reached)

Ranked by how much the chain depends on it:

1. **Single scattering albedo spectrum of Mars dust** ([R2645, R4395] Mars Pathfinder), and the refractive index of iron-oxide dust. Without it dust colour is a guess.
2. **Titan haze absolute optical depth and spectrum** [R1144, R4392]: the exponents are usable, the normalisation is not. Needed for any Titan-like haze world.
3. **Methane (and ammonia, water) absorption strengths in the visible** for the giants ([R2258, R2263], 2010). Band positions are known, strengths are not.
4. **Ozone cross sections with a clear licence**: [R4023] data (IUP Bremen) and the [R0599] files (MPI Mainz atlas). The Malicet/Brion numbers above are real but their file licence was not found. A small fitted table needs the licence only if redistributed.
5. **Refractivity and King factors for H2S, SO2 ([R0443] polarizabilities), NH3 and CO** from the original papers, and a clean check of the S and U CO2 value at 532 nm (computed 13.3e-27 cm2 against a snippet value of 12.4e-27; unresolved, 7.6 percent).
6. **Condensation curves for NH3, H2O and CH4** ([R2725, R2726] and 2010, [R4653]) for the cloud-base rule.
7. **[R0575: Table 3]** (A, B, C, D) and its optical depth table, to replace the recomputed Earth test value with the published one.
8. **[R3470]**: the H2 and He dispersion; the replacement is consistent at 1 micron but is a cross-section law.
9. **[R2279]** originals; **[R2047]**; **an Angstrom-exponent review**; the Venus cloud total.
10. Whether the proxy can be widened: the hosts to add are doi.org, ui.adsabs.harvard.edu, www.science.org, amt.copernicus.org, acp.copernicus.org, opg.optica.org, www.sciencedirect.com, ntrs.nasa.gov, and the Mainz and Bremen data hosts (a Project settings > environment change).

## 8. Direct answers for the lanes

- **GEN.417 (Rayleigh of the mixture):** use the section 2.2 table with per-gas reference temperature; sigma_mix = sum f_i sigma_i; tau = N_col sigma_mix with N_col = P/(m_bar g) from the stored partials. Tests from section 3.
- **GEN.413 (ozone):** Kozakis fit by T_eff, clamp 0.01-150 percent PAL, seeded spread; Chappuis cross sections from a fitted table; UV index from the Hartley band; keep the generator free of a fixed Earth column (Earth about 279-300 DU is a test value).
- **GEN.414 (haze, dust, cloud):** grey dust (alpha about 0, radius 1-2 micron, tau 0.5-1.5 with storms), steeper haze (lambda^-1 to -2.3), cloud tau from the Venus/Earth ranges with the cloud type from the condensation rule; the dust colour and haze absolute depths wait for Boss's research.
- **GEN.419 (giants):** Sudarsky class by T_eff, colour table from Cahoy Table 6, ice giants blue by methane bands plus Rayleigh plus pale aerosols; band strengths are missing.

## 9. How the numbers above were checked

Where I say COMPUTED, I re-ran it here with Python: sigma_N2(550) = 4.6346e-27 cm2, sigma_CO2(550) = 1.1637e-26 cm2 (ratio 2.511), Earth air column 2.148e26 cm^-2 and tau_R = 0.0969, Kozakis Sun fit 191.8, 199.4 and 108.0 DU at 100, 10, 1 percent PAL, Kasten and Young airmass at 0, 60, 85, 90 degrees. The rest of the computed figures (H2 and He ratios, the ozone sky table, the Malicet/Brion peaks, the O2 and CO numbers) come from three read-only reading passes whose scratch files were not kept in the project; they can be rerun from the cited arXiv PDFs and PyPI packages.
