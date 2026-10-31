# Boss's papers: what they fill, what stays open (2026-10-17)

Status: research result. Written 2026-10-17 by the paper-filing thread (about 12:45 PM EDT as that thread's own clock read; the coordinator relayed it at 12:34 PM EDT by this thread's clock, so the lane clock runs ahead) and filed by the TODO thread on 2026-10-17. Boss uploaded 29 files between 12:10 and 12:12 PM EDT and said at 12:13 PM EDT this is the last batch he expects; 28 are kept in `/mnt/project-files/papers/` (index: `INDEX.md` there). The reader files with the table, equation, figure or page of every number are in `/mnt/project-files/reports/boss-papers-2026-10-17-files/`. Where a paper contradicts a second-hand value in our files the paper wins and both values are stated (section 3). The items that use each number carry it in their notes in [TODO.md](../../TODO.md) and the gap plan has the summary in [section 4.10](gap-filling-plan.md).
This report was read-only work: no code and no edits to TODO.md, issues, docs or docs/plan; the TODO thread filed it.

Inputs: 29 uploaded files in five messages between 12:10 PM and 12:12 PM EDT, plus Boss's "Pay Only" list (12:11 PM EDT) and his note at 12:13 PM EDT that this is the last batch he expects. One file is a byte-identical second copy (Sneep and Ubachs), so 28 files are filed in `/mnt/project-files/papers/` (see `INDEX.md` there for citation, items and status of each). The gap items are those of `docs/design/research/gap-filling-plan.md`; the requests are R1-R15 of GEN.531 (GR1-GR15 in the plan).

Detail behind every number: five reader files in `/mnt/project-files/reports/boss-papers-2026-10-17-files/` (A optics, B1 sulfuric acid and forcing, B2 radiation, comets, outgassing, C body, D galaxy). Each lists the table, equation, figure or section and printed page of every number.

## 0. How the numbers were obtained

- Five reading workers each read their papers (text layer where it exists; page images for the scans, which are King, Bates, Tabazadeh, Zeleznik and Myhre). They marked **V** = read from the paper's own text or table image, **R** = recomputed by them (arithmetic in their file), **D** = digitized from the geometry of a figure, **FIG** = read off a plotted curve (about plus or minus 10 to 20 percent). I use the same marks below. Nothing is filled from memory except where a file says RECALLED.
- I re-read these myself against the PDFs and they agree with the workers: Hassler (0.210 mGy/day, 0.64 mSv/day, Q 3.05, 76 mGy/yr, 232 mSv/yr); Palme and O'Neill Table 3 oxides and total 98.41, Th 83.4 ppb, U 21.8 ppb; Aumer and Binney beta values and the Table 2 best-fit row; Bierson and Nimmo Table 1 (920, 3500, 0.60, 0.20, -0.1, -0.11); Sneep and Ubachs Table 2; Irwin 0.4-0.8 ppm and 4.4-5.0 times solar; Graus Eq. 1 with V50 = 18 km/s and sigma 2.5; Tabazadeh Eq. 1 and Table 1 rows 10 to 40 (page image); Mumma Table 4, SO2 row 0.2; the first page of Myhre. Numbers I did not re-read stay with the worker's V mark and are not independently double-checked.
- Standing rule applied throughout: where a paper contradicts a second-hand value in our files, the paper wins and both values are stated (section 3).

## 1. Bottom line

Counts are of the plan's gap ids and requests that these papers touch.

**Filled from the papers (numbers in section 2):**
- AD4 Rayleigh, for air, N2, O2, Ar, CO2 (and CO, CH4, N2O, SF6 at 532 nm only); the Earth tau_R test value 0.09721.
- AD7 water side of H2SO4-H2O (Tabazadeh) and the acid and water activities (Zeleznik).
- AD3 greenhouse forcing laws (Myhre).
- AD5 Mars surface dose (Hassler).
- BC4 crust rows (Rudnick and Gao) and mantle oxides, K, Th, U of the first-edition Palme and O'Neill.
- BC7 ice compaction law and the ice and rock densities (Bierson and Nimmo).
- BC3 mineral parameters of Part II (Stixrude and Lithgow-Bertelloni); the 100 GPa test passes.
- GR4 for the Milky Way (Aumer and Binney); GR6 and GR3 for the Milky Way (Haywood); GR15 form and width (Graus).

**Partly filled:**
- AD2 (Wolff has figures only, no tables).
- AD7 acid vapour pressure (no paper here prints it).
- BC3 (Part I and the 2024 Part III set, which our files use, are not here).
- BC4 mantle (2003 edition only).
- GR6 and GR3 beyond the Milky Way (Pinna 2019a covers one galaxy).

**Still open after these papers:**
- AD4 for SO2, H2S, H2, He, NH3 (needs Bogaard et al. 1978, not uploaded).
- AD1 giant-planet visible absorption (Karkoschka, paid, not uploaded).
- AD2 Titan haze (Tomasko 2008, paid, not uploaded).
- AD8 photochemistry rates (JPL 19-5, not uploaded; the Curchod review has no numbers).
- BC12 (Fray and Schmitt, paid).
- BC13 wR, wT, dose limits (the ICRP extract stops before the tables; the full report is paid).
- AT.k7 N2, O2, Ar, Kr, Xe in comets (Mumma has none of them).
- AT.a1, a2 mass law and later-mantle law (Gaillard and Pahlevan have none).
- GR4 for other galaxies, GR1, GR10, GR12, GR14, GX5 and every row Knapen's book or Gowanlock leaves untouched (GR5, GR7 to GR9, GR11, GR13).

**Will be filled with labelled stand-ins from the gap plan** (the plan's D, B, T, P fills stay as written, flagged "stand-in", replaceable when a source arrives): everything in the open list above, plus AD2 spectra if Boss accepts digitized values (section 2.2), and the acid vapour pressure until a pure-acid curve is verified.

## 2. The numbers the gap plan needs

### 2.1 AD4 Rayleigh (GEN.422)

King factor: F_K = (6 + 3 rho)/(6 - 7 rho), rho = natural-light depolarization ratio (King 1923 Eq. 78 and 90, p. 351 and 353; Sneep Eq. 7, p. 296; Bucholtz Eq. 6). [V] The mixture rule from King Eq. 91 and 94: F_mix - 1 = sum over gases of w_r (F_r - 1), with w_r = v_r (mu_r - 1)^2 / sum(v (mu - 1)^2), equivalently sigma_mix = sum v_r sigma_r. [R, worker A] Inverse: rho = 6 (F - 1)/(3 + 7 F). [R]

Air (dry, 0 C, 760 torr) from Bates 1984 Table 1 (printed p. 787) and Bucholtz 1995 Tables 1-3 (pp. 2767-2769): [V]

| lambda (micron) | Bates Q_R (cm^2) | Bates F'_K | Bucholtz sigma (cm^2, 15 C) |
|---|---|---|---|
| 0.30 | 5.656e-26 | 1.055 | 5.642e-26 |
| 0.40 | 1.673e-26 | 1.051 | 1.673e-26 |
| 0.50 | 6.663e-27 | 1.047 printed (probably misprint for 1.049) | 6.656e-27 |
| 0.55 | 4.513e-27 | 1.049 | 4.509e-27 |
| 0.70 | 1.695e-27 | 1.048 | 1.692e-27 |
| 1.00 | 4.022e-28 | 1.047 | 4.010e-28 |

Bucholtz analytic fit, Table 3 (p. 2769): sigma = A lambda^-(B + C lambda + D/lambda), lambda in micron, cm^2: A = 3.01577e-28 (0.2-0.5 micron) and 4.01061e-28 (above 0.5); B = 3.55212 and 3.99668; C = 1.35579 and 1.10298e-3; D = 0.11563 and 2.71393e-2; stated accuracy under 0.4 percent at 0.2 micron, under 0.1 percent above 0.5. [V] Refractivity (Bucholtz Eq. 4-5, p. 2766) and N2/O2 refractivity and King-factor formulas (Bates Eq. 6-14; N2 F = 1.034 + 3.17e-4/lambda^2; O2 F = 1.096 + 1.385e-3/lambda^2 + 1.448e-4/lambda^4, lambda in micron) are in detail file A sections 3-4. [V]

Earth test value, Bucholtz Table 4 (p. 2770): tau_R at 550 nm = 0.09721 for the 1962 US Standard atmosphere (0.0976 tropical, 0.09688 subarctic summer). [V] The plan's computed 0.0969 is 0.3 percent lower, inside the plan's 2 percent tolerance; the published value can replace the computed one as the test target. Bucholtz's own depolarization column is back-computed from Bates's rounded King factors, so Bates and Bucholtz are one source, not two. [R]

Measured cross sections, Sneep and Ubachs Table 2 (p. 304), 532.2 nm, 15 C, 1013 hPa, units 1e-27 cm^2: [V]

| Gas | measured | n-based | King factor used |
|---|---|---|---|
| Ar | 4.45 (0.3) | 4.56 | 1 |
| N2 | 5.10 (0.24) | 5.30 | 1.035 |
| CO | 6.19 (0.4) | 6.82 | 1.016 |
| CO2 | 12.4 (0.8) | 13.29 | 1.145 |
| CH4 | 12.47 (0.23) | 14.69 | 1.000 (assumed) |
| N2O | 15.90 (0.08) | 18.19 | 1.225 |
| SF6 | 32.3 (0.5) | 34.1 | 1.000 |

The paper's own statement (abstract and p. 306): the refractive-index method is confirmed "at the few percent level", with larger gaps for CO, N2O and CH4. The plan's "5 to 9 percent generic error" is therefore wider than the measured 2.5 to 4 percent for Ar and N2 and 7 percent for CO2, and narrower than the 10 to 18 percent for CO, CH4 and N2O. [R, worker A, from Table 2]

**Open: SO2, H2S, H2, He, NH3 (and water vapour).** None of the five Rayleigh papers has a depolarization ratio for them. King 1923 gives only air, O2, N2, Ar and obsolete white-light values. Bogaard, Buckingham, Pierens and White 1978 is cited by the plan and by Summa 2021 and is not uploaded. The +-5 to 9 percent generic Rayleigh error of the plan stays for these gases. Summa 2021 gives static King factors only for N2 (1.036) and CO (1.013) from theory, which agree with Bates (1.034) and Sneep (1.016). [R]

### 2.2 AD2 Mars dust (GEN.423)

Wolff 2009 prints **no table** of single scattering albedo, n or k. Printed in the text and Table 3 (p. 12): solar-band w0 = 0.922 to 0.926 (+-0.003) for effective radius 1.2 to 1.8 micron; 0.940 +-0.003 for the "1.8 refractive index, 1.4 scattering" case; literature 0.89 to 0.94. [V] The wavelength spectra of w0, n and k for the four radii (30 channels, 440 to 2920 nm) exist only in Figures 12, 13, 14. Worker A recovered them from the vector geometry of the PDF (the n = 1.5000 fixed point at 500 nm and the round CRISM channel wavelengths both reproduce, which checks the calibration). Accuracy about 0.0005 in w0, 0.0002 in n, 0.0001 in k. [D] That table is detail file A section 6.5. Sample, reff = 1.8 micron: 500 nm w0 0.8310, n 1.5000, k 0.00549; 700 nm w0 0.9712, n 1.4987, k 0.00095; 2400 nm w0 0.9659, n 1.4890, k 0.00453. The paper's other source for Mars dust indices (Tomasko et al. 1999, Table 8) and Ockert-Bell 1997 are not here. Titan haze (Tomasko 2008) is paid and not uploaded.

### 2.3 AD7 sulfuric acid and water (GEN.471, 491)

Tabazadeh et al. 1997 (scan, all four pages): [V]
- Eq. 1 (p. 1932): ln P0_H2O = 18.452406985 - 3505.1578807/T - 330918.55082/T^2 + 12725068.262/T^3, P0 in mbar, T in K, valid 185 to 260 K (supercooled liquid water). Re-read by me on the page image.
- Table 1 (p. 1932): water partial pressure P_H2O (mbar) = exp(a + b/T + c/T^2) for 15 compositions from 10 to 80 weight percent H2SO4 (for example 40 percent: a 20.379, b -4828.5, c -112550; 60 percent: 20.405, -5252.1, -100840), valid 185 to 260 K, linear interpolation between compositions, within 2 percent of Clegg and Brimblecombe. Rows 10 to 40 re-read by me; the rest are the worker's.
- Eq. 2-3 and Table 2: acid molality and weight percent as a function of water activity and temperature, three branches. The branches are not continuous at a_w = 0.05 and 0.85 (steps of about 0.2 to 1 weight percent). Test at 200 K, 60 weight percent: Eq. 1 gives P0 = 3.165e-3 mbar and Eqs. 2-3 give a_w = 0.0735. [R]
- **No acid vapour pressure is given.** The plan's line asking Tabazadeh for the acid vapour pressure by weight percent is not met by this paper.

Zeleznik 1991 (scan, 42 of 44 pages, printed 1180-1181 missing; they hold only solid-hydrate heat-capacity figures): [V]
- Gibbs energy function of the binary solution, Eq. 12-13 with the 50 parameters of Table 6A (pp. 1186-1188), valid 200 to 350 K at 0.1 MPa; the worker's transcription reproduces Table 7 to 4 decimals at 200, 298.15 and 350 K. Eq. 12 is printed with "ln x_i" and only works as ln x_j; Eq. 11's second line is printed without its minus sign. Both corrections are the worker's recomputation, not the paper's.
- Tables 7 and 8 (pp. 1189-1195) give activities and osmotic coefficients by mole fraction at 200, 250, 298.15, 300, 350 K and by molality at four temperatures. **No partial pressures**: the paper holds no vapour phase.
- Against Tabazadeh the water activity differs by 10 to 22 percent at 60 to 70 weight percent and 200 to 220 K (Zeleznik lower) and by under 5 percent at 40 to 50 percent and at 250 K. [R] Tabazadeh was built for the 185-260 K stratospheric range, so it is the natural choice there.
- The acid vapour pressure still needs a pure-acid curve times the acid activity; that curve is Ayers 1980 (paid, not uploaded), or the Kulmala-Laaksonen form our files already carry.

### 2.4 AD3 greenhouse forcing (GEN.413, 414)

Myhre et al. 1998 Table 3 (printed p. 2718): [V]

| Gas | Expression, dF in W m^-2 | alpha (IPCC 1990) | alpha (this paper) |
|---|---|---|---|
| CO2 | alpha ln(C/C0), C in ppmv | 6.3 | **5.35** |
| CH4 | alpha (sqrt M - sqrt M0) - (f(M,N0) - f(M0,N0)), M in ppbv | 0.036 | **0.036** |
| N2O | alpha (sqrt N - sqrt N0) - (f(M0,N) - f(M0,N0)), N in ppbv | 0.14 | **0.12** |
| CFC-11 | alpha (X - X0) | 0.22 | **0.25** |
| CFC-12 | alpha (X - X0) | 0.28 | **0.33** |

f(M,N) = 0.47 ln[1 + 2.01e-5 (MN)^0.75 + 5.31e-15 M (MN)^1.52]. Uncertainties assumed: 1 percent CO2, 10 percent CH4, 5 percent N2O and CFC-11, 3 percent CFC-12. CO2 doubling gives 3.71 W m^-2 [R]. Our parts file labels 5.35 and 0.12 as "IPCC 1990"; IPCC 1990 had 6.3 and 0.14 (printed in the same table), so those labels should read Myhre 1998. The paper states no validity range and prints no concentrations behind its forcing tables (they are in IPCC 1995). Table 2's printed total 2.247 does not equal its 15 rows (2.2537; it equals the sum without the CH3CCl3 row). Infrared opacity grids (correlated-k) are a download job, not in these papers.

### 2.5 AD5 and BC13 radiation (GEN.436, 444)

Hassler et al. 2014 Mars surface, GCR only, about 300 sols (Table 1, p. 1244797-3): dose rate 0.21 +-0.04 mGy/day, mean quality factor 3.05 +-0.26, dose-equivalent rate 0.64 +-0.12 mSv/day; MSL cruise 0.48 mGy/day, Q 3.82, 1.84 mSv/day. Table 3: 76 mGy/yr and 232 mSv/yr at the surface. One SEP event 0.025 mGy (Table 2) against about 50 microGy in the text (unexplained). The measured column at Gale is about 21 g/cm^2 (altitude -4.4 km). Our file's code point (16.4 g/cm^2, 237 mSv/yr) matches the value but not the column. [V] The paper uses the ICRP 60 quality function Q(L), not ICRP 103 weighting factors.

**BC13 is not closed.** Neither ICRP upload holds any wR, wT or dose limit: the free extract (35 pages) runs front matter, preface, executive summary, glossary and report pages 1 to 34, and the tables are in section 4 (printed pp. 61-78) and section 5.10 (from p. 98). The "preface, executive summary and glossary" upload is a strict subset of the extract. What both do hold: the Q(L) function (Q = 1 for L < 10 keV/micron; 0.32 L - 2.2 for 10 to 100; 300/sqrt(L) above 100), the definitions, "about 5 percent per Sv" and a 10 mSv radon reference level. [V, worker B2]

### 2.6 AT.k7 comets (GEN.488, 495)

Mumma and Charnley 2011, Table 4 and Figure 4: no N2, O2, Ar, Kr or Xe appears anywhere. SO2 is listed for Hale-Bopp only, 0.2 percent of water (re-read by me), giving a second comet next to 67P (0.127 percent in our files). Table 4 columns are Halley, Hyakutake, Hale-Bopp, C/2001 A2 and 73P/SW3 (re-read by me): NH3 is 1.5, 0.5, 0.7 and below 0.3 percent of water in the four comets that have it, and H2S 0.4, 0.8, 1.5, 1.15 and 0.25, so the repo's NH3 upper bound of 3.6 percent is not supported by this review. Three infrared taxonomy classes (organics-normal, enriched, depleted). Figure-read ranges are approximate. [V for table values, FIG for the figure]

### 2.7 AT.a1, a2 mantle redox (GEN.482, 498)

No mass law and no later-mantle law in either paper. Gaillard and Scaillet 2014 (HAL manuscript): equations (1) to (11) without coefficients; speciation against pressure only as curves; its text claims convergence to FMQ-2 +-0.5 at 10 to 0.01 bar while its own Figure 2 shows about FMQ-2.75 with 0.7 dex spread at 0.01 bar; the source fO2 is quoted three ways (FMQ-1.7, -1.5, -1.4). Pahlevan, Schaefer and Hirschmann (submitted manuscript): Earth only, last magma-ocean equilibration above IW+1 (H2/H2O below 0.3), above the Sossi Earth value (IW+0.5) our files use; its Figure 5 top axis contradicts its own Figure 2 by a factor of 6 to 9 (13 bar against about 116 bar of H2 at IW+1). The plan's fills for these rows stay stand-ins. Hirschmann 2012 and Gaillard et al. 2021/2022 are not uploaded. [V and FIG, worker B2]

### 2.8 BC4 oxides, K, Th, U (GEN.428)

Rudnick and Gao 2014 Table 11 (p. 43), oxides in weight percent, Th and U in ppm: [V]

| | Upper | Middle | Lower | Total |
|---|---|---|---|---|
| SiO2 | 66.6 | 63.5 | 53.4 | 60.6 |
| TiO2 | 0.64 | 0.69 | 0.82 | 0.72 |
| Al2O3 | 15.4 | 15.0 | 16.9 | 15.9 |
| FeOT | 5.04 | 6.02 | 8.57 | 6.71 |
| MgO | 2.48 | 3.59 | 7.24 | 4.66 |
| CaO | 3.59 | 5.25 | 9.59 | 6.41 |
| Na2O | 3.27 | 3.39 | 2.65 | 3.07 |
| K2O | 2.80 | 2.30 | 0.61 | 1.81 |
| Th ppm | 10.5 | 6.5 | 1.2 | 5.6 |
| U ppm | 2.7 | 1.3 | 0.2 | 1.3 |

K is not tabulated; K2O x 0.8301 gives 23,244 / 19,093 / 5,064 / 15,026 ppm [R]. The printed K/U row is 7 to 10 percent off K2O-derived K over printed U, so use K2O and U directly. The chapter is the 2003 text under a 2014 cover (no citation after 2004 found), so there is no 2003-against-2014 list. **Every crust value in our files matches this table exactly.** Palme and O'Neill's own crust column (Table 5: Th 3,500 ppb, U 910 ppb, K 9,100 ppm) disagrees with Rudnick and Gao; Rudnick's values win for the crust.

Palme and O'Neill 2003 (first edition, not the 2014 edition the plan named): Table 3 (p. 12), primitive mantle, weight percent: MgO 36.77 +-0.44, Al2O3 4.49 +-0.37, SiO2 45.40 +-0.30, CaO 3.65 +-0.31, FeOt 8.10 +-0.05, **total 98.41, not 100** (the paper's Eq. 6 sets the five to 98.41; the other 1.59 is minor oxides 1.34 plus under 0.15 trace); Mg# 0.890. Mg/Si 1.045 and Al/Si 0.112 are mass ratios. Table 4: K 260 ppm, Th 83.4 ppb, U 21.8 ppb. [V; re-read by me] Our files' "Palme and O'Neill 2014" Th 0.0849 and U 0.0229 ppm are 1.8 and 5.0 percent above these; this file cannot confirm them. Decision for TODO: either keep the 2014 numbers labelled second-hand until the 2014 edition arrives, or use these 2003 numbers labelled "2003 edition". The verified source wins under the standing rule.

### 2.9 BC3 mantle minerals (GEN.425)

Stixrude and Lithgow-Bertelloni 2011 (Part II): the only parameter table is A1 (pp. 1207-1208): 45 phases, with A2 (interaction) and A3 (Landau) parameters. The worker transcribed A1 in full, not A2 (A1 is in detail file C section 3). The paper has no density-at-pressure table; the BC3 test "density at 100 GPa within 2 percent" was recomputed from Table A1 with a cold Birch-Murnaghan isotherm and passes: perovskite within 0.05 percent of our copy, periclase within 0.6 percent. [R] Our mantle parts table is the 2024 Part III set, not Part II (for example Mg-post-perovskite K0 247.7 in ours against 231 here; Mg-perovskite q0 0.834 against 1.1). The authors say parameters at full precision are "available on request" and that rounded values do not reproduce high-pressure results. Part I (2005) and Part III are not here. BC3 stays "not a blocker".

### 2.10 BC7 ice compaction and density (GEN.442, 457)

Bierson and Nimmo (accepted manuscript), Table 1 and Eq. 1-2: ice density 920 and silicate density 3500 kg/m^3, initial porosity 0.60, strength-supported porosity 0.20 (ice), compaction parameters b1 = -0.1 per MPa (ice) and b2 = -0.11 (silicate). Ice porosity phi_i = max(0.60 exp(-0.1 P), 0.20), P in MPa; silicate phi_s = min(0.6 P^-0.11, 0.6). A constant 70 percent rock mass fraction matches 14 of 17 Kuiper Belt objects within 2 sigma; minimum formation time 4 Myr (rock fraction 0.5) or 5 Myr (0.7), from early Al-26. Predictions in one sentence each (p. 10): Makemake and 2007 OR10 similar to or slightly below Pluto (about 1800 kg/m^3), MU69 about 750. [V; the first three Table 1 values and the porosity parameters re-read by me] **Contradiction:** our ice compaction guess has a reference pressure P0 = 1e6 Pa, ten times lower than the paper's exp(-0.1 P) (P0 = 10 MPa) with a 0.2 floor; the paper wins. Our 0.08 scatter and 0.58-0.82 range are not in the paper. Eq. 3 (dphi/dt = phi P/eta) is printed with a positive sign, which would make porosity grow, and the file carries a stray note "error in equation 2b": check against Besserer et al. 2013 or the published version before coding Eq. 3. Fray and Schmitt (ice densities, vapour-pressure windows) is paid and not here.

### 2.11 Galaxy rows

**GR4 age-velocity relation, Milky Way only (R4).** Aumer and Binney 2009 Eq. 16 (p. 1292): sigma(t) = v10 ((t + tau1)/(10 Gyr + tau1))^beta, v10 the dispersion at 10 Gyr. Table 6 (p. 1296): [V]

| component | beta | tau1 (Gyr) | v10 (km/s) |
|---|---|---|---|
| sigma_U | 0.307 | 0.001 | 41.899 |
| sigma_V | 0.430 | 0.715 | 28.823 |
| sigma_W | 0.445 | 0.001 | 23.831 |
| 3D (Table 2, best fit) | 0.349 | 0.149 | 55.179 |

Star formation fitted with the best total fit: alpha -2.549, gamma 0.117 per Gyr, tau_max 12.602 Gyr. beta = 0.5 is excluded and "the very upper limit" (Table 3). Published and arXiv v1 numbers are identical; cite MNRAS 397, 1286. sigma_W/sigma_U runs from 0.41 at 1 Gyr to 0.57 at 10 Gyr, so the plan's constant 0.54 is an old-age value only. Table 5 row 3 prints sigma_final 61.82 km/s where the paper's own Eq. 16 gives 59.8 (a misprint, present in both versions). Aumer and Binney's v10_W of 23.8 compares with the plan's Sharma and Hayden 21.1 km/s at 10 Gyr, exponent 0.441; their beta_U 0.307 sits above the plan's radial exponent 0.251. **Nothing here covers other galaxies**: GR4 for other galaxies stays a plan derivation (Martinsson scaling with the 20 percent closure check).

**GR6 and GR3 thick disc, Milky Way.** Haywood et al. 2013: thick disc forms from about 13 to about 8 Gyr ago (4 to 5 Gyr long), the age scale uncertain by 1 to 1.5 Gyr; inner thin disc starts about 8 Gyr ago, outer thin disc 9 to 10 Gyr ago (so "thin-disc start" is two numbers); M_thick/M_thin at least 0.2, possibly up to about 1. [V] Aumer and Binney's disc age is 12.6 Gyr (tau_max). Gowanlock et al. 2011 gives disc density forms from which thick/thin ratios of 0.24 and 0.29 follow [R]; it has no star formation or metallicity formulas, no bulge, no halo. Pinna et al. 2019a covers one galaxy (FCC 170): star formation peaks at the 14 Gyr edge of its age grid (resolution about 3 Gyr), thick-disc history has a second bump near 10 Gyr and stops about 11 Gyr, and no thick/thin ratio is printed (the worker's digitized 0.45 is a bin-mass ratio, not a disc ratio; its galaxy mass is internally inconsistent). The plan's statements about FCC 153, FCC 177 and a "2 to 4 Gyr first phase" need Pinna 2019b (not uploaded) and Bensby 2014 (not uploaded) and stay unverified.

**GR15 occupation width (R15).** Graus et al. v1 preprint, Eq. 1: the occupation is an erf in peak velocity, not in mass: f(>Vpeak) = (1 + erf((Vpeak - V50)/(sqrt2 sigma)))/2 with V50 = 18 km/s and sigma = 2.5 km/s (alternative V50 = 10.5 km/s), sigma assumed from FIRE-2, not fitted. [V; re-read by me] Converting to mass with Mpeak proportional to Vpeak^3 anchored at 1e9 solar masses at 20 km/s gives widths of about 0.18 dex (V50 = 18) and 0.32 dex (V50 = 10.5) [R]. The conventional model's edge-to-edge span (25 to 10 km/s) is 1.19 dex against the 1.22 dex between the Benitez-Llambay and Frenk edges, so the plan's derived width is confirmed, not tightened. No reionization redshift is used, so GR12 is not served. The v1 numbers may differ from the published MNRAS 488, 4585.

**Knapen et al. book.** All 10 chapters present. Qualitative or second-hand for GR1, GR3, GR6, GR10, GX5; nothing for GR5, GR7 to GR9, GR11 to GR15 or the GX5 baseline. Two conflicts with our files: its Comeron 2014 statement (hot-to-cold mass ratio independent of galaxy mass) disagrees with the mass-dependent Comeron 2021 law the plan uses, and its Merritt 2016 halo-fraction scatter (1.01 dex, 8 galaxies) is about three times the Harmsen 0.33 dex in the plan. Minus signs are lost in the text layer; check any sign taken from it against the primary paper.

### 2.12 Ice-giant composition (not a listed gap)

Irwin et al. (author manuscript): H2S above the main Uranus cloud (1.2 to 3 bar) is 0.4 to 0.8 ppm at the cloud top, column (2 to 5) x 10^19 molecules per cm^2, polar cap not detected; deep S/N above 4.4 to 5.0 times solar; the cloud is probably H2S ice, not NH3. Deep CH4 4 percent and He:H2 = 0.131 are assumed inputs, not retrievals. Relative humidity depends on the temperature profile (113 percent with one profile, 16 percent with the Spitzer profile), so use the cloud-top mole fraction. Nothing in the plan carries a number to compare. [V; re-read by me]

## 3. Where the papers contradict our files

Standing rule applied: the verified source wins; both values stated.

| Our file or plan | Paper | What differs |
|---|---|---|
| Ice compaction P0 = 1e6 Pa (parts file) | Bierson and Nimmo | Paper: exp(-0.1 P) per MPa, floor 0.2; ten times different. Scatter 0.08 and range 0.58-0.82 are not in the paper |
| "Palme and O'Neill 2014" Th 0.0849, U 0.0229 ppm | Palme and O'Neill 2003 | Th 0.0834, U 0.0218 ppm in this edition (+1.8 and +5.0 percent). 2014 numbers cannot be confirmed from this file |
| Mantle EOS rows = Stixrude Part II (assumed) | Part II Table A1 | Our rows are the 2024 Part III set; q0, theta0 and post-perovskite K0 differ |
| Labels "IPCC 1990" on 5.35 and 0.12 | Myhre Table 3 | IPCC 1990 values were 6.3 and 0.14; 5.35 and 0.12 are Myhre's |
| "5 to 9 percent" generic Rayleigh error | Sneep and Ubachs | 2.5-4 percent for Ar and N2, 7 for CO2, 10-18 for CO, CH4, N2O |
| Earth tau_R computed 0.0969 | Bucholtz Table 4 | Published 0.09721 (0.3 percent higher) |
| Tabazadeh as source of acid vapour pressure | Tabazadeh | Gives water pressure only |
| Repo comet NH3 upper bound 3.6 percent | Mumma and Charnley | About 1.5 percent in Table 4 |
| Constant sigma_z/sigma_R 0.54 | Aumer and Binney | 0.41 at 1 Gyr rising to 0.57 at 10 Gyr |
| Plan expected Wolff tables of w0, n, k | Wolff | No tables; figures only |
| Gaillard FMQ-2 +-0.5 (own text) | Gaillard Fig. 2 | About FMQ-2.75, 0.7 dex spread at 0.01 bar |
| Code Mars point 16.4 g/cm^2, 237 mSv/yr | Hassler | Value matches; the measured column is 21 g/cm^2 |

Everything else compared agreed: all crust oxides, Th, U and K in the parts files match Rudnick and Gao exactly; Hassler confirms our 0.21 and 0.64; Bates and Bucholtz confirm the air numbers.

## 4. Printing faults in the papers (for whoever codes them)

- Sneep and Ubachs: Eq. 13 prefactor (CO2) is wrong by a factor 1000 as printed (1.1427e3 reproduces Table 2); Eq. 17 (CO) gives about half the paper's own Table 2 n-based value; Eq. 23 (O2) lacks the 10^12 and its constants are at 0 C while the N2 constants are at 15 C; Table 1 prints N2 delta as 0.974 (must be 0.0974) and sigma-bar in units 100 times off; the O2 cross section at 532 nm prints "4066 x 10^-27" (they recompute 4.66e-27 against 4.50 +-0.15e-27 measured; which was meant is unresolved). The working forms in `atmosphere-sources-check.md` already avoid most of these.
- Bates Table 1 prints F'_K = 1.047 at 0.500 micron where the cross section implies 1.049 (Bucholtz prints 1.049). Bucholtz Table 4, subarctic winter at 0.55 micron prints 9.761e-2, about 0.5 percent above its neighbours and the fit (a misprint, not stated as one).
- Zeleznik Eq. 12 ("ln x_i" should be ln x_j) and Eq. 11 (missing minus sign), as in 2.3.
- Tabazadeh Table 2 branches are discontinuous at a_w = 0.05 and 0.85.
- Myhre Table 2 total does not equal the sum of its rows.
- Aumer and Binney Table 5 row 3 sigma_final, as in 2.11.
- Bierson and Nimmo Eq. 3 sign and the "error in equation 2b" note, as in 2.10.
- Hassler: the "25 percent of GCR dose equivalent" SEP sentence holds only for absorbed dose (24 percent; 8 percent for dose equivalent); the mission total is 982 mSv from the tables against about 1.01 Sv in the text.
- Pahlevan Figure 5 top axis, as in 2.7.

## 5. Papers still missing

**Paid-only, from Boss's 12:11 PM EDT list, not uploaded, and what each gates** (Boss, 3:07 PM EDT 2026-10-18: he has already provided all the AD1 data he can; none of these is a request to him now, each gap stays a labelled stand-in until a measurement is found):

| Paper | Gates |
|---|---|
| Fray and Schmitt 2009 | BC12 (GEN.470, 497): temperature windows of each sublimation fit and the cold end; BC7 ice densities (GEN.442, 457), partly covered by Bierson's 920 |
| ICRP 103 (2007) full report | BC13 (GEN.436, 444): wR, wT, dose limits. The free extract does not hold them (printed pp. 61-78 and from 98 are missing) |
| Karkoschka 1994; Karkoschka and Tomasko 2010 and 2009 | AD1 (GEN.422, atmosphere colour): CH4 and NH3 visible absorption of the giants |
| Tomasko et al. 2008 | AD2 (GEN.423): Titan haze n, k and phase functions |
| Ayers, Gillett and Gras 1980 | AD7 (GEN.471, 491): pure-acid vapour pressure by weight percent. Water side is filled; acid side is partly filled |

**On the fetch list, not paid-only, and not uploaded** (Boss's message does not say whether he tried; this thread cannot tell): Bogaard et al. 1978 (AD4: SO2 and H2S depolarization, the main open Rayleigh item); JPL Evaluation 19 (AD8 GEN.496, 501); Stixrude and Lithgow-Bertelloni 2005 Part I (BC3); Palme and O'Neill 2014 (BC4); Lodders 2003 (BC1 condensation temperatures); Hirschmann 2012 and Gaillard et al. 2021/2022 (AT.a1, a2); Bensby et al. 2014 and Pinna et al. 2019b (GR6). Sharma and Hayden 2020 and Nadler 2020 (not uploaded, not on the list) are named by the readers as sources the plan leans on.

**Several items on the original list are now reachable without Boss**, per `papers-to-fetch-2026-10-17.md` section E (arXiv copies of Haywood, Aumer and Binney, Pinna, Graus; free publisher copies of Lodders and Stixrude Part I per OpenAlex, untested through our proxy).

## 6. Uploads that were not on the fetch list

Summa 2021 (cross-check for CO and N2 only), Irwin (Uranus H2S, section 2.12), Curchod and Orr-Ewing 2024 (a review, no numbers for AD8), Pahlevan et al. (section 2.7), Gowanlock et al. 2011 (section 2.11) and the Knapen et al. book. The Curchod review prints 2019 for JPL Evaluation 19 where the plan says 2020. Gowanlock and Knapen serve none of the 15 requests directly.

## 7. Surprises and questions for the TODO thread

1. Twelve entries of the section 3 table carry a source label or value that these papers do not support or support differently; the gap plan's registry should record each as "paper read, differs".
2. The 2014 editions of Palme and O'Neill and of Rudnick and Gao: Rudnick's 2014 chapter is the 2003 text, so only Palme's second edition could still differ. Whether it is worth a fetch is the TODO thread's and Boss's call.
3. AD4 for SO2 and H2S is the largest hole this batch leaves in the colour work (Venus-like and volcanic worlds); Bogaard 1978 is the single paper that closes it.
4. Boss's four batches never carried the 1998 Bucholtz comment and reply (Appl. Opt. 37, 428 and 430), so whether they change Tables 3 or 5 is unchecked.
5. If Boss accepts digitized Wolff spectra for GEN.423, they should be labelled "digitized from Figures 12-14" in the data module, never cited as the paper's table.
