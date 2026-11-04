# Bulk composition chain: sources check (GEN.425)

**Status (2026-10-17):** the result of GEN.425 (done), copied from the read-only thread's report so that every lane has it in the repository. It feeds GEN.426 to GEN.430, DB.69 and the surveys GEN.432 and GEN.433; the lane that builds each item reads its section first. Where it says a source was not reached, the input is open and is Boss's to research (section 5). The check scripts and parsed tables stay in /mnt/project-files/reports/bulk-composition-sources-check-2026-10-17/ and are not committed. Times in it are as the lane's clock reported.

Date: 2026-10-17 (task started about 11:53 PM EDT on 10-16 as the wake clock reported; lane clocks drift). Read-only research, no repo change. Chain: goal issue #1949, plan note `docs/design/composition/atmosphere-colour.md` section 5, trace `composition-density-link-2026-10-17.md`. Feeds GEN.426 (design), GEN.427 (materials table), GEN.428 (fractions drawn), GEN.429 (density, mass, radius derived).

Check scripts and the parsed tables are in `/mnt/project-files/reports/bulk-composition-sources-check-2026-10-17/` (the scripts read the `pdftotext` output of the named arXiv PDF from the working directory; they are not committed to the repo).

## 0. How far the check could go (read this first)

Labels: **READ** = the paper's own text was read (arXiv PDF or the ADS scan). **COMPUTED** = re-derived here with Python from READ inputs. **2ND-HAND** = seen only in a citing paper or in open-source code. **NOT REACHED** = no copy could be read.

Network, as it stood when this was written (Boss widened access just before 12:02 AM EDT; the first half of the work ran before that):

| Host | Result |
|---|---|
| `export.arxiv.org` PDFs | work (used for every paper below); the `/api/query` endpoint timed out at 60 s on all seven tries |
| `arxiv.org` (abs and pdf) | **now reachable** (200 for `/abs/1507.08081` and `/pdf/0707.2895`) |
| `doi.org` | answers, but only redirects to `adsabs.harvard.edu`, which does not resolve from this container |
| `articles.adsabs.harvard.edu` | **now reachable for old scans**: the first try for [R3932] timed out, the second returned the 29-page PDF with a text layer. [R2726] (`2003ApJ...591.1220L`) returned 403 (no free scan) |
| `ui.adsabs.harvard.edu` link gateway | redirect loop, nothing returned |
| `iopscience.iop.org` | returns a Radware bot-check page instead of the paper; I did not try to get around it |
| `cdsarc.cds.unistra.fr` | the server answers but refuses: 403 "unable to read htaccess" on `/ftp/J/ApJS/99/713/`, 500 on `viz-bin/cat`. `cdsarc.u-strasbg.fr` does not resolve |
| `ntrs.nasa.gov` (search API) | answered 200; nothing was needed from it |
| `water.lsbu.ac.uk` (Chaplin ice page) | no answer |
| `pypi.org`, `files.pythonhosted.org` | work (burnman 2.1.0 source used as an independent cross-check) |

A single command that probed a long list of other hosts at once was refused by the session's safety classifier, so I tested only hosts I needed. **Not tested:** ScienceDirect, OUP, Wikipedia, Zenodo, NIST, PMC, Crossref, Semantic Scholar, github.com.

## 1. Verdict and counts

12 items were on the list (the six papers named, plus equations of state, SCvH, the three stellar-chemistry papers, and the water fraction against the snow line).

| Result | Count | Items |
|---|---|---|
| READ and verified | 8 | [R1429, R2741, R3992, R4336, R4880, R4881, R4882], equation-of-state constants (Seager Table 1, cross-checked against two independent sources) |
| Verified with a limit that matters | 3 | [R3916] (stoichiometry reproduced, but methane column off by 15 percent and it needs stellar O, Mg, Si the project does not store); [R1146] (an inversion method, no closed-form law, so GEN.429 cannot take an equation from it); [R3932] (original read now, but the numerical tables were not obtained) |
| Assembled, not one law | 1 | Water fraction against the snow line: [R2318, R2865, R3193, R4880], plus a stoichiometric estimate computed here |

Also read: [R0709] and the abstract of [R0708] for what replaced SCvH. **Not reached:** [R2726], the SCvH numerical tables, Chaplin's low-pressure ice data.

**Corrections to what the notes or the papers say (9):**

1. [R3992: Table 3] labels the 8300 kg/m3 iron row "Fe(alpha)"; Table 1 gives 8.30 g/cm3 for Fe(epsilon) and 7.86 for Fe(alpha). The row is the epsilon phase.
2. [R2741]: the text says the envelope radius goes as `t^0.11` (solar) and `t^0.18` (enhanced opacity) with no sign; Eq. 3 correctly has exponent **-0.18** (planets shrink with age).
3. [R1429] pure water ice is **12 percent larger** in radius than [R4882] pure water (1.58 against 1.41 R_earth at 1 M_earth), because it uses ANEOS "warm ice"; [R4880] says older water equations of state are wrong at a few hundred GPa. Use the newer one for water.
4. [R4882]'s "100% Fe" column is not pure iron. It uses the PREM outer-core fit (rho0 7.05, K0 201 GPa), which carries Earth's light elements, so it is 4 to 8 percent larger than pure epsilon-iron [R1429, R3992].
5. [R4882]'s closed form is calibrated on Earth's seismic model (PREM). It is valid, but it is not "nothing copied from the solar system"; the constants 1.07 and 0.21 are fitted to a table built from Earth's mantle and core. See section 3.
6. [R4882]'s Table 2 water columns ("25%h2o", "50%h2o") do not say what the rest of the planet is made of, in the text I read. Do not use the table for water mixtures without that.
7. [R3992] MgO bulk modulus 177 GPa differs 9 percent from [R4232] (161.4 GPa, via burnman). Silicate and iron constants agree within 1 percent in density (section 2.7).
8. Water has two Vinet/BME parameter sets in [R3992]: BME for ice VII (K0 23.7 GPa, K0' 4.15) and a Vinet fit to DFT for ice VII, VIII and X (K0 14.3771 GPa, K0' 6.57972, rho0 1460). They are not interchangeable; the second covers the high-pressure range.
9. [R3916]'s Table 1 methane column (0.37 percent for the Sun) is not reproduced by its own stoichiometry (0.32 percent computed), everything else is. Cause not found.

## 2. Source by source

### 2.1 [R3992] (READ, verified)

**Law.** Spherical hydrostatic structure `dm/dr = 4 pi r^2 rho`, `dP/dr = -G m rho / r^2` with a cold equation of state per material; Vinet or fourth-order Birch-Murnaghan (BME) fitted to experiment up to roughly 200 GPa, Thomas-Fermi-Dirac (TFD, Salpeter and Zapolsky) above. Outer boundary `P(R) = 0`. Temperature matters little: a 3.5 percent density error is a 1.2 percent radius error (`R ~ rho^-1/3`).

**Constants (Table 1, copied as printed):**

| Material | K0 (GPa) | K0' | rho0 (g/cm3) | Fit |
|---|---|---|---|---|
| Fe (alpha) | 162.5 +/- 5 | 5.5 +/- 0.8 | 7.86 | BME |
| Fe (epsilon) | 156.2 +/- 1.8 | 6.08 +/- 0.12 | 8.30 | Vinet |
| FeS | 35 +/- 4 | 5 +/- 2 | 4.77 | BME |
| H2O ice VII | 23.7 +/- 0.9 | 4.15 +/- 0.07 | 1.46 | BME |
| MgO | 177.0 +/- 4 | 4.0 +/- 0.1 | 3.56 | BME |
| MgSiO3 enstatite | 125 | 5 (assumed) | 3.22 | BME |
| MgSiO3 perovskite | 247 +/- 4 | 3.97 | 4.10 | BME4 (K0'' = -0.016 /GPa) |
| (Mg0.88,Fe0.12)SiO3 perovskite | 266 +/- 6 | 3.9 +/- 0.04 | 4.26 | BME |
| C (graphite) | 33.8 +/- 3 | 8.9 +/- 1.0 | 2.25 | BME |
| SiC | 227 +/- 3 | 4.1 +/- 0.1 | 3.22 | BME |

**Modified polytrope (Table 3): `rho(P) = rho0 + c P^n`, P in Pa, valid below 1e16 Pa.** Fe 8300, 0.00349, 0.528; MgSiO3 perovskite 4100, 0.00161, 0.541; (Mg,Fe)SiO3 4260, 0.00127, 0.549; H2O 1460, 0.00311, 0.513; graphite 2250, 0.00350, 0.514; SiC 3220, 0.00172, 0.537. Fit error 2 to 5 percent for P < 5e9 and P > 3e13 Pa, up to 12 percent in between.

**Generic mass-radius form (Eq. 23):** `log10 Rs = k1 + (1/3) log10 Ms - k2 Ms^k3`, `k1 = -0.20945, k2 = 0.0804, k3 = 0.394`, good to 1 percent for Ms < 4; the paper says it holds up to about 20 M_earth. Effective scales per composition (Table 4, `R = r1 Rs`, `M = m1 Ms`): Earth-like Fe 0.325 / MgSiO3 0.675: m1 6.41 M_earth, r1 3.19 R_earth; Fe 4.34, 2.23; MgSiO3 7.38, 3.58; H2O 5.52 and 4.43 (from the polytropic rows).

**Range.** Solid planets, cold, up to about 20 M_earth for the generic law; the numerical curves go to 10 M_Jupiter but electron degeneracy then flattens the radius. The paper's proposal: a "super Earth" is a planet without a substantial gas envelope.

**Verdict: verified**, one label slip (correction 1). The constants agree with later sources (section 2.7).

### 2.2 [R4881] (READ, verified)

Three-layer grid, 0.1 to 100 M_earth: Fe (epsilon, [R3992] Vinet to 2.09e4 GPa then TFD), MgSiO3 (perovskite and post-perovskite from [R0650], post-perovskite at 122 GPa, about 1.5 percent denser, dissociation at 0.90 and 2.10 TPa from [R4465], TFD above 4.9 TPa), H2O (Chaplin ices up to 2.216 GPa, then Frank, Fei, Hu 2004 for ice [R1457] DFT 37.4 GPa to 8.893 TPa, TFD above). Surface pressure 1 bar.

One formula worth lifting (Frank, Fei and Hu 2004 as printed, ice VII, 2.216 to 37.4 GPa): `rho (mol/cm3) = 0.0805 + 0.0229 (1 - exp(-0.0743 p)) + 0.1573 (1 - exp(-0.0061 p))`, p in GPa.

There is **no closed-form radius law in this paper**; its product is the grid (Table 1 for six 2-layer curves, Table 3 for ternary diagrams) and contours of central pressure, core mass fraction and core radius fraction. Water: the paper follows ice along the melting curve; superionic phase included.

**Verdict: verified as a method and as the source of the three-layer equations of state; no equation for GEN.429 comes out of it directly.**

### 2.3 [R4882] (READ, verified)

**Law (Eq. 2):** `R/R_earth = (1.07 - 0.21 CMF) (M/M_earth)^(1/3.7)`, **for 0 <= CMF <= 0.4 and 1 <= M <= 8 M_earth**; agrees with the paper's own curves within about 0.01 R_earth, and inverted (Eq. 3, `CMF = (1/0.21)(1.07 - (R/R_e)/(M/M_e)^(1/3.7))`) gives CMF to about 0.02. Error propagation (Eq. 4): `|dCMF| ~ 5 sqrt((dr/r)^2 + ((1/3.7)(dm/m))^2)`.

**Equations of state (second-order Birch-Murnaghan, Eq. 1: `P = 1.5 K0 [(rho/rho0)^(7/3) - (rho/rho0)^(5/3)]`), fitted to PREM:** lower mantle rho0 3.98 g/cm3, K0 206 GPa; outer core rho0 7.05, K0 201; inner core rho0 7.85, K0 255; each good to about 1 percent. Mantle: linear PREM table to 23.83 GPa, BM2 to 3.5 TPa, TFD (MgSiO3) above. Core: BM2 to 12 TPa, TFD Fe above. Earth's CMF is 0.325 in the paper; TFD for MgO, SiO2, MgSiO3 and Mg2SiO4 are almost identical (average A = 20, Z = 10), so Mg/Si does not matter at high pressure.

**Table 2** gives radii for 100% Fe, 50, 30, 25, 20% Fe, rock, 25% H2O, 50% H2O, 100% H2O at 0.125 to 32 M_earth (copied into `compare_mass_radius_laws.py`).

**Checks (COMPUTED):** the closed form reproduces the paper's own 30% Fe column to 0.002, 0.006, 0.005 and 0.016 R_earth at 1, 2, 4 and 8 M_earth. The BM2 constants give back PREM at the core-mantle and inner-core boundaries (mantle 5.569 g/cm3 at 135.8 GPa, outer core 9.914 g/cm3 at 135.8 GPa and 12.17 g/cm3 at 328.85 GPa; the PREM values themselves were not read here, the paper only states 1 percent).

**Verdict: verified.** Caveats: calibrated on Earth's interior (correction 5); the inner core and light elements in the iron are baked in; mass range 1 to 8 only. [R4880] applies the same `M^(1/3.7)` scaling to water worlds (section 2.4).

### 2.4 [R4880] (READ, verified)

Argues that the 2 to 4 R_earth planets are partly water worlds. What is usable:

- `R = M^(1/3.7)` (Earth-like, R_earth and M_earth units) for rocky planets, and **`R = f M^(1/3.7)` with `f = 1 + 0.55 x - 0.14 x^2`**, x the ice (water) mass fraction; **`f = 1.41` for 100 percent water** (check: 1 + 0.55 - 0.14 = 1.41). Used there for 1 to 20 M_earth.
- Rock and ice grow together beyond the snow line: the rock to ice ratio "does not strongly depend on host star metallicity"; x between 1/2 and 2/3 is assumed for the icy planets.
- Radius reference: **1 millibar** surface; cosmic H2-He is 75 percent H2 and 25 percent He by mass; water in the interior stays solid along the melting curve; envelope curves are isothermal at 300, 500, 700, 1000 K.
- "A gas envelope of at most a few percent by mass" is enough for any planet between 2 and 4 R_earth.

**Checks (COMPUTED):** `f(x)` times the "rock" column of [R4882: Table 2] gives 25 percent water 1.204 against the table's 1.182 (1.9 percent) and 50 percent water 1.323 against 1.270 (4 percent) at 1 M_earth. So `f` is good to about 2 to 4 percent if the base is pure rock. [R1429] eq. 7 (section 2.6) agrees with [R4882] to about 5 percent for mixtures.

**Verdict: verified.** The 1:2 ratio of rocky to icy cores and the claim that they are water worlds are the paper's population argument, not a law to copy; the paper itself says gas dwarfs against water worlds cannot be told apart from mass and radius alone.

### 2.5 [R2741] (READ, verified)

**Core (Eq. 1):** `R_core = (M_core/M_earth)^0.25 R_earth` for an Earth-like rock/iron core (2:1 rock to iron, ANEOS olivine and SESAME 2140 iron) to within about 2 percent; about 10 percent if the iron fraction varies.

**Envelope (Eq. 3, enhanced-opacity fit, defaults 5 M_earth, 100 F_earth, 5 percent, 5 Gyr):**
`R_env = 2.06 R_earth (M/M_earth)^-0.21 (f_env/5%)^0.59 (F/F_earth)^0.044 (age/5 Gyr)^-0.18`, total `R = R_core + R_env + R_atm`, `R_atm ~ 9 H` (about 9 scale heights, usually near 0.1 R_earth, from 100 bar to the 20 mbar level where the transit radius is defined).

**Domain.** 1 to 20 M_earth; envelope mass fraction 0.01 to 20 percent; flux 0.1 to 1000 F_earth (tables at 0.1, 10, 1000); age 100 Myr to 10 Gyr; solar metallicity and 50x enhanced-opacity atmospheres. Hydrogen/helium by [R3932]. Radius defined at **20 mbar**. Rocky/gaseous divide: the paper suggests **1.75 R_earth**.

**Checks (COMPUTED, 24 rows x 11 columns x 6 tables parsed into `lopez_fortney_2014_tables2-7.csv`):** with `R_atm` left out, Eq. 3 plus Eq. 1 against the tables for envelope fractions of 0.5 percent and above: at 1 Gyr enhanced opacity mean +3.1 percent, RMS 6.3 percent, range -11 to +15; at 10 Gyr RMS 8.3 percent, range -22 to +35; solar metallicity at 1 Gyr bias +8 percent, RMS 11 percent. That matches the paper's own "within about 0.1 dex" and its warning that the fit is a rough guide.

**Verdict: verified** (correction 2). For GEN.429 interpolate the tables (CSV above) instead of using Eq. 3; the equation is a sanity check.

### 2.6 [R1429] (READ, verified)

Equations of state: iron SESAME 2140, olivine and water ANEOS, hydrogen and helium [R3932] (Y = 0.28). Warm ice adds a thermal pressure `P_T = 3.59e-5 rho T` (P in Mbar, rho in g/cm3, T in K) on a Uranus/Neptune adiabat. Rock and iron are cold (no thermal term; accurate to 1 to 2 percent). Giant cores: 50/50 ice/rock by mass.

**Closed forms (Eqs. 7 and 8, R in R_earth, M in M_earth, log base 10), fitted 0.01 to 100 M_earth:**

- Ice/rock: `R = (0.0912 imf + 0.1603)(log M)^2 + (0.3330 imf + 0.7387) log M + (0.4639 imf + 1.1193)`, imf = ice mass fraction (1 pure ice, 0 pure rock).
- Rock/iron: `R = (0.0592 rmf + 0.0975)(log M)^2 + (0.2337 rmf + 0.4938) log M + (0.3102 rmf + 0.7932)`, rmf = rock mass fraction (1 pure rock, 0 pure iron).

Authors' accuracy claim: 2.5 percent average for ice/rock, 1.5 percent for rock/iron, better than 0.5 percent near 1 M_earth, up to 10 percent for ice/rock at 0.01 M_earth.

**Checks (COMPUTED against the paper's own Table 1, M <= 100):** rock/iron with rmf 0.67 mean 1.1 percent, max 3.7; rmf 0.5 mean 1.6, max 3.2; pure iron mean 2.5, max 4.1; pure rock via Eq. 8 mean 2.3, max 10. Ice/rock: pure ice mean 3.4, max 17; 50/50 mean 3.0, max 10; pure rock via Eq. 7 mean 3.5, max 13. So the rock/iron equation meets its claim; the ice/rock one is worse than claimed at the ends (the claim "2.5 percent" is the average over its whole range, my numbers are for the ten tabulated masses, so the end errors dominate).

Radius reference is the 1 bar level; hot giants need the model atmosphere tables (Tables 2 and 3 there: 0.02 to 9.5 AU, 17 to 3594 M_earth, core 0 to 100 M_earth; not needed for the solids chain).

**Verdict: verified**, with the water warning in correction 3.

### 2.7 Equation-of-state constants: Vinet, Birch-Murnaghan (verified, cross-checked)

Seager Table 1 (section 2.1) is the verified source. Independent checks, COMPUTED:

- Epsilon-iron: Seager Vinet (8.30, 156.2, 6.08) against the calibrant in burnman 2.1.0 (V0 6.81706 cm3/mol so rho0 8.192, K0 148, K0' 5.86): density at 50, 100, 200, 330 GPa 9.987/9.947, 10.995/10.999, 12.41/12.48, 13.744/13.886 g/cm3, within 1.0 percent.
- MgSiO3 perovskite: Seager (4.10, 247, 3.97; third-order BM here) against [R4232] as shipped in burnman (V0 24.445 cm3/mol so rho0 4.107, K0 250.5, K0' 4.14): 25, 50, 100, 135 GPa density ratio 1.000, 1.001, 1.004, 1.006. [R4882]'s lower-mantle BM2 (3.98, 206) is within 1.6 percent of it at 25 GPa and 1.1 percent at 135 GPa.
- Periclase: [R4232] rho0 3.584 (agrees with 3.56) but K0 161.4 GPa against 177 (correction 7).
- Two-layer integration (COMPUTED, `two_layer_structure_check.py`): Fe Vinet (Seager) over MgSiO3 BM3 (Seager), CMF 0.3, no phase changes, no TFD: radii 0.792, 0.978, 1.198, 1.453, 1.744 R_earth at 0.5, 1, 2, 4, 8 M_earth. [R4882: Table 2] (30% Fe): 0.823, 1.005, 1.22, 1.47, 1.75 at the same masses. So **the Seager constants alone reproduce a published table to 1 to 4 percent** (a little low because the BM and Vinet forms are extrapolated cold and the Earth-like upper mantle is missing). Use BM only up to about 2.8e13 Pa; Seager Fig. 1 shows the MgSiO3 BME curve turns unphysical there and TFD must take over.

**Verdict: verified.** What the table lacks: a verified density for **low-pressure ice phases** (ice Ih and the others below 2.2 GPa). Seager's 1.46 g/cm3 is ice VII at zero pressure (high-pressure phase), the wrong number for a small icy body or moon. Ice Ih was not found in any source read (section 5).

### 2.8 Saumon, Chabrier, ApJS 99, 713 (original READ from the ADS scan; tables NOT obtained)

From the paper itself: hydrogen and helium equations of state for low-mass stars (M < 1 M_sun), brown dwarfs and giant planets, covering **2.10 < log T(K) < 7.06 and 4 < log P (dyn/cm2) < 19**; free-energy minimisation; **mixtures by the additive volume rule plus an ideal entropy of mixing** (so the density of a mixture with helium mass fraction Y follows from the two pure tables); a hypothetical plasma phase transition for hydrogen, with tables both with it and with it interpolated away; tables published on the AAS CD-ROM Series Vol. 5 and by anonymous FTP; the paper's own figure uses Y = 0.25. Radiation is not included.

What replaced it: [R0709] keep SCvH only for ρ <= 0.05 g/cm3 and add ab initio data in between; they state SCvH underestimates H2 pressure dissociation, is too stiff at low pressure, and overestimates the internal energy in the dissociation regime. Sizes quoted there: SCvH pressure about 10 percent above the simulations at 2 g/cm3 along an isentrope, and departures "as much as 25%" in one density-temperature regime (which quantity that 25 percent refers to is not clear in the extracted text). [R0708] ([R0708], abstract only) add H-He interactions, which matter at low temperature and high density (giants and brown dwarfs).

**Verdict: verified as a source, but the numbers are not in hand.** The CDS catalogue J/ApJS/99/713 refused every request from here (section 0). For a planet generator the envelope radius comes from Lopez and Fortney or [R1429] tables, which already used SCvH; the generator does not need to integrate SCvH. Only a lane that integrates the envelope itself needs the tables.

### 2.9 [R1146] (READ, verified as a method)

A Bayesian inversion, not a forward law. Assumptions to carry over: bulk **Fe/Si and Mg/Si mass ratios of the planet equal the star's** (cites [R0455, R2180, R4335], because Fe, Mg, Si condense together at small distances); a pure-iron core; a mantle in the Na2O-CaO-FeO-MgO-Al2O3-SiO2 system (more than 98 percent of Earth's mantle) with Na2O 2, CaO 2, Al2O3 4.8 weight percent fixed at chondritic values; mineral equations of state from [R4231]/2011 by Gibbs minimisation; core by [R0332] (Mie-Gruneisen-Debye with rho0 8.334, q 0.489, gamma0 2.434).

Relation worth keeping (Eq. 1): `Fe/Si_bulk = Fe/Si_mantle + M_core / M_Si,mantle`, so a larger core must come with less iron (FeO) in the mantle; Earth's mantle values in the paper's Table 2: Fe/Si 0.17, Mg/Si 0.83, Si 22 weight percent, core radius 3400 km. Prior ranges (Table 1): mantle Si 14 to 30 weight percent, mantle Fe/Si up to 2.23, Mg/Si 0.89 +/- 0.08 ([R2726] Sun: Fe/Si 1.69 +/- 0.18). Tested to 10 M_earth and 1.7 R_earth, plus Kepler-36b. Conclusions: mass and radius fix the core size; stellar Fe, Si, Mg reduce the degeneracy.

**Verdict: verified.** It gives no equation for GEN.429. The usable consequence is the bulk constraint, turned into a core mass fraction in section 3.3.

### 2.10 [R4336] (READ, verified)

Gibbs-energy minimisation (HSC Chemistry 7.1) of nebular condensation inside the Alibert planet-formation model, for 18 stars covering **Fe/Si 1.34 to 2.81, Mg/Si 0.66 to 1.85, C/O 0.12 to 0.75 (mass ratios)**, solar-mass and solar-luminosity hosts. Result: **Mg/Si and Fe/Si in rocky, icy and giant planets fall on the 1:1 line with the star**; deviations only for material that formed where the disc is too hot to condense Mg, Si, Fe (inside about 0.3 AU, roughly 1000 K, for the disc they show), C/O shows only a very weak link.

Thresholds as printed (mass ratios, pressure about 1e-4 bar): Mg/Si below 0.86 puts Mg in pyroxene (MgSiO3) with leftover Si in other silicates; 0.86 to 1.73 splits between pyroxene and olivine; above 1.73 all Si goes to olivine and leftover Mg to oxides; C/O above 0.6 makes Si carbides. **Checked (COMPUTED):** 0.86 and 1.73 are exactly the molar ratios 1 and 2 (24.305/28.0855 = 0.865, twice that 1.731); C/O 0.6 is molar 0.8 (0.8 x 12.011/15.999 = 0.601).

**Verdict: verified.** Limits: solar-mass stars only, one family of discs, solar-system-type formation.

### 2.11 [R3916], A&A 580, L13 (READ; stoichiometry reproduced)

The model: heavy elements as H2O, CH4, Fe, MgSiO3, Mg2SiO4 in equilibrium condensation; with stellar abundances N_X: `N_O = N_H2O + 3 N_MgSiO3 + 4 N_Mg2SiO4`, `N_Mg = N_MgSiO3 + 2 N_Mg2SiO4`, `N_Si = N_MgSiO3 + N_Mg2SiO4` (valid while N_Mg > N_Si). The iron fraction among refractories is `f_iron = m_Fe / (m_Fe + m_MgSiO3 + m_Mg2SiO4)`. The paper's results: C-7 31.6, K-93 34.7, K-10 27.5, Sun 33.2 percent ([R2726] condensation, 38.0 percent). The paper compares the 17 percent iron from mass and radius [R1156] with 28 to 35 percent from the stars and calls the agreement good "on a 10 to 20 percent level" (their words).

**Check (COMPUTED, `check_santos_2015_stoichiometry.py`):** with Asplund-2009-style solar number abundances (log eps C 8.43, O 8.69, Mg 7.60, Si 7.51, Fe 7.50, He 10.93; **typed from memory, not read from a source here**), the script gives MgSiO3 0.183, Mg2SiO4 0.077, H2O 0.509, Fe 0.130 mass percent, Z 1.216 and f_iron **33.2**, matching the paper's Sun column (0.19, 0.08, 0.50, 0.13, 1.26, 33.2) except methane (0.32 against 0.37; correction 9). Because 33.2 is reproduced exactly the abundances are consistent with the paper's, but a lane must cite its own Asplund source.

**Verdict: verified with a limit.** Needs oxygen, magnesium, silicon and carbon as well as iron from the host; and the project stores only [Fe/H] (section 5).

### 2.12 Water fraction against the snow line (assembled)

There is no single law. What the sources give:

| Piece | Source | Value |
|---|---|---|
| Ice/rock mass ratio beyond the water ice line | [R2865] (READ) | **1 +/- 0.5** (ice 45 +/- 15 weight percent of condensed mass) without refractory organics; **0.35 +/- 0.25** (30 +/- 10 percent) with them; across model variants 0.1 to 1.5. These are **all** ices (H2O, CO, CO2, CH3OH, NH3), not water alone; strong gradient near the water line |
| Rock and ice fractions for sub-Neptune cores | [R4880] (READ) | x = 1/2 to 2/3 ice |
| Stoichiometric water share (COMPUTED here from [R3916]'s equations and the Asplund-style Sun) | [R3916] | water 0.509, refractories 0.390 mass percent of the nebula, so **water / (water + refractories) = 0.566**; with methane counted as ice 0.68 |
| Snow-line position, dust and accretion | [R3193] (READ) | `R_SL = [3 mu m_p (G M*)^(3/2) Mdot^2 kappa_R / (128 pi^2 k_B sigma_SB T_ice f alpha)]^(2/9)` with `T_ice` 160 K (150 to 170 K), steady viscous disc. For 1.0, 0.8, 0.6, 0.4, 0.2 M_sun: small grains 3.3, 2.6, 1.9, 1.2, 0.5 AU (1-sigma 2.2 to 5.0 for 1.0 M_sun); larger grains 1.6, 1.2, 0.9, 0.6, 0.2 AU. Observed scatter in accretion rate (0.4 dex) gives about 0.2 dex in snow-line distance. For comparison the same table lists luminosity-based main-sequence values 2.5, 1.23, 0.61, 0.32, 0.16 AU |
| Snow line moves in time | [R2318] (READ in part) | located where T_mid = **170 K**; moves inward as accretion falls; the old fixed `2.7 (M*/M_sun) AU` line is only a main-sequence reference |

**Computed derivation of the project's constant.** The equilibrium temperature of a blackbody at 1 AU from 1 L_sun is `T = (L / (16 pi sigma a^2))^(1/4) = 278.3 K`; so `a_snow = (278.3 / T_ice)^2 sqrt(L/L_sun)` AU: **2.68 AU at 170 K**, 3.03 at 160 K, 3.44 at 150 K. The fixed `SNOW_LINE_AU_AT_1_LSUN = 2.7` [R1854] was this 170 K value; GEN.458 removed it: `SNOW_LINE_TEMPERATURE_K = 170` now gives the blackbody centre (`formation.snow_line_au`, 2.68 AU for the Sun) and each star draws and stores its own line around it. Mulders and Kennedy and Kenyon show accretion heating and dust size move the true line by a factor of 2 to 3 around it, so a drawn scatter of about 0.2 dex is justified by Mulders; the optically thin formula is the unscattered centre.

**Verdict: assembled; usable.** The water mass fraction inside the snow line is not in any source read (condensation models give essentially zero; delivery is a separate research item, section 5).

## 3. What the two tables should use

### 3.1 GEN.427 materials table

Per material: zero-pressure density, Vinet or BM constants, valid pressure range, source. Use Seager Table 1 and the Table 3 polytrope for the cold, high-pressure regime (verified, cross-checked in 2.7) and take the following **values to the 3 listed significant figures**:

| Material | rho0 (g/cm3) | K0 (GPa) | K0' | Form | Use |
|---|---|---|---|---|---|
| Fe (epsilon), core | 8.30 | 156.2 | 6.08 | Vinet | any pressure above about 100 GPa; TFD above 2.09e4 GPa |
| Fe (alpha), small bodies | 7.86 | 162.5 | 5.5 | BME | low pressure |
| MgSiO3 perovskite, mantle | 4.10 | 247 | 3.97 | BM4 (K0'' -0.016 /GPa) or BM3 | above 27 GPa; post-perovskite 1.5 percent denser above 122 GPa |
| MgSiO3 enstatite | 3.22 | 125 | 5 (assumed) | BME | small bodies, upper mantle (low-pressure silicate) |
| (Mg0.88,Fe0.12)SiO3 | 4.26 | 266 | 3.9 | BME | oxidised mantle |
| MgO | 3.56 | 177 (161.4 per [R4232]) | 4.0 | BME | choose and cite one |
| FeS | 4.77 | 35 | 5 | BME | sulfide core option |
| H2O ice VII | 1.46 | 23.7 | 4.15 | BME | 2.2 to 44 GPa (Seager); Frank et al. formula 2.2 to 37.4 GPa; DFT table above |
| graphite / SiC | 2.25 / 3.22 | 33.8 / 227 | 8.9 / 4.1 | BME | carbon-rich (C/O > 0.6 by mass) |
| H/He (cosmic) | SCvH or CMS19 tables | - | - | tables | envelopes; Y 0.25 to 0.28 in the sources |

**Not in any source read:** ice Ih and the low-pressure ices (Ih, III, V, VI) as bulk densities for small icy bodies and moons, porosity of small bodies, and CH4, NH3, CO2 ices. A zero-pressure density for them cannot be written down from the sources here; the "uniform 3.93 to 5.51 g/cm3 for every class" problem in the trace (class D "frozen hydrocarbons and ice") stays unsolved for low-mass bodies until that is researched.

### 3.2 GEN.429 density, mass, radius

Recommended order, by body mass M (M_earth), composition by mass fractions (Fe core CMF, silicate, water x, envelope f_env):

1. **Small bodies (below about 0.01 M_earth), no sources read for the regime:** pressure is negligible, so density is the harmonic mean of the zero-pressure densities `rho = 1 / sum(w_i / rho0_i)` (volume additivity; not from a paper; ignores porosity). This is a derivation, flag it as such. Seager shows the equation of state reduces to `rho = rho0` for scaled mass below 1 (Ms < 1, i.e. M below m1; m1 is 5.80 M_earth for the iron modified polytrope), so compression is small there.
2. **Solid planets 0.01 to about 100 M_earth (sources valid to 100):** integrate the two- or three-layer structure (Zeng's mass-as-coordinate form: `dr/dm = 1/(4 pi rho r^2)`, `dP/dm = -G m /(4 pi r^4)`) with the constants above, or interpolate a grid built once with them. That integration reproduces the published [R4882] table to 1 to 4 percent at 0.5 to 8 M_earth (2.7). It makes the composition the cause of the density, which is what Boss asked. The closed forms are checks and fallbacks:
   - rock plus iron: [R4882: Eq. 2] (1 to 8 M_earth, CMF 0 to 0.4, within 0.01 R_earth) or [R1429: Eq. 8] (0.01 to 100 M_earth, 1.5 percent);
   - rock plus water: [R4880] `f(x)` on the rock radius (2 to 4 percent) or [R1429: Eq. 7] (5 percent for mixtures against [R4882]; pure ice 12 percent big);
   - envelope on top: [R2741] tables (1 to 20 M_earth, 0.01 to 20 percent, radius at 20 mbar).
   Agreement of the closed forms with each other (COMPUTED, Earth-like CMF 0.325, R_earth): at 1, 2, 4, 8 M_earth [R4882] gives 1.002, 1.208, 1.457, 1.757; [R1429] eq. 8 gives 1.003, 1.211, 1.445, 1.703; Seager generic gives 0.970, 1.188, 1.443, 1.733; Lopez core `M^0.25` gives 1.000, 1.189, 1.414, 1.682. So **3 to 4 percent spread in radius, 10 to 12 percent in density**, is the honest uncertainty for rocky bodies from the literature.
3. **Over 20 M_earth with a gas envelope:** use the [R1429] giant tables (core 0 to 100 M_earth, 17 to 3594 M_earth, 0.02 to 9.5 AU, ages) rather than a closed form; the project already has a giant mass-radius relation (GEN.34) that this chain does not need to touch.
4. **Radius reference level differs:** Seager 0 pressure, [R4881] 1 bar, [R1429] 1 bar, Lopez and Fortney 20 mbar, [R4880] 1 millibar. The code's radius level has to be fixed once (1 bar is the common one) and the envelope tables corrected by the scale-height term.

### 3.3 GEN.428 composition fractions: what can be derived from the star

Core mass fraction from stellar Fe/Si and Mg/Si (mass ratios), assuming **all iron in the core and Mg and Si as oxides** (O/Si molar = 2 + Mg/Si molar for any Mg/Si, which follows from MgSiO3, Mg2SiO4, SiO2 and MgO): with silicate mass per gram of Si `1 + Mg/Si + (2 + r) 15.999/28.0855`, `r = (Mg/Si) 28.0855/24.305`, **CMF = (Fe/Si) / ((Fe/Si) + silicate mass per gram of Si)**.

Test (COMPUTED, `cmf_from_star_and_snow_line.py`): the Lodders-2003 Sun used by Dorn and Thiabaud (Fe/Si 1.69, Mg/Si 0.88) gives **CMF 0.320**; Earth's 0.325 ([R4882: Table 1]) is within 2 percent. Thiabaud's corners: Fe/Si 1.34, Mg/Si 0.66: 0.293; 2.81, 1.85: 0.351; 1.34, 1.85: 0.205; 2.81, 0.66: 0.465. So the stellar range gives CMF from about 0.2 to 0.47 (Mercury-like at the top is a relevant sanity check, not a target). It overstates the core when some iron sits as FeO in the mantle (Dorn: Earth's mantle Fe/Si 0.17) and ignores light elements in the core. The Santos route (f_iron 33.2 for the Sun) is the same idea with oxygen accounting.

## 4. Where it applies

| Law | Mass range | Composition range | Reference radius |
|---|---|---|---|
| [R3992] generic (Eq. 23) | up to ~20 M_earth (Ms < 4, extended to ~40 with Table 4 k) | homogeneous or layered Fe/silicate/water, carbon | P = 0 |
| [R4882: Eq. 2] | 1 to 8 M_earth | CMF 0 to 0.4, rock only | about 1 bar |
| [R4880] (f) | 1 to ~20 M_earth | x 0 to 1 water on rock | 1 mbar |
| [R1429] Eqs. 7 and 8 | 0.01 to 100 M_earth | ice/rock and rock/iron | 1 bar |
| [R2741] | 1 to 20 M_earth | f_env 0.01 to 20 percent, 0.1 to 1000 F_earth, 100 Myr to 10 Gyr | 20 mbar |
| [R4336] | solar-mass hosts | Fe/Si 1.34 to 2.81, Mg/Si 0.66 to 1.85, C/O 0.12 to 0.75 | n/a |
| [R3916] | rocky planets | Mg > Si, C as CH4 only | n/a |
| [R2865] | solar chemistry | beyond the water line, ice/rock 0.1 to 1.5 | n/a |

## 5. What Boss must research himself (or a lane must go and read)

1. **Low-pressure ice and silicate densities for small bodies and moons.** Ice Ih (and the other ices below 2.2 GPa), CH4, NH3 and CO2 ices, olivine and enstatite bulk densities with porosity. Not in any source read. Chaplin's water page was unreachable; [R0777] are the sources [R4881] name.
2. **Stellar Mg, Si, O and C against [Fe/H].** The project stores only [Fe/H] (`metallicity_fe_h`), but the core fraction and the water fraction need Fe/Si, Mg/Si, O and C/O. Needs a drawn relation with scatter; no value should be fixed.
3. **Solar abundance set.** The check above used values typed from memory (Asplund-style). The chain needs a cited abundance table to calibrate against; [R2726] (full text 403) and [R0169] were not read.
4. **Water inside the snow line** (rocky planets in the warm zone): condensation models say nearly zero; delivery of water by scattered or drifting bodies is a separate field. Not sourced; the existing separate `water_mass_fraction` draw needs a rule for how it relates.
5. **Iron oxidation and light elements.** FeO in the mantle and light elements in the core change CMF by ten percent or more [R1146, R4882]. No closed form exists; either accept "all iron in the core" with its stated bias or research a draw.
6. **SCvH or CMS19 tables** if any lane integrates an envelope itself. CDS refused from this container; the AAS CD-ROM and anonymous FTP routes are not available here.
7. **Hot, molten planets** (very close in): Seager says temperature matters little for solids; the water layers do depend on temperature (Fortney's warm ice, [R4880]'s 300 to 1000 K envelopes). Not needed for the first version.

## 6. Open points for the design (GEN.426)

- Which radius level the code uses (section 3.2, item 4).
- Whether the structure is integrated per planet, tabulated once, or a closed form with a stated error. Per planet costs a root-find per body; a grid in M, CMF, x and f_env is the usual answer.
- The closed forms disagree by 3 to 4 percent in radius (10 to 12 percent in density); a lane that derives density must say which source it follows and why.

## 7. Files

Scripts (all read the `pdftotext` text of the named arXiv paper from the working directory unless noted): `check_lopez_fortney_2014.py`, `check_fortney_2007_eq7_8.py`, `compare_mass_radius_laws.py`, `check_eos_constants.py` (self-contained), `check_santos_2015_stoichiometry.py` (self-contained), `cmf_from_star_and_snow_line.py` (self-contained), `two_layer_structure_check.py` (self-contained, needs scipy). Data: `lopez_fortney_2014_tables2-7.csv` (1584 rows: age, opacity, flux, mass, envelope fraction, radius). Papers read: [R3992], 1301.0818, 1512.08827, 1906.04253, 1311.0329, astro-ph/0612671, 1502.03605, 1507.01343, 1507.08081, 1407.7271, 1505.03516, 0710.1065, 1902.01852, 2107.04434 (abstract), and the ADS scan of ApJS 99, 713 (1995).
