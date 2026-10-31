# Group 6: evolved and frost atmospheres, trace gases, white dwarf and pulsar planets, host-star UV regime

Rows: j1 to j5, l1 to l3, m2 to m6, m11, p4, p6 of part 09. Labels (READ, 2ND-HAND, DATA, COMPUTED, SNIPPET, RECALLED, NOT REACHED) are those of `_worker-rules.md`. READ here means the arXiv PDF text was extracted and read by me (abstract and the passages quoted); where I read only the abstract it says "abstract".

## 1. Summary and verdicts

Most rows in this group are not probabilities. j1, j5, m2, l1, l3, m5, m6 and the abiotic part of m3 and m11 are outcomes of inventory (group 2), escape, photochemistry and speciation, so they should be computed, with the volatile inventory as the only random element. j2, j3, j4, m4 and the biotic parts of m3 and m11 have no observational ensemble at all and the only constraint is the Solar System, which breaks the NO FIXED CONTENT rule if used to tune a knob. The two rows that do have real data are p4 (tau_sat; N=55 stars in the early-M bin of Jackson 2012, N=219 volume-complete mid-to-late M dwarfs in Medina 2022) and, weakly, l3 and l1 (upper limits on transiting planets around 1148 K2 white dwarfs and 313 polluted white dwarfs, one microlensing and one transiting giant). The single `LU[0.162, 2] Gyr` for p4 mixes a mass dependence: 0.162 Gyr is the early-M X-ray value (READ), 2 Gyr is the fully convective value (READ, 2.3 +- 1.3 Gyr for 0.1 to 0.2 Msun). p6 can be replaced by a continuous percentile draw using the published Johnstone 2021 tracks (MIT-licensed code, integer percentiles 2 to 98, masses 0.1 to 1.25 Msun).

Verdict key: (a) published ensemble plus bracket suffices for a tunable seeded draw; (b) tunable only to a few stated constraints; (c) needs a source not in the public record. "Derive" means the quantity should be computed by physics, not drawn.

| Row | Quantity | Verdict | Draw or derive | Bracket-only seeded draw honest? Widest defensible bracket |
|---|---|---|---|---|
| j1 | Venus-like CO2-N2 | (b) | Derive: runaway gate (Kopparapu, READ) + water loss + inventory | Yes for the retention residual: [0, 1]; the 0.1 floor has no support (JWST TRAPPIST-1 b and c, N=2, allow bare rock) |
| j2 | Mars-like thin CO2 | (c) | Derive from inventory + escape (d table) | Yes: [0, 1]; Solar System only |
| j3 | Titan-like N2-CH4 moon | (c) | Derive: Jeans retention + CH4 source | Yes: [0, 1]; Solar System only, no exomoon ensemble |
| j4 | Pluto/Triton-like frost | (c) | Derive: vapour-pressure equilibrium + retention | Yes: [0, 1]; Solar System only (KBO retention model gives a physical bracket) |
| j5 | Abiotic O2 world | (b) | Derive: water loss, EUV, O2 sinks | Yes for the residual: [0, 1] (Harman 2018 says about 0 for wet HZ planets with trace N2; Luger and Barnes says hundreds of bar otherwise) |
| l1 | Silicate-vapour exosphere | (c) for the weight, derivable | Derive from P_eq vs exosphere criterion; threshold is a design constant | Threshold 1e-6 bar +- 0.5 dex has no data; Koester fraction is NOT the tuning constraint |
| l2 | Airless | complement of l1 | Consequence | n/a |
| l3 | Surviving giant at WD | (b) | Derive: engulfment gate (Veras 2016) x giant occurrence | Yes: [0, 0.5] (upper value is a prediction, 2ND-HAND); c[0, 0.1] is below the published prediction and above the transit upper limit for P < 27 d |
| m2 | Abiotic O3 | (b) | Derive from photochemistry grid; state-dependent flag | [0, 1] (conditional: about 0 for wet N2-bearing HZ planets, 1 for desiccated CO2 worlds) |
| m3 | N2O | (c) | Biotic part = design flag; abiotic SEP part = derived, model-only | [0, 1]; the c[0.05, 0.5] has no support |
| m4 | PH3 | (c) | Biotic design flag | [0, 0.05] is defensible; no confirmed detection on any rocky planet |
| m5 | HCl | (b) as a budget, (c) as a probability | Derive: halogen budget x speciation x solubility | Not a coin flip; see section 2.4 |
| m6 | HF | same as m5 | same | same |
| m11 | NO2, CH3Cl | (c) | CH3Cl = design flag; NO2 = derived from lightning rate | [0, 1] for CH3Cl; lightning rate is Solar-System-anchored |
| p4 | M dwarf saturated, tau_sat | (a) for early M, (b) for mid and late M | Derive from Rossby saturation + spin-down track, or tabulate tau_sat(M) | Replace LU[0.162, 2] Gyr with a mass-dependent value (section 2.5) |
| p6 | Rotation-track weights | (a) | Replace 3 weights by a percentile draw | Draw u in [0.02, 0.98]; 0.25/0.5/0.25 is defensible only as quartile bins |

Places where a verified source overrides a partly verified one (usability rule):
1. m2: part 09 says "O3 about 10 % of Earth column (Domagal-Goldman 2014, 2ND-HAND)". The paper's Table 3 (READ) gives abiotic O3 columns of 4.5e15 (Sun), 1.0e15 (K2V), 8.8e15 (M3.5V), 1.8e15 (M4V) and 1.3e18 (F2V) molecules per cm2 against 8.6e18 for present Earth: that is 0.05, 0.01, 0.10, 0.02 and 15 percent. Only the F-star case is near 10 percent. READ wins.
2. p4: the project stand-in "2 Gyr" for tau_sat is overridden for early M by Jackson 2012 Table 2 (READ, 162 Myr) and supported for 0.1 to 0.2 Msun by Medina 2022 (READ).
3. m4: "TEST 20 ppb Venus claim (Greaves 2020 READ)" is contradicted by two other READ papers (Villanueva 2021 upper limit < 1 ppbv at 3 sigma; Lincowski 2021 finds the line forms above 80 km and is fit by mesospheric SO2). Both sides are verified; the claim cannot be a calibration point.
4. j1: "runaway limit 1.04 to 1.11" is replaced by READ values 1.0512 (Kopparapu 2013), 1.107 at 1 Earth mass, 0.99 at 0.1 and 1.188 at 5 Earth masses (Kopparapu 2014).
5. j5: Harman 2018 (abstract READ from NTRS) conditions its "no false positives" result on the planet lying in the conventional HZ with no irrecoverable water loss, which does not contradict Luger and Barnes 2015 (READ) but restricts it. The 2015 arXiv preprint of Harman (1509.07863, READ abstract) has the opposite headline (K and M stars may build O2); the 2018 version supersedes it.

## 2. Per-row and per-subgroup detail

### 2.1 j1 to j5: Venus-, Mars-, Titan-, Pluto-like and abiotic-O2 atmospheres

**Are they probabilities?** No. Each is the result of (inventory delivered, group 2) x (escape history, tables d and p) x (photochemistry) x (redox, table a). The part 09 weight c[lo,hi] is a mixing weight between overlapping gates, not a measured fraction. The honest model: draw the volatile inventory and redox; compute whether the planet is above the runaway limit, loses its water, keeps CO2-N2; the row label is then an output.

**Ensemble (point 1).**

| Item | Measurement | N | Status | Source |
|---|---|---|---|---|
| Venus Zone planets | 43 Kepler candidates in the Venus Zone; occurrence eta 0.32 (+0.05, -0.07) for M dwarfs and 0.45 (+0.06, -0.09) for GK dwarfs | 43 | READ (abstract) | Kane, Kopparapu, Domagal-Goldman 2014, arXiv:1409.2886 |
| TRAPPIST-1 b dayside | T_B = 503 (+ -27) K at 15 um vs Teq 391 to 400 K; "little or no planetary atmosphere" | 1 planet | READ | Greene 2023, arXiv:2303.14849 |
| TRAPPIST-1 c dayside | f_p/f_star = 421 +- 94 ppm, T_B = 380 +- 31 K; rules out cloud-free O2/CO2 with Psurf >= 10 bar (10 ppm CO2) and pure CO2 above 0.1 bar; Venus analogue with H2SO4 clouds disfavoured at 2.6 sigma | 1 planet | READ | Zieba 2023, arXiv:2306.10150 |
| Hot Rocks Survey | LHS 1478 b eclipse 138 +- 53 ppm (2.8 sigma, one visit a non-detection); GJ 3473 b 186 +- 45 ppm; "may not uniquely distinguish bare rock from atmosphere" | 2 so far | SNIPPET (titles and ids verified: 2410.11048, 2604.02332) | August 2024; Holmberg 2026 |
| Rocky Worlds DDT | GJ 3929 b "likely a bare rock" (Xue et al., ApJL, Dec 2025); first complete data set announced July 2026 | 1 | SNIPPET | search result only; arXiv id not found |
| Evolved TRAPPIST-1 climates | model predictions only; no observation | 7 planets, models | READ (abstract) | Lincowski 2018, arXiv:1809.07498 |

The JWST data are the only exoplanet observations that bear on j1 and j5: two planets (b, c) with no thick CO2 or O2 atmosphere. With N=2, the data allow bare-rock fractions from 0 to 1 for the generator's whole zone-h population; they are an upper-limit check, not a fit. No exoplanet observation exists for Mars-like, Titan-like or Pluto-like atmospheres.

**Physical bracket (point 2).**

| Item | Value | Status | Source |
|---|---|---|---|
| Runaway greenhouse S_eff, Sun, 1-D | 1.0512 (Table 3, 2600 to 7200 K polynomial) | READ | Kopparapu 2013, arXiv:1301.6674 |
| Runaway S_eff vs planet mass | 0.99 (0.1 Me), 1.107 (1 Me), 1.188 (5 Me) | READ | Kopparapu 2014, arXiv:1404.5292 Table 1 |
| Recent Venus 1.776, maximum greenhouse 0.356, early Mars 0.32 (1 Me) | S_eff, same table | READ | Kopparapu 2014 |
| 3-D runaway threshold | about 375 W/m2, "significantly higher than previously thought" | READ (abstract) | Leconte 2013, arXiv:1312.3337 |
| 3-D threshold in S_eff | 375/340.25 = 1.10 if S0 = 1361 W/m2 | COMPUTED | inputs: 375 W/m2, S0/4 |
| Venus Zone inner edge | atmospheric erosion at about 25 times Earth flux | READ (text mention) | Kane 2014 |
| O2 build-up per terrestrial ocean lost | 120 bar (energy-limited) to 240 bar (diffusion-limited) for Venus with 1 TO; 1800 to 2400 bar for 10 TO; about 5 bar/Myr (Earth mass) to 25 bar/Myr (super-Earths); "inefficient oxygen sinks" required | READ | Luger and Barnes 2015, arXiv:1411.7412 |
| Water-lost depends on EUV flux, runaway duration, initial water | stated in abstract | READ | same |
| Abiotic O2 for lifeless HZ planets with weak cold trap | O2-dominated if low non-condensable inventory; ~0.15 bar partial pressure example | READ | Wordsworth and Pierrehumbert 2014, arXiv:1403.2713 |
| Desiccated CO2 worlds, M dwarf | up to about 40 to 50 % of CO2 converted; O2 up to 20 %, O3 above 1 ppm | READ | Gao 2015, arXiv:1501.06876 |
| GJ 1132 b: more than 5 wt % water needed for a steam atmosphere; magma ocean absorbs at most about 10 % of O2 | | SNIPPET | Schaefer 2016, arXiv:1607.03906 |
| Mars-like CO2-N2-Ar, Titan, Pluto: vapour pressure, Jeans retention | volatile retention vs size and temperature ("Schaller and Brown 2007") | SNIPPET; successors Johnson 2015 (arXiv:1503.05315) and Brown 2012 (arXiv:1112.2764) verified by title only | NOT READ |

**Verdict (point 3).** j1 (b): the gate is derived from READ limits; there is no ensemble for the atmosphere fraction, only the two JWST planets as an upper-limit check and the Venus Zone occurrence of planets (not atmospheres). j5 (b): the physics (Luger and Barnes, Wordsworth, Gao, Schaefer 2016) is published and outcome-dependent on the O2 sink assumption (surface absorption of about 100 bar is assumed by Luger and Barnes, READ). j2, j3, j4 (c): no extrasolar ensemble exists; the only constraint is the Solar System (Mars, Titan, Pluto and Triton).

**Knobs (point 4).** If the generator keeps a draw for the residual after the physics gate: knob w_retain in [0, 1] (log-uniform not appropriate because 0 is allowed; use uniform). For j1 the physical gate parameters are S_eff,runaway(M, Teff) with target 1.05 to 1.19 (READ, Kopparapu) and a 3-D alternative of about 1.10 (Leconte, abstract). The water-loss knob is the EUV exposure (group p) and initial water (group 2). The O2-sink knob (fraction of O2 absorbed by the surface) is unsourced except by assumption in Luger and Barnes; bracket [0, 1]. Over what range the targets hold: Kopparapu 2600 to 7200 K, 0.1 to 5 Me.

**Consistency test (point 5).** Simulated: fraction of zone-h planets with thick CO2 or O2 atmospheres at T-1 b and c insolation (Teq 391 to 400 K for b), compared with "no thick atmosphere" for two planets; tolerance is a binomial limit on N=2, so the test is "population must not forbid 0 of 2 with high probability". A selection model is needed: transit geometry and the JWST target list (nearby M dwarfs, short periods).

**No-fixed-content flags (point 6).** j2 to j4 are Solar System-only constraints (Mars, Titan, Pluto). The 97:3 CO2:N2 split (Sossi 2020, READ in part 08) is a magma-ocean calculation, not a Venus copy, but Venus is the test point. The "Schaefer 2012 crust composition" speciation (section 2.4) is Earth-only.

**Group for a builder (point 7).** Builder A: j1, j5, m2 (and the abiotic part of m3) share one model: runaway gate, water loss, O2 build-up, and a photochemistry look-up of O3 columns. j2 to j4 form a second small builder (retention and vapour pressure) with no tuning data.

### 2.2 l1 to l3: white dwarf and pulsar planet options

**Ensemble (point 1).**

| Item | Measurement | N | Status | Source |
|---|---|---|---|---|
| Metal pollution of DA white dwarfs | 48 of 85 (56 %) show metals; 23 (27 %) must be currently accreting; "at least 27 % ... could be as high as about 50 %"; cooling ages 20 to 200 Myr; Teff < 23000 K; 1 of 14 with Mwd > 0.8 Msun accreting; median progenitor about 2 Msun | 85 | READ | Koester, Gaensicke, Farihi 2014, arXiv:1404.2617 |
| Same, volume-limited sample | 45 +- 4 % with atmospheric metals | volume-limited | SNIPPET (via a 2025 paper); Wilson 2019 id not verified | Wilson 2019 |
| Same, earlier | "one quarter to one half" of white dwarfs; accreted mass about that of Vesta or Ceres over 200 Myr | | READ (abstract) | Zuckerman, Melis, Klein 2010, arXiv:1007.2252 |
| Transiting planets at WDs, K2 | 1148 high-probability white dwarfs; 10 eclipsing objects (1 disintegrating planet, 9 likely stellar); hot Jupiters < 1.5 %; Earth-size HZ planets < 28 % | 1148 | SNIPPET (abstract via search; PDF downloaded, text length 107 kB, not read in detail) | van Sluijs and Van Eylen 2018, arXiv:1711.09691 |
| Transiting planets at polluted WDs, TESS | 313 polluted WDs; two known irregular-transit systems recovered, no new; transit fraction 0.8 (+0.6, -0.4) % | 313 | READ (abstract) | Robert et al. 2024, arXiv:2407.21743 |
| Transiting giant at a WD | WD 1856+534 b: P = 1.4 d, Jupiter-sized, M < 13.8 MJ (95 %), < 11.7 MJ if thin disk (age < 10 Gyr); scattered inward without tidal disruption | 1 | READ | Vanderburg 2020, arXiv:2009.07282 |
| Microlensing giant at a WD | MOA-2010-BLG-477Lb: WD 0.53 +- 0.11 Msun, planet 1.4 +- 0.3 MJ, projected separation 2.8 +- 0.5 au (so a > 2.8 au); cites a prediction that over half of WDs have Jovian companions | 1 | READ | Blackman 2021, arXiv:2110.07934 (candidate id 2109.07533 was wrong) |
| Review of survival | Jovian planets must start at least 20 percent beyond the maximum stellar radius (2 Msun AGB example); engulfment depends mainly on maximum RGB radius; "fewer than 10 %" of giant-branch planet hosts have planets within 0.5 au | review | READ | Veras 2016, arXiv:1601.05419 |
| Pulsar planets | PSR B1257+12 (3 planets), PSR B1620-26 b (2.5 MJ, 23 au, circumbinary around a WD plus MSP), PSR J1719-1438 | few | READ (Veras text; Patruno text) | Veras 2016; Patruno and Kama 2017, arXiv:1705.07688 |
| Pulsar planet search | 151 young pulsars, no planets; excludes Mercury analogues in one third of the older sample and M > 0.4 Me with Pb < 1 yr in all but 5 % | 151 | READ (abstract) | Kerr 2015, arXiv:1507.06982 |
| Why pulsar planets are rare | needs companion destruction; dead zones rare | theory | READ (abstract) | Martin 2016, arXiv:1609.06409 |
| Pulsar planet atmospheres | atmospheres can survive X-ray irradiation and pulsar wind if the planet is a super-Earth with a large atmospheric mass fraction (at least 100 Myr) or has a moderately strong magnetosphere (several Gyr); X-ray spectrum of one pulsar | 1 system | READ | Patruno and Kama 2017 |
| Wolszczan and Frail 1992 | first pulsar planets | | NOT REACHED | Nature |

**Physical bracket (point 2).**
l1: P_eq(Ts) of the silicate vapour at the surface temperature (i6 formula from the part 09 table, Kite 2016 READ in part 09, not re-read here). The exosphere criterion (collision-free above a column) is a physics statement, not a probability. l3: survival requires an initial orbit outside the engulfment radius, with tides and mass loss (Veras 2016, READ). Pulsar planets: atmosphere survival from the energy-limited escape and wind stripping calculation of Patruno and Kama (READ).

**Verdict (point 3).**
- **Is the polluted-WD fraction the right tuning constraint for l1? No.** Koester 2014 measures the fraction of white dwarfs that accrete rocky debris (27 to about 50 percent for 20 to 200 Myr cooling ages, Teff below 23000 K, A-type progenitors). It tunes the upstream gate "this WD hosts rocky bodies that reach it" (and, with debris-disk lifetimes, the number of exposed rocky bodies). It says nothing on whether an exposed body keeps a silicate-vapour exosphere. l1 is a deterministic function of Ts through P_eq, so the draw should be removed or kept only as a width on the threshold (0.5 dex DESIGN), with no data to tune it. Verdict for l1 and l2: (c) as probabilities, derivable from physics.
- **l3: (b).** Constraints: one hot giant (WD 1856 b), one cold giant (Blackman), transit upper limits (hot Jupiters < 1.5 % of 1148 K2 WDs; 0.8 percent transit fraction of 313 polluted WDs, both short-period only), and a prediction of more than half of WDs hosting Jovian companions (2ND-HAND). The c[0, 0.1] prior is not supported as a number: the right quantity is P(survive) = P(initial a above the engulfment radius of the progenitor, after tides) x occurrence of giants around 1 to 8 Msun progenitors. This should be derived in the generator from its own giant occurrence (group 1 table f) and Veras' engulfment rule. Bracket-only draw: the widest defensible upper value is the published prediction (> 0.5, 2ND-HAND); the transit limits constrain only P < 27 d.
- **Pulsar planets: (c).** No observation of any pulsar planet atmosphere exists; the record has 3 to 5 planetary systems around about 3000 neutron stars (Patruno: roughly 3000 known; 2200 young pulsars, 400 millisecond pulsars, READ). The only physics is Patruno and Kama's survival calculation. A pulsar-planet atmosphere should be derived (escape against X-ray and wind), not drawn.

**Knobs (point 4).** l1: threshold log10 P_eq,0 = -6 (DESIGN), width 0.5 dex (DESIGN); no target statistic. l3: knob = occurrence of giants beyond the engulfment radius; targets: (i) transit fraction < 1.5 % (K2) and 0.8 (+0.6, -0.4) % for hot companions, (ii) one microlensing detection (cannot constrain a fraction without a survey completeness model). Engulfment knob: a_engulf = 1.2 x R_max(M_prog) for Jovian planets (the "20 percent" figure is an example for a 2 Msun AGB star, READ; use Veras' stellar-model-dependent value).

**Consistency test (point 5).** Simulated hot-giant fraction at WDs vs the K2 and TESS upper limits (inject-and-recover with transit probability a_R/R_WD and detection limits is needed). A population test for l1 does not exist.

**No-fixed-content flags (point 6).** Pulsar-planet atmospheres and white-dwarf exospheres have no Solar System anchor either; the only anchors are one pulsar system and one WD 1856 b / Blackman planet each. The existing recipe (dIW ~ N(-1, 0.5) for WD debris) is from six polluted WDs (Doyle 2019, READ in part 05), an extrasolar constraint, so no flag there.

**Group for a builder (point 7).** One builder: l1 to l3 with the WD pollution gate, engulfment rule and pulsar escape. It shares nothing with the other groups except the rocky-body list from table i.

### 2.3 m2 to m4 and m11: abiotic O3, N2O, PH3, NO2 and CH3Cl

| Row | Ensemble | Physics | Verdict |
|---|---|---|---|
| m2 abiotic O3 | none (no exoplanet O3 detection) | Domagal-Goldman 2014 Table 3 (READ): O3 and O2 columns in 5 % CO2 atmospheres with H2 flux 1e6 and Fe flux 1e10 cm-2 s-1; Gao 2015 (READ): O3 > 1 ppm, O2 up to 20 % for H < 1 ppm; Harman 2018 (abstract READ): lightning NO removes false positives on HZ planets with trace N2 | (b), derive |
| m3 N2O | none | biotic: design flag (life exists). Abiotic: SEP-driven N2O; Airapetian 2020 review (READ passage) says a recent model gives more than 300 times the earlier N2O output; present Earth is about 400 ppbv, biological. Segura 2005 (abstract READ): CH4, N2O, CH3Cl have longer lifetimes on M dwarf planets | (c); design flag plus derived abiotic part |
| m4 PH3 | Venus claim only (Solar System) | Sousa-Silva 2020 (READ): detectable at surface fluxes 1e10 to 1e14 cm-2 s-1, "no known abiotic false positives on terrestrial planets from any source that could generate the high fluxes required"; PH3 + OH lifetime 28 h night, 5 h day (READ) | (c); biotic design flag |
| m11 NO2, CH3Cl | none | CH3Cl biotic (Segura 2005, abstract READ); NO2 from lightning, Harman 2018 abstract (READ) gives no numbers; lightning rates from Hodosan 2016 (arXiv:1606.09172, "guided by Solar System data", SNIPPET) | (c); CH3Cl design flag; NO2 derived |

**Is any of this an observational probability?** No. There is no measurement of any of these gases on any extrasolar planet, so no ensemble can be built today. m2 is an output of a photochemistry grid (columns vary by four orders of magnitude between stellar types, from 1e15 to 1.3e18 cm-2 at the same CO2 and volcanic flux); m3, m4, m11 (CH3Cl) are mostly biosignature design flags whose probability equals "life exists and produces it". A bracket-only seeded draw is the honest answer for m3, m4, m11 with the widest defensible bracket [0, 1] (m4 narrowed to [0, 0.05] if kept as a design choice about biology). The part 09 priors c[0.05, 0.5], c[0, 0.05], c[0, 0.3] have no data behind them: they are choices. m2 should not be a presence probability: the generator can look up an O3 column (Domagal-Goldman Table 3 grid and Gao) as a function of star type, CO2, H2 flux, desiccation and a lightning flag.

**Disputed PH3 record (all READ).** Greaves 2020 (2009.06593): about 20 ppb, abiotic photochemistry short by more than four orders of magnitude. Villanueva 2021 (2010.14305): < 1 ppbv at 3 sigma (0.7 to 1.5 for line-width choices). Lincowski 2021 (2101.09837): the line forms above 80 km and fits mesospheric SO2. Do not tune m4 to Venus.

**Knobs and test.** None of the proposed draws can be tuned: no target statistic exists. Consistency test: none possible except "no rocky planet in the generator shows PH3 above 1 ppb unless life flag" if a future observation requires it.

**No-fixed-content flags.** Every constraint on m3 (400 ppbv), m4 (Venus 20 ppb) and m11 (lightning NO flux based on modern Earth; SNIPPET) is Solar System-only. The rule is broken if any knob is tuned to those.

**Group for a builder.** m2 with Builder A (section 2.1). m3, m4, m11 (CH3Cl) form a "design flags" group that belongs with the life-chemistry design, not with a fit.

### 2.4 Halogens: HCl and HF (m5, m6) and the derived-budget question

**Can a derived halogen budget replace c[0.1, 0.5]? Yes for the stellar input, only partly for the planet.** HCl and HF are equilibrium species of the Cl and F in the crust or mantle, not a coin flip. Schaefer, Lodders and Fegley 2012 (arXiv:1108.4660, READ) show that HCl and HF are "the dominant Cl- and F-bearing gases at lower temperatures", appear at percent level only for the continental-crust composition (Table 2: F 0.053 wt %, Cl 0.047 wt %) and are below the 0.1 % plot threshold for BSE composition (F 0.002 wt %, Cl 0.004 wt %); the part 09 values (HCl 1.4 % at 1000 K 100 bar; HF 4.5 % at 1000 K 1e-6 bar) are crust-composition results. COMPUTED from Table 2: crust over BSE enrichment is 26.5 for F and 11.8 for Cl; molar Cl/F is 0.48 in the crust and 1.07 in the BSE; the solar-photosphere Cl/F number ratio is 10^(5.31 - 4.40) = 8.1 (meteoritic 10^(5.23 - 4.42) = 6.5). So the planet's halogen mixing ratio is set mainly by differentiation (crust enrichment) and volatile depletion, not by the star.

**Stellar side (point 1).**

| Quantity | Value | N | Status | Source |
|---|---|---|---|---|
| Solar F | log eps 4.40 +- 0.25 (sunspot HF) and 4.42 +- 0.06 (meteoritic) | 1 star | READ | Asplund, Amarsi, Grevesse 2021, arXiv:2105.01661 |
| Solar Cl | log eps 5.31 +- 0.20 (sunspot HCl) and 5.23 +- 0.06 (meteoritic) | 1 star | READ | same |
| Cl in other stars | 15 giants plus 1 M dwarf, [35Cl/Fe] = -0.10 +- 0.15 over -0.72 < [Fe/H] < 0.20 | 16 | READ (abstract) | Maas, Pilachowski, Hinkle 2016, arXiv:1609.01626 |
| Cl in M giants | 52 M giants, -0.5 < [Fe/H] < 0.12; [Cl/Ca] about solar; Cl from core-collapse supernovae plus some Type Ia | 52 | READ (abstract) | Maas and Pilachowski 2021, arXiv:2101.11068 |
| F in K giants | 49 nearby K giants; no evidence for the neutrino process | 49 | READ (abstract) | Jonsson et al. 2017, arXiv:1612.04384 |
| F in bulge giants | chemical evolution of F in three bulge fields | few | READ (abstract) | Jonsson et al. 2014, arXiv:1403.2594 |
| F in disk giants | 12 thin and thick disk giants; larger dispersion than Na | 12 | READ (abstract) | Brady et al. 2024, arXiv:2406.02855 |

No F or Cl abundance exists for main-sequence FGK stars from direct lines; both come from molecular lines (HF, HCl) in cool giants or sunspots, so the generator's host F and Cl must be scaled from [Fe/H] and [alpha/Fe] by a chemical-evolution relation. The published scatter (0.15 dex for Cl) is smaller than the planetary processing uncertainty.

**Planet side.** Halogen depletion in terrestrial bodies: Clay et al. 2017 (Nature 551, 614; SNIPPET only: Cl in chondrites about six times lower than earlier estimates, Br/Cl and I/Cl matching BSE) says most Earth halogens sit in surface reservoirs. Volcanic degassing on Earth: HCl 4.3 +- 1 Tg/a, HF 0.5 +- 0.2 Tg/a, S-Cl and S-F partitioning contrasts of about 9 and 36 (SNIPPET from a search result). Every constraint on the star-to-planet transfer (condensation, depletion, crust enrichment) comes from Earth, chondrites and Venus: a Solar System-only calibration.

**Verdict.** (b) for a derived HCl/HF budget (stellar scaling is data-supported, speciation is published), but the planet-side depletion and crust-enrichment factors have no extrasolar data. As a probability it is (c). Recommend: compute HCl and HF mixing ratios from (a) host [Cl/Fe]-like scaling, (b) a volatility depletion curve vs formation temperature (condensation temperatures are not verified here: RECALLED), (c) a crust-enrichment factor tied to the differentiation draw (table a, groups 2 and 3), (d) Schaefer speciation, (e) removal by an ocean or by condensation as salts (derived from the j and c tables). The only random elements are inventory and differentiation. If a seeded presence draw is kept, the widest defensible bracket is [0, 1] for "above 1 ppm"; c[0.1, 0.5] has no support. Knobs for the derived approach: depletion exponent (unsourced), crust enrichment 11.8 (Cl) and 26.5 (F) (Earth-only, COMPUTED from READ Table 2). Test: the Venus HCl (about 0.5 ppm) and HF (about 5 ppb) values are RECALLED and Solar System-only; no exoplanet test exists.

**Group for a builder.** m5 and m6 alone, plus the speciation call (Schaefer fit exists in part 05).

### 2.5 p4 and p6: M-dwarf saturation time and the rotation-track weights

**Question: do the activity-rotation papers fix tau_sat and the initial rotation spread by mass with enough N?**

**Published saturation times (point 1).**

| Quantity | Value | N | Calibration of age | Status | Source |
|---|---|---|---|---|---|
| X-ray tau_sat, early M (B-V 1.275 to 1.41) | log tau_sat = 8.21 +- 0.04 (stat) +- 0.31 (sys) yr, i.e. 162 Myr (stat range 148 to 178; sys range 79 to 332 Myr, COMPUTED) | 55 stars (717 stars in all 7 bins from 13 clusters) | open clusters 5 to 740 Myr; field stars at assumed 4.5 Gyr shown but not fitted | READ (Table 2) | Jackson, Davis, Wheatley 2012, arXiv:1111.0031 |
| Other Jackson bins | log tau_sat 7.84 to 8.35 for 0.29 <= B-V < 1.275; "no obvious trend with colour; scatter around about 100 Myr" | 67 to 220 per bin | same | READ | same |
| Saturation level and slope, early M | log(LX/Lbol)sat = -3.14 +- 0.02 +- 0.35; alpha = 1.18 +- 0.31 | 55 | same | READ | same |
| Rossby saturation | Ro_sat = 0.13 +- 0.02 (beta = -2) or 0.16 +- 0.03 (beta = -2.7); log(LX/Lbol)sat = -3.13; beta = -2.70 +- 0.13 | 824 stars | Rossby convective times, rotation periods | READ | Wright et al. 2011, arXiv:1109.4634 |
| Fully convective stars | Ro_sat = 0.14 (+0.08, -0.04); log(LX/Lbol)sat = -3.05 (+0.05, -0.06); unsaturated fully convective stars follow the same relation (beta = -2.3) | 19 new stars plus earlier sample | MEarth periods, Chandra, ROSAT | READ | Wright et al. 2018, arXiv:1807.03304  |
| Johnstone Ro_sat | Ro_sat = 0.0605 +- 0.0033 (with convective times); R_X,sat = 5.135e-4 +- 3.3e-5; slopes -0.135 +- 0.030 and -1.889 +- 0.079 | Wright 2011 sample | Ro scale depends on tau_c | READ | Johnstone 2021, arXiv:2009.07695 |
| eROSITA M dwarfs | 704 stars K5 to M7; rotation for 180; saturated regime Prot < 8.5 d adopted; mass dependence of the saturated LX/Lbol examined | 704 (180 with Prot) | TESS periods | READ (abstract, passages) | Magaudda et al. 2022, arXiv:2106.14548 |
| Saturated to unsaturated transition, mid-to-late M | average transition age 2.4 +- 0.3 Gyr (assuming constant thin-disk SFH over 8 Gyr); Prot < 10 d stars 2.0 +- 1.2 Gyr, 10 to 90 d 5.6 +- 2.7 Gyr, > 90 d 12.9 +- 3.5 Gyr; fast rotators 0.6 +- 0.3 Gyr (0.2 to 0.3 Msun), 1.9 +- 0.9 (0.1 to 0.15), 2.7 +- 1.0 (0.15 to 0.2), combined 2.3 +- 1.3 Gyr (0.1 to 0.2 Msun) | 219 single stars within 15 pc, 0.1 to 0.3 Msun | Galactic kinematics (UVW, vertical velocity dispersion vs age of Lu 2021); 64 % thin disk, 5 % thick, 31 % unassigned | READ | Medina et al. 2022, arXiv:2205.02331 |
| MEarth rotation | 387 rotation periods (0.1 to 140 d); P < 10 d younger than 2 Gyr on average, P > 70 d about 5 Gyr | 387 | kinematics | READ (abstract) | Newton et al. 2016, arXiv:1511.00957 (candidate "Newton 2017" is the H-alpha paper, 1611.03509, SNIPPET) |
| H-alpha activity lifetime by spectral type | M0 0.8, M1 0.4, M2 1.2, M3 2.0, M4 4.5, M5 7.0, M6 7.0, M7 8.0 Gyr (uncertainties 0.4 to 1.0 Gyr) | 38000 low-mass stars | 1-D dynamical heating simulations | READ (Table 2) | West et al. 2008, arXiv:0712.1590 |
| UV transition region | constant level for 240 +- 30 Myr, then declines 2.1 dex to 10 Gyr (M0 to M2.5, 0.3 to 0.6 Msun) | young groups, Hyades, field | clusters and field | READ (abstract) | Loyd et al. 2021 (HAZMAT VII), arXiv:2011.10158 (not Richey-Yowell lead) |

**Answer for p4.** For 0.5 to 0.6 Msun (B-V 1.3 to 1.4, K7 to M1) the X-ray tau_sat is measured: 162 Myr, with a 0.31 dex systematic range, from N=55 stars with cluster ages. For 0.2 to 0.3 Msun the activity age of fast rotators is 0.6 +- 0.3 Gyr and for 0.1 to 0.2 Msun it is 2.3 +- 1.3 Gyr, from N=219 stars with kinematic ages that assume constant star formation (a lower limit, Medina says). The overall M-dwarf mean is 2.4 +- 0.3 Gyr. The Halpha activity lifetimes of West 2008 are longer (M4 4.5, M5 to M6 7, M7 8 Gyr) and use different physics (Halpha, not X-ray; heating simulation ages). So a single LU[0.162, 2] Gyr mixes mass: the lower edge is the measured early-M X-ray value, the upper edge is the fully convective value from Medina, and late M can exceed it (West 2008). Mid-M (0.3 to 0.5 Msun) has no direct X-ray tau_sat in the record.

**Recommended replacement.** tau_sat(M) piecewise in log mass: 0.16 Gyr (sys 0.08 to 0.33) at 0.5 to 0.6 Msun, 0.6 +- 0.3 Gyr at 0.2 to 0.3 Msun, 2.3 +- 1.3 Gyr at 0.1 to 0.2 Msun, interpolated; spread drawn as log-normal with the published uncertainty. Or compute it: R_X = R_X,sat for Ro <= Ro_sat on the drawn spin-down track (Mors code, MIT licence; README states a function that returns the age a star stays saturated; masses 0.1 to 1.25 Msun). Caveat: Johnstone says the tracks do not describe fully convective M dwarfs well after about 2 Gyr and that fully convective stars may not spin down much. The Ro_sat constants are convention dependent (0.0605 with times against 0.13 to 0.16 with Wright 2011 times): never mix them.

**Knobs for p4 and targets.** (i) tau_sat for 0.55 Msun: target log tau_sat = 8.21 +- 0.04 +- 0.31; (ii) tau_sat for 0.25 Msun: 0.6 +- 0.3 Gyr; (iii) tau_sat for 0.15 Msun: 2.3 +- 1.3 Gyr; (iv) R_X,sat: log -3.14 +- 0.02 +- 0.35 (Jackson), -3.05 (+0.05, -0.06) (Wright 2018), 5.135e-4 +- 3.3e-5 (Johnstone), -3.13 (Wright 2011); (v) post-saturation slope alpha 1.18 +- 0.31 (Jackson) with 1.22 +- 0.10 mean over bins. Validity: clusters 5 to 740 Myr; Medina 15 pc volume-complete 0.1 to 0.3 Msun.

**Consistency test for p4 and p6.** Simulated statistic: fraction of M dwarfs in the saturated regime vs age and mass, compared with (a) Medina's 2.4 +- 0.3 Gyr and mass trend, (b) the age distribution of rapid rotators (Prot < 10 d: 2.0 +- 1.2 Gyr), (c) the fraction of active stars (74 +- 4 %, 19 +- 4 %, 8 +- 3 % in the three period bins), (d) Jackson's 162 Myr. Tolerance: the published uncertainties above. Selection models are needed: volume-complete (Medina and Pass use 15 pc), detection bias in X-ray catalogues (Magaudda includes only X-ray-detected stars) and in rotation catalogues (Johnstone: Kepler misses the slowest rotators). Kinematic age assignments carry thick-disk contamination (5 % thick, 31 % unassigned).

**p6: weights.** Tu 2015 (arXiv:1504.04546, READ) uses the 10th, 50th and 90th percentile tracks of the initial rotation distribution for 0.9 to 1.1 Msun only, with t_sat of 5.7, 23 and 226 Myr as in part 09; the constraints are clusters at 150, 550 and 1000 Myr (230, 134, 36 stars) plus NGC 6530 (28 stars) and h Per (117 stars); the abstract gives saturation lasting "from about 10 Myr to about 300 Myr". Gallet and Bouvier 2013 (arXiv:1306.2130, READ) use the 25th, 50th and 90th percentiles as slow, median, fast with initial periods of 10, 7 and 1.4 d from 13 cluster distributions. Johnstone 2021 (READ) uses the 5th, 50th and 95th percentiles at 150 Myr for four masses, 12 clusters in 7 age bins, and provides tracks at every integer percentile 2 to 98 for masses 0.1 to 1.25 Msun (code MIT-licensed, `raw.githubusercontent.com/ColinPhilipJohnstone/Mors/master`, README READ; the track data were not downloaded). No paper assigns weights to three tracks. The 0.25/0.5/0.25 weights can be justified (COMPUTED, reasoning only): tracks at the 10th/50th/90th percentile are close to the medians (12.5, 50, 87.5) of quartile bins [0, 0.25], [0.25, 0.75], [0.75, 1], which have probability 0.25, 0.5, 0.25. The best replacement is a uniform percentile u in [0.02, 0.98] mapped onto Johnstone's integer-percentile tracks, which also makes the mass dependence automatic (tracks for 0.1 to 1.25 Msun) and removes the solar-mass restriction of Tu 2015. The 200-fold spread in L_X at 100 Myr (Tu, part 04) is real and must be drawn. Knob: percentile u (no tuning needed: it is the observed 150 Myr distribution); target: the cluster rotation percentiles (READ). Test: reproduce Johnstone's modelled X-ray distributions against cluster data within the published scatter.

**No-fixed-content flags.** Johnstone's tracks are fitted to open clusters and field stars of all masses, not the Sun alone; the Sun appears only as a minor weight in fits. No Solar-System-only knob in p4 or p6.

**Group for a builder.** p4, p6 with p3 and p5 (same Mors tracks and Rossby relation). One worker can replace tau_sat and the track weights together and build the consistency tests against Jackson, Medina, Pass and Wright.

## 3. Sources actually verified (arXiv PDF text extracted; no licence statement found in the extracted text unless noted)

Read in the body text (READ): Kopparapu 2013 (1301.6674), Kopparapu 2014 (1404.5292), Kane 2014 (1409.2886, N=43), Luger and Barnes (1411.7412), Wordsworth and Pierrehumbert (1403.2713), Gao (1501.06876), Domagal-Goldman (1407.2622), Greene (2303.14849), Zieba (2306.10150), Jackson (1111.0031, 717 stars), Wright 2018 (1807.03304), Medina 2022 (2205.02331, N=219), West 2008 (0712.1590, 38000 stars), Gallet and Bouvier (1306.2130), Tu (1504.04546).07695), Koester (1404.2617, N=85).07282), Blackman (2110.07934), Veras (1601.05419), Patruno and Kama (1705.07688), Schaefer 2012 (1108.4660), Asplund (2105.01661).06593).14305), Lincowski 2021 (2101.09837).05224), Segura 2005 (astro-ph/0510224), Airapetian 2020 (1905.05093), Harman 2018 (NTRS 20180007112 abstract), Mors README and LICENSE (github raw, MIT).
Abstract only (READ abstract): Leconte (1312.3337), Lincowski 2018 (1809.07498), Wright 2011 (1109.4634, N=824), Newton 2016 (1511.00957, N=387), Pass 2023 (2306.00799, N=123), Magaudda (2106.14548, N=704), Loyd 2021 (2011.10158), Pineda (2102.12485), Richey-Yowell 2023 (2305.06561), Zuckerman (1007.2252), Robert (2407.21743, N=313), Martin (1609.06409), Kerr (1507.06982, N=151), Maas 2016 and 2021 (1609.01626, 2101.11068), Jonsson (1612.04384, 1403.2594), Brady (2406.02855), Harman 2015 preprint (1509.07863).
Title and id only: France 2016 (1602.09142), Shkolnik and Barman (1407.1344, 2ND-HAND), van Sluijs and Van Eylen (1711.09691, SNIPPET), Hot Rocks I and V (2410.11048, 2604.02332), Johnson 2015 (1503.05315), Brown 2012 (1112.2764), Meadows 2018 (1705.07560 and the Proxima paper 1608.08620).

Candidate corrections: Kane 2014 is 1409.2886; Blackman is 2110.07934; Martin 2016 is "Why are pulsar planets rare?"; "Newton 2017" is the Halpha paper (1611.03509), the MEarth rotation paper is Newton 2016; "Jonsson 2014" is 1403.2594 (bulge), the solar neighbourhood paper is Jonsson 2017; HAZMAT VII is led by Loyd; "Ostberg and Kane 2019" not found (Ostberg et al. 2023, AJ 165, 168, arXiv 2302.03055, SNIPPET); Harman 2018 has no arXiv version (1509.07863 is the 2015 preprint with a different thesis); Loyd 2016 and Youngblood 2016 MUSCLES ids not located; Goldblatt 2013 has no arXiv version; Wright 2011 beyond the abstract and Wilson 2021 not checked.

## 4. Boss must supply

1. Schaller and Brown 2007 (ApJ 659, L61), full text or the retention-diagram table (no arXiv, IOP not reachable); Johnson 2015 (1503.05315) and Brown 2012 (1112.2764) are on arXiv and can be read next.
2. Harman 2018 (ApJ 866, 56) full text for the lightning NO flux and N2 threshold (only the abstract was read).
3. Wilson 2019 (volume-limited fraction of polluted WDs, 45 +- 4 %) arXiv id for the l3 upper bracket.
4. Rocky Worlds DDT and Hot Rocks papers with arXiv ids (GJ 3929 b, LHS 1478 b, GJ 3473 b, others) for a proper N of bare-rock versus atmosphere planets in the generator's zone-h range.
5. A F and Cl condensation and depletion source for the halogen budget (volatility of Cl and F in the nebula, Lodders 2003 or Clay 2017 full text) and the Venus HCl/HF measurements (only for test points).
6. Mid-M (0.3 to 0.5 Msun) X-ray saturation time: no paper found; Magaudda 2022 full text (Psat, mass dependence) or Reiners 2022 (arXiv:2204.00342, downloaded, not read) would fix it.
7. The Mors percentile track data (download from GitHub) to confirm the integer percentiles and the saturation-lifetime function.

## 5. Recalled / to verify

- Venus HCl about 0.5 ppm and HF about 5 ppb (RECALLED, not used except as a flagged test point).
- Cl and F 50 percent condensation temperatures (RECALLED).
- Wilson 2019 polluted fraction 45 +- 4 percent (SNIPPET).
- Hot Rocks and Rocky Worlds results (SNIPPET).
- Goldblatt 2013 numbers (282 and 294 W/m2) (SNIPPET).
- Airapetian 2016 N2O model and Hodosan 2016 lightning rates (SNIPPET); Barth 2024 summary of the Harman NO flux (SNIPPET).
- Clay 2017 halogen depletion factors; Earth degassing fluxes HCl 4.3 and HF 0.5 Tg/a (SNIPPET).
- Ostberg et al. 2023 arXiv id 2302.03055 (SNIPPET).
- Kite 2016 P_eq formula for l1 (READ in part 09, not re-read here).
