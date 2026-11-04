# Part 05: Mixing ratios by body class

Research lane 05 for the v9.0 atmosphere composition work. Status labels follow `_worker-rules.md`: READ, 2ND-HAND, DATA, COMPUTED, SNIPPET, RECALLED, NOT REACHED. Solar System bodies appear only in rows marked TEST (calibration points, never presets). Where sources disagree, both values and the spread are given.

Scratch files used for COMPUTED rows: thermodynamic data from the cantera 3.2.0 wheel (`cantera/data/nasa_gas.yaml`, `nasa_condensed.yaml`, BSD-3, originating from McBride, Gordon and Reno, NASA TM-4513, 1993); small Python scripts for equilibrium constants and a speciation solver (ideal gas, no graphite, no solubility, no S2).

Network note for the coordinator: reachable and used were export.arxiv.org and arxiv.org PDFs, ntrs.nasa.gov API downloads, nature.com PDF, par.nsf.gov PDF, faculty.washington.edu PDF, en.wikipedia.org raw, pypi.org wheels. Failed: the arXiv API (HTTP 429), arxiv.org/search (400), science.org (403), pmc.ncbi.nlm.nih.gov (bot stub), hal-insu.archives-ouvertes.fr (bot check).

---

## 1. Scope and how the numbers feed the generator

Two generator paths are supported by this part.

**Path A, draw bulk properties first, then fill composition.** Inputs: surface or reference pressure P, temperature T, surface gravity g, body mass, an oxygen fugacity offset dIW (log10 units relative to the iron-wuestite buffer), host spectral type, age, and a class letter. Output: a mole-fraction vector over the species list in section 13, normalised to 1, whose mean molecular weight then fixes scale height and density. Sections 3 to 10 each end with a recipe in pseudo-formulas of the form `mixture = f(P, T, g, dIW, star, age)`.

**Path B, build up from basics.** Inputs: stellar elemental abundances (section 2), a volatile inventory, mantle redox, accretion zone, and an escape history. Output: the same mole-fraction vector, derived by (i) condensing refractory and ice phases at the formation temperature, (ii) degassing the mantle at the magma ocean stage (section 3), (iii) re-equilibrating at the final T and P, and (iv) removing light species by escape. Path B is the physically honest route; Path A can be fitted to it by sampling Path B offline and storing the result as class templates with the quoted spreads.

Rules the generator should respect everywhere:
- Number fractions of heavy elements scale together with a metallicity factor m (solar units) while H and He are only weakly changed.
- Scatter should be drawn log-normal around the median, with the log10 range given in each table, because the quoted spreads in the literature are mostly one-sigma in dex.
- Condensation removes species where T is below the condensation temperature at the local partial pressure (section 12 gives the vapour-pressure relations that were reached). The vapour-pressure relations must be applied before mixing ratios are used for optics.
- Redox state (dIW) is a mantle property and is not the same as the atmospheric oxygen fugacity at the surface; the buffer fits in section 3.2 give fO2 at T and P.

---

## 2. Elemental budget and stellar to planetary links

### 2.1 Solar reference abundances

Number abundances relative to H are 10^(log epsilon - 12), COMPUTED from the READ log epsilon values of [R0168] ([R0168], photosphere column unless noted). [R0169] values are listed where they differ by more than 0.05 dex.

| Element | log eps (A21 photosphere) | n/nH (COMPUTED) | Other value | Source | Status |
|---|---|---|---|---|---|
| H | 12.00 | 1 | CI chondrite 8.22 | [R0168] | READ |
| He | 10.914 | 0.0821 | 10.93 (A09) | [R0168, R0169] [R0169] | READ |
| C | 8.46 | 2.88e-4 | 8.43 (A09), CI 7.39 | [R0168] | READ |
| N | 7.83 | 6.76e-5 | CI 6.26 | [R0168] | READ |
| O | 8.69 | 4.90e-4 | CI 8.39 | [R0168] | READ |
| Ne | 8.06 | 1.15e-4 | 7.93 (A09) | [R0168] | READ |
| Na | 6.22 | 1.66e-6 | 6.24 (A09) | [R0168] | READ |
| Mg | 7.55 | 3.55e-5 | 7.60 (A09) | [R0168] | READ |
| Al | 6.43 | 2.69e-6 | | [R0168] | READ |
| Si | 7.51 | 3.24e-5 | | [R0168] | READ |
| P | 5.41 | 2.57e-7 | CI 5.43 | [R0168] | READ |
| S | 7.12 | 1.32e-5 | | [R0168] | READ |
| Cl | 5.31 | 2.04e-7 | 5.50 (A09), CI 5.23 | [R0168] | READ |
| Ar | 6.38 | 2.40e-6 | 6.40 (A09) | [R0168] | READ |
| K | 5.07 | 1.17e-7 | 5.03 (A09) | [R0168] | READ |
| Ca | 6.30 | 2.00e-6 | | [R0168] | READ |
| Ti | 4.97 | 9.33e-8 | | [R0168] | READ |
| Fe | 7.46 | 2.88e-5 | 7.50 (A09) | [R0168] | READ |
| Ni | 6.20 | 1.58e-6 | | [R0168] | READ |
| Kr | 3.12 | 1.32e-9 | 3.25 (A09) | [R0168] | READ |
| Xe | 2.22 | 1.66e-10 | 2.24 (A09) | [R0168] | READ |

Derived ratios (COMPUTED): C/O = 0.589, N/O = 0.138, S/O = 0.027, Mg/Si = 1.10, Fe/Si = 0.89, Ar/Ne = 0.021. Photospheric mass fractions X = 0.7438, Y = 0.2423, Z = 0.0139, Z/X = 0.0187 +/- 0.0009 [R0168]; protosolar X = 0.7121, Y = 0.2725, Z = 0.0154 (READ, same). [R0169]: X = 0.7381, Y = 0.2485, Z = 0.0134, Z/X = 0.0181 (READ). [R2733] give Z/X = 0.0191 (READ, [R2733] table). The spread in solar Z is therefore 0.0134 to 0.0154 depending on photospheric versus protosolar and the author. A generator that scales other stars should state which Z_sun it uses.

### 2.2 Cool solar gas, all heavy elements locked in simple hydrides

COMPUTED from the table above: assume all O in H2O (or 80% of O if 20% is held by silicates), all C as CH4, N as NH3, S as H2S, and H closed as H2.

| Species | x (all O in gas) | x (80% O in gas) | Status |
|---|---|---|---|
| H2 | 0.857 | 0.857 | COMPUTED |
| He | 0.141 | 0.141 | COMPUTED |
| H2O | 8.4e-4 | 6.7e-4 | COMPUTED |
| CH4 | 4.96e-4 | 4.96e-4 | COMPUTED |
| NH3 | 1.16e-4 | 1.16e-4 | COMPUTED |
| Ne | 1.97e-4 | 1.97e-4 | COMPUTED |
| H2S | 2.3e-5 | 2.3e-5 | COMPUTED |
| Ar | 4.1e-6 | 4.1e-6 | COMPUTED |

For comparison, the [R4863] equilibrium transition (section 7) says CH4 holds the carbon below about 500 K at planet pressures, so this is the correct low-temperature end member.

### 2.3 How host-star composition maps to planet composition

| Finding | Value | Source | Status |
|---|---|---|---|
| Planetary Mg/Si and Fe/Si track the stellar values; C/O of solids depends only weakly on stellar C/O; rocky planets formed inside the ice line carry very little carbon | qualitative | [R4336] abstract and conclusions | 2ND-HAND (abstract read, body not fully read) |
| Refractory material is more than 50 wt% of accreted solids, up to 30 wt% refractory organics | >50 wt% | [R4335] | READ |
| In planet-host dwarfs, C/O > 0.8 in 34% and <= 0.8 in 66%; Mg/Si < 1 in 56%, 1 to 2 in 44%, > 2 in 0% | 34% | Delgado [R1060] | READ |
| [R0535]: carbon-rich dwarfs are under 0.13% and Mg/Si peaks near 1.0 with about 60% between 1 and 2 | <0.13% | [R0535, R0536] | SNIPPET (PDF downloaded, abstract only used) |
| Disagreement: C/O > 0.8 fraction is 34% (Delgado Mena) versus under 0.13% (Brewer and Fischer). The spread is two orders of magnitude, so a generator should use a small tail (about 0.1 to a few percent) and treat the 34% as an upper bound from a heterogeneous sample | | both above | READ and SNIPPET |
| Compositions of terrestrial planets range from Earth-like to graphite and SiC dominated for high C/O stars (threshold C/O = 0.8) | | [R0455] | READ |
| Mg/Si and the oxygen budget set mantle mineralogy | | [R4470] | READ |

Recipe (Path B start): draw stellar [Fe/H], C/O, Mg/Si from the Delgado Mena distribution (sigma roughly 0.2 dex in each), then set planetary Mg/Si and Fe/Si equal to stellar; draw volatile delivery from the accretion zone (section 3.1).

### 2.4 Giant-planet bulk metallicity: mass-metallicity draw rule

| Relation | Value | Validity | Source | Status |
|---|---|---|---|---|
| Heavy-element mass Mz = (57.9 +/- 7.03) (M/MJ)^(0.61 +/- 0.08) Earth masses; 1 sigma intrinsic scatter factor 1.82 (10^sigma, 0.26 dex) | | 47 transiting giants, roughly 0.1 to 10 MJ | [R4352] | READ |
| Zpl/Zstar = (9.7 +/- 1.28) (M/MJ)^(-0.45 +/- 0.09), Zstar = 0.014 x 10^[Fe/H] | | same | [R4352] | READ |
| Mz = Mcore + fZ (Mp - Mcore), Mcore = 14.7 (+1.8/-1.6) Earth masses, fZ = 0.09 +/- 0.01, astrophysical scatter 0.66 +/- 0.08 x MZ | | 147 planets | [R0710] | READ |
| No scatter was found for the Zpl/Zstar fit itself; use the 0.26 dex of the Mz relation as the proxy | | | [R4352] | NOT REACHED (see Boss list) |

COMPUTED comparison, with Z_star = 0.0139 ([R0168] photosphere):

| M/MJ | Mz (Earth masses) | Z by Mz fit | Z by ratio fit | Ratio fit in solar units | Chachan Mz | Chachan Z |
|---|---|---|---|---|---|---|
| 0.03 | 6.8 | 0.715 | 0.653 | 47 | (all mass) | 1.0 |
| 0.1 | 14.2 | 0.447 | 0.380 | 27 | 16.2 | 0.51 |
| 0.3 | 27.8 | 0.291 | 0.232 | 17 | 22.0 | 0.23 |
| 1 | 57.9 | 0.182 | 0.135 | 9.7 | 42.0 | 0.13 |
| 3 | 113 | 0.119 | 0.082 | 5.9 | 99 | 0.10 |
| 10 | 236 | 0.074 | 0.048 | 3.4 | 299 | 0.094 |

The two Thorngren fits differ by about 30% in Z at Jupiter mass and the Chachan fit is lower for sub-Saturns. The bulk Z is not the atmospheric Z: atmospheric retrievals give a spread around bulk values (see 2.5), so the generator should draw atmospheric enrichment as `log10(m_atm) = log10(Z_bulk/Z_star) + N(0, 0.3)` with the 0.3 dex width RECALLED from the retrieval scatter below, not a sourced number.

### 2.5 Atmospheric metallicity and C/O trends

| Result | Value | Source | Status |
|---|---|---|---|
| log(H2O/H) relative to host expectation = -1.09 log(M/MJ) - 0.95 across 19 planets (water abundance falls with mass) | slope -1.09 | [R4663] | READ |
| Solar-system giants: log(CH4/H) = -1.11 log(M/MJ) + 0.38 | slope -1.11 | [R4663] | READ (this is a TEST relation) |
| H2O in hot giants is sub-stellar while Na and K are stellar to super-stellar; roughly 100x H2O depletion is possible from high C/O | | [R4663] | READ |
| A later revisit finds the data cannot distinguish a linear mass relation from independence | | [R4258] | 2ND-HAND (abstract) |
| C/O between 0.5 and 2 changes H2O and CH4 by orders of magnitude; HCN and C2H2 appear at C/O >= 1; TiO and VO are underabundant at C/O >= 1 | | [R2801, R2803], [R2801] | READ |
| Madhusudhan classes O1, O2 (C/O < 1), C1, C2 (C/O >= 1, split at Teq near 1200 K) | | same | READ |
| TEST: WASP-77Ab [C/H] = -0.46 (+0.17/-0.16), [O/H] = -0.49 (+0.14/-0.12), [(C+O)/H] = -0.48 (+0.15/-0.13), C/O = 0.59 +/- 0.08 (they quote solar C/O as 0.55; [R0168] gives 0.589, a 0.04 spread) | | Line et al. 2021, [R2701] | READ |
| TEST: WASP-107b Z = 43 +/- 8 x solar, CH4 1.0 +/- 0.5 ppm (4.2 sigma; about 3 orders of magnitude below equilibrium because of vigorous mixing, Kzz = 10^(11.6 +/- 0.1) cm2/s), Tint = 460 +/- 40 K | | [R4106] | READ |
| TEST: WASP-39b best-fit models 3 to 10 x solar metallicity with sub-solar to solar C/O; SO2 is photochemical (below 1e-12 in thermal equilibrium at 10 x solar) | | [R0064, R4437] | READ |

Recipe for giants (any T): `m = 10^N(log10 Zratio(M), 0.3)`; `C/O = 0.55 * 10^N(0, 0.15)`, truncated at 0.2 to 1.5, with a 2 to 5% chance of C/O >= 0.8 from section 2.3. These widths are RECALLED design choices bracketed by the Madhusudhan and Line cases above, not published fits.

---

## 3. Outgassed secondary atmospheres of rocky planets

### 3.1 Volatile inventory

| Input | Value | Source | Status |
|---|---|---|---|
| Volatile contents in wt% (delivered, stored in mantle) for Low: H2O 0.05, CO2 0.01 delivered; 0.015, 0.0022 stored | | [R3378], Sci. Rep. 10, 10907 (nature.com PDF) Table 1 | READ |
| High: delivered H2O 0.5, CO2 0.1; stored 0.045, 0.005 | | same | READ |
| Dry: delivered H2O 0, CO2 0.6; stored 0.005, 0.018 | | same | READ |
| Partition D_H2O = 0.01; extrusive melt fraction 10%; maximum outgassing at 2 to 4 Earth masses | | same | READ |
| Mantle at 4000 K, 0.05 wt% H2O (550 ppm) and 0.01 wt% CO2 (130 ppm) | | [R2286] | READ |
| Example 1 Earth mass planet with XCO2 = 22 ppm and XH2O = 150 ppm | | [R3378] | READ |
| Reduced mantles give lower surface pressure and a thicker atmosphere (more H2/CO per unit gas mass); CH4 excluded from their model | | [R3378] | READ |

The spread across Low, High and Dry is a factor of 10 in water and 60 in CO2; the generator should draw the delivered volatile fraction log-uniformly between these cases.

### 3.2 Redox: oxygen fugacity buffers

[R1989] buffers as quoted in [R2286] (2ND-HAND): `log10 fO2 = A - B/T + C (P - 1)/T`, T in K, P in bar.

| Buffer | A | B | C | Status |
|---|---|---|---|---|
| IW | 6.899 | 27714 | 0.05 | 2ND-HAND |
| QIF | 7.679 | 29673 | 0.05 | 2ND-HAND |
| QFM | 8.555 | 24014 | 0.092 | 2ND-HAND |
| NiNiO | 8.951 | 24556 | 0.046 | 2ND-HAND |

Cross-check: the NASA stoichiometric FeO data (DATA, cantera nasa_condensed.yaml) give log fO2(IW) = -12.20 at 1500 K and -7.74 at 2000 K (COMPUTED). The RECALLED [R1464] fit (-27489/T + 6.702) gives -11.62 and -7.04. Spread between the two thermodynamic routes: 0.6 to 1.2 log units at 1500 to 2000 K; non-stoichiometric wuestite explains part of it. [R3378] labels IW "oxidised" and QFM "reduced" in the running text, which is backwards (IW is the reducing buffer); the generator should follow the physical order IW < QIF < QFM < NiNiO.

Meteorite-class outgassing redox, [R3944] (ApJ 843:120, Table 2; PDF read, arXiv id not recorded): `log10 fO2 = a + b 10^3/T + c 10^6/T^2 + d 10^9/T^3 + f 10^12/T^4` at 1 bar, 300 to 2000 K.

| Composition | a | b | c | d | f | Status |
|---|---|---|---|---|---|---|
| CI | 2.4976 | -9.8605 | -17.0701 | 7.5220 | -1.0404 | READ |
| CV | 9.0621 | -31.193 | 5.1092 | -1.8475 | 0.2000 | READ |
| H | 5.0743 | -22.906 | -5.6610 | 2.0634 | -0.2618 | READ |
| EH | 4.9495 | -24.024 | -4.6236 | 1.7177 | -0.2332 | READ |
| Eucrite | 5.4856 | -25.127 | -3.6580 | 1.3014 | -0.1650 | READ |

Redox ranking is CI > CV > H > EH; EH is about IW-1 and eucrite IW-1.4 to IW-0.8 (READ, same paper). The existing project redox draw (`REDOX_EARTH_DELTA_IW = 3.5`, slope 2.0 per dex of mass, sigma 1.0) is therefore far more oxidised than meteoritic building blocks (dIW about -1 to 0) and corresponds to a modern Earth-like upper mantle; [R4154] (NOT REACHED, snippet: Fe3+/sumFe 0.037 +/- 0.005, magma ocean near IW+0.5, CO-rich then CO2-N2) suggests a magma-ocean stage near IW to IW+1 before oxidation. See Boss list.

### 3.3 Speciation equations

Fits of log10 K(T) = a + b/T, 600 to 2000 K, bar standard state, from NASA polynomials (DATA) fitted by COMPUTED least squares:

| Reaction | a | b (K) | Max residual (log10) | Status |
|---|---|---|---|---|
| H2 + 1/2 O2 = H2O | -2.913 | 12952 | 0.04 | COMPUTED |
| CO + 1/2 O2 = CO2 | -4.514 | 14742 | 0.03 | COMPUTED |
| CO + 3 H2 = CH4 + H2O | -13.146 | 11697 | 0.07 | COMPUTED |
| H2S + 3/2 O2 = SO2 + H2O | -4.166 | 27156 | 0.01 | COMPUTED |
| N2 + 3 H2 = 2 NH3 | -12.234 | 5732 | 0.08 | COMPUTED |
| CO + H2S = COS + H2 | -1.477 | 4 | 0.06 | COMPUTED |
| CO + NH3 = HCN + H2O | 0.296 | -2742 | 0.02 | COMPUTED |

Cross-check against [R4863]: K(CH4 CO) = pCH4 pH2O/(pCO pH2^3) = 5.24e-14 exp(27285/T), atm (READ), versus my fit ln K = -30.27 + 26934/T in bar (COMPUTED); these agree to the stated accuracy once atm versus bar is accounted for.

Ortenzi's own relations: (H2O/H2)^2 = K3 fO2 and CO2/CO = K4 fO2^(1/2) (their Eq 16 and 18, READ); carbonate melt fraction X_CO3 = K1 K2 fO2 / (1 + K1 K2 fO2), with log10 K1 = 40.07639 - 2.53932e-2 T + 5.27096e-6 T^2 + 0.0267 (p - 1)/T and log10 K2 = -6.24763 - 282.56/T - 0.119242 (p - 1000)/T (READ; degassing at about 2000 K).

### 3.4 Speciation results (COMPUTED)

Solver inputs: ideal gas, element budget per H atom C = 0.1, N = 0.02, S = 0.02 (a design choice for illustration, not a sourced mantle inventory), species H2, H2O, CO, CO2, CH4, N2, NH3, H2S, SO2, COS, HCN, no graphite, no solubility.

| Condition | Result (mole fractions) | Status |
|---|---|---|
| 1200 K, 1 bar, dIW = 0 | H2 0.48, H2O 0.30, CO 0.11, CO2 0.051, H2S 0.033, CH4 1.7e-5, SO2 1.6e-7 | COMPUTED |
| 1200 K, 1 bar, dIW = +4 | H2O 0.78, CO2 0.16, SO2 2.1e-2, H2S 1.1e-2 | COMPUTED |
| 1000 K, 100 bar, dIW <= +1 | CH4 about 0.24; CO falls to 1e-6 to 1e-3 | COMPUTED |
| 800 K, 100 bar | NH3 about 1.5e-2 to 2e-2 | COMPUTED |
| H2S/SO2 crossover | near dIW = +3 at 1200 K | COMPUTED |

Direct comparison with the published bulk-silicate-Earth speciation ([R3945: Table 3], READ):

| Condition | Continental crust composition | BSE composition |
|---|---|---|
| 1000 K, 1e-6 bar | 49% H2O, 39% CO2, 4.5% HF | 61% H2O, 20% SO2, 12% CO2 |
| 1000 K, 100 bar | 55% H2O, 42% CO2, 1.4% HCl | 80% H2O, 16% CO2, 2.3% H2 |
| 2000 K, 1e-6 bar | 52% SiO, 25% O, 15% O2 | 28% Mg, 26% O, 24% SiO |
| 2000 K, 100 bar | 48% H2O, 38% CO2, 4% HF | 62% H2O, 20% SO2, 12% CO2 |

Further READ points from the same paper: H2 is never the major H-bearing gas; CO2 holds about 30% (molar) of the gas up to 1150, 1300, 1550, 1950 and 2300 K at 1e-6, 1e-4, 1e-2, 1 and 100 bar; N2 dominates at 500 to 700 K; NH3 appears for the BSE at low T; H2S is the main S gas for the BSE at low T. Spread note: my COMPUTED table at dIW = 0 is H2-rich, the Schaefer composition (a more oxidised, crust-like or BSE inventory) is H2O/CO2-rich; this is the same dIW dependence, not a conflict.

### 3.5 Generator recipe for rocky outgassed atmospheres

```
inputs: P_surf, T_surf, dIW, delivered volatile class (Low/High/Dry), star
1. fO2 = 10^(A - B/T + C (P-1)/T)   (Holloway buffers, 3.2) + dIW
2. H2O/H2 = sqrt(K_H2O fO2) ;  CO2/CO = K_CO2 sqrt(fO2)   (3.3)
3. If dIW < 0.5 and T < 900 K and P > 10 bar: CH4/CO via K(CH4 CO), N2->NH3 via K(NH3) (3.3)
4. S: H2S below dIW about +2.5 (T = 1200 K), SO2 above ~ +3  (3.4)
5. Ar from 40K decay (3.6); N2 as the remainder after volatile closure
6. Remove condensates: H2O below T_sat(P), NH3 etc. via section 12
```

Validity: 600 to 2000 K for K(T); ideal gas holds to about 100 bar for the CO2/H2O mixture and fails above that (use fugacity coefficients; NOT REACHED here).

### 3.6 Noble gases and radiogenic argon

Decay law for 40Ar: `Ar40(t) = f_Ar x K40_now x (exp(lambda t) - 1)`, with lambda = ln 2 / 1.248 Gyr = 0.5554 per Gyr (2ND-HAND: half-life 1.248 Gyr from NUBASE2020 via Wikipedia raw), branch fraction f_Ar = lambda_EC/lambda about 0.104 (10.3% plus about 0.1% plus 0.001%, 2ND-HAND from Wikipedia), 40K/K = 117 ppm (2ND-HAND, same). COMPUTED: Ar40 per K40 present is 0.077 at 1 Gyr, 0.212 at 2 Gyr, 1.162 at 4.5 Gyr; the K40 at 4.5 Gyr ago was 12.2 times the present amount.

TEST check (Earth; inputs RECALLED): BSE K about 240 ppm (McDonough and Sun), mantle mass 4.0e27 g gives about 3.3e18 mol of 40Ar produced in 4.5 Gyr (3.4e-6 mol per gram K, COMPUTED); the atmosphere holds about 1.7e18 mol of 40Ar (0.934% Ar at 99.6% 40Ar in 1.8e20 mol of air, RECALLED); the ratio is 0.5, consistent with the RECALLED 50 to 66% degassed fraction. Status: COMPUTED on RECALLED inputs.

Generator rule: Ar40_atm (mol) = f_degas x f_Ar x K40_now (mol) x (exp(lambda t) - 1), where K40_now = 117e-6 x total planetary K (mol) and f_degas is drawn between 0.05 and 0.7 depending on tectonic regime (the 0.5 to 0.66 Earth value is RECALLED). This is a design range anchored on one TEST planet; the real distribution for other planets is NOT REACHED. Primordial 36Ar and 38Ar plus Ne, Kr, Xe: solar ratios in section 2.1 apply for captured nebular gas (Ar/Ne = 0.021, Kr and Xe about 1e-9 and 1.7e-10 relative to H). For secondary atmospheres the noble gases are strongly depleted relative to solar; the depletion factor is RECALLED (Earth Ne 18.18 ppm, Kr 1.14 ppm, Xe 0.087 ppm in air, RECALLED TEST values).

### 3.7 Related outgassing results (context rows)

| Result | Source | Status |
|---|---|---|
| Volcanic CH4 fluxes comparable to biological ones are unlikely and would come with abundant CO | [R2490, R4729] | READ |
| 1 to 10 Earth-mass planets with 10 to 100 day periods keep water-dominated atmospheres; FeO reduction creates H2O that is soluble in magma | [R2383] | READ |
| CO2 atmospheres suppress H2O outgassing; a molten mantle increases radius by about 5% | [R0502] | READ |
| Atmodeller: CO-dominated at magma-ocean stage evolving to CO2-rich atmospheres of a few hundred bar at 280 K | [R0503] | READ |
| Abiotic O2 build-up near 5 bar per Myr on Earth-mass planets, up to 25 bar per Myr on super-Earths, reaching hundreds to thousands of bar (GJ 667Cc up to about 2000 bar, up to 10 oceans lost) | [R2769] as quoted in [R2286] | 2ND-HAND |
| Liggins classes at 800 K: Class O (oxidised: SO2, sulfur allotropes), Class I (intermediate: CO2, CH4, CO, COS), Class R (reduced: H2, NH3, CH4) | [R2686] | READ |
| H/C and age set whether the main species is H2, H2O or CO2 | [R2686] | READ |

Atmospheric class by mantle redox (generator mapping): dIW < 0 gives Class R (H2, CH4/CO, NH3); dIW 0 to +2 Class I; dIW > +2.5 Class O. The 0 and 2.5 cuts are my interpolation of Liggins and the section 3.4 crossovers and are labelled COMPUTED/design, not published numbers.

---

## 4. Lava worlds and silicate vapour

| Item | Value | Source | Status |
|---|---|---|---|
| Dominant gases above the BSE or komatiite melt: Na, O2, O, SiO, Fe, Mg | | [R3059] | READ |
| Of 193 Kepler candidates studied: 189 Type I, 3 Type II, 1 Type IV | | same | READ |
| Type boundaries in the text (Tp in K): I < 2033, II 2033 to 2588, III 2588 to 2890, IV 2890 to 3168, V > 3168 | | same text | READ |
| Type boundaries in the figure caption lines: 2096, 2460, 2735, 2974, 3170 K | | same figure captions | READ |
| Spread of the two boundary lists: up to 130 K on the I/II boundary (2033 vs 2096) and 128 K on II/III (2588 vs 2460) | | COMPUTED | COMPUTED |
| Vapour pressure of BSE melt: log10(Peq) = k1 exp(k2 Ts), 1400 < Ts < 2800 K, k1 = {-34.7, -60.1, -60.2}, k2 = {-0.00112, -0.00122, -0.00101} for small (< 0.01), intermediate (about 0.5), large (> 0.9) vaporised fraction | | [R2380] | READ |
| Peq about 2e-4 bar at 2000 K and 8e-3 bar at 2500 K (small fraction) | | evaluated from Kite eq. | COMPUTED |
| Atmosphere exchanges quickly when bulk-silicate FeO is low; pools are well mixed below about 2400 K substellar temperature | | [R2380] | READ |
| Main species Na, K, Fe, Si, SiO, O, O2; SiO UV absorption causes a thermal inversion | | [R2105] | READ |
| Tabulated Type I to V vapour compositions and the original silicate-vapour chemistry papers ([R3059, R3943] numeric tables) | | | NOT REACHED (only figures read; see Boss list) |
| TEST: 55 Cnc e has a volatile atmosphere likely rich in CO2 or CO; best fits are CO2-rich, surface pressure 0.01 to 100 bar with CO2 volume mixing ratio above 1e-5; a CO-only scenario is allowed | | Hu et al. 2024, [R2023] | READ (TEST) |

Recipe (lava world, Ts = substellar surface temperature, composition BSE):
```
logP_eq(Ts) = k1 * exp(k2 * Ts)       ([R2380], 1400-2800 K)
x_Na  : x_SiO : x_Mg : x_Fe : x_O2 : x_O  ~ taken from the Miguel/Schaefer table by Type (I-V)
    Type I (Tp < ~2033-2096 K): Na dominant
    Type II-III: SiO and O2 grow
    Type IV-V (Tp > ~2890 K): Fe, Mg, Si, SiO dominate
```
Mole-fraction values by type were not extracted in this lane (figures only); the [R3945] table in 3.4 (2000 K, 1e-6 bar: 28% Mg, 26% O, 24% SiO for the BSE) is the only quantitative anchor. See Boss list for tabulated Type I to V compositions.

Tracers for optics: SiO, Na, K, Mg, Fe, TiO, VO (the last two follow the Burrows and Sharp sequence in section 7: TiO exists above about 2000 to 2300 K, VO above about 1700 to 1900 K).

---

## 5. Sub-Neptunes, super-Earths, hycean and water worlds

| Item | Value | Source | Status |
|---|---|---|---|
| Four atmosphere classes for thick atmospheres: hydrogen-rich (X_H >= 0.7), water-rich, oxygen-rich, hydrocarbon-rich | | Hu and [R2025] | READ |
| Water-rich atmospheres require carbon as CO2 and exist only for C/O well below solar; loss of light elements from a C-rich, O-poor atmosphere leaves C2H2 and C2H4 | | Hu and [R2025] | READ |
| Hot Neptunes with metallicity 10 to 10,000 x solar span H2-dominated to CO2-, CO-, H2O- or O2-dominated compositions | | [R3176] | READ |
| Interior-atmosphere coupling: H about 1 to 2 wt% and O about 10 wt% of the metal core; about 0.1 wt% H in the mantle; more than 93% of H sits in the atmosphere; the envelope is oxidised, mantle and core reduced | | [R3961] | READ |
| Hycean: H2 envelope mass fraction <= 0.1%, minimum water mass fraction 10%, ocean-bottom T about 300 to 400 K, pressure at the habitable-surface boundary 1 to 1000 bar, radii up to 2.6 R_Earth at 10 Earth masses | | [R2804] | READ |
| Water mixing ratio in a hycean photosphere: H2O volume mixing ratio of 10% corresponds to the 10% mass fraction case | | same | READ |
| TEST: K2-18b (HST) H2O 0.033% to 8.9% (1 sigma), N2 < 10.9%, CO < 7.45%, CO2 < 2.4%, NH3 < 13.5%, CH4 < 0.248% (2 sigma) | | [R0340] | READ |
| TEST: K2-18b (JWST) CH4 about 1% and CO2 about 1% in an H2-rich atmosphere (5 sigma and 3 sigma), NH3 not detected, log X_CH4 = -2.04 (+0.61/-0.72) | | [R2805] | READ |
| Spread note: the HST CH4 upper limit (0.248%) lies below the JWST value (about 1%, log X = -2.04 gives 0.9%), so the two datasets disagree by a factor of about 4 at face value | | both above | COMPUTED |
| TEST: TOI-270 d mean molecular weight 5.47 (+1.25/-1.14), atmospheric metal mass fraction 58 (+8/-12)%, C/O 0.47 (+0.16/-0.19), interior about 90 wt% rock and iron; fully miscible supercritical metal-rich envelope | | [R0339] | READ |
| TEST: GJ 9827 d steam world, O/H about 4 by mass, H2O volume mixing ratio above 31% as the background gas | | [R3517] | READ |
| TEST: GJ 1214b reflective metal-rich atmosphere, metallicity above 100 x solar (models 1 to 3000 x solar considered) with thick haze | | [R2316] | READ |

Recipe (sub-Neptune or water world):
```
inputs: M, R (or bulk density), Teq, envelope mass fraction f_env, dIW, star
1. Draw envelope metal mass fraction Zenv: log-uniform in [0.01, 0.6]
   (bracketed by GJ 1214b >100 x solar, TOI-270 d 58 wt%, K2-18b near solar to tens of x)
2. mean molecular weight mu: Zenv = 0.01 -> ~2.4 ; 0.58 -> ~5.5 (TOI-270 d value) ; beyond -> H2O background (GJ 9827 d, >31% H2O)
3. Oxygen carrier: H2O set by H2 + 1/2 O2 = H2O at the magma/envelope interface fO2 (3.2-3.3)
   (Schlichting and Young: envelope oxidised even though mantle is reduced)
4. C: CH4 below ~500-700 K, CO2/CO above; N: N2 (NH3 not seen in K2-18b)
5. If hycean branch (f_env <= 0.001, ocean mass >= 10%): H2O vapour limited by T_surface sat. pressure (12.2)
```
The Zenv range is a design bracket from TEST planets, not a fitted distribution; an occurrence-weighted distribution is NOT REACHED (see Boss list).

---

## 6. Evolved secondary and frost-equilibrium atmospheres (Venus-, Mars-, Titan-, Pluto-like classes)

These are class recipes. The Solar System bodies are TEST rows only, used to confirm that a recipe can reproduce a known endpoint.

### 6.1 Recipes

```
Oxidised terrestrial (CO2-N2 type, Venus/Mars-like):
   take the outgassed mixture from 3.5 at high dIW (CO2, H2O, N2, SO2) ; remove H2O by escape/condensation ;
   N2 from nitrogen inventory ; Ar from 3.6 ;   O2 from H2O photolysis + escape (Luger and Barnes: up to ~5 bar/Myr)
Reduced N2-CH4 (Titan-like, T < ~150 K, N2 dominant):
   N2 from NH3 photolysis or outgassing ; CH4 mole fraction set by surface-ocean / interior reservoir ; H2 from CH4 photolysis (escape-limited)
Frost-equilibrium (Pluto/Triton-like): p_i = chi_i * p_sat,i(T_surface)   (ideal solution, Raoult)
   N2 sets total P ; CH4 and CO ride on top at the ice-mole-fraction scaled ratio
```

### 6.2 TEST rows (Solar System calibration only)

| TEST body | Values | Source | Status |
|---|---|---|---|
| Titan | stratospheric CH4 1.48 +/- 0.09%, near-surface CH4 5.65 +/- 0.18%, H2 (1.01 +/- 0.16)e-3 | [R3289] (NTRS file) | READ |
| Mars | CO2 0.960, Ar 0.0193, N2 0.0189, O2 1.45e-3, CO below 1e-3 | no source | READ |
| Mars (later) | CO2 0.951, N2 0.0259, Ar 0.0194, O2 1.61e-3, CO 5.8e-4 | no source | SNIPPET |
| Mars spread between the two sets | N2 0.0189 vs 0.0259 (30% apart) | both above | COMPUTED |
| Venus | CO2 96.5%, N2 below 3.5%, SO2 150 ppm, Ar 70 ppm, H2O 20 ppm, CO 17 ppm | search summary (Wikipedia-type compilation) | SNIPPET |
| Pluto | N2 above 99%, CH4 about 0.30%, CO about 0.05% (cited); surface CH4 (3 to 3.6)e-3 in N2 ice; near-surface CO/N2 about 4e-3 and CO/CH4 about 1.7e-3 reported by one set of authors | [R3186] review, [R3186] | 2ND-HAND |
| Triton | CO/N2 about 0.015% in the atmosphere (ground-based surface ice gave about 0.1%, predicted atmospheric mixing ratio 1.5e-4); CH4 partial pressure several times larger than Voyager; CO abundance similar to surface abundance, controlled by a thin CO-rich veneer | [R2642] (abstract READ); [R3186] | READ and 2ND-HAND |

Two values for Pluto CO (about 0.05% atmospheric; CO/N2 about 4e-3 near surface) are not the same quantity and lie a factor of 8 apart; the generator should state whether it means a column-integrated or a near-surface ratio. The vapour-pressure relations for N2, CO, CH4 and CO2 ices are NOT REACHED ([R1445] not available; see 12.2 and the Boss list), so the frost-equilibrium recipe can only be used with the project's existing two-point table until those relations are obtained.

---

## 7. Gas giants and brown dwarfs

### 7.1 Condensation and chemistry sequence versus temperature (cloud-free, solar composition)

[R0608] (astro-ph/9807055) sequence, READ:

| Temperature (K) | Event |
|---|---|
| 2300 to 2000 | TiO disappears; silicates and Fe appear |
| 1900 | Mg2SiO4 forms |
| 1900 to 1700 | VO disappears |
| 1700 | MgSiO3 forms |
| 1400 | CrH disappears; Li goes to LiCl below this |
| 1500 to 1200 | CO converts to CH4 (pressure dependent) |
| below 1200 | K, Rb, Cs go to chlorides |
| 1100 to 700 | Sulfide, NaCl and KCl clouds |
| about 700 (600 to 800) | N2 converts to NH3 |
| about 350 | water cloud |
| about 200 | NH3 and NH4SH clouds |

CH4 is the main carbon carrier except in the youngest and hottest objects (Teff >= 1300 K) (READ, same paper). The transition temperatures are pressure dependent; the generator must use the equations below at the local P.

### 7.2 Condensation temperatures (T in K, P_T in bar, [X/H] in dex)

[R4561] and 2006 (astro-ph/0511136), READ:

| Condensate | Relation |
|---|---|
| Fe | 10^4/T_cond = 5.44 - 0.48 log P_T - 0.48 [Fe/H] |
| Mg2SiO4 | 10^4/T_cond = 5.89 - 0.37 log P_T - 0.73 [Fe/H] |
| MgSiO3 | 10^4/T_cond = 6.26 - 0.35 log P_T - 0.70 [Fe/H] |
| MnS | 10^4/T_cond = 7.45 - 0.42 (log P_T + [Mn/H] + [S/H]) |
| Na2S | 10^4/T_cond = 10.05 - 0.72 (log P_T + [Na/H] + 0.5 [S/H]) |
| ZnS | 10^4/T_cond = 12.52 - 0.63 (log P_T + [Zn/H] + [S/H]) |
| NH4SH | 10^4/T_cond = 48.91 - 4.15 (log P_T + 0.5 [N/H] + 0.5 [S/H]) |
| H2S (solid) | 10^4/T_cond = 86.49 - 8.54 log P_T - 8.54 [S/H] |
| Fe / Fe(OH)2 boundary | log P_T = 11.94 - 12088/T - 2 [Fe/H] |
| FeH gas | log X = -6.13 - 1905/T + 0.5 log P_T + [Fe/H] |

Quartz forms at about 1550 K and enstatite at about 1600 K (READ). Monatomic Fe, Mg and SiO are the dominant refractory gases above the clouds (READ). Validity: the fits are for solar-like H2/He gas at the pressures of brown-dwarf and giant atmospheres (roughly 1e-3 to 1e2 bar); their metallicity range is stated in the paper and is not extrapolated here.

### 7.3 CO / CH4 / NH3 equilibrium and quench

[R4863], READ:
- Equilibrium: K(CH4 CO) = pCH4 pH2O / (pCO pH2^3) = 5.24e-14 exp(27285/T), atm.
- Quench times: t_q1 = 1.5e-6 p^-1 m^-0.7 exp(42000/T) s, t_q2 = 40 p^-2 exp(25000/T) s, 1/t_CO = 1/t_q1 + 1/t_q2, with p in bar and m the metallicity in solar units; the quench level is where t_mix = H^2 / Kzz equals t_CO.
- The CH4 to CO transition lies near 500 K for planets (up to 600 K with sluggish mixing, 400 K with clouds or strong mixing) and 1100 +/- 100 K for brown dwarfs.
- The NH3/N2 ratio is insensitive to mixing in their grid (Teff 500 to 1100 K).

### 7.4 Mixing ratios by temperature regime

Median and range columns are COMPUTED from the solar budget in 2.2 scaled by metallicity m (CH4, H2O, NH3, H2S scale with m while H2 and He are weakly changed), then corrected by the condensation relations above; the log10 range column combines m = 0.3 to 30 with the C/O range 0.3 to 1.0 (design brackets from [R2801, R2803] and Line 2021). They are not measured values.

| Class (Teff or Teq, K) | H2 | He | H2O | CH4 | CO | NH3 | H2S | Na+K | Controls | Status |
|---|---|---|---|---|---|---|---|---|---|---|
| Cold giant, ice-giant envelope (Teq < 200) | 0.86 | 0.14 | cloud-depleted (0 to 8e-4 m) | 5e-4 m (log range -4.0 to -1.5) | ~0 | cloud-depleted | cloud-depleted | ~0 (condensed) | m, cloud depth, interior T | COMPUTED |
| Warm giant (200 to 700) | 0.86 | 0.14 | 8e-4 m | 5e-4 m | quench-set, 1e-9 to 1e-5 | 1e-4 m | 2e-5 m | alkali chlorides/sulfides depleted | m, C/O, Kzz | COMPUTED |
| T dwarf (700 to 1300) | 0.86 | 0.14 | 8e-4 m | 5e-4 m | 1e-6 to 1e-4 (quench) | 1e-4 m | 2e-5 m | Na, K present above clouds | Kzz, gravity | COMPUTED |
| L dwarf (1300 to 2200) | 0.86 | 0.14 | 8e-4 m | trace | 5e-4 m | N2 | 2e-5 m | Na, K strong; FeH, CrH | silicate and Fe clouds | COMPUTED |
| Y dwarf (< 500) | 0.86 | 0.14 | cloud-depleted | 5e-4 m | quench-set | 1e-4 m | cloud-depleted | ~0 | water cloud near 350 K | COMPUTED |

Recipe:
```
given Teff (or T(P)), g, m, C/O, Kzz:
   x_H2 = 0.857*(1 - Zscale), x_He = 0.141*(1 - Zscale)   [Zscale = metal mole fraction]
   O budget = m * 4.9e-4 minus the O taken by silicate/Fe clouds (~20%, [R4561] rock-forming elements)
   C budget = m * 2.9e-4 * (C/O)/(0.589) ; CH4/CO by K(CH4 CO) with the quench of 7.3
   N budget: NH3 at T < 700 K (Burrows and Sharp), N2 above ; S budget: H2S [R4562]
   condense Fe, Mg2SiO4, MgSiO3, Na2S, NH4SH, H2O, H2S by 7.2 and 12.2 before computing gas fractions
```
The ice-giant interior (above 1 kbar, water-rich) has no sourced mixing ratio in this lane; only the envelope recipe above applies. A mass-radius-consistent water/ammonia/methane interior mixture is NOT REACHED (Boss list). TEST: Jupiter-class retrievals follow log(CH4/H) = -1.11 log(M/MJ) + 0.38 (2.5), which the Solar System giants satisfy by construction.

### 7.5 Brown dwarf spreads

Zahnle and Marley give the CH4 to CO crossover temperature as 1100 +/- 100 K for brown dwarfs and 500 K (400 to 600 K) for planets, a 600 K spread between classes at the same chemistry because of surface gravity and Kzz. The generator should not interpolate between the two numbers; it should evaluate the quench formula in 7.3 with the object's own gravity, metallicity and Kzz.

---

## 8. Hot and ultra-hot Jupiters

| Item | Value | Source | Status |
|---|---|---|---|
| For HD 189733b and HD 209458b, CH4, NH3, HCN and C2H2 are enhanced by quenching and photochemistry; CO, H2O, N2 and CO2 follow equilibrium except below about 1 microbar | | [R3177] | READ |
| Photochemical and kinetic models for hot Jupiters (HD 209458b, HD 189733b) | not extracted | [R4527] | NOT REACHED (PDF downloaded, body not read) |
| TiO/VO as thermal-inversion absorbers: cold-trapping depletes them unless Kzz is at least a few times 1e7 cm2/s (0.1 micron grains) up to about 1e11 cm2/s (10 micron grains); VO unlikely to matter | | [R4172] | READ |
| Ultra-hot dayside (Tday about 2900 K): H- opacity, thermal dissociation of H2O and other major molecules, and thermal ionisation of metals must be included; WASP-18b [M/H] = -0.01 +/- 0.35, C/O < 0.85 (3 sigma), inverted temperature profile | | [R0134] | READ |
| H2O in hot giants is sub-stellar while Na and K are stellar to super-stellar | | [R4663] | READ |
| C/O >= 1 gives HCN and C2H2, with TiO and VO underabundant | | [R2801, R2803] | READ |
| Teq of TiO survival: Burrows and Sharp: TiO exists above 2000 to 2300 K, VO above 1700 to 1900 K | | astro-ph/9807055 | READ |

Recipe:
```
given Teq, Tday, m, C/O, Kzz:
   Tday < 1500 K : as warm giant (7.4), with CH4/NH3/HCN enhanced by Kzz quench
   1500 - 2200 K : CO + H2O + H2 ; CO2 rises with m ; TiO/VO only if cold trap leaky (Kzz > 1e8 and a < 1 micron)
   Tday > 2200 K : thermally dissociate H2O (to OH + H, O) , H2 (to H) ; add H-, Na, K, Fe, Mg, SiO, TiO (Section 4 list)
Mole-fraction values for dissociation depend on P and T and are not tabulated in this lane (T-P-dependent; use FastChem/GGchem data, NOT REACHED here).
```

---

## 9. Comets and comae

### 9.1 Abundances relative to H2O

| Species | 67P percent of H2O (TEST) | Cross-comet range, percent of H2O | Source | Status |
|---|---|---|---|---|
| CO2 | 4.7 +/- 1.4 | median about 12 (25 comets; about half CO2-dominated, a third CO-dominated, a tenth comparable) | [R3835]; Harrington [R3535] | READ; SNIPPET |
| CO | 3.1 +/- 0.9 | 0.5 to 23 (Hale-Bopp 23.0, Hyakutake 18.2); median about 3 | [R3835, R4704] (Tables 2 and 3); Harrington [R3535] | READ; READ; SNIPPET |
| O2 | 3.1 +/- 1.1 | not tabulated | [R3835] | READ |
| CH4 | 0.34 +/- 0.07 | not extracted | [R3835] | READ |
| C2H6 | 0.29 +/- 0.06 | not extracted | [R3835] | READ |
| CH3OH | 0.21 +/- 0.06 | 0.2 to 2.6 | [R3835, R4704] | READ |
| H2CO | 0.32 +/- 0.10 | not extracted | [R3835] | READ |
| NH3 | 0.67 +/- 0.20 | 0.1 to 0.73 | [R3835, R4704] | READ |
| HCN | 0.14 +/- 0.04 | 0.09 to 0.29 | [R3835, R4704] | READ |
| H2S | 1.10 +/- 0.46 | not extracted | [R3835] | READ |
| N2 | 0.089 +/- 0.024 | not extracted | [R3835] | READ |
| OCS | 0.041 | not extracted | [R3835] | READ |
| SO2 | 0.127 | not extracted | [R3835] | READ |
| Ar | 5.8e-4 | not extracted | [R3835] | READ |
| Kr | 4.9e-5 | not extracted | [R3835] | READ |
| Xe | 2.4e-5 | not extracted | [R3835] | READ |
| Ne | < 5e-6 | not extracted | [R3835] | READ |

The 67P CO2 value (4.7) and the 25-comet CO2 median (about 12) differ by a factor 2.6 and the 67P CO (3.1) sits at the cross-comet median (about 3). Spread for CO alone across comets is a factor of 46 (0.5 to 23). Generator rule: draw `CO2/H2O` log-normal with median 0.12 and sigma about 0.5 dex, and `CO/H2O` log-normal with median 0.03 and sigma about 0.6 dex (sigma values are COMPUTED by eye from the quoted ranges, not published fits). Per-comet tabulations of [R0048] (NTRS, OCR damaged) and [R3206] / [R0436] are NOT REACHED; a protoplanetary-disk midplane comparison is in [R1234] and an ice-origin review in [R3833] (both downloaded, tables not extracted).

### 9.2 Production-rate scaling (COMPUTED)

Energy balance of a sublimating water-ice surface: `Q(r) = (1 - A) S0 pi R^2 / (r^2 L)` molecules per second, with L the latent heat per molecule, with S0 = 1361 W/m2 at 1 AU, A = 0.04, L = 51 kJ/mol (RECALLED, to verify), r in AU. This gives a flux limit of 1.54e22 molecules/s/m2 for the subsolar point at 1 AU and 4.8e28 molecules/s for R = 1 km. Values by distance (flux per m2): 6.2e22 at 0.5 AU, 1.5e22 at 1, 3.9e21 at 2, 1.7e21 at 3, 6.2e20 at 5. The r^-2 law holds inside about 2 to 2.5 AU; beyond it radiative loss and the water ice line (about 2.5 AU, RECALLED) matter and activity is driven by CO and CO2 (RECALLED, consistent with CO and CO2 dominance in the Harrington Pinto sample above). An empirical production-rate law (Combi, Meech and Svoren family) is NOT REACHED.

---

## 10. White dwarf and pulsar planets

| Item | Value | Source | Status |
|---|---|---|---|
| Debris accreted by polluted white dwarfs is volatile-depleted and rock-forming (O, Mg, Si, Ca, Fe dominated); the best extrapolation from calcium to total heavy-element accretion uses the bulk terrestrial Ca mass fraction A about 0.016 | | [R1325] | READ |
| Upper limit on the H mass fraction of the debris is 1%, consistent with volatile-depleted parent bodies; carbon is consistent with the terrestrial value, not chondritic | | [R1325] (citing Hartmann et al.) | 2ND-HAND |
| Six polluted white dwarfs give fO2 consistent with Earth, Mars and typical asteroids; most oxidised bodies about dIW = -0.6, WD 1536+520 dIW = -1.37; Solar System bodies span about 6 log units in dIW | | [R1154] | READ |
| Planets that survive post main sequence evolution and their orbital fates | | [R4532] | NOT REACHED (PDF downloaded, body not read) |
| Atmospheres of planets around white dwarfs and pulsars (composition, escape) | | none | NOT REACHED |

Recipe: for an exposed rocky body around a white dwarf, composition is bulk-silicate-Earth like (3.4 table, 2000 K row) with dIW drawn from N(-1, 0.5) (spread bracketed by Doyle values), volatiles near zero (H < 1% by mass of debris), so only a thin outgassed or vaporised atmosphere (Na, SiO, Mg, Fe) is expected when heated. For pulsar planets no sourced composition exists; treat as bare rock with vapour only if heated, and flag as a gap.

---

## 11. Host-star and M-dwarf dependence

| Finding | Source | Status |
|---|---|---|
| On Earth-like planets around M dwarfs the biogenic gases CH4, N2O and CH3Cl have substantially longer lifetimes and higher mixing ratios than on Earth (UV spectrum differs from the Sun); active M-star planets still develop an Earth-like ozone layer | [R3999], astro-ph/0510224 | READ (abstract) |
| N2O builds up to observable concentration around M dwarfs with low UV flux; CH3Cl may become detectable; the O2 0.76 micron feature is harder to detect for later M dwarfs | [R3845, R3846], [R3845] | READ (abstract) |
| Abiotic O3 appears without detectable O2 in prebiotic-Earth-like atmospheres | [R1137] | READ (abstract) |
| A catalytic H2O2 cycle sustains CO2; below about 1 ppm H2, about 40% of CO2 converts to CO and O2 in about 1 Myr; the loss is capped near 50% | [R1511] | READ |
| F and G stars do not accumulate abiotic O2; K and especially M stars can | [R1814] | READ |
| Extreme water loss in the M dwarf habitable zone gives abiotic O2 build-up (hundreds to thousands of bar in the cases quoted by [R2286]) | [R2769] | READ (title); numbers 2ND-HAND via [R2286] |

Recipe: the host star changes the atmosphere mainly through the FUV/UV to visible ratio and flare history. For an M dwarf host and a rocky planet with an oxidised mantle, move oxygen from the CO2/H2O reservoir to O2 on a timescale set by water loss (up to about 5 bar/Myr on Earth-mass planets, 2ND-HAND), and move CO2 to CO plus O2 by the factor in the Gao row (cap 50%) when H2 is below 1 ppm. Quantitative flare-driven rates and an F/G/K/M tabulation of mixing ratios are NOT REACHED beyond these abstracts.

---

## 12. Vapour pressure and condensation relations that were reached

### 12.1 Water ice (COMPUTED from DATA)

From NASA condensed thermodynamic data for H2O(s) and H2O gas (cantera 3.2.0 `nasa_condensed.yaml` and `nasa_gas.yaml`, range 200 to 273.15 K): `log10 p_sat(bar) = 7.5622 - 2670.16/T`, maximum fit residual 6e-4 in log10. Check values: 6.118e-3 bar at 273.15 K (triple point 6.117e-3 bar, RECALLED, agrees), 2.73e-4 bar at 240 K, 1.63e-6 bar at 200 K; extrapolated below 200 K to 5.4e-8 bar at 180 K and 6.1e-11 bar at 150 K (outside the data range, flag EXTRAPOLATED).

### 12.2 NH3 and beta-O2 ices (READ, COMPUTED check)

[R0418]: `ln p_vap(bar) = A - B/T + C ln T + D2 T + D3 T^2 + D4 T^3`.

| Species | A | B (K) | C | D2 (1/K) | D3 (1/K^2) | D4 (1/K^3) |
|---|---|---|---|---|---|---|
| NH3 ice | -5.55 +/- 0.45 | 3605 +/- 31 | 4.82792 | -0.024895 | 2.1669e-5 | -2.3575e-8 |
| beta-O2 ice | 15.29 +/- 0.29 | 1166.2 +/- 8.6 | -0.75587 | 0.14188 | -1.8665e-3 | 7.582e-6 |

COMPUTED check: NH3 gives 6.17e-2 bar at 195 K (triple point about 6.1e-2 bar, RECALLED) and 3.7e-2 bar at 190 K, 3.4e-3 bar at 170 K, 1.6e-4 bar at 150 K, 2.5e-7 bar at 120 K; O2 gives 2.7e-3 bar at 54 K, 4.9e-4 bar at 50 K, 3.7e-5 bar at 45 K, 1.4e-6 bar at 40 K. Their experimental vapour pressures lie below the [R1445] extrapolation by a factor 1.5 for NH3 and 1.6 for beta-O2 (READ), a direct measure of the spread to carry. The project's two-point table (`CONDENSATION_K` in `atmosphere.py`) can be validated against these.

### 12.3 Not reached

N2, CO, CH4, CO2, SO2, H2S and NH4SH(s) sublimation relations ([R1445, R1708] updates) are NOT REACHED; only the NH4SH and H2S condensation relations of 7.2 exist.

---

## 13. Species by class (final table)

Entries are mole fractions (gas phase, dry unless noted) with the evidence status of the entry. Column key: RR = reduced rocky outgassing (dIW <= +1, 800 to 1200 K, 1 to 100 bar); RO = oxidised rocky (dIW about +4); LV = lava-world vapour (2000 K or more); SN = H2-rich sub-Neptune or hycean; WW = steam or water world; GG = giant 200 to 1300 K; HJ = hot or ultra-hot Jupiter; CM = comet coma (percent of H2O, 67P TEST). A dash means negligible by chemistry in that regime; "?" means no sourced number.

| Species | RR | RO | LV | SN | WW | GG | HJ | CM (%H2O) | Status |
|---|---|---|---|---|---|---|---|---|---|
| H2 | 0.48 | - | ? | ~0.8 | 0.02 to 0.7 (H2O background) | 0.86 | ~0.9 | - | COMPUTED/READ |
| He | - | - | - | 0.14 | - | 0.14 | 0.14 | - | COMPUTED |
| N2 | 0.01 (1/3 N budget) | main | - | ? | ? | - | - | 0.089 | COMPUTED/READ |
| O2 | trace | 0.15 (crust, 2000 K, 1e-6 bar); abiotic build-up possible | 0.15 | - | ? | - | - | 3.1 | READ |
| O3 | - | trace, abiotic possible | - | - | - | - | - | - | READ (abstract) |
| CO2 | 0.05 | 0.16 (COMPUTED), 0.12 to 0.42 (BSE/crust, Schaefer) | ? | ~0.01 | up to dominant | - | rises with m | 4.7 | COMPUTED/READ |
| CO | 0.11 | ~0 | ? | < 7e-2 (K2-18b limit) | ? | quench 1e-9 to 1e-4 | 5e-4 m | 3.1 | COMPUTED/READ |
| CH4 | 1.7e-5 (1200 K), 0.24 (1000 K, 100 bar) | - | - | ~0.01 | ? | 5e-4 m | quench-enhanced | 0.34 | COMPUTED/READ |
| H2O | 0.30 | 0.78 | trace | 0.03 to 0.09 (K2-18b bound) | > 0.31 | 8e-4 m (cloud-limited) | dissociated above ~2200 K | 100 | COMPUTED/READ |
| NH3 | 1.5e-2 to 2e-2 (800 K, 100 bar) | - | - | not detected | ? | 1.2e-4 m | quench-enhanced | 0.67 | COMPUTED/READ |
| H2S | 0.033 | 0.011 | ? | ? | ? | 2.3e-5 m | ? | 1.10 | COMPUTED/READ |
| SO2 | 1.6e-7 | 0.021 | ? | photochemical | ? | - | photochemical (WASP-39b TEST) | 0.127 | COMPUTED/READ |
| Ar | 40Ar radiogenic, section 3.6 | same | - | 4.1e-6 (solar) | ? | 4.1e-6 | 4.1e-6 | 5.8e-4 | COMPUTED/READ |
| Ne | - | - | - | 2.0e-4 (solar) | ? | 2.0e-4 | 2.0e-4 | < 5e-6 | COMPUTED/READ |
| Kr | - | - | - | 1.3e-9 per H (solar) | ? | same | same | 4.9e-5 | COMPUTED/READ |
| Xe | - | - | - | 1.7e-10 per H (solar) | ? | same | same | 2.4e-5 | COMPUTED/READ |
| NH4SH | - | - | - | - | - | cloud (7.2) | - | - | READ |
| HCN | via CO + NH3 (3.3) | - | - | ? | ? | trace | quench- and C/O-enhanced | 0.14 | READ |
| C2H6 | - | - | - | ? | ? | trace | ? | 0.29 | READ |
| C2H2 | - | - | - | ? | ? | trace | C/O >= 1 | ? | READ |
| OCS (COS) | via CO + H2S (3.3) | - | - | ? | ? | ? | ? | 0.041 | READ |
| HCl | 1.4% (crust composition, 1000 K, 100 bar) | ? | ? | ? | ? | KCl, NaCl condensates | ? | ? | READ |
| HF | 4.5% (crust composition, 1000 K, 1e-6 bar) | ? | ? | ? | ? | - | ? | ? | READ |
| N2O | - | abiotic ? | - | - | - | - | - | - | NOT REACHED |
| NO2 | - | lightning-driven, ? [R1813] | - | - | - | - | - | - | NOT REACHED |
| Na | - | - | Type I dominant | - | - | condensed below ~1100 K | 1.7e-6 per H | - | READ/COMPUTED |
| K | - | - | present | - | - | KCl below ~1200 K | 1.2e-7 per H | - | READ/COMPUTED |
| SiO | - | - | 0.24 (BSE at 2000 K, 1e-6 bar) | - | - | - | refractory gas | - | READ |
| Mg | - | - | 0.28 (BSE at 2000 K, 1e-6 bar) | - | - | - | refractory gas | - | READ |
| Fe | - | - | present | - | - | cloud (7.2) | refractory gas | - | READ |
| TiO | - | - | ? | - | - | - | cold-trapped unless Kzz large | - | READ |
| VO | - | - | ? | - | - | - | same (less likely) | - | READ |
| PH3 | - | - | - | - | - | disequilibrium | - | - | NOT REACHED (qualitative only) |
| CH3OH | - | - | - | - | - | - | - | 0.21 | READ |
| H2CO | - | - | - | - | - | - | - | 0.32 | READ |

---

## 14. Boss must research

1. [R1445] (Planet. Space Sci. 57) and [R1708] for N2, CO, CH4, CO2, SO2, H2S, NH3 sublimation coefficients: needed for the Pluto/Triton frost-equilibrium recipe and validation of `CONDENSATION_K`. Not reachable here (journal; no arXiv copy found).
2. [R4154] (Sci. Adv.) and [R1491] and 2021: magma-ocean redox and CO to CO2-N2 evolution; snippet only. Needed to fix the dIW prior of the generator.
3. [R1953] (graphite, C solubility) and carbon solubility in melt: our speciation has no graphite and no solubility, so CH4/CO at the surface is biased.
4. [R0436, R3206], Dello [R1063], Harrington [R3535] full tables: per-comet abundances of 30 or more species; [R4704] has only ranges. Also [R0048] with usable OCR.
5. [R2725] and 2006, [R2726]: brown dwarf chemistry and the original condensation tables; [R3943] for silicate vapour.
6. [R3059] numeric tables (Type I to V compositions) for lava worlds; only figures were read.
7. The [R4352] scatter on the Zpl/Zstar fit and the width of the atmosphere-to-bulk offset in giants.
8. An occurrence-weighted distribution of sub-Neptune envelope metallicity (Zenv) for the generator; used here only as a bracket of TEST planets.
9. Hot Jupiter thermochemistry tables (FastChem or GGchem output) for T above 2200 K, and [R4527] / Line 2021 numbers.
10. Mass-radius-consistent water, ammonia and methane interior mixture of ice giants (no mixing ratio found).
11. Pulsar planet and white-dwarf planet atmospheres: no sourced composition.
12. Fugacity coefficients for CO2-H2O at above 100 bar (ideal gas fails).
13. Combi, Meech and Svoren production-rate law for comets and the CO, CO2 activity onset distance.
14. N2O, NO2, O3 abiotic abundances ([R1137, R1813] numeric values), and the Atmospheric Toxicity document item [R1813] (ApJ 866, 56, real).

## 15. Recalled / to verify

All values below are RECALLED and not verified:
- [R1464] iron-wuestite fit log fO2 = -27489/T + 6.702.
- BSE K about 240 ppm (McDonough and Sun); K/U about 12,000 (Arevalo); Earth 40Ar degassed fraction about 50 to 66%.
- Earth air: Ne 18.18 ppm, Kr 1.14 ppm, Xe 0.087 ppm, Ar 0.934%; 40Ar 99.6% of Ar; atmosphere 1.8e20 mol; mantle mass 4.0e27 g.
- Latent heat of water ice sublimation 51 kJ/mol; water ice line near 2.5 AU; triple points of H2O (6.117 mbar) and NH3 (about 61 mbar).
- Marti and Mauersberger ice vapour pressure fit (not used in the tables).
- Design widths 0.3 dex (atmospheric metallicity), C/O 0.55 x 10^N(0, 0.15), CO2/H2O and CO/H2O sigma values, dIW prior N(-1, 0.5) for white dwarf debris.
- Mars N2 0.0259, Venus compositions (SNIPPET status; keep as TEST only).

## 16. Cross-links

- Part on thermodynamics and gas tables: use the same K(T) fits (3.3) and vapour pressures (12) for density and scale height.
- Part on clouds and hazes: condensation relations 7.2 and 12; NH4SH, Fe, silicate and water clouds.
- Part on photochemistry and escape: SO2 and O3 photochemistry (7, 11), abiotic O2, H2 loss, Jeans escape for Ar/Ne/He fractionation.
- Part on radiation and shielding: M-dwarf UV (11).
- Project code: `src/planetgen/physics/atmosphere.py` `SPECIES`, `CLASS_MIXES`, `CONDENSATION_K`, `REDOX_*` constants; the report `composition-density-link-2026-10-17.md` lists the balance-gas and aerosol gaps this part is meant to close. The redox prior (dIW = 3.5 + 2 log10(M/M_Earth)) is more oxidised than the meteoritic and magma-ocean values in 3.2 and should be revisited.
