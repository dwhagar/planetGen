# Atmosphere colour: one consistent atmosphere, then its colour

**Status (2026-10-17, drafted 12:19 PM EDT, revised 12:50 PM EDT for Boss's papers and Boss's conflict rules of 12:39 PM EDT):** the checked design (GEN.409), Boss approved section 18 as recommended (9:40 to 9:55 PM EDT). Nothing is built, and this note changes no code, stored value, schema or generator epoch; the build items GEN.410 to GEN.419 do (section 13). It replaces the plan note of 2026-10-16. Boss (11:24 PM EDT and 11:25 PM EDT) made the colour of the air of any world, worked out from the atmosphere it really has, the seminal feature of 9.0. Sources: the master [atmosphere-composition.md](atmosphere-composition.md) and its parts (cited "part 01" to "part 10"), [atmosphere-sources-check.md](atmosphere-sources-check.md) (GEN.408, "sources check"), [sky-colour-analysis.md](sky-colour-analysis.md), [composition-model.md](composition-model.md) sections 4a and 4b, [tunable-gaps.md](tunable-gaps.md) section 6, [boss-papers-2026-10-17.md](../research/boss-papers-2026-10-17.md) (on origin/main; sections 2.1 to 2.4, 2.6, 2.7, 3, 4; a paper wins over a second-hand value in our files and both values are stated), and the code of section 3.

**Verification states** (Boss, 3:02 AM EDT: verified and partly verified are usable unless contradicted; the verified source wins): **verified** (read, data file, or computed; wins every contradiction); **partly verified** (second-hand, or true inside a stated limit; usable inside it); **open** (no usable source; the entry names the seeded physical range used meanwhile, never a fixed value); **stand-in** (a measured Solar System shape or an unvalidated engineering formula, a flagged temporary row with the item that replaces it). Recalled or snippet values are placeholders, never sources. **Conflict rules (Boss, 12:39 PM EDT):** when in doubt verified wins; when two verified sources are in doubt the newest wins; when two sources are both in question the most permissive wins. Third rule, Decided (Boss, 2026-10-17, 1:12 PM EDT): the widest value or range either source allows, never a value outside both. Each place that settles a disagreement carries 'Rule applied' and section 17 lists them.

## 1. What Boss decided

- Store the albedo and the greenhouse values used in generation; the migration works the formula back from stored values, or draws from the planet class where that fails.
- Break the molar density into per-gas partial pressures and calculate it from them, so both match and the light meets a known mixture.
- Build ozone in from real research; if the research is missing, a generator from the oxygen and the star's UV; Boss researches on request.
- Sky colour and atmosphere colour are the same; the wording is "atmosphere colour" and it replaces the planet colour where that is what it meant.
- The colour is absolute, not adapted to the star (card, evening of 2026-10-16); gas giants are in the first version.
- Haze, dust and similar scatterers get reasonable draws or derivations.
- One table of every gas and compound, including every compound the descriptions mention.
- Nothing is fixed per class or copied from the Solar System; every input is a stored physical quantity or a seeded draw. Solar System bodies are test targets.
- Later rulings: value-changing chain stays in 8.4 (12:13 AM EDT, 2026-10-17); radius reference level 1 bar (1:24 AM EDT); gated possibility tables (1:40 to 1:51 AM EDT); the usability rule (3:02 AM EDT); Tunable gaps defaults D4 and D5 (3:36 AM EDT relay).

## 2. The chain

- **8.4 phase 1:** GEN.408 sources check (done); GEN.412 albedo and greenhouse recorded, with the inverse.
- **8.1:** GEN.409 this design; GEN.410 gas table; GEN.472 gated table loader; GEN.470, GEN.497 saturation curves; GEN.496 photolysis lifetimes.
- **8.4 phase 2:** GEN.411 molar density from partials; GEN.452 condensation curves; GEN.422 no recalculation redraw; GEN.423 density, scale height, greenhouse follow the mixture; GEN.424 named balance gas; GEN.454 giant mixtures; GEN.455 exospheres; DB.67, DB.71.
- **8.4 phase 3:** GEN.413 ozone; GEN.414 haze, dust, cloud; GEN.456 lava-world vapour; GEN.462 ground albedo by band.
- **GEN.462 wiring (ground under the air):** `body_colour.ground_albedo` hands the area-weighted band reflectance of the stored surface mix (GEN.441, `surface_composition.band_albedo`) to the colour model in place of the 0.15 stand-in; the open ocean, the metal and a body with no stored mix keep the stand-in for their share. The Bond albedo for the temperature and the space weathering of an airless surface are GEN.684.
- **8.4 phase 4:** DB.68; GEN.415 text from the mixture; TEST.209 sweep; DOC.63.
- **9.0:** phase 1 GEN.416 star colour, GEN.417 Rayleigh depth, GEN.418 solid-world colour, GEN.419 giant colour; phase 3 GEN.420, MAP.244, DOC.64; phase 4 GEN.421 headline.

The foundations sit in 8.4 because the classes are reworked there (GEN.90). VIEW.11 extends GEN.416 and stays in 9.1. The trace of 2026-10-17 found that only the pressure flows from the stored mixture; the drawn molar density, scale height and multiplier never see it, the recalculation redraws albedo and mix, the freeze-out lowers the pressure after the temperature is set, and habitability uses 28 g/mol for every balance gas (GEN.411, 422, 423, 424).

## 3. What exists today

- `physics/atmosphere.py`: ten stored species (line 38; `p_<gas>_kpa` columns, `db/models.py` 765 and 1246); a mole-fraction template per class (`CLASS_MIXES`, 86 to 105; C, D, S have no air); helium, sodium, ammonia are unstored balance gases; the cold side caps each partial at a two-point curve (`CONDENSATION_K`, 77 to 84, the H2S row copies SO2) and lowers `atmospheric_pressure` after the temperature was set (241 to 245).
- `physics/planets.py` `calculate_atmospheric_conditions` (630 to 716): albedo from the class `albedo_range` (E, F, G, H, K, L, M, N, O, P, V) or (0.12, 0.35) (A, B, C, D, I, J, Q, S, T; 657, 658); T_e the airless temperature (662); scale height R T_e / (M g) with `atm_molar_density` (671 to 677, stored as `scale_height_km`); pressure effective density times g times H (692 to 694); greenhouse factor min(500, `atm_molar_density` / 0.04345 times a class multiplier) (710 to 713); T_s = max(2.725 K, ((1 - A) S (1 + G) / (4 sigma))^(1/4)) (714, 715). Albedo and multiplier are discarded. `atm_molar_density` is a mean molar mass in kg/mol drawn independently of the mixture. A moon passes its parent's distance (`generation/validation.py` 287, 715).
- `physics/radiation.py` 216, 223: `ozone_ratio` (square root of pO2, cap 1.2) and a `uv_index` that is not the CIE index. No aerosol quantity, no colour code.
- Giants (classes I, J, T) run the same code; `constants.py` 143 to 162 already defines their air density at a 1 bar level.
- Database: schema v106 (migration 0106), no albedo or greenhouse column; generator epoch 11 (`src/tests/golden/epochs.json`, first release 8.0.1708).

## 4. The decisions

| # | Decision | State |
|---|---|---|
| 1 | Store albedo and greenhouse; inverse for old rows, class draw where it fails (section 8) | Decided, Boss, 2026-10-16, 11:24 to 11:25 PM EDT (the plan note does not time each line) |
| 2 | Partial pressures first; molar density calculated from them (sections 6, 7) | Decided, same |
| 3 | Ozone from research, else a generator from O2 and star UV (section 9) | Decided, same |
| 4 | "Atmosphere colour" replaces "sky colour" and the planet colour | Decided, same |
| 5 | Absolute colour, no adaptation; giants included | Decided, card of the evening of 2026-10-16 (minute not recorded) |
| 6 | One gas and compound table with every compound the descriptions mention | Decided, 11:24 to 11:25 PM EDT |
| 7 | Nothing fixed per class, no presets; Solar System rows live in tests only (TEST.215) | Decided, standing rule; 1:43 AM EDT on stand-ins |
| 8 | Value-changing chain stays in 8.4, outside the 8.1 stamp | Decided, 12:13 AM EDT, 2026-10-17 |
| 9 | Radius reference level 1 bar | Decided, 1:24 AM EDT |
| 10 | Gated possibility tables, stand-in marker, one shared component table | Decided, 1:40 to 1:51 AM EDT |
| 11 | Verified wins; partly verified usable inside its limit | Decided, 3:02 AM EDT |
| 12 | Life-linked gases (N2O, PH3, NO2, CH3Cl) default off, never tuned; Kzz and haze yield are nuisance parameters | Default taken, 3:36 AM EDT relay, stands unless Boss overrules |
| 13 | Stored greenhouse quantity is G = (T_s / T_e)^4 - 1, column `greenhouse_factor` (the code's class draw `greenhouse_multiplier` stays internal) | Default taken (master finding 5) |
| 14 | Old rows: albedo and G from the stored scale height, never the pressure; failures take a seeded class draw; `albedo_basis` flag | Default taken |
| 15 | Stored versus derived as in section 5 | Default taken |
| 16 | `atm_molar_density` keeps name and unit (kg/mol); new rows hold the sum of x_i M_i; old rows are rewritten after the inverse has used the old value | Default taken |
| 17 | Gas table columns and add list (section 6) | Default taken |
| 18 | Rayleigh data per gas (section 6.2): refractivity values centre N2, O2, Ar, CO2 (Bates King factors, Bucholtz fit as check); measured 532 nm gaps of Sneep and Ubachs set the per-gas error (2.5 to 18 percent); generic 5 to 9 percent only for SO2, H2S, H2, He, NH3, H2O | Partly verified; Bogaard 1978 not uploaded, so those six stay open. Rule applied: verified, then newest, then most permissive (17) |
| 19 | Ozone: Kozakis polynomial by host times a seeded structure factor LU[0.85, 1.71]; abiotic branch below 0.01 percent O2; column stored | Default taken; host T_eff map open. Rule applied: verified, widest model spread (17) |
| 20 | Cloud draw from the condensation level and vapour budget (10.1) | Default taken; efficiency, cover open |
| 21 | Haze: threshold, diffusion-limited supply times yield LU[0.01, 1] (10.2) | Stand-in; absolute depth open |
| 22 | Dust gated by the lofting threshold, depth LU[0.5, 1.5] meanwhile; albedo spectrum from Wolff 2009 (digitized from Figures 12 to 14) (10.3) | Stand-in; no Wolff table, Titan haze unread. Rule applied: verified (17) |
| 23 | Giant reference level P_ref = min(1 bar, cloud top) (section 11) | Default taken |
| 24 | Giant colour derived from optical columns; Sudarsky classes and the Cahoy table are labels and test targets | Default taken; band strengths open |
| 25 | Single scattering to tau(550 nm) 0.3, two-stream above (supersedes the plan note's 0.5) | Default taken; verified basis. Rule applied: verified (17) |
| 26 | Planck spectrum from 4500 K; cooler stars use BT-Settl chromaticity with a seeded spread | Default taken; licence open |
| 27 | Output: hue swatch (brightest channel 1), CIE xy, sky-to-star luminance ratio, a word; none stored | Default taken |
| 28 | Three migration revisions, numbers taken at build time | Default taken |
| 29 | This design bumps no epoch; the build items bump per phase | Default taken |
| 30 | The open inputs of section 16 run on their seeded ranges | Open |

## 5. Stored versus derived

| Quantity | Kind | Note |
|---|---|---|
| Partial pressures (ten columns; added gases per DB.71) | stored | the mixture is the cause; the sum equals the stored pressure |
| `atmospheric_pressure_pa` | stored | after the freeze-out |
| `atm_molar_density`, kg/mol | stored, derived at generation | sum of x_i M_i (decision 16); the scale height and old readers use it |
| `bond_albedo`, `greenhouse_factor` | stored | draws; not recoverable from the new pressure; the temperature is exactly what they give |
| `albedo_basis` (0 generated, 1 exact inverse, 2 class draw) | stored | the migration writes 1 and 2 |
| `surface_temperature`, `scale_height_km`, `atm_density` | stored (exist) | follow from the above (GEN.423) |
| `ozone_column_du` | stored | the seeded structure factor makes it not derivable |
| `tau_haze`, `tau_dust`, `tau_cloud` (550 nm), `cloud_cover`, `aerosol_draw` | stored | `aerosol_draw` is one seeded unit number that picks every within-bracket optical property |
| T_e, Rayleigh depth, ozone absorption, air mass, aerosol material and constants, giant reference level, Sudarsky label | derived | from stored values and the tables |
| Atmosphere colour (hex, xy, luminance ratio, word) | derived | a pure function; the star's colour from its stored temperature |

## 6. The gas and compound table (GEN.410) and the mixture

### 6.1 Columns

One row per species; each column group carries its own verification state, source and licence (the GEN.472 loader refuses a row without them).

| Group | Columns | Source |
|---|---|---|
| Identity | key (formula, never a class letter), name, formula, CAS, molar mass (`periodictable`), isotopologue flag, flags (infrared, volatility class, near-UV fragile, toxic, `column_only` for O3) | part 01 T1 |
| Phase | triple and critical points, boiling or sublimation point, saturation curve (id, coefficients, validity window, blend point), latent heats | part 01 sections 3, 6, 7; part 06 |
| Thermal | heat capacity at three temperatures, gamma, viscosity, Henry constant, critical constants for real-gas Z | part 01 sections 4, 8 to 10 |
| Rayleigh | refractivity formula and the temperature it is defined at (288.15 or 273.15 K), or polarizability and depolarization; sigma(550 nm) check; error spread | sources check 2.1, 2.2; part 02 2.1 |
| Absorption | UV and visible bins with temperature trend; O3 Chappuis, Hartley, Huggins; alkali and oxide data for hot giants | part 02 2.3 to 2.5 |
| Infrared | band centres, strength class, forcing law (Myhre 1998, 6.3), collision-induced pairs | part 03 section 2; Myhre Table 3 |
| Photochemistry, escape | photolysis threshold, lifetime by host star, reaction rows, escape data | part 04; GEN.496 |
| High energy | EUV and X-ray absorption, attenuation, proton stopping, W value | part 07 |
| Distribution | the gated possibility rows (GEN.472): median mole fraction and log range by gate; comet rows per 6.3 | parts 05, 09; Mumma and Charnley 2011 |

Condensates (water, ammonia, silicates, salts, tholins, dust minerals) are a second registry: density, saturation curve, n and k by wavelength, single scattering albedo, exponent (part 06 2.8). The saturation curves replace `CONDENSATION_K` (GEN.452): the Wagner liquid fit above the triple point, the Fray and Schmitt sublimation fit below, and the Grundy 2024 and Blakley 2024 laboratory fits only inside their windows (N2 23 to 31 K, CO 26 to 35 K, CH4 31 to 43 K, O2 29 to 37 K, NH3 88 to 105 K; part 01 7.1, 7.3). The code's H2S row is 20 to 50 K wrong.

**Species.** The ten stored species stay. The add list is He, Ne, Kr, Xe, NH3, NH4SH, HCN, C2H6, C2H2, OCS, HCl, HF, N2O, NO2, O3, Na, K, SiO, Mg, Fe, TiO, VO, PH3 (master finding 10), plus the class-text words with no row today: "mineral dust" (an aerosol row, 10.3), "trace elements" (lava-world vapour rows), "hydrocarbons" (CH4, C2H6, C2H2, C2H4).

### 6.2 Rayleigh data

sigma = 24 pi^3 (n^2 - 1)^2 F_K / (lambda^4 N^2 (n^2 + 2)^2), N the density at which (n - 1) is defined; F_K = (6 + 3 rho) / (6 - 7 rho), so rho = 6 (F - 1) / (3 + 7 F) (King 1923 Eq. 78, 90; Sneep and Ubachs Eq. 7; read). Mixtures add by number fraction, equivalently F_mix - 1 = sum w_r (F_r - 1) with weights v_r (mu_r - 1)^2 (King Eq. 91, 94; read). The King factors of N2 and O2 are Bates 1984 (N2 F = 1.034 + 3.17e-4 / lambda^2; O2 F = 1.096 + 1.385e-3 / lambda^2 + 1.448e-4 / lambda^4, lambda in micron; read from the paper's page images). Bucholtz's depolarization column is back-computed from Bates's rounded factors, so Bates and Bucholtz are one source, not two (recomputed in Boss's papers 2.1).

**Air shortcut and check.** Bucholtz 1995 Table 3 gives sigma = A lambda^-(B + C lambda + D / lambda), lambda in micron, cm^2: A 3.01577e-28 below 0.5 micron and 4.01061e-28 above; B 3.55212 and 3.99668; C 1.35579 and 1.10298e-3; D 0.11563 and 2.71393e-2; stated accuracy under 0.4 percent at 0.2 micron and under 0.1 percent above 0.5 micron (read). It is the regression check for Earth-like air at 15 C. At 550 nm the sigma of air is 4.509e-27 cm^2 (Bucholtz Table 2), 4.513e-27 (Bates Table 1, 0 C), 4.507e-27 (the fit, computed here) and 4.5135e-27 (this note's per-gas sum, computed here).

**Per-gas cross sections at 550 nm,** cm^2 (ratio to N2), from the refractivity form with the factors above: N2 4.635e-27; O2 4.072e-27 (0.879); Ar 3.991e-27 (0.861); CO2 1.164e-26 (2.51). Gases that carry a literature ratio: CH4 9.905e-27 (2.14, He 2021); N2O 1.587e-26 (3.42, He 2021); CO 5.63e-27 (1.21, the Sneep and Ubachs measured ratio; Allen and Cox form gave 1.24, partly verified, dropped); H2 9.38e-28 (0.202, Dalgarno and Williams, partly verified) to 1.008e-27 (0.217, Peck and Huang constants, data), so 0.20 to 0.22; He 6.0e-29 (0.013); NH3 7.2e-27 to 7.4e-27 (1.55 to 1.60, partly verified); Ne 2.19e-28, Kr 9.21e-27, Xe 2.40e-26 (data). H2S 1.98e-26 and SO2 2.41e-26 are the polarizability route. H2O vapour has no known reference temperature: generic 3.38e-27 with a spread to the n-based 2.8e-27 to 3.1e-27. None of the five Rayleigh papers holds a depolarization ratio for SO2, H2S, H2, He, NH3 or water; Bogaard, Buckingham, Pierens and White 1978 is not uploaded, so these six stay open at the generic 5 to 9 percent error (a 5 percent log-normal draw, 13 gases, ratios 0.907 to 1.060; part 02 3.1).

**Measured error per gas.** Sneep and Ubachs 2005 Table 2 (read from page images, 532.2 nm, 15 C, 1013 hPa, 1e-27 cm^2); the gap is the refractivity value over the measured value:

| Gas | measured | refractivity-based | gap | measured ratio to N2 |
|---|---|---|---|---|
| Ar | 4.45 (0.3) | 4.56 | 2.5 percent | 0.87 |
| N2 | 5.10 (0.24) | 5.30 | 3.9 percent | 1 |
| CO2 | 12.4 (0.8) | 13.29 | 7.2 percent | 2.43 |
| CO | 6.19 (0.4) | 6.82 | 10.2 percent | 1.21 |
| N2O | 15.90 (0.08) | 18.19 | 14.4 percent | 3.12 |
| CH4 | 12.47 (0.23) | 14.69 | 17.8 percent | 2.45 |
| SF6 | 32.3 (0.5) | 34.1 | 5.6 percent | 6.33 |

This note's formulas give 5.30 (N2), 4.56 (Ar) and 13.32 (CO2) at 532.2 nm (computed here), the refractivity column. The paper's own statement is that the refractivity method holds "at the few percent level", with larger gaps for CO, N2O and CH4. Rule applied: verified (Sneep Table 2 replaces the generic error). Decision 18 takes the refractivity value as the centre for N2, O2, Ar and CO2, because it is the value behind the published Earth target (14.1, 15) and the measured gap is inside the stated few percent, and takes the full gap down to the measured value as the seeded range (one-sided log-uniform factor from measured / refractivity to 1; tests use the centre). Rule applied: verified for the centre, most permissive for the range (17). For CH4, N2O and CO2 the ratios come from He 2021 (2.14, 3.42, 2.51) where both are measured and Sneep and Ubachs give 2.45, 3.12, 2.43: two verified measurements, newest (He 2021) centres, and the older one bounds the draw. Rule applied: newest, range most permissive (17). The coded Sneep and Ubachs CO form and the older CH4 ratio 2.76 (the refractivity ratio, 2.77 in Table 2) are dropped.

Phase function p = 3 / (4 (1 + 2g)) [(1 + 3g) + (1 - g) cos^2 theta], g = rho / (2 - rho) (Monte Carlo verified, part 02 2.8).

### 6.3 Forcing laws, comets and mantle redox (Boss's papers)

- **Forcing (Myhre et al. 1998 Table 3, printed p. 2718, read).** dF = 5.35 ln(C/C0) W m^-2 for CO2 (C in ppmv; IPCC 1990 had 6.3); CH4 0.036 (sqrt M - sqrt M0) minus an overlap term; N2O 0.12 (IPCC 1990 0.14) likewise; overlap f(M,N) = 0.47 ln[1 + 2.01e-5 (MN)^0.75 + 5.31e-15 M (MN)^1.52]; CFC-11 0.25, CFC-12 0.33. Assumed uncertainties: 1 percent CO2, 10 percent CH4, 5 percent N2O. CO2 doubling is 3.71 W m^-2 (5.35 ln 2 = 3.708). The "IPCC 1990" labels on 5.35 and 0.12 in the parts file are wrong: they are Myhre. The paper states no validity range, so the laws are Earth-range only; correlated-k grids stay a download job. Stored G (section 8) is not computed from them; they are a test target.
- **Comets (Mumma and Charnley 2011 Table 4, read).** NH3 is 1.5, 0.5, 0.7 and below 0.3 percent of water in the four comets that have it, so the repo bound of 3.6 percent is unsupported: about 1.5 percent is the upper end. SO2 is 0.2 percent of water in Hale-Bopp, a second comet beside 67P (0.127 percent). N2, O2, Ar, Kr and Xe in comets have no source in the review and stay open.
- **Mantle redox (Gaillard and Scaillet 2014; Pahlevan, Schaefer and Hirschmann).** No mass law and no later-mantle law (gap rows AT.a1, a2): they stay stand-ins. The sulfur outgassing rule (10.4) and any H2 against H2O outgassing branch inherit that flag. Gaillard's text says FMQ-2 +- 0.5 while its Figure 2 shows about FMQ-2.75 with 0.7 dex spread at 0.01 bar; neither is used.

## 7. Mixture, pressure and molar density agree (tests in TEST.209)

(i) Partial pressures sum to the stored pressure to 1e-9 relative. (ii) The molar mass is the sum of x_i M_i. (iii) No partial exceeds the saturation curve at the surface temperature unless the body names a liquid or solid reservoir. (iv) The text names only table gases above the GEN.415 thresholds (`atmosphere.describe`: the four largest gases from 1 percent of the air, then the two largest traces from 0.01 percent, each under the gas table's name; a class carries no gas phrase, `tuning.PLANET_CLASSES[...]["air"]` says only whether it has air, and `atmosphere.PENDING` stands in until the air is solved; tested for every class and zone in `test_atmosphere_text_consistency.py`). (v) T_s is the formula of the stored A and G. (vi) The Rayleigh depth reads the column P / (m g) with the mixture's mean molar mass (master finding 5), never the drawn one.

## 8. Albedo and greenhouse: the inverse for old rows

The old generator discards two draws but its stored scale height keeps T_e. With H = `scale_height_km` times 1000, M = `atm_molar_density`, g = `gravity` times 9.80665, S = L / (4 pi d^2) (stored luminosity and distance; a moon uses its parent's distance):

    T_e = H M g / R        1 - A = 4 sigma T_e^4 / S        G = (T_s / T_e)^4 - 1        (airless: T_e = T_s)

This is exact (14.2). The pressure route is not: the freeze-out lowers the pressure after the temperature is set, so a body that lost 20 percent of its air gives T_e of 204 K and A of 0.71 where the truth is 0.29. Part 03 6.3 is right that one equation has two unknowns for a physical greenhouse model; the old generator has a third stored relation, so its rows are solvable.

**It fails (take the class draw) when:** the scale height, molar mass, gravity, luminosity or distance is missing or zero for a body with air; T_e or T_s sits on the 2.725 K floor (the equation becomes an inequality); A falls outside [0, 1) or G outside [0, 500]; or the stored pressure exceeds `atm_density` times g times H times the retention factor (gravity in g) by more than 1e-6, which marks a row not built by this code or edited after generation. The ratio of G to the molar-mass ratio against the class multiplier range is logged, not enforced, because class ranges changed between releases.

**As built (GEN.412).** `planetgen/physics/climate_draws.py` holds the formulas (`equilibrium_temperature_k`, `surface_temperature_k`, `greenhouse_from_temperature`), the exact inverse (`from_stored`, which returns the reason when it has no solution) and the seeded class draw (`class_draw`, keyed on `albedo-backfill:<table>:<id>`). A new body records `bond_albedo`, `greenhouse_factor` (the capped G, 0 with no air) and `albedo_basis` 0; migration 0142 (schema v142) adds the columns to `planets` and `moons` and fills older rows with basis 1 or 2 by the rules below (the moon's starlight is its planet's, a circumbinary planet's is the pair's effective luminosity). In 200 sample systems every body inverted exactly. The columns are recorded draws, so the generator makes the same draws as before and takes no epoch.

**As built (GEN.423, DB.67).** Both items were met by the pieces above: `atmosphere_state.air_state` gives the molar mass `sum x_i M_i`, the density `P M / (R T)` and the scale height `R T / (M g)` from the stored partial pressures, `atmosphere_state.greenhouse_warming` gives `G` from the same mixture at the absorbed starlight (`surface_air.solve` settles the two together), and a body stores the albedo, `G` and the molar density it was built with (migration 0142 for the first two, the older `atm_molar_density` column for the third). Nothing in generation or the schema was left to change, so there is no epoch and no migration for them; `test_surface_air.py::test_every_generated_body_stores_the_density_scale_height_and_greenhouse_its_mixture_gives` holds the stored values to the mixture on generated planets, giants and moons (agreement to 1e-6, found at 1e-15).

**Class draw.** A is uniform on the class `albedo_range` (0.12 to 0.35 where none), limited to A at least 1 - 4 sigma T_s^4 / S so that G is not negative, from the one random wrapper keyed on the body seed and the label `albedo-backfill`; G comes from the stored T_s, so the stored temperature never moves. Giants have no class range: the fallback is open, bracket 0.02 to 0.91 (Sudarsky 2000 whole sample about 0.02 to 0.8; Cahoy 2010 Bond up to 0.914). New rows record the A and G they were built with (GEN.412); once aerosol and Rayleigh draws exist, A should follow from them and the ground by A = R + (1 - R)^2 A_s / (1 - R A_s) (part 03 7.2, about 20 percent accurate, partly verified), a question for Boss (section 18, item 3).

## 9. Ozone (GEN.413)

**Built (GEN.413, first slice):** `physics/ozone.py` has the median, the host mix, the F-star scaling, the structure spread, the abiotic branch; `ozone_column_du` is stored on planets and moons (migration 0135) and `radiation.uv_index` and `ozone_loss_flag` read it (the `ozone_ratio` square-root proxy is gone). The oxygen comes from `p_o2_kpa` until lane I's mixture accessor (GEN.411) lands. The Chappuis band is a table of cell means in `physics/data/ozone_chappuis.json`, read by `physics/atmosphere_colour.py` (section 12, GEN.418). Not built yet: the ozone as a stored gas of the mixture, and the GEN.189 loss events reducing the column (the flag is kept). **GEN.668 (removed by TEST.524):** a helper that scaled the stored column by the share of the air a loss event keeps was built and nothing called it, since no generation step removes air after the ozone is drawn (escape and the water loss act on the air before it, and the final-age refresh, GEN.839, works the dose and the ozone again from the air that is left); it was deleted.

**Model.** Stored in Dobson units (1 DU = 2.687e16 cm^-2). With x = log10 of the O2 mixing ratio at 1 bar and Earth gravity (the fit's variable), the median is the Kozakis 2022 polynomial N = a x^5 + b x^4 + c x^3 + d x^2 + e x + f (sources check 4.2, Table 3, read; valid 0.01 to 150 percent of the present O2 level, clamped outside). Evaluated here, the Sun host gives 183, 192, 210, 199, 108 DU at 150, 100, 25, 10, 1 percent (the paper's 5.64e18 cm^-2 = 210 DU maximum at 25 percent is reproduced); G0V 263 and 169 DU at 100 and 1 percent; K2V 160 and 63; K5V 179 and 38; M5V 162 and 17. Hot hosts peak at lower O2, cooler hosts rise to the highest. For other pressures and gravities the O2 column p_O2 / g replaces the mixing ratio in x (stand-in until the Chapman layer depth, part 04 3.3, replaces it).

**Spread.** The sources check's plus or minus 30 percent around the fit cannot reach the published columns: the Sun host gives 192 DU at 100 percent, the 3-D model 279 DU (Cooke 2022, read), 1-D models 296 to 327 DU (Grenfell 2007) and 313 DU (Segura 2005); at 10 percent the ratios to the fit run 0.85 (Cooke) to 1.34 (Segura 2003). Decision 19 multiplied the median by a structure factor log-uniform on [0.85, 1.71], which spans every sampled model and level; GEN.748 narrows it where the models agree: the range at an oxygen level is the lowest to the highest ratio the models give there (1.45 to 1.70 at 100 percent, 0.85 to 1.34 at 10 percent, 0.85 to 1.71 at 1 percent and below, linear in log PAL between; `ozone.STRUCTURE_KNOTS`, stand-in OZ3); row m1 of part 09 (LU[192, 279] DU) is too narrow. The Earth test range 279 to 300 DU lies inside the support (163 to 328 DU) and is a test, not a preset.

**Host map (open).** The fit has five hosts (G0V, Sun, K2V, K5V, M5V); their T_eff are in Kozakis Table 1, not in the notes. Once read, the column interpolates on the stored T_eff between the nearest two; F stars have no fit (Rugheimer 2015: 8.6e18 and 1.15e19 cm^-2 for F7V and F0V, read) and use that bracket; cooler than M5V keeps the M5V row.

**Without O2.** Below the fit's 0.01 percent limit (part 09 gates at 0.1 percent; the fit is verified to 0.01 percent, so the gate moves there; rule applied: verified, section 17) the abiotic branch gives 0.01 to 0.1 percent of the same host's 100 percent column for G, K, M hosts and 15 percent for an F host (Domagal-Goldman 2014 Table 3, read via Tunable gaps); "about 10 percent of Earth" is contradicted and dropped. Abiotic O2 outcomes are scenario flags (part 04 3.6).

**Optics.** Chappuis absorption, 410 to 760 nm, from the Burkholder and Talukdar 1994 files (Mainz, CC BY 4.0 snapshot; peak 5.163e-21 cm^2 at 602 nm; 220 K to 298 K ratios 0.68 at 420 nm, 0.84 at 450, about 1.0 from 500 to 650), tabulated at build time. It is Beer-Lambert absorption along the sun path through a thin shell at the Chapman peak pressure (part 04 3.3: 63 hPa at 100 percent O2, 630 hPa at 10 percent); at 300 DU the vertical depth is 0.042 at 600 nm. It matters at twilight, not at noon (5 to 8 percent of luminance, xy shift 0.005 to 0.01). Hartley and Huggins serve the UV index (GEN.445) from part 02's coarse bins until the Bremen licence (none stated) is settled. The `radiation.py` proxies are replaced.

## 10. Haze, dust and cloud (GEN.414)

Depths are at 550 nm; the wavelength law and colour come from the material (part 06 2.8); the gated rows are n1 to n9 of part 09.

**10.1 Cloud.** A condensate forms where x_i P(z) reaches the checked saturation curve along the temperature profile (part 06 3.3). The in-cloud depth is tau_c = 3 Q_ext m_c / (4 rho_p r_eff), Q_ext about 2 above 1 micron (part 06 3.5, read), m_c limited by the vapour above the condensation level; the draw is log-uniform from 0.1 to that limit, r_eff open (1 to 20 micron; Earth model droplets 10 to 17.5 micron, Groot 2020). Comparison depths: Earth bins 0, 5, 10, 20 (grey); Venus total 20 to 30 at 0.63 micron (second-hand); giant water clouds up to about 10,000 (Cahoy 2010). `cloud_cover` is 0 without a condensate, else a seeded draw on [0.1, 1] (open). Clouds are grey; the plane albedo is the Monte Carlo table of part 06 3.10 (g = 0.85, omega = 1: 0.419 at tau 1 and mu0 0.2, 0.602 at tau 10 and mu0 0.5); an absorbing cloud saturates (0.60 at omega 0.99).

**10.2 Haze (stand-in).** Onset at CH4/CO2 about 0.1, thick above 0.2, production cut about tenfold in H2-rich gas (master finding 28). Mass flux is the diffusion-limited CH4 supply N K_zz f_CH4 / H_a times a yield LU[0.01, 1] (laboratory rates span 250 times); the steady depth is 3 Q_ext Mdot / (4 rho_p r_eff v_f) (engineering, not validated, part 06 3.7). K_zz is a nuisance parameter (decision 12). The exponent is 0.97 to 2.34 for aggregates (3 to 4 only for tiny spheres); colour from tholin or soot constants (Khare k falls from 0.057 to 0.012 between 450 and 650 nm). Absolute depth is open until Tomasko 2008 and Doose 2016 are read. The confirmed Titan laws (above 80 km 1.012e7 lambda^-2.339, 30 to 80 km 2.029e4 lambda^-1.409, below 30 km 6.27e2 lambda^-0.9706, lambda in nm; total 9.6 at 500 nm; team slides) are a test target, not a preset.

**10.3 Dust (stand-in).** Dust exists on a dry surface when the wind friction velocity exceeds the lofting threshold (Shao and Lu or Iversen and White; the two differ by 15 to 60 percent, draw inside the span; part 06 3.8). The wind has no source: u*/u*_t is log-uniform on [0.3, 3] (open). Depth meanwhile LU[0.5, 1.5] (sources check 8; Mars-based, hence stand-in), grey (exponent near 0; 7 to 11 percent less opacity at 440 than at 880 nm), r_eff 1.2 to 1.8 micron, asymmetry 0.673 +- 0.081. Single scattering albedo: Wolff 2009 prints no table of albedo, n or k. It prints the solar-band w0 0.922 to 0.926 (+-0.003) for effective radius 1.2 to 1.8 micron (Table 3, p. 12, read; 0.940 for the 1.8 index, 1.4 scattering case; literature 0.89 to 0.94). The spectra at 30 channels from 440 to 2920 nm exist only in Figures 12 to 14 and are **digitized from Figures 12 to 14** (recovered from the vector geometry, accuracy about 0.0005 to 0.001 in w0; never called the paper's table; reader file A section 6.5). Digitized w0 for radius 1.2 and 1.8 micron: 0.765 and 0.772 at 440 nm; 0.780 and 0.788 at 460; 0.825 and 0.831 at 500 (k 0.0087 and 0.0055); 0.941 and 0.945 at 600; 0.968 and 0.971 at 700 nm; about 0.965 to 0.972 from 800 to 2480 nm. At 450 nm the digitized value is 0.772 to 0.780. This replaces the second-hand range 0.74 to 0.89 for the 450 nm value (the UV band 0.67 to 0.74 is below the digitized range and stays partly verified). Rule applied: verified (17). Stand-in until a table is read; Ockert-Bell 1997 and Tomasko 1999 Table 8 are not uploaded. `aerosol_draw` picks the radius in [1.2, 1.8] micron, which interpolates between the two columns. Hematite-bearing dust absorbs blue in every dataset; the hematite share comes from the ground composition (GEN.441, GEN.462). The Fe3+ attribution has no source and is not used.

**10.4 Sulfur aerosol.** Row n9: S8 or H2SO4 when the outgassed sulfur flux is at least 100 times the Earth volcanic flux in N2-CO2 air, or at the Earth flux in H2 air (Hu 2013, read); the acid share is partly open. Tabazadeh et al. 1997 gives the water side only: Eq. 1, ln P0 = 18.452406985 - 3505.1578807/T - 330918.55082/T^2 + 12725068.262/T^3 (mbar, 185 to 260 K, supercooled water), and Table 1, P_H2O = exp(a + b/T + c/T^2) for 10 to 80 weight percent H2SO4 (read; test at 200 K, 60 weight percent: P0 = 3.165e-3 mbar, a_w = 0.0735, recomputed). It gives no acid vapour pressure. Zeleznik 1991 gives the Gibbs energy of the binary solution (200 to 350 K, 0.1 MPa, Table 6A) and activities in Tables 7 and 8, and no partial pressures (it holds no vapour phase). The two disagree on water activity by 10 to 22 percent at 60 to 70 weight percent and 200 to 220 K and by under 5 percent at 40 to 50 percent and at 250 K. Rule applied: newest and built for the stratospheric range, Tabazadeh 1997 for water, Zeleznik as the bracket (17). The acid vapour pressure needs the pure-acid curve of Ayers 1980 (not uploaded) or the Kulmala-Laaksonen form in our files, times the acid activity (stand-in). Tabazadeh Table 2 branches are discontinuous at a_w 0.05 and 0.85, and Zeleznik Eq. 12 (ln x_j) and Eq. 11 (a minus sign) are misprinted.

## 11. Gas and ice giants: the reference level (GEN.419)

A giant has no ground, so the colour is read at a pressure level. The radius level is 1 bar (decision 9) and the code already defines giant air density there. Decision 23: P_ref = min(1 bar, cloud top). The cloud base is the first crossing of T(P) with the condensation curve going up (part 06 3.3), T(P) from the Guillot 2010 semi-grey profile (part 03 5.1; opacity link open). Condensate decays upward as exp(-f_sed z / L) (Ackerman and Marley 2001, read) while air thins as exp(-z / H), so the pressure where the cloud depth from the top reaches 1 is P_top = P_base tau_c^(-1 / (1 + f_sed H / L)) (algebra derived here; H/L up to 10, f_sed 0.01 to 10, second-hand). A cloud-free hot giant keeps 1 bar. The sky is the Rayleigh column above P_ref (0.089 at 550 nm for 1 bar of 86/14 H2/He at 24.8 m/s^2, part 03 7.1), the absorbers, and the cloud or deeper gas as the lower boundary.

Sudarsky classes (I below about 150 K, II about 250 K, III above about 350 K, IV 900 to 1500 K, V above about 1500 K; approximate; the 2003 paper uses T_eq 350 to 800 K for III) are a text label only. The colour comes from Rayleigh, condensate optics and cloud albedo (part 06), CH4, NH3 and H2O bands, and Na and K wings (S = 2e-14 and 3.5e-14 cm^2 nm). Above the Ti condensation curve (about 1700 K at 1 mbar, Fortney 2008) TiO and VO make a giant dark red-brown, 1400 to 1700 K is mixed, below it alkali wings and water rule (part 02 3.7). Band positions are read (Cahoy Table 3: CH4 0.46 to 0.99 micron, Na 0.59, K 0.40 and 0.77, H2O 0.94); strengths are open (Karkoschka 1994; Karkoschka and Tomasko 2009, 2010), so methane-bearing giants carry `colour_quality = partial`. Test targets: about 4 percent methane with a thin aerosol is the blue of an ice giant, none is pearly white (Irwin 2022, read); the Cahoy albedos of section 15.

**Built (GEN.419 first slice; `physics/giant_colour.py`):** `giant_atmosphere_colour` cuts the stored mixture at `P_ref` and runs the section 12 colour model with the cloud deck's diffuse reflectance (`deck_albedo`, the two-stream reflectance of an `Aerosol` over a black base) as the ground, or the dark `DEEP_ALBEDO` of 0.03 (stand-in) for a clear giant. The cloud-top algebra above is not coded (no cloud deck is stored on a giant, so `reference_pressure_pa`, decision 23, is given none and holds 1 bar). The Rayleigh depth above 1 bar of 86/14 H2/He at 24.8 m/s^2 is 0.0927 at 550 nm against the 0.089 of part 03 7.1 (4 percent; the same code gives 0.097 for Earth air against 0.099 there). Methane, ammonia and water use the gas table's coarse bins until the Oxford k-tables (GEN.422) and Karkoschka strengths are read, so any of them above 1e-4 in mole fraction gives `colour_quality = "partial"`. **Not built:** the temperature profile and condensation curves that give the cloud base (Guillot, part 06 3.3), the alkali wings and TiO/VO of hot giants, and the Cahoy albedo comparison (a geometric albedo seen from space, a different quantity from the colour at the level).

## 12. The colour model and the thick-air limit (GEN.416 to GEN.418)

**Built (GEN.418, solid worlds; `physics/atmosphere_colour.py`):** a pure function of the air column, the star's temperature and the stored inputs (ozone column, haze, dust, cloud and its cover, the ground's albedo, per wavelength if the caller has it), returning the sky hex, xy and word, the sky's luminance as a share of the star's, the sun disc's colour and luminance, the regime (`thin`, `blend`, `thick`) and, when the aerosols move the colour by 0.03 in xy or more against the same air without them, the label "(aerosol-dominated)". The pieces below are used as follows.

- **One spectrum from two models.** Items 5 and 6 are joined, not switched: the hue is the two-stream (delta-Eddington) diffuse flux at the ground times `(R / R_550)^(1 - w)`, where `R` is the ratio of the single-scatter zenith radiance of item 5 to the two-stream diffuse flux of the same air with no ozone and a black ground, and `w` is a smoothstep of the total depth at 550 nm from 0.2 (w = 0, the single-scatter hue) to 0.4 (w = 1, two-stream only). The ozone and the ground's bounce enter as the two-stream run with them over the run without, so they act the same way in both regimes and the colour has no step (the test scans 400 depths from 0.03 to 3.9 at 550 nm and no channel moves by more than 3 per step of 1.2 percent in pressure). The aerosols are in the single-scatter radiance as extra extinction with their own phase function (Henyey-Greenstein, mixed by scattering depth).
- **Reproduces section 14.1.** Case A #99bfff (xy 0.2452, 0.2522), case C #8bb1ff, case D #ffbb79 (xy 0.4259, 0.3973); case B is #fffbf0 because the 3300 K star takes the BT-Settl offset (the worked example used a blackbody, #fffbec). Case C2 (dust 0.9) is #ffdec2 with the digitized albedos: the direction (blue to warm) is the result, as the note says, the hex is not.
- **The sky's luminance** is the two-stream diffuse flux over the star's flux on a horizontal surface at the top, in both regimes (0.052 for the Earth-like air at 30 degrees from the zenith with a black ground; 0.0015 for the thin CO2 case; 0.3 for 1 MPa of nitrogen). The old item "1 MPa nitrogen is white-cyan with an orange sun" holds (#e7ffff, sun #ffb84a). The old "9 MPa CO2 world is white" of the GEN.418 text is replaced by this section's result (orange, case D), which is the newer.
- **Ozone** takes its Chappuis cross section from `physics/data/ozone_chappuis.json`, the 295 K Brion-team table of the Mainz atlas as 5 nm cell means (peak 5.11e-21 cm^2 at 605 nm; built by `scripts/build_ozone_chappuis.py` from the `ozone_uv` download, so a page that shows a colour needs no download). Its temperature dependence is the ratio of the Burkholder and Talukdar tables at 220, 230, 250, 273 and 298 K to the 298 K table (0.68 at 420 nm for 220 K), interpolated to the temperature of the body's ozone layer (GEN.748: the skin temperature of the body, 2 to the power of -1/4 of its equilibrium temperature, 213.6 K for Earth). The layer sits at the Chapman peak pressure, as in the radiation model. At 300 DU the Earth-like sky loses 2.9 percent of its luminance with the sun 30 degrees from the zenith (xy moves 0.003 and 0.004) and 11 percent with it 78 degrees away (xy moves 0.011 and 0.014), the twilight effect of the section 9 estimate.
- **The cool star** is a blackbody bent by `exp(a t + b t^2)` until its xy is the blackbody's plus the `photometry.cool_offset` plus the seeded spread, so the spectrum stays smooth and the colour does not step at 4500 K.
- **Not in it:** the wavelength-dependent ground (GEN.441 and GEN.462 supply the array), the aerosol optics (GEN.414 supplies `Aerosol` depths, albedos and asymmetry), the giants' reference level and band absorbers (GEN.419), the twilight sun path and Kasten and Young air mass (plane-parallel, sun within 78 degrees of the zenith), and the Monte Carlo correction table of item 5 (the blend keeps the two-stream brightness, so the single-scatter shortfall of 10 to 21 percent does not enter the luminance ratio).

1. **Light.** Planck spectrum of the stored T_eff from 4500 K (xy error at most 0.001 against the CIE tables of part 02 2.9). Cooler: BT-Settl chromaticity (offset from a blackbody 0.002 at 3300 K, 0.02 at 3000 K, 0.05 at 2500 K) with a seeded spread of 0.02 in x and 0.01 in y (STScI file licence not stated).
2. **Rayleigh depth** tau(lambda) = [sum x_i sigma_i(lambda)] P / (m g) (section 6.2; the Bucholtz fit is the regression check for Earth-like air).
3. **Absorbers and aerosols** add extinction and absorption along the path: ozone (9), haze, dust, cloud (10), giant bands (11).
4. **Air mass** by the Chapman integral with X = R / H, never sec(z) beyond 75 degrees; Kasten and Young is an Earth regression target only (10 percent high at the horizon; part 02 2.8).
5. **Thin air, tau(550) up to 0.3** (3.1 bar of Earth-like air): single scatter L = F0 p mu0 [exp(-tau/mu0) - exp(-tau/mu)] / (4 pi (mu0 - mu)). It falls short of the Monte Carlo flux by 10 percent (tau 0.05, mu0 1) to 21 percent (tau 0.1, mu0 0.2); a correction table from a Monte Carlo run is a GEN.418 task (open), and until then the swatch is the single-scatter hue.
6. **Thick air, tau above 0.3: two-stream.** Eddington T = [(2/3 + mu0) + (2/3 - mu0) exp(-tau/mu0)] / (4/3 + tau), R = 1 - T, spherical albedo s = integral of 2 mu (1 - T) dmu (numerical below tau 2), ground albedo A_g: R_tot = R + T A_g (1 - s) / (1 - A_g s), ground flux T / (1 - A_g s); within 1.3 percent of Monte Carlo at mu0 of 0.5 and 1, up to 3.5 percent low at 0.2 (verified, part 02 2.8). **The limit:** T goes as 1/tau and tau as lambda^-4, so the ground light reddens without bound, the direct beam dies as exp(-tau/mu0), and the sky is the star's colour pushed red (14.1, case D). Single scatter saturates to #ff0000 there and must not be used. Absorbing aerosols over a Rayleigh layer use layer adding with part 06 3.10 (open: no delta-Eddington form in the notes); two-stream overestimates reflection for CO2 ice clouds (Kitzmann 2013).
7. **Ground bounce:** divide the downward flux by (1 - A_g s); the ground colour by band is GEN.462.
8. **Sun disc and twilight:** exp(-m tau_total(lambda)), Chapman m, Chappuis included.
9. **Colour:** X, Y, Z with the Wyman, Sloan and Shirley 2013 fit of the CIE 1931 curves (`sky-view.md` section 4), linear sRGB (D65 matrix), negatives clipped, brightest channel 1, sRGB transfer. Absolute (decision 5). The JSON also carries xy and the sky-to-star luminance ratio, since the swatch hides brightness.

## 13. Migrations and the generator epoch

Revisions take the next free numbers at build time (0106 is the highest present; the plan note's "0102 to 0105 stack" is stale), on the planets and moons tables:

- **DB.67 (8.4 phase 2):** adds `bond_albedo`, `greenhouse_factor`, `albedo_basis`; backfill by section 8 for every row with air (airless rows from T_s alone). It runs before the revision that rewrites `atm_molar_density`.
- **DB.71 (8.4 phase 2):** storage for the added gases (layout is DB.71's; any table gas must be storable). Rows saved before GEN.85 with no partials get one seeded draw from the gated mixture tables, flagged. Then `atm_molar_density` is rewritten from the partials and `scale_height_km` from R T_e / (M g) with T_e from the stored albedo.
- **DB.68 (8.4 phase 4):** adds `ozone_column_du`, `tau_haze`, `tau_dust`, `tau_cloud`, `cloud_cover`, `aerosol_draw`, backfilled by the GEN.413 and GEN.414 functions on their own labelled streams, so a rerun repeats and no other draw shifts.

**Epoch statement.** This design changes nothing: epoch 11 stays. The build items do, and each bump is decided there with `scripts/bump_version.py --check-pr`. GEN.410 is a table with no behaviour. GEN.412 only records draws the code already makes (no bump if the check agrees). GEN.411, 422, 423, 424, 452, 454, 455 change values and share one bump in 8.4 phase 2; GEN.413, 414, 456 add stored values and share one in phase 3 (decision 29). GEN.416 to GEN.419 are pure functions of stored values in 9.0 phase 1: no column, no migration, no bump. Galaxies made under an older epoch keep their generated values and receive the new columns by backfill.

## 14. Worked example (computed 2026-10-17, 12:19 PM EDT; dust case 12:45 PM EDT)

### 14.1 Optics and colour

Inputs: test air x = N2 0.78084, O2 0.209476, Ar 0.00934, CO2 0.0004 (normalised), P = 101325 Pa, g = 9.80665 m/s^2; sun 30 degrees from the zenith, view at the zenith (scattering angle 30 degrees), Rayleigh phase 3/4 (1 + cos^2); stars as Planck spectra. `sig(gas, nm)` is the formulas of sources check 2.2 with the reference densities of 2.1 (N2, Ar, CO2 at 2.546899e19, O2 at 2.68678e19 cm^-3) and the listed King factors. Other constants: molar masses (part 01 T1: N2 28.014, O2 31.998, Ar 39.95, CO2 44.009 g/mol), Avogadro 6.02214076e23 (the script; part 04 line 494 quotes 6.022e23, which changes tau by under 0.01 percent), Planck c2 1.438777e-2 m K, the Wyman fit and D65 matrix of `sky-view.md` section 4.

    m = sum(f * M)                                    # g/mol
    tau(nm) = P * NA / (m * 1e-3 * g) * 1e-4 * sum(f * sig(gas, nm))
    sky = planck(T) * p * mu0 * (exp(-t/mu0) - exp(-t/mu)) / (4*pi*(mu0 - mu))     # t = tau(lambda), mu = 1
    edd = planck(T) * ((2/3 + mu0) + (2/3 - mu0)*exp(-t/mu0)) / (4/3 + t)         # thick air, ground albedo 0
    swatch = sRGB(D65 @ XYZ(spectrum)), brightest channel 1                        # absolute, no adaptation

| Case | Inputs | tau 450 / 550 / 650 nm | Colour (hue-only swatch, CIE xy) |
|---|---|---|---|
| A Earth-like air | 5772 K | 0.2207 / 0.0970 / 0.0492 | **#99bfff**, xy 0.2452, 0.2521: a blue sky (the star alone is #fff1ea, xy 0.3265, 0.3360) |
| B Red dwarf world | same air, 3300 K | same | **#fffbec**, xy 0.3243, 0.3441: near-white cream (the star alone is #ffc280) |
| C Thin CO2 | 600 Pa, g 3.7, 5772 K | 0.00591 / 0.00258 / 0.00130 | clear: **#8bb1ff**, a blue hue at 0.030 of case A's brightness |
| C2 with dust | grey dust tau 0.9, g 0.673, omega from Wolff digitized (0.78 at 450 nm, 0.963 at 650 nm for radius 1.8 micron) | | radius 1.8 micron **#ffe7d3**, xy 0.3429, 0.3490; radius 1.2 micron **#ffe6d2** (direction only) |
| D Dense CO2 | 9.2e6 Pa, g 8.87, 5772 K | 37.8 / 16.5 / 8.35 | transmitted light **#ffbb79**, xy 0.4259, 0.3972 |

**A:** tau(550) = 0.0970 against the published 0.09721 (Bucholtz 1995 Table 4, 1962 US Standard atmosphere; 0.2 percent below) and the computed check 0.0969 (0.1 percent above); the Bucholtz fit on the same column gives 0.0968, the difference being the real profile and 330 ppm CO2 behind the published value; the 450 to 650 nm exponent is 4.08, the sky analysis value; tau(450) / tau(650) = 4.5 makes the sky blue. **B:** the same air under a 3300 K star is cream in absolute terms, not blue; the sky analysis found blue (#8cb6ff) only after adapting the eye to the star, which decision 5 does not do. A 3300 K blackbody differs from BT-Settl by 0.002 to 0.003 in xy (part 02 2.9). **C:** the physics of a clear thin CO2 sky is a blue hue at 3 percent of Earth's brightness (tau ratio 0.027); with dust of depth 0.9 the hue turns peach because dust is grey in extinction and brighter in the red (the digitized albedo of Wolff 2009 replaces the earlier assumed 0.74 and 0.89 at 450 nm, which gave #ffe4d0 and #ffece1; the change is within the direction). Case C2 lies beyond single-scatter validity (tau above 0.3), so only the direction (blue to warm) is a result, not the hex. **D:** single scatter returns #ff0000 (the plan base's clipping artefact); the two-stream light at the ground is orange, tau(550) 16.5 against 16.4 in part 03 7.1, and the direct beam is about 0.7 percent of the star's at its brightest wavelength.

### 14.2 The inverse (`planets.py` lines 649 to 714 reproduced)

Inputs (a test body, not a preset): L = 3.828e26 W, d = 1 AU = 149597870.7 km, A = 0.29, M = 0.0290 kg/mol, gravity 1 g, atm_density 1.0 kg/m^3, multiplier 1.7, R = 8.31446261815324, sigma = 5.6703744191844314e-8, CO2 base 0.04345, floor 2.725 K, cap 500.

    S = L / (4 pi d^2) / 1e6 = 1361.17 W/m2     T_e = ((1 - A) S / (4 sigma))^(1/4) = 255.490 K     H = R T_e / (M g) = 7.469 km
    G = min(500, M / 0.04345 * 1.7) = 1.134638  T_s = ((1 - A) S (1 + G) / (4 sigma))^(1/4) = 308.820 K  P = 73250 Pa
    inverse: T_e = H M g / R = 255.490 K, A = 1 - 4 sigma T_e^4 / S = 0.290000000000, G = (T_s / T_e)^4 - 1 = 1.134637514384
    pressure route after a 20 percent freeze-out: T_e = 204.39 K, A = 0.709 (wrong)

G for the three test sets of part 03 5.4 (READ inputs): Earth 0.625, Venus 99.60, Titan 0.694.

## 15. Test targets (physical targets and test rows, never presets)

| Target | Value | State |
|---|---|---|
| tau_R(550 nm), 1 bar Earth-like air | **0.09721** (Bucholtz 1995 Table 4, 1962 US Standard; subarctic summer 0.09688, tropical 0.0976), tolerance 2 percent. Checks: computed 0.0969 (earlier; 0.0980 with Penndorf's F), this note 0.0970, Bucholtz fit on the same column 0.0968 | verified (read); the computed value stays the check. Rule applied: verified (17) |
| tau_R at 350, 443, 550, 670, 865 nm | 0.62, 0.23, 0.10, 0.043, 0.015 (Groot 2020 Table 1); this note 0.629, 0.235, 0.097, 0.0435, 0.0155 | verified (read) |
| sigma_N2(550) 4.635e-27 cm^2; ratios CO2 2.51, H2 0.20 to 0.22, He 0.013, CH4 2.14 (Sneep and Ubachs measured 2.45); exponent about -4.08; doubling P doubles tau | as listed | verified / partly verified |
| Sneep and Ubachs 532.2 nm cross sections (Ar 4.45, N2 5.10, CO 6.19, CO2 12.4, CH4 12.47, N2O 15.90, SF6 32.3 e-27 cm^2); this note's refractivity values 4.56, 5.30, 13.32 for Ar, N2, CO2 | tolerance: the Table 2 gap per gas (section 6.2) | verified (read) |
| Myhre CO2 forcing | 5.35 ln(C/C0); doubling 3.71 W m^-2 | verified (read, recomputed) |
| Mars dust w0 (Wolff 2009) | solar band 0.922 to 0.926 (Table 3); digitized 440 nm 0.765 to 0.772, 500 nm 0.825 to 0.831, 700 nm 0.968 to 0.971 | solar band verified (read); spectrum digitized from Figures 12 to 14 |
| Tabazadeh water pressure at 200 K, 60 weight percent | P0 3.165e-3 mbar, a_w 0.0735 | verified (read, recomputed) |
| Blue sky under 5772 K; cream under 3300 K | #99bfff; #fffbec | verified (computed) |
| Thin CO2 (600 Pa, g 3.7) tau(550) | 0.0026 (part 03 7.1: 0.0027 at 640 Pa, g 3.71) | verified (computed) |
| Star xy at 5772, 6504, 3000 K | 0.3266, 0.3358; 0.3135, 0.3236; 0.4369, 0.4041 (tolerance 0.001) | verified (CIE tables) |
| Kozakis Sun polynomial | 191.8, 199.4, 108.0 DU at 100, 10, 1 percent O2 | verified (computed) |
| Ozone support at the Sun host, 100 percent O2 | contains 279 to 300 DU | partly verified |
| Chappuis depth at 300 DU | 0.042 vertical at 600 nm, 0.8 at air mass 20 | verified (computed) |
| G and the inverse | Earth 0.625, Venus 99.6, Titan 0.694; round trip to 1e-12 | verified (computed) |
| Venus cloud 20 to 30; Earth cloud about 10; Titan haze laws and 9.6; Mars dust exponent near 0 and omega bands | as in section 10 | partly verified |
| Giant albedos (Cahoy; geometric B, V, R, I and Bond) | 0.8 AU 0.556, 0.322, 0.156, 0.047, 0.408; 2 AU 0.733, 0.742, 0.727, 0.673, 0.914; 5 AU 0.609, 0.567, 0.531, 0.453, 0.675; 10 AU 0.450, 0.386, 0.358, 0.313, 0.489 | verified (read) |
| Thick-air limit | tau(550) 16.5 gives orange (xy 0.4259, 0.3972), never #ff0000 | verified (computed) |

Solar System rows (Earth, Venus, Titan, Mars, Jupiter) live in test files only; a test fails when a generator table reads one (TEST.215).

## 16. Open inputs and seeded ranges

| Open input | Range meanwhile | Boss's research |
|---|---|---|
| Mars dust albedo table (w0, n, k) | Wolff 2009 digitized from Figures 12 to 14; 450 nm w0 0.772 to 0.780 | a printed table; Ockert-Bell 1997, Tomasko 1999 Table 8 |
| Titan haze absolute depth | yield LU[0.01, 1] times the supply | Tomasko 2008, Doose 2016 |
| CH4, NH3, H2O visible strengths for giants | methane: free data, partly verified (2026-10-17, about 1:55 PM EDT); ammonia and water: labelled stand-ins; `colour_quality = partial` | Karkoschka 1994, 2009, 2010 (no longer a request to Boss, 2026-10-18, 3:07 PM EDT) |
| Rayleigh of SO2, H2S, H2, He, NH3, H2O | generic 5 to 9 percent (SO2 and H2S polarizability route, plus or minus 8 percent) | Bogaard et al. 1978 (not uploaded) |
| Acid vapour pressure of H2SO4-H2O | Kulmala-Laaksonen form times Zeleznik activity (stand-in) | Ayers 1980 |
| N2, O2, Ar, Kr, Xe in comets | seeded log range from the gated rows | none in Mumma and Charnley |
| Mantle redox rows (AT.a1, a2) | stand-ins | Hirschmann 2012, Gaillard 2021/2022 |
| Bremen ozone cross sections; Huggins fit | coarse bins | fetched at build time (no licence gate, Boss 2:35 PM EDT) |
| Kozakis host T_eff; F-star ozone | bracket of section 9 | Kozakis Table 1 |
| Infrared depth from line absorbers | grey tau0 of part 03 5.4 | correlated-k grids |
| Wind, cloud efficiency, cover, r_eff | u*/u*_t LU[0.3, 3]; cover [0.1, 1]; r_eff 1 to 20 micron | none in the notes |
| Giant albedo fallback | 0.02 to 0.91 | Sudarsky, Cahoy |
| Item numbers | the master chain's GEN numbers are used here; the gap plan and Boss's papers report give the same gaps other numbers (GEN.422 is AD4 Rayleigh there, GEN.471 and 491 acid) | the TODO thread reconciles them |

## 17. Source conflicts and the rule that decided each

Rules (Boss, 12:39 PM EDT): verified wins; two verified sources in doubt, the newest wins; two sources both in question, the most permissive wins. The third rule means the widest value or range either source allows, never a value outside both: Boss confirmed at 1:12 PM EDT that it means the widest range (Decided), as composition-model.md section 4a item 12 records it. Where a paper contradicts a second-hand value in our files, the paper wins and both values are stated.

| # | Conflict | Decision | Rule applied |
|---|---|---|---|
| 1 | Two-stream switch: plan note tau 0.5, part 02 Monte Carlo 0.3 | 0.3 | verified |
| 2 | Earth tau_R(550): computed 0.0969, Bucholtz Table 4 0.09721 | 0.09721 is the test target, 0.0969 the check | verified |
| 3 | Rayleigh error: generic 5 to 9 percent, Sneep and Ubachs gaps 2.5 to 17.8 percent | per-gas gaps; generic only for the six open gases | verified |
| 4 | N2, O2, Ar, CO2: refractivity value (Bates, Bucholtz) against Sneep measured, 2.5 to 7.2 percent lower | refractivity centre, range down to the measured value | verified (centre), most permissive (range) |
| 5 | CH4, N2O, CO2 ratios: He 2021 2.14, 3.42, 2.51; Sneep measured 2.45, 3.12, 2.43; older coded CH4 2.76 | He 2021 centre, Sneep bounds the draw, 2.76 dropped | newest, most permissive (range), verified (2.76) |
| 6 | CO ratio: Allen and Cox 1.24 (partly verified), Sneep 1.21 | 1.21 | verified |
| 7 | H2 ratio: Dalgarno and Williams 0.202, Peck and Huang 0.217 | 0.20 to 0.22 | most permissive |
| 8 | Dust albedo at 450 nm: second-hand 0.74 to 0.89, Wolff digitized 0.772 to 0.780 | digitized | verified |
| 9 | Dust lofting threshold: Shao and Lu against Iversen and White, 15 to 60 percent | draw inside the span | most permissive |
| 10 | Ozone spread: plan note plus or minus 30 percent, published models 279 to 327 DU against the 192 DU fit | structure factor LU[0.85, 1.71] | verified (support), most permissive (width) |
| 11 | Abiotic ozone gate: part 09 at 0.1 percent O2, Kozakis fit valid to 0.01 percent | fit used to 0.01 percent, abiotic branch below | verified |
| 12 | Abiotic O3 "about 10 percent of Earth" (row m2) against Domagal-Goldman 2014 (0.01 to 0.1 percent) | the latter | verified |
| 13 | Water activity of H2SO4-H2O: Tabazadeh against Zeleznik, 10 to 22 percent at 60 to 70 weight percent | Tabazadeh centre, Zeleznik as bracket | newest |
| 14 | Forcing labels: parts file "IPCC 1990" on 5.35 and 0.12, Myhre Table 3 (IPCC 1990 6.3, 0.14) | Myhre 1998 | verified |
| 15 | Comet NH3 bound 3.6 percent, Mumma and Charnley about 1.5 percent | about 1.5 percent | verified |
| 16 | Mantle redox: Gaillard text FMQ-2 +- 0.5, Figure 2 about FMQ-2.75 +- 0.7 | not used; stand-in. If used, the union, FMQ-3.45 to FMQ-1.5 (derived here) | both in question, most permissive |
| 17 | Chain phase of GEN.420, MAP.244, DOC.64: plan note 9.0 phase 2, composition-model phase 3 (GEN.421 phase 4) | composition-model | newest |
| 18 | Giant reference level: composition-model "cloud tops", plan note open | P_ref = min(1 bar, cloud top), which allows both | newest, most permissive |
| 19 | Class albedo default for giants 0.12 to 0.35, Cahoy Bond albedos 0.49 to 0.91 | fallback 0.02 to 0.91 | verified, widest |
| 20 | `greenhouse_multiplier` (code class draw) against `greenhouse_factor` (Boss and part 03) | store `greenhouse_factor` | verified |
| 21 | Migration stack 0102 to 0105 (plan note), 0106 in the repository | 0106 is the highest | verified |
| 22 | Venus cloud depth "30 to 40" (unverified) against 20 to 30 (second-hand, Venus total at 0.63 micron) | 20 to 30 | verified |
| 23 | Sky analysis recommended an adapted colour; Boss's card chose absolute colour | absolute | Boss's ruling, not a source conflict |

## 18. Approved by Boss

**Decided (Boss, 2026-10-17, 9:40 to 9:55 PM EDT, cards answered in the project chat, the rest closed by his typed "I approve all your defaults" at 9:54 PM EDT, relayed by the coordinator): every ruling below is approved as recommended.** By tap: row 1 yes, row 3 derive, row 4 yes, row 5 yes, row 6 yes, row 7 derive, row 8 the switch at tau(550) 0.3, rows 9 to 11 yes, row 12 the research order, row 13 yes, row 15 yes; row 2 (old rows inverted from the stored scale height, class draw on failure, flagged) closed by the typed approval; row 14 needed none. Rulings and the recommendation each card carried:

1. **Greenhouse quantity (13).** Store G = (T_s / T_e)^4 - 1 as `greenhouse_factor`. Recommend yes: it reproduces Earth, Venus and Titan and makes the inverse exact.
2. **Old rows (14).** Invert from the stored scale height, class draw on failure, flagged. Recommend yes; the stored pressure is wrong for any body that froze out gas.
3. **Albedo of new rows.** Derive A from the Rayleigh, aerosol and ground draws, or keep the class-range draw. Recommend derive, so that a cloudy world is not dark; the class range stays the fallback for old rows.
4. **Molar density (16).** Rewrite old rows from their stored partials after the inverse. Recommend yes, in the DB.71 revision.
5. **Ozone spread (19).** Structure factor LU[0.85, 1.71] instead of plus or minus 30 percent, and widen row m1; move the abiotic gate to 0.01 percent. Recommend yes.
6. **Giant reference level (23).** P_ref = min(1 bar, cloud top). Recommend yes.
7. **Giant colour (24).** Derive from optical columns and test against the Cahoy table. Recommend derive; interpolating the published grid is the fallback if band strengths stay unavailable.
8. **Thick-air switch at tau(550) 0.3 (25).** Recommend yes (verified basis).
9. **Dust and haze stand-ins (21, 22).** Accept the flagged ranges, and accept the Wolff spectrum digitized from Figures 12 to 14 as a labelled stand-in until a table is read. Recommend yes.
10. **Cool stars (26).** BT-Settl chromaticity below 4500 K; no licence gate (Boss, 2:35 PM EDT), the licence column keeps `open`. Recommend yes.
11. **Epoch grouping (29).** One bump for 8.4 phase 2, one for phase 3. Recommend yes.
12. **Research order (30).** Recommend: Bogaard 1978 (six open Rayleigh gases), giant band strengths, Tomasko 2008 and Doose 2016, a printed Mars dust table, then the Bremen ozone data (fetched at build time, no licence gate).
13. **Rayleigh centre and range (18).** Refractivity centre for N2, O2, Ar, CO2 with a one-sided draw down to the Sneep and Ubachs value; He 2021 centre for CH4, N2O. Recommend yes; the alternative, centring on the measured values, lowers Earth-like tau by about 4 percent and fails the 2 percent test against the published 0.09721.
14. **Third conflict rule.** The reading in section 17 (the widest value or range either source allows, never outside both) is the one recorded in composition-model.md section 4a item 12, and Boss confirmed it at 1:12 PM EDT (Decided); nothing to approve.
15. **Item numbers.** Ask the TODO thread to reconcile the GEN numbers of the master chain with those of the gap plan before the build items are written. Recommend yes.
