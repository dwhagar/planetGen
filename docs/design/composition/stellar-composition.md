# Stellar composition: the checked design of what every star and remnant is made of

**Status (2026-10-17):** design for GEN.435; Boss approved section 14 as recommended (9:40 to 9:45 PM EDT). Built: GEN.467 (2026-10-17), GEN.439 and GEN.440 (see composition-model.md) and GEN.590 (Boss's "Require the data", 2026-10-18, 12:55 PM EDT); see the "Built" and "As built" paragraphs below. It settles how the composition of every kind of star, brown dwarf and remnant follows from the values the generator already stores, so that GEN.439 (the relations), GEN.440 (the derived function), GEN.464 (the solids law), GEN.428 (the planets' core fraction), GEN.463 and GEN.467 (the two replaced figures) and TEST.211 (the sweep) can be built from it. Master plan: [composition-model.md](composition-model.md). Research: [body-composition-research.md](body-composition-research.md) and [part A](body-composition-parts/A-stars.md); the other sources are named where used.

**Checked against Boss's papers** (`docs/design/research/boss-papers-2026-10-17.md` on origin/main, sections 2.8, 2.11, 3 and 4). Where a paper contradicts a value in our files the paper wins and both values are stated. Changes: the age gate of D10 (3.0 to 8.0 Gy, Haywood 2013), the age-metallicity answer Q3, the test targets, section 12 and the new section 13. **Source rule** (Boss, 12:39 PM EDT on 2026-10-17): when in doubt, verified wins; when two verified sources are in doubt, the newest wins; when two sources are both in question, "most permission take wins". **Assumption for Boss to confirm:** the third rule is read as most permissive, meaning the widest value or range is kept (or both options kept and the spread carried). Each place that resolves a disagreement says "Rule applied", and section 13 lists them all.

**Verification states** (Boss's standing rule, 3:02 AM EDT on 2026-10-17: verified and partly verified values are usable unless a verified source contradicts them): **verified** (a source was read, or the number was recomputed here from read inputs); **partly verified** (ships as sourced, with the limit named); **open** (no source; the range used meanwhile is stated); **stand-in** (a flagged temporary shape with a replacement item). Every number below is cited to a repo file and section, or was computed in this session with the code of the appendix.

## 1. What exists in the code (read on 2026-10-17)

- **Stored per star.** `metallicity_fe_h` on `stars` and `bright_stars` (migration 0095) and on `galaxy_components` (clusters and satellites, which also carry a `seed`), `age_gy`, `lifespan_gy`, `mass_kg`, `yerkes_class`, temperature, luminosity, radius and `uid` (80 bits in the first layout, 96 bits by default since Boss, 2026-10-15; nullable). `bright_stars` also holds `population`, `seed`, `initial_mass_sol` and `phase_end_age_gy`. **The `stars` table has no population and no initial-mass column**; the population is a column only for the bright stars and otherwise lives in the sector's `POPULATION` setting at generation. Neutron stars and black holes have their own tables (mass, radius or horizon, spin) and no [Fe/H].
- **How [Fe/H] is drawn** (`physics/metallicity.py`, `tuning.py` 756-768): normal about a population mean, then clamped to (-3.5, +0.6). Means and sigmas: young 0.00 and 0.10, intermediate 0.00 and 0.15, old -0.20 and 0.25, bulge -0.05 and 0.40, nuclear cluster +0.20 and 0.40, nuclear disc +0.20 and 0.30, halo -1.6 inside 25 kpc and -2.2 outside with 0.5. The disc gradient is 0.06 dex/kpc (Hayden 2015 gives 0.053, part A 2.6; GEN.405 owns it). The stream is `metallicity:{mass}:{age}:{distance}` or `metallicity:bright:{seed}`, separate from every other draw. A binary's secondary copies its primary's value (`generation/binary.py` 136); cluster and satellite stars take their component's peak with sigma 0.15 and 0.30 (`COMPONENT_FE_H_SIGMA`).
- **What reads it:** planet occurrence, the disc's solids (the star's metal mass fraction over the Sun-like reference's, Z/Z0, GEN.464; it was 10^[Fe/H], which is 1.7 to 3 times too low for a metal-poor high-alpha star), giant factor 10^(2 [Fe/H]), the close-binary odds and the horizontal-branch colour shares. There is no X, Y or Z, no element ratio, no alpha enhancement and no age-metallicity link. Age and [Fe/H] are drawn independently.
- **Remnant figures:** main-sequence lifetime 10 Gy times M^-2.5; white dwarfs use the Kalirai IFMR (0.109 M + 0.394, clamped 0.5 to 1.35 Msun) and a radius of 5800 km times M^(-1/3); neutron stars draw mass uniformly in 1.1 to 2.2 Msun and radius uniformly in 10 to 13 km; black holes 5 to 20 Msun with spin 0 to 0.998.

## 2. Decisions

| # | Decision | State | Who and when (EDT) | Range or default used meanwhile |
|---|---|---|---|---|
| D1 | Reference scale is Asplund 2021 proto-solar log eps (X, Y, Z 0.7121, 0.2725, 0.0154); the spread among Asplund 2021, Lodders 2009 and Lodders 2019 is the uncertainty on ratios | Decided | Boss, 2026-10-17, 1:24 AM | Spread on the Sun: Fe/Si 1.69 to 1.94, Mg/Si 0.89 to 1.07, C/O 0.34 to 0.44 by mass (part A 2.3) |
| D2 | Neutron-star radius and main-sequence lifetime are derived, not set (GEN.463, GEN.467) | Decided | Boss, 1:22 AM, "Replace both" | This note gives the inputs (X, Z) and the radius band, section 7 |
| D3 | Bodies draw composition first and keep the inverse for editing; stars keep the stored [Fe/H] and get a composition filled in (Boss, 2026-10-16 11:58 PM: either route works) | Decided for bodies (1:21 AM); Default taken for stars | Boss; this draft | Route 2 for stars would move every stored [Fe/H] and force a rescatter |
| D4 | Composition categories are gated possibility tables, weights add to 100 percent, Dirichlet where only a range is known, Solar System shapes only as stand-ins, one shared component vocabulary, weights are odds | Decided | Boss, 1:40, 1:43, 1:44, 1:47 and 1:51 AM | Loader GEN.472; the star rows of section 3 are rows of it |
| D5 | The verification column carries verified, partly verified, open or stand-in | Decided | Boss, 3:02 AM | Section 12 lists each row's state |
| D6 | Brewer 2016 and Hypatia are downloaded at build time and no derived table is stored; tests run offline on fixtures | Decided | Boss, between 11:19 and 11:23 AM | The fits of part A stay in `tuning.py` as constants with provenance; the download only feeds the distribution test |
| D7 | No new column for the 75 million stars unless a value cannot be derived | Default taken (item text; Core 4 tapped at 12:13 AM) | GEN.435, GEN.440 | Section 4: no value fails this test |
| D8 | White-dwarf debris redox uses the core-formation scale (GEN.498) | Default taken (stands unless Boss overrules) | Relayed 3:36 AM | Not used here beyond the polluted flag |
| D9 | The branch and the four [X/Fe] offsets are **derived from a seeded stream keyed on stored values**, not stored (answers A-S1) | Default taken in this draft | Boss to rule | Section 5, Q1 |
| D10 | No stored population is needed: stars younger than `THICK_START_GY` = 8.0 Gy are thin-disc; older stars use the logistic in [Fe/H]; a bright star's column is used when present | Default taken; the 8.0 is partly verified (Haywood 2013, read in Boss's papers; age scale uncertain by 1 to 1.5 Gy) | Boss to rule | Range 6.5 to 9.5 Gy. The first draft used 3.0 Gy, the project's intermediate/old boundary (`tuning.py` 718); the paper contradicts it, both values stated. Rule applied: verified. Consequences in Q3 |
| D11 | Oxygen from the Nissen 2014 NLTE fit; neon follows oxygen; the switch to magnesium stays | Default taken; Open on the physics | Boss to rule | Z/Z0 spread between the two options 8 percent at most, section 5 Q2. Rule applied: most permissive (both kept, spread carried) |
| D12 | Ages and metallicities stay independent; only the age gate of D10 links them | Default taken; stand-in (Haywood's old thick sequence implies a real link) | Boss to rule | Section 5 Q3 |
| D13 | Remnants: keep Kalirai IFMR, replace the white-dwarf radius law by the integrated Chandrasekhar table, neutron-star mass from the observed mixture and radius from the band, black holes carry no composition | Default taken | Boss to rule | Section 7 |
| D14 | The function returns the birth composition; a photospheric correction is not applied | Default taken; Open | Boss to rule | Solar value 0.064 dex (metals) and 0.070 dex (He) over 4.57 Gy (part A 2.1); no age law read |
| D15 | No table of GEN.407 changes in this design; the tip by metallicity is GEN.502 and reads the derived Z | Default taken | Boss to rule | Section 6 |
| D16 | Members of a multi-star system, and the stars of a cluster, share the group's [Fe/H], branch and offsets | Default taken | Boss to rule | Section 9; within-group scatter 0 (open, none sourced) |
| D17 | The disk solids law reads Z/Z0 (GEN.464) | Default taken (a filed finding, not a Boss ruling) | TODO thread, from the lane draft of about 1:10 AM | Z/Z0 values of section 3 |

## 3. The composition of a star from the stored values

Inputs, all stored or fixed constants: [Fe/H] (call it F), age, kind, mass, the optional population, and a key. Model constants are in the appendix and trace to part A 2.2, 2.4 and 2.7.

1. **Key.** `abundance:{uid}` (the star's `uid`; where it is null, `abundance:{mass}:{age}:{distance}`, the quantities that already seed the [Fe/H] stream; a bright star uses `abundance:bright:{seed}`). A stream of its own, so no earlier draw moves. Draws in fixed order: u (branch), four normals (C, O, Mg, Si), u (CEMP), u (CEMP carbon). All are always consumed.
2. **Branch.** thin if age < `THICK_START_GY` (8.0 Gy, D10); otherwise high-alpha if u < P(F) = 1/(1+exp((F+0.350)/0.132)) (Bensby 2014 fit, part A 2.6), else thin. Because u is fixed per star, the branch flips at most once as F falls, so an edited [Fe/H] never flickers. P(F) is 0.50 at -0.35, 0.757 at -0.5, 0.968 at -0.8 and 0.9999 at -1.6, so the halo is high-alpha without any population label. The fit comes from part A; Bensby 2014 is not among Boss's uploads, so it could not be checked against the paper (partly verified).
3. **Means** [X/Fe] for C, O, Mg, Si, Ca from the hinge fits of part A 2.7 (N14 NLTE for C and O), with F clamped to each branch's sampled range (thin -0.9 to +0.5, high-alpha -1.6 to +0.1).
4. **Residuals:** the four normals through the Cholesky factor of the BF16 correlation matrix (C-O 0.56, C-Mg 0.69, C-Si 0.71, O-Mg 0.36, O-Si 0.37, Mg-Si 0.77), times the branch sigma (thin 0.044, 0.052, 0.066, 0.042; high 0.050, 0.065, 0.073, 0.073 for C, O, Mg, Si).
5. **Others:** S, Ar follow Ca, Ne follows O, everything else [X/Fe] = 0 (N not fitted).
6. **CEMP:** for F < -2, with probability 0.25, [C/Fe] is replaced by 0.7 + 1.3 u (Frebel and Norris: 20 to 30 percent, [C/Fe] > 0.7). Part A's rise to probability 1 below -4 is unreachable because the stored [Fe/H] stops at -3.5.
7. **X, Y, Z and elements.** r = Z/X = (sum over metals of PROTO_i 10^(F + [X_i/Fe]) + remainder 10^F) / X0; X = (1 - Y_P)/(1 + (1+g) r); Z = r X; Y = 1 - X - Z; each element w_i = X (PROTO_i/X0) 10^(F + [X_i/Fe]). Y_P = 0.24672 (Planck 2018, READ). **g = 1.674**, derived as (1 - Y_P - X0)/Z0 - 1, not the rounded 1.68 of part A 2.4: with 1.68 the Sun reproduces X to 1e-4 only (0.71201 against 0.7121), with the derived g exactly. The observed dY/dZ is 2.1 +- 0.9 (Casagrande 2007, READ), so 1.674 is 0.5 sigma low and inside it. The 14 listed elements cover 99.7 percent of Z0; the remaining 4.9e-5 scales with Fe (partly verified, flagged).
8. **Ratios:** log eps_i = A21 + F + [X_i/Fe]; Fe/Si, Mg/Si, C/O, O/Si as number and mass ratios (mass factors 1.9884, 0.8654, 0.7507 and 0.5697).
9. **Phase shifts** (additive on the birth values, applied by class): subgiant IV none; red giant III and II first dredge-up (section 10, example F); supergiants IB to 0 open (the grid stops at 2.5 Msun); Wolf-Rayet and carbon stars from the flags of GEN.243 and GEN.249; Pop III only if a kind is added (X = 0.75328, Y = 0.24672, Z = 0, no draw).
10. **Radioactive birth abundances** for GEN.459: Eu from the Eu fit, Th = Eu x 0.405 x 10^[Th/Eu] with [Th/Eu] = 0.014 +- 0.045 (Botelho 2019, READ), U/Th birth 0.572; K not fitted. Limit: Th itself is measured near solar [Fe/H] only (slope -0.338 +- 0.049 inside +-0.1 dex, master 15.1), so Th elsewhere is inferred through Eu (partly verified). Cross-check on U/Th with Palme and O'Neill 2003 (first edition, read in Boss's papers: Th 83.4 ppb, U 21.8 ppb): decayed back 4.568 Gy (appendix `birth_u_th`) the mass ratio is 0.554; the Sun's atomic 0.572 is 0.587 by mass, so the two sit 6 percent apart. Our files' "Palme and O'Neill 2014" values (Th 0.0849, U 0.0229 ppm; 1.8 and 5.0 percent above) cannot be confirmed and give 0.572, so the label here is "2003 edition" and the 6 percent is carried as the uncertainty on U/Th. Rule applied: verified for the label; most permissive for the spread.

**Z/Z0 against 10^[Fe/H]** (recomputed here; matches part A 2.5 to three digits in Z/Z0): thin +0.2 gives 0.862 (neon on O) or 0.884 (neon on Mg); thin 0.0 gives 0.981 or 0.988; high-alpha 0.0 gives 1.298 or 1.274; thin -0.5 gives 1.380 or 1.386; high-alpha -0.5 gives 2.043 or 1.962; high-alpha -1.0 gives 2.982 or 2.760; -2.0 gives 3.027 or 2.798; -3.0 gives 3.031 or 2.802. GEN.464's test reads this table.

## 4. What is derived on demand and what is stored

| Quantity | Stored or derived | Note |
|---|---|---|
| [Fe/H], age, lifespan, mass, class, temperature, luminosity, radius, `uid` | stored (unchanged) | The only inputs; no new column |
| Population | stored on `bright_stars` only; otherwise replaced by the age gate (D10) | Not needed on the 75 million |
| Branch, [Mg/Fe], [Si/Fe], [O/Fe], [C/Fe], CEMP flag | derived from the key and stored values | Reproducible at any worker count; the same key gives the same draw |
| X, Y, Z, element mass fractions, Fe/Si, Mg/Si, C/O, O/Si, Z/Z0 | derived | A cache is allowed; a test asserts cache equals recomputation |
| Surface shifts after dredge-up, carbon-star C/O, Wolf-Rayet fractions | derived from class, mass, Z and flags | Flags belong to GEN.249 and GEN.243 |
| Central H fraction, fuel factor | derived | Section 6 |
| Th, U, K birth abundances | derived | Decay to the planet's age is GEN.459's |
| White-dwarf core type and radius, neutron-star layers, black-hole horizon and spin quantities, brown-dwarf radius | derived from stored mass (and age for the brown dwarf) | Section 7 |
| An admin's wish to set one star's [Mg/Fe] by hand | the only case that cannot be derived | Not offered; if Boss wants it, a sparse override record, never a column on 75 million rows (open) |

## 5. The open questions of the item, answered

**Q1. Are the branch and the four offsets stored or derived from the seed?** Derive. Part A 3.3 and part B list five new columns; GEN.440 and the item rule them out, and nothing needs them stored: every value is a function of the stored [Fe/H], age and key. The branch draw u is fixed per star so edits are stable; changing a star's [Fe/H], age or mass in the admin (edit-everything.md) re-derives the composition with no stored state to repair. Cost: one 4 by 4 Cholesky factor (a constant) and seven draws per call.

**Q2. Neon and oxygen.** Take oxygen from the Nissen 2014 NLTE triplet fit (plateau +0.606 on the high-alpha branch; the Bensby 2014 own-O plateau is +0.582, the two differ by up to 0.1 dex elsewhere, part A 2.7 note b) and let neon follow oxygen, with a `NE_FOLLOWS` switch to magnesium. Neon is not measured in photospheres (set values range 7.93 to 8.15, part A 2.1), so it is the least certain entry of Z. The two options differ in Z/Z0 by at most 8 percent (2.982 against 2.760 at -1.0); the generic single-[alpha/Fe] model gives about 2.0, so both stand 1.4 to 1.5 times above it. The physical reason for neon with oxygen is recalled, not read: open. Part A 3.2 step 4 says neon follows magnesium while the filing of GEN.464 says the lane preferred oxygen; this draft follows GEN.464 and keeps both numbers as the filing requires.

**Q3. Age-metallicity link (A-D3).** [Fe/H] stays unlinked to age (D12). The age gate carries the one link, and Boss's papers move it. Haywood 2013 (read): the thick disc formed from about 13 to about 8 Gy ago, with the age scale uncertain by 1 to 1.5 Gy; the inner thin disc starts about 8 Gy ago and the outer thin disc 9 to 10 Gy ago; thick/thin mass ratio at least 0.2, possibly up to about 1. The first draft gated at 3.0 Gy, the project's intermediate/old boundary, which would make 3 to 8 Gy stars thick-disc. The paper contradicts that, so the paper wins and the gate is 8.0 Gy (range 6.5 to 9.5). Rule applied: verified.

The cost, computed here (appendix P(F), 3e5 draws, ages uniform): for the project's old disc (3 to 10 Gy, [Fe/H] N(-0.20, 0.25)) the high-alpha share is 0.329 with the gate at 3.0, 0.165 at 6.5, 0.095 at 8.0 and 0.023 at 9.5. For the solar-neighbourhood mix N(-0.10, 0.22) with ages uniform 0 to 10 Gy it is 0.216 with no gate, 0.151 at 3.0 and 0.043 at 8.0. P(F) is Bensby's fraction over all ages, so applying it only above 8 Gy under-produces thick-disc stars: 0.095 of old-disc stars, against a thick disc of at least 0.167 of disc mass (ratio 0.2; a number share is not a mass share, and the whole-disc ratio belongs to the density lane, GEN.134). A stand-in correction multiplies P by (10 - 3)/(10 - gate) = 3.5 and caps it at 1; it gives 0.189 for the old disc, but it jumps at 10 Gy and gives 0.23 at [Fe/H] = 0, so it is not adopted. The honest fix needs the age-resolved thick fraction from Bensby 2014 (not among Boss's uploads, so every Bensby-based value here is part A's, unchecked against the paper) or Pinna 2019b (not uploaded; Pinna 2019a covers one galaxy, FCC 170, with about 3 Gy age resolution and no thick/thin ratio, so it cannot test the gate). Open. At [Fe/H] = -0.5, where the branches differ most, Z/Z0 is 1.380 thin against 2.043 high-alpha (section 3), so a thin-branch star of 3 to 8 Gy there has 32 percent less Z than a thick-disc star would.

Also open, for GEN.134: the project's old window ends at 10 Gy and the halo starts at 10 Gy (`tuning.py` 718), while Haywood has the thick disc forming back to 13 Gy and Aumer and Binney fit a disc age of 12.6 Gyr (tau_max, Table 6 fit), so the old disc window stops 2.6 to 3 Gy short of the disc's start. Not changed here.

Independence of [Fe/H] and age otherwise loses little. Bensby 2014 gives -0.027 dex/Gy with scatter 0.265 dex (second-hand through part A, partly verified), which for ages uniform over 0 to 10 Gy is an explained sd of 0.078 dex and a correlation of 0.28. A young star (sigma 0.10) is below -0.5 with probability 3e-7 and an intermediate one (sigma 0.15) 4e-4. Linking would move every stored [Fe/H]. Bensby's kinematically thick fraction of 0.08 at 0 to 2 Gy (shape only) is not reproduced; Haywood's age floor says it is a kinematic misclassification, so a young alpha-rich star is not produced on purpose.

**Q4. Remnant limits.** Section 7.

## 6. How the abundances change the star, and whether the tables of GEN.407 change

GEN.407 is built (PR #1985, epoch 4); GEN.502 holds what it left. What composition does, with what the repo research can support:

- **Lifetime.** Fuel is the hydrogen fraction X. Against the proto-solar 0.7121: example A 1.014, B 1.036, C 1.053, Pop III 1.058 (section 10). Burning rate depends on luminosity, which depends on opacity and mean molecular weight; no read source in the repo gives dL/dZ or dT/dZ. GEN.467 should take (mass, X, Z); the project's 10 Gy M^-2.5 is replaced there. **Built (GEN.467, 2026-10-17):** `physics/ms_lifetime.py` is the fit of Hurley, Pols and Tout 2000 (MNRAS 315, 543, eqs. 4 to 7 and the Appendix coefficients, read from the paper), a function of mass and [Fe/H] (`Z = 0.0154 x 10^[Fe/H]`, the proto-solar Z above); a solar mass lasts 10.2 Gy, 20 solar masses 9.0 Myr (5.6 before), 150 solar masses 3.2 Myr (0.04 before). Every caller passes the default [Fe/H] of 0 until the composition is drawn first (9.0). For a check, the Sun's central hydrogen from three models is 0.708 at zero age, 0.340 at 4.57 Gy, 0.145 at 6.5 and 0.008 at 8 (Bahcall 2001, READ); the fit X_c = X0 (1 - (t/8)^1.15) returns 0.336, 0.150 and 0 (error under 0.01).
- **Opacity.** Z sets the Rosseland opacity. The only opacity law in the repo research is the hydrogen-burning edge M_edge = 0.0865 (kappa/1e-2)^(-1/9) Msun (Burrows 2001, READ): a 10 percent opacity change moves the edge 1 percent. Nothing else in the generator reads it.
- **Luminosity, radius, temperature, colour on the main sequence.** Boss's cluster text says metal-poor stars are hotter and brighter at the same mass (globular-clusters.md 34, an unchecked note) and puts the 12 Gy turn-off at 0.85 Msun against 0.93 in the code; no tracks were read (MIST and BaSTI were unreachable, stellar tables check). Open; the tables stay Z-independent as a stand-in. Colour is a blackbody of the stored temperature (GEN.416); composition does not change it in this design.

| GEN.407 value (`tuning.py`) | Change? | Reason |
|---|---|---|
| `GIANT_BRANCH_TIP_LUMINOSITY_SOL` 2500 | Yes, by GEN.502, from Z | Serenelli 2017 Table 2 (READ): log L 3.2335, 3.326, 3.386 at Z 1e-4, 1e-3, 4e-3 (midpoints of the printed ranges). Interpolated in log Z and held outside: example C gives 2183 Lsun, D 1712, B 2432; above Z 4e-3 the table ends, so solar-Z stars stay at 2432 (flagged) |
| `GIANT_BRANCH_TIME_SLOPE` 0.78 | No | Core-mass law is GEN.502's, not composition |
| `HORIZONTAL_BRANCH_COLOURS` rows, keyed on [Fe/H] | No | Calibrated on 57 clusters in [Fe/H] (Zinn and West scale); a high-alpha star has Z/Z0 up to 3 at the same [Fe/H], but no source tells whether colour follows Z or [Fe/H] for field stars. Open; the scale of the stored [Fe/H] against Zinn and West is unstated in the code |
| HB age shift 0.1 dex/Gy, HB luminosity 35 to 50 | No | Placeholder and Z-independent in the sources read |
| Main-sequence L(M), R(M), `SUN_EFFECTIVE_TEMPERATURE_K` | No | No read source |
| WD IFMR, cooling law | No | IFMR metallicity dependence not read (Cummings 2018 PARSEC and MIST differ by 0.03 to 0.07 Msun) |

## 7. Brown dwarfs and remnants

| Kind | Composition derived from the stored values | State |
|---|---|---|
| Brown dwarf | Natal gas X, Y, Z as a star; fully convective. Deuterium limit M_D = 13.55 - 0.9 log10(Z/0.0189) - 0.2 (Y-0.25)/0.01 MJ (linear fit of Spiegel 2011 Table 1; their solar Z is 0.0189, not 0.0154; spread 11.0 to 16.3 MJ); lithium burned above 0.06 Msun (63 MJ); hydrogen burning above 0.07 to 0.074 Msun at solar Z, 0.092 at Z = 0 (Burrows 2001); radius and Teff from the BHAC15 grid (solar Z only), Li depletion from the same grid | D, Li, H limits verified at solar Z; Z dependence of H limit open (log-linear is a parametrisation); radius at other Z open |
| White dwarf | Mass stored. Core He below about 0.45 (recalled, open), C/O above, O/Ne above 1.05 Msun (Althaus 2010, READ). Atmosphere DA with p 0.80 (Xu 2024, READ), else DB, DC, DQ, DZ; 25 to 50 percent of cool ones (under 10,000 K, the floor of the hot white dwarfs) polluted: `detail["polluted"]`, drawn from the star's own stream at the middle share 0.375 and never stored (GEN.253). Layer masses and the C:O split: open. Radius from the integrated table of the appendix (Sirius B 1.018 Msun gives 0.00803 Rsun against 0.00803 +- 0.00011 measured) | IFMR, DA share, pollution verified; core C:O and layers open |
| Neutron star | Mass stored; draw from the observed mixture (slow 1.49 +- 0.19, recycled 1.54 +- 0.23, floor 1.17, heaviest 2.35 +- 0.17; M_TOV 2.01 to 2.16). Radius band 11.6 to 13.1 km, nearly mass independent (Miller 2021, hull of the 1-sigma intervals; J0740 alone is 13.7 +2.6/-1.5, so the band is not a cap); the function is GEN.463's. Layers by density: 56Fe to 62Ni (8.02e6 to 2.71e8 g/cm3) through 118Kr (4.27e11), neutron drip near 4e11, inner crust Zr and Sn clusters, uniform n, p, e core (Chamel and Haensel 2008, READ). Crust thickness 1 km and 0.01 Msun are recalled, open | Mass, band, crust sequence verified; radial structure needs an EOS: open |
| Black hole | No composition (no-hair theorem, READ). Mass, spin, charge 0 stored or fixed; r_s, Hawking T, horizon, ISCO and efficiency by formula (recomputed: 10 Msun gives r_s 29.53 km, T 6.17e-9 K, and at spin 0.9 ISCO 2.321 and efficiency 0.156) | Verified; spin cap 0.998 recalled |
| Protostar, pre-main sequence | Natal gas of the cloud, deuterium not burned until the limit above | Partly verified |
| Stripped-helium, Thorne-Zytkow, Wolf-Rayet | Wolf-Rayet: WN X_N 1 percent and H-free early, late WN X_H 15 +- 10 percent, WC carbon 20 to 55 percent by mass (Crowther 2007, READ); WO, which the paper names without fractions, carbon 30 to 55 and oxygen 10 to 40 percent by mass, a stand-in (registry row BC24, GEN.831). Stripped-helium and Thorne-Zytkow: no source read | WR partly verified; others open |

**Findings on the existing figures** (computed here): the code's white-dwarf radius 5800 km times M^(-1/3) is 22 percent below the integrated radius at 0.6 Msun (6877 against 8842 km at 0.6 Msun), equal within 1.5 percent at 1.0 Msun and 31 and 64 percent above it at 1.2 and 1.3 Msun; `mathcheck.py` tests the 0.6 Msun radius against Earth's 6371 km within 15 percent, which hides this. The uniform neutron-star mass draw of 1.1 to 2.2 Msun puts as much weight at 2.0 as at 1.4, while the observed mixture peaks near 1.3 to 1.5.

## 8. How a star's iron, magnesium and silicon set the rocky planets' core fraction (GEN.428)

The planet draw reads the host's number ratios from section 3. With mu the standard atomic weights (Fe 55.845, Mg 24.305, Si 28.0855, Ca 40.078, Al 26.9815, O 15.999), the core mass fraction with all iron in the core is, from Unterborn 2023 Eq. 8 as read in part B 2.4:

CMF = (Fe/Mg) mu_Fe / [(Fe/Mg) mu_Fe + mu_bar], mu_bar = (mu_Mg + mu_O) + (Ca/Mg)(mu_Ca + mu_O) + (Al/Mg)(mu_Al + 1.5 mu_O) + (Si/Mg)(mu_Si + 2 mu_O).

Check: Fe/Mg 0.71, Si/Mg 0.79, Ca/Mg 0.07, Al/Mg 0.09 give 0.2917 (part B table: 0.292). Ca and Al come from [Ca/Fe] (section 3) and the Bensby [Al/Fe] fit (thin 0.032 - 0.029 min(F - 0.06, 0); high 0.225 - 0.298 max(F + 0.48, 0)). Results:

| Case | CMF |
|---|---|
| Asplund 2021, Lodders 2009, Asplund 2009 reference ratios | 0.307, 0.302, 0.315 (set spread 0.013) |
| Lodders 2003 (GEN.428 text) | 0.320 |
| Draw, solar-neighbourhood mix N(-0.10, 0.22), 6e4 stars: median and 2.5 to 97.5 percent | 0.284 [0.199, 0.337] |
| Thin disc only N(0, 0.2) | 0.293 [0.242, 0.341] |
| Halo N(-1.6, 0.5) | 0.184 [0.144, 0.237] |
| Observed BF16 Fe/Si 1.17, 1.81, 2.39 at Mg/Si 0.93 | 0.228, 0.314, 0.376 |
| Thiabaud box corners (Fe/Si 1.34 and Mg/Si 1.85; 2.81 and 0.66) | 0.189, 0.450 |

Metal-poor high-alpha stars give smaller cores (Fe/Si falls to about 1.0 while Mg/Si rises to 1.2). GEN.428 quotes a stellar range of 0.2 to 0.47; that is the Thiabaud compilation envelope. The drawn distributions reach 0.34 at 97.5 percent in the solar neighbourhood, so the upper part of 0.2 to 0.47 is a plausibility box, not a drawn range (master 7). Mantle FeO and core light elements move the planet's CMF by about 0.04 per 0.1 (part B 2.3) and belong to GEN.428, not to the star.

For comparison only: Palme and O'Neill 2003 (first edition, read in Boss's papers, Table 3) give a primitive mantle of MgO 36.77, Al2O3 4.49, SiO2 45.40, CaO 3.65, FeOt 8.10 percent, total 98.41 (not 100; the rest is minor oxides 1.34 and trace), so Mg/Si is 1.045 by mass (recomputed here: 22.17 over 21.22). The Sun's by-mass Mg/Si in this note is 0.949, and a planet's mantle after core formation is not its star's bulk ratio, so this is not a contradiction and no target is taken from it. The label in our files, "Palme and O'Neill 2014" (Th 0.0849, U 0.0229 ppm), cannot be confirmed from the papers and is not used here.

## 9. Clusters and multi-star systems (GEN.269)

- **Multi-star systems.** All members form from one cloud: they share the group's [Fe/H] (already true for a binary, `binary.py` 136; for the member tree this is the storage half of GEN.269, Core 2's finding of 2026-10-16, 4:58 AM EDT), and the branch and four offsets come from the key of the group's primary (`abundance:{primary uid}`), not the member's. Birth X, Y, Z and ratios are identical across members. Each member then takes its own phase shift (a giant sibling differs from a dwarf one) and its own remnant layers; the pollution and carbon-star flags are per member. A merger product keeps the real system age for the age gate. Mass-transfer pollution of an accretor is open. Within-group scatter: 0, open; no source was read (Hinkel 2014's 0.11 dex is between catalogues, not between siblings).
- **Clusters.** The cluster record holds its [Fe/H] (peaks -1.55 and -0.55, `CLUSTER_FE_H_GROUPS`); its stars draw [Fe/H] about it with sigma 0.15 (recalled in the code, open). Branch and residuals are keyed on the cluster seed (`abundance:cluster:{seed}`), so the residuals are shared and only the F-dependent means move with each star. A cluster of 10 to 13 Gy ages passes the age gate and the logistic: at -1.55 it is high-alpha with probability 0.9999, at -0.55 with 0.82. Satellite galaxies follow the same rule as a stand-in; their real alpha pattern is not in the research (open).

## 10. Worked examples

Inputs are illustrative, not presets and not copied from the Solar System. Keys carry a suffix chosen as the first that keeps all four normals inside +-1.6 for readability. Code in the appendix; X, Y, Z sum to 1 in every case.

| Kind and inputs | Steps | Outputs |
|---|---|---|
| **A, main sequence, thin.** F -0.20, age 2.0 Gy, key `example-A1` | age 2.0 < 8.0: thin. Normals -0.15, 1.30, 0.005, 0.42, 0.0; residuals C -0.014, O 0.110, Mg -0.004, Si 0.026, Ca 0.012; [X/Fe] C -0.012, O 0.194, Mg 0.003, Si 0.017, Ca 0.072, Ne 0.194 (refitted to Brewer and Fischer 2016, GEN.988 to GEN.991; the first fit gave C -0.004, O 0.115, Mg 0.076, Si 0.043, Ca 0.06) | X 0.7190, Y 0.2682, Z 0.01282 (Z/(Z0 10^F) 1.319). Mass fractions: O 6.38e-3, C 1.76e-3, Ne 1.89e-3, Fe 8.4e-4, Mg 4.5e-4, Si 4.9e-4, N 4.9e-4. Fe/Si 1.703, Mg/Si 0.917, C/O 0.276 (mass). CMF and fuel factor not recomputed |
| **B, old, high-alpha.** F -0.80, age 9.0 Gy, key `example-B2` | age 9.0 >= 8.0 passes the gate; P(F) 0.968, u 0.533: high. O +0.503, Mg +0.312, Si +0.175, C +0.131 | X 0.7375, Y 0.2566, Z 0.00591 (Z/(Z0 10^F) 2.42). Fe/Si 1.186, Mg/Si 1.303, C/O 0.188. CMF 0.205. Tip 2432 Lsun |
| **C, halo.** F -1.60, age 12 Gy, key `example-C1` | P 0.9999, high. O +0.681, Mg +0.424, Si +0.240, C +0.216 | X 0.7497, Y 0.2490, Z 0.00135 (3.49). Fe/Si 1.019, Mg/Si 1.448, C/O 0.152. CMF 0.175. Tip 2183 |
| **D, carbon-enhanced metal-poor.** F -3.0, age 12.5, key `example-D8` | F < -2: CEMP draw 0.069 < 0.25, [C/Fe] = 0.7 + 1.3 x 0.396 = 1.215 (a key suffix that gives a CEMP star; most keys do not) | X 0.7530, Y 0.2469, Z 9e-5. C/O 1.80 mass (2.40 number), Fe/Si 0.99. A rocky planet here could be carbon-rich |
| **E, Pop III** (only if a kind is added) | Z = 0 exact, no draw | X 0.75328, Y 0.24672, Z 0 |
| **F, red giant.** Mass 1.5 Msun, F -0.10, age 4.1 Gy, key `example-E4` (code lifetimes: giant between 3.99 and 4.35 Gy) | Branch thin (age 4.1 < 8.0; the logistic alone gives P 0.131 against u 0.393, thin either way). Birth X 0.7197, Y 0.2678, Z 0.01257, C/Z 0.177, N/Z 0.050. Dredge-up (Boothroyd 1999, READ): weight in log Z between Z 0.001 and 0.02 is 0.845; dY = +0.0165; C factor 0.7257; N factor 2.02 (Z 0.02 only; low-Z N not extracted) | Surface X 0.7032, Y 0.2842, Z 0.01257; C/Z 0.129, N/Z 0.100; C/N number 4.17 to 1.50; C/O 0.409 (number). Tip 2432 Lsun |
| **G, AGB.** Mass 2.0 Msun, Z 0.0125 | Carbon-star minimum mass 1.4 at Z 0.01 and 1.5 at Z 0.02 (Karakas 2016, 2ND-HAND), linear in Z gives 1.425 | 2.0 > 1.425: C/O above 1 is possible once thermal pulses have run. The flag is GEN.249's; the C/O value is open (yield tables not extracted) |
| **H, Wolf-Rayet WC.** key `example-W1` | C = 0.20 + 0.35 u, u 0.767 (the refined 20 to 55 percent) | X_H 0, C 0.469, the rest He and O, split open |
| **H2, Wolf-Rayet WO.** Subtype from the spectral text ("WO2"), a drawn position u of 0.5 | C = 0.30 + 0.25 u, O = 0.10 + 0.30 u, both at one u (BC24 stand-in); a WN star is early (X_H 0) up to subclass 5 and late (X_H 0.05 to 0.25) from 6 | At u 0.5: C 0.425, O 0.25, He 0.325; the three always add to 1 |
| **I, brown dwarf.** 0.050 Msun (52.4 MJ), age 1 Gy, gas of A | M_D = 13.55 + 0.189 - 0.326 = 13.41 MJ; 52.4 > 13.4 so deuterium burned; below 62.9 MJ so lithium kept; below 0.07 Msun so no hydrogen burning | R 0.0960 Rsun (66,787 km), Teff 1383 K (BHAC15, solar Z), Li depletion 0 to -0.03 dex |
| **J, white dwarf.** Mass 0.61 Msun (progenitor 1.98 Msun by the inverse IFMR), key `example-WD1` | Core C/O (0.45 to 1.05). Atmosphere u 0.391 < 0.80: DA | R 8752 km (0.01258 Rsun) against 6839 km from the code's law |
| **K, neutron star.** 1.36 Msun, R 12.45 km (Miller 2021 R(1.4)) | Mean density 3M/(4 pi R^3) | 3.35e14 g/cm3 = 1.22 times saturation (2.7 to 2.8e14); compactness 2GM/Rc^2 0.323; g 1.16e14 cm/s2. The layer sequence of section 7 applies from the surface inward; thicknesses need an EOS (open) |
| **L, black hole.** 12 Msun, spin 0.7 | Formulas of the appendix | r_s 35.44 km, T_H 5.14e-9 K, r_+ 30.37 km, ISCO 3.393 GM/c^2 = 60.12 km, efficiency 0.1036 |

## 11. Test targets (the Sun's composition is a target, never a preset)

| Target | Expected | Tolerance | State |
|---|---|---|---|
| Scale reproduces itself: F = 0, all [X/Fe] = 0 gives (X, Y, Z) = (0.7121, 0.2725, 0.0154) | exact | 1e-6 with g derived from the reference; 1e-4 with 1.68 | verified |
| Thin branch means only at F = 0 (the Sun's role as a target) | X 0.7128, Y 0.2720, Z 0.01512 | Z within 3 percent of 0.0154 | verified |
| Z/Z0 against 10^F at F 0, -0.5, -1, -2 for thin and high-alpha | table of section 3 | 0.01 | verified (GEN.464) |
| X + Y + Z = 1; Y >= Y_P; Z monotone in F on a branch for fixed residuals; branch flips at most once | invariants | 1e-9 | verified (checked in the session) |
| 6e3 draws, mix N(-0.10, 0.22): Fe/Si 1.82, Mg/Si 0.98, C/O 0.38, high share 0.224 (age gate off, as part A 2.8 drew it; part A had Fe/Si 1.62, Mg/Si 1.02 and C/O 0.40 before the refit to Brewer and Fischer, GEN.988 to GEN.991) | the Brewer and Fischer medians below | as stated | verified |
| Age gate: no star younger than `THICK_START_GY` is high-alpha; old disc (ages uniform 3 to 10 Gy, [Fe/H] N(-0.20, 0.25), 3e5 draws) high-alpha share 0.095 at 8.0 Gy (0.329 at 3.0) | invariant; 0.095 | exact; 0.01 | stand-in (the share depends on the open age-resolved fraction, Q3) |
| Against BF16 medians (Fe/Si 1.81, Mg/Si 0.93, C/O 0.385) | within 0.2, 0.1, 0.05 | as stated | verified (the draw sat 0.19 below on Fe/Si before the thin Mg and Si refit, GEN.989; the table's C/O median over all its stars is 0.368) |
| Mg-Si residual correlation 0.68 (0.77 of the 844-dwarf table), thin Mg/Si ratio residual 0.054 dex | within 0.03 and 20 percent | | verified |
| Sun's central hydrogen at 4.57 Gy, 1 Msun | 0.340 | 0.03 | verified |
| First dredge-up, 1 Msun, Z 0.02: Y +0.0225, C/Z x0.88, N/Z x1.45 | | 10 percent | verified |
| White-dwarf radius at 1.018 Msun | 0.00803 Rsun | 5 percent | verified |
| DA peak mass | 0.59 Msun | 0.05 | verified |
| Heaviest neutron star reachable by the draw tail | 2.35 +- 0.17 | | verified |
| Black hole r_s at 10 Msun | 29.53 km | 1e-3 | verified |
| Brown dwarf radius, 0.07 Msun at 1 Gy | 0.0944 Rsun | 5 percent | verified |
| Same star, same composition at one and four workers; no table row lacks a state; the loader counts each state | | | TEST.211, TEST.214 |
| Solar values live in test files only; no generator table reads them | | | TEST.215 |

## 12. Contradictions and gaps found between sources

1. **Storage.** Part A 3.3 and part B ask for five new columns; GEN.440 and GEN.435 forbid them. Resolved by D9; the verified rule of the item wins.
2. **Population.** Part A treats the population as stored for every star; the code stores it on `bright_stars` only. Resolved by D10.
3. **Neon.** Part A 3.2 step 4 (magnesium) against the GEN.464 filing (oxygen). Resolved by a switch, D11.
4. **Rounded g.** Part A 2.4 prints 1.68 while its invariant 2 asks 1e-6; the exact value is 1.674.
5. **Second abundance fit.** gap-filling-plan 4.7 (Brewer 2016, 1592 stars, linear, no quality cut, sd 0.111, 0.127, 0.060, 0.087) differs from part A 2.7 (hinge fits, cuts, sd 0.04 to 0.07). At F = -0.5 the mixture mean [Mg/Fe] is 0.251 against 0.100, and at +0.3 the thin value 0.006 against -0.098 (the hinge holds a plateau). Differences up to 0.1 dex sit inside the 0.11 dex catalogue-to-catalogue spread (Hinkel 2014), so both ship with that limit; the part A fit is used because it separates the branches. The C/O median 0.47 (tunable-gaps f7, number) against 0.51 here (0.385 mass) is consistent with a different solar reference, not verified.
6. **CMF range.** GEN.428's 0.2 to 0.47 is the Thiabaud envelope; the drawn range is narrower (section 8).
7. **White-dwarf radius and neutron-star mass law** in the code (section 7).
8. **Brown dwarf limits.** Spiegel's 13.0 +- 0.8 MJ (abstract) and 13.55 MJ (Table 1, solar Z, Y 0.25) are two values of one source; the formula takes the Table 1 centre. Its solar Z is 0.0189 and must not be replaced by 0.0154.
9. **CEMP above probability 1 below -4** is unreachable under the clamp at -3.5.
10. **Radiogenic Th.** Part A step 7 derives Th from Eu over the whole range; master 15.1 supports Th/Fe trends only within +-0.1 dex of solar. Shipped with that limit.
11. **Lifetime.** The Sun reaches zero central hydrogen at 8 Gy in models, the code gives 10 Gy; decided (D2).
12. **Age gate.** The first draft's 3.0 Gy (the project's intermediate/old boundary) against Haywood 2013 (thick disc 13 to 8 Gy ago, uncertain by 1 to 1.5 Gy). The paper wins, both values stated: D10, Q3.
13. **Palme and O'Neill.** Our files' "2014" Th 0.0849 and U 0.0229 ppm against the 2003 edition read in Boss's papers (Th 83.4 ppb, U 21.8 ppb, primitive-mantle oxide total 98.41). The paper wins; label "2003 edition". This note takes nothing from either except the U/Th cross-check of step 10.
14. **Age-velocity relation.** Aumer and Binney 2009 (Eq. 16, Table 6: beta_U 0.307, beta_V 0.430, beta_W 0.445, 3D best fit beta 0.349, tau1 0.149, v10 55.179 km/s; disc age 12.6 Gyr): sigma_W/sigma_U runs 0.41 at 1 Gyr to 0.57 at 10 Gyr, which I recomputed from the table (0.414 and 0.569), not the constant 0.54 of the plan. This note does not use the relation, so nothing changes here; the thick-disc classification by velocity (Bensby's) is not generated.
15. **Not uploaded.** Bensby 2014 and Pinna 2019b are not among Boss's uploads and Pinna 2019a covers one galaxy; every Bensby-based value (the branch logistic, the Mg, Si, Ca, Al hinge fits) stays part A's and partly verified.

Not found anywhere in the research (open, for Boss's list): mass and metallicity dependence of main-sequence luminosity and temperature (isochrone grids), the white-dwarf C:O split and layer masses, the neutron-star EOS and crust depth, age dependence of photospheric diffusion, within-group abundance scatter, the age-resolved thick-disc fraction, the alpha pattern of satellites and the bulge, and yields for carbon stars and supergiants.

## 13. Source conflicts and the rule that decided each

Rules (Boss, 12:39 PM EDT on 2026-10-17): verified wins; if two verified sources are in doubt, the newest wins; if both are in question, most permissive wins. The third means the widest value or range, or both options kept with the spread carried: Boss confirmed this at 1:12 PM EDT (Decided), as composition-model.md section 4a item 12 records it.

| Conflict | Rule applied | Outcome |
|---|---|---|
| Age gate 3.0 Gy (project) against Haywood 2013 (thick disc 13 to 8 Gy ago) | Verified | 8.0 Gy, range 6.5 to 9.5 (D10, Q3) |
| Palme and O'Neill "2014" Th 0.0849, U 0.0229 ppm against the 2003 edition | Verified | Label "2003 edition"; 6 percent U/Th spread carried (step 10) |
| Brewer 2016 linear fits against part A's hinge fits (up to 0.1 dex on [Mg/Fe]) | Most permissive | Both ship with the 0.1 dex spread as the stated uncertainty; the hinge fit stays the default because only it separates the branches. A literal widest-scatter reading would use Brewer's sd of 0.06 to 0.13 instead of 0.04 to 0.07: not done, Boss to say. Newest would favour Brewer (2016 against 2014) |
| Neon follows O (GEN.464 filing) or Mg (part A 3.2 step 4) | Most permissive | Switch kept, both Z/Z0 columns carried (8 percent at most); default O (D11) |
| g = 1.68 (part A) against 1.674 (recomputed from X0, Z0, Y_P) | Verified | 1.674; Casagrande's 2.1 +- 0.9 is wider and contains it |
| WD radius law in the code against the integrated table (Sirius B) | Verified | Integrated table (D13) |
| NS radius 10 to 13 km uniform in the code against Miller 2021 hull 11.6 to 13.1 km | Verified | The band; J0740 alone (13.7 +2.6/-1.5) is wider and not a cap |
| Lifetime 10 Gy (code) against the Sun's 8 Gy hydrogen exhaustion (models) | Verified | Derived, not set (D2, decided) |
| Brown-dwarf deuterium limit 13.0 +- 0.8 MJ (Spiegel abstract) against 13.55 MJ (Table 1) | Verified | Table 1 |
| CEMP probability 1 below -4 (part A) against the clamp at -3.5 (code) | Verified | 0.25 over the reachable range |
| GEN.428 CMF 0.2 to 0.47 against the drawn 0.20 to 0.34 | Most permissive | Envelope kept as the bound, draw reported separately |
| Th from Eu across the range (part A) against Th measured only within +-0.1 dex of solar (master 15.1) | Most permissive | Shipped with the limit (partly verified) |
| C/O median 0.47 (tunable-gaps f7) against 0.51 here | Most permissive | Both kept; the solar reference differs (unverified) |
| Part A 3.3 and B's new columns against the item text of GEN.435 and GEN.440 | None: an instruction from the item, not a source value | No columns (D9) |
| Solar reference scale: Asplund 2021 against Lodders 2009 and 2019 | None: Boss's ruling, 1:24 AM EDT | Asplund 2021, spread as uncertainty (D1) |

## 14. Approved by Boss

**Decided (Boss, 2026-10-17, 9:40 to 9:45 PM EDT, cards answered in the project chat): every row below is approved as recommended, with the research row approved as a build-time download (GEN.548).** The table keeps the recommendation each card carried.

| # | Rule on | Recommendation (approved) |
|---|---|---|
| D9 | Derive the branch and four offsets from a seeded key, store nothing | Approve: no value fails the derivation test, and no migration is needed |
| D10 | Age gate at 8.0 Gy (range 6.5 to 9.5, Haywood 2013) instead of a stored population | Approve 8.0 with the logistic unscaled; the thick-disc share is then low (0.095 of old-disc stars) until Bensby 2014 or Pinna 2019b gives an age-resolved fraction; do not adopt the 3.5 times stand-in |
| D11 | Oxygen from the NLTE fit, neon follows oxygen, magnesium as a switch | Approve the default; neon's physics is open |
| D12 | Keep [Fe/H] independent of age | Approve; a link moves every stored value for a 0.28 correlation |
| D13 | Replace the white-dwarf radius law; observed neutron-star mass mixture; no black-hole composition | Approve; the white-dwarf change is a generator epoch change joined to the one rescatter |
| D14 | Return birth composition, no photospheric correction | Approve until an age law is read |
| D15 | No GEN.407 table changes here; tip by Z through GEN.502 | Approve |
| D16 | Shared [Fe/H], branch and offsets in groups and clusters, zero within-group scatter | Approve as a flagged default |
| Rule | Third source rule: the widest value or range, or both kept (Boss confirmed at 1:12 PM EDT that 'most permission take wins' means the widest range of the two sources: Decided) | Nothing to confirm; section 13 shows where it was used |
| Brewer scatter | Use Brewer 2016's wider sd (0.06 to 0.13) in place of part A's (0.04 to 0.07)? | No: keep part A, carry the 0.1 dex mean difference as uncertainty |
| Research | Isochrone grids for L, T and lifetime against Z; white-dwarf layers; neutron-star EOS | Boss approved a build-time download as in D6 (GEN.548) |

### 14.1 The three downloads (GEN.548, as built)

Boss's rule (12:47 AM EDT, 2026-10-18): static files are downloaded at setup, never at run time. `update.sh` and `install.sh` run `python3 -m planetgen.cli.datasets fetch` (gap-filling-plan.md 6.5); the readers only look in the cache folder and say so when they are on the stand-in. Every citation below was opened on 2026-10-18 (Crossref, arXiv, Zenodo records); the files are not in the repository.

| Fill | Data | Source (verified) | Setup download | Readers (`planetgen.datasets`) | What the fill replaced |
|---|---|---|---|---|---|
| SC1 | Main-sequence luminosity, temperature and radius against mass, age and [Fe/H]; main-sequence lifetime | MIST v1.2 basic isochrones, v/vcrit 0.4, solar-scaled, 15 values of [Fe/H] from -4.00 to +0.50, 107 ages (Choi et al. 2016, ApJ 823, 102; Dotter 2016, ApJS 222, 8; mist.science) | one 221.5 MB archive, checked by sha256, reduced to the 5 MB `main_sequence.npz` of main-sequence rows, then deleted | `isochrones.main_sequence_point`, `main_sequence_lifetime_gy`, `abundances`, `metallicities` | the Baraffe et al. 2015 track at solar metallicity and the fitted relations of section 6 |
| SC2 | White-dwarf radius against mass and temperature, and the hydrogen, helium, carbon and oxygen mass fractions | Bédard, Bergeron, Brassard and Fontaine 2020, ApJ 901, 93 (arXiv:2008.07469): 23 masses 0.20 to 1.30, thick and thin hydrogen layer, homogeneous equimassic C/O core; Warwick page of A. Bédard | `sequences.tar`, 2.5 MB | `white_dwarf_models.radius_solar`, `layer_mass_fractions`, `sequence`, `masses_sol` | the radius law of D13 and a fixed layer split |
| SC3 | 21 unified neutron-star equations of state; the ground-state outer-crust nuclei for 20 nuclear mass models | Malik, Pais and Providência 2024, A&A 689, A242 (Zenodo 10.5281/zenodo.10950421); Servais and Chamel 2026, A&A 710, A322 (Zenodo 10.5281/zenodo.20085135); both CC BY 4.0 | 21 EOS tables (10 KB each) and `compositions.zip`, 1.8 MB | `neutron_star_matter.equations_of_state`, `crust_models`, `outer_crust`, `crust_nucleus` | the nickel-to-krypton crust table and the radius floor of D13 |

What the data does not close: the white-dwarf sequences are carbon-oxygen models with a 1:1 core split by assumption (a helium-core or oxygen-neon white dwarf, and a measured C:O split, remain open); the isochrone grid is solar-scaled (no [alpha/Fe] enhancement, which the thick-disc branch of D9 needs) and its rotation is a model choice; the equations of state are model families, so the maximum mass and radius they give are spreads of models, not measurements. **As built (GEN.590), Boss's answer 'Require the data' (12:55 PM EDT on 2026-10-18).** `stellar_composition.star_makeup` reads the three readers and has no fallback: a main-sequence star gets `detail["main_sequence"]` (log luminosity, log Teff and log radius from the MIST track at its mass, age and [Fe/H], `None` when the track does not reach it) and `main_sequence_lifetime_gy`; a white dwarf needs its `teff_k` and gets the sequences' radius and, for a carbon-oxygen core, the hydrogen, helium, carbon and oxygen layers (a DA star takes the thick-layer case, any other the thin one: a choice of this build); a neutron star gets the outer crust of the mass model `u_crust_model` picks among the 20 (a draw appended after the earlier ones, so none moves), ending at the drip density that model gives. Without the files these kinds raise `DataMissing`, and bulk generation refuses to start (a "Check the data" step before the reset in the Generate page's jobs, `run_galaxy.require_datasets` in the CLI and the neighbourhood button). The Chandrasekhar radius table of GEN.439 and the Chamel and Haensel crust list are gone. Kept as stand-ins, because the data cannot replace them: the brown-dwarf radius and hydrogen limit (the MIST grid starts at 0.1 solar masses) and the inner crust and core of a neutron star, which need an equation of state integrated (nothing does yet; the radius is the measured band of `physics/neutron_star.py`). The helium-core and oxygen-neon white dwarfs get a radius from the carbon-oxygen sequences and no layers; both are listed in `open_items`.

### 14.2 The star's page (UX.254, as built)

The star's sheet (`describe/sheet.py`, the page and `GET /api/objects/<ref>/sheet`) shows its chemistry in the Composition group, worked out when the sheet is made from the stored class, state, mass, [Fe/H] and age by `star_makeup`: the **disc population** (thin disc, or thick disc or halo, alpha-rich), **[alpha/Fe]** (`stellar_abundances.alpha_enhancement`, the mean [X/Fe] of O, Mg, Si and Ca, in dex), the **metal mass fraction** Z the star formed with, and, for a red giant, the **helium added by the first dredge-up**. No column and no stored value: the same star reads the same on every request. The draws are keyed as the planets' host ratios are (`bulk_draw.star_key`: the stored mass in kg and the distance of its sector from the galactic centre, light-years, which `sectors.galactic_radius_pc` gives; the galaxy's default for a star in no placed sector), so the page and the planets agree. A star with no [Fe/H], a black hole or neutron star, or a kind whose downloaded data are missing shows none of the rows. A remnant whose mass changed after its planets were drawn is keyed on the mass it has now, and the planets on the mass they saw (an open item; the main-sequence and giant stars, which keep their mass, agree).

## Appendix: the code behind the numbers

```python
import math, random
EPS = dict(C=8.46, O=8.69, Mg=7.55, Si=7.51, Fe=7.46, Ca=6.30, Al=6.43)        # A21 log eps, part A 2.2
AW  = dict(Fe=55.845, Mg=24.305, Si=28.0855, C=12.011, O=15.999, Ca=40.078, Al=26.9815)
PROTO = dict(C=2.8356e-3, N=7.7519e-4, O=6.4144e-3, Ne=1.8966e-3, Na=3.1232e-5, Mg=7.0592e-4, Al=5.9447e-5,
             Si=7.4394e-4, S=3.4596e-4, Ar=7.8448e-5, K=3.7603e-6, Ca=6.5459e-5, Fe=1.3184e-3, Ni=7.6147e-5)
X0, Z0, Y_P = 0.7121, 0.0154, 0.24672
THICK_START_GY = 8.0                                                           # Haywood 2013, range 6.5 to 9.5 (D10)
G = (1 - Y_P - X0) / Z0 - 1                                                    # 1.674026
REST = Z0 - sum(PROTO.values())                                                # 4.95e-5, scaled with Fe
SD = dict(thin=dict(C=.044, O=.052, Mg=.066, Si=.042), high=dict(C=.050, O=.065, Mg=.073, Si=.073))
CORR = [[1,.56,.69,.71],[.56,1,.36,.37],[.69,.36,1,.77],[.71,.37,.77,1]]      # C, O, Mg, Si (BF16)

def xfe(el, F, br):                                                            # part A 2.7
    F = min(max(F, -.9), .5) if br == "thin" else min(max(F, -1.6), .1)
    t = br == "thin"
    return {"Mg": (.006 - .343*min(F-.03, 0)) if t else (.357 - .373*max(F+.73, 0)),
            "Si": (.017 - .166*min(F+.07, 0)) if t else (.253 - .199*max(F+.99, 0)),
            "Ca": (-.041 - .163*min(F-.42, 0)) if t else (.282 - .251*max(F+.99, 0)),
            "O":  (-.011 - .373*F) if t else (.606 - .495*max(F+.91, 0)),
            "C":  (-.030 - .163*F) if t else (.191 - .521*max(F+.10, 0))}[el]

def chol(a):
    n = len(a); L = [[0.0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1):
            s = sum(L[i][k]*L[j][k] for k in range(j))
            L[i][j] = math.sqrt(a[i][i]-s) if i == j else (a[i][j]-s)/L[j][j]
    return L

def derive(F, age, key, ne="O"):
    rng = random.Random(key); u = rng.random(); z = [rng.gauss(0, 1) for _ in range(4)]
    u_cemp, u_c = rng.random(), rng.random()
    P = 0.0 if age < THICK_START_GY else 1/(1+math.exp((F+0.350)/0.132))
    br = "high" if u < P else "thin"
    L = chol(CORR); cz = [sum(L[i][k]*z[k] for k in range(4)) for i in range(4)]
    d = {e: xfe(e, F, br) + cz[i]*SD[br][e] for i, e in enumerate(("C", "O", "Mg", "Si"))}
    if F < -2 and u_cemp < (0.25 if F > -4 else 1.0): d["C"] = 0.7 + 1.3*u_c    # CEMP
    d["Ca"] = xfe("Ca", F, br); d["S"] = d["Ar"] = d["Ca"]
    d["Ne"] = d["O"] if ne == "O" else d["Mg"]
    for e in ("N", "Na", "Al", "K", "Fe", "Ni"): d[e] = 0.0
    r = (sum(PROTO[e]*10**(F+d[e]) for e in PROTO) + REST*10**F) / X0          # Z/X
    X = (1-Y_P)/(1+(1+G)*r); Z = r*X
    return dict(branch=br, d=d, X=X, Y=1-X-Z, Z=Z)

def ratio(d, a, b):                                                             # mass ratio; F cancels, d[Fe] = 0
    return 10**(EPS[a]-EPS[b]+d[a]-d[b]) * AW[a]/AW[b]

def cmf(F, d, al):                                                              # al = [Al/Fe]; section 8
    n = lambda a, x: 10**(EPS[a]+F+x)
    mg = n("Mg", d["Mg"]); fe, si, ca, alu = n("Fe", 0)/mg, n("Si", d["Si"])/mg, n("Ca", d["Ca"])/mg, n("Al", al)/mg
    mubar = (AW["Mg"]+AW["O"]) + ca*(AW["Ca"]+AW["O"]) + alu*(AW["Al"]+1.5*AW["O"]) + si*(AW["Si"]+2*AW["O"])
    return fe*AW["Fe"] / (fe*AW["Fe"] + mubar)

# Remnants (CODATA 2018; Msun 1.98841e30 kg, Rsun 695.7e6 m as in `constants.py`)
G_, c_, hbar, h, me, mu = 6.6743e-11, 299792458.0, 1.054571817e-34, 2*math.pi*1.054571817e-34, 9.1093837015e-31, 1.66053906660e-27
def white_dwarf(xc):                      # zero-temperature Chandrasekhar star, mu_e = 2; integrate dx/dr and dm/dr to x = 1e-4
    from scipy.integrate import solve_ivp
    A = math.pi*me**4*c_**5/(3*h**3); B = 8*math.pi*mu*2*(me*c_)**3/(3*h**3)
    f = lambda r, y: [-G_*y[1]*B*y[0]**3/r**2 / (8*A*y[0]**4/math.sqrt(1+y[0]**2)), 4*math.pi*r*r*B*y[0]**3]
    ev = lambda r, y: y[0]-1e-4; ev.terminal = True
    s = solve_ivp(f, [1, 5e7], [xc, 4/3*math.pi*B*xc**3], events=ev, rtol=1e-10, atol=[1e-12, 1e10], method="LSODA")
    return s.y[1][-1]/1.98841e30, s.t[-1]/695.7e6             # Msun, Rsun; a grid of 80 xc in 0.3..12 is interpolated
def black_hole(m, a):                     # m in Msun; returns r_s km, T_H K, r_+ km, ISCO in GM/c^2, efficiency
    M = m*1.98841e30; rg = G_*M/c_**2
    z1 = 1+(1-a*a)**(1/3)*((1+a)**(1/3)+(1-a)**(1/3)); z2 = math.sqrt(3*a*a+z1*z1)
    isco = 3+z2-math.sqrt((3-z1)*(3+z1+2*z2))
    return 2*rg/1e3, hbar*c_**3/(8*math.pi*G_*M*1.380649e-23), rg*(1+math.sqrt(1-a*a))/1e3, isco, 1-math.sqrt(1-2/(3*isco))

def birth_u_th(th_ppm, u_ppm, t=4.568):                                         # mass ratio U/Th t Gy back; 235U atom fraction 0.7204 percent
    l = lambda h: math.log(2)/h
    a8, a5, m8, m5 = 0.992742, 0.007204, 238.0508, 235.0439
    fu = (a8*math.exp(l(4.468)*t)*m8 + a5*math.exp(l(0.7038)*t)*m5) / (a8*m8 + a5*m5)
    return u_ppm*fu / (th_ppm*math.exp(l(14.05)*t))                             # 2003 edition 0.0834, 0.0218: 0.554; "2014" labels: 0.572
```
