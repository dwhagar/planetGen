# Giant planet, ice giant and brown dwarf interiors: the layers, metallic hydrogen, helium rain and the core (GEN.641)

**Status (2026-10-23):** research note for GEN.641 and the design GEN.642 (draw and describe) and DB.89 (store) build from. Boss, 6:01 PM EDT 2026-10-23: gas giants get a full composition, including liquid or metallic hydrogen and helium, whether there is a rocky core, and what it is made of. It changes no stored value and no draw. Sibling notes: [bulk-composition.md](bulk-composition.md) (decisions 20 to 24, 34: the envelope, the heavy-element mass law, the radius tables), [composition-model.md](composition-model.md), [body-composition-parts/C-volatile-and-giants.md](body-composition-parts/C-volatile-and-giants.md) (the radius laws and the Fortney tables). The whole-body totals stay as GEN.428 stores them (`envelope_mass_fraction`, `envelope_helium_fraction`, `envelope_z_fraction`); this note adds the interior inside them.

Verification states follow the rule of 3:02 AM EDT 2026-10-17: **READ** (the paper's text was read for the number), **READ (abstract)**, **2ND-HAND** (a number the read papers quote from another paper, which was not opened), **STAND-IN** (a Solar System fit or a placeholder, flagged, with the item that replaces it), **OPEN** (no source; a seeded range is used and flagged). Conflicts follow Boss's rules (bulk-composition.md section 11): verified wins; two verified in doubt, the newest wins; both in question, the widest range wins. Nothing is fixed per planet and nothing is copied from the Solar System as a preset: Jupiter, Saturn, Uranus and Neptune are test targets (section 6).

## 1. What exists at main 8.0.2042

- A giant has a class letter and a composition text; no layer, no core mass, no helium split. `generation/planet.py` writes flavour text about a "liquid metallic state" using the n = 1 polytrope central pressure P_c = pi G M^2 / (8 R^4) (`physics/constants.py` line 195) and one fixed `HYDROGEN_METALLIZATION_PRESSURE_PA = 1e11` ("a narrative milestone, not a claim of precision"). Nothing stores it.
- The shared materials table (`physics/data/materials.json`, GEN.427) has the rows `hydrogen`, `helium`, `rock`, `iron`, `iron_epsilon`, `water_liquid`, `ice_vii`, `ice_x` (the row named "ice X (and superionic water above about 100 GPa and 2000 K)") and the other ices. It has no row for metallic hydrogen; section 5 adds the layer states, not a new material, because metallic hydrogen is hydrogen at a pressure.
- Decisions 20 to 24 and 34 of bulk-composition.md fix the envelope fraction, the heavy-element mass M_Z law (Chachan mean, Thorngren selectable), a core drawn log-uniform 5 to 20 Earth masses (decision 22, OPEN) and the Fortney radius tables. They treat a giant as a core plus one homogeneous envelope. The sources below show that picture is out of date; section 3 replaces decision 22.

## 2. What the sources say

All numbers are for the four Solar System giants, which are the only ones measured, so each is a constraint on a law, not a preset.

### 2.1 Hydrogen: molecular to metallic

| Quantity | Value | Source, state |
|---|---|---|
| Insulating molecular liquid to metallic liquid, above about 1 Mbar at high temperature | stated as the general picture | Helled and Howard 2024 (arXiv 2407.05853), section 2.0.1: READ |
| Crossover in Jupiter's conditions, pure H and H-He mixture of proto-solar composition, 6000 K isotherm, quantum Monte Carlo | 0.4 to 0.6 Mbar (40 to 60 GPa); helium delays metallization by about 70 GPa at 1500 K and 30 GPa at 6000 K | Mazzola, Helled and Sorella 2018 (arXiv 1709.08648), abstract and section 3; Helled, Mazzola and Redmer 2020 (arXiv 2006.12219) section 3: READ |
| Older and other estimates | liquid-liquid transition at 1 to 2 Mbar below a critical temperature of about 2000 to 5000 K; experiment 1.4 Mbar (Weir 1996) | Debras and Chabrier 2019 (arXiv 1901.05697) section 3: READ; Weir via Helled 2020: 2ND-HAND |
| Magnetic field generation in Jupiter | starts near 140 GPa | Helled 2020 section 1: READ |

Rule 3 (widest range) gives the **molecular-to-metallic pressure P_MM drawn log-uniform from 0.4 to 2 Mbar** once per body (decision D1). The fixed 1 Mbar constant of `constants.py` goes (no shim).

### 2.2 Helium rain

| Quantity | Value | Source, state |
|---|---|---|
| Helium is immiscible in metallic hydrogen below a temperature that rises with pressure; droplets fall and release gravitational energy | the mechanism | Stevenson and Salpeter 1977, via Helled and Howard 2024 section 2.0.1 and Mankovich and Fortney 2020: 2ND-HAND |
| Onset for a proto-solar adiabat (Y = 0.27) with the Schottler and Redmer 2018 phase curve | about 5500 K at 2 Mbar, for a one-bar temperature T1 = 149 K | Mankovich and Fortney 2020 (arXiv 1912.01009) Figure 1: READ |
| That curve taken at face value | Jupiter (T1 = 166 K) is 500 K or more too warm to demix, yet its atmosphere is helium-depleted; good solutions need the phase curve shifted warmer (delta T_phase up to about 400 K explored) | the same paper, section 4: READ |
| Jupiter's atmosphere | Y / (X + Y) = 0.234 (Galileo) against proto-solar 0.2745; "modest" differentiation, begins at about 3.5 Gyr | Helled and Howard 2024 Figure 4 caption: READ; Galileo value 2ND-HAND |
| Saturn's atmosphere | predicted Y = 0.07 +/- 0.01 (He/H2 = 0.036 +/- 0.006); begins at about 1.3 Gyr; the luminosity from rain extends the cooling time to the solar age | Mankovich and Fortney 2020 abstract: READ |
| Where it happens | from the 2 Mbar level, at about 80 percent of Jupiter's radius and 55 percent of Saturn's, with a continuous helium gradient below | Mankovich and Fortney 2020 Figure 3 caption: READ |
| Saturn's deep interior | a helium-rich shell (Y near 1 at the bottom) above the heavy-element core, inner boundary near 35 percent of Saturn's radius | the same, Figure 3 and section 2: READ; Helled and Howard 2024 section 3.1.2: READ |
| Neon follows helium | neon differentiation may have added luminosity | Mankovich and Fortney 2020 abstract: READ (abstract) |

Cooler planets differentiate earlier and deeper: a Saturn-mass planet at a given age rains more than a Jupiter-mass one. The phase diagram is the weakest input: Helled and Stevenson 2024 say it is "very difficult to determine the critical line of a mixture at megabar pressures and many thousands of degrees".

### 2.3 The core: compact, dilute, or both

| Quantity | Value | Source, state |
|---|---|---|
| The old picture (pure heavy-element core under one homogeneous envelope) | "no longer valid"; gravity data need a deep region enriched in heavy elements with a gradient, a "fuzzy" or "dilute" core | Helled and Stevenson 2024 (arXiv 2403.11657) sections 1 and 2: READ |
| Jupiter's dilute core | mass about 100 Earth masses (a third of Jupiter), H-He-dominated, heavy-element fraction about 20 percent; small-core solutions extend to about 20 percent of Jupiter's mass | Helled and Stevenson 2024 section 2.4 and Figure 2: READ; Howard et al. 2023 (arXiv 2302.09082) abstract: READ (abstract) |
| Jupiter, other models | dilute core out to 60 to 75 percent of the radius; region between 0.1 and 2 Mbar with an entropy and composition change; outer envelope near solar Z | Debras and Chabrier 2019 abstract and Figures 3 to 4: READ |
| Jupiter's total heavy elements | 18 to 45 Earth masses (adiabatic about 20, non-adiabatic up to 60) | Helled and Howard 2024 section 3.1.1; Helled and Stevenson 2024 section 2.1: READ |
| Saturn's dilute core | stable region out to 0.59 +/- 0.01 R_S, 55.1 +/- 1.7 Earth masses of which 17.4 +/- 1.2 are rock and ice (Z about 0.32 there); central density about 6 g/cm3; up to 80 percent heavy elements in the inner part | Mankovich and Fuller 2021 (arXiv 2104.13385) text: READ; Helled and Howard 2024 section 3.1.2: READ |
| Saturn, compact-core models | pure core of 15 to 18 Earth masses with under 5 Earth masses in the envelope; 4-layer fit 12 to 18 Earth masses total | Militzer et al. 2019, Ni 2020, via Helled and Howard 2024: 2ND-HAND |
| A compact pure heavy-element core | "a few Earth masses at most" for Jupiter; may not exist; invisible to gravity (a 5 Earth-mass core adds about 1e-4 of the moment of inertia); formation models predict it, and convective mixing or a giant impact could erase it | Helled and Stevenson 2024 section 2.3: READ |
| Why the core is dilute | four candidate origins: evaporated accreted solids mixing in (confined to about the inner 20 Earth masses by runaway gas accretion), giant impacts, convective erosion, and delayed runaway accretion at about 100 Earth masses | Helled and Stevenson 2024 section 3: READ |
| What it is made of | water and rock mix with metallic hydrogen in all proportions at these temperatures; gravity cannot tell water from rock; rocks may sink deeper than water during accretion, so a compact core is rock-and-iron rich and the ice mixes higher up | Helled and Stevenson 2024 section 4 ("Differentiation", "The chemical composition of heavy elements"): READ; Lozovsky 2017 via it: 2ND-HAND |
| Temperature | central temperature about 20,000 K adiabatic, 35,000 K for superadiabatic (Jupiter); 10 Mbar region 9,500 to 13,000 K; boundary 1 to 2 Mbar at 5,000 to 9,000 K; heavy-element mass and central temperature trade off | Helled and Howard 2024 section 3.1.1; Debras and Chabrier 2019 Figure 4: READ |

Conflict (rule 3): compact core 0 to a few Earth masses (Jupiter, Helled 2024) against 15 to 18 (Saturn, Militzer 2019, 2ND-HAND). **The compact core is drawn between 0 and 20 Earth masses** (the older decision 22 range of 5 to 20 reaches the Saturn value; the lower end now includes none), D3.

### 2.4 Ice giants

| Quantity | Value | Source, state |
|---|---|---|
| H/He share of the mass | about 10 to 20 percent (Helled and Fortney 2020); up to about 2 Earth masses in three-layer models | Helled and Howard 2024 section 3.2: READ |
| Heavy-element share | Uranus 76 to 90 percent, Neptune 77 to 90 percent, depending on whether the heavy element is rock (SiO2) or water | Helled et al. 2011 (arXiv 1010.5546) via Helled and Howard 2024: READ |
| Structure | three layers (rock core, water-rich mantle, H/He envelope) or gradients with no distinct layers; outer envelope Z up to 0.17 (Voyager rotation) or 0.08 (revised), inner envelope Z about 0.9; solar-Z envelopes also fit | Nettelmann et al. 2013 via Helled and Howard 2024: 2ND-HAND |
| Rock to water ratio | unknown: at most about 2 (Podolak 2022); 0.05 to 0.38 maximum for non-adiabatic empirical models (Neuenschwander 2024); solar 1/2; carbon-rich planetesimals give more methane and water than rock (Malamud 2024: methane planets) | Helled and Howard 2024 sections 3.2.1 and 3.2.2: READ; each paper 2ND-HAND |
| Is the hydrogen metallic | the adiabats cross the molecular-to-metallic line, but "it is unclear whether H-He exist in the planetary deep interiors" | Helled and Howard 2024 Figure 1 caption: READ. Estimate here: the pressure at the base of a 2 Earth-mass envelope on a 13 Earth-mass core is about 21 GPa (0.2 Mbar), below P_MM, so no metallic layer (a thin-shell estimate, `G M_core M_env / (4 pi R_core^4)`, COMPUTED) |
| Water in the deep interior | superionic ice (ice XVIII) at 100 to 400 GPa and 2000 to 3000 K; solid oxygen lattice, mobile protons; the generator needs its density and that it is solid | C-volatile-and-giants.md section 2.3 (Millot 2019 abstract, French 2009 via Zeng and Sasselov): READ (abstract), 2ND-HAND |
| Low luminosity | composition gradients and thermal boundary layers; central temperatures 5,000 to 6,500 K adiabatic, above 20,000 K non-adiabatic | Helled and Howard 2024 section 3.2.1: READ |

### 2.5 Brown dwarfs

A brown dwarf is fully convective and partially degenerate (Chabrier and Baraffe 2000, astro-ph/0006383, sections 2 and 3: READ for the degeneracy parameter psi running from about 2 to 0.05 and the central density rising as mass falls toward the hydrogen-burning limit). Its central pressure is hundreds to thousands of times P_MM, so only a thin outer layer is molecular; it forms from a collapsing cloud with no solid core, so its heavy-element mass fraction is the star's metallicity Z with no compact core (D6). Whether helium rain occurs inside cool T and Y dwarfs is OPEN (no source was read); the same rule as section 2.2 applies if it does, with T1 replaced by the effective temperature, and is **not** switched on for brown dwarfs until a source is read (GEN.641 follow-up).

## 3. The layers a giant gets

Layers are listed outside in. A layer has a kind, a pressure at its top and bottom, a share of the whole body's mass, the mass fractions of its materials (rows of the shared materials table) and a state.

| Layer | Kind key | Material rows | State | Present when |
|---|---|---|---|---|
| Molecular envelope | `molecular_envelope` | `hydrogen`, `helium` (Y_out), heavy elements Z_out | molecular fluid | always for a giant |
| Metallic envelope | `metallic_envelope` | `hydrogen`, `helium` (Y_in), heavy elements Z_mid | metallic fluid | envelope base pressure above P_MM |
| Helium rain zone | `helium_gradient` | as above with Y rising inward | metallic fluid, gradient | T1 below the onset (section 4.2) |
| Helium shell | `helium_shell` | `helium`, `hydrogen` | metallic fluid | strong depletion only (section 4.2) |
| Dilute core | `dilute_core` | `hydrogen`, `helium`, `rock`, `iron`, `water_liquid` or `ice_x`, other ices | metallic fluid with heavy elements, gradient | mass above the compact core and M_Z large enough (section 4.3) |
| Compact core | `compact_core` | `rock` (silicate), `iron`, with the ice share the draw gives | molten or solid | drawn (D3) |
| Ice mantle (ice giants) | `ice_mantle` | `water_liquid` or `ice_x` (superionic), `ammonia`, `methane` per the ice set of GEN.428, `rock` share | ionic or superionic fluid | envelope fraction below one half |
| Rock core (ice giants) | `rock_core` | `rock`, `iron` | molten or solid | always for an ice giant |

An ice giant uses `molecular_envelope`, `ice_mantle` and `rock_core` and has a `metallic_envelope` only if its envelope base pressure passes P_MM. A gas giant uses `molecular_envelope`, `metallic_envelope`, optionally the helium layers, `dilute_core` and optionally `compact_core`. Which of the two a body is follows from the envelope fraction f_env (a gas giant above one half, an ice giant below; Helled and Fortney's 10 to 20 percent for the Solar System ice giants) and not from a class letter (GEN.645).

## 4. The draws

Every number below is OPEN or STAND-IN unless it names a source; each draw is a gated GEN.472 table row set, registered as `giant-*`, with the bracket it comes from. The draw order is fixed (a reorder is an epoch change) and follows GEN.428 and GEN.460: M, M_Z, f_env and Y_proto come first (decisions 20, 21, 24).

### 4.1 Where the layers sit

The hydrostatic profile of the layers is the n = 1 polytrope that the repo already uses for the central pressure: P(xi)/P_c = (sin xi / xi)^2 with xi = pi r / R, and the mass inside xi is M (sin xi - xi cos xi)/pi, P_c = pi G M^2 / (8 R^4). It is a rough shape for a hydrogen-helium sphere, good enough to place layer boundaries to a few percent of the radius (COMPUTED here: for Jupiter's mass and radius the 2 Mbar level falls at 0.81 R_J and 16 percent of the mass lies above it, against 0.8 to 0.9 R_J in Debras and Chabrier; for Saturn 0.59 R_S against 0.55 in Mankovich and Fortney). It is a stand-in until the radius tables of GEN.460 supply real density profiles. For an envelope under 30 percent of the mass (ice giants, sub-Saturns) the shell formula `P_base = G M_core M_env / (4 pi R_core^4)` replaces it, because the mass is in the core. The radius R is the one GEN.429 and GEN.460 derive, so the layers never change mass or radius.

| Level | Pressure | Mass above, Jupiter | Mass above, Saturn | Radius, Jupiter | Radius, Saturn |
|---|---|---|---|---|---|
| 0.4 Mbar | lowest P_MM | 0.039 | 0.165 | 0.907 R | 0.801 R |
| 1 Mbar | middle | 0.088 | 0.332 | 0.859 R | 0.700 R |
| 2 Mbar | highest P_MM and the rain onset | 0.156 | 0.534 | 0.807 R | 0.588 R |

### 4.2 Metallic hydrogen and helium

- **P_MM** log-uniform 0.4 to 2 Mbar, drawn once (D1). The molecular envelope is the mass above it; the metallic envelope is the rest of the envelope down to the dilute core.
- **Helium split.** Rain is active when the planet's one-bar temperature T1 (the stored equilibrium plus internal heat; the body's own value, never a Solar System number) is below the onset T1_on (STAND-IN, D2). The onset is the one number the sources cannot give: the phase curve says 149 K unshifted and Jupiter shows depletion at 166 K, so the phase curve must sit warmer (Mankovich and Fortney section 4). Two measured points give a straight line between them: Y_out / Y_proto = 0.85 at T1 = 166 K (Jupiter) and 0.26 at T1 = 140 K (Saturn), reaching 1 at 172.6 K and held at 0.26 below 140 K. These two points are the only calibration; they are flagged in the gap plan registry as a stand-in that GEN.641's follow-up replaces with the Schottler and Redmer curve and the stored adiabat (the curve is Figure 1 of Mankovich and Fortney; the paper itself is to be fetched).
- **Mass balance.** Helium is conserved: Y_proto M_env = Y_out M_mol + Y_in M_met + M_He,shell, so the helium removed from the molecular envelope reappears below it. The metallic envelope gets the removed helium as a gradient (Y_in rising inward, the `helium_gradient` layer) while Y_out / Y_proto stays above 0.5; below 0.5 half of the removed helium goes into a `helium_shell` between the rain zone and the dilute core (the Saturn picture, "inevitably forms a helium-rich shell", STAND-IN D2).
- **Neon** follows helium in the sources; it is not a layer here (OPEN, no number).

### 4.3 The core

- **Compact core mass M_cc**, log-uniform from 0.5 to 20 Earth masses with probability 0.5 that none exists, capped at M_Z (D3, OPEN: the probability is a placeholder; Helled and Stevenson say it is unknown and that gravity cannot see it).
- **Dilute core mass M_dil**: log-uniform from 20 to 100 Earth masses (the evaporated-solid mixing zone of Muller 2020 at the low end, the delayed-runaway mass of Venturini and Helled at the high end, both 2ND-HAND through Helled and Stevenson), capped at 0.65 M (Saturn's 55 of 95 Earth masses is 0.58, READ). Below M = 30 Earth masses there is no separate dilute core and the heavy elements are a gradient through the whole body (the ice-giant picture).
- **Heavy-element budget.** M_Z (GEN.460) is split: Z_out M_mol in the outer envelope (Z_out = Z_env from decision 21), M_cc in the compact core, and the rest in the dilute core and the metallic envelope. The dilute core's mean Z_dil is that rest over M_dil and must fall between Z_out and 0.8 (Saturn's inner part, READ); if it does not, M_dil is moved inside its bracket until it does, and when no value fits the planet has no dilute core and the heavy elements are spread through the metallic envelope. Z rises inward through the dilute core; the shape of the gradient is OPEN (a linear rise from Z_out at its top is the placeholder) and only the mean, which the budget fixes, is a constraint.
- **What the heavy elements are made of.** The rock-to-ice split of the heavy elements comes from the formation branch GEN.428 stored: a body born beyond the snow line takes the water share W of its solids (decision 2, 0.14 to 0.39 of the solids) and the ice set of decision 18; a body born inside takes none. The compact core takes the rock and iron first (rocks sink deeper than water, Helled and Stevenson section 4), and the ice goes outward first: the molecular envelope (which holds it as vapour), then the metallic envelope, then the dilute core. Rock and iron that the compact core cannot hold fill the dilute core, then the deep metallic envelope. Heavy elements the dilute core cannot hold (it takes at most 80 percent of its mass) enrich the deep envelope, never the atmosphere, whose metallicity the caller fixed. The iron-to-rock split of the refractory part is the star's, as for a terrestrial planet (GEN.428's `metal_mass_fraction` over `rock_mass_fraction`, CMF 0.2917 for solar abundances).
- **Temperature** is not drawn or stored (OPEN): the boundary temperatures of Debras and Chabrier (5,000 to 9,000 K at 1 to 2 Mbar, 9,500 to 13,000 K at 3 to 7 Mbar, 19,000 to 30,000 K at 40 Mbar) are test targets, not inputs.

### 4.4 Ice giants

- f_env is GEN.428's; the H/He share is 0.05 to 0.20 of the mass (Helled and Fortney 10 to 20 percent for Uranus and Neptune; widest range with the three-layer fits that give 1 to 2 Earth masses).
- **Rock fraction of the heavy elements** 0.05 to 0.67 (Neuenschwander 2024 low, Podolak 2022 and Zeng 2019 high; rule 3 widest, OPEN). Water, ammonia and methane share the rest by the ice set (GEN.428), which already follows Marboeuf's table; the carbon-rich option of Malamud 2024 (methane above water) is a table row with weight set by the share of refractory organics, OPEN.
- Mantle water is `water_liquid` above a drawn depth and `ice_x` (superionic) when the layer's pressure exceeds 100 GPa; the radius does not change (the radius tables already carry it, decision 12).
- A gradient (no distinct layers) is the default representation of Z; the three layers are the rows a description names.

### 4.5 Brown dwarfs

One envelope: molecular above P_MM, metallic below, Z the star's metallicity, no compact core, no dilute core (D6). The molecular layer is thin, because the central pressure is far above P_MM; its mass share is computed when GEN.642 is built (an n = 3/2 polytrope for the degenerate range), not in this note. Helium is the proto-solar Y of the natal gas (stellar-composition.md) until a source for cool-dwarf helium rain is read.

## 5. What gets stored and described

- **Summary columns** on the planet and moon tables (DB.89, lane I's BulkFractions extends): `core_mass_earth` (compact core, 0 if none), `dilute_core_mass_earth`, `metallization_pressure_gpa` (P_MM), `helium_depletion` (Y_out / Y_proto, 1 when no rain), `has_helium_shell` (bool). All numeric columns are declared in `db/column_units.py` in the same PR.
- **Layer rows**: the child table that stores the detailed mixture (`planet_composition` / `moon_composition`, GEN.428) gains a `layer` index and a `layer_kind`; layer 0 is the whole-body mixture that GEN.428 writes today. Each layer row also stores `top_pressure_gpa`, `bottom_pressure_gpa`, `mass_fraction` (of the body) and `state`. A recalculation never redraws a stored layer (F10).
- **Description** (GEN.642, one writer with GEN.430; the numbers are an illustration of the wording, not a draw): "A gas giant of 1.0 Jupiter masses. Its outer 6 percent is molecular hydrogen with helium at 85 percent of the natal share; below 0.9 Mbar the hydrogen is metallic. About 40 Earth masses of heavier elements are spread through a dilute core that reaches 25 percent of the mass; the heavy elements there are rock and iron with 5 percent water ice. No compact core." Thresholds by fixed rules, as GEN.430.

## 6. Test targets (never presets)

Each target is met by some draw inside the allowed distributions, not by one draw, and lives in test files only.

| Target | Tolerance | Source |
|---|---|---|
| Jupiter-class mass and age: P_MM level at 0.80 to 0.93 R, dilute core 15 to 35 percent of the mass, Z_out under 0.1, total heavy elements 18 to 60 Earth masses | inside the bracket | Debras and Chabrier 2019; Helled and Stevenson 2024 |
| Saturn-class (95 Earth masses, T1 near 140 K): Y_out / Y_proto 0.15 to 0.4, a helium shell present in some draws, dilute core out to 0.4 to 0.65 R | inside the bracket | Mankovich and Fortney 2020; Mankovich and Fuller 2021 |
| Uranus and Neptune class: H/He 0.05 to 0.20, no metallic layer, Z 0.75 to 0.95, rock fraction of Z 0.05 to 0.67 | inside the bracket | Helled and Fortney 2020 |
| Layer masses sum to 1 and layer pressures are monotone; layer heavy-element masses sum to M_Z (1e-9) | exact | invariants |
| Helium is conserved across the layers (1e-9) | exact | invariants |
| A giant with envelope under one half has no `metallic_envelope` unless its base pressure passes P_MM | exact | section 3 |
| Same seed, same layers; a recalculation changes none (F10) | exact | bulk-composition.md section 8 |

## 7. Decisions

| # | Decision | State | Used meanwhile, source |
|---|---|---|---|
| D1 | The molecular-to-metallic pressure is drawn per body, log-uniform 0.4 to 2 Mbar; the fixed 1 Mbar goes | Default taken | Mazzola 2018 (0.4 to 0.6), Debras and Chabrier 2019 (1 to 2), rule 3 widest |
| D2 | Helium rain by the one-bar temperature: a line through two measured points, a helium shell at strong depletion | STAND-IN, flagged in the gap plan registry | Mankovich and Fortney 2020 give the mechanism and the 149 K unshifted onset; the line and the shell rule are placeholders until the Schottler and Redmer curve and the adiabat replace them |
| D3 | Compact core 0 to 20 Earth masses, probability 0.5 of none, replaces decision 22's 5 to 20 | Open | Helled and Stevenson 2024 (a few Earth masses, may not exist); Militzer 2019 (15 to 18) |
| D4 | A dilute core of 20 to 100 Earth masses, at most 0.65 of the mass, gradient in Z | Default taken | Debras and Chabrier 2019, Mankovich and Fuller 2021, Helled and Stevenson 2024 |
| D5 | Rock and iron to the compact core first, then the dilute core; ice outward first (molecular envelope, metallic envelope, dilute core); split from GEN.428's formation branch | Default taken | Helled and Stevenson 2024 section 4 |
| D6 | A brown dwarf has no compact or dilute core; metallic envelope below P_MM; no helium rain until a source is read | Default taken, rain open | Chabrier and Baraffe 2000 |
| D7 | The layer pressures come from the n = 1 polytrope (shell formula under 30 percent envelope) until GEN.460's tables supply real profiles | STAND-IN | `constants.py` line 195; COMPUTED checks in 4.1 |
| D8 | Layer rows live in the composition child table with a layer index; summary columns extend BulkFractions | Default taken | contracts file of lane I, section 1.1 |
| D9 | Whole-body accounting: heavy elements in envelope layers count as `envelope_z` of the envelope; heavy elements in core layers count as metal, rock, water and ice; all hydrogen and helium (the dilute core's too) is envelope; the helium share is Y/(X+Y) of the hydrogen and helium | Default taken | Lane I's `BulkFractions` contract |

## 8. Gaps

1. **The Schottler and Redmer 2018 phase curve** and the Morales 2013 curve are not in hand; their data (Figure 1 of Mankovich and Fortney is a plot of them) would replace D2's line. To fetch.
2. **Equation-of-state tables** for hydrogen-helium (the CMS and Saumon-Chabrier-van Horn tables) are not downloaded; the polytrope is the stand-in (D7).
3. **Cool brown dwarf helium rain**, neon rain, and the helium content of Y dwarfs: no source read.
4. **The probability of a compact core** and the gradient shape in the dilute core: no source; the placeholders are flagged.
5. **Mazzola 2018 crossover** is at 6000 K; the pressure at 3000 to 5000 K has no number in the read papers, so the bracket 0.4 to 2 Mbar is temperature-independent.
6. **Exoplanet giants**: inflated hot Jupiters may have a radiative layer or a different entropy; none of the sources read covers them (GEN.460's inflation row stays the radius authority).
7. **Hot and warm Jupiters' T1 for the rain test**: the equilibrium temperature with internal heat, not the Solar System's 166 K; the stored `temperature` is a surface value and the one-bar value needs the greenhouse and internal-heat terms of the atmosphere lane.

## 9. Sources fetched for this note

arXiv 2407.05853 (Helled and Howard 2024, full), 2403.11657 (Helled and Stevenson 2024, full), 1912.01009 (Mankovich and Fortney 2020, abstract, Figures 1 and 3, sections 2 to 4), 1901.05697 (Debras and Chabrier 2019, abstract and structure figures), 1709.08648 (Mazzola, Helled and Sorella 2018, abstract and metallization text), 2006.12219 (Helled, Mazzola and Redmer 2020, sections 3 and 4), 2104.13385 (Mankovich and Fuller 2021, abstract and mass numbers), 2302.09082 (Howard et al. 2023, abstract), 1010.5546 (Helled et al. 2011, summary), astro-ph/0006383 (Chabrier and Baraffe 2000, degeneracy sections). Raw papers stay out of the repository (Boss, 2026-10-17).
