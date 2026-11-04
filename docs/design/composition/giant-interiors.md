# Giant planet, ice giant and brown dwarf interiors: the layers, metallic hydrogen, helium rain and the core (GEN.641)

**Status (2026-10-23):** research note for GEN.641 and the design GEN.642 (draw and describe) and DB.89 (store) build from. Boss, 6:01 PM EDT 2026-10-23: gas giants get a full composition, including liquid or metallic hydrogen and helium, whether there is a rocky core, and what it is made of. The kernel is wired into the generator since generator epoch 80 (GEN.460, GEN.642; section 10); DB.89 (the layer table) is next. Sibling notes: [bulk-composition.md](bulk-composition.md) (decisions 20 to 24, 34: the envelope, the heavy-element mass law, the radius tables), [composition-model.md](composition-model.md), [body-composition-parts/C-volatile-and-giants.md](body-composition-parts/C-volatile-and-giants.md) (the radius laws and the Fortney tables). The whole-body totals stay as GEN.428 stores them (`envelope_mass_fraction`, `envelope_helium_fraction`, `envelope_z_fraction`); this note adds the interior inside them.

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
| Insulating molecular liquid to metallic liquid, above about 1 Mbar at high temperature | stated as the general picture | [R1893], section 2.0.1: READ |
| Crossover in Jupiter's conditions, pure H and H-He mixture of proto-solar composition, 6000 K isotherm, quantum Monte Carlo | 0.4 to 0.6 Mbar (40 to 60 GPa); helium delays metallization by about 70 GPa at 1500 K and 30 GPa at 6000 K | [R2966], abstract and section 3; [R1894] section 3: READ |
| Older and other estimates | liquid-liquid transition at 1 to 2 Mbar below a critical temperature of about 2000 to 5000 K; experiment 1.4 Mbar | [R1038] section 3: READ; Weir via [R1892, R1894]: 2ND-HAND |
| Magnetic field generation in Jupiter | starts near 140 GPa | [R1892, R1894] section 1: READ |

Rule 3 (widest range) gives the **molecular-to-metallic pressure P_MM drawn log-uniform from 0.4 to 2 Mbar** once per body (decision D1). The fixed 1 Mbar constant of `constants.py` goes (no shim).

### 2.2 Helium rain

| Quantity | Value | Source, state |
|---|---|---|
| Helium is immiscible in metallic hydrogen below a temperature that rises with pressure; droplets fall and release gravitational energy | the mechanism | [R4220], via [R1893] section 2.0.1 and [R2852]: 2ND-HAND |
| Onset for a proto-solar adiabat (Y = 0.27) with the [R3973] phase curve | about 5500 K at 2 Mbar, for a one-bar temperature T1 = 149 K | [R2852: Figure 1]: READ |
| That curve taken at face value | Jupiter (T1 = 166 K) is 500 K or more too warm to demix, yet its atmosphere is helium-depleted; good solutions need the phase curve shifted warmer (delta T_phase up to about 400 K explored) | the same paper, section 4: READ |
| Jupiter's atmosphere | Y / (X + Y) = 0.234 (Galileo) against proto-solar 0.2745; "modest" differentiation, begins at about 3.5 Gyr | [R1893: Figure 4] caption: READ; Galileo value 2ND-HAND |
| Saturn's atmosphere | predicted Y = 0.07 +/- 0.01 (He/H2 = 0.036 +/- 0.006); begins at about 1.3 Gyr; the luminosity from rain extends the cooling time to the solar age | [R2852] abstract: READ |
| Where it happens | from the 2 Mbar level, at about 80 percent of Jupiter's radius and 55 percent of Saturn's, with a continuous helium gradient below | [R2852: Figure 3] caption: READ |
| Saturn's deep interior | a helium-rich shell (Y near 1 at the bottom) above the heavy-element core, inner boundary near 35 percent of Saturn's radius | the same, Figure 3 and section 2: READ; [R1893] section 3.1.2: READ |
| Neon follows helium | neon differentiation may have added luminosity | [R2852] abstract: READ (abstract) |

Cooler planets differentiate earlier and deeper: a Saturn-mass planet at a given age rains more than a Jupiter-mass one. The phase diagram is the weakest input: [R1895] say it is "very difficult to determine the critical line of a mixture at megabar pressures and many thousands of degrees".

### 2.3 The core: compact, dilute, or both

| Quantity | Value | Source, state |
|---|---|---|
| The old picture (pure heavy-element core under one homogeneous envelope) | "no longer valid"; gravity data need a deep region enriched in heavy elements with a gradient, a "fuzzy" or "dilute" core | [R1895] sections 1 and 2: READ |
| Jupiter's dilute core | mass about 100 Earth masses (a third of Jupiter), H-He-dominated, heavy-element fraction about 20 percent; small-core solutions extend to about 20 percent of Jupiter's mass | [R1895] section 2.4 and Figure 2: READ; [R2013] abstract: READ (abstract) |
| Jupiter, other models | dilute core out to 60 to 75 percent of the radius; region between 0.1 and 2 Mbar with an entropy and composition change; outer envelope near solar Z | [R1038] abstract and Figures 3 to 4: READ |
| Jupiter's total heavy elements | 18 to 45 Earth masses (adiabatic about 20, non-adiabatic up to 60) | [R1893] section 3.1.1; [R1895] section 2.1: READ |
| Saturn's dilute core | stable region out to 0.59 +/- 0.01 R_S, 55.1 +/- 1.7 Earth masses of which 17.4 +/- 1.2 are rock and ice (Z about 0.32 there); central density about 6 g/cm3; up to 80 percent heavy elements in the inner part | [R2851] text: READ; [R1893] section 3.1.2: READ |
| Saturn, compact-core models | pure core of 15 to 18 Earth masses with under 5 Earth masses in the envelope; 4-layer fit 12 to 18 Earth masses total | Ni 2020, via [R1893]: 2ND-HAND |
| A compact pure heavy-element core | "a few Earth masses at most" for Jupiter; may not exist; invisible to gravity (a 5 Earth-mass core adds about 1e-4 of the moment of inertia); formation models predict it, and convective mixing or a giant impact could erase it | [R1895] section 2.3: READ |
| Why the core is dilute | four candidate origins: evaporated accreted solids mixing in (confined to about the inner 20 Earth masses by runaway gas accretion), giant impacts, convective erosion, and delayed runaway accretion at about 100 Earth masses | [R1895] section 3: READ |
| What it is made of | water and rock mix with metallic hydrogen in all proportions at these temperatures; gravity cannot tell water from rock; rocks may sink deeper than water during accretion, so a compact core is rock-and-iron rich and the ice mixes higher up | [R1895] section 4 ("Differentiation", "The chemical composition of heavy elements"): READ; 2ND-HAND |
| Temperature | central temperature about 20,000 K adiabatic, 35,000 K for superadiabatic (Jupiter); 10 Mbar region 9,500 to 13,000 K; boundary 1 to 2 Mbar at 5,000 to 9,000 K; heavy-element mass and central temperature trade off | [R1893] section 3.1.1; [R1038: Figure 4]: READ |

Conflict (rule 3): compact core 0 to a few Earth masses [R1895] against 15 to 18 (Saturn, 2ND-HAND). **The compact core is drawn between 0 and 20 Earth masses** (the older decision 22 range of 5 to 20 reaches the Saturn value; the lower end now includes none), D3.

### 2.4 Ice giants

| Quantity | Value | Source, state |
|---|---|---|
| H/He share of the mass | about 10 to 20 percent [R1892]; up to about 2 Earth masses in three-layer models | [R1893] section 3.2: READ |
| Heavy-element share | Uranus 76 to 90 percent, Neptune 77 to 90 percent, depending on whether the heavy element is rock (SiO2) or water | [R1889, R1890, R1891] via [R1893]: READ |
| Structure | three layers (rock core, water-rich mantle, H/He envelope) or gradients with no distinct layers; outer envelope Z up to 0.17 (Voyager rotation) or 0.08 (revised), inner envelope Z about 0.9; solar-Z envelopes also fit | [R1893]: 2ND-HAND |
| Rock to water ratio | unknown: at most about 2; 0.05 to 0.38 maximum for non-adiabatic empirical models; solar 1/2; carbon-rich planetesimals give more methane and water than rock ([R2827]: methane planets) | [R1893] sections 3.2.1 and 3.2.2: READ; each paper 2ND-HAND |
| Is the hydrogen metallic | the adiabats cross the molecular-to-metallic line, but "it is unclear whether H-He exist in the planetary deep interiors" | [R1893: Figure 1] caption: READ. Estimate here: the pressure at the base of a 2 Earth-mass envelope on a 13 Earth-mass core is about 21 GPa (0.2 Mbar), below P_MM, so no metallic layer (a thin-shell estimate, `G M_core M_env / (4 pi R_core^4)`, COMPUTED) |
| Water in the deep interior | superionic ice (ice XVIII) at 100 to 400 GPa and 2000 to 3000 K; solid oxygen lattice, mobile protons; the generator needs its density and that it is solid | C-volatile-and-giants.md section 2.3 ([R1457] via Zeng and Sasselov): READ (abstract), 2ND-HAND |
| Low luminosity | composition gradients and thermal boundary layers; central temperatures 5,000 to 6,500 K adiabatic, above 20,000 K non-adiabatic | [R1893] section 3.2.1: READ |

### 2.5 Brown dwarfs

A brown dwarf is fully convective and partially degenerate ([R0707], astro-ph/0006383, sections 2 and 3: READ for the degeneracy parameter psi running from about 2 to 0.05 and the central density rising as mass falls toward the hydrogen-burning limit). Its central pressure is hundreds to thousands of times P_MM, so only a thin outer layer is molecular; it forms from a collapsing cloud with no solid core, so its heavy-element mass fraction is the star's metallicity Z with no compact core (D6). Whether helium rain occurs inside cool T and Y dwarfs is OPEN (no source was read); the same rule as section 2.2 applies if it does, with T1 replaced by the effective temperature, and is **not** switched on for brown dwarfs until a source is read (GEN.641 follow-up).

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

The hydrostatic profile of the layers is the n = 1 polytrope that the repo already uses for the central pressure: P(xi)/P_c = (sin xi / xi)^2 with xi = pi r / R, and the mass inside xi is M (sin xi - xi cos xi)/pi, P_c = pi G M^2 / (8 R^4). It is a rough shape for a hydrogen-helium sphere, good enough to place layer boundaries to a few percent of the radius (COMPUTED here: for Jupiter's mass and radius the 2 Mbar level falls at 0.81 R_J and 16 percent of the mass lies above it, against 0.8 to 0.9 R_J in Debras and Chabrier; for Saturn 0.59 R_S against 0.55 in Mankovich and Fortney). It is the stand-in of a run without the downloaded tables; with them the profile is the structure of section 11. For an envelope under 30 percent of the mass (ice giants, sub-Saturns) the shell formula `P_base = G M_core M_env / (4 pi R_core^4)` replaces it, because the mass is in the core. The radius R is the one GEN.429 and GEN.460 derive, so the layers never change mass or radius.

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
- **Dilute core mass M_dil**: log-uniform from 20 to 100 Earth masses (the delayed-runaway mass of Venturini and Helled at the high end, both 2ND-HAND through Helled and Stevenson), capped at 0.65 M (Saturn's 55 of 95 Earth masses is 0.58, READ). Below M = 30 Earth masses there is no separate dilute core and the heavy elements are a gradient through the whole body (the ice-giant picture).
- **Heavy-element budget.** M_Z (GEN.460) is split: Z_out M_mol in the outer envelope (Z_out = Z_env from decision 21), M_cc in the compact core, and the rest in the dilute core and the metallic envelope. The dilute core's mean Z_dil is that rest over M_dil and must fall between Z_out and 0.8 (Saturn's inner part, READ); if it does not, M_dil is moved inside its bracket until it does, and when no value fits the planet has no dilute core and the heavy elements are spread through the metallic envelope. Z rises inward through the dilute core; the shape of the gradient is OPEN (a linear rise from Z_out at its top is the placeholder) and only the mean, which the budget fixes, is a constraint.
- **What the heavy elements are made of.** The rock-to-ice split of the heavy elements comes from the formation branch GEN.428 stored: a body born beyond the snow line takes the water share W of its solids (decision 2, 0.14 to 0.39 of the solids) and the ice set of decision 18; a body born inside takes none. The compact core takes the rock and iron first (rocks sink deeper than water, Helled and Stevenson section 4), and the ice goes outward first: the molecular envelope (which holds it as vapour), then the metallic envelope, then the dilute core. Rock and iron that the compact core cannot hold fill the dilute core, then the deep metallic envelope. Heavy elements the dilute core cannot hold (it takes at most 80 percent of its mass) enrich the deep envelope, never the atmosphere, whose metallicity the caller fixed. The iron-to-rock split of the refractory part is the star's, as for a terrestrial planet (GEN.428's `metal_mass_fraction` over `rock_mass_fraction`, CMF 0.2917 for solar abundances).
- **Temperature** is not drawn or stored (OPEN): the boundary temperatures of Debras and Chabrier (5,000 to 9,000 K at 1 to 2 Mbar, 9,500 to 13,000 K at 3 to 7 Mbar, 19,000 to 30,000 K at 40 Mbar) are test targets, not inputs.

### 4.4 Ice giants

- f_env is GEN.428's; the H/He share is 0.05 to 0.20 of the mass (Helled and Fortney 10 to 20 percent for Uranus and Neptune; widest range with the three-layer fits that give 1 to 2 Earth masses).
- **Rock fraction of the heavy elements** 0.05 to 0.67 ([R4880] high; rule 3 widest, OPEN). Water, ammonia and methane share the rest by the ice set (GEN.428), which already follows Marboeuf's table; the carbon-rich option of [R2827] (methane above water) is a table row with weight set by the share of refractory organics, OPEN.
- Mantle water is `water_liquid` above a drawn depth and `ice_x` (superionic) when the layer's pressure exceeds 100 GPa; the radius does not change (the radius tables already carry it, decision 12).
- A gradient (no distinct layers) is the default representation of Z; the three layers are the rows a description names.

### 4.5 Brown dwarfs

One envelope: molecular above P_MM, metallic below, Z the star's metallicity, no compact core, no dilute core (D6). The molecular layer is thin, because the central pressure is far above P_MM; its mass share is computed when GEN.642 is built (an n = 3/2 polytrope for the degenerate range), not in this note. Helium is the proto-solar Y of the natal gas (stellar-composition.md) until a source for cool-dwarf helium rain is read.

## 5. What gets stored and described

- **Summary columns** on the planet and moon tables (DB.89, lane I's BulkFractions extends): `core_mass_earth` (compact core, 0 if none), `dilute_core_mass_earth`, `metallization_pressure_gpa` (P_MM), `helium_depletion` (Y_out / Y_proto, 1 when no rain), `has_helium_shell` (bool). All numeric columns are declared in `db/column_units.py` in the same PR.
- **Layer rows** (built, schema v154): `planet_interior_layers` and `moon_interior_layers`, one row per layer outside in, with `layer_index`, `kind`, `state`, `top_pressure_gpa`, `bottom_pressure_gpa`, `top_radius_fraction`, `bottom_radius_fraction`, `mass_earth` and the mass of each material in Earth masses (`hydrogen`, `helium`, `iron`, `rock`, `water_liquid`, `ice_x`, `nh3_ice`, `ch4_ice`). The tables are separate from the mixture child table of GEN.428, which was not extended. The body carries them as `interior_layers` (a list of flat dicts, `giant_interior.body_values`) and `giant_interior.interior_of` rebuilds the `GiantInterior`. A recalculation never redraws a stored layer (F10).
- **Description** (GEN.642, one writer with GEN.430; the numbers are an illustration of the wording, not a draw): "A gas giant of 1.0 Jupiter masses. Its outer 6 percent is molecular hydrogen with helium at 85 percent of the natal share; below 0.9 Mbar the hydrogen is metallic. About 40 Earth masses of heavier elements are spread through a dilute core that reaches 25 percent of the mass; the heavy elements there are rock and iron with 5 percent water ice. No compact core." Thresholds by fixed rules, as GEN.430.

## 6. Test targets (never presets)

Each target is met by some draw inside the allowed distributions, not by one draw, and lives in test files only.

| Target | Tolerance | Source |
|---|---|---|
| Jupiter-class mass and age: P_MM level at 0.80 to 0.93 R, dilute core 15 to 35 percent of the mass, Z_out under 0.1, total heavy elements 18 to 60 Earth masses | inside the bracket | [R1038, R1895] |
| Saturn-class (95 Earth masses, T1 near 140 K): Y_out / Y_proto 0.15 to 0.4, a helium shell present in some draws, dilute core out to 0.4 to 0.65 R | inside the bracket | [R2851, R2852] |
| Uranus and Neptune class: H/He 0.05 to 0.20, no metallic layer, Z 0.75 to 0.95, rock fraction of Z 0.05 to 0.67 | inside the bracket | [R1892] |
| Layer masses sum to 1 and layer pressures are monotone; layer heavy-element masses sum to M_Z (1e-9) | exact | invariants |
| Helium is conserved across the layers (1e-9) | exact | invariants |
| A giant with envelope under one half has no `metallic_envelope` unless its base pressure passes P_MM | exact | section 3 |
| Same seed, same layers; a recalculation changes none (F10) | exact | bulk-composition.md section 8 |

## 7. Decisions

| # | Decision | State | Used meanwhile, source |
|---|---|---|---|
| D1 | The molecular-to-metallic pressure is drawn per body, log-uniform 0.4 to 2 Mbar; the fixed 1 Mbar goes | Default taken | [R2966] (0.4 to 0.6), [R1038] (1 to 2), rule 3 widest |
| D2 | Helium rain by the one-bar temperature: a line through two measured points, a helium shell at strong depletion | STAND-IN, flagged in the gap plan registry | [R2852] give the mechanism and the 149 K unshifted onset; the line and the shell rule are placeholders until the Schottler and Redmer curve and the adiabat replace them |
| D3 | Compact core 0 to 20 Earth masses, probability 0.5 of none, replaces decision 22's 5 to 20 | Open | [R1895] (a few Earth masses, may not exist) |
| D4 | A dilute core of 20 to 100 Earth masses, at most 0.65 of the mass, gradient in Z | Default taken | [R1038, R1895, R2851] |
| D5 | Rock and iron to the compact core first, then the dilute core; ice outward first (molecular envelope, metallic envelope, dilute core); split from GEN.428's formation branch | Default taken | [R1895] section 4 |
| D6 | A brown dwarf has no compact or dilute core; metallic envelope below P_MM; no helium rain until a source is read | Default taken, rain open | [R0707] |
| D7 | The layer pressures come from the n = 1 polytrope (shell formula under 30 percent envelope) until GEN.460's tables supply real profiles; with the hydrogen-helium download the structure of section 11 replaces it (GEN.733), and the polytrope is the fallback of a run that disregards missing data | STAND-IN | `constants.py` line 195; COMPUTED checks in 4.1 |
| D8 | Layer rows live in the composition child table with a layer index; summary columns extend BulkFractions | Default taken | contracts file of lane I, section 1.1 |
| D9 | Whole-body accounting: heavy elements in envelope layers count as `envelope_z` of the envelope; heavy elements in core layers count as metal, rock, water and ice; all hydrogen and helium (the dilute core's too) is envelope; the helium share is Y/(X+Y) of the hydrogen and helium | Default taken | Lane I's `BulkFractions` contract |

## 8. Gaps

1. **The [R3973] phase curve** are not in hand; their data (Figure 1 of Mankovich and Fortney is a plot of them) would replace D2's line. To fetch.
2. **Equation-of-state tables** for hydrogen-helium: the [R0708] tables are downloaded since GEN.733 (section 11); the polytrope stays as the fallback of a run that disregards missing data (D7).
3. **Cool brown dwarf helium rain**, neon rain, and the helium content of Y dwarfs: no source read.
4. **The probability of a compact core** and the gradient shape in the dilute core: no source; the placeholders are flagged.
5. **[R2966] crossover** is at 6000 K; the pressure at 3000 to 5000 K has no number in the read papers, so the bracket 0.4 to 2 Mbar is temperature-independent.
6. **Exoplanet giants**: inflated hot Jupiters may have a radiative layer or a different entropy; none of the sources read covers them (GEN.460's inflation row stays the radius authority).
7. **Hot and warm Jupiters' T1 for the rain test**: the equilibrium temperature with internal heat, not the Solar System's 166 K; the stored `temperature` is a surface value and the one-bar value needs the greenhouse and internal-heat terms of the atmosphere lane.

## 9. Sources fetched for this note

[R1893] ([R1893], full), 2403.11657 ([R1895], full), 1912.01009 ([R2852], abstract, Figures 1 and 3, sections 2 to 4), 1901.05697 ([R1038], abstract and structure figures), 1709.08648 ([R2966], abstract and metallization text), 2006.12219 ([R1894], sections 3 and 4), 2104.13385 ([R2851], abstract and mass numbers), 2302.09082 ([R2013], abstract), 1010.5546 ([R1889, R1891], summary), astro-ph/0006383 ([R0707], degeneracy sections). Raw papers stay out of the repository (Boss, 2026-10-17).

## 10. Wired into the generator (GEN.460, GEN.642; generator epoch 80)

`physics/giant_draw.py::draw_giant` runs once a giant's mass is settled (`planets.apply_giant`). The order is part of the epoch:
the envelope's metallicity (a multiple of the star's Z, log-uniform, `ENVELOPE_Z`, stand-in), the heavy-element mass M_Z (class J:
[R0710] with its scatter; classes I and T: a uniform share of the mass, `HEAVY_FRACTION`, stand-in), the radius (the
Fortney tables where `giant_radius.in_range` and `covers` hold, plus the Sestovic inflation, kept in the class's radius range; else
the anchored relation of GEN.34), then the interior (`draw_interior`) whose whole-body totals become the stored bulk composition.

Stand-ins, each named in `giant_draw`: the envelope metallicity range, the I and T heavy-element shares, the water share of the
heavy elements beyond the snow line and the ammonia and methane proportions to it, the single Bond albedo (Li et al. 2018's 0.503),
the own heat flux (the cooling law through Jupiter and Neptune of `rogue_surface.giant_internal_flux_w_m2`, `giant_draw.own_flux_w_m2`, GEN.816: L ~ M^1.65 below a Jupiter mass, t^-1.3 in age; it has no helium-rain share until GEN.810, so a Saturn-like giant comes out near 115 K against 134 K) and the one-bar over effective temperature (1.33, from Jupiter's 166 K over
125 K). A giant moon (the host's distance is not carried) takes the anchored relation and a 130 K one-bar temperature. The interior
is stored (DB.89, `giant_interior.body_values`); the old stand-in envelope and core shares (`bulk_draw.draw_giant`) are removed.

**Settled once the age and the orbit are final (GEN.961, GEN.976).** The draw above reads the star's age and the starlight of the circle at the
semi-major axis, and the star's age rises when the planets are in (`Star.adjust_age_for_planets`; the eccentricity is drawn after the
radius, which it needs). Left alone, 28 percent of the table-sized giants of 138 measured had a stored radius over 1 percent too large and
the worst 9.3 percent. After the age is final and the placed bodies are trimmed, the "giant settle pass" of `StarSystem`
(`planets.settle_giant`, `giant_draw.settle`) takes the table radius again at the final age and the mean light of the orbit
(`orbit_elements.flux_distance_au`, no draw), and a giant that it moves by more than `SETTLE_TOLERANCE` (0.1 percent) has its volume,
density, gravity, interior (layers, bulk composition and text, one-bar temperature at the mean light) and air made again. The star's ratios,
the envelope's Z and the heavy-element mass are the draws it keeps (`DrawnGiant`); the interior and the air are drawn from streams of the giant's own
(`giant_draw.settle_stream`, keyed by its mass, distance, star and the age and radius it settles from), so the run's stream does not move and
a system whose planets the age did not move is bit-identical. Giants the tables do not size (ice giants, giants under 17 Earth masses, a
given radius, a moon) are not touched. A moon that the smaller giant can no longer hold (larger than `planet.radius / 10**(1/3)`, or outside the range of orbits the new radius and scale height leave) is dropped (`drop_unstable_moons`). A moon's tide does not read its parent's radius (`tidal_heating.forcing_of`: the host's mass, the moon's distance and eccentricity; GEN.987, tested), so the moons' heat needs nothing more. The tides of every planet's and captured moon's orbit are taken again at the final age in the same pass (`planets.settle_orbit`, GEN.985; the birth eccentricity is kept with the body, the position and velocity are made again on the new ellipse, a planet's moons are told). Not settled, and filed: the air, class, dose and spin of a rocky planet or moon (GEN.992 and the notes of GEN.955; `reports/gen961-stale-age-audit-2026-10-30.md`, `reports/gen986-final-age-rocky-planets-2026-10-30.md`).

## 11. The structure from the equation of state (GEN.733)

The pressure at a layer boundary is read from a hydrostatic profile (section 4.1), which was the n = 1 polytrope for a gas giant, the
n = 3/2 polytrope for a brown dwarf and three layers of constant density for an ice giant (D7). A polytrope has one index for the
whole body where the real fluid goes from an ideal molecular gas to a degenerate metal, so the profile now comes from the equations
of hydrostatic balance, `dP/dm = -G m / (4 pi r^4)` and `dr/dm = 1 / (4 pi r^2 rho)`, integrated inward from the one-bar level in ln P
(`physics/giant_structure.py`).

**The density.** The hydrogen and helium of the envelope come from the (T, P) tables of [R0708] (DirEOS2021,
registry fill GS1; `datasets/hhe_eos.py`, downloaded by `update.sh`): the mixture of helium share 0.275 with the helium-hydrogen
interaction, moved to the body's helium share by the ideal volume law (GS3). The temperature at each pressure is on an adiabat, a
line of the table's entropy (the tables' adiabatic gradient carries numerical noise, GS6). The heavy elements are the repository's
own curves (iron's Vinet fit, the silicate mantle table, the water isotherm of `bulk_structure`, GS4), as a core, an ice mantle and
a share of the envelope mixed in by the ideal volume law (GS2).

**What fixes the entropy.** The radius is GEN.460's and is kept. The generator's one-bar temperature comes from the starlight, and
an irradiated or inflated giant is not on the adiabat that temperature starts, so the model takes the one-bar temperature of the
adiabat that closes it on the radius (`fit_adiabat`, about ten integrations, 10 to 40 ms a giant) and reports it (`Profile.adiabat_1bar_k`).
The generator's own one-bar temperature is only the fit's starting point. A radius no adiabat in 30 to 8000 K closes (a body too
light or too dense for it) takes the nearest adiabat, says it did not close (`Profile.closed`) and the layers use the fractions of
the radius it did close at.

**Checks against the sky** (the real tables, `tests/test_datasets_real.py`): Uranus and Neptune close on one-bar adiabats of 50 to
100 K for the heavy-element masses and ice shares of the interior models (measured 76 and 72 K); Jupiter and Saturn close with their
radii; a 40 Jupiter-mass brown dwarf reaches the degenerate regime at the centre. Where the tables matter most is a Saturn-like body,
whose envelope is far from n = 1: the mass above 100 GPa is 31 percent of the body against the polytrope's 37 percent and the
central pressure is ten times higher.

**Without the download.** A run told to disregard missing data (GEN.599) takes the profiles of section 4.1 and records it
(`fetch.stand_ins_used`); any other run is refused (GEN.590), like the other datasets.

**Limits** (registry 6.134): the heavy elements are a core and a uniformly mixed envelope with no dilute core or helium gradient
inside the structure (GS2); one natal helium share for the body (GS3); ammonia and methane ice take water's density (GS4); the
adiabat of a body under 100 K is carried up the molecular gas's own (GS5); the tables' unphysical regions are read as they stand
(GS6); the centre is closed as a sphere (GS7).

## 12. The internal heat (GEN.810)

A giant's own heat (`physics/giant_heat.py`, registry 6.135) is stored as its `internal_heat_flux_w_m2`, the heat it radiates per square metre of
its 1 bar surface, and has two parts. The **base** is the cooling law of `rogue_surface.giant_internal_flux_w_m2` (GEN.816: luminosity as the mass to the 1.65 below a Jupiter mass, the age to
the -1.3), the same law the one-bar temperature in `giant_draw` is made from before the interior exists. The **helium rain** is read back from what the interior stored: every layer's
helium against the natal share of its hydrogen and helium says which layers lost helium (the envelope above) and which gained it (the
helium-enriched zone, the shell, the dilute core); the potential at each layer's half-mass radius follows from the layers' stored radii
and masses (each a shell of uniform density, exactly the potential of a uniform sphere for one layer); the energy is the mass that sank
times the potential it fell through; half of it leaves as heat (GH2), over the star's age (GH3). Nothing is added when the helium did not
rain (`helium_depletion` is 1), and a giant whose helium went to a shell radiates more than one whose helium stayed in a gradient, because
the shell is deeper.

A Saturn-like giant (95 Earth masses, 140 K at one bar) lands at about 2.4 W/m2 against the 1.88 W/m2 measured (the base 1.1 and the rain 1.3: with
the share at 0.3, the low end of GH2's knob, it is 1.85), a Jupiter-like one at 5.7 against 5.4, and a Neptune-like one at 0.43 against 0.43; none is fitted. The rain's heat is not fed back to the one-bar temperature that decided the rain (GH4). An ice giant
radiates the base only, which is what the law through Neptune gives it.
