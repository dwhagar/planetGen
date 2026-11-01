# Surface pressure supply law: how much gas a rocky body has delivered to its surface (GEN.726)

Status: design, 2026-10-24, lane I (Blocker I). Boss, 3:50 PM EDT 2026-10-24: "Hold 4c until the supply law is done"; earlier, 3:30 PM EDT: "I want to be as accurate to known physics as possible". This note defines `P_supply`, the one term of the retention design `P = min(P_supply, P_retention, P_vapour)` (`habitability/atmospheres-retention-and-classes.md`) that had no definition, names the source of every term and the tier of every number that is not read from a source. GEN.411 and GEN.422 (the wiring of the mixture draw and the order of the pressure sum) wait for it.

## 1. What is being defined

The pressure of a hydrostatic atmosphere is its mass over the area, `P = g M_air / (4 pi R^2)`. The supply law gives the mass `M_air` that the body's interior has delivered to the surface by its age and that no surface sink has taken back:

    M_air = sum over gases of  (atoms of the gas's elements degassed) - (what the sinks hold)

It does not subtract escape. The loss of a body that is too small or too hot to keep what it has is the retention row's job (`d1` to `d5`, the cosmic shoreline of [R4862], GEN.511) and the Jeans filter that already removes the light gases from the mixture (`mixture_escape`, GEN.693); subtracting it here as well would count it twice. The vapour cap is the third term of the minimum. So the law answers: if nothing escaped, how much gas would be here?

What it replaces: on main the pressure is a class density range times `gravity^1` ("a starting point, not derived from a retention model"). GEN.422 deletes the class `atm_density` and `atm_molar_density` ranges; this law is where the number comes from instead.

## 2. The terms

### 2.1 The element inventory (read)

`mixture_rocky.rocky_inventory` already gives the delivered budget of the outgassing stage: H, C, N, S and radiogenic Ar, in atoms per 100 g of the mantle that holds them. Water and carbon dioxide are [R3378: Table 1] (water 0.05 to 0.5 weight percent with CO2/H2O 0.2, the dry class CO2 0.018 to 0.1 weight percent (Ortenzi's stored 0.018 to the wet High case; his delivered 0.6 is 16 times the bulk silicate Earth's carbon and is not used, GEN.813), the water-rich class 10 to 50 weight percent); nitrogen and sulfur per carbon atom are brackets that contain the bulk silicate Earth and the CI chondrites; argon is made by the rock's potassium. These are kept as they are (tier 1 and 2 as the gap plan registry rows AM13 and AM14 say). What was missing is the mass they refer to and the share that has come out.

**The mass is the mantle's**, `M_mantle = M * rock_mass_fraction` (the silicate share of `bulk_composition`; the core and the ice do not hold the delivered volatiles). The inventory's per-100-g basis was a scale-free device for the mixture, and the Ortenzi values are mantle contents, so the code reads them per 100 g of mantle.

[R2923] (EPSL 313-314, 56-66; [R2923], read) puts the water of the bulk Earth at 1000 to 3000 ppm and the carbon at up to 500 ppm, "both largely sequestrated in the solid Earth", and says up to 90 percent of the water and up to 97 percent of the carbon may be in the mantle. Ortenzi's 500 to 5000 ppm water and 100 to 1000 ppm CO2 are the same order of magnitude and are the part a planet can outgas, so the inventory above is an upper bound on what reaches the surface.

### 2.2 The share degassed by the body's age (the part that is calibrated)

The share of an element that has reached the surface is

    phi_e = eps_e * (1 - exp(-Theta)),     Theta = E_release / (M_mantle * e_proc)

`E_release` is the heat the interior has carried out since the body formed (J), a new output of the interior heat model (`interior_heat.Heat.released_j`, the integral of the heat flow over the age). The physical picture is [R4413] ("The habitability of a stagnant-lid Earth", [R4413], read, sections 2.2 and 2.3): melt is produced where the mantle is hotter than its solidus, the melt rises and erupts, and the gas the melt cannot hold dissolved comes out. The melt production is driven by the same convective heat loss, so the mantle processed through melting grows with the heat carried out. `e_proc` (J per kilogram of mantle) is the heat released per kilogram of mantle for one e-folding of the undegassed share. The form `1 - exp(-Theta)` is the first-order depletion of a well-mixed reservoir, the standard form of mantle degassing models.

Alternatives considered and left out: integrating the melt production of [R4413] equations 15 to 19 (it needs the solidus, the melt fraction and the crust growth, none of which the interior model carries), and the modern outgassing flux of the carbonate-silicate models (a flux today, no history). The heat release is the history the project already computes and stores the inputs of.

`eps_e`, the share of the element the melt cannot hold:

| element | eps | basis |
|---|---|---|
| C, N, Ar | 1 | CO2 and N2 are far less soluble than water in basaltic melt, "the low solubility of CO2 causes most of the carbon to be outgassed" ([R4413], abstract); argon is a noble gas |
| S | 0.5 | sulfur stays partly in the melt as sulfide (stand-in, tier 6) |
| H | drawn 0.05 to 0.5, log | water is highly soluble and "limited to a few tens of bars" [R4413]; the magma ocean retains much of it ([R0501]: "much of the water delivered during accretion can be retained in the interior"); [R2923]: up to 90 percent in the mantle. The bracket's top is Marty's 10 percent surface share doubled for the water that never was in the mantle (stand-in, tier 6) |

**Calibration of `e_proc`.** The Earth and Venus both hold about 1e22 mol of carbon at the surface: Earth 7 to 11e21 mol in carbonates and organic matter (a review's figure met in a search, not yet checked against the paper: tier 3) and Venus 1.05e22 mol of CO2 in its air (9.2 MPa over 4 pi R^2 / g = 4.77e20 kg of air, 96.5 percent CO2). With the log-middle of the Ortenzi CO2 range, 0.032 weight percent of a 4.0e24 kg mantle (1.3e21 kg), Earth's 4.6e20 kg of CO2 is a share of 0.36, so `Theta_Earth = 0.45`. The interior model gives the Earth twin its heat release, and `e_proc` is set so that the twin reaches 0.45. Venus is then a check, not a fit: its heat release per kilogram of mantle (stagnant lid in the model) should give a share inside a factor of 1.5 of the 0.43 its air needs (4.6e20 kg of CO2 over 3.3e24 kg of mantle at the same inventory).

### 2.3 The sinks

**The ocean (water, sulfur).** When the surface holds liquid water, the vapour in the air is the saturation pressure of water at the surface temperature times the relative humidity the habitability code already uses; the rest of the water is the ocean (`hydrosphere`). So the air's water is `min(supply, saturation)`: the vapour term of the minimum. Sulfur gases are taken up by an ocean (sulfate and sulfide; `SULFUR_SINK` of `mixture_evolved`, already a registry row AM2).

**The carbonate-silicate cycle (carbon dioxide).** On a body with an ocean and rock to weather, silicate weathering removes CO2 at a rate that grows with temperature and with the CO2 pressure, and volcanic outgassing puts it back; over geological time the two balance ([R2628, R4596], "Carbonate-silicate cycle predictions of Earth-like planetary climates and testing the habitable zone concept", PMC7708846, read). The continental weathering law in the form of that paper (Eq. 4 of its simplified model; Table 1 for the ranges):

    F_w = rho * f_land * (pCO2 / pCO2_mod)^alpha * exp((Ts - Ts_mod) / T_e)

with `pCO2_mod = 288e-6 bar` and `Ts_mod = 289 K` (preindustrial Earth), `alpha` 0.1 to 0.5 (geologic constraints over the last 100 Myr), `T_e` 10 to 40 K (laboratory and Phanerozoic constraints), `f_land` the land fraction (Earth 1 in their scaling, so here `f_land / f_land_Earth`), and outgassing `F_out = F_out_mod (q / q_Earth)^m` with `m` 1 to 2 and `F_out_mod` 6 to 10 Tmol C per year. In balance `F_out = F_w`, so

    pCO2* = pCO2_mod * [ (q/q_E)^m * (F_out_mod / rho_mod) / (f_land / f_land_E) ]^(1/alpha) * exp(-(Ts - Ts_mod) / (alpha T_e))

with `rho_mod f_land_E = F_out_mod` for the modern Earth (a steady state at 288 ppm by construction). The weathering and the outgassing both scale with the surface area, so it cancels. `Ts` depends on `pCO2` through the greenhouse, so the balance is solved with the temperature function of the body (`atmosphere_state.greenhouse_warming`, the fixed point of GEN.422). The CO2 left in the air is `min(supply, pCO2*)`: a world with little CO2 has all of it taken; one with a lot is held at the balance and the rest is in carbonates. The unstable branch [R1661] report in the outer habitable zone (weathering falling with pCO2) is not modelled: the solution takes the stable root.

The cycle needs liquid water at the surface (ocean fraction above zero, `Ts` between the freezing point and the critical point of water) and rock to weather. Without liquid water (Venus, Mars today) nothing takes the CO2 back and the air holds all the supply, which is what makes Venus a 90 bar CO2 world with the same carbon as the Earth. A frozen world keeps its CO2 until the greenhouse thaws it (the balance's temperature factor lets `pCO2*` climb, capped by the vapour term when the CO2 would condense).

Nitrogen and the noble gases have no surface sink (N2 is almost insoluble); the Earth's mantle holds more nitrogen than its air ([R2923]: 0.84 +/- 0.43 ppm of the silicate Earth in the mantle against 3.9e18 kg in the air), which the share `phi_N` already carries.

### 2.4 From atoms to a pressure

The mixture draw takes the element budget of the gas and returns mole fractions in equilibrium with the mantle. The supply law changes only the budget it is given: the inventory times `phi_e` (what has been degassed), then the sinks take atoms out (ocean: H and S; carbonate: C and the O that went with it). The moles of gas follow from the atoms and the atoms per mole of the mixture; the mass from the molar mass; the pressure from the mass. In steps, in the order GEN.422 fixes: mixture, pressure with freeze-out, density, scale height, greenhouse, temperature (the fixed point), with the sinks inside the loop because the carbonate balance needs the temperature.

## 3. Anchors the tests keep

| body | what is checked | value | basis |
|---|---|---|---|
| Earth twin | `Theta` from the interior model gives the calibrated share | 0.45 (set) | surface carbon 1e22 mol |
| Earth twin | `pCO2*` at the modern flux and 289 K | 288e-6 bar (identity) | [R2628] |
| Venus twin | degassed CO2 inside a factor of 1.5 of the air | 4.6e20 kg | 9.2 MPa, 96.5 percent CO2 |
| Venus twin | no ocean: no sink, the air is the supply | | |
| Mars twin | the supply is above the observed 0.6 kPa: the observed air is the retention row's loss, not the supply's | 2.5e16 kg observed | non-thermal escape (not modelled here) |
| Titan | out of scope: the `titan` and `frost` families set their own pressure | | |

## 4. Registry rows this adds (source tier, tweak level; GEN.725 rule)

| name | value | source | tier | tweak |
|---|---|---|---|---|
| `PROCESSING_ENERGY_J_PER_KG` | set by the Earth calibration | Earth and Venus surface carbon, [R4413] picture | 5 confident educated guess based on physics | 1 nudged within the range of two anchors |
| `DEGAS_EFFICIENCY_H` | 0.05 to 0.5, log | [R0501, R2923, R4413] | 6 marginal educated guess | 1 |
| `DEGAS_EFFICIENCY_S` | 0.5 | sulfide retention in melt | 6 | 0 |
| `WEATHERING_ALPHA` | 0.1 to 0.5, uniform | [R2628: Table 1] | 1 verified (range) | 0 |
| `WEATHERING_T_E_K` | 10 to 40 K, uniform | same | 1 | 0 |
| `OUTGASSING_EXPONENT` | 1 to 2, uniform | same | 1 | 0 |
| `LAND_FRACTION_REFERENCE` | 0.29 (Earth) | the Earth's land share | 1 | 0 |
| `CO2_MODERN_BAR`, `TS_MODERN_K` | 288e-6, 289 | same | 1 | 0 |

The draws go through `util/draw.py` after the existing mixture draws, in the order above, so a change of the order is a generator-epoch change (and the whole law is one).

## 5. What it does not do

- The magma-ocean epoch: the inventory is taken as already split between mantle and surface at the start of the geological stage ([R0501] treat the first million years; the water share `eps_H` stands for it).
- Impact erosion and delivery after the main accretion (the physics gap register's A5), and non-thermal escape (the retention rows).
- Seafloor weathering, organic carbon burial, the biological weathering factor `f_bio` (taken as 1) and carbonate metamorphism of a stagnant lid [R1403]: the cycle is one box, the stable root.
- A tectonic or stagnant-lid switch on the outgassing: the heat flow `q` carries the regime only through its value, because the regime is drawn after the hydrosphere. This is the largest simplification.
- Giant and ice-world envelopes, Titan-like and frost air: their pressures come from their own families.

## 6. In the generator (GEN.411, GEN.422)

Built 2026-10-24 (lane I, Blocker I): `physics/surface_air.py` solves the order of section 2.4 for every rocky planet and moon and every giant, and `physics/planets.calculate_atmospheric_conditions` stores the result (`atmosphere.store_air`). The class density ranges `atm_density_range` and `atm_molar_density_range`, the greenhouse multiplier ranges, `CLASS_MIXES` and the old mixing-ratio draw are gone (GEN.683); no body is given an air by its class.

**One body, one fixed point.** Draws first (design constants of the law: the efficiencies, the cycle's three, the hydrate skin; the seed of the replay), then repeated steps at a temperature and a pressure:

1. The recipe subject of the body at (`Ts`, `P`, ocean, geology) (`mixture_subject.subject_of_body`, GEN.756) and the mixture draw over the inventory times the degassed share, less the carbonate sink at that temperature (`surface_pressure.carbonate_pco2_bar`). The draw replays from one seed at every step, so only the temperature and pressure move it. The regimes of the recipe depend on them, so the draw is kept from step 8 on (`FREEZE_AFTER`); a regime on the edge of its gate cannot flip the solution back and forth.
2. The moles of the mixture (`MixtureDraw.moles`, per 100 g of mantle, times the mantle mass) and the light-gas escape (`mixture_escape`) give the column.
3. The freeze-out (`freeze_out`): a gas above its saturation pressure condenses until it stands on the curve, and water over an ocean or ice stands at the relative humidity 0.77 ([R2848], the Earth's mean) while the water lasts; the quadratic of each gas is solved exactly, the mean molar mass changes with each, and the result is repeated until none moves. An ocean dissolves the sulfur gases by Henry's law (`mixture_sinks.sulfur_remaining_fraction`, GEN.757). The retention row `d6` caps the total. An icy cold body then holds the guest gases above their dissociation pressure in the cages of its ice skin (GEN.762, `clathrate_brine.hydrate_hold_pa`).
4. The pressure is the sum of what is left (`P = g M / (4 pi R^2)`, `surface_pressure.pressure_pa`); the density, mean molar mass and scale height follow (`atmosphere_state.air_state`).
5. The greenhouse and the temperature of that air (`atmosphere_state.greenhouse_warming`).

The step is damped in ln P and in T and halves its damping when it overshoots; a creep toward the solution (the water vapour and the temperature feeding each other) is taken by a secant step of at most `MAX_GAIN` times the damped one. It stops when ln P moves by less than `TOLERANCE_LN_P` (1e-3, the noise of the greenhouse solver) and T by less than 0.1 K, or after `MAX_STEPS`; `Air.converged` says which. A world whose ocean freezes and thaws more than `MAX_OCEAN_FLIPS` times has no state that is both (the ocean weathers the CO2 away and freezes, the ice leaves the CO2 to build up and thaws): it is solved as one that keeps its ocean, with the weathering law taking itself away where it is cold.

**Giants, brown dwarfs, envelopes.** There is no surface: the reference level is the 1 bar level, the mixture there does not depend on the surface temperature, and the greenhouse temperature there is of the radiated flux, the absorbed starlight plus the giant's own heat flux (`giant_draw.own_flux_w_m2`, the cooling law through Jupiter and Neptune, GEN.816; Jupiter 154 K, Uranus 83 K, Neptune 76 K against 166, 76 and 72 measured, Saturn 115 K against 134 until its helium-rain heat is added, GEN.810). One step is the solution.

**Where the old values went.** The pressure of a class-A to V rocky world is now the supply of its interior less its sinks; a class whose retention row says "none kept" (d1) is airless; the stored density is the mean molar mass over the air's own state, not a class range.

**Registry rows this adds** (GEN.725 rule, in `research/gap-filling-plan.md`):

| name | value | source | tier | tweak |
|---|---|---|---|---|
| `RELATIVE_HUMIDITY` | 0.77 | [R2848] | 1 verified | 0 |
| `LONG_LIVED_FLUX_FRACTION` | 0.5 | share of the Earth twin's heat flow that keeps a mantle geologically alive | 7 educated guess from inferred physics | 0 |

| `MAX_OCEAN_FLIPS`, `MAX_GAIN`, `MAX_STEPS`, tolerances | | numerical | not stand-ins | |

Removed stand-ins: the weathering sink bracket of the rocky inventory (now the carbonate balance and Henry's law), the class density ranges, the greenhouse multiplier ranges.
