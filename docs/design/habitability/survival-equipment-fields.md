# Survival equipment fields: a combinable Body + Filter + Supply system for the habitability index (GEN.655, DB.93)

Research thread "Survival equipment fields". Written for Boss's message of 6:56 PM EDT 2026-10-23: "change the equipment system to be a +/- system where we just combine the different things a person could need and then derive a name from that. Categories of equipment are: Body Protection (like a suit, the lowest kind of suit that fits the environment) / Air Filter or Scrubber / Air Supply. Research first and give me a suggestion for the equipment fields for a world's habitability index to provide an accurate idea of survival equipment for the planet based on the possible needs."

**Status.** Built in 8.1 (GEN.655 the equipment module and its inputs, DB.93 the columns, API object and migration 0136; GEN.656 the regional temperatures, migration 0137, section 5.5). This page was the design draft of 2026-10-23, revised with Boss's answers (section 2); where the build differs from the draft, the difference is in section 14. Still to come, each its own item: the chips (UX.154), the wiring of the toxic-species rows to the stored gases (GEN.657, the limits are verified; the wiring waits for the generation lock, section 14) and the dust input's toxicity half (GEN.658 built the loading half; toxicity waits for the surface composition to be stored, GEN.441). The numbers in sections 8 and 9 come from a scratch prototype run against the repo's own physics on 17 named worlds and on 2,261 planets and moons from 160 generated systems; the module's own tests (`tests/test_equipment.py`) hold the named worlds. Values marked DESIGN are planetGen defaults meant to be tuned by a sweep; values marked RECALLED are from memory and must be checked against the source before they are leaned on; section 12 says which is which.

## 1. The short version

1. **Replace the one number (`equipment_tier`, 0 to 4) with items that add up**, in four slots: **Body Protection** (the lowest suit that fits, 0 to 6, with the stress that set it, plus a powered frame for high gravity), **Air Filter or Scrubber** (dust, a gas filter of 3 grades, a CO2 scrubber), **Air Supply** (0 to 4) and **Radiation Protection**. Radiation climbs a five-step ladder modelled on real gear: **shielded garment** (the apron and Demron class), **storm vest** (the AstroRad class flown on Artemis I), **hardened suit** (game technology), then **shielded shelter** and **hardened habitat** (Boss, A and F).
2. **Only what is needed is listed, and each item says how badly.** No "none" ever appears. Every item is **required** (past the physical limit) or **recommended** (the fuzzy band where a person survives but is impaired or at risk, or where only a region or the night asks for it). A world that needs nothing shows no equipment line at all (Boss, E).
3. **The physics always decides the gear** (Boss, B). Where the PHI-6 bands are looser than the physics (a filter in air past IDLH, the 10 MPa hyperbaric edge against the gas-density ceiling, H2S's IDLH), the equipment follows the physics and the scores stay as they are.
4. **Average temperature sets the habitability; regions and the night go in the text** (Boss, C). Each world gets four regional temperatures (Polar, Northern, Equatorial, Southern) and a day and night swing. The score and the required gear read the average; the detail text names any region or hour that needs more, and that heavier gear is listed as recommended.
5. **The name is derived, never stored**, as two short lists in Boss's order (Body, Filter, Supply, Shelter): Mars reads "Required: pressure suit + dust-tight seals. Recommended: hardened habitat." A short code (`B4p F100 S2 [R2]`, recommended items in brackets) is what the API and search carry (Boss, D).
6. **One source of truth.** The new `physics/equipment.py` returns a list of needs once; the stored levels, the describer sentences and the species fit all read it, replacing the parallel `habitability.equipment()` and `protection.requirements()`.
7. **What the sample shows (160 systems, 2,261 rocky bodies).** 2,009 need a pressure suit, 100 need armour or a vehicle, 14 need no body protection. Radiation, required: nothing 1,025, shielded garment 214, storm vest 84, hardened suit 118, shielded shelter 301, hardened habitat 519. Every Red-dose body in the sample is a storm-particle world (a flaring star), which is exactly where suits help; on cosmic-ray worlds like Mars no suit helps and the draft lists only the habitat.

## 2. Boss's decisions (7:37 PM EDT 2026-10-23)

Boss answered in text and then cleared the cards; the text answer is the decision ("my text answer stays", 7:39 PM EDT).

| | Question | Boss's answer | What it changed in this draft |
|---|---|---|---|
| A | A fourth slot for radiation? | "1 and 2, so radiation shelter is a higher level than just hardened suit, within the game environment such suits exist." | Both: a **hardened** modifier on Body Protection (a radiation-hardened suit exists in the game) and the **Radiation Shelter** slot above it. Ladder: hardened suit, shielded shelter, hardened habitat (5.1, 5.4) |
| B | How far does a filter reach? | "1 and always follow the physics" | Gear follows OSHA and NIOSH and the physical limits everywhere they are stricter than the built table; PHI-6 bands and scores unchanged, no epoch (5.2, 5.3, 8.1) |
| C | Day and night temperature | "1 & Day / Night and regions (Polar, Northern, Equatorial, Southern) habitability is based on the average surface temperature with details about regions within the detail text if needed." | Swing estimate plus four regional temperatures, stored; score and required gear from the average; regions and night in the detail text, heavier gear they need listed as recommended (5.0, 5.5) |
| D | Name and code format | "1, Yes" | Words joined with plus signs, plus a short code |
| F | Radiation suit levels (7:46 PM EDT) | "let's add to the radiation part different levels so look up levels of actual radiation suits and incorporate that before hardened suit." | Two levels from real gear before the hardened suit: **shielded garment** and **storm vest**. Each level is a layer of shielding (g/cm2) run through the repo's own dose curves, so the gear needed follows the physics of each world's dose (4.5, 5.4) |
| E | One sortable number? | "1, We never want to list what isn't required ... it lists only things needed such as 'pressure suit recommended' or required depending on severity, we can use that for edge cases where it's a bit fuzzy." | No single number; no "none" anywhere; every listed item carries **required** or **recommended** (6.4) |

How A and F now work, to confirm or correct: the **required** radiation gear is the lowest step of the ladder that brings the year's dose under 1 Sv (out of Red); the **recommended** gear is the lowest step that brings it under 0.1 Sv (Green). Where a shelter is needed, the suit for going outside is judged on a working year outside (2,000 hours). One change from my earlier reading, because the physics decides (B): a suit does almost nothing against cosmic rays, so a Yellow cosmic-ray world (Mars, the Moon) no longer lists a hardened suit. It lists the hardened habitat as recommended, and the Moon a shielded garment as well.

Small choices I made and will keep unless you object: oxygen from 50 to 160 kPa is a gas-mix mask, not an environment suit (5.3); ambient CO2 from 0.93 to 5 kPa is a scrubber on the intake, above that own air (5.2); an open mask and a filter cannot be combined, the sealed supply replaces them (6.1); a locked world's regions read day side, terminator and night side instead of the four latitude bands (5.5); UV, gravity and armour thresholds are DESIGN values (5.1).

## 3. What exists today, and what it cannot say

The built line (`habitability.equipment`, stored as `equipment_tier`, section 8 and 9 of `docs/design/habitability-index.md`) is one of five names: ideal, breathing mask, mask with scrubber, environment suit, full life support with radiation hardening. It is the heaviest need of Pressure, Temperature, Air and Radiation.

What a single tier cannot do, with the cases that show it:

- **It cannot combine.** A cold world with a toxic gas and a thin sky is "environment suit"; so is a vacuum, a Titan and a Venus. A visitor learns nothing about whether they need a heater, a pressure shell, or only air.
- **It ignores what it has no row for.** Only CO2, CO, H2S and SO2 count as toxic. Ammonia, hydrogen chloride, hydrogen fluoride, chlorine, nitrogen dioxide, hydrogen cyanide and ozone are in the research notes and not in the code. Dust, UV and gravity are not read at all.
- **It misses real cold and heat.** The cold edge of the suit is -50 C (an environment suit), but -35 C on Earth still reads "ideal" and so does a 40 C day (Green, but the day wants cooling clothing).
- **"Mask with scrubber" is the wrong object.** A scrubber removes CO2 (lithium hydroxide, soda lime, amine); a cartridge filter removes SO2, H2S, ammonia and similar by adsorption. They are different devices with different limits.
- **"Environment suit" spans the whole range** from a liquid-cooled suit at 60 C to an armoured suit at 460 C.
- **Two implementations.** `habitability.equipment()` decides the tier; `protection.requirements()` decides the gear named in the sentence. They agree today (a test checks it) but each has its own copy of the thresholds' logic.
- **The filter rule is looser than the physics.** The PHI-6 Air bands let a cartridge filter work up to ten times the acute limit (`FILTER_FACTOR`, "NIOSH APF 10"). OSHA, NIOSH and the repo's research notes forbid air-purifying respirators in atmospheres Immediately Dangerous to Life or Health; an assigned protection factor of 10 multiplies the exposure limit, not the IDLH. For H2S the band allows 0.5 kPa (500 times the long-term limit).

## 4. The science, by category

What each piece of gear does and where it stops. Status: **V** checked on the web this session; **N** from the research notes in `docs/design/` (Boss's imported research, not re-checked); **R** recalled, check at build time; **D** a planetGen design value.

### 4.1 Body protection (what shields the skin and the body from pressure, heat, cold and chemistry)

| Stress | Limit and what fails | Gear | Status |
|---|---|---|---|
| Low pressure | Below 6.3 kPa (water's vapour pressure at 37 C, the Armstrong limit) body fluids boil: ebullism, circulatory arrest in 10 to 15 s | Full gas-pressurised suit | N |
| Low pressure | A mask on pure oxygen holds a safe alveolar oxygen down to about 20 kPa; below that the lungs need positive pressure, and the chest tolerates only about 4 kPa of it (6.7 briefly) without counterpressure, so the whole body must be squeezed | Counterpressure or gas-pressurised suit | N |
| Suit pressure | The ISS spacesuit runs at 29.6 kPa (4.3 psi) of pure oxygen; mechanical counterpressure suits aim for the same 30 kPa | A suit is a 30 kPa shell | V (npj Microgravity 2021; MIT BioSuit) |
| High pressure | Nitrogen narcosis from 300 to 400 kPa of N2 (argon about twice as narcotic); high-pressure nervous syndrome from 1.5 to 3 MPa; breathing-gas density above 6 g/L halves work capacity, above 10 g/L the lungs cannot move gas | Mixed gas up to about 5 to 7 MPa; beyond it a rigid suit that keeps 101 kPa inside | N |
| Heat | The wet-bulb limit for a resting young adult is about 31 C, not 35 C: the Penn State lab measured 30.55 plus or minus 0.98 C | Cooling garments; a liquid-cooled garment above 45 C in sealed gear | V (Vecellio et al. 2022); N |
| Heat | A mask alone survives about 120 to 140 C of dry air and 50 to 55 C of steam before skin at the seal burns; above that a fully enclosed suit with its own cooling loop | Heat-shielded suit with a liquid cooling garment | N |
| Extreme heat | Venus's 735 K: a suit is a short-lived thermal capacitor (refrigeration into molten-salt storage), sorties of 1 to 2 hours | Armoured thermal suit | N |
| Cold | Exposed face skin freezes in 30 to 60 s below about -45 to -50 C with wind; mask seals and valves fail at -50 to -60 C; cold gas chills lung blood (breathing gas must be warmed on Titan) | Insulated suit to -50 C; heated, aerogel-insulated suit below | N |
| Corrosive gases | Hydrogen fluoride, strong acids, ammonia burn and are absorbed through skin; a mask leaves the skin exposed | Fully encapsulating chemical-barrier suit | N |
| UV | The stored index is DNA-weighted UV against Earth's 1; Earth's skin burns in minutes at the noon peak, so a world at 30 times that burns in seconds | Cover-up clothing, then an opaque suit | D |
| Gravity | Not in the notes. Sustained loads above about 1.5 g strain the body; above about 3 g no person can work | A powered frame; a vehicle | D |

### 4.2 Air filter or scrubber (what is taken out of the air before it is breathed)

| Contaminant | Device | Limit | Status |
|---|---|---|---|
| Dust and aerosol | Particle filter (P100); dust-tight seals and a suitport for a sealed suit | Mars dust: 0.1 mg/m3 for 30 days (particles under 10 micrometres), ceiling 10 mg/m3; lunar silicate is razor sharp and wrecks seals | N (NASA-STD-3001) |
| Toxic gases (SO2, H2S, NH3, HCl, HF, Cl2, NO2, HCN, O3) | Chemical cartridge on a half mask (assigned protection factor 10), a full facepiece (50) or a powered unit (up to 1,000) | Useful from the exposure limit up to ten, fifty or a thousand times it, never into IDLH; needs 19.5 percent oxygen at sea level (about 19.8 kPa) | N; the numbers are OSHA's and I could not re-open them (R) |
| Carbon monoxide | Hopcalite catalyst: escape devices only, about 2 hours, exothermic | Not a working filter | V (MSA W65 sheet) |
| Carbon dioxide | A cartridge does not take it out. A chemical scrubber does (2 LiOH + CO2 gives Li2CO3 + H2O; soda lime; amine beds), as the loop of a rebreather or backpack | Long-term 0.93 kPa; NASA's one-hour spaceflight limit 0.4 kPa; 4 kPa is IDLH; above 5 kPa acute; above 10 kPa a narcotic. An intake scrubber is workable up to about 5 kPa, a mass of absorbent per hour that grows with the partial pressure | N, D |

### 4.3 Air supply (where the breathing gas comes from)

| Need | Device | Limit | Status |
|---|---|---|---|
| Too little oxygen | Mask with supplemental oxygen | Inspired O2 below 8 kPa gives collapse in 10 to 15 s; below 16 kPa judgement fails (the notes' "hypoxic alert"); the built table uses 8 kPa (the more permissive) | N |
| Too much oxygen | Dilute the gas, or breathe a mix | Above 50 kPa for days damages the lungs; 140 to 160 kPa causes convulsions | N |
| Gas that no filter may meet | Own gas: supplied air through a hose, a self-contained bottle (30 to 60 min), or a closed loop with its own scrubber and oxygen make-up (6 to 8 h) | An air-purifying respirator is forbidden in IDLH and in an oxygen-deficient atmosphere | N |
| High pressure | Heliox or trimix at ambient pressure (saturation) | Narcosis at 300 to 400 kPa of N2; density ceiling about 5 to 7 MPa | N; the density figure is checked in 4.4 |
| Beyond the ceiling | A rigid suit or vehicle that breathes at 101 kPa | Hydra X record 7.11 MPa; atmospheric diving suits at 1 atm inside | N |
| Breath temperature | Warm or cool the inspired gas | Below about -50 C or above about 140 C the air itself cannot be breathed | N |

### 4.4 A check I made on the density ceiling

For a helium and oxygen mix holding 40 kPa of oxygen at 300 K, the breathing-gas density reaches 8.5 g/L at 5 MPa, 10.1 g/L at 6 MPa and 11.7 g/L at 7 MPa, so the notes' "5 to 7 MPa" is right. It depends on temperature (density goes as 1/T), and the built table's 10 MPa hyperbaric edge (the more permissive of 5 and 10) is looser than the physics. Following the physics (decision B), the gear reads the computed density: an armoured suit is recommended from 6 g/L and required from 10 g/L, whatever the table's pressure edge says. The score keeps its 10 MPa edge.

### 4.5 Radiation (what a person cannot wear)

ICRP: 20 mSv per year averaged over five years (occupational), 50 mSv in any one year (the table's Blue edge and single-year limit); 100 mSv per year a visit; 1 Sv per year is the threshold of radiation sickness. NASA's career limit is 600 mSv, with 250 mSv per solar storm event (V, NASA OCHMO technical brief). Mars gets about 0.7 mSv per day from galactic cosmic rays (N). Wearable shielding cuts the blood-forming-organ dose in a storm by 20 to 30 percent and does nothing against cosmic rays; heavy metal in a suit makes it worse (secondary neutrons); real shielding is hydrogen-rich (water, polyethylene, regolith) and thick (N). Burying a base cuts the cosmic-ray and storm doses and not the dose from the ground (the stored `dose_ground_msv_yr` is the part burial cannot remove).

**Real radiation suits, from least to most (decision F).**

| Real gear | Shielding | What it stops | What it does not stop | Game level |
|---|---|---|---|---|
| Anti-contamination coverall (Tyvek class, worn with a respirator) | none to speak of | radioactive dust on the skin and, with the respirator, in the lungs | any penetrating radiation | not a radiation level: dust-tight seals and the filter (5.2) do this job |
| Spacesuit (EMU class) | thin layers: helmet shell about 0.18 g/cm2, visor shells to 0.38; the water cooling tubes cover only 40 percent of the body (V, NASA) | electrons and low-energy protons | storm cores, cosmic rays | the baseline of every sealed suit, not a level |
| Lead or lead-free apron, 0.25 to 0.5 mm lead equivalent; Demron full suit | about 0.3 to 0.6 g/cm2 | 90 to 97 percent of 100 kVp X-rays (V, apron studies); Demron at least 50 percent of gamma up to 130 keV (vendor claim) | about 6 percent of 662 keV gamma (DERIVED, mass attenuation about 0.11 cm2/g, R); protons above about 25 MeV; heavy metal makes secondaries in particle fields | **1 shielded garment**, a 1 g/cm2 hydrogen-rich layer in the game (lead is the wrong material for particles) |
| Partial-body gamma shield (StemRad 360 Gamma) | heavy shielding over the pelvis only | enough bone marrow to recover from an acute gamma dose (R, vendor) | the rest of the body | the idea behind level 2: shield the organs that matter most |
| Storm vest (AstroRad, high-density polyethylene, 26 kg, flown on Artemis I) | several g/cm2 over the torso | an estimated 40 to 60 percent of the effective dose of a great solar storm (V, Artemis I results, extrapolated to the 1972 and 1989 events) | cosmic rays | **2 storm vest** |
| Hardened suit | game technology | | | **3 hardened suit** |
| Storm shelter (water, regolith, polyethylene) | 20 to 50 g/cm2 | almost all of a storm | half of the cosmic-ray dose | **shielded shelter** |
| Buried habitat | about 300 g/cm2 (about 2 m of regolith) | cosmic rays to the Green band | the ground's own dose | **hardened habitat** |

## 5. The slots: levels, thresholds and the fields that feed them

### 5.0 The conditions every slot reads (all from the stored row unless marked)

| Condition | Source | Stored today? |
|---|---|---|
| P, total pressure, kPa | `atmospheric_pressure_pa` | yes |
| T mean, C; coldest and hottest over the orbit | `surface_temperature_k`; `habitability.temperature_extremes_c` (from `flux_variation`) | yes |
| **Day and night swing** | New: damped-oscillator estimate from the rotation period, air column heat capacity (c_p P / g) and its radiative time constant, regolith inertia if airless (5.5) | no (decision C) |
| **Regional temperatures** (Polar, Northern, Equatorial, Southern) | New: annual-mean insolation by latitude from the obliquity, flattened by the air's heat transport (5.5) | no (decision C) |
| Wet bulb | Stull formula from T and relative humidity (the H2O partial pressure over saturation) | derived |
| Partial pressure of each gas | the child tables `planet_gases` and `moon_gases` (DB.71, GEN.424): every gas of the table, the balance gas included | yes |
| Inspired oxygen | O2 fraction times (P minus 6.3 kPa) | derived |
| Breathing-gas density of the best mix | heliox with 40 kPa O2: P M / (R T) | derived (`atm_molar_density` is the mean molar mass) |
| Dose, Sv/yr, and its parts | `surface_dose_msv_yr`, `dose_gcr`, `dose_sep`, `dose_ground`, `dose_helio_mult` | yes |
| UV | `uv_surface_index`, `ozone_loss_flag` | yes |
| Gravity | `gravity_g` | yes |
| Dust loading and toxicity | `tau_dust` (GEN.414) and `scale_height_km`, as the airborne loading `4 rho r tau / (3 Q H)` (`aerosols.dust_loading_mg_m3`, GEN.658: 2,500 kg/m3 grains, 1.5 micrometres, Q 2.5, DESIGN); the toxicity of the dust waits for the stored surface composition; a rocky body under 10 kPa still has regolith dust | loading yes |
| Corrosive surface liquid | `ocean_ph` (acid below 2, alkaline above 12) | yes |

### 5.1 Slot 1: Body Protection (`eq_body`, 0 to 6, and `eq_body_kind`)

The levels nest: each one does everything the one below it does, so "the lowest suit that fits" is the lowest level that meets every stress. The level is the **maximum** over the sub-needs; the sub-need that set it is the **kind** (cold, heat, pressure, chemical, diving, UV, gravity), which gives the name its word. The required level reads the average surface temperature; a recommended level (5.5, 6.4) can sit above it.

| Lv | Name | What it is | Triggers (any one; limits-table names in brackets) |
|---|---|---|---|
| 0 | none | ordinary clothes | none of the rows below |
| 1 | garment | weather clothing: warm, cooling or sun-cover | T from -20 up to 0 C (cold_clothing, comfort_low); T from 31 to 45 C (comfort_high, heat_garment); wet bulb at or above 31 C; UV index 3 to 30 (DESIGN) |
| 2 | thermal suit | insulated or liquid-conditioned suit, not sealed; mask or hood to breathe | T from -50 to -20 C (cold_suit); T from 45 to 140 C (heat_garment, heat_shielded); 400 to 10,000 kPa (the heated diving suit helium needs); UV index over 30 |
| 3 | sealed suit | fully encapsulated, own microclimate, ambient pressure | T below -50 C (heated, aerogel); T from 140 to 400 C (heat-shielded, aluminised, liquid cooling); a corrosive or skin-absorbed gas (HF, HCl, NH3, Cl2) over its IDLH; a surface liquid below pH 2 or above 12 |
| 4 | pressure suit | counterpressure or gas-pressurised suit (a 30 kPa shell) with thermal control | P below 14.3 kPa (mask_min); the form is counterpressure from 6.3 to 14.3 and full gas-pressurised below 6.3 (armstrong) |
| 5 | armoured suit | rigid suit that keeps 1 atm inside, thermal-capacitor class | P over 10,000 kPa (pressure_hyperbaric_max); T over 400 C (Venus class, sorties of hours) (DESIGN) |
| 6 | vehicle only | no suit exists for it | T over 500 C; gravity over 3 g; P over 100 MPa (DESIGN) |

**Severity bands (physics first, decision B).** Required past the physical limit; recommended in the band before it:

| Stress | Recommended | Required | Source of the band |
|---|---|---|---|
| Low pressure | pressure suit from 20 down to 16 kPa | pressure suit below 16 kPa | a mask on pure oxygen is unusable below 16 to 20 kPa (notes; index design 7.2). The score keeps its 14.3 kPa Green edge |
| Cold | heated sealed suit from -40 to -50 C | insulated suit from -20 C, heated sealed suit below -50 C | notes give -40, the table -50 |
| Humid heat | cooling clothing at a wet bulb of 30 to 31 C | at 31 C or more | Penn State's 30.55 plus or minus 0.98 C |
| High pressure | armoured suit at a best-mix density of 6 to 10 g/L | over 10 g/L | work capacity halves at 6 g/L, ventilation fails at 10 |
| Gravity | powered frame from 1.25 to 1.5 g | from 1.5 g (vehicle only over 3 g) | DESIGN |

**Radiation suit levels** (shielded garment, storm vest, hardened suit) sit in their own slot (5.4); a hardened suit hardens whatever suit the world needs ("hardened pressure suit").

Add-on `eq_mobility`: a powered frame (recommended from 1.25 g, required from 1.5 to 3 g, DESIGN).

New limits-table rows: `armoured_heat_c` 400, `vehicle_heat_c` 500, `frame_gravity_g` 1.5, `vehicle_gravity_g` 3.0, `uv_cover_index` 3, `uv_suit_index` 30, `corrosive_ph_low` 2, `corrosive_ph_high` 12, `vehicle_pressure_kpa` 100,000.

### 5.2 Slot 2: Air Filter or Scrubber (three parts that combine)

| Part | Field | Values | Triggers |
|---|---|---|---|
| Dust | `eq_dust` | 0 or 1 | Airborne dust over 0.1 mg/m3 recommended and over 10 mg/m3 required (the loading from `tau_dust`, GEN.658); a rocky body under 10 kPa has regolith dust and the filter is required. Renders as "dust filter" with an open mask and as "dust-tight seals" when the supply is sealed (a sealed suit needs a seal and suitport, not a filter) |
| Gas | `eq_gas_filter` | 0 to 3 | The worst filterable gas, as a ratio to its long-term limit (HR = partial pressure over chronic): 0 none; 1 gas filter (HR up to 10, half mask); 2 full-face gas filter (HR up to 50); 3 powered gas filter (HR up to 1,000). Never over IDLH: that is Supply 2 |
| CO2 | `eq_co2_scrub` | 0 or 1 | Ambient CO2 above 0.93 kPa and up to 5 kPa (a bulk absorbent on the intake). Above 5 kPa the person breathes own gas |

Toxic table (new rows in `protection_limits`, long-term limit and IDLH, kPa partial pressure; every row read from the NIOSH Pocket Guide for GEN.657, see section 14; the draft's recalled NO2 IDLH of 20 ppm is 13 ppm in the newest NIOSH document):

| Gas | Long-term | IDLH | Gas | Long-term | IDLH |
|---|---|---|---|---|---|
| CO | 0.005 | 0.12 | HCl | 0.0005 | 0.0046 (NIOSH 2025: 45 ppm) |
| H2S | 0.001 | 0.0101 (built acute 0.05) | Cl2 | 0.0001 | 0.0010 |
| SO2 | 0.0005 | 0.0101 | NO2 | 0.0001 | 0.0020 |
| NH3 | 0.0025 | 0.030 | HCN | 0.0005 | 0.0051 |
| HF | 0.0003 | 0.0030 | O3 | 0.00001 | 0.0005 |

The built table's acute values (0.05 kPa for H2S, 0.1 for CO) are looser than NIOSH's IDLH for H2S (100 ppm, 0.0101 kPa). Boss's rule of 2026-10-14 (the more permissive value wins) keeps 0.05; the table above shows both so the choice is visible.

### 5.3 Slot 3: Air Supply (`eq_supply`, 0 to 4)

| Lv | Name | What it is | Triggers (any one) |
|---|---|---|---|
| 0 | ambient air | breathe the air | none below |
| 1 | gas-mix mask | open mask that adds oxygen or dilutes it | inspired O2 under 8 kPa with P at or over 14.3 kPa; inspired O2 over 50 kPa (to the 160 acute limit); P from 250 to 400 kPa (the Green band's mix mask) |
| 2 | own gas | sealed full mask or suit on its own gas: hose, bottle or closed loop with its scrubber | P under 14.3 kPa (a suit's loop); O2 over 160 kPa; CO2 over 5 kPa; any filterable gas over IDLH, or over 1,000 times its long-term limit; breathing air colder than -50 C or hotter than 140 C; and any time level 1 and a filter are both wanted (the open mask cannot be sealed for a filter) |
| 3 | mixed gas | heliox or trimix at ambient pressure, density-managed | P from 400 to 10,000 kPa |
| 4 | isolated | the wearer breathes 101 kPa inside a rigid suit or vehicle | P over 10,000 kPa; Body Protection 5 or 6 |

Body Protection 4 carries a level-2 loop with it (the suit's own); 5 and 6 carry level 4. So the supply is not named again beside a pressure suit unless it is mixed gas (hyperbaric inside a suit never happens in the sample).

### 5.4 Slot 4: Radiation Protection (`eq_rad_suit` 1 to 3, then `eq_shelter` 1 or 2)

A five-step ladder (decisions A and F). Each step is a layer of shielding, in g/cm2 of hydrogen-rich material, laid on top of the world's own air column and run through the repo's own dose curves (`radiation.gcr_dose_msv_yr`, `radiation.sep_dose_msv_yr`) against the stored parts of the dose (`dose_gcr_msv_yr`, `dose_sep_msv_yr`, `dose_ground_msv_yr`). Nothing is drawn: the gear follows each world's mix of cosmic rays and storm particles.

| Step | Name | Shielding (DESIGN, from the real gear in 4.5) | Storm cut, airless body | Cosmic-ray cut, airless body | Storm cut under Mars's air |
|---|---|---|---|---|---|
| 1 | shielded garment | 1 g/cm2 over the whole body | 71 percent | 17 percent | 5 percent |
| 2 | storm vest | the garment plus 5 g/cm2 over the torso (60 percent of the effective dose) | 81 percent | 31 percent | 17 percent |
| 3 | hardened suit | 5 g/cm2 over the body plus 5 more over the torso | 89 percent | 43 percent | 30 percent |
| 4 | shielded shelter | 50 g/cm2 | 98 percent | 70 percent | 79 percent |
| 5 | hardened habitat | 300 g/cm2, or as deep as the world needs (the derived depth below) | all | 94 percent | 99 percent |

The cuts are computed, not chosen: storm particles are soft, so on an airless body even a garment stops most of them; under an air column that has already stopped the soft part, thin layers do little. Cosmic rays need a habitat.

**The rule.**

- **Required:** the lowest step that brings the year's dose to 1 Sv or less (out of Red).
- **Recommended:** the lowest step that brings it to 0.1 Sv or less (Green), when it is above the required step.
- When no shelter is needed, a suit step is judged on the whole year (an unshielded base gives no more cover than a suit). When a shelter or habitat is needed, it is judged on the whole year, and the gear for going outside is the lowest suit step that keeps a **working year outside**, 2,000 hours (DESIGN), under the same limit; when even a hardened suit cannot, the hardened suit is listed and the text gives the **hours outside a year** it allows ("Even in a hardened suit, 1 Sv is reached after 800 hours outside a year").
- A step that does not lower the dose enough is never listed, so no world lists a suit that would not help (decision E).

Stored beside the levels: `shelter_column_g_cm2`, the shielding a shelter needs to bring the dose to the limit (it tells a player how much regolith or water to pile on), and `outside_hours`, the hours outside a year the listed suit allows. Burial does not reduce the ground dose.

### 5.5 Regions and the day and night swing (decision C)

Habitability, the PHI-6 score and the **required** gear read the average surface temperature. Four regional temperatures and a day and night swing are worked out and stored for the detail text:

- **Regions.** Polar (latitude 60 degrees and over, 13 percent of the surface), Northern (20 to 60 N, 26 percent), Equatorial (20 S to 20 N, 34 percent), Southern (20 to 60 S, 26 percent). Each takes the annual-mean insolation of its band, which depends on the obliquity (from the stored spin axis); a thick atmosphere carries heat poleward and flattens the differences (a factor from the air column, DESIGN, fitted so Earth's poles sit about 40 K under its mean and an airless body keeps the full radiative contrast). Northern and Southern differ when the orbit is eccentric and perihelion falls near a solstice; their seasonal extremes come from the obliquity.
- **Day and night.** A damped-oscillator estimate: the day length against the air's radiative time constant (column heat capacity c_p P / g over 4 sigma T^3), and the regolith's thermal inertia for an airless body. It gives each region a day high and a night low.
- **Locked worlds.** A tidally locked world has no day and night; its regions read day side, terminator and night side (small choice, section 2).
- **In the text.** "The equator averages 12 C; polar nights fall to -70 C, where a heated sealed suit is recommended." When a region or the night crosses a limit the average does not, the heavier gear is listed as **recommended** with the place named.

Stored (GEN.656; since GEN.1108, schema v196, as the seven bands of `physics/latitude_bands.py`: `band1_mean_k` to `band7_state`, 49 columns, and `temp_mean_swing_k`, in place of the four named regions and the three zones; see database-schema.md) on planets and moons (their own columns, so search and the text can read them). Derived from stored values only: no new draws.

**As built (GEN.656, `physics/regional_temperature.py`).** The annual-mean insolation of each band is integrated over a circular orbit for the stored axial tilt, turned into a temperature by the fourth-root law and scaled so the area-weighted bands average the stored mean. The air then removes the share `x / (x + 20)` of the contrast, `x` the air column (pressure over gravity) in Earth's: Earth keeps its poles about 35 K under the mean (observed 35 to 40 K), Titan loses 30 percent and Venus 80 percent of the contrast. The orbits are circular, so Northern and Southern are equal. A tidally locked world (rotation period equal to the orbital period) has no bands: its day side, terminator ring and night side take the radiative multiples 1.131, 0.5 and 0 of the mean, pulled toward the mean by `x / (x + 0.1)` (a bar of air removes 90 percent of the contrast; the night floor is 20 K). Columns `temp_day_side_k`, `temp_terminator_k` and `temp_night_side_k` hold them, so the latitude group and the locked group are never both set. The swing is `sqrt(2) T / sqrt(1 + (omega tau)^2)` with `tau = (c_p P / g + I sqrt(2 / omega)) / (4 sigma T^3)`, ground inertia `I` 300 J/(m2 K s^0.5) and `c_p` 1000 (both DESIGN stand-ins), times one minus the transport share. Checked in `tests/test_regional_temperature.py` against Earth, Mars, the Moon, Mercury and Venus. Not built yet: the seasonal extremes of the eccentric orbit, and naming the heavier gear for a night that crosses a limit the mean does not (the Temperature row names a place whose temperature needs a heavier class than the mean).

## 6. How the combination is built

### 6.1 The algorithm

1. Each stress rule returns a **need**: `(slot, level or part, kind, domain, reason)`. The reason is the sentence the describer prints ("the air pressure is 0.6 kPa"). The domain is the PHI-6 domain the need comes from.
2. Per slot, the level is the maximum of its needs. The Filter slot keeps three independent parts.
3. Cross-slot rules, applied in this order: (a) Body 5 or 6 sets Supply 4; Body 4 sets Supply to at least 2. (b) If Supply is 1 and a gas filter or the CO2 scrubber is also needed, Supply becomes 2 (an open mask cannot also be sealed for a filter). (c) If Supply is 2 or more, the gas filter and the CO2 scrubber drop out (the person is not breathing the ambient air; the loop scrubs its own CO2); the dust part stays and renders as seals. (d) A shielded garment or storm vest is worn over whatever Body the world needs, even ordinary clothes; a hardened suit raises Body to at least 3 (a sealed suit). (e) Shelter is independent. (f) Each item keeps its severity; an item both required and recommended at different levels lists the required one and the heavier one as recommended.
4. Name and code are formed from the final slots (6.2). Nothing but the levels is stored.

### 6.2 The name

Two lists, each joined with " + " in Boss's order (Body, Filter parts, Supply, Shelter): **Required: ...** and **Recommended: ...**. A list with nothing in it is left out, and a world with nothing in either shows no equipment line (decision E). A powered frame goes after the body word; a radiation suit step joins the body word: "pressure suit + storm vest", or "hardened pressure suit" (or "hardened suit" with no other suit).

- Body words by level and kind: 1 cold "cold-weather clothing", heat "cooling clothing", UV "sun-cover clothing"; 2 cold "insulated suit", heat "liquid-cooled suit", diving "heated diving suit", UV "full-cover suit"; 3 cold "heated sealed suit", heat "heat-shielded suit", chemical "chemical-barrier suit"; 4 "pressure suit"; 5 "armoured suit"; 6 "vehicle only".
- Filter words: "dust filter" (or "dust-tight seals" when Supply is 2 or more); "gas filter", "full-face gas filter", "powered gas filter"; "CO2 scrubber".
- Supply words: "gas-mix mask", "self-contained air", "mixed-gas rig", "isolated air". Not named when the body already carries it (pressure suit, armoured suit).
- Radiation words: "shielded garment", "storm vest", "hardened suit"; "shielded shelter", "hardened habitat".

The code lists required items, then recommended ones in brackets: `B<level><kind> F<dust><gas><co2> S<level> H<suit step> R<shelter level>`, then `[...]` with only what is recommended. `H1` shielded garment, `H2` storm vest, `H3` hardened suit; `R1` shielded shelter, `R2` hardened habitat. Kind letters c cold, h heat, p pressure, x chemical, d diving, u UV, g gravity. Mars: `B4p F100 S2 [R2]`.

### 6.3 The strip

Only listed items get a chip: a solid chip for required, an outlined chip for recommended. Nothing is drawn for a slot that needs nothing.

```
Mars        Body: pressure suit (required)   Filter: dust-tight seals (required)   Radiation: hardened habitat (recommended)
Arctic day  Body: insulated suit (required)
Earth       (no equipment line)
```

### 6.4 Required or recommended

An item is **required** when its condition passes the physical limit in 5.1 to 5.4 at the average surface temperature. It is **recommended** when the condition sits in a severity band before that limit (5.1), when it reaches the Green dose band (5.4), or when only a region or the night needs it (5.5). The describer says why in one short sentence each ("Recommended: a pressure suit, because at 18 kPa a mask on pure oxygen is marginal.").

## 7. Fields to store, and how they are served

Replace `planets.equipment_tier` and `moons.equipment_tier` (TINYINT, 0 to 4) with, on both tables (NULL for a gas giant and for an unscored body, like the rest of the score):

| Column | Type | Meaning |
|---|---|---|
| `eq_body`, `eq_body_kind` | TINYINT 0 to 6, 0 to 7 | required Body Protection level and what set it (none, cold, heat, pressure, chemical, diving, UV, gravity) |
| `eq_body_rec`, `eq_body_rec_kind` | TINYINT | recommended Body level and kind, stored only above the required one, else NULL |
| `eq_rad_suit`, `eq_rad_suit_rec` | TINYINT 0 to 3 | required and recommended radiation suit step: 1 shielded garment, 2 storm vest, 3 hardened suit |
| `eq_mobility` | TINYINT 0 to 2 | powered frame, same encoding |
| `eq_dust` | TINYINT 0 to 2 | dust filter or seals, same encoding |
| `eq_gas_filter` | TINYINT 0 to 3 | required gas filter grade |
| `eq_co2_scrub` | TINYINT 0 to 2 | CO2 scrubber on the intake, same encoding |
| `eq_supply`, `eq_supply_rec` | TINYINT 0 to 4 | required and recommended Air Supply level |
| `eq_shelter`, `eq_shelter_rec` | TINYINT 0 to 2 | required and recommended shelter: 1 shielded shelter, 2 hardened habitat |
| `shelter_column_g_cm2`, `outside_hours` | DOUBLE | shielding a shelter needs; hours outside a year the listed suit allows (NULL when unlimited) |
| `band1_mean_k` to `band7_mean_k` and the other `bandN_*` columns, `temp_mean_swing_k` | DOUBLE, K, h, d and a state word | the annual means, the day and night swing, the hours of day and night, the longest polar day and night and the state of each of seven bands of latitude (5.5, GEN.1108); replaced `temp_polar_k`, `temp_north_k`, `temp_equator_k`, `temp_south_k` and `temp_swing_k` |

Each in its own column so the search can filter it (the rule of section 4 of the index design). New numeric columns go in `db/column_units.py` in the same PR; these are dimensionless levels, like the tier columns they replace. The name and code are computed by `physics/equipment.py` at read time (a pure function of the levels), like `hab_note`.

**API.** `GET /api/planets|moons/<id>` and the system document gain one object, and `equipment_tier` is removed (no backward-compatibility layer):

```json
"equipment": {
  "required": ["pressure suit", "dust-tight seals"],
  "recommended": ["hardened habitat"],
  "code": "B4p F100 S2 [R2]",
  "items": [
    {"slot": "body", "item": "pressure suit", "severity": "required", "domain": "Pressure", "text": "The air pressure is 0.6 kPa, so a pressure suit is required."},
    {"slot": "filter", "item": "dust-tight seals", "severity": "required", "domain": "Air", "text": "Fine regolith dust wrecks open seals."},
    {"slot": "radiation", "item": "hardened habitat", "severity": "recommended", "domain": "Radiation", "text": "The dose is 236 mSv a year, nearly all cosmic rays, which no suit stops; about 250 g/cm2 of regolith (the shower maximum of the cosmic-ray curve) brings a long stay into the Green band."}
  ],
  "regions": {"polar_c": -110, "northern_c": -65, "equatorial_c": -45, "southern_c": -70, "swing_k": 90}
}
```

**Search and lists.** The single "Equipment a Human Needs" facet becomes "Needs" ticks per item (garment, thermal suit, sealed suit, pressure suit, armoured suit, vehicle, shielded garment, storm vest, hardened suit, dust, gas filter, CO2 scrubber, gas-mix mask, own air, mixed gas, shielded shelter, hardened habitat), each with "required" or "required or recommended"; there is no "none" option (decision E), and "worlds that need no gear" is the search with every tick off and "needs nothing" set. The list column "Equipment" shows the required list, with the recommended list in the tooltip; the symbol cluster beside the PHI-6 symbols (UX.101, UX.152) can carry the chips as small icons.

## 8. Relation to the PHI-6 domains and to the old line

### 8.1 Domains

| PHI-6 domain | Slots it drives |
|---|---|
| Pressure | Body (pressure suit under 14.3 kPa, diving suit, armoured suit); Supply (own gas under 14.3, mixed gas over 400, isolated over 10 MPa) |
| Temperature | Body (cold, heat); Supply (breathing air colder than -50 or hotter than 140 C) |
| Air | Filter (dust, gas, CO2); Supply (thin or rich oxygen, a gas past IDLH, CO2 over 5 kPa) |
| Radiation | Radiation Protection (suit steps, shelter); Body for UV |
| Water | nothing worn: water is treated or made from ice, which the describer already says |
| Suitability | nothing worn: a base brings its own food and shelter |

The same condition gives the domain its colour and the slot its level, so I propose a **consistency test** instead of rebuilding the score: a domain at Blue requires nothing from its slots; Green requires at most a level-1 item (garment, gas-mix mask, one filter part); Yellow at most levels 2 and 3 and radiation gear only as recommended; Red may require anything. Declared exceptions, because the gear follows the physics and the scores stay (decision B): a gas between its IDLH and ten times the acute limit (Yellow in the score, own air in the gear); pressure from 14.3 to 20 kPa (Green in the score, pressure suit recommended or required); best-mix density over 10 g/L below 10 MPa (Yellow, armoured suit required). A later item can bring the score bands to the physics if Boss wants them to agree; this draft does not change them.

### 8.2 The old line, mapped and retired

| Old tier | What the new slots say |
|---|---|
| 0 ideal | nothing listed; plus the cases the old line missed (a cold day, a toxic gas it had no row for) |
| 1 breathing mask | Supply 1 alone, or with a garment |
| 2 mask with scrubber | a gas filter or the CO2 scrubber, supply 0 |
| 3 environment suit | any of Body 2 and up, or Supply 2 and up |
| 4 full life support with radiation hardening | a required radiation step: shielded garment, storm vest, hardened suit, shielded shelter or hardened habitat, whichever the dose needs (and whatever else the world needs) |

Removed in the same PR, with no shim: `habitability.equipment()`, `EQUIPMENT_NAMES`, `EQUIPMENT_LABELS`, `EQUIPMENT_RULE` and `EQUIPMENT_MEANINGS`, the stored `equipment_tier`, the `equipment` search facet and its help link, the page chips that name a tier, and `tuning.SPECIES_EQUIPMENT_MIN_TECH` (replaced by a table of the least tech level per slot level, section 9.3).

## 9. Worked examples and what the sample shows

### 9.1 Named worlds, old line against new

Inputs as in the index design's worked examples, or the body's textbook values. Mars and the Moon use their day and night extremes (-125 to +20 C and -173 to +120 C).

| World | Old line | Required | Recommended | Code |
|---|---|---|---|---|
| Earth, 15 C | ideal | (nothing listed) | | |
| Earth, Arctic day at -35 C | ideal | insulated suit | | `B2c` |
| Everest summit (33.7 kPa, -30 C) | breathing mask | insulated suit + gas-mix mask | | `B2c S1` |
| Mars (0.24 Sv/yr, cosmic rays) | environment suit | pressure suit + dust-tight seals | hardened habitat | `B4p F100 S2 [R2]` |
| Moon (0.5 Sv/yr, cosmic rays) | environment suit | pressure suit + dust-tight seals | shielded garment + hardened habitat | `B4p F100 S2 [H1 R2]` |
| Airless moon of a flare star (5 Sv/yr, storm particles) | full life support | pressure suit + dust-tight seals + storm vest | hardened suit + shielded shelter | `B4p F100 S2 H2 [H3 R1]` |
| Close planet of a flare star (100 Sv/yr, storm particles) | full life support | hardened pressure suit + hardened habitat (about 800 hours outside a year) | | `B4p F100 S2 H3 R2` |
| Venus surface (9.2 MPa, 735 K) | environment suit | armoured suit | | `B5h S4` |
| Titan (147 kPa, 94 K) | environment suit | heated sealed suit + self-contained air | | `B3c S2` |
| Europa surface (2,000 Sv/yr) | full life support | hardened pressure suit + hardened habitat (minutes to hours outside) | | `B4p S2 H3 R2` |
| Dense CO2 (300 kPa, 40 C, 95 percent CO2) | environment suit | cooling clothing + self-contained air | | `B1h S2` |
| Hycean (2 MPa, 60 C, H2 and He) | environment suit | liquid-cooled suit + mixed-gas rig | | `B2h S3` |
| Photochemical CO (100 kPa, 5 percent CO) | environment suit | self-contained air | | `S2` |
| Dune world (90 kPa, 45 C) | ideal | cooling clothing | | `B1h` |
| Valley with SO2 and H2S at 10 times the limit | mask with scrubber | full-face gas filter | | `F020` |
| Cold air with 2 percent ammonia | ideal | chemical-barrier suit + self-contained air | | `B3x S2` |
| Air with 50 ppm HCl at 100 kPa | ideal | chemical-barrier suit + self-contained air | | `B3x S2` |
| 2 bar air with 40 percent oxygen | environment suit | gas-mix mask | | `S1` |
| 10 bar air, 50 percent O2 | environment suit | heated diving suit + mixed-gas rig | | `B2d S3` |
| Earth air under 200 times Earth's UV | ideal | full-cover suit | | `B2u` |
| Thin world at 18 kPa, 21 percent O2, 10 C | breathing mask | gas-mix mask | pressure suit | `S1 [B4p]` |

Zero slots are left out of the code as well as the name. Mars and the Moon would also carry regional notes in their text (polar nights).

Where the new line differs from the old on purpose: the photochemical CO world needs air, not a suit (a pure gas hazard, no body hazard); the oxygen-rich 2 bar world needs a mix mask (the 7.3 alternative of the index design, tier 1 diluent mask); and several worlds the old line called "ideal" were missing a row (cold, ammonia, HCl, UV).

### 9.2 What 160 generated systems give (2,261 rocky planets and moons)

| Slot | Distribution |
|---|---|
| Body | 0: 14, garment 1: 20, thermal 2: 59, sealed 3: 59, pressure suit 4: 2,009, armoured 5: 40, vehicle only 6: 60 |
| Filter | dust or seals 2,024 (every body under 10 kPa, the interim rule); gas filter and scrubber 0: a body that needs own gas drops them, and none of the 61 bodies that did not need own gas had a toxic gas or CO2 in range |
| Supply | 0: 27, mask 1: 34, own gas 2: 2,097, mixed gas 3: 3, isolated 4: 100 |
| Radiation, required | nothing 1,025; shielded garment 214; storm vest 84; hardened suit 118; shielded shelter 301 (with a garment 100, storm vest 77, hardened suit 76, hardened suit with an hours limit 36, no suit 12); hardened habitat 519 (503 with a hardened suit and an hours limit) |
| Radiation, recommended | 274 list nothing; of the rest nearly all get the hardened habitat (1,555), with a garment 406, storm vest 241, hardened suit 180 |
| Pressure 14.3 to 20 kPa | 37 bodies move from "mask" to "pressure suit recommended or required" under the physics |
| Old tier | ideal 21, mask 24, environment suit 980, full life support 1,236 |

The old tier 4 (1,236 bodies) is exactly the Red dose case, and all 1,236 are storm-particle worlds (over half their dose from a flaring star): the new ladder splits them across every step from the garment (the soft storms of airless bodies) to the habitat. (The name counts below are from the first prototype run, before the severity split.) Of the 237 bodies with 10 kPa or more of air, the new line names at least 15 distinct combinations (the commonest: heat-shielded suit + self-contained air 45, vehicle only 25, pressure suit 24, liquid-cooled suit + self-contained air 14, armoured suit 14, liquid-cooled suit + gas-mix mask 13, insulated suit 12, gas-mix mask 11, full-cover suit + self-contained air 11, cold-weather clothing 9).

Caveats: the sample's atmospheres contain only the ten stored gases, so the toxic list beyond CO, H2S and SO2 cannot fire in it (it will matter when GEN.424 lands); the dust rule is the interim one; the regions and swing are not in the prototype yet (they only add recommended items and text, since the required gear reads the average); and the thresholds of UV and gravity are DESIGN values. A sweep of the real build, with a table of every combination, is the first thing the build item should produce.

### 9.3 Species and tech level (POP.12)

A species could live somewhere with equipment when its tech level reaches the heaviest slot level the world needs. Replace the five-value tuple by a table (DESIGN, from the old tuple 0, 2, 2.5, 3, 4): Body 1 at 0, 2 at 2, 3 at 3, 4 at 3, 5 at 4, 6 never; dust 2; gas filter 1, 2, 3 at 2.5, 2.5, 3; scrubber 2.5; Supply 1 at 2, 2 at 3, 3 at 3.5, 4 at 4; shielded garment 1, storm vest 2.5, hardened suit 3.5; shielded shelter 2, hardened habitat 4. A species needs the tech for every required item; recommended items do not gate the fit. The slot thresholds read the species' limits (`habitability.bands`) as the old tier did.

## 10. What it touches

**Code.** `physics/habitability.py` (equipment removed), `physics/habitability_world.py` (`equipment_tier` in `SCORE_FIELDS`, `EQUIPMENT_*`), `physics/protection.py` (`requirements`, `kit_sentence`, `SUIT_RANK`, `OVERALL_GEAR`, `SHORT_GEAR` become readers of the needs list), `physics/protection_limits.py` (new rows), `physics/species_limits.py` and `db/species_fit.py` (assess), `db/query.py` (search facet and clause), `web/searchpage.py`, `web/lib/habitability_ref.py` and `web/class_pages.py` (the levels page, UX.90), `web/lib/systempage.py` (equipment chip), `web/templates/habitability.html`, `system.html`, `classes.html`, `cli/generate.py` (its mention), `tuning.py` (`SPECIES_EQUIPMENT_MIN_TECH`), `db/column_units.py`.

**Docs.** `habitability-index.md` sections 2 (Equipment), 7.3 (the five tiers, which this supersedes and implements in part), 8, 9 and 10; `database-schema.md`, `api.md`, `html-interface.md`.

**Tests.** `test_habitability.py`, `test_habitability_world.py`, `test_protection.py`, `test_species_limits.py`, `test_class_pages.py`, `test_web_search.py`, `test_wording_agrees_with_the_index.py`, `test_api.py`, `test_systempage.py`, and the epochs golden only if a draw changes (none does).

**Items.** GEN.89 (score storage), GEN.148 and GEN.207 (dose bands, gates), GEN.208 (the one limits table), GEN.216 and GEN.218 (the describer), UX.90 (levels page), UX.100, UX.101, UX.150 and UX.152 (lists, symbols, filters, score display), GEN.643 (habitable), GEN.649 and GEN.653 (life scores, food), GEN.87 (dose and UV inputs), GEN.414 (haze, dust columns), GEN.424 and DB.71 (all gases stored), POP.12, POP.13 and POP.14 (species limits and placement), MAP.143 and MAP.220 (map colours and markers, which may want the strip). The classes of GEN.645 to GEN.648 do not read equipment. The rogue-planet score (found in the model-gaps report, missing item 3) should carry the same columns when it lands.

**Epoch and migration.** The levels are a pure function of the stored row, so no generated value changes: one Alembic migration (0136) that adds the columns and drops `equipment_tier`; it does not backfill, `planetgen rescore-habitability` (`store.refresh_habitability`) fills the levels of a galaxy made before it, and a body not re-scored shows no equipment line until then. No generator epoch. The three new inputs arrive with their own lanes' epochs and are read when present.

## 11. Proposed items (numbers for TODO to assign)

1. Equipment items: `physics/equipment.py` (needs with severity, levels, the two lists, code), new limit rows (severity bands, the five radiation steps with their shielding and the 2,000-hour working year, toxic table), the radiation ladder computed with the dose functions of `physics/radiation.py`, the consistency test with its declared exceptions, and removal of `habitability.equipment` and the tier constants.
2. Equipment columns: migration with exact backfill, `column_units`, API object (required, recommended, items, regions), removal of `equipment_tier`.
3. Search "Needs" ticks, the levels help page (UX.90) rewritten from the same limit rows, the describer sentences from the needs list (no "none" anywhere).
4. Chips for the lists and system page, solid for required and outlined for recommended (with UX.101 and UX.152).
5. Regional temperatures (Polar, Northern, Equatorial, Southern; day side, terminator, night side for locked worlds) and the day and night swing: module, columns, the regional sentences of the detail text, tests against Earth, the Moon, Mars and Mercury.
6. Toxic-species rows in the limits table, wired to the gas rows of GEN.424; verify the RECALLED limits against NIOSH and ACGIH.
7. Dust loading and dust toxicity input (with GEN.414 and lane K's surface composition).
8. Per-item tech-level table for the species fit (POP.12).

Order: 1 and 2 together (the cut-over), then 3, 4, 5 and 8; 6 and 7 follow the lanes that supply their inputs. No generator epoch: everything is derived from stored values.

## 12. Evidence notes

**Checked on the web this session (V).** NASA OCHMO technical brief: career limit 600 mSv, 250 mSv per solar event. The ISS spacesuit at 4.3 psi (29.6 kPa) of pure oxygen and the 30 kPa target of mechanical counterpressure suits (npj Microgravity 2021; MIT BioSuit). The wet-bulb limit of 30.55 plus or minus 0.98 C (Vecellio et al. 2022, Penn State). NIOSH's 2025 IDLH for hydrogen chloride, 45 ppm. Hopcalite escape respirators exist and last about two hours (MSA W65 sheet). AstroRad: high-density polyethylene, 26 kg, flown on Artemis I, an estimated 40 to 60 percent effective-dose cut for great storms, no help against cosmic rays (Science News, Space.com, StemRad). Lead aprons: 0.5 mm lead equivalent passes about 5 percent at 100 kVp, 0.25 mm about 10 percent (PMC 6038217 and apron studies). Shuttle-era suit helmet and visor layers 0.01 to 0.38 g/cm2; the cooling tubes cover 40 percent of the body (NASA NTRS 20040065979 series). Demron: at least 50 percent of gamma up to 130 keV (vendor claim).

**From Boss's research notes in `docs/design/` (N), not re-checked.** The pressure, hyperbaric, thermal, mask-limit and respirator-hierarchy tables of `Environmental.Extremes.and.Protection.Research.md` and `Thermal.Respiratory.Protection.Research.md`; the gas limits of `Atmospheric Toxicity.md`; the radiation figures of `Naturally Occuring Ionizing Radiation.md` and `activity-magnetism-radiation-hydrosphere.md`.

**Recalled (R), check at build time.** The assigned protection factors (10, 50, 25 to 1,000) and the 19.5 percent oxygen rule of OSHA 29 CFR 1910.134 (the web search could not show the numbers); the long-term limits and IDLH of NH3, HF, Cl2, NO2, HCN and O3; the 100 ppm IDLH of H2S (the built table uses 500 ppm for acute); lead's mass attenuation at 662 keV (about 0.11 cm2/g); the StemRad 360 Gamma marrow figure.

**Design values (D).** The UV thresholds 3 and 30; gravity 1.5 and 3; armoured suit above 400 C and vehicle above 500 C; the interim dust rule under 10 kPa; the shielding of each radiation step (1, 1 + 5 torso, 5 + 5 torso, 50 and 300 g/cm2), the 60 percent torso share and the 2,000-hour working year; the severity band for gravity; the polar heat-transport factor; the tech-level table; the idea that an intake scrubber is workable up to 5 kPa of CO2.

**A caveat for the build.** The repo's cosmic-ray curve is log-interpolated from the Moon (0 g/cm2) to Mars (16.4), so it credits a 1 g/cm2 garment with a 17 percent cut on an airless body; measurements show thin layers barely change the cosmic-ray dose. The storm part, which decides every Red world in the sample, is not affected; the Moon's recommended garment does depend on it. The owner of `radiation.py` may want to hold the curve flat below a few g/cm2.

**Out of scope.** Wind and storms, micrometeoroids, lightning, flammable atmospheres (an oxygen-enriched suit in hydrogen or methane air), pathogens and local life, and how long gear lasts (sortie length). Each could become a further slot or a flag; none changes the three Boss named.

## 13. Files

The prototype (`eq_proto.py`, `sweep.py`) lived in the session's scratch folder and is not part of the repo. The build is `physics/equipment.py` with the limit rows of `physics/protection_limits.py`.

## 14. Built notes (GEN.655, DB.93)

- **The hardened suit is a word, not a body level.** A hardened suit is the suit the world asks for, hardened: "hardened pressure suit" replaces "pressure suit" in the list. `eq_body` keeps the level of the suit alone; `eq_rad_suit` = 3 says it is hardened. A recommended hardened suit over a required plain one lists as the hardened name under Recommended.
- **Inspired oxygen has a recommended band too.** Between 8 and 16 kPa of inspired oxygen a gas-mix mask is recommended (`inspired_o2_rec_kpa`); below it, required. The pressure suit is required under `suit_required_kpa` (16 kPa) and recommended under `suit_recommended_kpa` (20 kPa), both scaled for another species by its pressure limit over the human 14.3.
- **Radiation.** A step's shielding runs through `radiation.gcr_dose_msv_yr` and `radiation.sep_dose_msv_yr` at the added column mass, as ratios of the stored dose parts (`dose_gcr_sv_yr`, `dose_sep_sv_yr`, `dose_ground_sv_yr`, filled by `habitability_world`); the ground dose is not reduced. A suit is judged over the whole year when no shelter is needed, else over the working year outside (`work_hours_outside_yr`, 2,000 h); `outside_hours` is what the dose allows. The heliosphere multiplier scales the cosmic-ray part.
- **Species.** `tuning.SPECIES_EQUIPMENT_TECH` gives a tech level per item name; a species can live there with equipment when its tech level reaches the highest of the required items (`equipment.tech_needed`); vehicle only is 99, never within reach.
- **Corrosive liquid.** A surface liquid below pH 2 or above 12 asks for a chemical-barrier suit; that is a body hazard, not a water score, and the PHI-6 water score is unchanged.
- **Migration.** Each table is one `ALGORITHM=COPY` statement (`planets` and `moons` have a FULLTEXT index, which refuses an instant column drop).
- **Tests.** `tests/test_equipment.py`, `tests/test_migration_survival_equipment.py`.

- **Toxic rows verified (GEN.657).** Each of the ten rows of `protection_limits.TOXIC_GASES` was read from its NIOSH Pocket Guide page (cdc.gov/niosh/npg, 2026-10-29): CO REL TWA 35 ppm, ceiling 200, OSHA PEL 50, IDLH 1,200; H2S REL ceiling 10, IDLH 100; SO2 REL TWA 2, ST 5, OSHA PEL 5, IDLH 100; NH3 REL TWA 25, STEL 35, IDLH 300; HF REL TWA 3, ceiling 6, IDLH 30; HCl ceiling 5, IDLH 50 in the pocket guide and 45 in NIOSH 2025-110 (the newest, used); Cl2 REL ceiling 0.5, OSHA ceiling 1 (the more permissive, used), IDLH 10; NO2 REL ST 1, OSHA ceiling 5, IDLH 13 (NIOSH 2017-202, the pocket guide's pointer; the 1994 value was 20, which the draft recalled); HCN REL ST 4.7, IDLH 50; O3 REL ceiling 0.1, IDLH 5. ACGIH values are not freely published and were not consulted. Only the NO2 IDLH changed from the draft, and no world reads NO2 yet. `tests/test_equipment.py` holds the published IDLH of every row and checks, for every gas and a ladder of ratios to its long-term limit, that the filter grade is 0 to 3 by the ratio (10, 50, 1,000), that past the IDLH or 1,000 times the limit the grade drops out for own air (Supply 2), and that a corrosive gas past its IDLH asks for the chemical-barrier suit.
- **The wiring waits for the lock.** `habitability_world.world_from` builds `World.gases` from five of the stored species (`_GAS_NAMES`: O2, CO2, CO, H2S, SO2), so NH3 (75 of 623 planets in a 60-system sample, up to 0.9 kPa, thirty times its IDLH), HCN and the rest never reach the filter and suit rules even though `planet_gases` stores them. Reading them changes stored `eq_*` columns and the rosetta golden, so it is a generation change (epoch, golden) and rides the generation merge lock.

## 15. Boss's follow-up of 11:25 PM EDT 2026-10-23 (GEN.671, UX.163)

- **GEN.671, a hardened suit against cosmic rays.** In the game a hardened suit works against cosmic rays, so section 2's "a suit does almost nothing against cosmic rays" no longer decides the gear. Where the dose is past the Green edge and no ordinary garment, vest or shelter below the habitat brings the year under it (a cosmic-ray world: Mars, the Moon), the **hardened suit is required** and the **hardened habitat is recommended**. Mars reads "Required: hardened pressure suit + dust-tight seals. Recommended: hardened habitat." (`B4p F100 S2 H3 [R2]`); the Moon the same. The recommended shielded garment and the recommended-only habitat of the draft's table are gone. The physics ladder itself is unchanged; the rule is in `equipment._radiation`.
- **UX.163, a Required line and a Recommended line.** A scored planet or moon gives a Required line and a Recommended line (`equipment.lines`) at the head of its PHI-6 rating in the details, the text and the API (`required_line`, `recommended_line`). A line with nothing in it is left off, never shown as "none" (Boss, 2026-10-24), so Earth has neither and the search column and `line` agree. Decision E (never list an item that is not needed) is unchanged.
- **Per-item specs in the details.** The body's description carries a "Survival equipment" facts table: one row per item with how badly it is needed, the PHI-6 category it comes from and the tech level that makes it, then the shelter's shielding (g/cm2) and the hours outside a year the suit allows.
- A galaxy rescored before this change keeps the old radiation levels until `planetgen rescore-habitability` is run again.

## 16. Which row an item sits under (GEN.702)

An item carries every PHI-6 domain it can come from (`Item.domains`, the likeliest first; `domain` is the first), because the stored levels do not record the cause: an air supply of level 1 or 2 may be for the air, the pressure or the cold; a hardened suit is for the suit's own stress and the radiation; a body level of 4 or more (a pressure suit, an armoured suit) also keeps out the pressure, the temperature and the air, and a sealed suit the air. A category's row lists the items whose domains include it and whose own colour is not Blue (`Item.rows`), the first domain when all of them are Blue. The API's items gain `domains`.

## 17. Hours outside in a day (GEN.768, Boss 2026-10-24)

The property sheet used to label the yearly `outside_hours` "Hours of daylight a person can be outside", which it is not: it is the
hours of a working year that the listed radiation suit allows before the dose limit, and `None` unless that is a limit. Boss asked
for the hours of daylight a person can be outside in one local day, with the daylight hours of the planet beside it. The first
part is built (`physics/outside_hours.py`) and the answer is stored (schema v187): `outside_hours_day`, the hours of one local day, and
`outside_hours_limit`, the name of the limit that sets them. `habitability_world.compute` works both out with the rest of the score
from the equipment levels, the stored `daylight_hours` (GEN.767) and the day and night swing, so a dose refresh or a re-score keeps them
current. A world whose daylight (or, for a locked world, terminator belt temperature) is not stored has neither (NULL). The property
sheet shows the day's hours in the glance and in the rotation group with what limits them in plain words (`outside_hours.LIMIT_WORDS`);
the yearly dose-limited hours keep their own row, "Hours a person can work outside per year (dose limit)".

**The answer is the smallest of five limits and the one that binds is named** (`outside_hours.LIMITS`: daylight, dose, heat,
cold, air, uv, and vehicle for a body that needs one, which has 0 hours):

| Limit | What it is |
|---|---|
| daylight | the daylight hours of the local day (GEN.767; a locked world's terminator belt, or none) |
| dose | the yearly hours the listed suit allows (`outside_hours`) over the work days of the year (`work_days_per_year`: the local solar days, 24-hour shifts for a locked world); no limit when the yearly hours are `None` |
| heat, cold | the time before the core temperature leaves 35 to 38.5 C (`exposure_hours`) |
| air | the gas of the listed Air Supply against a shift (`protection_limits.supply_hours`): a charge lasts a share of the 8-hour shift (mask 1, own gas 0.5, mixed gas 0.75, isolated air 1); ambient air has none |
| uv | the day's ultraviolet against what the worn garment is rated for: the band edge it stands for (bare 3, sun-cover clothing 30, full-cover suit 300 times Earth's, a sealed suit any), the dose linear in the index, so past the edge the hours are the daylight times edge over index |

**Heat and cold** are a one-node balance of the 70 kg, 1.8 m2 reference body at a light to moderate working rate (150 W): its heat
capacity (3470 J/(kg K)) against the metabolism, the dry flow through the garment's insulation and the surface film (radiation plus
convection that goes as the square root of the gas density), the evaporation of sweat (ISO 7933: the vapour-pressure difference
between the wet skin and the air over the clothing's and the film's vapour resistance, at most 500 W, none in a sealed suit) and the
garment's cooling or heating power (a cooling vest, a liquid cooling garment, a suit's thermal control, a heated suit;
`protection_limits.GARMENTS`). The garment for the hot end is the one the ladder lists at the day's hottest temperature and for the
cold end the one at the coldest (`equipment.body_at`), so a person wears the clothing the ladder lists for the weather they
are in; an unsealed garment is worn open (30 percent of its insulation) in air colder than the skin. The day's temperatures are the
world's cold and hot extremes (`habitability.temperature_extremes_c`) widened by half the day-night swing of the equator band (`band4_swing_k`).

**What it shows.** The ladder's gear is chosen to cope, so heat and cold bind only near the hot and cold ends of a garment's band (a
liquid-cooled suit at 120 C with a 40 K swing: about 6 h; an insulated suit at -49 C: about 11 h), the air supply binds in a sealed
suit (4 h on own gas in a vacuum or the deep cold), the dose binds on a world that needs a shelter, and the ultraviolet binds only
past a garment's rating. Earth has the full daylight.

Every number is a `protection_limits` row or table with a source tier and a tweak level in section 6.8 of
[gap-filling-plan.md](../research/gap-filling-plan.md) (the rows OH1 to OH10; list in
[stand-in-values.md](../research/stand-in-values.md)). Nothing is drawn.
