# The composition of a surface liquid (GEN.1217 design)

Status: Boss answered the six questions of section 9 on 2026-10-31 (11:56 and 11:57 AM EDT); the answers are in section 9 and the design follows them. GEN.1215 (the acid and ammonia columns) is built as a branch; the rest (GEN.1218) is not. It follows the research of GEN.1216 (`/mnt/project-files/reports/gen1216-surface-liquid-composition-2026-10-31.md`, called "the report" below) and the census of GEN.1214 (`/mnt/project-files/reports/tides1/gen1214/census-report-2026-10-31.md`). It decides what a body stores about its liquid, so that GEN.1215 (the acid and ammonia columns) and GEN.1218 (the build) fit one plan.

## 1. Where the code is today

`physics/surface_liquids.py` holds a table of liquids (one row each: window, density, viscosity, tier and tweak, `surface_liquids.json`). `physics/surface_liquid_world.py` (GEN.1035) gives a body at most one `surface_liquid` id with a `surface_liquid_fraction` of the surface and a mean depth. Water takes the id when the hydrosphere is a liquid; otherwise lava, otherwise the volatile with the best score among those stable at the surface and supplied (air at saturation, or the body's ice).

What is not stored, and is worked out each time the life classes ask:

- **Acid.** `acid_medium.level` (GEN.888) works the strength of sulfuric acid out from the air's SO2, an oxidiser, and the water activity. Nothing is stored; a pool and a cloud are one number, 0, 1 or 2.
- **Ammonia.** `strange_solvents.ammonia_water_mass_fraction` (GEN.894) is `nh3_share_of_ice * ice / water`, the ammonia as a share of the water's mass. The Class 8 gate is that value at 0.1 or more. It is not a property of any liquid the body holds.
- **Water chemistry.** `ocean_class`, `ocean_ph`, `water_activity` are stored for a water ocean (GEN.88), with no salinity and no ion mix.
- **Hydrocarbon seas.** The liquid id is a single volatile (`methane`), though Titan's seas are methane, ethane and nitrogen together. `surface_liquids.ideal_sea` gives a Raoult mix with a 0.3 nitrogen factor, but the world step does not use it.

## 2. What the research gives (and does not)

The report's headline: there are anchor numbers for every family and no published law that turns a planet's conditions into a composition. The generator therefore builds each law from saturation, solubility and phase equilibrium, with the anchors as checks, and the results are stand-ins of tier T6 to T9 in the registry (GEN.725). Physics drives; no game rates.

| Family | What can be a law | What is only an anchor |
|---|---|---|
| Water seas | Freezing point and activity from the salt (eutectics: NaCl 251.9 K, MgSO4 269.6 K, Mg and Ca perchlorate about 205 and 198 K, Na perchlorate 239 K); carbonate-system pH from CO2 pressure and alkalinity | Salinity against age (use modern 35 g/kg, range 20 to 50), ion mix against rock type |
| Sulfuric acid | Strength from the water activity of the air (exists, `acid_medium`); freezing point from the two eutectics (211 K at 36 wt percent, 199 K at 38 wt percent) | A full liquidus; HCl and HNO3 (no table, stays a stand-in) |
| Ammonia-water | Liquid down to the eutectic 176 K at 33 wt percent; seven measured liquidus points 0 to 27 wt percent | Anything above 27 wt percent; viscosity (sources differ 2 to 20 times) |
| Liquid CO2 | Span-Wagner saturation table (216.6 K triple, 304.1 K critical) | Density above 304 K |
| Methane-ethane-nitrogen | Raoult or regular-solution equilibrium with the N2 solubility table (Kumar and Chevrier) | Model spread of a factor of 10 in methane; the observed seas (methane 71 to 80 mol percent) |
| Sulfur | Window 392 to about 450 K, viscosity peak at 432 to 458 K | Liquid density |

Life: water has numbers (pH -0.06 to 12.5, water activity down to 0.605). Acid, ammonia-water, dense CO2 and hydrocarbon life have none, so the gates of Classes 5 to 9 stay at tier T8 or T9, tweak L3 or L4, and the registry says so.

## 3. The model: a main liquid and what is dissolved in it

A body's surface holds a **main liquid** (what the existing columns say) and a **mix**: the mass or mole fractions of the few things dissolved in it that something else reads. The mix is stored as typed columns, one per quantity, because the repo gives every numeric column a unit (`db/column_units.py`) and every column a sheet row and a rule.

Stored quantities (all nullable, all NULL when the body has no such liquid):

| Column | Unit | Meaning | Read by |
|---|---|---|---|
| `surface_liquid_acid_wt`, `surface_liquid_acid_layer` | percent of mass, text (`pool` or `cloud`) | sulfuric acid in the acid liquid: its strength, and whether it is a pool at the ground or the cloud level (built, GEN.1215) | Classes 5 and 9 gate, the sheet |
| `surface_liquid_ammonia_wt` | percent of mass | ammonia in the water-ammonia liquid (built, GEN.1215, for a body with a liquid water ocean) | Class 8 gate, the sheet, the freezing point of the sea |
| `surface_liquid_salinity_g_kg` | g/kg | dissolved salt of a water sea | the sheet, the freezing point and water activity of the sea |
| `surface_liquid_methane_mol`, `_ethane_mol`, `_nitrogen_mol` | mole fraction | the hydrocarbon sea's makeup | Class 7, the sheet, the density of the sea |

GEN.1215 stores the first two (the acid and ammonia fractions) and makes the life classes read them instead of recomputing; GEN.1218 adds the other four and the liquids of section 4. Naming them all here means the GEN.1215 migration does not have to be redone. The layer column keeps the cloud case apart from the pool case; both are stored (question 3). Built: `surface_liquid_world.generate` stores the three columns, and the life classes read them through `surface_liquid_world.acid_level` and `ammonia_water_mass_fraction`, falling back to working them out for a row without them. The Class 8 gate line is still a tenth of the water's mass (`life_classes._solvent_term`); question 4's answer moves it to about 5 percent of the liquid's mass, which is Aux 2's habitability change, tuned by the GEN.1225 research.

## 4. Liquids and supply to add

1. **Ethane** as a liquid and as an ice (today it hides in the unnamed "other" ice). Titan's lakes need it.
2. **Nitrogen ice** of cold bodies. The census found the comet volatile table has no N2, so every body's `ice_n2_mass_fraction` is 0 and Triton and Pluto kinds cannot hold N2 ice. The research (GEN.1225, `/mnt/project-files/reports/gen1225-ammonia-gate-and-nitrogen-ice-2026-10-31.md`, section 3.5) gives the rule: (a) a body whose ice formed in a zone colder than the N2 condensation temperature gets N2 ice, and that temperature is 22 to 30 K at the disc's partial pressure (28 K at 1e-4 bar), added to `ICE_CONDENSATION_K` (CS1, T1 L0 for the vapour-pressure curve, T4 for the partial-pressure rule); (b) where it forms, the N2 mass per water mass is the solar ceiling (0.107 to 0.21) times a trapping efficiency of 0.001 to 0.013, which is 1e-4 to 1.4e-3 of the water mass with a middle of about 4e-4 (new row GN1, T3 L2, from comet 67P and Pluto; the two 67P N2/CO values disagree by 5 times); (c) the body keeps it only if its Jeans parameter allows (new row GN2, retention, T8 L4: the sources give the regime but no radius cutoff, so it ties to the escape machinery the repo has). Warmer zones get no N2 ice (their nitrogen is the NH3 of the water ice, CS2).
3. **Sulfuric acid solution** as a liquid with a strength (not one fixed 75 percent row) and a freezing curve from the two eutectics.
4. **Ammonia-water solution** as a liquid whose melting point follows the stored ammonia fraction (176 K at 33 wt percent, 267.7 K at 5.2 wt percent, from the seven measured points), replacing pure ammonia at 195 K as the only ammonia liquid.
5. A water sea gains a **salt** that sets its freezing point (NaCl, MgSO4, perchlorates as the brine rows already do).

## 5. How a body gets its liquid

The order stays: water if the hydrosphere is liquid, lava if molten, otherwise the volatiles. What changes:

1. **Windows follow the mix.** The window of a solution is its own (the ammonia-water liquid is liquid from 176 K up; sulfuric acid at 36 wt percent from 211 K), not the pure liquid's. This could let more cold bodies hold a liquid without any invented supply: the 310 ammonia ice holders the census calls "too cold for pure ammonia" would be in the window of a water-ammonia sea at 176 to 270 K where they also hold the water for it (not yet counted).
2. **A second liquid.** If two liquids are stable and supplied (28 of 37 CO2-window bodies in the census show water only), the one with the larger cover is the main liquid and the other is kept as a **second liquid**: two more columns (`surface_liquid_2`, `surface_liquid_2_fraction`) and no more. See question 1.
3. **The pick** among volatiles stays by supplied amount, but a body cold enough for methane to be liquid (above 90.7 K) keeps methane over carbon monoxide; the census found 26 CO seas at 70 to 90 K that a Titan-like pick would not choose. This is a tie-break rule, not a new number.
4. **Hydrocarbon mix.** The mole fractions come from Raoult with the N2 solubility table (a table, not the single 0.3 factor), at the sea's temperature and the air's methane partial pressure; ethane from the sea's accumulated ethane (a supply, section 4 item 1). Checks against the report's observed seas: Ligeia 71 methane, 17 nitrogen, 12 ethane (mol percent), Punga 80, 20, 0.
5. **Acid.** The strength and the place (pool or cloud) now come from the same `acid_medium` solver but are stored; the sheet shows them.

## 6. What the life classes read

Class 5 and 9: `surface_liquid_acid_wt` at or over 81 (pool) or the cloud strength; Class 8: `surface_liquid_ammonia_wt` at or above the liquidus share for the body's temperature (about 5 percent at 267 K, 14 percent at 251 K, 20 percent at 235 K, 27 percent at 207 K, from the seven measured points of [R0794]) and at most the 33 percent peritectic share, with a 3 and a 5 percent shoulder (registry LF44 rewritten, T3 L1; GEN.1225 section 2.5). The research found no growth or activity above about 0.9 percent ammonia and no source for a 10 percent requirement, so the class is a speculation at T9 L4 in what it says about life, and the gate is a statement about the liquid being liquid; Class 7: a hydrocarbon sea with an ethane or methane mole fraction; Class 6: unchanged (liquid or supercritical CO2, `co2_phase`). Water classes: pH, activity and salinity inside the sourced limits. The gates are registered as stand-ins (T8 to T9, L3 to L4) because no source gives a growth limit; the registry rows say what paper would replace each (the report's upload list).

## 7. Build plan

- **GEN.1215** (columns, lock): migration for `surface_liquid_acid_wt` and `surface_liquid_ammonia_wt` on `planets` and `moons`; unit registry; sheet rows and rules; `acid_medium` and `strange_solvents` write them in `surface_liquid_world.generate`; the life classes read the stored values; the ammonia definition changes from a share of the water's mass to a share of the liquid's mass (question 4). Shares one epoch with GEN.1214 and GEN.1222.
- **GEN.1218a** (liquids and supply, lock): items 1, 3, 4 of section 4 and the windows of section 5.1, then the N2 ice of item 2.
- **GEN.1218b** (mix, lock): the hydrocarbon mix columns, salinity, the second liquid.
- Each part adds registry rows (stand-ins with tier and tweak), a wiring-register line and a census run (this census script, `census_surface_liquid.py`, re-run on the same seeds).

## 8. Checks

After each part: the census counts methane, CO, CO2, ammonia, nitrogen seas and the water-under-150-K count; the acid and ammonia columns must be NULL on every body with no such liquid and the sheet rule must say why; Ligeia and Punga-like seas inside the report's ranges; the nine-class counts do not fall for the water classes.

## 9. Questions for Boss, and his answers (11:56 and 11:57 AM EDT 2026-10-31)

1. A second liquid beside the main one: **yes, one extra liquid.**
2. The makeup of a sea as separate columns: **yes, separate columns.**
3. Sulfuric-acid clouds: **store them too, the strongest level, with a pool-or-cloud flag** (`surface_liquid_acid_layer`).
4. The Class 8 ammonia gate: **near 5 percent of the liquid's mass as a starting guess, but tuned from real physics and chemistry** (Sources 1, GEN.1225). Tuned result: the gate follows the liquidus share for the body's temperature, with the 33 percent peritectic as the upper bound (section 6).
5. Nitrogen ice on far cold bodies: **invent its amount from where nitrogen freezes, as a registered stand-in, also tuned from physics** (GEN.1225). Tuned result: section 4 item 2.
6. Methane and carbon monoxide ice together, methane able to be liquid: **a methane sea.**
