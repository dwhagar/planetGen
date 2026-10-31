# The composition of a surface liquid (GEN.1217 design)

Status: draft for Boss's answers to the questions in section 9. Nothing here is built. It follows the research of GEN.1216 (`/mnt/project-files/reports/gen1216-surface-liquid-composition-2026-10-31.md`, called "the report" below) and the census of GEN.1214 (`/mnt/project-files/reports/tides1/gen1214/census-report-2026-10-31.md`). It decides what a body stores about its liquid, so that GEN.1215 (the acid and ammonia columns) and GEN.1218 (the build) fit one plan.

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
| `surface_liquid_acid_wt` | percent of mass | sulfuric acid in the acid liquid (a pool at the ground, or the strength of the cloud level) | Classes 5 and 9 gate, the sheet |
| `surface_liquid_ammonia_wt` | percent of mass | ammonia in the water-ammonia liquid | Class 8 gate, the sheet, the freezing point of the sea |
| `surface_liquid_salinity_g_kg` | g/kg | dissolved salt of a water sea | the sheet, the freezing point and water activity of the sea |
| `surface_liquid_methane_mol`, `_ethane_mol`, `_nitrogen_mol` | mole fraction | the hydrocarbon sea's makeup | Class 7, the sheet, the density of the sea |

GEN.1215 stores the first two (the acid and ammonia fractions) and makes the life classes read them instead of recomputing; GEN.1218 adds the other four and the liquids of section 4. Naming them all here means the GEN.1215 migration does not have to be redone. The acid column keeps the cloud case apart from the pool case by the existing `acid_medium` level, which stays a derived number from the stored acid fraction and a stored place (see question 3).

## 4. Liquids and supply to add

1. **Ethane** as a liquid and as an ice (today it hides in the unnamed "other" ice). Titan's lakes need it.
2. **Nitrogen ice** of cold bodies. The census found the comet volatile table has no N2, so every body's `ice_n2_mass_fraction` is 0 and Triton and Pluto kinds cannot hold N2 ice. The share must come from the N2 snow line physics (N2 condenses far out, about 25 K), not the comet table; a stand-in with tier "from physics" and tweak L3 until a source gives the abundance.
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

Class 5 and 9: `surface_liquid_acid_wt` at or over 81 (pool) or the cloud strength; Class 8: `surface_liquid_ammonia_wt` over a line (see question 4); Class 7: a hydrocarbon sea with an ethane or methane mole fraction; Class 6: unchanged (liquid or supercritical CO2, `co2_phase`). Water classes: pH, activity and salinity inside the sourced limits. The gates are registered as stand-ins (T8 to T9, L3 to L4) because no source gives a growth limit; the registry rows say what paper would replace each (the report's upload list).

## 7. Build plan

- **GEN.1215** (columns, lock): migration for `surface_liquid_acid_wt` and `surface_liquid_ammonia_wt` on `planets` and `moons`; unit registry; sheet rows and rules; `acid_medium` and `strange_solvents` write them in `surface_liquid_world.generate`; the life classes read the stored values; the ammonia definition changes from a share of the water's mass to a share of the liquid's mass (question 4). Shares one epoch with GEN.1214 and GEN.1222.
- **GEN.1218a** (liquids and supply, lock): items 1, 3, 4 of section 4 and the windows of section 5.1, then the N2 ice of item 2.
- **GEN.1218b** (mix, lock): the hydrocarbon mix columns, salinity, the second liquid.
- Each part adds registry rows (stand-ins with tier and tweak), a wiring-register line and a census run (this census script, `census_surface_liquid.py`, re-run on the same seeds).

## 8. Checks

After each part: the census counts methane, CO, CO2, ammonia, nitrogen seas and the water-under-150-K count; the acid and ammonia columns must be NULL on every body with no such liquid and the sheet rule must say why; Ligeia and Punga-like seas inside the report's ranges; the nine-class counts do not fall for the water classes.

## 9. Questions for Boss

1. Should a body be able to hold a second liquid beside its main one (for example a water ocean with a carbon dioxide pool, or a methane sea with an ammonia-water lake)? Recommended: yes, one second liquid. Alternative: keep one liquid per body.
2. Should the makeup of a sea (acid strength, ammonia share, salt, methane-ethane-nitrogen shares) be stored as separate columns, or as one text column? Recommended: separate columns, as everywhere else in the database.
3. Should the sulfuric-acid clouds of a thick-aired world also be stored, or only acid pools at the ground? Recommended: store the strongest acid level (pool if there is one, else the cloud) with a flag for which it is.
4. Class 8 life today needs ammonia at a tenth of the water's mass. The research says an ammonia-water liquid is a third ammonia at its coldest, and the one life measurement we have (methane-making microbes) stops at one to one and a half percent. Should the gate be set at a share of the liquid (recommended: 5 percent of the liquid's mass, marked as a guess with no data), or kept as it is?
5. Nitrogen ice on far cold bodies has no source in our tables. May we invent its amount from the physics of where nitrogen freezes, as a registered stand-in that research can replace later (recommended), or leave nitrogen out?
6. Should a cold body that holds both methane and carbon monoxide ice, and is warm enough for methane to be liquid, get a methane sea rather than a carbon monoxide one (recommended)?
