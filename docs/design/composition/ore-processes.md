# Ore-forming processes for rocky planet and moon crusts (GEN.1379)

Status: design (docs only, lock-free). Written 2026-11-07 for the real-material build (Boss, 5:24 PM EST on 2026-11-07: all non-fictional generation in the real-material report is built in 8.1; 5:20 PM EST: GEN.798 is widened to rocky crusts). This is N10 part 1 of `/mnt/project-files/research/real-material-8.1-build-plan-2026-11-07.md`. It decides, for the 19 deposit types of `/mnt/project-files/research/real-material-data/ore_deposit_types_and_parameters.csv`, which stored state lets a body have the deposit, how much and how rich, and which fields 8.1 stores so that the 8.2 resource scan (GEN.798, GEN.799) can place and score deposits without more physics. It builds nothing; the build items are GEN.1380 to GEN.1383 (the mechanism slices) and GEN.1378 (the sedimentary and biogenic deposits).

**Source rule.** Every value below that came from the data file is **RECALLED** (the file names its sources and says none was opened). The builder of each slice opens the primary source, or an accepted secondary source, of every value it uses before the value enters GEN.725 (item rule); until then the value is a stand-in at the tier written here. Numbers I added are marked **DESIGN** (a modelling choice, tweak level L0 to L3 as stated).

## 1. Scope and what is not decided here

- **Bodies.** Rocky planets and moons (type `t`, with a crust). Giants and icy-shell bodies have no ore crust; belts, fields and comets are scored by GEN.798 from their composition (GEN.473, GEN.477), and the one extraterrestrial type in the file (metal cores exposed by collisions) belongs to that path, not to this one.
- **What is stored (8.1).** The *potential* of a body: which deposit types it admits, how much ore each holds in total, and how rich. Not individual deposits. A deposit has a place only when the 8.2 scan asks for one (GEN.799), from the stored potential and the object's own seed, so the 8.1 census does not pay for map points nobody looks at.
- **Not decided.** The resource score and the recovery difficulty (GEN.798); the table that stores locations (DB.106); the fictional materials (8.2).

## 2. The generator state the ores read

What the generator stores (checked in the code on 2026-11-07), what it will store once the upstream real-material items land, and what this design adds.

| State | Field | Status | Used by |
|---|---|---|---|
| Tectonic regime | `lid_regime` (`mobile` or `stagnant`, GEN.721) | stored | every hydrothermal and orogenic type |
| Mantle redox | `mantle_redox` (reduced, intermediate, oxidized) and `mantle_delta_iw` | stored | porphyry (oxidised), chromitite and diamond (reduced), sulfide saturation |
| Water | `hydrosphere`, `water_mass_fraction`, ocean cover and state (liquid, ice shell, none) | stored | VMS, SEDEX, weathering, sedimentary types, all hydrothermal |
| Heat | `internal_heat_flux_w_m2`, `urey_ratio` | stored | magma volume, hydrothermal vigour |
| Age | star age; body age (the planet's formed at the star's) | stored | age gates (BIF, unconformity U, the share of a type's lifetime already run) |
| Climate | surface temperature, pressure, liquid-water presence, evaporation (GEN.1319 brine) | stored | laterite and bauxite, evaporites, placers |
| Life history | the GEN.750 timeline and `life_stage` (oxygenic photosynthesis time, ocean oxygenation) | stored | BIF, phosphorite |
| Bulk and mantle oxides | `rock_mass_fraction`, `metal_mass_fraction`, FeO and Al2O3; Cr2O3, TiO2, P2O5, Na2O, K2O, MnO, NiO (N5) | partly; N5 completes it | magmatic, pegmatite, carbonatite, phosphorite |
| Star element set | Ti, V, Cr, Mn, Co, Cu, Zn, P, Li (N1); REE, PGE, Au, Pb and the siderophile split (N2) | planned | every grade ceiling |
| Mantle carbon and water | carbonates, diamond or graphite by redox, hydrous phases (N7) | planned | carbonatite, kimberlite |
| Salt sequence | evaporite sequence from the brine (N9) | built (GEN.1366, `physics/evaporite_sequence.py`) | evaporite |
| **Crust thickness** | `crust_thickness_km` | **new, derived** (section 4.1) | porphyry, orogenic, pegmatite, craton |
| **Craton** | `craton_fraction` (share of the crust that is thick, old and cold) | **new, derived** (section 4.2) | kimberlite, unconformity U, orogenic |
| **Sulfur supply** | `crust_sulfur_index` (0 to 1) | **new, derived** (section 4.3) | magmatic sulfide, VMS, SEDEX |

Planned fields are inputs from other real-material items; a slice that needs one that has not merged waits for it (the build plan orders them: GEN.1380 needs only N5 and N8).

## 3. Three rules that decide whether a type exists

These follow the data file's own grouping and are the gate of every row in section 5.

1. **Magmatic** (chromitite, Ni-Cu-PGE sulfide, Fe-Ti-V oxide): needs a large volume of mafic or ultramafic magma, and for the sulfide type a crustal sulfur source. They exist with or without plate tectonics: a stagnant-lid body gets them from plume or impact melt, so the gate is magma volume (heat flow and mass), not the lid.
2. **Hydrothermal** (porphyry, VMS, orogenic, epithermal, SEDEX and MVT): need liquid water in the crust and a heat engine. Porphyry, orogenic and epithermal need a **mobile lid** (arcs, collision, shear zones); VMS needs seafloor spreading, so a mobile lid *and* an ocean; SEDEX and MVT need a sedimentary basin with brines, so liquid water and a surface. A stagnant-lid or dry body has none.
3. **Sedimentary and weathering** (BIF, phosphorite, evaporite, laterite and bauxite, placer, roll-front uranium): need surface liquid water and time; BIF needs an Fe(II) ocean that is partly oxidised (an oxygen source); phosphorite needs a productive sea (life); laterite and bauxite need a warm wet climate and a stable surface.

Pegmatite, carbonatite, kimberlite and unconformity uranium have their own gates (section 5): continental crust and granite production; a carbonated mantle and rifting; a thick old craton with carbon at depth and low oxygen fugacity; an old basin over reduced basement.

## 4. The three derived fields

Each is a **DESIGN** stand-in (registered in GEN.725 by the slice that builds it, tier T6 or T7 until a source is opened, tweak L2 where it is tuned to the Earth, Mars and Venus values). They are deterministic from stored state plus the object's seed through `draw.py`, so a body made again gets the same record.

### 4.1 Crust thickness

`crust_thickness_km` is the mean thickness of the crust. The wanted behaviour: a small hot body that has stayed stagnant has a thick lid crust, a mobile-lid body recycles and has a thin oceanic crust with thicker continents, and a body whose mantle is very hot has a thin crust. Earth is 7 km oceanic and about 35 km continental, Mars about 50 km and Venus about 30 km (all RECALLED, to open before use). The stand-in is a monotone function of the surface gravity and the mantle temperature already computed by the interior model (`interior_heat`), fitted to those three bodies: two parameters, L2. Used by porphyry (depth of the intrusion), pegmatite and orogenic gold (continental crust), and to scale the volume of ore-bearing crust.

### 4.2 Craton fraction

`craton_fraction` (0 to 1) is the share of the crust that is thick, old and cold enough to hold a keel deeper than 150 km (kimberlite) and a Proterozoic basin over basement (unconformity uranium). Gate: mobile lid, body age at least 2 Gyr, and enough continental crust (a non-zero share of felsic crust, which `surface_composition` already tracks as a rock type); zero otherwise. Magnitude: a DESIGN fraction around 0.3 for an Earth-like body, scaled by age (a younger body has had less time to stabilise keels), drawn once per body.

### 4.3 Crustal sulfur index

`crust_sulfur_index` (0 to 1) says how much sulfur the crust offers a rising mafic magma to saturate on. The data file names the sources of the sulfur: sulfide-bearing sedimentary rocks (needs a surface and an ocean with sulfate reduction) and evaporites (sulfate). The stand-in: a function of the star's S/Fe (an element the star set holds), the water inventory, and whether sedimentary sulfide has had time to build (age). A body with no liquid water has the index of its primordial sulfur only (low); an Earth-like body has a high one. Impact melt sheets (Sudbury) are the other sulfur route: a body with a recent large impact in its record can have the sulfide type without sedimentary sulfur (uses the collision count the planet architecture already records).

## 5. The 19 deposit types

Columns: the **gate** is a condition on the stored state (an AND of its parts; "mobile" and "stagnant" are `lid_regime`); the grade and tonnage ranges are the data file's (RECALLED); the **slice** is the build item that stores the type. Sources named are the data file's, unopened.

| # | Type | Process | Gate on the generator state | Grade (RECALLED) | Tonnage (RECALLED) | Key minerals | Slice |
|---|---|---|---|---|---|---|---|
| 1 | Porphyry Cu-Mo-Au | magmatic-hydrothermal | mobile, liquid water, `mantle_redox` oxidised (fO2 at or above FMQ+1 to +2), crust thickness enough for a 1 to 4 km intrusion, star Cu | Cu 0.3 to 1.0 percent, Au 0.1 to 1 g/t, Mo 0.01 to 0.05 percent | 1e8 to 1e10 t ore | chalcopyrite, bornite, molybdenite, pyrite | GEN.1381 |
| 2 | Volcanogenic massive sulfide | seafloor hydrothermal | mobile, ocean (liquid water), extensional crust, heat flow above a floor | Cu 1 to 4, Zn 3 to 10, Pb 0 to 3 percent, Ag, Au | 1e6 to 1e8 t | pyrite, chalcopyrite, sphalerite, galena | GEN.1381 |
| 3 | SEDEX and MVT Zn-Pb | basinal brine | liquid water, a basin (surface, sediment), brines 100 to 250 C, sulfate reduction (sulfur index) | Zn 5 to 12, Pb 2 to 8 percent | 1e6 to 1e8 t | sphalerite, galena, pyrite, barite, dolomite | GEN.1381 |
| 4 | Orogenic gold | metamorphic fluid | mobile, collision or accretion (continental crust), metamorphism, quartz veins | Au 3 to 10 g/t | 1e6 to 1e8 t ore, giants over 100 t Au | native gold, pyrite, arsenopyrite | GEN.1381 |
| 5 | Epithermal Au-Ag | near-surface hydrothermal | mobile, volcanism, shallow water table | Au 1 to 10 g/t | 1e5 to 1e7 t | electrum, sulfides, quartz, adularia | GEN.1381 |
| 6 | Magmatic Ni-Cu-PGE sulfide | magmatic segregation | large mafic or ultramafic magma volume (plume or impact melt), sulfur index high, PGE from the star | Ni 0.5 to 3, Cu 0.3 to 2 percent, PGE 0.5 to 10 g/t | 1e6 to 1e9 t | pentlandite, pyrrhotite, chalcopyrite, PGM | GEN.1380 |
| 7 | Chromitite | magmatic | large layered mafic intrusion or ophiolite, reduced to moderate fO2, Cr2O3 in the mantle (N5) | stratiform 30 to 45 percent Cr2O3, podiform 40 to 60 | 1e7 to 1e10 t | chromite | GEN.1380 |
| 8 | Fe-Ti-V oxide | magmatic | late Fe-rich liquid in a large layered intrusion, or an anorthosite massif (mostly Proterozoic: an age window), TiO2 and V from the mantle (N5, N1) | magnetite 20 to 30 percent, V2O5 0.3 to 1; ilmenite TiO2 about 32 | 1e8 to 1e9 t | magnetite, ilmenite, titanomagnetite | GEN.1380 |
| 9 | Banded iron formation | chemical sediment | ocean with dissolved Fe(II), a rising oxygen source (photosynthesis or abiotic), age window, no sulfidic ocean | Fe 25 to 35 percent, enriched 55 to 65 | 1e9 to 1e11 t | hematite, magnetite, chert, siderite | GEN.1378 |
| 10 | Laterite Ni and bauxite | weathering | warm wet climate, stable surface, liquid water, millions of years; ultramafic (Ni) or aluminous (Al) source rock | Ni 1 to 2 percent; Al2O3 40 to 60 percent | 1e8 to 1e10 t | garnierite, goethite, gibbsite, boehmite, diaspore | GEN.1382 |
| 11 | LCT and NYF pegmatite | late magmatic | continental crust, granite production, orogenic crust, water, Li, B, F, P | Li2O 1 to 2.5 percent, Ta2O5 0.01 to 0.05 | 1e5 to 1e8 t | spodumene, lepidolite, beryl, columbite-tantalite, cassiterite | GEN.1382 |
| 12 | Carbonatite REE-Nb-P | mantle carbonate melt | rifting or plume, carbonated mantle (mantle carbon, N7), REE and Nb in the star set (N2) | REO 3 to 10 percent, Nb2O5 1 to 3 | 1e7 to 1e9 t | bastnasite, monazite, pyrochlore, apatite | GEN.1383 |
| 13 | Kimberlite diamond | xenocryst delivery | craton (keel over 150 km), carbon in the mantle, reduced fO2, kimberlite magmatism | 0.1 to 1.0 carat per tonne | 1e6 to 1e8 t ore | diamond, garnet, olivine, ilmenite | GEN.1383 |
| 14 | Unconformity-related U | basinal brine, redox boundary | old basin over reduced basement (craton share and age window), oxidising brines | U3O8 1 to 20 percent | 1e4 to 1e6 t U | uraninite, coffinite | GEN.1383 |
| 15 | Roll-front U | groundwater redox front | oxidising groundwater, reductant (organics or sulfide), uranium-bearing volcanic or granitic source | U3O8 0.05 to 0.3 percent | 1e3 to 1e5 t U | uraninite, coffinite | GEN.1381 |
| 16 | Phosphorite | marine sedimentary | productive ocean with upwelling, organic matter (life), phosphorus supply (N18, star P) | P2O5 18 to 35 percent | 1e8 to 1e10 t | francolite, apatite | GEN.1378 |
| 17 | Evaporite | evaporation | arid climate or restricted basin, brine chemistry from weathering inputs (GEN.1319, N9) | halite and potash beds, gypsum | 1e8 to 1e12 t | calcite, gypsum, anhydrite, halite, sylvite, carnallite | GEN.1378 and N9 |
| 18 | Placer | mechanical concentration | running water, long weathering, uplift | ilmenite, zircon, rutile 1 to 10 percent heavy minerals; Au 0.5 to 5 g/t | 1e6 to 1e9 t | ilmenite, zircon, rutile, monazite, cassiterite, gold | GEN.1382 |
| 19 | Asteroid M-type metal | differentiation and impact | not a crust deposit: a disrupted differentiated parent body | Fe-Ni, Ni 5 to 25 percent, PGE tens of ppm | n/a | kamacite, taenite, schreibersite, troilite | GEN.473, GEN.477 and GEN.798 |

Row 19 is listed so the 19 types of the file all have an owner; this design stores nothing for it. Rows 9, 16 and 17 are the sedimentary and biogenic deposits (GEN.1378); this design gives their gate and the fields they store so the same record serves the scan.

## 6. How much, and how rich (the endowment and grade model)

For each type a body admits (its gate holds), 8.1 stores two numbers and one flag.

**Endowment** `E_k`, the total ore mass of type *k* on the body, kg. The model has three factors:

1. **Opportunity**: the area (or volume) of crust that satisfies the type's setting. Hydrothermal and orogenic types scale with the crust of a mobile-lid body, so with surface area times the continental share; magmatic types with the volume of magma the body's heat flow and mass support; weathering and sedimentary types with the area under liquid water or a wet climate. DESIGN, L1.
2. **Time**: the share of the body's life during which the setting existed (an age window, and for types that need life or oxygen, the time since that event on the GEN.750 timeline). Types with a window (BIF, anorthosite Fe-Ti, unconformity uranium) are zero outside it.
3. **A lognormal draw** for the scatter between bodies of the same state. The data file's tonnage range for the type (its 10^a to 10^b) gives the median and the spread of the *per-deposit* tonnage, and the deposit count per unit opportunity gives the sum; the Earth's known totals are the L2 anchor (porphyry copper at roughly half of terrestrial Cu, the file says, to be checked against its source).

**Grade** `g_k`, the mass fraction of the key element or oxide in the ore. Drawn inside the file's range, then **capped by what the host supplies**: the concentration factor over crustal abundance (the file's: Cu about 140, Zn 750, Pb 3000, Ni 320, Cr 3300, U 740, Li 750, Au 2300, Fe 11, Al 3, all RECALLED with the crustal abundances the file cites) multiplies the body's crustal mean of that element, which comes from the star set and the oxides (N1, N2, N5). So a metal-poor star gives poorer ore, a Cr-poor mantle gives a lean chromitite, and the ordering the file states (Fe and Al easy, Cr, Pb and Au hard) comes out of the factors, not out of a table. A type whose capped grade falls below the bottom of its range is not stored (the body has no ore of that type, only the dispersed element).

**Flag** `status`: `present`, `absent_by_state` (a gate failed) or `absent_by_chance` (the gate holds but the endowment draw gave no deposit), so the page and the census can tell a dry planet from an unlucky one.

**Draws.** One stream per body and type keyed on the object ID and the type id (the same pattern as the oddity streams, GEN.1036), so adding a type moves no other draw; the gate and the factors are deterministic from the stored state. Every constant (the opportunity scalings, the age windows, the deposit count per unit opportunity, the lognormal widths) is a GEN.725 stand-in with its tier and tweak level, replaced when its source is opened.

## 7. What 8.1 stores (the fields)

One row per body and admitted type, in a child table (the pattern of `field_members` and the other child tables), plus two columns on the body.

| Where | Field | Unit | Meaning |
|---|---|---|---|
| `planets`, `moons` | `crust_thickness_km` | km | section 4.1 |
| `planets`, `moons` | `craton_fraction` | fraction | section 4.2 |
| `planets`, `moons` | `crust_sulfur_index` | fraction | section 4.3 |
| `body_deposit_types` (new) | `object_id`, `deposit_type` (id of the 19), `status` | text | one row per body and type that was assessed |
| | `endowment_kg` | kg | total ore mass, `E_k` |
| | `grade` | mass fraction | `g_k` of the key element or oxide |
| | `key_element` | text | the element or oxide the grade is of |
| | `age_window_gyr_start`, `age_window_gyr_end` | Gyr | when the type formed (for the scan's age line and for the census) |
| | `capped_by_host` | 0 or 1 | the grade was held down by the host supply |

All numeric columns get a `db/column_units.py` entry and a unit registry (GEN.568) entry; the new table and the three columns need a migration (renumbered at merge) and an epoch (`--affects fill`). The scan reads rows, never recomputes them; the object page shows the types with their status in plain words ("Ore: porphyry copper, rich" and so on) by a page change that is not part of this design.

## 8. Tests and the census (what the build slices must prove)

- **Gates (invariants, per slice).** No hydrothermal row on a stagnant-lid or dry body; no VMS without an ocean; no BIF without an oxygen source and an Fe(II) ocean inside its age window; no kimberlite without a craton and mantle carbon; no phosphorite without life; no unconformity uranium before its basin age. Each as a test on generated bodies, and each with a forced counter-example that the gate refuses (the pattern of the oddity forcing).
- **Anchors.** An Earth-like body (mobile lid, oceans, 4.5 Gyr, oxygenated) admits all of rows 1 to 18; a Mars-like body (stagnant, cold, dry) admits only the magmatic rows, weathering of its past, and none of the hydrothermal ones; a Venus-like body (stagnant, hot, dry at its age) admits the magmatic rows. Totals within an order of magnitude of the Earth's known endowments for the types where the source gives one.
- **Isolation.** Dial-style test: a body made again with the same seed has the same rows; turning a type off (a registry flag) changes no other draw.
- **Census (N26, the closing gate of the build plan).** Over a generated sample, the share of bodies that admit each type by body class, and the share of admitting bodies that have a `present` row, against the state rates (not against a literature rate, since no survey says how many planets have porphyry copper). New types appear only where the physics says (the reachability pattern of TEST.455).
- **Registry.** Every constant has a GEN.725 row; `scripts/wiring_guard.py --check` and the burn-down register pass (a slice that stores a field the scan does not read yet gets a wiring-register line naming GEN.799).

## 9. Slice order and what each needs

| Slice | Types | Needs first |
|---|---|---|
| GEN.1380 (lock, wave 4) | rows 6, 7, 8 | N5 (mantle oxides), N8 (minerals), the three derived fields of section 4 (crust thickness and sulfur index land here) |
| GEN.1381 | rows 1 to 5, 15 | GEN.1380's table and fields, N1 (Cu, Au), N9 for the brines (SEDEX) |
| GEN.1382 | rows 10, 11, 18 | N1 (Li), climate fields, GEN.1380's table |
| GEN.1383 | rows 12, 13, 14 | N7 (mantle carbon), N2 (REE, Nb), `craton_fraction` |
| GEN.1378 | rows 9, 16, 17 | GEN.750 timeline, N9, N18 |

The first slice creates the table and the three body columns; later slices add rows. GEN.798 (the resource design) reads this record and is written before GEN.1380 closes so the stored fields are enough (the build plan's rule).

## 10. Open points for the builders (not for Boss)

1. **Fe(II) ocean and the oxygen rise** for BIF: the GEN.750 timeline gives the oxygenic photosynthesis time; the abiotic oxygen route (hydrogen escape, photolysis) is a second source the file mentions. GEN.1378 decides whether the abiotic route gets its own gate.
2. **The oxygen fugacity scale.** The file states fO2 against FMQ; the generator stores ΔIW. The offset between the two buffers (about 3.5 log units at magmatic temperatures) is RECALLED and depends on temperature: GEN.1380 opens a source (the buffer paper is on the list of papers Boss cannot access, so use a second-hand citation) before gating porphyry or chromitite on it.
3. **Tonnage scatter versus total.** The file gives per-deposit ranges; the endowment sums a count of deposits. The count per unit opportunity is the least anchored number; the Earth totals are the only check, so GEN.1380 should state its tier plainly (L2 at best) in the registry row.
4. **Moons.** A moon's crust follows the same gates with its own state; a tidally heated moon (Io) gets magmatic rows from its heat flow, an icy-shell moon none. The slices test one rocky moon of each kind.
