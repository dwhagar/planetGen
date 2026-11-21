# Adaptations of life to cold, heat, dryness, salt and dim light (GEN.1725)

Status: the table and the traits cold and heat are built; nothing reads a limit yet (the refuge gate of GEN.1723 will). Code: `physics/data/life_adaptations.json`,
`physics/life_adaptations.py`, `physics/adaptation_traits.py` and `adaptation_draw.py` (GEN.876's traits, extended). Tests: `tests/test_life_adaptations.py`. Research:
`reports/life-adaptations-climate-2026-11-14.md` (GEN.1724) and `reports/ice-light-and-lichen-cold-edge-2026-11-14.md` (GEN.1727). Boss, 9:17 PM EST 2026-11-14: "implement it, make it happen".

## 1. What it is

Oxygen-making life lives far outside the warm, wet, lit sea. It lives in sea ice, under lake ice, in melt holes on glaciers, in rock, in desert crusts, in salt, and in hot springs, and each place needs
something the life has become: a pigment that works in almost no light, a sugar that keeps cells from freezing, a crust that wakes on a few hours of water. The table lists those adaptations, once, as data,
so that the generator can say how far a body's life goes and the page can say what it is adapted to.

## 2. The table

`life_adaptations.json` holds 39 rows (the adaptations of the GEN.1724 tables C, H, D, S and L, and the three of the chemical classes, GEN.1811), the eight axes, the stages and the habitats.

* **Axes.** `t_cold` and `t_hot` (K, the local temperature of the water the life works in), `light_floor` (umol photons m-2 s-1), `water_hours` (hours of liquid water a year), `water_activity`, `ground_water_depth` (m),
  `insulation` (K) and, for the classes of a solvent that is not water, `acid_strength` (percent sulfuric acid by mass) and `ammonia_fraction` (ammonia as a share of the water's mass). Each axis has **one ordinary limit** (the life with no adaptation) with its own tier and tweak, so no two rows disagree about it, and a direction (wider is lower or higher).
* **Rows.** An adaptation: the climate it serves (the codes of the research), the kinds of life that carry it, the life classes (`all`, or any of the classes of `life_classes`: the oxygen-making `1a`, `1b`, `1c`, `2` and the chemical `3a` to `9`, or none), the axis, the limit at the tolerant
  (grade 2) and extreme (grade 3) grades (each wider than the one before), whether it keeps the life **active** or only **dormant** (freeze tolerance, desiccation tolerance and spores let life survive a limit, not work
  past it), the stage the life must have reached, the habitat it opens, its source rows, tier and tweak. A row that moves no number (the snow algae's sunscreen) has no axis and opens a habitat.
* **Stages.** The delay after the first life before a life can have an adaptation, from the one example (GEN.1724 LA8): cyanobacteria 0.8 to 1.1 Gyr, land mats 0.6, algae 2.2, mosses 3.3, rooted plants 3.4, C4 3.8. The animal stage has
  no delay of its own (the life timeline decides).
* **Habitats.** The refuge kinds the gate of GEN.1723 counts (open water, water under ice, melt film, crust land, seep, brine, hot pool, and `chemical_medium`: the medium of a class that is not water) with the gate each one needs and the adaptations that open it.

`life_adaptations.limit(axis, grade, life_types=..., class_id=..., life_age_gyr=..., animals=...)` gives the widest limit of the rows that apply at a grade, or the ordinary limit. At the extreme grade for class 1a the table
reaches the measured edges: 253 K (lichen, field, -20 C), 346 K (cyanobacteria, 73 C), 0.04 umol m-2 s-1 (under the Arctic pack), 75 hours of water a year (hypoliths), 68 m of root.

## 3. The traits cold and heat

GEN.876's traits widen a limit the habitability index holds, by a factor. Cold and heat widen a limit in the units of an axis (the temperature the life is active down to and up to), so a grade means a limit from the table,
not a factor (`adaptation_traits.limit`). They are drawn once per body from its own stress, like the three before them, and stored in `life_traits` (`cold:3,heat:2`, shown as "cold (extreme)"):

* the stress is the shortfall in K of the **habitat temperature** (`life_supply.habitat_temperature_k`, the temperature of the water the life works in: the surface's under an open ocean or on a dry surface, the melting point of ice under an ice lid, the surface's plus the rise with depth in an aquifer; TEST.800) past the ordinary limit (273.15 K for cold, 313.15 K for heat), scored on a straight line from 0 to the width of the table for class 1a (20 K for cold,
  33 K for heat); the chance and the grade follow the same law as the other traits (`P_CAP` 0.8 times the score; extreme with chance `GRADE3_SHARE` times the score squared);
* they come **last** in the draw order, so the radiation, sulfide and carbon dioxide traits of every body are what they were (a test pins it); a body now draws ten uniforms, two for each of five stresses;
* **nothing reads them yet.** The scores of the habitability index and the life classes are unchanged, and a trait never makes a world habitable to people. The refuge gate of GEN.1723 reads `adaptation_traits.limit` to widen where
  a class 1 life can work; widening the class 1 gate in `life_classes.json` waits for it (the same file).
* the **light floor** under ice has no trait: nothing stores the ice a world's ocean wears (GEN.1728), so its stress cannot be read. The same goes for dryness (no stored hours of liquid water) and brine (the water activity
  has its column, but no body has been scored for it). They are rows of the table with no draw.

The habitat temperature is one value for the body, not the local one of each habitat (GEN.1723 stores that): a world with a cold mean and a warm tropics still draws the cold trait a little too often, and one with a warm mean and cold poles too seldom. Before TEST.800 the trait read the surface mean, which gave a cold trait to 95 of 143 life worlds in a 40-system sample (66 percent; 94 of the 137 chemical-class worlds, whose crust is warm under a frozen surface); with the habitat temperature it is 13 of 144 (9 percent), and 23 of 144 carry any trait.

## 4. Reuse

The traits extend GEN.876 (`adaptation_traits`, `adaptation_draw`) and add no second draw; the table reads the same tier and tweak ladders as the registry (`docs/design/research/stand-in-values.md`); the refuge limits of GEN.1722 (RF1 to RF11)
are rows, not new constants. Habitats are read by GEN.1723. Light under ice is `ice_light` (GEN.1726).

## 4a. The record of a body (GEN.1808)

`physics/adaptation_record` puts the table and the refuge gate together for one finished body (Boss, 5:54 PM EST 2026-11-20: apply to all life types, with "information to each about the adaptations that were necessary to make it work"). It reads the refuge (`surface_refuge.refuge_of`: the cold and hot ends of the temperatures of its wet lit steps, the hours of liquid water, the least light under the ice, the habitats counted) and the body's water table and water activity, and for each axis whose ordinary limit does not reach the demand finds the lowest grade of the table (`life_adaptations.limit`, for the body's class, life age and kinds of life) that does, and the rows that give it. A demand no grade reaches is kept as `grade` `None`: the life there is past the table. Habitat rows with no axis are listed for each habitat the refuge counted and the class may use. Light is converted at the gate's own floor (stand-in AE1, 5 micromoles per watt) and the temperature range holds the middle 95 percent of the wet lit time (AE2).

The record is `planet.life_adaptations`, a list of dicts (`life_class`, `axis`, `habitat`, `demand`, `grade`, `rows`), `None` for none, the entries of each class on the body one class after the other. It is made once by the life pass after the air history, moves no number of the body and is not stored yet (GEN.1809). Every class of `life_classes` has a record.

### The chemical classes (GEN.1811)

The classes 3a to 9 eat no light and have no refuge, so their needs are not read from the refuge. `adaptation_record.chemical_demands` reads the one habitat temperature of the class gate (the cloud layer's for 5a) and the medium of the class, and puts them to the table the same way. What the research of GEN.1810 (`reports/life-refuges-classes-3a-9-gen1810-2026-11-21.md`) found is in the rows and the record:

* **3a, 3b, 4a, 4b, 4c (water).** They carry the cold-active (C16) and the hyperthermophilic (H3) microbe rows and one new cold row, C17: metabolism in a brine to 248.15 K (-25 C, the cold edge of a permafrost brine microbe, Marginal). The hot end is the record (122 C, 395 K), which the gate already holds. A cold water habitat is recorded in `brine`, a hot one in `hot_pool`.
* **5a, 5b, 9 (sulfuric acid).** Row A1 reaches 4.7 percent acid (growth around pH 0) and 11 percent (*Picrophilus* adapted to 1.2 M); the classes sit at 81 percent and more and the water activity there is 0.004 or less (the lowest a cell has divided at is 0.585). The record says **past the table** for the acid strength, the water activity and any temperature: no measured life lives in the medium.
* **6 (dense CO2) and 7 (a hydrocarbon sea).** The solvent holds no water, so the demand is a water activity of 0, past every row; no life has been measured in either.
* **8 (ammonia-water).** Row N1 reaches 0.43 percent ammonia of the water (growth directly in 0.25 M ammonia); the gate asks 10 percent, 25 times that. The record says past the table. The gate stays (TODO's answer): the gap is recorded, not closed.

No gate moved and no generated number changes: the record reads the finished body. Habitat pressure (the 140 MPa record of a piezophile) is a separate item (GEN.1819).

## 5. What it leaves out

Production factors (how much each adaptation costs the life in growth) are empty except where the research gives one (the under-ice shade adaptation and cryoconite, 1.0): the research found no adaptation that raises refuge production
per area by more than a small factor, and no source for the costs. The build treats an empty factor as no change.
