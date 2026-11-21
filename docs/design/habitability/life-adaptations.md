# Adaptations of life to cold, heat, dryness, salt and dim light (GEN.1725)

Status: the table and the traits cold and heat are built; nothing reads a limit yet (the refuge gate of GEN.1723 will). Code: `physics/data/life_adaptations.json`,
`physics/life_adaptations.py`, `physics/adaptation_traits.py` and `adaptation_draw.py` (GEN.876's traits, extended). Tests: `tests/test_life_adaptations.py`. Research:
`reports/life-adaptations-climate-2026-11-14.md` (GEN.1724) and `reports/ice-light-and-lichen-cold-edge-2026-11-14.md` (GEN.1727). Boss, 9:17 PM EST 2026-11-14: "implement it, make it happen".

## 1. What it is

Oxygen-making life lives far outside the warm, wet, lit sea. It lives in sea ice, under lake ice, in melt holes on glaciers, in rock, in desert crusts, in salt, and in hot springs, and each place needs
something the life has become: a pigment that works in almost no light, a sugar that keeps cells from freezing, a crust that wakes on a few hours of water. The table lists those adaptations, once, as data,
so that the generator can say how far a body's life goes and the page can say what it is adapted to.

## 2. The table

`life_adaptations.json` holds 36 rows (the adaptations of the GEN.1724 tables C, H, D, S and L), the six axes, the stages and the habitats.

* **Axes.** `t_cold` and `t_hot` (K, the local temperature of the water the life works in), `light_floor` (umol photons m-2 s-1), `water_hours` (hours of liquid water a year), `water_activity`, `ground_water_depth` (m)
  and `insulation` (K). Each axis has **one ordinary limit** (the life with no adaptation) with its own tier and tweak, so no two rows disagree about it, and a direction (wider is lower or higher).
* **Rows.** An adaptation: the climate it serves (the codes of the research), the kinds of life that carry it, the life classes (`all`, or the oxygen-making `1a`, `1b`, `1c`, or none), the axis, the limit at the tolerant
  (grade 2) and extreme (grade 3) grades (each wider than the one before), whether it keeps the life **active** or only **dormant** (freeze tolerance, desiccation tolerance and spores let life survive a limit, not work
  past it), the stage the life must have reached, the habitat it opens, its source rows, tier and tweak. A row that moves no number (the snow algae's sunscreen) has no axis and opens a habitat.
* **Stages.** The delay after the first life before a life can have an adaptation, from the one example (GEN.1724 LA8): cyanobacteria 0.8 to 1.1 Gyr, land mats 0.6, algae 2.2, mosses 3.3, rooted plants 3.4, C4 3.8. The animal stage has
  no delay of its own (the life timeline decides).
* **Habitats.** The refuge kinds the gate of GEN.1723 counts (open water, water under ice, melt film, crust land, seep, brine, hot pool) with the gate each one needs and the adaptations that open it.

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

## 5. What it leaves out

Production factors (how much each adaptation costs the life in growth) are empty except where the research gives one (the under-ice shade adaptation and cryoconite, 1.0): the research found no adaptation that raises refuge production
per area by more than a small factor, and no source for the costs. The build treats an empty factor as no change.
