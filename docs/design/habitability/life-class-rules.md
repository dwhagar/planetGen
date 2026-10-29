# Life class rules

Generated from `src/planetgen/physics/data/life_classes.json` by `scripts/build_life_class_rules.py`; do not edit by hand. The classifier (`planetgen.physics.life_classes`, GEN.886) reads that file and nothing else; the design is [life-classes.md](life-classes.md). Every number in the file carries a source tier (T1 to T9) and a tweak level (L0 to L4) (see [stand-in-values.md](../research/stand-in-values.md)).

Rule file version 1. Reference star 5772 K, reference band 400 to 700 nm.

## The classes

| Class | Name | Drinks | Eats | Breathes | Gives off | Reaches | Tier |
|---|---|---|---|---|---|---|---|
| 1a | Oxygenic light-eaters: visible pigment (Chlorophyll a) | water | light 400 to 700 nm, carbon dioxide and water, plus nitrogen, phosphorus and iron | oxygen once it exists (aerobes); before that nothing | oxygen (net, after the sinks fill) and buried organic carbon; methane only through a small food web | Animals | T1 |
| 1b | Oxygenic light-eaters: blue and UV tuned (Blue-Optimized Porphyrins) | water | light 300 to 500 nm, carbon dioxide and water | oxygen once it exists | oxygen and buried organic carbon | Animals | T6 |
| 1c | Oxygenic light-eaters: far-red (chlorophyll d and f) | water | light 700 to 800 nm and beyond, carbon dioxide and water | oxygen once it exists | oxygen and buried organic carbon | Animals | T5 |
| 2 | Anoxygenic light-eaters (Bacteriochlorophylls) | water | near-infrared light 800 to 1050 nm, carbon dioxide and a donor: hydrogen, ferrous iron or hydrogen sulfide | nothing (anaerobes); free oxygen poisons them | sulfate, sulfur, Fe(III), organic carbon; methane from their food web | Microbial mats; animals only if another oxidant source exists | T1 |
| 3a | Hydrogen eaters: geochemical hydrogen (serpentinisation, vents) | water | hydrogen from rock and water, and carbon dioxide | carbon dioxide (methanogens, acetogens) | methane and acetate | Microbial | T2 |
| 3b | Hydrogen eaters: radiolytic hydrogen and sulfate (deep crust, aquifers) | water | hydrogen from the radiolysis of water, and carbon dioxide or sulfate | carbon dioxide or sulfate | methane, acetate, sulfide | Microbial | T2 |
| 4a | Sulfur, iron and nitrogen eaters: oxygen as the oxidant | water | hydrogen sulfide, ferrous iron, ammonium or methane | oxygen | sulfate, Fe(III), nitrogen, carbon dioxide | Animals (chemosymbiotic tubeworm, clam and mussel communities) | T2 |
| 4b | Sulfur, iron and nitrogen eaters: nitrate and nitrogen oxides as the oxidant | water | hydrogen sulfide, ferrous iron, ammonium or methane | nitrate and nitrogen oxides | sulfate, Fe(III), nitrogen, carbon dioxide | Animals (small) | T2 |
| 4c | Sulfur, iron and nitrogen eaters: chlorate and perchlorate as the oxidant | water | hydrogen sulfide, ferrous iron, ammonium or methane | chlorate and perchlorate | sulfate, Fe(III), nitrogen, carbon dioxide, chloride | Microbial; animals are speculative (T8) | T8 |
| 5a | Acid-world life: cloud-deck plankton | concentrated sulfuric acid (81 to 100 percent) | reduced sulfur, organic carbon, blue and UV light | the acid itself (reduced to SO2 and sulfur); traces of oxygen | SO2, sulfur, carbon dioxide, water | Mats and aerial plankton; at most simple colonial forms | T6 |
| 5b | Acid-world life: pool and crust mats | concentrated sulfuric acid | reduced sulfur, organic carbon | the acid itself | SO2, sulfur, carbon dioxide, water | Mats | T6 |
| 6 | Dense-CO2 life | liquid or supercritical carbon dioxide | hydrogen as fuel | carbon dioxide (the solvent is the oxidant) | methane, water (which separates out) | Microbial, at most mats | T7 |
| 7 | Cold hydrocarbon-sea life | liquid methane and ethane | hydrogen and acetylene (hydrogenation) | nothing | methane | Microbial mats | T8 |
| 8 | Ammonia-water life | water with ammonia | hydrogen and chemical gradients in cold alkaline brine | carbon dioxide or carbonate | methane, nitrogen | Microbial | T8 |
| 9 | Silicon-bearing life in acid | concentrated sulfuric acid | as Class 5, with organosilicon chemistry | the acid itself | SO2, silica (solid SiO2) | Microbial | T9 |

## What each class needs

| Class | Solvent | Energy | Gates (all must hold) | Windows | Supply |
|---|---|---|---|---|---|
| 1a | water (factor 1) | light, 400 to 700 nm, efficiency 1 | hydrosphere in surface ocean, hycean; pressure_pa >= 10.0; co2_partial_pa >= 5.0; band_light_vs_earth >= 0.01 | temperature_k 255 to 395 (falls off over 15 below, 15 above) | none |
| 1b | water (factor 1) | light, 300 to 500 nm, efficiency 1 | hydrosphere in surface ocean, hycean; pressure_pa >= 10.0; co2_partial_pa >= 5.0; band_light_vs_earth >= 0.01 | temperature_k 255 to 395 (falls off over 15 below, 15 above) | none |
| 1c | water (factor 1) | light, 700 to 800 nm, efficiency 0.8 | hydrosphere in surface ocean, hycean; pressure_pa >= 10.0; co2_partial_pa >= 5.0; band_light_vs_earth >= 0.01 | temperature_k 255 to 395 (falls off over 15 below, 15 above) | none |
| 2 | water (factor 1) | light, 800 to 1050 nm, efficiency 0.8 | hydrosphere in surface ocean, hycean; pressure_pa >= 10.0; co2_partial_pa >= 5.0; band_light_vs_earth >= 0.01; share_o2 < 0.03 | temperature_k 255 to 370 (falls off over 15 below, 15 above) | donor_supply (placeholder 1 until it is stored) |
| 3a | water (factor 1) | chemical (internal and tidal heat) | pressure_pa >= 0.0 | temperature_k 255 to 395 (falls off over 15 below, 15 above) | h2_supply (placeholder 1 until it is stored) |
| 3b | water (factor 1) | chemical (internal and tidal heat) | pressure_pa >= 0.0 | temperature_k 255 to 395 (falls off over 15 below, 15 above) | radiolytic_h2_supply (placeholder 1 until it is stored) |
| 4a | water (factor 1) | chemical (internal and tidal heat) | share_o2 >= 0.001 | temperature_k 255 to 395 (falls off over 15 below, 15 above) | donor_supply (placeholder 1 until it is stored) |
| 4b | water (factor 1) | chemical (internal and tidal heat) | pressure_pa >= 0.0 | temperature_k 255 to 395 (falls off over 15 below, 15 above) | oxidant_supply (placeholder 0.3 until it is stored) |
| 4c | water (factor 1) | chemical (internal and tidal heat) | pressure_pa >= 0.0 | temperature_k 255 to 395 (falls off over 15 below, 15 above) | oxidant_supply (placeholder 0.1 until it is stored) |
| 5a | sulfuric (factor 0.1) | chemical (internal and tidal heat) | acid_medium >= 1.0 | none | none |
| 5b | sulfuric (factor 0.1) | chemical (internal and tidal heat) | acid_medium >= 1.0 | temperature_k 283 to 373 (falls off over 15 below, 15 above) | none |
| 6 | co2 (factor 0.1) | chemical (internal and tidal heat) | share_co2 >= 0.5; pressure_pa >= 518000.0; co2_phase in liquid, supercritical | temperature_k 217 to 400 (falls off over 20 below, 40 above) | h2_supply (placeholder 0.5 until it is stored) |
| 7 | hydrocarbon (factor 0.2) | chemical (internal and tidal heat) | none | temperature_k 90 to 190 (falls off over 10 below, 10 above); pressure_pa 50000 to 500000 (falls off over 30000 below, 300000 above) | hydrogen_acetylene_supply (placeholder 0.5 until it is stored) |
| 8 | ammonia_water (factor 0.1) | chemical (internal and tidal heat) | ammonia_water_mass_fraction >= 0.1 | temperature_k 176 to 273 (falls off over 10 below, 10 above) | h2_supply (placeholder 0.5 until it is stored) |
| 9 | sulfuric (factor 0.1) | chemical (internal and tidal heat) | acid_medium >= 1.0; silicon_rich_crust >= 1.0 | temperature_k 283 to 363 (falls off over 15 below, 10 above) | none |

## Why each number is what it is

### Class 1a: Oxygenic light-eaters: visible pigment (Chlorophyll a)

- solvent water: T9, L3. Liquid water; the index's own L_solv (water activity and chaotropicity) applies
- energy light: T1, L0. A star spectrum rich in 400 to 700 nm: F, G, K and bright early M; the Sun's share of photons in the band is the reference (efficiency 1)
- gate hydrosphere: T9, L3. Sunlit water: an open surface ocean or a hycean sea (life-classes.md 3.1, Class 1 and 2 window); under ice or in an aquifer there is no light
- gate pressure_pa: T9, L3. The pigments that harvest light need air (10 Pa or more), the gate the old LIFE_CHEMICALS table used (GEN.92, GEN.647)
- gate co2_partial_pa: T7, L3. Carbon dioxide at least about 5 Pa for a carbon fixer (life-classes.md 3.1, Class 1 window)
- gate band_light_vs_earth: T7, L3. Light in the pigment's band at least 1 percent of Earth's (life-classes.md 3.1); the 0.01 percent allowance for deep microbes waits for the habitat-depth model
- window temperature_k: T7, L3. Mean habitat 255 to 395 K for microbes (life-classes.md 3.1, Class 1 window); smooth falloff outside like the index's own terms
- note: Class 1a on an Earth-like world is exactly today's PHI_bio (life-classes.md 8.1). The old oxygen gate of Chlorophyll a is gone: Class 1 makes the oxygen (GEN.855).

### Class 1b: Oxygenic light-eaters: blue and UV tuned (Blue-Optimized Porphyrins)

- solvent water: T9, L3. Liquid water
- energy light: T7, L3. Hot bright stars: late B to F; a UV shield from fluorescent proteins (life-classes.md 3.1); tier T6 to T7: a reasoned invention
- gate hydrosphere: T9, L3. Sunlit water: an open surface ocean or a hycean sea (life-classes.md 3.1, Class 1 and 2 window); under ice or in an aquifer there is no light
- gate pressure_pa: T9, L3. The pigments that harvest light need air (10 Pa or more), the gate the old LIFE_CHEMICALS table used (GEN.92, GEN.647)
- gate co2_partial_pa: T7, L3. Carbon dioxide at least about 5 Pa for a carbon fixer (life-classes.md 3.1, Class 1 window)
- gate band_light_vs_earth: T7, L3. Light in the pigment's band at least 1 percent of Earth's (life-classes.md 3.1); the 0.01 percent allowance for deep microbes waits for the habitat-depth model
- window temperature_k: T7, L3. Mean habitat 255 to 395 K for microbes (life-classes.md 3.1, Class 1 window); smooth falloff outside like the index's own terms

### Class 1c: Oxygenic light-eaters: far-red (chlorophyll d and f)

- solvent water: T9, L3. Liquid water
- energy light: T7, L3. K and M stars; an efficiency penalty for the lower photon energy: the chemistry is T5, the penalty 0.8 is T7 (life-classes.md 3.1)
- gate hydrosphere: T9, L3. Sunlit water: an open surface ocean or a hycean sea (life-classes.md 3.1, Class 1 and 2 window); under ice or in an aquifer there is no light
- gate pressure_pa: T9, L3. The pigments that harvest light need air (10 Pa or more), the gate the old LIFE_CHEMICALS table used (GEN.92, GEN.647)
- gate co2_partial_pa: T7, L3. Carbon dioxide at least about 5 Pa for a carbon fixer (life-classes.md 3.1, Class 1 window)
- gate band_light_vs_earth: T7, L3. Light in the pigment's band at least 1 percent of Earth's (life-classes.md 3.1); the 0.01 percent allowance for deep microbes waits for the habitat-depth model
- window temperature_k: T7, L3. Mean habitat 255 to 395 K for microbes (life-classes.md 3.1, Class 1 window); smooth falloff outside like the index's own terms

### Class 2: Anoxygenic light-eaters (Bacteriochlorophylls)

- solvent water: T9, L3. Liquid water
- energy light: T7, L3. Near-infrared light, 800 to 1050 nm; the lower photon energy costs efficiency (T7)
- gate hydrosphere: T9, L3. Sunlit water: an open surface ocean or a hycean sea (life-classes.md 3.1, Class 1 and 2 window); under ice or in an aquifer there is no light
- gate pressure_pa: T9, L3. The pigments that harvest light need air (10 Pa or more), the gate the old LIFE_CHEMICALS table used (GEN.92, GEN.647)
- gate co2_partial_pa: T7, L3. Carbon dioxide at least about 5 Pa for a carbon fixer (life-classes.md 3.1, Class 1 window)
- gate band_light_vs_earth: T7, L3. Light in the pigment's band at least 1 percent of Earth's (life-classes.md 3.1); the 0.01 percent allowance for deep microbes waits for the habitat-depth model
- gate share_o2: T9, L3. Free oxygen below 3 percent poisons an anoxygenic or hydrogen-eating microbe: the existing gate (life-classes.md 3.1, Class 2 window)
- window temperature_k: T7, L3. Habitat 255 to 370 K (life-classes.md 3.1, Class 2 window)
- supply donor_supply: T9, L4. The donor supply (hydrogen, Fe2+ or H2S delivered) is not stored yet (GEN.887); the placeholder 1 assumes an anoxic ocean has its donors

### Class 3a: Hydrogen eaters: geochemical hydrogen (serpentinisation, vents)

- solvent water: T9, L3. Liquid water in contact with rock (ocean, ice-covered ocean or aquifer)
- energy chemical: T9, L3. The share of the body's heat flow that life can use (CHEMICAL_SHARE), as the index already scores sub-ice life; no light needed
- gate pressure_pa: T9, L3. Any pressure: a deep ocean or aquifer is not exposed to the air
- window temperature_k: T7, L3. Mean habitat 255 to 395 K for microbes (life-classes.md 3.1, Class 1 window); smooth falloff outside like the index's own terms
- supply h2_supply: T9, L4. The hydrogen supply from serpentinisation (rock type, water, heat) is not stored yet (GEN.887); the placeholder 1 assumes ultramafic rock is present

### Class 3b: Hydrogen eaters: radiolytic hydrogen and sulfate (deep crust, aquifers)

- solvent water: T9, L3. Liquid water in contact with rock
- energy chemical: T9, L3. The share of the body's heat flow that life can use; radiogenic heat and dose make the hydrogen
- gate pressure_pa: T9, L3. Any pressure: a deep ocean or aquifer is not exposed to the air
- window temperature_k: T7, L3. Mean habitat 255 to 395 K for microbes (life-classes.md 3.1, Class 1 window); smooth falloff outside like the index's own terms
- supply radiolytic_h2_supply: T9, L4. Radiolytic hydrogen from the rock's uranium, thorium and potassium and the water is not stored yet (GEN.887); the placeholder 1 assumes enough

### Class 4a: Sulfur, iron and nitrogen eaters: oxygen as the oxidant

- solvent water: T9, L3. Liquid water
- energy chemical: T9, L3. The energy of the fuel meeting the oxidant, scored on the share of the heat flow life can use
- gate share_o2: T9, L3. An oxygen sea: oxygen at least 0.1 percent of the air for the oxidant to reach the vents (life-classes.md 3.1, Class 4 window)
- window temperature_k: T7, L3. Mean habitat 255 to 395 K for microbes (life-classes.md 3.1, Class 1 window); smooth falloff outside like the index's own terms
- supply donor_supply: T9, L4. The reduced fuel supply (H2S, Fe2+, NH4+, CH4) is not stored yet (GEN.887); the placeholder 1 assumes vents or seeps supply it

### Class 4b: Sulfur, iron and nitrogen eaters: nitrate and nitrogen oxides as the oxidant

- solvent water: T9, L3. Liquid water
- energy chemical: T9, L3. The energy of the fuel meeting the oxidant
- gate pressure_pa: T9, L3. Any pressure: a deep ocean or aquifer is not exposed to the air
- window temperature_k: T7, L3. Mean habitat 255 to 395 K for microbes (life-classes.md 3.1, Class 1 window); smooth falloff outside like the index's own terms
- supply oxidant_supply: T9, L4. Nitrate and nitrogen-oxide delivery is not stored yet (GEN.887); the placeholder 0.3 stands for a slow oxidant supply

### Class 4c: Sulfur, iron and nitrogen eaters: chlorate and perchlorate as the oxidant

- solvent water: T9, L3. Liquid water
- energy chemical: T9, L3. The energy of the fuel meeting the oxidant
- gate pressure_pa: T9, L3. Any pressure: a deep ocean or aquifer is not exposed to the air
- window temperature_k: T7, L3. Mean habitat 255 to 395 K for microbes (life-classes.md 3.1, Class 1 window); smooth falloff outside like the index's own terms
- supply oxidant_supply: T9, L4. Perchlorate delivery is not stored yet (GEN.887); the placeholder 0.1 stands for a very slow supply

### Class 5a: Acid-world life: cloud-deck plankton

- solvent sulfuric: T9, L3. Concentrated sulfuric acid as the liquid; the existing SOLVENT_FACTOR of sulfuric acid, 0.1 (planetGen's choice)
- energy chemical: T9, L3. Chemical energy of the acid medium; light too in the cloud layer
- gate acid_medium: T6, L4. A liquid sulfuric acid layer is present: nothing in the generator models it yet (GEN.888), so this class scores 0 until then

### Class 5b: Acid-world life: pool and crust mats

- solvent sulfuric: T9, L3. Concentrated sulfuric acid as the liquid; SOLVENT_FACTOR 0.1
- energy chemical: T9, L3. Chemical energy of the acid medium
- gate acid_medium: T6, L4. A liquid sulfuric acid pool is present: not modelled yet (GEN.888)
- window temperature_k: T6, L3. Warm dry ground, acid pools 283 to 373 K (life-classes.md 3.1, Class 5b window)

### Class 6: Dense-CO2 life

- solvent co2: T9, L3. Liquid or supercritical CO2 as the solvent; no existing factor, 0.1 like the other strange solvents (planetGen's choice)
- energy chemical: T9, L3. Chemical energy of hydrogen fuel, scored on the share of the heat flow life can use
- gate share_co2: T7, L3. A dense CO2 air or sea: carbon dioxide at least half of the air (life-classes.md 3.1, Class 6 window)
- gate pressure_pa: T2, L0. Liquid CO2 needs above 5.2 bar (the triple point is 5.18 bar)
- gate co2_phase: T2, L0. The CO2 is liquid or supercritical at the surface
- window temperature_k: T8, L3. Liquid below 304 K, supercritical above that and about 73 bar (life-classes.md 3.1)
- supply h2_supply: T9, L4. The hydrogen supply is not stored yet (GEN.887); 0.5 stands for a half of Class 3's

### Class 7: Cold hydrocarbon-sea life

- solvent hydrocarbon: T9, L3. The existing SOLVENT_FACTOR of a hydrocarbon sea, 0.2 (planetGen's choice); the index's `solvent()` rule says when a methane sea exists
- energy chemical: T9, L3. Chemical energy of hydrogen and acetylene, scored on the share of the heat flow life can use
- window temperature_k: T8, L3. 90 to 190 K (life-classes.md 3.1, Class 7 window)
- window pressure_pa: T8, L3. N2 and CH4 air around 1.5 bar (life-classes.md 3.1, Class 7 window)
- supply hydrogen_acetylene_supply: T9, L4. The hydrogen and acetylene supply (haze organics) is not stored yet (GEN.887); 0.5 is a stand-in

### Class 8: Ammonia-water life

- solvent ammonia_water: T9, L3. Water with ammonia as the liquid; 0.1 like the other strange solvents (planetGen's choice)
- energy chemical: T9, L3. Chemical energy of hydrogen and gradients, scored on the share of the heat flow life can use
- gate ammonia_water_mass_fraction: T8, L3. Ammonia at least 10 percent of the water mass (life-classes.md 3.1, Class 8 window); the ammonia inventory is not stored yet (GEN.887), so this class scores 0 until then
- window temperature_k: T8, L3. 176 to 273 K: the water-ammonia eutectic is about 176 K (recalled)
- supply h2_supply: T9, L4. The hydrogen supply is not stored yet (GEN.887); 0.5 is a stand-in

### Class 9: Silicon-bearing life in acid

- solvent sulfuric: T9, L3. Concentrated sulfuric acid as the liquid; SOLVENT_FACTOR 0.1
- energy chemical: T9, L3. Chemical energy of the acid medium
- gate acid_medium: T6, L4. A liquid sulfuric acid layer or pool is present: not modelled yet (GEN.888)
- gate silicon_rich_crust: T9, L4. A silicon-rich crust (GEN.441 crust composition); not exposed as a quantity yet
- window temperature_k: T9, L3. As Class 5; silicones resist 98 percent acid at 90 C (life-classes.md 3.1, read in the climates report)

## Classes that do not drink water

The hydrogen inventory of the chemistry term is held at 1 Earth units: T9, L4. The index's chemistry term needs carbon, nitrogen and hydrogen in Earth units, and hydrogen comes from the water; a class that does not drink water gets its hydrogen from the supply term (GEN.887), so its inventory is held at this floor, and phosphate is the index's NON_AQUEOUS_PHOSPHATE_UMOL_L.
