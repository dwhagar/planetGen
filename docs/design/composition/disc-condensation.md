# Disc condensation, species by species (GEN.1363)

Item GEN.1363, Lane 14. Part A (merged): `physics/disc_condensation.py`, the budget, the lines and the C/O across a line. Part B (built): the water share of the solids an icy body can reach follows the budget at its formation temperature, with its star's carbon and oxygen. Part C (built): a gas giant's C/O comes from the gas and solids of the disc where it formed. One epoch for the three.

## 1. What was missing

`condensate_shares` has three ices (water 182 K, ammonia 131 K, methane 41 K) and puts every carbon atom in methane and every nitrogen atom in ammonia, the equilibrium of [R2726: Table 10 and 11] at 10^-4 bar. The disc picture of [R3321, R3322] is different: about half the carbon is refractory grains, the volatile carbon is mostly CO and CO2, most nitrogen is N2, sulfur is mostly troilite, and the noble gases come out as clathrates in water ice. The repo's comets already follow a disc budget (`comet_ice_budget`, [R3322: Table 1]); the bodies of a system did not.

## 2. The model (part A)

Per 10^4 hydrogen atoms (the unit of [R3322: Table 1]), scaled by the star's element ratios over the rock-forming metals:

| Species | Line (K) | Budget | Source and tier |
|---|---|---|---|
| H2O | 182 | the oxygen the rock (22.75 percent of the Sun's) and the CO and CO2 leave | [R2726], T1 |
| NH3 (hydrate) | 131 | 30 percent of the nitrogen | [R2726] line, T1; split DSC5, T3 |
| NH4 salts | 131 | 10 percent of the nitrogen | DSC5, T3 recalled |
| H2S | 80 | 30 percent of the sulfur | line recalled DSC2, T3; share from 67P DSC6, T3 |
| Xe, Kr, Ar clathrate | 69, 53, 48 | all (equilibrium) | [R2726: Table 8], T1; amounts [R0168], T1 |
| CO2 | 47 | volatile carbon, 1 part in 6 of what is not methane | [R3322: Table 1], T1 |
| CH4 | 41 | 4.2 percent of the volatile carbon | line [R2726]; share from 67P [R3835], DSC4, T3 |
| N2 | 22 | 60 percent of the nitrogen, as molecules | line recalled DSC2, T3; share DSC5, T3 |
| CO | 20 | volatile carbon, 5 parts in 6 of what is not methane | [R3322: Table 1], T1 |
| refractory carbon | 500 | 50 percent of the carbon | [R3321, R3322], T2 |

The oxygen of the CO and CO2 comes out of the water because a gas takes oxygen whether or not it is frozen. For the Sun that leaves 2.17 oxygen as water per 10^4 hydrogen against Lodders's 3.78 (all the oxygen the rock leaves): the water mass of the nebula falls from 0.571 percent to 0.33 percent, and the water share of the solids between the water and ammonia lines from 0.539 to 0.34 (the grains add 0.145 percent). Part B decides how much of that reaches the icy-body water range; that is an output change, which is why part A alone changes nothing.

The mass of a species in percent of the nebula is its moles times its molar mass over water's, times the percent one mole of water makes in Table 11 (`PERCENT_PER_H2O_MOLE`), so the water of a disc with no carbon is exactly Table 11's.

The gas C/O (`gas_and_solid_c_over_o`) rises outward: 0.76 hot (the rock has the oxygen), 0.89 once water is frozen, about 1.05 once CO2 is frozen, then the gas is used up at the CO line; the solids stay oxygen-rich (0.44 with water and grains, 0.59 with everything). This is the lever for a giant's C/O (part C).

## 3. What parts B and C wire

**Part B.** `bulk_draw._icy_water_range` takes the top of an icy body's water from `disc_condensation.water_share_of_solids` (water ice over rock, carbon grains and every ice) at the body's formation temperature, with the carbon and oxygen of its star (`bulk_draw.ice_forming_ratios`). For the Sun that is 0.34 beyond the water line and 0.23 beyond the CO line, where Lodders's table gave 0.54 and 0.38; a star with more oxygen over its magnesium than the Sun's (the [Fe/H] 0 star of the tests, 1.25) gets 0.43. `condensate_shares.water_share_of_solids` is deleted. The icy-body water is a draw between a floor and this top, so every icy body's water ceiling falls by about a third.

Not wired: the ice over rock that sets the jump of the surface density at the snow line (`condensate_shares.WATER_ICE_SOLID_BOOST`, 2.17) stays on Lodders's table. The disc's solids would give about 1.6 against the grains-and-rock inside the line, which would cut every core beyond the snow line by a quarter; that needs the minimum-mass nebula weighed first (open point below).

**Part C.** `envelope_c_over_o(temperature, enrichment, carbon, oxygen)` mixes the gas of the disc with `sigma` times its solids until the envelope holds `enrichment` times the star's carbon and oxygen: a lightly enriched giant takes the gas C/O of its place, a heavily enriched one the oxygen-rich solids' (0.44 for the Sun's disc with every grain delivered). `mixture_subject` gives a giant that is not a moon its formation temperature (`182 K * sqrt(snow line / formation distance)`) and its star's carbon and oxygen; `mixture_giant` multiplies the result by the recipe's drawn `f7` over 0.57 (the center of f7's core, stand-in DSC7) so the observed spread stays. The sample of G2V giants has a disc median of 0.34 (0.28 to 0.40), where the retrieved giants' median is 0.57: the recipe draws enrichments of 3 to 80 times the star's, most of which is solids.

## 4. Open points


- Water mass: the repo's Lodders equilibrium had all spare oxygen as water; the disc measurement has CO and CO2 holding a large share. Part B takes the disc side, so the icy bodies get less water than before; comets (ice over rock near 1, [R2865]) are the check.
- N2 22 K and H2S 80 K are recalled; both are registered at T3 until a source is opened.
- Equilibrium supply; no trapping in water ice, no drift, no cooling disc.
- The solid boost at the snow line stays on Lodders's table (see section 3); a disc-based boost changes every giant core and is a separate decision.
- The giant C/O median (0.34 against 0.57) says either the recipe's enrichments are high for this model or the solids a giant swallows are less oxygen-rich than the table holds; DSC7 keeps the observed spread and lets the disc set the center.
