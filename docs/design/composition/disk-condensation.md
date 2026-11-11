# Disc condensation, species by species (GEN.1363)

Item GEN.1363, Lane 14. Part A (merged): `physics/disk_condensation.py`, the budget, the lines and the C/O across a line. Part B (built): the water share of the solids an icy body can reach follows the budget at its formation temperature, with its star's carbon and oxygen. Part C (built): a gas giant's C/O comes from the gas and solids of the disc where it formed. One epoch for the three.

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

**Part B.** `bulk_draw._icy_water_range` takes the top of an icy body's water from `disk_condensation.water_share_of_solids` (water ice over rock, carbon grains and every ice) at the body's formation temperature, with the carbon and oxygen of its star (`bulk_draw.ice_forming_ratios`). For the Sun that is 0.34 beyond the water line and 0.23 beyond the CO line, where Lodders's table gave 0.54 and 0.38; a star with more oxygen over its magnesium than the Sun's (the [Fe/H] 0 star of the tests, 1.25) gets 0.43. `condensate_shares.water_share_of_solids` is deleted. The icy-body water is a draw between a floor and this top, so every icy body's water ceiling falls by about a third.

Not wired: the ice over rock that sets the jump of the surface density at the snow line (`condensate_shares.WATER_ICE_SOLID_BOOST`, 2.17) stays on Lodders's table. The disc's solids would give about 1.6 against the grains-and-rock inside the line, which would cut every core beyond the snow line by a quarter; that needs the minimum-mass nebula weighed first (open point below).

**Part C.** `envelope_c_over_o(temperature, enrichment, carbon, oxygen)` mixes the gas of the disc with `sigma` times its solids until the envelope holds `enrichment` times the star's carbon and oxygen: a lightly enriched giant takes the gas C/O of its place, a heavily enriched one the oxygen-rich solids' (0.44 for the Sun's disc with every grain delivered). `mixture_subject` gives a giant that is not a moon its formation temperature (`182 K * sqrt(snow line / formation distance)`) and its star's carbon and oxygen; `mixture_giant` multiplies the result by the recipe's drawn `f7` over the median of `f7` (0.511, stand-in DSC7), to the power 0.5 (DSC8), so the observed spread stays without moving the center. GEN.1471 (measured in GEN.1395, reports/gen1395-giant-co-measure-2026-11-08.md): the first build divided by the core's center 0.57 instead of the draw's median, which lowered every giant by a sixth, and counted the recipe's spread in full on top of the disc's own. The sample of G2V giants then had a median of 0.33 where the retrieved giants' median is 0.57, 53 percent under 0.35 where LExACoM has 17. The enrichment range was not the lever: the recipe's enrichment is the planet's bulk heavy-element mass (18 times the star's for the median giant), and above about 5 times the envelope's C/O is 0.45 to 0.47 whatever the amount. The envelope's C/O is now set with a share of that bulk, `ENVELOPE_METAL_SHARE` 0.06 (DSC9): the bulk metallicity is an upper limit on the atmosphere's [R4350: Abstract], and the water of gas giants is substellar [R4663: Abstract]. The generated median is 0.57 and 16 percent are under 0.35 (0.59 and 14 on 331 fresh giants); 22 to 26 percent are above 0.8 where 4 to 9 percent are observed.

## 4. Open points


- Water mass: the repo's Lodders equilibrium had all spare oxygen as water; the disc measurement has CO and CO2 holding a large share. Part B takes the disc side, so the icy bodies get less water than before; comets (ice over rock near 1, [R2865]) are the check.
- N2 22 K and H2S 80 K are recalled; both are registered at T3 until a source is opened.
- Equilibrium supply; no trapping in water ice, no cooling disc. The drift of pebbles across the lines is `physics/disk_drift.py` (GEN.1473, `docs/design/composition/disk-drift-evaporation.md`): a giant's gas is weighted by the ratio of the pebble flux to the gas flux at the age it gathers it; the static budget here is the case where every weight is 1.
- The solid boost at the snow line is read from the disk budget (section 6, GEN.1396).
- The giant C/O upper tail (GEN.1471): a giant that took mostly gas has the oxygen-poor gas C/O of the disc, 0.89 for the Sun's (1.05 below 45 K), because the water is frozen out; 22 to 26 percent of the generated giants are above 0.8 where LExACoM has 4 to 9. The model has no drift and evaporation of water between the lines, which would spread the gas C/O between the star's and 1, and no intermediate state, so no setting of DSC8 or DSC9 closes it. It needs that model (open research).

## 5. The nitrogen, sulfur and noble gases of a giant (GEN.1372)

A giant's envelope was one scale for every metal but the carbon and oxygen. `envelope_enrichments(temperature_k, enrichment, ...)` gives the nitrogen, sulfur, argon, krypton and xenon their own factor over the star's amount: `(gas + sigma * solid) / (gas + solid)` of the element's atoms at the formation temperature, with the same `sigma` that sets the carbon and oxygen (`envelope_c_over_o`). Nitrogen is N2 gas until the ammonia line (131 K), then 40 percent ice; the sulfur is 30 percent H2S (frozen below 80 K) and 70 percent troilite, which is rock and takes the rock's scale; argon, krypton and xenon freeze below their lines (48, 53, 69 K) and are gas above. `mixture_subject` gives a giant its star's nitrogen, sulfur and argon over magnesium over the Sun's (`bulk_draw.disk_volatile_ratios`; krypton and xenon take the argon ratio), `mixture_giant._disk_enrichments` passes the factors to `element_budget`, and the constants record them as `enrichment_n`, `enrichment_s`, `enrichment_ar`, `enrichment_kr`, `enrichment_xe`. Phosphorus, neon and every rock-forming metal keep the one scale; a giant with no disk (a brown dwarf, a rogue, a hand-built subject) keeps it for every element.

Against the Galileo probe's Jupiter ([R5163: Table 2]: NH3 2.96 +/- 1.13, H2S 2.75 +/- 0.66 times solar) the model gives nitrogen 3.0 and sulfur 3.1 for a giant that formed at 110 K with four times its star's heavy elements. It misses the probe's roughly three times solar argon, krypton and xenon: a giant that formed warmer than 48 K holds them at the star's own amount, because the cold trapping or pebble drift that would carry them in is not modelled. Registry rows GE1 to GE3 (section 6.221).

## 6. The jump of the solids at the snow line (GEN.1396)

The planet mass law needs how much more solid a ring of the disk holds beyond the snow line than inside it (`formation.mmsn_surface_density_gcm2`, the isolation mass goes as the density to the 3/2). Until GEN.1396 that was `condensate_shares.WATER_ICE_SOLID_BOOST`, 1 plus the equilibrium water over the rock (2.17, [R2726: Table 11]). That counted water that the disk budget gives to CO and CO2, and left out the ammonia and the carbon grains. The jump is now `disk_condensation.SNOW_LINE_SOLID_JUMP`, from the budget of part A at the middle of the water and ammonia lines (156.5 K, `snow_line_solid_jump`): all the solids standing there over the rock. For the Sun that is 1.967 (rock 0.489, water 0.328, ammonia and salts, refractory carbon 0.145, other ices; percent of the nebula).

**Which divisor (decided by Lane 11 with Lane 14).** The divisor is the rock alone, not the rock and the grains (1.52). The solid density it multiplies is the minimum-mass nebula's rock normalization (`MMSN_SOLID_SURFACE_DENSITY_SOL_GCM2`, 7.0 g/cm^2, the silicates and metal that [R1854] counts), and the planet mass law then scales that by `PLANET_DISC_SOLIDS_ENHANCEMENT`, tuned to the observed occurrence (stand-in PA1). A jump relative to the rock keeps the inner disk as tuned and puts everything that freezes out beyond the line (and the carbon grains, which the normalization did not count) on the outside, as the ratio of solids that the budget gives. Dividing by rock and grains would hide the carbon, which the normalization leaves out inside the line too, in the enhancement that was tuned without it; it is the smaller number only because it assumes the normalization already holds the grains. Checked against the icy bodies: the budget gives a rock share of 0.51 of the solids beyond the line (1 over 1.967); icy dwarf planets are rockier, about 0.7 (Pluto's 30 to 70 ice and rock mix, [R0552]; a rock mass fraction of about 72 percent from the densities of the large Kuiper belt objects, [R0269]), because a body can lose ice by impacts and heat but cannot gain it. So the budget's rock share is the lower bound, in the right direction, and the old boost (0.46) was lower still.
