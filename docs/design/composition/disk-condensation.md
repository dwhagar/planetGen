# Disc condensation, species by species (GEN.1363)

Item GEN.1363, Lane 14. Part A (built): `physics/disk_condensation.py`, the budget, the lines and the C/O across a line, wired to nothing. Part B (to build): the ice shares of a body follow the budget at the body's formation temperature. Part C (to build): the C/O of a giant follows the gas C/O at its formation distance.

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

## 3. Open points

- Water mass: the repo's Lodders equilibrium has all spare oxygen as water; the disc measurement has CO and CO2 holding a large share. Part A takes the disc side for the budget, so the icy bodies of part B get less water than before; comets (ice over rock near 1, [R2865]) are the check.
- N2 22 K and H2S 80 K are recalled; both are registered at T3 until a source is opened.
- Equilibrium supply; no trapping in water ice, no drift, no cooling disc.
