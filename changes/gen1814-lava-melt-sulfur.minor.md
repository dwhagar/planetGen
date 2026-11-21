**Headline:** A lava world's melt holds sulfur: a third sulfur sink, sized by the sulfide saturation of the molten layer, sits beside the crust sulfate and sulfide sinks.

### Added

- `physics/melt_sulfur.py`: the sulfur content at sulfide saturation (SCSS) of a basalt melt, ln S[ppm] = 11.79 - 7073/T - 92 P/T (Steenstra et al. 2022, R5340; fitted from 42 values, 1 atm to 40 GPa, 1,573 to 4,605 K), turned into a capacity in mol per square metre: the melt column reaches the depth where the solidus (`tidal_heating.solidus_k`, inverted) has risen to the surface temperature, at most the body's magma-ocean base, times the column's mean SCSS. `sulfur_balance.balance` takes it as `melt_capacity_mol_m2`, after the sulfate and sulfide crust sinks and before the air takes the rest (`SulfurBalance.melt_mol_m2`); `surface_air._shares` passes it. A surface under the melting point has no melt and no capacity, so cooler bodies are untouched. Stand-ins MS1 to MS4 in section 6.319 (Confident; melt depth and capacity-as-upper-bound registered with their tweak levels). Of 200 test systems only two lava worlds change (GEN.1814).

generation-output: changed
generation-affects: fill
