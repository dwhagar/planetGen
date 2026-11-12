**Headline:** the 2,000 K cap on a surface is gone: an air's hottest surface is now the temperature its main gas is half ionized at its pressure (12,500 K for 1 bar of carbon dioxide, tens of thousands of K for a deep air), and the greenhouse table no longer stops at 4,000 K.

### Changed

- `surface_air.MAX_SURFACE_TEMPERATURE_K` (2,000 K, a recalled stand-in) is replaced by `ionization.ceiling_k`: the Saha half-ionization temperature [R5238] of the air's main gas (the largest partial pressure, its atom of the lowest first ionization energy [R5228], `physics/data/ionization_energies.json`) at the air's pressure. Registry row CC7 is rewritten (tier T5, tweak L0) and IR18 added for the extrapolation below.
- `greenhouse_table.warming` no longer clamps its surface at 4,000 K (`FITTED_UP_TO_K` is where the fit ends): above it the same log ratio gives the warming, an extrapolation. Bodies the cap and the clamp held (86 census bodies sat at 2,000 K) now settle where their air puts them; the before and after of the 62 census systems is in the PR.
- A body held at the ceiling is a mistake of the model to explain: it logs one normal line (body, uncapped value, ceiling, air) and the test suite sets `PLANETGEN_CEILING_IS_ERROR`, so any seeded census or fuzz run that makes one fails.
- `trace_air.apply` now stores the molar mass, density and scale height of the air it leaves in every case (it skipped them when it removed equilibrium NO2 or N2O from a hot air without drawing a trace gas, which the 4,000 K airs exposed).
- With no 2,000 K cap a runaway steam or carbon dioxide air goes on to 4,000 K and thousands of bar, where the ideal-gas equilibrium stops converging (4 of 62 census systems failed with `EquilibriumError` between 3,900 and 9,300 bar above 3,980 K); the fixed point now keeps the last air the equilibrium solved, unsettled (`Air.converged` false), after its step-backs, where it raised. Finding for TODO: the equilibrium solver has an edge near 4,000 bar and 4,000 K that wants its own item (non-ideal gas).
- Tests that read the old constant now read the ceiling of the body's air.
- Changes the temperature, pressure and class of every body the cap or the clamp held, so it is a lock item: the golden and epoch are the train's.

generation-output: changed
generation-affects: fill, population
