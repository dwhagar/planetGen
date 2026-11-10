**Headline:** A molten world whose buffered oxygen fills the gas pressure at 6000 K now generates, with the carbon gas solved at a slightly higher pressure, instead of failing the whole system.

### Changed

- GEN.1513: `mixture_lava.vapor_mixture` solved the carbon gas of a volatile-bearing melt at the surface pressure with the oxygen held at the buffer's fugacity. At 6000 K the buffer's fugacity can be 10^5.13 bar against 1.4e5 bar of gas (94 percent), the solve then has no room for the carbon and raised "no convergence at 6000 K and 143497 bar", which failed the generation of the zeta Ophiuchi example about one run in eight. New `carbon_gas` tries again at a gas pressure higher by `CARBON_GAS_PRESSURE_STEP` (1.25) up to `CARBON_GAS_PRESSURE_STEPS` (8) times, as the GEN.1296 share rule already raises it for the melt's own oxygen; it converges within three steps from 4500 to 6000 K and fugacities of 1e3 to 1e7 bar, and ends in a `MixtureError` otherwise. A gas that converged is left at its pressure, so nothing that generated before changes.
- Check: 40 seeds of the zeta Ophiuchi example built with no error (seed 9 was the failing one before).

generation-affects: none
