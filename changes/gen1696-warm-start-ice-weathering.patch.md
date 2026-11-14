**Headline:** A world that can be warm or frozen now comes out warm, and a frozen world's weathering stops, so more worlds with water keep their oceans.

### Changed
- **A body with the ice-albedo law enters the air solve from the warm side (GEN.1696).** The solve started from the airless temperature, so a world with two possible climates (ice-free albedo 0.30 above 305 K, ice-covered 0.62 below 240 K) always found the frozen one. A planet cools from its formation, so a body that holds surface water and is not a giant now starts from `air_branches.WARM_START_K` (320 K) when its airless temperature is lower; the steam start of GEN.1238 and bodies without the law are unchanged. Oddity W20 (the other branch of a bistable climate) now finds the frozen alternative of such a body, where it used to find the warm one.
- **Continental weathering works only on the ice-free land (GEN.1697).** `surface_air._weathering_land` takes the land share `1 - ocean_cover` times one minus `IceAlbedoLaw.ice_fraction` of the body's own law at the surface temperature, so under a global ice sheet the weathering all but stops and the volcanic CO2 builds up until it melts the ice (Hoffman et al. 2017). A body with no ice law keeps its whole land. No new law or constant: the ice share is the one that already sets the albedo.
- On 800 paired Earth copies (the system's own step order, fresh body per trial) the surface temperature under 250 K falls from 176 to 131 of 800 and the median surface temperature rises from 293.7 K to 294.6 K; Class M goes 14 to 12 (small on its own; the methane haze GEN.1695 and the galactic ozone flag GEN.1694 gate the rest).

generation-output: changed
generation-affects: fill
