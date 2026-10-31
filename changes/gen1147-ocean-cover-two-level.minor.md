**Headline:** The share of a world's surface its ocean covers now follows a two-level surface (basins and continents) instead of a square-root law, and a mantle stores 2 to 8 Earth oceans by mass.

### Changed
- **Ocean cover (GEN.1147).** `hydrosphere.ocean_fraction` fills a surface of two Gaussians of height (ocean crust 60 percent, continental crust 40 percent, relief as 1 / g) with the body's water layer and returns the share under the water, in place of `min(1, sqrt(2 D / L))`. The four numbers are stand-ins tuned so Earth's 2.74 km layer covers 0.71 (row HY2). Icy bodies keep a relief of at most Earth's. The cover is now 0.27 at 0.2 km of water, 0.50 at 0.64 km and 0.59 at 1.3 km, so a world with a little less water than Earth has more land than before.
- **Mantle store (GEN.1147).** The most water a mantle holds is 2 Earth oceans at 0.5 Earth masses up to 8 at 5 (Guimond et al. 2023, row HY1), where it was a flat 10 oceans. `surface_water_fraction` takes the body's mass.
- **Not changed.** The Bern delivery law for the water a body is born with.

generation-output: changed
