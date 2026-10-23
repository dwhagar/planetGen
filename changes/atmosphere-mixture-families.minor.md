**Headline:** The atmosphere mixture draw covers every class: Venus-like and Mars-like air, abiotic oxygen, Titan-like and Pluto-like air, comet comae, and the vapour over a lava world.

### Added
- **`draw_mixture` families (GEN.479, GEN.456, GEN.499).** `physics/mixture_evolved.py` (a Venus-like or Mars-like body is the outgassed gas without its water and most of its sulfur over a crust that buffers the oxygen; an abiotic-oxygen world keeps the oxygen of its lost water; Titan-like air is nitrogen and methane at the methane's saturation; a Pluto-like air stands over its ices by Raoult's law), `physics/mixture_coma.py` (a comet's coma from the drawn ratios to water) and `physics/mixture_lava.py` (the vapour over an ideal-solution melt of the bulk silicate Earth's oxides at the held oxygen fugacity, from the melt's element potentials: 1.8e-4 bar at 2000 K and 1.1e-2 bar at 2500 K against Kite et al. 2016's 2e-4 and 8e-3, sodium and potassium first, SiO above O2 beyond about 2900 K). A white-dwarf planet with a vapour exosphere takes the lava chemistry; one without has no gas. Not called by the generator yet: no generator epoch.
- **Two alkali oxide liquids in the species table** (`Na2O(L)`, `K2O(L)`, from the same NASA Glenn database) and `Species.supercooled_g_over_rt`, the liquid's free energy below its melting point.

### Changed
- **The planned rocky classes U, Y and Z have the rocky rows of the possibility tables** (`ROCKY` in `scripts/atmosphere_recipes_data.py`, the data file rebuilt): the rows of table j are gated on them and table j3 and j4 need a retention row. No existing class changes, so no generated body changes.
- **`mixture_draw.family_of(recipe, subject)` takes the subject**, because a hot-zone body is a lava world only where its surface is molten.
