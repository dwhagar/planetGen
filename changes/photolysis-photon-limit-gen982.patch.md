**Headline:** The photolysis table can now say how much of a gas the light breaks in a thick air, and a star's far-ultraviolet light at a body can be asked for; nothing reads them yet.

### Changed
- **GEN.982, parts 1 and 2.** `photolysis_lifetimes.json` gains the mean UV cross sections of water, ammonia, methane, carbon dioxide and hydrogen sulfide (Hu, Seager and Bains 2012, table 6), and `photolysis.shielded_column_mol_m2` and `photon_limited_flux_mol_m2_s` give the column that absorbs the light and the most of a gas the light can break per square metre per second (for methane over the Earth's area about 78 Tmol a year, the order of the anoxic methane fluxes). `stellar_library.far_uv_power_w` and `activity.far_uv_flux_ratio` give a star's 91.2 to 200 nm flux at a distance relative to the quiet Sun at 1 AU, the number the photolysis rates take. Both are for part 3 (the methane haze supply), which is a lock branch; no generated value changes.
generation-output: unchanged
