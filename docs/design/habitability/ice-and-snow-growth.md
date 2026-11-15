# How thick a world's ice and snow are (GEN.1728)

Status: built (the law); not yet read by anything. Code: `physics/ice_growth.py`. Tests: `tests/test_ice_growth.py`. Registry: section 6.295 of the gap plan (IG1 to IG8). Filed by GEN.1726, whose light law
([light-through-ice.md](light-through-ice.md)) takes a thickness from its caller; GEN.1723's refuge term is the caller.

## 1. What it gives

`lid_of(mean_k, amplitude_k, absorbed_flux_w_m2, ocean_flux_w_m2, freezing_k)` returns the ice thickness and the snow depth on it, in metres, and whether the lid is perennial. The inputs are what the
generator already has for any world: the mean surface temperature, the half-range of the year's cycle (`seasonal_swing`), the starlight the surface absorbs, the heat the water gives up from below and the
freezing point of its water. Nothing is stored; the caller asks when it needs the lid, so there is no new column and no epoch.

## 2. The ice

* **A season's growth.** The ice grows from below as the cold above carries heat out through it. Stefan's law with the snow and the air in the way ([R5269], equations 18 and 21):
  `H = sqrt((k he)^2 + 2 k u / (rho L)) - k he`, where `u` is the freezing degree days (kelvin seconds), `he = h_s / k_s + 1 / h_air` is the resistance of the snow and the air, `k` the ice's conductivity
  (`A / T`, the hydrosphere's ice law, 2.07 W per m and K at melting), `rho` 917 kg per m3 and `L` 3.352e5 J per kg. The freezing degree days of a sine year about the mean come out in closed form
  (`freezing_degree_days`), checked against a direct sum. Two thousand degree days make 1.4 m with no snow, the first-year Arctic lid.
* **A lid that never thaws.** When the warmest month is under freezing there is no summer to start over. The ice settles where the heat it conducts down balances the water's flux and the wear of the surface:
  `H = k ((Tf - T) / (F + rho L a) - h_s / k_s)` ([R5271]'s balance, annual mean). With the standard Arctic case (2 W per m2 from the water, 0.4 m of ablation a year, 0.34 m of snow, a 256 K mean) it gives
  2.9 m, the multiyear thickness the model of [R5271] reaches (2.88 m). With no ablation it is the linear form of the conduction limit of an ice shell.
* **The cap.** The seasonal thickness is held under the conduction limit of the coldest month (the hydrosphere's `ice_shell_thickness_km`, reused), which only binds for an icy body with a feeble interior.

## 3. The snow

The snow of a year is the Earth's share of the absorbed starlight that goes into the water cycle (about 0.32, the weathering law's own energy limit [R1660]), over the latent heat of the water at the mean
temperature. It falls in the freezing share of the year, about two thirds of it stays on the ice (11 cm of water standing in May of the 17 cm that fell on the Arctic ice, [R5270]) and it lies at 300 kg per m3
([R5270]). A world with a faint star, a very cold air or a short freezing season gets little snow.

**Check.** With 75 W per m2 absorbed at 255 K the law gives 0.65 m of snow against the 0.34 m the Soviet stations measured on multiyear ice ([R5270]), a factor 1.9 high: it has no dry cold air to thin the
snowfall. It errs toward more snow, so toward less light under the ice. The registry row IG6 says so and what would replace it.

## 4. What it leaves out

* The melt of the summer. A lid whose summer thaws part of it keeps the seasonal thickness, whichever is the case; where the lid thaws completely it is the season's growth, which is the right answer.
* Snow-ice, the slush that forms where the weight of the snow pushes the ice under the water; ridging and drift; the salt of the ice (the freezing point is the caller's `freezing_k`).
* The wear of the surface of a lid under another air. One measured case (40 cm a year, [R5271]) is used for every perennial lid (IG8).
* A snowfall that depends on the dryness of the air (IG6).

## 5. Reuse

`hydrosphere.ICE_CONDUCTIVITY_A_W_M`, `ICE_DENSITY_KG_M3` and `ice_shell_thickness_km` (the ice's conduction), and the weathering law's latent heat (`continental_weathering.latent_j_m3`, factored out of
`runoff_m_yr` with no change of value), its Earth precipitation and absorbed-flux constants. `ice_light.transmission` takes the thickness it gives.
