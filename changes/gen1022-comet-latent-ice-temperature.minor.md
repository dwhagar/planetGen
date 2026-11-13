**Headline:** The N2, argon, krypton and xenon of a comet's ice now rise and fall with its CO, as the comets measured do

### Changed
- GEN.1022 (option B of the comet species model): a comet's CO/H2O (row k5) now implies the temperature its ice formed at (the laboratory law through Halley's CO, 12.5 K per decade), and row k7 draws N2 and argon about the value that temperature gives between the warm end (67P, the median comet) and the cold end (C/2016 R2 for N2 per CO, Hale-Bopp for argon per water), with a residual of their own (N2 0.05 dex, argon 0.4 dex). Krypton follows the argon drawn for the same comet and xenon its krypton, by the ratios of 67P. Before, the seven species were drawn independently of each other and of the CO, which allowed a CO-poor, argon-rich comet that the data do not show.
- The CO/H2O of every comet is the k5 draw unchanged (median 2.5 percent, 0.65 dex). OCS, SO2 and O2 keep their GEN.1021 brackets. The draw count is unchanged. Row k5 stores its draw as `co_over_h2o`, and a species of a multi-species row can read the ones drawn before it (`draw_<name>`).
- Stored comet argon, nitrogen, krypton and xenon columns move for every comet, on top of GEN.1021's move.

generation-output: changed
generation-affects: fill
TODO items resolved: GEN.1022
