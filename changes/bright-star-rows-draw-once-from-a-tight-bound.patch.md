**Headline:** The bright-star backfill draws each sky row once from a bound that follows the stars it is drawing, so a sky share walks 5 to 10 times faster.

### Fixed

- PERF.280: the bright-star candidate draw (GEN.308/GEN.309) takes one Poisson count and one stream a row and luminosity band (it took one a block of 256 wedges), and its bound is the most each population's density can be in the row times its share of the band (`galaxy.density.population_bounds_array`), where it was the whole disk's and bulge's largest density times the young stars' share, a hundred times what a cell expects. A deep sky tier listed 340 candidates for each star it kept; it now lists about one. The bounds of a couple of thousand rows are worked out together from arrays. The walk of one 1/16 share of a tier takes 0.1 to 1.8 s where it took 1.1 to 8.5 s. The stars are drawn from the same exact thinning, so counts and distributions are unchanged (tier 28, three shares, 64 seeds: 83.6 stars, 82.7 expected), but the random streams differ, so every candidate-path star moves: a lock item.

generation-output: changed
generation-affects: scatter, fill
