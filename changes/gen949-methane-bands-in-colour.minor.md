**Headline:** The colour of a giant, and of any air with methane in it, now uses Karkoschka's measured methane absorption at 1 nm resolution instead of 35 coarse bins.

### Changed

- The gas table's coarse bins give the mean methane absorption of a 10 nm cell, which hides the structure that sets a giant's colour: the strong bands (727, 790, 840 and 890 nm) are many times deeper than their cell mean and the windows between them nearly clear. `uv_cross_sections.methane_optical_depth` now reads the new setup dataset `methane_visible` ([R2259], NASA PDS volume GBAT_0001, `data/1995low.tab`, 101 kB; `update.sh` downloads it, fill AD12) and takes the depth of each cell between 400 and 1050 nm as the optical depth of the mean transmission of its samples. A run that disregards missing data takes the old bins and says so.
- A giant whose only band-bearing gas is methane now reports `colour_quality = "full"`; ammonia and water above 1e-4 still make it `partial`.
- Run `python3 -m planetgen.cli.datasets fetch` (or `update.sh`) once to get the new file; generation will not start without it unless told to disregard missing data.

generation-output: changed
generation-affects: fill
