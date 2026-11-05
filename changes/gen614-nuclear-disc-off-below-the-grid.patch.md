**Headline:** A galaxy whose centre is smaller than the grid no longer packs the nuclear stellar disc into one cell: the disc is switched off when it is shorter than 2 pc.

### Changed
- `galaxy.density.model_terms` (GEN.614): the nuclear stellar disc's length is 0.057 times the bulge radius, so under a bulge radius of about 35 pc it was smaller than the grid's 2 pc core ball (`tuning.NUCLEAR_RESOLUTION_PC`) while the nuclear cluster was softened to it, and the disc put about five times too many systems into one cell. The disc is now off there: its mass and density are zero, and the galaxy's other components are unchanged. Chosen over softening because the cluster is one Hernquist cusp that can be held flat to the ball's mean, while a disc has two lengths and no single flat core that holds its mass. A galaxy with a bulge radius of 35 pc or more, the Milky Way's 1.58 kpc included, is unchanged.
- GEN.615 (no thick disc, bar, stellar halo or nuclear structures under 3e8 solar masses of stars) was already built with the Structure dial and has its tests in `tests/test_halo_threshold.py`; nothing to add.

### Added
- `tests/test_gen614_nuclear_disc_grid.py`: bulge radii from 0.7 pc to 1.58 kpc, the cut at 35 pc, no division by zero, a switched-off disc adds no density or mass in the model, the potential components or the array form.

generation-output: changed
TODO items resolved: GEN.614 (golden and epoch ride train 10; GEN.615 was already built)
