### Fixed
- **The galaxy had almost no bulge, so edge-on views showed a flat disk
  (GEN.118), and no thick disk (GEN.119).** The density model's bulge was
  a 200 pc sphere holding 0.6% of the stars; the Milky Way's holds about
  31%. The model is now three published Milky Way components:
  - a thin disk (scale length 2.6 kpc, scale height 300 pc) carrying the
    spiral arms. `--disk-scale-height-pc` is now the exponential scale
    height surveys quote, so the disk is twice as thick as before;
  - a thick disk (2.0 kpc, 900 pc, 4% of the plane's density at the Sun);
  - a boxy bar bulge, Dwek et al. 1995's fit to the COBE/DIRBE image,
    angled 27 degrees from the Sun-center line.

  New `planetgen plan` defaults: `--disk-scale-length-pc 2600`,
  `--disk-scale-height-pc 300`, `--bulge-scale-radius-pc 1580` (along the
  bar), `--bulge-amplitude 3.11`, and the calibration point at 8.2 kpc
  (3.15 scale lengths). The young and intermediate populations' scale
  heights are now 50 and 150 pc. The Galaxy Map's prisms draw the same
  model. `tests/test_galaxy_milky_way.py` and two new math checks compare
  it with Bland-Hawthorn & Gerhard 2016.

  At the default shape a full galaxy has about 23.9 billion candidate
  sectors (was 12.3), 21.0 billion qualifying (was 10.3) and 299 billion
  systems (was 88), reaching 4.1 kpc above the plane (was 1.3). An
  existing galaxy keeps its stored shape numbers but the new formulas read
  them differently, so plan and generate it again.
- **The Galaxy Map's slab lines could cross after the view turned.** The
  buttons were ordered by each slab's middle but each line ends where the
  slab's outline comes nearest its button; when those ends come out of
  order the buttons now follow them.
