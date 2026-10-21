**Headline:** The galaxy's disks are exponential in height, the form their scale heights are measured with, so a galaxy's midplane stars are twice as dense for the same mass.

### Changed
- **The thin and thick disks' vertical profile is `exp(-|z| / h)`, not the isothermal sheet `sech^2(z / 2h)`** (`galaxy/density.py` `_vertical`, `generation/object_first.py`, `galaxy/skeleton.py`, `html/static/galaxyprisms.js`). The scale heights the model uses (BHG16's 300 pc and 900 pc) are exponential fits to star counts; a `sech^2` sheet with the same tail has a column of `4 h` times its midplane density where the fit's is `2 h`, so for one stellar mass the midplane was half the fits' (GEN.404 made the density follow the mass and exposed it). The mass of each disk is unchanged (its column is `2 h` times the midplane now, `2 pi L^2 x 2 h` in `density._disk_mass`, `shape_draw.disk_masses` and `density.bulge_to_total`), so at the reference radius between two arms the density is 0.750 of the density the earlier calibration pinned there (it was 0.376 under GEN.404) and the column through that point is 524.26 pc of relative density before and after. The nuclear stellar disc keeps its `sech^2` height (`density._sech_vertical`).
  - The open clusters' heights (`generation/open_clusters.py`) and the neutron-star and black-hole height factor (`galaxy/remnant_distribution.py`) follow the stars' profile.
  - The reference galaxy's bulge amplitude is 3.11 / 2 (`mathcheck.REFERENCE_BULGE_AMPLITUDE`): BHG16's 3.11 was measured against the `sech^2` disk, and the bulge keeps its 31% of the mass.
  - `MAGNETAR_YOUNG_STAR_SHARE` is 0.0016 (it was 0.0019), the young population's share of the stars over the new profile.
  - The reference galaxy's outline reaches shell 3,596 (it was 3,147 under GEN.404, 3,762 before it); the rest of the gap to the pinned density is the mean system mass (GEN.545).
- The generator epoch is 40.
generation-output: changed

### Added
- `tests/test_galaxy_milky_way.py`: the thin disk's column is `2 h` times its midplane, the density at the reference radius is 0.75, and its column is 524.26 as before.

TODO items resolved: GEN.544. Added: none.
