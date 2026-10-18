**Headline:** Every galaxy draws its own thick disk, bar and stellar halo from its seed, and one heating law sets the heights and motions of all its stars.

### Changed
- **The thick disk, the bar and the stellar halo are drawn from the galaxy seed, not copied from the Milky Way** (`galaxy/structure.py`, `galaxy/structure_draw.py`). The structure is rebuilt from the seed on every load and stored nowhere (`galaxy_shape` gains no column; DB.66 keeps storing the drawn values in 8.3). A shape built without a seed uses `REFERENCE_STRUCTURE`, the published Milky Way values, as a check target and not as a generator default.
  - Circular speed comes from the baryonic Tully-Fisher relation (McGaugh 2012, 0.0375 dex) until the potential gives it.
  - Thick disk: the mass ratio sits between Yoachim and Dalcanton and Comeron 2018 with 0.268 dex of scatter, the length ratio follows Comeron 2018, the height ratio is 2.5 to 4.2, and the amplitude is worked out from the three, replacing the fixed 0.04.
  - Stellar halo: the share of the stars is `0.05 (M*/5e10)^0.79` with 1.01 dex of scatter (a stand-in until the halo research lands), cut to 0.1% to 30%.
  - Bar (GEN.533, sources in `docs/design/gen533-bar-sources.md`): presence from Erwin 2018 by stellar mass, the in-plane axis ratio from Diaz-Garcia 2016, a boxy or thin vertical ratio, an angle uniform on 0 to 180 degrees, and the unbarred bulge's axis ratios from Costantin 2018. A galaxy without a bar has no bar gap in its star formation.
- **One law sets the scale heights and the velocity dispersions of the disc populations** (the dispersion of the gas the stars are born from, 6.5 +/- 1 km/s as a molecular stand-in until GEN.536 reads the CO values, a central disc dispersion of 0.248 Vc, an age exponent of 0.44), replacing the fixed height ratios and dispersions of the young, intermediate and old discs, the halo and the bulge. The reference check puts the old disc 32% under Sharma and Hayden's 21.1 km/s at the Sun (inside the plan's 35%) and the young disc at 0.22 of the old disc's height (65 pc for 300 pc; observed about 50 pc). The atomic gas's 10 km/s (Tamburro 2009) gave 0.40 and twice the molecular clouds the catalogues hold, so it is not used.
- Removed from `tuning.py`: the thick-disk, bulge axis, bar angle, halo mass and system-fraction constants and the fixed height ratios and dispersions the draw replaces.
- The generator epoch is 20.
generation-output: changed

### Added
- `tests/test_structure_draw.py` (13 tests): determinism, load equals plan, the Tully-Fisher draw, the thick-disk bracket and bound, the halo law, the bar's ranges, barred against unbarred, the one-law heights and dispersions, and that the stored row gains nothing. The skeleton tests in `tests/test_galaxy_gen.py` store a fixed seed and build the toy shape with its drawn structure.

TODO items resolved: GEN.532, GEN.533. Added: none.
