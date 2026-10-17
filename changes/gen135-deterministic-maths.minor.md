**Headline:** Distances, radii and rotations no longer depend on the Python release or the BLAS library, so a seed gives the same galaxy everywhere.

### Changed
- New `planetgen.util.detmath` (`hypot`, `dist`, `fsum`, `matmul`): square roots of squares added left to right, the exactly rounded sum, and a matrix product without BLAS. `math.hypot`, `math.dist` and `@` at about 70 call sites in `galaxy/`, `generation/`, `physics/`, `population/`, `db/` and `api/` now use them. `test_reproducible_draws.py` fails on `math.hypot`, `math.dist`, `@` and `numpy.dot` in the package outside `web/` (drawing code) and the migrations.
- Many stored positions and radii differ in the last bit, and each nebula's shape differs a little; this is generator epoch 2.

generation-output: changed

- The epoch rule needs the golden files, the epoch and the note's line together only once a golden file is checked in; until the golden test lands, the epoch and the note go together alone.

### Fixed
- TEST.205: the SQL portability lint no longer reads `ALGORITHM=INPLACE, LOCK=NONE` in an ALTER TABLE as a column named `lock` (broken on main since DB.18).

TODO items resolved: GEN.135, TEST.205 / Added: none
