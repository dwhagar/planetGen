**Headline:** the table of how far each class may be fudged is now rebuilt by the update script whenever the stand-in rows or the class fit change, and a galaxy is not drawn on a table that no longer matches them.

### Changed

- `physics/fudge_table.py` (GEN.1317, new): says whether the derived table in the data directory is current (its header holds the generator epoch, the hash of the stand-in rows, the hash of the class fit and the matrix it was built from), stale or missing; `class_fit` reads that table when it is current and the committed fixture `physics/data/fudge_table.json` otherwise (registry rows FT6 and FT7).
- `scripts/build_fudge_table.py`: `--update` rebuilds the table from the stored sensitivity matrix in 0.3 s and says when the matrix itself is more than 50 epochs old; `--force` rebuilds a current table; `--status` prints the state. A full matrix refresh is still a separate run: about 1.5 hours on four cores, 6 hours on one.
- `update.sh` and `install.sh` (`scripts/deploy-common.sh`) call the update after the datasets step and warn at the end if it failed; `generation/run_galaxy.py` refuses a bulk draw on a stale or missing table the way it refuses missing datasets.
- Tests: the datasets-gate fixtures install a current table; `tests/test_fudge_table_update.py` (new, 8 tests) covers the build, the header, a stale epoch, an edited stand-in row, reproducibility, a broken file and the refusal.

generation-affects: none
