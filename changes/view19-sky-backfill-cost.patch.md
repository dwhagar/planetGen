**Headline:** A read-only script measures what a sky of a given limiting magnitude costs the star backfill.

### Changed

- `scripts/measure_sky_backfill.py` builds a galaxy shape from a seed in memory, integrates the model's own star densities to count the stars a 6.5 or a 7.5 magnitude sky shows (with and without a flat dust bound), draws the luminosity bands and the mass rings over spheres of cells to time them, and draws the two ledgers in both orders to show they never overlap. It writes nothing to the database. The measured numbers are in reports/system-sky-backfill-cost-2026-11-05.md (VIEW.19).

generation-output: unchanged
