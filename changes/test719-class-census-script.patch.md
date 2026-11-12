**Headline:** The class census is now a repeatable script, `scripts/class_census.py`, with the same tables on every run.

### Changed

- `scripts/class_census.py collect` draws the census samples (natural, hab, stars, young, youngforced, capsoff) on fixed seeds through the host's memory and load gate, with `--workers`, `--forced-young` (every random primary star aged 0.5 to 50 Myr, for the disk-by-age row of GEN.1553) and `--caps-off all|d34`; a stopped run continues. `report` prints tables A to K of reports/class-census-2026-11-11.md in one markdown format from a directory of collected files, with numpy only (no pandas). On the stored 11-11 data it gives the report's class counts, funnel, flux-distance mass, pressure pile-up, star and value-check tables exactly.
- `tests/test_class_census_script.py` holds the table format, the funnel, the crosstab, the age bins, the fixed sample seeds and one real system drawn into a record the report reads.

generation-output: unchanged
