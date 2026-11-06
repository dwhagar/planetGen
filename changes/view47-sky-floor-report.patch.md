**Headline:** A run that charts the sky now says how faint a star its stored sky shows at each distance, and the unused whole-sphere draw wrappers are gone.

### Changed
- Sky charting already runs by default (`queue.sky_limit` 7.5, set by VIEW.25), so a generation run, the API job and the setup call all draw the sky tiers of the observer cells they touch. The run log now ends the "Sky tiers" line with the sky's floor at 10, 50, 250, 1000 and 3000 parsecs (`sky_tiers.floor_table`): the dimmest star the stored sky shows there through the least dust of any direction, and whether the table draws it.
- `bright_stars.backfill_sphere` and `backfill_mass_sphere` are removed. Charting draws each batch of a sphere itself (`luminosity_sphere_batches`, `mass_sphere_batches` with `backfill_cells` and `backfill_mass_cells`) because each cell has its own ceiling; the tests keep the same two helpers.

generation-output: unchanged
