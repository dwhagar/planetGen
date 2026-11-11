**Headline:** The one-off system page can now run with missing data files, as the Generate page can.

### Added
- **"Disregard missing data files" on the one-off system page (UX.283).** `/admin/generate/system` now has the same tick box the Generate page has; ticked, the page passes `--disregard-missing-data` to `planetgen system`, which then takes a labelled fallback for each missing dataset instead of refusing to run.
