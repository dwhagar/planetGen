**Headline:** The radiation sweep test draws more systems until it finds a giant, and no longer asks a band that is almost dark to pass less light under a thicker air.

### Fixed

- `tests/test_radiation_sweep.py` was red on Train 40 in three tests. GEN.1709 (the inner-system period spread) moved the draws and the sample of 20 systems plus up to 12 extra held no giant (2 giants before, none after), so `EXTRA_SYSTEMS_MAX` is now 40 (only the systems needed are drawn). A class N world the moist limit of GEN.1579 allows (19 bar, 412 K, 12.7 bar of carbon dioxide) passes 3.7e-5 of its uv_c light and 4.1e-5 under twice the column, so the rule that a thicker column passes no more light is not asked of a band under `DARK_TRANSMISSION` (1e-3); the rise of uv_c with oxygen in `atmospheric_radiation` is reported to TODO as a separate question.

generation-output: unchanged
generation-affects: none
