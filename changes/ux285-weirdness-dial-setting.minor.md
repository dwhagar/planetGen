**Headline:** The settings page has a "Galaxy and world generation" group with the oddity rate, the dial that says how often generated worlds get rare oddities (1 is the natural rate, 0 turns them off).

### Changed
- New option `generation.weirdness_dial` (default 1, 0 to 10) in `config.json`, the settings page (group "Galaxy and world generation", shared with the cell cap) and `planetgen config`: the weirdness dial of the oddity table. A run that names its own dial (`SystemConfig.WEIRDNESS_DIAL`, the generator page, an administrator's force) uses that; otherwise a dial the settings name; otherwise the tuning's 1. It applies to worlds generated afterwards.
- `oddities.site_dial()` reads the option live (the settings are re-read only when the file changes); a settings file that does not fit the model leaves the tuning's dial.
- The settings page's administrator groups are eight now (the new group sits after "Speed and disk"); `config.schema.json`, `config.json.example` and `docs/reference/config.md` are regenerated.
- New tests in `test_oddities_wiring.py` (the file's dial wins when the run names none, a run's own dial wins over the file, a broken file falls back).

generation-output: unchanged
