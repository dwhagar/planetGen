**Headline:** Two test files again build their stand-in planets with a build record, so main is green again after GEN.1651.

### Fixed
- `test_tidal_spin` and `test_radiation_dose_field` made stand-in planets without a `build` record, which `generate_spin` and `update_exposure` now write to (GEN.1651); the stand-ins carry a `BuildRecord()`.

generation-output: unchanged
