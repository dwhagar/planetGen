**Headline:** Two tests that went red when the star disk columns and the primordial heat model landed now match them, and a release note that failed the version check is fixed.

### Fixed
- **`test_every_column_is_written` listed four star columns as never written (TEST.352).** `disk_kind`, `disk_mass_solar`, `disk_inner_au` and `disk_outer_au` are written now, so they leave the test's "NULL in this galaxy" list.
- **`test_every_generated_rocky_body_stores_its_draws_and_an_urey_ratio_that_is_an_output` capped the Urey ratio at 1.5 (TEST.353).** A tiny moon (0.0004 Earth masses, 4.66 Gy old) starts cold and is heated by its own decay, so less heat flows out than the decay makes and the ratio is 3.0. That is what the GEN.699 thermal history says for such a body, so the cap is gone and the test still checks the ratio is positive and finite.
- **`changes/map-unresolved-glow.minor.md` had a `### Tests` heading (TEST.357),** which `bump_version.py --check-pr` rejects. Its one line is now a bullet under `### Added`.
TODO items resolved: TEST.352, TEST.353, TEST.357. Added: none.
