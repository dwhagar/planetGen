**Headline:** The browser test site builds again, so the map browser tests run.

### Fixed
- **Every `site_app` browser test errored in setup with "Column count doesn't match value count" (TEST.474).** The bright stars `tests/test_web_a11y.py` writes lacked the `rotation_period_hours` column GEN.846 added to `bright_stars`; the rows now carry it (no rotation). Test-only.
