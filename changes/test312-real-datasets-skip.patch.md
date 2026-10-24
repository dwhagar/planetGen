**Headline:** The real-data checks of the downloaded data sets now wait for every data set, so a machine that updated before the ozone file was added no longer shows a failing test.

### Changed
- **TEST.312.** `test_datasets_real` skipped only when the first three data sets were missing, so a machine whose cache held those but not the ozone cross sections (AD9) failed `test_the_catalogue_reports_nothing_missing`. It now skips unless the whole catalogue is in the cache.
