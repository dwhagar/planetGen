**Headline:** The TODO list retires the three scatter-limit settings columns, which Core 3 built.

### Changed
- **DB.54 is done (#1887, migration 0101).** The scatter limits by place are three galaxy settings. The migration order is now DB.54 0101, DB.55 0102, GEN.336 0103, DB.53 0104.

### Added
- **TEST.197 and TEST.198, two browser tests that fail under load (bugs).** `test_page_is_accessible_and_fits[web.sector-phone-light]` and `test_galaxy_map_does_not_prefetch_on_a_data_saving_connection` each failed once in a full run and passed alone; goal issue #1694.

TODO items resolved: DB.54 / Added: TEST.197, TEST.198
