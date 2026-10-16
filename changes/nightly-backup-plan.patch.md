**Headline:** A nightly full backup from the daily maintenance run is planned for 9.0, with its own schedule and retention settings.

### Added

- OPS.55: a backup step in the daily maintenance run that takes every part into one archive; options `settings.nightly_backup`, `settings.nightly_backup_hour` (default 4) and `settings.nightly_backups_kept` (default 7).

### Changed

- ADM.87 records each archive's kind (manual, settings-save, nightly) and prunes nightly archives by their own count.

TODO items resolved: none / Added: OPS.55
