**Headline:** The backup archive (settings now, control database and parity data for 9.0) and the restore and repair chain are planned, with the settings page staying in 8.1.

### Changed

- ADM.83 and ADM.85 take their backups from one archive mechanism (ADM.87): one bz2 tarball per backup, `settings.backups_kept` archives kept (default 20, 1 to 200).

### Added

- ADM.87 (the archive and the settings part, 8.1), and for 9.0 the umbrella ADM.91 with ADM.88, DB.40 (control and parity parts), ADM.89, DB.41 (their restores), DB.42 (repair from an archive), ADM.90 (the repair page), TEST.180 and DOC.56, in 9.0 phases 3 and 4.
- The design note `docs/design/backups-and-repair.md`.

TODO items resolved: none / Added: ADM.87, ADM.88, ADM.89, ADM.90, ADM.91, DB.40, DB.41, DB.42, TEST.180, DOC.56
