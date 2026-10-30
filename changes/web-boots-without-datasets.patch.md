**Headline:** The web app starts even when a downloaded data set is not there yet, so an update no longer turns every page into a 500 while it downloads.

### Changed
- `update.sh` pulls the new code in stage 1 but downloads data sets in stage 3 (52 minutes on the live run). Apache's daemon reloaded on the new code in between, and `create_app()` died at boot on the missing `mamajek_dwarf_sequence` file (the class reference catalog is built at start). The catalog now builds at start when it can and is built by the first /classes request after the download when it cannot; the missing file is logged as a warning.
