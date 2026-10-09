### Added
- `planetgen fingerprint` (GEN.58): a canonical SHA-256 digest of each sector's generated content and one for the region (`--ring`, `--sector`, or the whole galaxy with its plan), skipping row ids, timestamps, update clocks and what is rebuilt from positions, so two builds of one galaxy can be compared sector by sector.
