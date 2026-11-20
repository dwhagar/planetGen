**Headline:** A release's version now carries the database schema, generator epoch and API version, like `8.2.2100+db134.gen62.api5` (OPS.64).

### Added

- The release stamp writes a PEP 440 local suffix, `+db<schema>.gen<epoch>.api<api>`, after the release number in `_version.py`, the README badge and the changelog heading. The three numbers are read from the newest migration revision, `epochs.json` and `API_VERSION` at stamp time, so a database, a galaxy or a client can be matched to a release by its version alone. The user agent, the pages and the API show the full code.

### Changed

- BUILD, REVISION, the first-release names in `epochs.json` and `api-versions.json` and the release commit name still use the plain `MAJOR.REVISION.BUILD`. The PR check refuses a hand-set suffix as it refuses a hand-bumped version. The galaxy version key reads the release parts and ignores the suffix.
