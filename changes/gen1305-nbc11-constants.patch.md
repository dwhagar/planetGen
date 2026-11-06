**Headline:** The stand-in row for the molecular-fraction law now names its constants, so the wiring guard sees they are read.

### Fixed
- Stand-in row NBC11 (section 6.179) named `nebula_composition.KMT_*`, which the wiring guard read as the name `KMT_` and found unread. It now lists the five constants (`KMT_CHI_COEFFICIENT`, `KMT_CHI_METALLICITY`, `KMT_LOG_TERMS`, `KMT_S_COEFFICIENT`, `KMT_FRACTION_TERMS`) that `molecular_fraction` reads, and its register line is gone. No code changes.

generation-output: unchanged
