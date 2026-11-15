**Headline:** `update.sh` no longer stops at the math check on a database that has no job-ledger tables yet.

### Fixed
- A run whose control database is missing or has no ledger tables yet (the update's math check runs at stage 4, before stage 5 migrates) keeps its plan in its own Redis run record and writes no ledger row; so does a run with no database at all (`generate system --output`). The store is chosen once per run and kept in its record. A run on a database that has the tables is unchanged (OPS.263, TEST.769).
