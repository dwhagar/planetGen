**Headline:** `generate system --output` runs with no database reachable: a run that names no database keeps its plan in its own Redis record.

### Fixed
- A run whose steps name no database (`system --output`, which touches none) keeps its plan, its end marks and its unit list in its run record in Redis, behind the same `planrows` functions the ledger answers through (`planrows.ledgered` chooses the store once), instead of opening the control database; its units have no rows.
- `run_galaxy._sector_args_of` takes a unit that names no database.
- Three web-generate tests keep up with the pages and classes of today: the estimate page's wording ("1.25 minutes", "1 TB"), and "Q", which is a planet class now, so the bad class in the refusal test is "AB".

generation-output: unchanged
