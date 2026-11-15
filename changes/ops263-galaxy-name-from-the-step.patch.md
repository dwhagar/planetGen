**Headline:** A run's ledger row takes the database name from its step as written.

### Fixed
- The ledger run row names the galaxy from the step's `database` parameter (a name, or a dict's `database`) without building a connection config, so a step that gives the database alone no longer fails with a missing `host` (the stage hooks test red on main).
