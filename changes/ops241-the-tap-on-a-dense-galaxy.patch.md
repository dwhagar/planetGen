**Headline:** The queue commands' declared reads and writes are now held to a second tap run on a dense galaxy, and the declarations that run showed to be wrong are corrected.

### Changed
- `tests/test_command_touches.py` runs a second tap test on the galaxy of the bright star scatter tests, dense enough that every scatter draws rows and the sectors hold systems, with `--all` for the habitability and ledger re-scores and `--since` for the link, settle and backfill runs. The first galaxy is too sparse for a command that sends no statement without rows.
- `queue/touches.py`: `rescore-habitability` and `refill-effect-ledger` declared one table read and none written; they read a whole system (`SYSTEM_RELOAD_READS`) and write `planets`, `moons`, `system_interests` and `system_summary` (the habitability) and `effect_ledger` and `system_summary` (the ledger). `empty-nebula-field` also reads the skeleton and writes `scatter_epochs`. Output of generation is unchanged.
