**Headline:** `rescore-habitability` now re-scores only the systems whose stored numbers were not worked out by the running code, so a rerun after no change touches nothing.

### Changed
- **A system carries the stamp of the code that scored it (GEN.1479).** New narrow table `system_summary` (migration 0236) holds, for each star system, `habitability_stamp`: a hash of the physics source the habitability score is made from (`db/derived_stamp.py`). The system save and `store.refresh_habitability` write it; every system already stored starts with none, so the first run after the migration scores them all once, as the command always did.
- **`planetgen rescore-habitability` takes only the stale systems.** A system whose stamp is the running code's is left alone, the units are cut from the first stale system to the last, and a run that stopped carries on where it stopped. `--all` re-scores every system as before. The plan line says how many were left alone. A side table, not columns of `star_systems`, because that table has a FULLTEXT key and a new column would copy every row.

generation-affects: none
