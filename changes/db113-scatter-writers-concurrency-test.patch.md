**Headline:** The scatter write fix for the MySQL deadlock is now tested with eight real writers at once.

### Added

- **DB.113: a test with eight writer threads on their own connections, writing bright-star and phenomenon scatter batches into the same tables with the Object ID step.** The writers meet at a barrier between their inserts and their ID step, so the clash is certain: with the old write at least one writer is killed with MySQL error 1213, and with `write_committed` all of them land, every row has an ID, no ID is shared and each sector is numbered with no gap.

generation-output: unchanged
