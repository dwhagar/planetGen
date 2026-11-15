**Headline:** the end of a stage is three jobs in the ledger, and the Redis marks of how far it got are gone.

### Changed

- OPS.218 (slice 1): the cut of a stage also writes its `conclude`, `finish` and `release` jobs, each with an edge to the one before. The end of a stage marks each done as it goes, and a retry or the reconciler reads the rows instead of `runrecord.ended_key`, which is removed. Opening the next stage stays an RQ job until OPS.225, which removes it.
