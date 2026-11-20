**Headline:** The ledger rows of a run Redis says is over are ended too, so no run stays `running` after `finish-run` carried it on (OPS.284).

### Fixed

- OPS.284 (bug): `planetgen finish-run` carried a run on, but the run it carried on stayed `running` in the ledger with its unfinished jobs `waiting` (found by the OPS.247 trial after a run whose process died on a full disk). The supervisor's `LostRuns` duty now also ends the ledger rows (jobs `canceled`, the run row `canceled`) of every run whose Redis record is `done`, `failed` or `canceled` while its ledger run is `running` or `waiting`. A live or held run is not touched.

generation-output: unchanged
generation-affects: none
