**Headline:** The job engine's ledger spells the state of a run or job that was stopped on purpose `canceled`, as the rest of the project does.

### Changed
- **The ledger's canceled state is spelled the American way (OPS.215).** The spelling check flagged the British spelling in `queue/ledger.py` and `control_schema.sql`, and the run state on Redis is already `canceled`, so the state stored in `engine_runs.state` and `engine_jobs.state` is now `canceled` too. The control schema is v24: a control database at v23 has the three ledger tables dropped and made again by `store._drop_engine_ledger_v23`, since the v23 check allowed only the old spelling and nothing has written a row. No generator output changes.

generation-output: unchanged
