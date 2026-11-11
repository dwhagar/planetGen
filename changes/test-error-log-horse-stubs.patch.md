**Headline:** The error-log tests for a killed work horse run again after the memory-gate and killed-flag changes.

### Changed
- `test_error_log.py`: the stub job of the three horse-killed tests carries a Redis stand-in with no memory-gate claim (PERF.254 reads the claims table when a horse is killed), and the stub for `runrecord.finish_run` takes the `killed` flag OPS.144 added. Tests only; no code under `src/planetgen` changes.

generation-output: unchanged
generation-affects: none
