**Headline:** The tests' in-process work queue is undone in the right order, so a patch a test makes to `ensure_workers` no longer outlives it.

### Changed
- **The inline work queue patches through the test's own `monkeypatch` again.** PERF.106's shared context (`tests/inline_queue.py`) undid its patches with a context of its own, so a test that patched `workers.ensure_workers` after it, with a fixture older than the context, had the stand-in put back after the context had restored the real function; the stand-in then stayed for the rest of that worker's tests, and `test_workers` and `test_update_guard_work` failed on later runs. The autouse fixture passes its `monkeypatch` to the context; a fixture of wider scope (`run_cli`) still uses a context of its own.
- No product code changes.

generation-output: unchanged

TODO items resolved: none (a failure found by Core 5 in PERF.106's gate, TEST.250 in the findings list). Added: none.
