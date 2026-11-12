**Headline:** The warm-up log tests read the one duration form that OPS.173 introduced, which turns main green again.

### Changed

- `test_warm_log.py` matches the closing lines of a warm-up run ("finished after 0 ms", "ended by Ctrl+C after 12.3 s") in any of the one duration forms instead of only `m:ss`; OPS.173 had left it on the old form.

generation-output: unchanged
