**Headline:** Every place that quietly goes on past an error now leaves a debug line saying what it went past.

### Changed
- **`log.swallowed` in every broad `except` (OPS.168).** 130 handlers in `queue/`, `web/`, `api/`, `db/`, `util/`, `cli/`, `generation/`, `admin/` and others caught every exception and wrote nothing, so a run with debugging on was blind exactly where a fault hid. Each now writes one debug line (the place, what it was doing, the error) and nothing at normal level; behavior is unchanged. A new test fails on any broad `except` that neither re-raises nor logs. generation-output: unchanged
