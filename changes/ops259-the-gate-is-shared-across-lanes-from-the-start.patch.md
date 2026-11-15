**Headline:** The train gate runs on every ready lane from the first file, with no split to agree on first, and its file times now come from real database runs.

### Changed
- `train_gate.py run --claims DIR` has each lane take the next unclaimed test file, longest first, so an idle or fast lane takes more and the lanes finish together (OPS.259). `train_gate.py status` counts what is passed, failed, running and waiting, and `--release-older SECONDS` frees the claim of a lane that died.
- `train_gate.py learn --replace` sets each file's time to its longest pass, and `scripts/data/gate_durations.json` is refreshed from the train 27 and 28 runs: database files such as the fingerprint, old-schema and round-trip tests had been listed at their skipped seconds, not the 20 to 45 minutes they take.

generation-output: unchanged
