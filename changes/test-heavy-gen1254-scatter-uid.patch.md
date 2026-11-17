**Headline:** The scatter-id test file is listed as heavy, so the quick test slice leaves it out and a run on a busy box is not ended by the hang guard.

### Changed
- **`test_gen1254_scatter_uid.py` joins `HEAVY_FILES`** (`src/tests/areas.py`). It builds charted sectors, and a system build is now about 2.4 times slower (GEN.1779, PERF.276): Lane 14 measured the file at 1,495 s alone on a quiet core against 617 s on record, over the 1,140 s hang guard. It is not a hang. The test and the guard stay as they are.

generation-output: unchanged
generation-affects: none
