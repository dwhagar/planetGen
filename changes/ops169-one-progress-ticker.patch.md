**Headline:** Every long loop that says where it is now does it through one progress ticker at one cadence (a line every 30 seconds), with the same wording everywhere, instead of five hand-made clocks that each said it their own way.

### Changed
- **One ticker (OPS.169).** New `util/ticker.py`: `Ticker(label, total)` says `Label: 1,240 of 4,000 done (31%), 41.3 per s, 0:30 elapsed.` at most every `tuning.PROGRESS_LOG_SECONDS`; `Ticker(..., percent=True)` is for a counter on a bar's scale; `Pulse` is a thread that speaks on the cadence for work that goes quiet between counts; `stamp()` is the one form of the hand `[12 s into the stage]` stamp, now `(0:12 elapsed)`.
- **The copies are gone.** `block_lights.PROGRESS_SAY_SECONDS`, `uncharted_index.SAY_SECONDS` and the `refresh_block_lights` and `BlockLightsFollower` hand throttles and thread are replaced by the ticker. The download lines (`datasets.fetch.PROGRESS_SECONDS`, was 10 s) and the warm-up lines (`warm_progress.PROGRESS_SECONDS`, was 60 s, still held back by ten percentage points) now use the same `tuning.PROGRESS_LOG_SECONDS`.
- The lines read `Block light totals, units: 3 of 40 done (8%), 12,000 rows written, 0:41 elapsed.` and `Uncharted sectors, shells: ...` in place of the old per-place wording; the block light stage's pulse reads `Making the block light totals: 37% done, still working, 4:07 elapsed.`

generation-output: unchanged
generation-affects: none
