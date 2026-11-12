**Headline:** The line for a run that stopped early now says why and at which stage, and a cancel is no longer logged as an error.

### Changed
- **The run-failure line names the cause and the stage (OPS.163).** `Run … stopped at stage 3 of 5: ValueError: …`, with the run id, stage and cause as fields. A cancel (Ctrl+C or SIGTERM) is logged as a notice. generation-output: unchanged
