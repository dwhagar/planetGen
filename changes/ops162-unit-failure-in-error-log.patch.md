**Headline:** A failed unit now reaches error.log with its run, stage, unit and job kind beside the traceback.

### Changed
- **A failed unit's line in error.log carries its context (OPS.162).** It names the unit and the stage and has the run id, stage, unit and job kind as fields, so a failure can be found by run. generation-output: unchanged
