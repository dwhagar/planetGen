**Headline:** The test of the update script's References pages step reads the script's eight stages, so `test_references_pages` is green again.

### Fixed

- `test_update_script_builds_the_pages_in_step_3_and_repeats_a_failure_at_the_end` still looked for `log_stage 3 9` after the update dropped the math check as a stage (OPS.286), which made it 8 stages. It reads `log_stage 3 8` and `log_stage 4 8` now (TEST.815). The `Stage 4 of 9` lines in other tests and docs are examples of the log format, not the script, and are left as they are.

generation-output: unchanged
