**Headline:** Three tests that were red on main pass again, from the repository root and from `src`.

### Changed
- `test_set_permissions_keeps_both_config_files_when_both_exist` ran the permissions script against the machine's real default jobs folder, so a machine that already had `/var/lib/planetGen/jobs` (root's, not empty) got a warning the test forbids. The script runs now name a jobs folder under the test's own folder unless the test names one in its `config.json` (TEST.683).
- `test_the_code_reads_no_variable_the_model_does_not_name` lists `PLANETGEN_WARM_OUTCOME_FILE` (update.sh tells the warm-up where to write that the galaxy was empty) and `PLANETGEN_QUEUE_EXTRA_FUNCTIONS` (a queue worker's own addition to the allowed functions), each with its reason.
- `test_update_sh_reports_a_skipped_warm_up_not_a_failure` read `scripts/deploy-common.sh` relative to the current folder and failed from `src`; it reads it relative to the test file.

generation-output: unchanged
