**Headline:** The shell log mask no longer takes the word "pass" for a password name, so "Skipping the population pass: the math check failed." prints whole.

### Fixed
- `log_mask` in `scripts/logrecord.sh` (SEC.48) matched the bare word `pass` followed by a colon and replaced the word after it with X's, in the shell gate and in the no-Python fallback. It now matches only `password`, `passwd` and `passphrase`, the same names `util/secrets_mask.py` holds. `test_math_gate::test_update_sh_check_math_failure_only_warns` passes again, and a new test pins that ordinary log text goes through unchanged.

generation-output: unchanged
generation-affects: none
