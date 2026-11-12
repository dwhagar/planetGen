**Headline:** Keys, tokens, passwords and credentialed addresses no longer reach any log: every line is masked to its first 4 characters, its last 2 and X's between (`abcdXXXXXXyz`).

### Changed
- **One masking helper, `util/secrets_mask.py` (SEC.48).** It finds a secret by the name beside it (`api_key=`, `"password":`, `--mysql-password`, `Authorization: Bearer`, a log field called `token`), by its place in an address (`mysql://user:password@host`, `redis://:password@host`) and by its value (every secret option of `config.json` and the environment, every variable named like a secret). A value of 6 characters or fewer is all X's.
- **It is applied at every way a line leaves the logging layer:** each record any logger makes (the program's and a library's, message, exception and stack), `logrecord.Record`, every formatter (so `error.log`, the job logs and the Queue stream), the live-log callback, the relay, `record_error` and `write_line`.
- **`tests/test_log_masking.py` fails if** a key, token, password or credentialed address appears unmasked in any of them.

generation-output: unchanged
generation-affects: none
