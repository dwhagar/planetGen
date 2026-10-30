**Headline:** The help text of `planetgen run-job` and the replacement shown for the five ground-type rows no longer carry item codes.

### Fixed
- **TODO codes in text a user sees (TEST.487).** The `run-job` description named PERF.203 and the five stand-in rows of the surface ground types named GEN.474 as the replacement. The help now says what the command does without the code, and the rows say "ground types derived from the bulk composition and the melting of the mantle"; the item code stays in a comment. `test_no_todo_codes_shown` passes. Nothing generated changes.

generation-output: unchanged

TODO items resolved: TEST.487.
