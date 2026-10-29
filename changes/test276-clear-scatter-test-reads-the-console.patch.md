**Headline:** The clear-scatter test no longer fails when an earlier test left the logger silent.

### Changed
- **The test that reads the console sets the console first.** `test_clearing_a_scatter_with_nothing_stored_says_so_and_changes_nothing` read `capsys.out` but the logger belongs to the test process: a test that ran `main()` with `--quiet`, `--silent` or `--output -` (or the web app, which has no console handler) left it silent or writing to stderr, so under a parallel gate the message came back as an empty string. The module now sets the logger to write to stdout before each test. Tests only; the product is unchanged.
