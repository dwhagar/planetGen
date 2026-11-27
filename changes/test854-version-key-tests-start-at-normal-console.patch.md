**Headline:** The version-key tests start each test with the console log at normal severity, so a test that asked for `--debug` can no longer put SQL lines in front of the next test's run line

### Fixed

- TEST.854: `tests/test_version_key.py::test_a_plan_names_the_seed_it_was_given_or_none_yet` failed once in a three-worker gate run of Train 64 (the first line it read was a `DEBUG ... [sql] ... SELECT DATABASE()` line, not the run line) and passed alone. A run configures the process's logger, so `test_the_line_is_first_in_a_debug_file_too` (which runs `--debug FILE`) leaves the console at debug severity with the SQL trace on in the worker that ran it. The next test of that worker that reads `capsys` then sees the lines written before its own run configures the logger (the helper `_generate` opens the control schema first, TEST.850), at debug severity, ahead of the run line. With one worker the order of the file keeps the debug test last of the readers, so it passed alone. Each test of the file now starts and ends with `log.configure(log.NORMAL)`, which clears the handlers and the debug file. Nothing in the product changes: the run line is first because the command's own `configure` runs before it (`generate.main`, behind `log.untraced()`). Tests only.

generation-output: unchanged
