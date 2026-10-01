### Changed
- Tests: a test that runs past 10 minutes now fails with a stack trace of where it was stuck (pytest-timeout, in the `test` extra), and CI's test jobs stop after an hour and list the 25 slowest tests.
