**Headline:** A new test keeps silent sleeps out of the code that jobs run.

### Changed
- **No sleep in command code (OPS.197).** A job that sleeps while it waits shows no activity, so it looks stuck. A syntax-tree test (`tests/test_no_sleep_in_commands.py`) now lists every `sleep` call in the package by its function and fails on any that is not named with its reason: a blocking primitive that declares its wait, a pause under a second between two tries of a transaction, or code that does not run in a unit. An entry that no longer matches a sleep fails as well.

generation-output: unchanged
generation-affects: none
