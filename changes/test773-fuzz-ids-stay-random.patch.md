**Headline:** Two galaxy fuzz tests no longer fail in a whole-file run: the test harness stopped pinning the random ids of runs and tickets.

### Fixed
- `test_fuzz_cli` pinned `secrets.token_bytes` for every length to make a new galaxy's seed repeatable; `secrets.token_hex` calls it too, so a run's id and its ticket in the line repeated from one `main` call to the next and a new run was handed the last run's finished verdict ("the run done before stage 1 was queued"). Only the galaxy seed's length is pinned now.

generation-output: unchanged
