**Headline:** A gate run now fails, instead of skipping, when the database tests cannot reach MySQL, and a merge conflict in a generated file has a one-command repair.

### Changed

- TEST.723: with `PLANETGEN_REQUIRE_MYSQL` set (every `scripts/test.py` command except `quick` sets it), a test skipped with "No MySQL test server reachable" is reported as a failure carrying the server's message, as the browser tests already are (TEST.693). A refused login or a dead server can no longer come back as a run of skips that looks like a pass.
- TEST.721: new `scripts/regenerate_generated.py` (`--list`, `--check`, rebuild all, `--only`, and `--resolve`, which rebuilds and `git add`s the conflicted generated files after a merge and refuses while a conflict in a real source is left). It covers the stand-in list, the references pages and the data tables built from scripts; the global documents stay the TODO thread's. The procedure is in `docs/reference/testing.md`.

generation-output: unchanged
generation-affects: none
