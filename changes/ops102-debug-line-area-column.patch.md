**Headline:** The debug log line names its area in brackets, so a debug log can be filtered by it.

### Changed

- The debug line is now `2026-10-23T17:29:42.729-0400 DEBUG    worker.py[2142] [draw] random.py:688 random | roll draw.random() -> 0.31` (OPS.102): the area (`[draw]`, `[sql]`, `[queue]`, `[web]`, `[decision]`, `[cache]`, `[branch]`) replaces the logger-name column, which was almost always `planetgen` (or `planetgen.trace`). The format is the only one written and read: `test_debug_mode.py` reads the new columns, and `docs/reference/cli.md` and the logging standard show them.
- The debug form of the formatter derives the area of a record no filter stamped, so a handler made outside `configure` cannot fail on a missing field.

generation-output: unchanged
