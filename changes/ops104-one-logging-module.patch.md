**Headline:** The program's log lines are now made and written by one module, with one formatter, and the shell's log helper takes its words from it.

### Changed

- `util/log.py` is the only place a log line is made and written: its `LineFormatter` replaces the five formatter classes that were spread over `util/log.py` and `util/logrecord.py` (the command line's plain text, the debug line, the file line, the JSON line and the display line), and `log.get(component)` replaces `logrecord.get_logger`. What the lines look like does not change.
- `scripts/logrecord.sh` takes its severities, the ones that go to standard error, the stage word and the divider from the module (a block `scripts/build_log_vocabulary.py` writes), and a test fails when the block is out of date.
- The test that watches for hand-made log lines also fails on a formatter class or a handler made outside `util/log.py`; the activity log's handlers and the queue run log's handler are named in it as stragglers for the next part.
