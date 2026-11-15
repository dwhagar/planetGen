**Headline:** A unit run again after a passing failure goes on after the pieces it finished, and a test now proves it.

### Changed
- `commands.py` says that `pieces(args, cursor)` takes the ledger cursor, not 0.
- `test_failure_kinds.py` fails a piece at the start, middle and end, and checks each finished piece is done once and the result is whole.

generation-output: unchanged
