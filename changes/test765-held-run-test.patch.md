**Headline:** A queue-commands test now expects the held run that a unit failing for good leaves since OPS.221.

### Fixed
- `test_queue_commands.py` expected the run of a failed unit to be `failed`; OPS.221 holds it for a person (`held`, with the unit in its error), so the test says that.

generation-output: unchanged
