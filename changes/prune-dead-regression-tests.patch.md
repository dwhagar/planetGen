**Headline:** Delete regression tests that guard nothing and the unused galaxy shape visualizer.

### Removed
- A skipped test module for a function that no longer exists, and the 18 tests that only check what the first two already show.
- About 20 tombstone tests that assert a retired constant, option or page is gone, and trimmed four that also checked a live behaviour.
- A duplicate test, five unused test helpers, and the unused imports and variables across the test suite.
- The `galaxy_shape_visualizer_cli.py` test tool and the matplotlib test extra it needed.
