**Headline:** The web app no longer pulls in numba through the post-main-sequence module, and the belt members table has its documented precision floor.

### Fixed
- **The web app imported numba (TEST.608).** `physics/post_main_sequence.py` imported the numba bright-piece kernel at the top, and `planetgen.web.app` reaches that module, so Apache's web user (which cannot write numba's cache) would fail on it (OPS.42). The kernel is now imported where it is used, in `bright_pieces`, which only the workers call.

### Changed
- **The belt members table is in the precision floor table (TEST.609).** `belt_members` places a body to 2.4 m like `asteroid_belts` (km columns at 10,000 AU); `tests/test_precision.py` and `unified-unit-quantity.md` section 6.1 list it.
