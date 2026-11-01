**Headline:** The life sweep follows the class registry instead of a fixed count of twenty letters.

### Fixed
- **TEST.542: the life sweep expected 20 class letters and the rules give 19, since Class Q was retired (GEN.1135).** `test_sweep_life.py` now checks that the letters tried are exactly the classes of the registry, and its text and the audit's no longer name a literal count.
