**Headline:** A rare failure of the example narrative test, on a subdwarf that rightly carries the main-sequence lifespan note, no longer happens.

### Fixed
- **The main-sequence note test (TEST.500 follow-up).** `test_example_narrative_never_misapplies_main_sequence_notes` skipped class V stars only, so a Yerkes VI subdwarf companion in the van Maanen's star example (drawn rarely) failed it, though `describe/build.py` writes the main-sequence note for classes V and VI by design and the physics treats a subdwarf as main sequence. The test now skips both. No generator output changes.

generation-output: unchanged

TODO items resolved: TEST.500.
