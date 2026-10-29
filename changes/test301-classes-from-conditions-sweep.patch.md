**Headline:** A sweep now checks that every planet and moon class can be derived from the body's conditions, and that a request for a habitable world is answered.

### Added
- **The classes-from-conditions sweep (TEST.301).** `src/tests/audit/sweep_classes.py` checks, on the shared body sweep, that every planet and moon has a class letter of its own kind, that `classify` is total, repeatable and leaves the body alone, and that the stored letter is the one the conditions give. It also generates 12 systems with `+habitable_world` and reports how many answer the request. Two results are known and tagged in `src/tests/data/class_sweep_expectations.json`: the stored letter is still drawn before the body is built (GEN.646), and most habitable-world requests answer unmet because the air solve of 4c leaves no free oxygen (GEN.855). Both print in a report and fail the moment a tag is removed and the result still falls short.

TODO items resolved: TEST.301 (the reclassification after an edit is left for 8.6). Added: none.
