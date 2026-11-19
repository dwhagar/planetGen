**Headline:** Kill test for the object index warm-up (OPS.225).

### Added

- A kill test for `warm-object-index`: the unit is killed once before or after a write or sweep of the index at drawn steps, and the rerun leaves an index with every object once and nothing missing, different or extra (`objectindex.check`). No double write was found.
