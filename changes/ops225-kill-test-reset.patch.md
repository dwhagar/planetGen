**Headline:** Kill test for the galaxy reset (OPS.225).

### Added

- A kill test for `reset-galaxy`: the wipe is killed once before or after a table or the tidying at drawn steps, and the rerun leaves every content table empty with the Object ID counters kept, as an undisturbed reset does. No double write was found.
