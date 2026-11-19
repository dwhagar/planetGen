**Headline:** Kill test for the population and territories passes (OPS.225).

### Added

- A kill test for `run-population` and `run-territories`: units are killed once each at random commits, the same galaxy is populated again from its own stored ids, and the species, polities and owners match an undisturbed pass. No double write was found.
