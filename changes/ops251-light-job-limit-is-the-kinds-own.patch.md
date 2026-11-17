**Headline:** A site job is now ended by the clock at a multiple of how long its own kind has taken before, not at a typed ten times the estimate.

### Changed
- `lightjobs`: once a kind has eight recorded runs, its time limit is three times its longest (never under 30 s); until then the typed limit holds. What was read of a kind is kept five minutes, and a database that does not answer leaves the typed limit.
- A light job that ends on its own now records its seconds (`jobprofiles.record`); one the clock or a stop request ended does not (it would teach the limit its own value).

### Added
- `jobprofiles.seconds_of`: the wall seconds of a kind's newest finished runs.

generation-output: unchanged
