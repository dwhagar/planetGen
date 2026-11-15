**Headline:** The job engine's fence and the end of a stage are proved by seeded scenarios.

### Changed
- `tests/test_engine_proof.py` draws interleavings of attempts and kills of the end of a stage, ten seeds each, and checks that two attempts never both commit, that a job is done only by its newest attempt, and that a finished stage concludes, finishes and opens the next once however often its end is killed (OPS.225, the first slice).
