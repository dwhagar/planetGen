**Headline:** Dice can now be keyed to a body's birth facts, and a new test names the end-stage streams that still move when a body is nudged.

### Changed

- GEN.1547 (built ahead for Lane 05): `draw.body_stream` and `draw.body_stream_of` (mass, formation distance, star mass); `tidal_q_prime` uses them with the same string, so no draw changes; `test_end_stage_replay.py` replays the end stages on a body moved by a hair (strict xfail until the call sites are switched).

generation-output: unchanged
generation-affects: none
