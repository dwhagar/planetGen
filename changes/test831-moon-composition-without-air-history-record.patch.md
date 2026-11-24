**Headline:** The moon tests build their systems without the record of the air history, which was most of a system's cost.

### Fixed
- `tests/test_moon_composition.py`: the 60-system origin and orbit test and the seeded-moons fixture build their systems with `air_record.record_of` giving no record (`tests/air_history_off.py`). The record is made at the end of a body's build, from the finished body, on streams of its own: the moons of the seeded systems (60 systems and 978 moons for the origin test; 40 systems and 1,199 moons for the fixture) and the next draw are identical with and without it. The origin test fell from over 900 s (the limit) to 485 s alone.
