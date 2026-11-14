**Headline:** The Earth analog test now runs the generator's own steps in the generator's order, on a fresh body each trial, so its Class M rate is the generator's.

### Changed

- TEST.739: `scripts/earth_analog.py` ran the air, the class, the life and the class again on one body reused for every trial. The life step reads the light the dose and score steps write, so the first trial read none and each later one read what the trial before left (one trial ran on 0.0006 W/m2 left by a frozen world). Each trial now takes a fresh copy of the Earth planet and runs the air, the prepared class, day and night, field, dose, colour, score, the thermal history, the final-age refresh, the class, the star's induced field and the life, as `Planet.__init__`, the tail of `planets.generate_spin` and `StarSystem._finish_planet` do.
- `test_earth_analog.py` reads the calls of those three out of the source and holds the script to the same order, so a step added to the generator without the script fails the test. The stored baseline (`scripts/data/earth_analog_baseline.json`) is measured again.
- No generated value moves: this is a measuring script and its test.

generation-output: unchanged
