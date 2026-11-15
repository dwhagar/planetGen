**Headline:** The Earth analog now draws a different air in every trial, and its stored baseline is renewed on main (Class M 30 of 800).

### Fixed

- TEST.752: `scripts/earth_analog.py` moved only `planet.distance` per trial, but the air keys on `formation_distance_au` (GEN.1547), so every trial shared one air draw. It now sets `formation_distance_au` beside `distance`. The stored baseline `scripts/data/earth_analog_baseline.json` is measured again on main: 800 trials, Class M 30 (it was 20), no errors.

generation-output: unchanged
