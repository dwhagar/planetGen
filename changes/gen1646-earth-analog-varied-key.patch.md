**Headline:** The Earth analog now moves its planet's orbit by a billionth of an AU each trial, so the oddity streams differ trial to trial; its stored baseline is renewed.

### Changed
- **GEN.1646, the analog (TEST.726).** `scripts/earth_analog.py` reused one planet whose oddity key (`oddities.body_key`: star, distance, mass) never changed, so the second climate branch (W20) drew its air from one fixed stream in every trial and pinned a dry, nitrogen-rich inventory (19.4 ppm N, 6.7 times the bulk silicate Earth) on every trial where it fired. Each trial now sets `planet.distance = 1 + i x 1e-9` AU, and `scripts/data/earth_analog_baseline.json` is measured again (800 trials). Nothing in the generator changes.
generation-output: unchanged
