**Headline:** The generation-mix test draws its 250 systems from a fixed seed and bounds the giants a sample holds per star, so it no longer fails on main by chance.

### Changed

- `tests/test_generation_mix.py`: the sample is a pure function of `SAMPLE_SEED` (`deterministic_entropy`). Eleven unseeded samples of the same code on 2026-11-09 held 26 to 58 giants (0.08 to 0.17 a star), and the flat bound of 50 giants and the 5 percent super-Jupiter floor failed on main about one run in three. Giants made per 200 systems were 40 before Train 18 and 45 after, so no generation change lies behind it.
- The giant bound is per star, as the registry states the rate (PA2: 0.045 to 0.16 a star): over 0.045 and under a third more than 0.16 of the sample's stars.
- The median density ceiling is Neptune's 1.64 plus 0.4: the median of a sample of 40 to 50 giants swings by about 0.2 (1.41, 1.69 and 1.84 on three samples, the same before Train 18 as after). The generated median sitting above the Solar System's 1.3 is left to the planet-type owners.

generation-affects: none
