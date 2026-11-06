**Headline:** The rogue yield per star is now checked against the microlensing census and the terrestrial scale.

### Added
- **TEST.499: a calibration test for the rogue yield.** `tests/test_rogue_yield_calibration.py` measures the bodies the planet architecture ejects per star (about 0.11 of 0.33 Earth masses or more, 0.085 of 0.1 to 10 Earth masses) against the microlensing band of 8 to 42 per star and the terrestrial band of 0.8 to 7.5, and checks that the Rogues level ceilings reach the bands. The donor check is a strict expected failure: the architecture alone is far below both bands and the top-up makes up the rest. The numbers are in `reports/rogue-yield-calibration-2026-11-06.md`.
