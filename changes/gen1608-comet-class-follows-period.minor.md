**Headline:** A comet's period class now agrees with its orbit across the 200-year line that makes it periodic.

### Changed
- GEN.1608 (bug): the period class of an elliptical comet is drawn before its orbit, and the orbit is drawn from eccentricity ranges (cut by the star's orbit ceiling), so a comet drawn `long_period` could have a period of tens of years and get the periodic `P/` designation, which reads the real period. `comet.period_class_of` now moves the class to follow the period across `PERIODIC_COMET_MAX_PERIOD_YEARS` (a `long_period` comet under 200 years is `halley_type`; a `jupiter_family` or `halley_type` comet at 200 years or over is `long_period`), before the nucleus makeup is drawn, so the ices follow the class it ends with. No draw is added or removed. In a run of 1,500 elliptical comets about 10 percent change class (148 for a Sun-like star with no ceiling, 322 for a 0.2 solar-mass host with a 20 AU ceiling); their stored class and nucleus makeup move.
- `tests/test_comet_period_class.py`: the rule for each case, and a census that a drawn comet is `long_period` exactly when it is not periodic and that its designation agrees.

generation-output: changed
generation-affects: fill
TODO items resolved: GEN.1608
