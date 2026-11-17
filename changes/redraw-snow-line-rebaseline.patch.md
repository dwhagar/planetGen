**Headline:** The redraw-rate test now expects the snow line's scatter to redraw at the rate its 2-sigma cut gives, not the 3-sigma rate of before.

### Fixed
- **`snow_line.scatter` was held to a ceiling measured at the 3-sigma cut (GEN.1611).** GEN.1611 cut the snow line's scatter at 2 sigma, so a draw is thrown away 4.55 percent of the time instead of 0.27, 17.6 times as often per star; `test_no_redraw_rate_is_above_its_registered_ceiling` read the old ceilings and failed (0.4167 redraws per system in `stars-5` against 0.025). The cause is the cut and not a fault: 150 default systems threw away 13 of 210 snow-line draws (0.062 per accepted draw; 0.048 expected from the cut). The nine groups' baselines for the kind are the old ones times 17.6 and the ceilings twice them (`tests/data/redraw_rate_expectations.json`). Test data only.

generation-affects: none

generation-output: unchanged
