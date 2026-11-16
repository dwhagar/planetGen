**Headline:** The first planet of a Kepler-like inner system now sits anywhere from a few days to a few hundred days, so planets exist from 30 to 400 days where there were almost none, and K dwarfs have planets in their habitable zone.

### Changed
- GEN.1709: the first planet of an inner system is log-normal about 25 days with a width of 0.7 dex (`INNER_FIRST_PERIOD_MEDIAN_DAYS`, `INNER_FIRST_PERIOD_SIGMA_DEX`, registry row PA19), where it was 5 days with 0.2 dex: small planets are flat per factor of 2 in period from 8 to 128 days ([R2020] Table 2) and the inner system ended by 130 days. `first_distance_au` takes the median and width, so the draw is reused. The starting count of the inner system is 9.0 (5.619 before, `INNER_INITIAL_COUNT_EXTRA`, row PA3) because the chain now ends at the snow line or 400 days (`inner_chain(limit_au=)`; past the snow line a core may take gas, which doubled the M dwarfs' giants within 2.5 AU).
- The giant-impact stage of GEN.1560 starts at the outer chain's start, not at 0.45 AU: it merged the inner system's planets beyond 128 days.
- Measured on 230 Sun-like systems (F7V to K2V, same seeds before and after): see the PR body.
- Generation output changes: every planetary system with an inner system.
