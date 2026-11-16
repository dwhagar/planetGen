**Headline:** The first planet of a Kepler-like inner system now sits anywhere from a few days to a few hundred days, and the giant-impact stage no longer merges inner-system planets, so planets of 100 to 400 days exist where there were none.

### Changed
- GEN.1709: the first planet of an inner system is log-normal about 20 days with a width of 0.7 dex (`INNER_FIRST_PERIOD_MEDIAN_DAYS`, `INNER_FIRST_PERIOD_SIGMA_DEX`, registry row PA19), where it was 5 days with 0.2 dex; small planets are flat per factor of 2 in period from 8 to 128 days (Hsu 2019 Table 2). The chain ends at the snow line or 400 days, whichever is nearer (`inner_chain(limit_au=)`): past the snow line a core may take gas, which doubled the M dwarfs' giants within 2.5 AU. `first_distance_au` takes the median and width, so the draw is reused.
- The giant-impact stage of GEN.1560 starts at the outer chain's start, not at 0.45 AU: it merged the inner system's planets beyond 128 days.
- Measured on 230 Sun-like systems (F7V to K2V, same seeds): planets of 0.5 to 1.5 Earth radii in 0.36 to 1.1 sunlight per host star 0.222 to 0.230 (K0V to K2V 0.01 to 0.04; F7V to G8V 0.32 to 0.32); planets of 1 to 20 Earth radii inside 400 days per single Sun-like star 1.17 to 0.93 (observed 1.23). Not enough to reach the survey's 0.37 to 0.88.
- Generation output changes: every planet system with an inner system.
