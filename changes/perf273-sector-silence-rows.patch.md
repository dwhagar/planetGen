**Headline:** the last recorded silences of the sector commands are gone: a sector's final passes over its bodies no longer stand 4 to 6 s with nothing moving.

### Fixed

- PERF.273: `generate-start-sector`, `generate-neighborhood` and `generate-sectors` stood 4 to 6.6 s on some runs of the silence test. No one solve was long (a sweep of 400 random bodies: median 11 ms, worst 1.2 s for an airless world; a world's seasonal solve in a system 0.25 to 0.7 s); the stretch was a world's final pass (its seasonal temperatures, class and life data, 0.25 to 0.7 s a body) run over every body of a system with nothing telling the unit a body had ended. OPS.265 had the pass start the quiet clock again as each body ends (`_finish_planet`, `_finish_moon`); with that, three repeat runs of the three commands show no stretch over 2 s, and their three 12 s rows are deleted from `RECORDED_SILENCE`, which is now empty. The galaxy is the same.
