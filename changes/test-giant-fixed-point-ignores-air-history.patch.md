**Headline:** The giant fixed-point test no longer counts the solves of the air history as the solves of a giant.

### Fixed
- GEN.1745: `test_a_giants_temperature_is_a_fixed_point_and_does_not_depend_on_where_the_solve_began` records every `surface_air.solve` of a giant setting. Since the air record (slice 4b) the history solves the air of a transitional world again at the light and heat of 12 ages (`air_history._solved_temperature`, a giant setting that keeps a surface), and those solves were counted as giants of the run: their 1500 K restart reached another fixed point than the first solve on the Train 51 head (1217.6 K against 1190.5 K). The test now leaves the history's solves out and keeps the final air of each giant, which is what it states. The history's own temperatures are checked by `tests/test_air_history.py`.

generation-output: unchanged
generation-affects: none
