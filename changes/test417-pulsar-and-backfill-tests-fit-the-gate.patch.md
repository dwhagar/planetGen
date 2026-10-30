**Headline:** Two of the slowest tests finish in minutes (the pulsar fill test from 93 minutes to 73 seconds, the single-sector backfill test from 15 minutes to 2).

### Fixed

- **The pulsar fill test no longer draws the toy galaxy's thousands of first-fill stars (TEST.417).** `test_a_saved_fill_stores_the_pulsar_systems_and_marks_every_scatter_row_built` filled the busiest cell with the first fill (GEN.633) on, as the two other scatter tests already skip it (`first_fill.without_the_deepest_floor`); it took 5586 s and now takes 73 s. It still checks that the pulsar systems are stored and that every scatter row is marked built.
- **The single-sector backfill test fills shell 30, not the dense shell 4 (TEST.417).** The two sector fills took 15 minutes near the centre and take two minutes there. It still checks that a sector filled with `backfill=False` leaves the backfill to the end of the run and that the next one calls it once.
