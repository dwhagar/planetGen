**Headline:** The air solver's start dependence and unsettled airs are measured on main and bounded by tests: 2.5 percent of air bodies depend on the start (15 percent before), 0 of 35 airs in the test sample are unsettled (10 of 132 before).

### Changed

- GEN.1555: re-measured on main after Train 36 (G2V, 60 systems, 242 bodies with 10 kPa of air; the census script of Lane 14's branch): 6 bodies (2.5 percent) have an answer that depends on the start, against 42 of 286 (15 percent) when the item was filed; 2 sit on the 375 K moist limit (GEN.1579), 1 near water's boiling point. The fixed point already iterates a giant's temperature (GEN.1716) and settles slow creeps (GEN.1718), so Lane 14's open PR 6103 (giant steps and the recipe freeze near the answer) is superseded; nothing in the solver changes here.
- New `tests/test_surface_air_start_independence.py`: the same setting, start and draw state give the same answer twice; a giant ends at one temperature from nine starts (before GEN.1716 a hot giant ended at 1,381 K from 120 K and 1,633 K from 600 K); at most 15 percent of air bodies of the first twelve G2V seeds depend on the start (0 of 17 today). Real bistability (a snowball and a warm answer) stays, as TODO decided on 2026-11-12.
- `test_every_generated_body_stores_the_density_scale_height_and_greenhouse_its_mixture_gives`: the cap on airs it cannot read because they did not settle goes from a fifth to a twentieth of the bodies (0 of 35 today; it was 10 of 132).

generation-output: unchanged
generation-affects: none
