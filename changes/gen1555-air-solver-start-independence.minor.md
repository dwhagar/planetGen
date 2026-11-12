**Headline:** the air solver no longer lets the start decide a giant's temperature or an air's hydrogen-rich or nitrogen-rich regime; the bodies whose answer depends on the start fall from 15 to 6 percent (the rest truly have two answers).

### Changed

- A giant's temperature is iterated to the solver's tolerance (up to `MAX_GIANT_STEPS` 24) instead of taking one step, so a hot giant no longer ends at 1,381 K from a 120 K start and 1,633 K from a 600 K one.
- The recipe regimes (hydrogen-rich, nitrogen-rich, and the rest) are kept once a step is near its answer (`FREEZE_NEAR_K`, `FREEZE_NEAR_LN_P`, earliest step 8, latest 30) and not at a fixed step 8, which drew them from a start-dependent transient.
- G2V census of 60 systems: bodies whose answer depends on the start 42 of 286 (15 percent) before, 22 of 285 (8 percent) after. Bodies with two real answers (snowball and warm) stay two-valued: the production start is fixed. The 6 bodies at the 375 K moist limit are the clamped fixed point under their own drawn constants, not a fault of the solver (a probe with other design draws finds cold answers); they stay.
- Registry rows AS2 and AS3 in `gap-filling-plan.md` (GEN.725); the giant test of `test_surface_air` reads the haze at the temperature the giant ends at.
- Changes the temperature of giants and the air regime of some rocky bodies, so it is a lock item: the golden and epoch are the train's.

generation-output: changed
generation-affects: fill, population
