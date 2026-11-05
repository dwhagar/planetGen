**Headline:** The usual shares the Generate page and the prevalence options start from now match what generation gives after train 10: asteroid belts in 0.65 of systems, not 0.04.

### Changed
- **Prevalence base shares remeasured (TEST.579).** `tuning.PREVALENCE_BASE_SHARES` was measured before the belt mass cut was lifted (GEN.1068, GEN.1156), so it said an asteroid belt is in 0.039 of systems where generation now gives 0.653 (a real change, not a stale draw: four seeds of 1,000 systems each read 0.635 to 0.663). All four were remeasured over 4,000 systems after train 10: habitable_world 0.584 (was 0.598), asteroid_belt 0.653 (was 0.039), large_star 0.056 (was 0.050), planets 0.822 (was 0.809). A run with `--prevalence asteroid_belt=+30` now scales the belt chance from 0.653, and the Generate page shows the new usual shares. A run with no prevalence option is unchanged. `test_measured_base_shares_still_hold` and the belt cases of `test_prevalence.py` pass; the case `habitable_world -50` fails until GEN.1275 merges (a system that bars habitable worlds still holds a life-class body, which GEN.1275 removes). This branch includes PERF.236's sample sizes.

generation-output: unchanged

TODO items resolved: TEST.579.
