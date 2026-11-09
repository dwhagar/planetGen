**Headline:** A hypergiant system is bounded by a test, the belt-and-habitable-world cases of the heaviest stars no longer pass the 900 s limit, and the seasonal swing's cheaper fixed point has a byte-for-byte proof.

### Changed

- **PERF.249 part 2: tests that bound the hypergiant cost.** Part 1 (PR 5634) removed the 152 million `array_split` calls and PERF.255 (PR 5787) keeps the Anderson Gram matrix between passes; with both, the K00 hypergiant system of the first bug hunt (seed 0, every orbit flag forced) takes 28 s and 0.35 GB, where it took 311 s and 8.4 GB. `tests/test_hypergiant_time.py` bounds that system at 240 s, and `tests/test_seasonal_swing.py` compares the Anderson fixed point, its passes and its convergence byte for byte with the old form that built the whole Gram matrix every pass, on six stiff maps that wrap the six slots.
- **The M3IA+ and M50 belt-and-habitable-world cases build two systems, like O5V.** Both took more than 900 s at fifteen systems (about 100 s a system beside other jobs). The test also lost its assertion that fifteen requests gave a habitable world at least once: since GEN.1397 a request has one in 0.069 of cases (`tuning.PREVALENCE_FORCED_HIT_RATES`), so fifteen had none a third of the time (M5V failed). The odds are the arithmetic of `tests/test_prevalence.py` and the count of `measure_habitable_count`. No code changes.

generation-output: unchanged
