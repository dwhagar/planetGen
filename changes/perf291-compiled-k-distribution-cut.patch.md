**Headline:** The cut of a k-distribution at the g edges, one fifth of the time a wide binary's planets take to generate, is now compiled, with the same bits as before.

### Changed

- PERF.291: `datasets.ir_ktables.reduce_distribution` (read by every infrared optical depth: `line_opacity.layer_depth`, `resort`, the air history's starlight aloft) hands the cut to `physics/ktable_kernel.py`, a compiled loop over the order numpy's own `argsort` gives, adding in the same order in float64 with no fast-math. On 400 real calls from a generated wide binary the answers are equal to the last bit (`array_equal`), and `tests/test_ktable_kernel.py` holds that on random rows with ties, runs of zeros, float32 rows and any leading shape. The cut took 17 s of 54 s for one wide G2V binary with comets on (1,416 calls, 12 ms each) and now takes 2.8 s; the system takes 39.8 s. The comets were never the slow step (their draws cost 0.01 s); the time is the GEN.1745 air history's anchor solves. The kernel module follows `radiative_kernel`: it is imported when the cut is called, so the web app never loads numba (`test_the_web_app_never_imports_numba`).

generation-output: unchanged
