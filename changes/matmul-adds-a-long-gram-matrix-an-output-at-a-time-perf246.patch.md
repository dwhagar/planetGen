**Headline:** A planet draw that took 71 s (a tiny body at 0.9 AU of an A0V star, 0.00047 Earth masses with its moons) takes 16 s, because the sum behind the seasonal swing's Anderson mixing is added faster with the same bits.

### Changed

- `src/planetgen/util/detmath.py` (PERF.246): `matmul` of a product too large to form at once (more than `_ACCUMULATE_LIMIT` products) went through a loop of `inner` array operations. The seasonal swing's `_anderson` takes the Gram matrix of 6 series of 138,240 to 190,080 values every pass (31 solves of 30 passes in the draw), 0.26 s a call, 55 s of a 64 s draw. When the outputs (rows times columns) are fewer than the inner terms, each output is now added alone: its products in one buffer and `numpy.add.accumulate` over it, left to right as the loop adds them, so the bits are the loop's. The call takes 0.027 s. Census draw 72 of A0V hot: 71 s to 16 s; the draws the class census abandons at 45 s were this (A0V hot draws 72, 136; K0III hot 176, 181, 231, 232 took 17 to 71 s).
- `src/tests/test_detmath_matmul.py`: 38 tests that the output-at-a-time path is bit for bit the loop's (random shapes and limits, a tiny term lost to the sum before it, a 6 by 138,240 Gram matrix).

generation-output: unchanged
