**Headline:** A hypergiant system is made in 28 s instead of 42 s, and the seasonal swing of every planet is found faster, with the same temperatures bit for bit.

### Changed

- **The seasonal swing keeps the Gram matrix of its fixed point row by row (PERF.249 part 2).** `seasonal_swing._anderson` formed the whole 6 by 6 Gram matrix of the held differences, 36 sums over every series, at every pass; a pass writes one slot, so only that slot's row and column change. It now forms that row (at most 6 sums) and keeps the rest, which has the same bits (an entry is the sum of the products of two rows left to right, whichever other rows are held). On the K00 hypergiant system (seed 0, every orbit flag forced) the system takes 28 s in place of 42 s on a four-core worker, with the same planets and values. `tests/test_seasonal_swing.py` compares the fixed point, passes and result byte for byte with the old whole-matrix form on six stiff maps; `tests/test_hypergiant_time.py` bounds a K00 system at 240 s.

generation-output: unchanged
