**Headline:** The exobase of an air, solved 17,000 times in the loss histories of one wide binary's planets, is now compiled with the same bits as the Python solve.

### Changed

- PERF.291: `thermosphere._exobase_solved` hands its arithmetic to `physics/exobase_kernel.py`, the same lines compiled. A square written `x ** 2` is CPython's library `pow`, which the compiler folds to `x * x` (a different last bit for about one value in a thousand), so the squares go through `math.pow(x, two)` with the exponent passed in. The answers are equal (`==`) on all 20,326 calls recorded from a generated wide binary, on 6,000 random airs, on a hot air that never closes and a cold one that closes at once (`tests/test_exobase_kernel.py` keeps the Python solve as the reference), and `kept_shares` over 300 real airs is equal. The loss history of an air takes 7.7 ms instead of 10.3 ms. The kernel is imported when the solve is called, so the web app never loads numba.

generation-output: unchanged
