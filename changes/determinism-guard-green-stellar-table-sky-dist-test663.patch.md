**Headline:** The check that no maths in the package depends on the library version passes again: the stellar table and the sky files use the package's own distance and matrix product.

### Changed

- `src/planetgen/galaxy/stellar_table.py` (TEST.663): the four `@` products (the Legendre projection of the residual and the three table products) are `detmath.matmul`, which adds left to right with no BLAS blocking.
- `src/planetgen/db/sky_tier_ledger.py`, `src/planetgen/generation/sky_charting.py`: the three `math.dist` calls are `detmath.dist`. `test_no_version_or_library_dependent_maths_in_the_package` was red on main with these seven lines.
- Checked: `test_reproducible_draws.py`, `test_stellar_table.py`, `test_sky_charting.py`, `test_sky_tiers.py`, `test_system_sky_route.py` and `test_golden_seed.py` pass (113 passed, 1 skipped, 12 min 36 s on 4 vCPUs).

generation-output: unchanged
