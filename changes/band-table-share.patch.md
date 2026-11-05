**Headline:** The slow seed test takes seconds, and lane machines can share the star-table cache

### Changed

- `test_a_new_seed_replaces_the_old_only_before_any_sector` puts one sector in the galaxy grid directly instead of generating one, so it takes about 13 s where it took over four minutes on a fresh machine; it still checks that a different seed is refused once a sector exists, that the stored seed stays, and that the same seed plans again (PERF.232).
- New `scripts/band_table_share.py` copies the kept star sampling tables to a shared folder under a key of the code's contents and back, so a fresh lane machine skips the 100 s table build when the code is the same; a copy made by other code is never read. Production does not use it (PERF.232).

generation-output: unchanged
