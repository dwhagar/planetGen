**Headline:** The order-independence test compares a sector's scatter-born systems without their order-dependent Object ID serials, and passes again.

### Changed

- `test_a_sector_comes_out_the_same_whether_the_backfill_reached_it_first_or_not` failed on main (GEN.1302): sector 1 0 5 held the same four scatter-born systems in either order, but a scatter row's serial is its brightness rank within the write that made it (`db/scatter_uids`), so the serials of the stars an earlier backfill wrote (0 and 1) and the ones the first fill added (2 and 3) came out in another order than when the first fill wrote all four at once, and every Object ID below them differed. That is by design (`docs/design/architecture/object-ids.md`, "Scatter rows", now says so). `tests/test_order_independence.py` puts each sector's scatter serials in brightness order (`_canonical_scatter_serials`) before it fingerprints, so everything else a sector holds is still compared. Run-order-independent scatter IDs would be a new item (a serial block for each mass ring, with an epoch). No generated value changes.

generation-output: unchanged
