**Headline:** A test now proves a sector comes out the same whatever order sectors were generated in and at any worker count.

### Changed
- **`tests/test_order_independence.py` (TEST.181).** The same four sectors of shell 1 are generated in two orders (one `galaxy` run a sector, ascending and descending), after a first run that backfills the bright stars around one of them, and at one worker and four, and their fingerprints (GEN.58 over the GEN.304 encoder, a failure naming the sectors and sections that differ) are compared; three sectors whose objects share an ID are saved in two orders and must hold the same rows and the same collision numbers (GEN.307). Reverting the address ranking of GEN.307 fails the shared-ID test.
- A sector generated before the backfill reaches its cell does not get the stars the backfill places there, so it differs from one the backfill reached first; that is not covered by the test and is reported to TODO under GEN.57's backfill clause.

TODO items resolved: TEST.181. Added: none.
