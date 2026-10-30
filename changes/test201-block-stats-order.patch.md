**Headline:** New tests show the stored block totals do not depend on the order the sectors are added in.

### Changed

- TEST.201: tests/test_block_stats_order.py adds the same generated sectors to the block sums in several random orders (`block_stats.accumulate`), through a rebuild on a database whose sectors were inserted in each order, and as one block of 3000 terms spanning six decades. The counts must be identical and the real sums (expected systems, age, luminosity, mass) within a relative 1e-12, the error of a 64-bit sum of this many terms. All pass, so the totals are order independent within that bound; compensated sums stay with GEN.109. Tests only: no stored value changes.
