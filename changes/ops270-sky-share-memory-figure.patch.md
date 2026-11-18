**Headline:** The bright-star backfill (Stage 9) now runs four sky shares at once on a 16 GB host instead of three, because its memory figure was measured.

### Fixed

- OPS.270: the backfill claimed 2,100 MB for each unit, the figure of the other scatters, which was measured when a worker built its sampling tables cold. A share forked from a warm template reaches 1,502 MB of its own memory in about 2.5 minutes and stays there (the template holds 175 MB, the sampling tables read back from the cache are 7 MB); the memory is the stellar-evolution caches filling to their size. The claim is now 1,650 MB, so half of a 16 GB host holds four shares, not three.
