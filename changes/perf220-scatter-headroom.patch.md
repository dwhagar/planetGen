**Headline:** The star scatter stages now leave half the machine's memory free, so a big run slows down instead of taking the database down with it.

### Changed

- The cap on the workers of the massive-star, bright-star, backfill and band scatter stages (PERF.220) is sized for a worker at its peak, 900 MB (Boss's kernel log showed 450 to 900 MB resident each), and for half the host's memory, leaving the rest to MySQL, Redis and the web server: 8 workers on a 16 GB host, was 15. MySQL was the first process the kernel killed on 2026-10-31. generation-output: unchanged
