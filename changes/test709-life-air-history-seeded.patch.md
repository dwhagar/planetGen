**Headline:** The life and air history API test draws the same systems every run, so it no longer fails now and then.

### Changed
- **The life and air history test is seeded (TEST.709).** Its helper drew up to forty unseeded systems to find an airy rocky planet, and about one run in five came out with a planet whose air history has no methane row. Each attempt now draws on a stream of its own, so every run reads the same planet. Only test code changes; it passed 50 runs in a row.

generation-output: unchanged
