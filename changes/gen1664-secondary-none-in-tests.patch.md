**Headline:** Tests that read the second star's bodies no longer fail on single stars and close binaries, where there is no second star list.

### Fixed

- After the hosts rename (GEN.1664), `StarSystem.secondary` is `None` unless a wide pair has one; the tests that read `system.secondary.orbiters` now add those bodies only when there is a secondary host. No assertion is loosened and no generated value changes.

generation-output: unchanged
