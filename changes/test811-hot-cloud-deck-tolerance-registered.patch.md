**Headline:** The stopping test of the giant cloud-deck bisection is a registered stand-in, so `test_constants_completeness` is green again.

### Fixed

- `hot_cloud_deck.TOLERANCE` (1e-3 of a decade in pressure, the stopping test of the bisection of a giant's profile, GEN.1577) had no registry row, tag or place in the constants module, which turned `test_constants_completeness` red on main from Train 51. It is row DA3 of section 6.317 (a numerical stopping test, tier T9, tweak L0, like belt_ledger's BR21 tolerance) and its docstring carries the tag (TEST.811).

generation-output: unchanged
