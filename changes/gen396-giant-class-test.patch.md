**Headline:** A test now holds that the giants of a one to two solar mass star are class III whatever their light and metallicity, closing GEN.396 (the old class III luminosity range is long gone).

### Changed

- GEN.396 asked for the class III range to follow the evolved star's own luminosity and mass instead of a hand-set 50 to 1,000 solar luminosities. The range was already gone: a giant's class comes from its track phase (`stellar_evolution._evolved_state`), and `YERKES_LUMINOSITY_RANGES` is removed. What was missing was the item's test: `test_post_main_sequence` now sweeps ten cases (1.0 to 2.0 solar masses, Fe/H from -1.0 to +0.3) along the whole track and asserts the class runs main sequence, subgiant, giant, and that a giant lighter than 50 solar luminosities is still class III. No generation output changes.
