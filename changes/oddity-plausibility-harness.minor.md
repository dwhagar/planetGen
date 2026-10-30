**Headline:** A harness forces each weird-world oddity on 200 bodies and checks that every outcome stays physical, including rings; it runs today on stub rows and fails on deliberately broken ones.

### Added
- **The oddity plausibility harness (TEST.490).** `src/tests/oddity_plausibility_cli.py` draws a gated body for a row, applies the row's effect with the row's own stream, and checks the six general invariants of the weird-worlds design (energy, structure, stability, air, class derived and total, no NaN) and the row's own limits. Stub rows pass and nine deliberately broken rows each fail at the invariant they break; the isolation tests show that rolling a row reads nothing from the run's stream at any dial. A module that implements a row registers its effect in `EFFECTS`, and the shipped table is held to the same checks as rows arrive.
- **The ring list (TEST.494).** Ring mass under 1e-4 of the planet, the outer edge inside the Roche limit unless the body names the R4 resonance, an optical depth that is the formula's output, no ring inside 0.1 AU, the ring plane in the equator, and no ring field on a body with no ring; five broken ring rows show each can fail.

No generator output changes.

generation-output: unchanged

TODO items resolved: TEST.490, TEST.494.
