**Headline:** The roundness engine now knows that a pile of rubble does not turn round just by being big: it rounds only if it took its shape hot, or for an icy pile at the strength of its ice at that temperature.

### Changed
- **GEN.1176.** `physics/roundness.py` gains a `rubble` structure: a Mohr-Coulomb failure stress from a friction angle of 35 degrees and a cohesion of about 100 Pa (the studies of 1950 DA and other fast rotators), a hot strength of 2.5 MPa for a pile that took its shape hot or partly molten, `rubble_strength_pa`, `rubble_round_radius_m`, and `structure` and `heated` options on `roundness_ratio`. A cold pile stays near a roundness ratio of 0.28 at any size. The numbers live in `physics/data/roundness_strengths.json`. Nothing reads the engine yet, so generated worlds do not change.
generation-output: unchanged
