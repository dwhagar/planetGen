**Headline:** The roundness engine can now take a per-body strength drawn from its range, so the line between a round body and a potato is a band instead of a single size.

### Changed
- **GEN.1178.** Every function of `physics/roundness.py` takes an optional `draw` (0 to 1): the quantile of a body in the range of each strength (ice 0.4 to 3 MPa, rock 4 to 30, metal 50 to 250, and the friction angle, cohesion and hot strength of rubble), with the stand-in at the median. `drawn_strengths_pa` gives the strengths at a draw. The round radius moves about 1.5 times either way. The draw itself is not made yet: the call from `draw.py` rides GEN.1159. Nothing reads the engine yet, so generated worlds do not change.
generation-output: unchanged
