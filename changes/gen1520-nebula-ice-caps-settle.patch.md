**Headline:** A giant of a planetary nebula's white dwarf no longer stops the build when a trace of vapor cycles over its ice cap.

### Fixed
- `mixture_giant.cap_ices` raised "the ice caps did not settle" on a sulfur-rich layer at 91 K: its S2 is capped at about 1e-39 of the layer, and each of 1,000 rounds took 3e-39 of it while the equilibrium made it again from the rest of the sulfur. The settle test asked only whether the budget could make more than `MIN_CONDENSED` of the species; it now also asks whether the excess over the cap is more than that. An excess under it is not a condensate worth recording, which is what the test already says for a species the budget cannot make. Found by the Train 45 gate in `test_census_nebula_wholesphere` (B0VII central star, seed 2 of the hunt).

generation-output: unchanged
