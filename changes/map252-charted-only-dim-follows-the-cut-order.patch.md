**Headline:** "Charted only" now dims uncharted objects in the order they would be cut: satellite-galaxy stars first, then the halo, then the small and dim up to the large and bright, and the farther from the focus the more.

### Changed
- Galaxy Map: under "Charted only" an uncharted star or scattered phenomenon keeps a share of its light by the cut order (Boss, 7:45 PM EDT; MAP.250's tiers and MAP.251's order): a satellite-galaxy star 0.03, a halo star 0.08, then two ramps at the same time from 0.12 to 0.5, stars by luminosity and phenomena by mass (nebulae are left out), so the smallest and dimmest dim hardest and the most massive and the brightest least. The farther an object is from the camera than its focus, the less it keeps (down to 0.35 of its share at four times as far). A charted object is untouched.
- Wire: the tile's flag byte (`charted`) carries the star's tier in bits 1 and 2 (halo 1, satellite galaxy 2); older cached tiles read as tier 0.
- Docs: `docs/design/zoom-star-visibility.md` section 12.
- Tests: `starpoint.test.mjs` (the order, the ramps, the distance rule, the shader wiring), `test_tilepack.py` and `tilepack.test.mjs` (the tier rides the flag byte) and a browser test that reads each object's share from the drawn Galaxy Map on the halo fixture.

TODO items resolved: MAP.252 / Added: none
