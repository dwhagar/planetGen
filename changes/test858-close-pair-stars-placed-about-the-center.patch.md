**Headline:** A close pair's two stars are placed about the pair's center in the galaxy again, and the test that checks it reads fixed seeds (TEST.858).

### Fixed
- A close pair's primary and secondary stars were both placed at the system's center since the pair's orbit moved into its node (GEN.640 slice c): `system_position` still read the dropped `binary_*_position_*` fields of the pair's host and fell back to zero. It now reads each star's place from the pair's node, as it already did for a wide pair, so the stars straddle the center, mass weighted.
- `test_spatial_position.py` drew random systems until one was a close pair, so it passed or failed by the draw. It now draws a close and a wide pair on fixed seeds and checks that the stars sit about the center.

generation-affects: none
generation-output: unchanged
