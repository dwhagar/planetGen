**Headline:** The wheel-anchor browser test puts the pointer next to the star in the same call that asks for the anchor.

### Changed
- **The wheel-anchor test no longer leaves a gap between reading the star's spot and asking for the anchor (TEST.237).** The test now reads where the star is on the screen and asks the wheel's anchor a few pixels from it in one browser call, so a camera or layout still settling under load cannot move the star in between, and it reads the spot again before the zoom steps. It proves the same things; a failure now prints the star's spot and the camera.

TODO items resolved: none. Added: none.
