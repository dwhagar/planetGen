**Headline:** The habitable colour browser test reads what the Galaxy Map drew for each sector instead of a few faint screenshot pixels.

### Changed
- **The sector colour test (TEST.448).** Choosing "Habitable worlds" fills a sector in the ramp's colours faintly and low in the map under the legend lines, so the patch of screenshot round each sector gave the same few levels for both, and for either as often as not. The Galaxy Map now exposes `canvas.galaxyBlockLooks()` (a hook for the browser tests, like `galaxyLines`: each charted block's shell, band, wedge, worlds, `noFill` and opacity) and the test asserts the sector with worlds is filled and the one without is left unfilled.
