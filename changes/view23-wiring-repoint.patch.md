**Headline:** The wiring register no longer names VIEW.23 as the item that wires the sky photometry and light-time modules; VIEW.31 (the sky route) does.

### Changed

- `docs/design/research/wiring-register.md` (VIEW.23): the lines for `planetgen.physics.light_travel`, `planetgen.physics.sky_photometry` and rows SKP1 to SKP4 of 6.171 point to VIEW.31, which is where the sky catalogue they feed is read by product code. generation-output: unchanged
