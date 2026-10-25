**Headline:** The Galaxy Map's sector blocks draw again: the sectors' colour (such as Color by Habitable worlds) was missing since the one-star-buffer change, and the sector page map failed at phone width.

### Fixed
- **The block prism shader compiles again (MAP.293, MAP.294).** The relative-to-eye position lines of the star buffer (MAP.267) had gone into the block prism vertex shader, which has no such uniforms, instead of the star vertex shader: the blocks' shader did not compile, so no sector was filled or coloured, and the sector page map showed the same error at phone width. The star shader now takes the lines, as the star buffer meant, and a browser test draws the filled blocks at desktop and phone widths and fails on any shader error.
