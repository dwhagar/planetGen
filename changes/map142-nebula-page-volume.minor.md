**Headline:** A nebula's own page now draws it as a soft, irregular cloud with fading edges, like the maps do.

### Changed
- The nebula page's 3D view draws the nebula from its own cloud of points (the same points and layer as the Galaxy Map and the Sector Map) instead of a translucent surface, so it fades unevenly toward its edge and reads as gas. It still turns slowly, zooms and re-centers, with the brightest stars of the galaxy dimmed round it. A nebula that was never placed in the galaxy keeps its still drawing.
- `GET /api/nebulae/<id>/surroundings` now also says where the nebula is (`center_pc`, parsecs in the galaxy frame; null for one never placed).

generation-output: unchanged
