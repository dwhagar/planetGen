**Headline:** A nebula's shape mesh never reaches past the nebula's radius, as the shape's documentation promises.

### Fixed
- Galaxy Map: the scale that makes a nebula's farthest surface point 1.0 is fitted on a 24-cell grid, and the mesh is cut on a finer one (28 cells in full detail), so a thin spike of the surface that the fit grid steps over reached past the bounding sphere (1.006 radii in the draw with seed 14, 1.029 in the one `test_nebula_shape_endpoint_serves_a_mesh` met once in Core 4's gate). `contains` already ends the nebula at the sphere and `radius_ly` is promised to bound it, so the few vertices out there are pulled back onto the sphere. Shapes and scales are unchanged; only a mesh vertex past the sphere moves. A new test pins the seed-14 spike and the API test allows for the five decimals the route rounds to.

TODO items resolved: TEST.253 / Added: none
