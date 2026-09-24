### Changed
- **The 3D Galaxy Map shades predicted density with cylindrical segment
  prisms instead of spheres.** Space is cut into rings, wedges and layers
  on the same grid as cylindrical sectors, a power-of-two number of
  sector widths across so the prisms scale with the view (one prism is
  one sector at full zoom). Each prism is solid and lit, colored and
  sized inside its cell by its mean density. The density is computed in
  the browser from the galaxy's own shape (`static/galaxyprisms.js`), so
  the map no longer asks the server for a density point cloud.
