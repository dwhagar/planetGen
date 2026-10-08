### Fixed
- The Galaxy Map's browser tests wait until the map has drawn the stage its breadcrumb names (new `galaxyReady` on the canvas) instead of a clock, so the camera, scale line, slab button and hover tests no longer fail under load (TEST.96, TEST.98, TEST.99, TEST.100).
