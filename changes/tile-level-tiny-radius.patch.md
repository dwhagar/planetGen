### Fixed

- The Galaxy Map's tile-level helper (`galaxyViewport.tile_level_for_view_radius`
  and its copy in `lib/galaxymap3d.py`) no longer crashes with
  `OverflowError` on a subnormal view radius; any tiny positive radius gets
  the finest level (MAP.90).
