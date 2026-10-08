### Fixed

- A sector page (and the Galaxy page) no longer embeds the tile cache's
  count of cached tiles in its first-frame data, so the page's text stays the
  same between two loads until something changes. It started to differ when
  the sector page took the Galaxy Map's tiles (MAP.68), which failed
  `test_sector_page_and_tiles_show_a_phenomenon_the_cli_added`.
