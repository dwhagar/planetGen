### Changed

- The sector page's Sector Map is now the Galaxy Map's own engine locked to
  that sector (MAP.68, with MAP.67's one URL scheme): the same stars,
  clouds and bodies, hover, tooltip, ring, info panel and NAV links, with
  zoom buttons, Reset view, "Mark rogue planets" and the Contents table's
  "Show on map" buttons as before. It has no steps or history of its own, and a
  click on a neighboring sector opens that sector's page. Because it is
  the Galaxy Map, its scale line now reads in sectors and parsecs, the
  neighbors are the galaxy's own, and the browser's tile cache is shared
  with the Galaxy Map.
- The old Sector Map code (`sectormap.js`, `render_map_panel` and its
  no-script link list) is gone; `starmap.py` now only builds the scene data,
  `GET /sector/<id>/scene`. A sector with no place in the galaxy has no map.
