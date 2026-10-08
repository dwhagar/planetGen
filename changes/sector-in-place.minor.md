### Changed

- On the Galaxy Map, clicking a generated sector no longer loads its
  page: the sector opens in place as the drill-down's last stage (MAP.66).
  Its stars, black holes, nebulae and other bodies are drawn where the
  sector is in the galaxy, the camera flies to fit it, and hovering or
  clicking one gives the same tooltip, ring, info panel, ☆ and NAV links
  as on the sector page. The breadcrumb, Back, Forward, Up and Reset view
  work as on any other stage, and the stage has its own URL,
  `/galaxy?sector=<designation>&open=1`. Clicking a neighboring sector
  steps sideways into it. A sector's panel keeps its "View sector →"
  link to the page.
- The Sector Map's scene is built by a new module (`static/sectorscene.js`)
  that both the sector page and the Galaxy Map use, and `starmap.py` hands
  the sector's place and size in the galaxy to the script. The scene is
  also available as JSON at `/sector/<id>/scene`.
