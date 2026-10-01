### Changed
- **The Galaxy Map zooms smoothly.** The blocks for each view are now built
  in a Web Worker (`static/galaxyblocks.js`), so a zoom step no longer
  freezes the page while hundreds of milliseconds of block listing runs.
  Zoom steps glide over 160 ms instead of jumping (they still jump with
  "reduce motion" turned on), and a change of block size crossfades
  instead of popping.
- **Zoom steps you've already seen are instant.** Built views are kept (up
  to about 4 million vertices), and while the map is idle it prepares the
  views and fetches the tiles one zoom step in and out. The first frame is
  still built on the page, and so is everything in a browser where the
  worker can't start.
