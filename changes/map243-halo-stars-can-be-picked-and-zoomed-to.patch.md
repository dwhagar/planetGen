**Headline:** Stars and objects beyond the galaxy's outline (the halo, the satellite galaxies) can now be picked on the whole-galaxy view and zoomed into with the wheel, down to a window of their own.

### Changed
- A window of the map and a sector's address now hold every cell the grid names (out to the designation's last shell, all bands), not only those inside the galaxy's density outline (`REACH_OUTLINE` in `galaxyprisms.js`; `windowView`, `sectorStage` and `sectorAtPoint` in `galaxystages.js` lose their outline argument). The whole-galaxy view still draws the real outline's blocks and the block totals still come from the stored outline. Before, the wheel stopped at 6,912 pc over a halo star and a click read "Block 0·0·0, 0 sectors charted", because no block or window existed there.
- The wheel on the whole-galaxy view hands off, beyond the outline, to the window round the cell that holds the star or phenomenon in front (else the point being zoomed toward), and on through the windows (81, 27, 9, 3 cells) as it does for a block of the disc.
- `docs/design/fly-through-view-distance.md` says what is built.

### Added
- A `halo` pick layer on the whole-galaxy view: a star beyond the outline is picked within 10 px (its position, brightness and distance from the core shown), a phenomenon there as inside a window.
- `canvasEl.galaxyOutsideOutline` for the browser tests, the `halo_site` fixture (a loose group of stars and a black hole beyond the outline) and two browser tests: a black hole in the halo is picked, and a star in the halo is picked and the wheel carries the view down through four windows round it.

TODO items resolved: MAP.243 (and MAP.246 folded in) / Added: none
