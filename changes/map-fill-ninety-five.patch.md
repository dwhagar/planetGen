**Headline:** The Galaxy Map is now always 95% of the browser window's height, recomputed on every resize.

### Changed

- MAP.350 follow-up: the Galaxy Map's height was the window less 40 pixels; it is now 95% of the window's height (less a header held at the top of the window), with the rest split above and below, and the page still scrolls past or before the map. `mapfit.js` `fillHeight` takes `FILL_SHARE` and `fillBuffer` gives the margin. No generated output changes.
