**Headline:** The Galaxy Map's window now takes the whole width of the page and the height the browser window leaves, so its shape follows the window's.

### Changed
- **The map fills the page (MAP.291).** The Galaxy Map's window was a fixed 4:3 box (square on a phone) no wider than the window was tall, which left bands of empty page either side of it on a wide screen. It is now as wide as the page's column and as tall as the browser window leaves room for under the header and the search and above the toolbar (never under 320 px), with no fixed aspect ratio, and it follows the browser window as it is resized. The picture, the labels and the scale bar all follow the new shape. The sizing is its own small script (`mapfit.js`) so the other viewers can use it when they are made to scale to their windows (UX.142).

TODO items resolved: MAP.291. Added: none.
