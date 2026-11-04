**Headline:** The nav chart and the System Map diagram stay inside the browser window instead of running past its bottom.

### Changed
- The nav chart (5:3) and the System Map diagram (1:1) were as tall as their column was wide, so in a wide or short window they ran past the bottom of the window and the page had to be scrolled to see all of the map. They now take at most 95 percent of what the window leaves under the header (never under the least height, 320 px), the way the Galaxy Map, the Sector Map and the 3D system view already do, keep their shape (the nav chart's labels are placed in percent of it), sit in the middle of their column, and re-fit when the window is resized, with no reload. The diagram is fitted when it is brought up from behind the Diagram/3D switch too. New `static/boxfit.js` (`mapfit.fitWidth`, `mapfit.watchWindowBox`); a box asks for it with `data-fit-aspect="W/H"`.

generation-output: unchanged
TODO items resolved: none (UX.142 slice 1 of 2: the nav chart and the System Map diagram; the admin shell is slice 2)
