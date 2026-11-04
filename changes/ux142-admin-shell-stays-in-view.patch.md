**Headline:** The admin shell (the rail, the search box and the density switch) stays in view under the site header while a long admin page scrolls.

### Changed
- On a window at least 700 px wide and 600 px tall, the admin shell is held just under the site header, so the rail, the search and the density switch are no longer scrolled away on a long admin page. A phone or a short window keeps it in the page flow, where it would otherwise take too much of the window. `adminshell.js` measures the header and hands its height to the CSS as `--site-header-h`, so the hold follows the header when it wraps or the window is resized.

generation-output: unchanged
TODO items resolved: UX.142 (slice 2 of 2; slice 1, the nav chart and the System Map diagram, is in the earlier change note for the same item)
