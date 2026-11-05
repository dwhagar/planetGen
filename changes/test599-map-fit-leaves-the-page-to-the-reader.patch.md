**Headline:** The System and Galaxy views no longer scroll the page away from a key you have just pressed: a Tab to the skip link at the top of the page stays on screen.

### Changed

- The map-fit of a map that fills the window (`mapfit.js`, MAP.350) scrolls the page so the map sits in view until the reader scrolls it. A key press, a mouse or finger press and the wheel now count as the reader taking the page over too, for every map on the page, also one built later (the 3D system view is built when its data has come, which on a slow machine is after the reader's first Tab). Before, a scroll down to the map (544 to 616 px on a phone) followed the first Tab and left the focused skip link above the window (TEST.599: `web.system-phone-light` and `-dark` failed about one run in four).
- Two JavaScript tests for it (they fail without the change), and the skip-link failure in `test_web_a11y.py` now says where the page was.
