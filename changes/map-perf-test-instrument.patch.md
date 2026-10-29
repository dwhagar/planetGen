**Headline:** The Galaxy Map's performance test measures a climb step after its fetch has begun and counts the icon files the map now reads.

### Fixed
- **The node-path performance test read a URL as a request.** The map reads each point-of-interest icon file with `fetch(new URL(...))` (MAP.320); the test's counter took `.url` of it, threw inside the map's marker code and left the markers half built, which showed as `indexOf` and `scale` errors on the page. The counter takes either (TEST.439).
- **A climb step could end before the map had asked for its nodes.** The map waits 300 ms (and reads its stored nodes) before it fetches, and the test called a step settled after 400 ms of stillness, so on a busy machine "up 4" ended with a part-filled buffer and the next step was blamed for the rest. A step now waits 1 s, still under the 1.5 s prefetch idle (TEST.439).

TODO items resolved: TEST.439.
