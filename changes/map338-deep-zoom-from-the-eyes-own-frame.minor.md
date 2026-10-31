**Headline:** The Galaxy Map zooms in to 30 km from an opened system (it stopped at 3,000 km), drawing the system from the camera's own frame so it holds steady.

### Changed
- **The nearest view is 1e-12 pc (MAP.338).** `GLOBAL_FLOOR_PC` goes from 1e-10 to 1e-12. The zoom slider, the ruler and the wheel reach it, and the ruler has a Surface stop at it.
- **The open system is drawn from the eye's own frame once the camera is within 1e-6 pc of its target.** `galaxymap3d.js` keeps the camera's exact offset from its target (`eyeBack`) beside the position (a double whose step is 1.8e-12 pc at 8 kpc), draws the map without the system, then draws the system over it with the eye at the origin and everything scaled so the view is about a unit across (the log depth buffer keeps its resolution there). `galaxysystem.js` gains `frame()` and `galaxystageview.js` `systemFrame()` for it; `galaxyDeepFrame()` is the test hook.
- **The address drops `sector=` once the camera is far from the sector (TEST.496).** After the wheel carried the camera out of an opened sector the address kept `?sector=...` with the camera at the whole galaxy; the address now names the sector only while the camera is within 81 sector edges of its target.

No generator output changes.

generation-output: unchanged

TODO items resolved: TEST.496 (the wheel-out leg), MAP.338 (near and far count as done earlier; this is the frame plus offset and the 1e-12 pc floor; it cannot be proven on a real GPU from the lane, only on the headless browser).
