**Headline:** Every way of going somewhere on the Galaxy Map (Go to, Up, Home, steps, locate, Back and Forward) is now a camera flight, and the camera's distance alone sets the zoom.

### Changed
- **Navigation runs on the camera (MAP.337, step 2).** `galaxystageview.js` no longer switches the stage the map draws when you go somewhere: `go` is gone, replaced by `flyToStage`, which flies the camera to the place (and loads its data, opens its sector or system) and lands there with the same picture a stage hop gave. The cross-fade between stages (`crossTo`, `stepCross`, `finishCross`) is deleted. The zoom slider, plus and minus and the pinch carry the camera in or out across the whole range with no stage hop, and out of an opened sector or system they go on to the whole galaxy.
- The browser hook `galaxyZoomLevel` names the zoom level the camera's distance reaches, for tests.

No generator output changes.

generation-output: unchanged

TODO items resolved: none (MAP.337 stays open for steps 3 and 4).
