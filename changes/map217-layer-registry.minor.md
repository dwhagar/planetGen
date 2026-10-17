**Headline:** The three maps share one layer registry, so a map's layers have one menu, one set of URL keys and one legend.

### Changed
- New `maplayers.js`: a layer (id, label, modes, group, default on or off, URL key for hiding and for showing, legend line, picker layer, availability, title) is registered once, and the registry gives the menu buttons, the `hide=`/`show=` URL keys, the legend block and the blocker fade for it. The Galaxy and Sector stage view registers its sector kinds as layers; the System Map uses it for its own (course, facilities and the rest); the 3D galaxy map has the legend block (`galaxymap3d-layer-legend`). A kind absent from a sector is left out of the menu but keeps its URL entry.

### Added
- Tests: the registry in `maplayers.test.mjs` (16), the System Map layer menu, and a browser test that a registered Sector Map layer has a menu button, a URL key and a legend line.

TODO items resolved: MAP.217 / Added: none
