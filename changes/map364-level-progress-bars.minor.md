**Headline:** The Galaxy Map shows a progress bar for each level of map data while it loads, so a slow connection shows how far along it is.

### Changed
- **MAP.364.** While the nodes a view asks for come in from the server, each level of the block ladder with nodes still to come shows a bar in the top left of the map: nodes received over nodes requested for that level (`static/levelprogress.js`), filled as the response is read (`tilepack.readNodeBody` reports each node as its bytes arrive). A load that ends within 400 ms, which is every level served from this browser's store or the server's cache, shows no bar, and the bars go when every level is in. A failed load, or a node that comes back incomplete, replaces its level's bar with a line in plain words that says it could not load and when it will try again. Each bar is a `role="progressbar"` named for its level and takes the map viewport's colors in both themes. The tests hold the bars against a throttled connection, a revisit that finds its nodes stored, and a refused request. The map has twelve block levels, not the six steps (Z0 to Z5) of the zoom plan, so the bars are by block level.

generation-output: unchanged

TODO items resolved: MAP.364.
