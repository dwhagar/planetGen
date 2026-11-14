**Headline:** A sector page now keeps one wireframe fixed round its own sector with an icon on the middle of each edge that opens the sector beside it, and the stars of one wide binary or multiple system are joined by dashed links and named once on the Galaxy Map.

### Changed

- **The sector page's wireframe stays on its sector (MAP.427).** The sector page (the Galaxy Map locked to one sector) draws the sector's own locked wireframe, as an opened sector in the Galaxy Map does (MAP.327), and no longer moves it with the camera. On the middle of each edge that has a charted sector beside it there is a point-of-interest icon named for that sector ("Sector ..."); a click opens that sector's page as the new focus. No icon is drawn on a galaxy edge or where the neighbouring sector has no page. `galaxymap3d.js` builds the frame for a pinned map too (`data.pinned`), and `galaxystageview.js` gives it `neighbourPageId` and `openNeighbourPage`.
- **The stars of one system are shown as one system (MAP.428).** `static/systemgroups.js` (pure, `tests/js/systemgroups.test.mjs`) finds the systems with two or more stars among the stars the Galaxy Map draws. The members of such a system are joined by thin dashed links (`.map-system-links`, an SVG over the map) from the brightest to each of the others once they stand at least 4 px apart on the screen, and the system is named once, by its name without the component letter ("Ofriggor"), at the middle of its members, instead of each member being named. The links and the shared name follow the "Star names" layer.
- Tests: a browser test holds the sector page's frame fixed under drag, key and wheel and follows an edge icon to the neighbouring sector; a second one checks the link between the stars of the dense fixture's close pair (now one system) and its absence once the pair is one dot.

generation-output: unchanged
generation-affects: none
