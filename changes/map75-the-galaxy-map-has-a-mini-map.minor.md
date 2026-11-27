**Headline:** The Galaxy Map has a mini map showing the whole galaxy and where you are looking from

### Added
- The Galaxy Map has a mini map MAP.75: a small square in the bottom right corner showing the whole galaxy from the direction you are looking, with a cross where the view is centred and a line to your eye. It is a second, locked camera on the map's own renderer, so it never zooms and costs no second WebGL context. It is left out on a canvas too small for it and on the Sector Map.

generation-output: unchanged
