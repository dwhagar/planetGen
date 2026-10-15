### Changed
- The API, the web pages and the Generate page speak the spherical sector grid: sectors are addressed by shell, band and wedge (`/api/galaxy/cell`, `/api/galaxy/bright-stars` and `/api/galaxy/uncharted` take `?shell=&band=&wedge=`, `/api/galaxy/shape` returns `shells`), the Generate page's modes are random start, a whole shell, around a sector, one address and a span of shells, bands or wedges, and the sector page and nav chart show the new address.
- The Galaxy Map and the Sector Map show a notice that the map is being rebuilt for the new grid until the map work lands.

### Removed
- The Generate page's whole-ring, column and Galaxy Map block modes, and the radius-in-sectors cylinder fields (the ball radius replaces them).

TODO items resolved: UX.111 / Added: none
