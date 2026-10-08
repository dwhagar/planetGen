### Changed
- A nebula's page now shows the nebula itself: its irregular shape in 3D, with the galaxy's brightest stars dimmed around it for reference (drag to turn, scroll to zoom), instead of the flat AU diagram (MAP.105). New `GET /api/nebulae/<id>/surroundings` supplies the stars.

### Fixed
- The "-" button on a supernova remnant's AU diagram works from the first view: a remnant too big to open inside the 1 ly limit may zoom out to twice its opening view (UX.38).
