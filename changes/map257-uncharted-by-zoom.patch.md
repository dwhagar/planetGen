**Headline:** Zooming in on a cell nothing was generated in, when the scatters left objects in it, now opens it as Uncharted, as it opens a charted sector.

### Fixed

- The wheel opened only sectors that had been generated: the sector the camera came up to was chosen among the charted ones near the eye, so a cell holding a bright star, a black hole or a nebula but no generated sector never opened by zooming to it, and the charted neighbour opened instead (MAP.257). At the end of the zoom range the cell under the camera's middle now opens as "Uncharted" by the same path a click uses (MAP.162: the scene of `GET /galaxy/uncharted/<shell>/<band>/<wedge>/scene`, `enterSector` with no record); a cell the server says is empty stays unopened and the nearest charted sector opens as before.

generation-output: unchanged
