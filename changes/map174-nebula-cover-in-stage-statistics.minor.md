**Headline:** The Galaxy Map's stage and window data now say how much of each cell lies inside a nebula, generated or not.

### Added
- **Nebula cover in the stage-view cell statistics (MAP.174).** `/galaxy/stage` and `/galaxy/window` carry `nebula_cover`, the `[shell, band, wedge, cover]` of every cell whose share inside a nebula is above zero, for cells that hold a generated sector and cells that do not, and `stats.nebula_cover` on each generated child. The figures are MAP.172's stored block aggregates, made from the seed's nebula field and the stored nebulae. They are read on every request and not kept with the cached stages, which only change when a sector does. The Color by choice that draws them is MAP.177.
- **`GET /api/galaxy/nebula-cover?at=...[&g=...]`** returns that cover for a stage or a window.
