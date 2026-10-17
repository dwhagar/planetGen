**Headline:** Every scatter and fill limit on the Generate page is a slider with a number box at its end, in segments for the core, the halo and the crowded cells.

### Changed
- **The luminosity floor is a slider, not a dropdown (UX.116).** The galaxy-wide mass limit and luminosity floor, on Plan, New galaxy and Redo scatters, are sliders that move a number box at their end; the box shows the value and takes a typed one. The luminosity floor keeps its presets as the slider's stops and now takes any number between them.

### Added
- **Three segments under the galaxy-wide limits on Plan and New galaxy.** *Core*: its luminosity floor and mass limit, which follow the galaxy-wide ones (110% of them) until one is set, with a warning if the floor goes below the galaxy-wide one. *Halo and intergalactic space*: its luminosity floor (default 500, from 50). *Dense cells*: the cell star ceiling (default 5,000, with what a full cell costs in time and storage) and the fill luminosity cap (default 100). Each segment has a Reset to the defaults. They set `planetgen plan --core-min-luminosity`, `--core-mass-limit`, `--halo-min-luminosity`, `--cell-star-ceiling` and `--fill-luminosity-cap`, and a number outside a slider's range is refused with the field named.

TODO items resolved: UX.116. Added: none.
