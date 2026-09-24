### Fixed
- **Galaxy Map: a generated sector could show up as a huge sphere
  filling the view.** Each placed sector's dot is kept a constant size
  on screen, but its two sprites were being sized from the camera's
  distance to the galactic center instead of to the sector, so a sector
  far from the core grew into a giant disc once zoomed in on it. The
  size now comes from the sector's own position, and the dots are
  smaller (3-7 px instead of 4-16 px).
- **System Map: orbit lines were missing.** The orbits layer's
  `z-index: -1` put it behind the map's own background, because the map
  box didn't form its own stacking context. The viewport now isolates
  its stacking, so orbits draw again, still hidden behind each body's
  sphere.
- **Star and atmosphere glows rendered as flat, opaque discs** on the
  System Map and Sector Map. The glow shader measured the rim the
  front-face way on a back-face-only shell, which clamps to full
  brightness everywhere. It now fades from full brightness at the body's
  limb to nothing at the shell's edge (`bodyRendering.makeGlowMaterial`
  takes the shell's scale for this), and the star and planet glow
  strengths were retuned for the new falloff.
- **The "nearest" systems on a system page weren't links.** They were
  parsed out of the stored `location` text, whose names are written
  before each neighbor's name is made unique and never follow a rename,
  so they rarely matched a real system. `GET /api/systems/<id>` now
  returns `nearest_neighbors` (id, current name, distance) computed from
  positions, and the page links every one of them.

### Added
- **`tests/test_skeleton_shape.py`** confirms the unfilled-sector
  skeleton's slot centers form a thin disk in galaxy-frame parsecs; the
  sphere the Galaxy Map draws comes from its 200 pc planned-tier radius
  cap, which the tile-based fetching work replaces.
