### Fixed

- A sector save no longer fails with "existing_count must be >= 0, got
  -1" (TEST.85). When a system name with no decoration left was drawn
  again during name reservation, the new name was also counted as a
  holder of the registry row it matched later in the same pass, leaving
  that row one holder short. Each name's key is now taken once per pass,
  so a redrawn name only counts in the next pass. Names that never
  collide come out exactly as before.
