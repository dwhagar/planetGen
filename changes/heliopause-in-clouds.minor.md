### Added

- **A nebula squeezes the heliosphere of every system inside it.** The
  gas of a nebula or supernova remnant pushes on a star's wind bubble
  much harder than open space does, so a system inside a dense cloud
  now reports a much smaller heliosphere (a Sun-like star's drops from
  about 85 AU to well under 1 AU in a dense cold cloud). The system
  page text says so, and the system API returns `inside`,
  `heliopause_au` and `heliopause_open_space_au`; navigation uses the
  squeezed heliopause as the edge of a system's local frame.
