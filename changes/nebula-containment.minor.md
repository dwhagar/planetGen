### Added

- **Schema v39: what sits inside a nebula.** Star systems, rogue planets,
  interstellar comets, black holes, neutron stars, asteroid fields and
  nebulae now record the innermost nebula or supernova remnant they sit
  inside (`inside_nebula_id` / `inside_remnant_id`). This is set by a 3D
  distance test when a sector is generated and whenever a nebula or
  remnant is placed, so sectors generated later inside an existing cloud
  see it. The migration fills it for existing data, and the sector detail
  query reports each system's cloud.
