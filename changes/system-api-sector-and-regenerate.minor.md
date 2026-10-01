### Added

- **The API can add a system to an existing sector and regenerate a
  system in place.** `POST /api/systems` takes an optional `sector_id`
  (and `position`): the new system is placed clear of the sector's other
  systems' Hill spheres, with its location, containment and nearest
  systems filled in. `PATCH /api/systems/<id>` takes `{"regenerate":
  recipe}` to replace a system's stars, planets, moons, belts and comets
  while keeping its id, name, place and links.
