### Added
- **The drill-down's stage API.** `GET /api/galaxy/stage?at=m.ring.wedge.slab`
  (and the site's cached `/galaxy/stage`) returns how many generated
  sectors each block inside a drill-down block holds, and the sectors
  themselves at the smallest level. It feeds the Galaxy Map's coming
  drill-down navigation; nothing on the map changes yet.
