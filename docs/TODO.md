# planetGen Roadmap

## How to use this document

This file tracks **open work only**. Finished work is not listed here:
`CHANGELOG.md` and git history are the record. For how things work today,
read the code and the reference docs (`src/stellarObjects/schema.sql` and
`docs/database-schema.md` for the schema, `docs/api.md` and
`docs/html-interface.md` for the web interface,
`docs/design/galaxy-coordinate-system.md` for galaxy geometry). When an
item ships, delete it here and describe it in the PR's `changes/` note.

Each open item says what's wrong (the observed symptom), where to start
looking, and what "done" means, so it can be picked up cold.

## Background

The long-term goal is a fully populated galaxy (every sector, system,
planet, moon and belt) stored in MySQL and browsed through a web
interface served from Apache. The original roadmap phases are all done:
object-graph serialization, the relational schema and migrations, CLI
tools writing to the database, lazy galaxy-scale generation from a
density skeleton, and the Flask API (`src/html/api/`) behind a thin CGI
frontend (`src/html/`) with Galaxy, Sector and System maps, search, NAV,
admin auth and wiki publishing.

## Open items

### Galaxy Map (`src/html/galaxy.py`, `lib/galaxymap3d.py`, `static/galaxymap3d.js`, `queryDb.galaxy_tiles`)

- [ ] **Rework the Galaxy Map; it isn't useful in its current form.**
  Investigate a representation driven by the real galaxy/sector geometry
  (`stellarObjects/galaxyGeometry.py`: concentric Fibonacci-sphere shells
  of Voronoi sector slots) instead of per-sector sprites plus an
  illustrative density cloud. Idea to evaluate first: group neighboring
  sectors and draw them as arc segments of shells (a band of shells by an
  angular range), shaded by density, so the map shows the galaxy's
  structure at every zoom without enumerating individual sectors. Done
  means a written comparison of options, then an implementation that
  loads quickly at full-galaxy zoom.
- [ ] **The unfilled-sector skeleton forms a sphere at the center
  instead of a spiral galaxy.** The "planned" tier (real,
  not-yet-generated sector addresses,
  `galaxyViewport.planned_slots_in_view`) draws as a ball around the
  galactic center, while the illustrative density cloud
  (`density_sample_points`) does look like a spiral. Verify that the
  centers of the qualifying unfilled sectors trace the disk/arms/bulge
  the density skeleton (`galaxy_shape`, `galaxyDensity.predicted_star_count`)
  describes. Likely cause: the old view query listed the 4,000 slots
  nearest the view center, which at the full-galaxy starting view is a
  ball around the core. Since the map moved to cube tiles, planned slots
  are only fetched within 20 pc of the camera target once zoomed in
  (`lib/galaxymap3d.py`'s `PLANNED_*` settings), and the ball no longer
  shows at full-galaxy zoom. Still to do: confirm with real data that
  qualification against the shape is applied and the slot centers are
  galaxy-frame parsecs. Done means a test that samples slot centers and
  asserts they follow the disk, plus any fix.
- [ ] **Remove the large sphere marker drawn for a star in a sector.**
  The Galaxy Map still shows a large sphere for a star inside a placed
  sector, left over from the idea of showing the brightest stars at
  galaxy scale, which proved unhelpful. Remove the marker and whatever
  API/query fields exist only to feed it.

### Performance

- [ ] **Extend caching beyond the Galaxy Map's tiles.** The 3D map's
  cube tiles are now cached on disk by the web layer
  (`html/lib/tilecache.py`, keyed by `GET /api/galaxy/stamp`) and in the
  browser. Every other CGI page still calls the API fresh, including
  results that rarely change (`/api/galaxy/sectors`, which the Galaxy
  Map page loads in full for its tables on every visit, `/api/galaxy/shape`,
  sector and system detail). Extend the same pattern (a disk cache in
  the web layer keyed by a cheap content stamp) where it pays off, and
  decide how the orbit updater's writes invalidate system detail.

### Sector Map (`src/html/sector.py`, `lib/starmap.py`, `static/sectormap.js`)

- [ ] **Show every kind of stellar phenomenon on the Sector Map, all
  clickable.** Only nebulae, asteroid fields, black holes and neutron
  stars appear (`queryDb._PHENOMENON_TABLES`). Supernova remnants, rogue
  planets and interstellar comets have no galaxy-frame position columns
  (`_UNPLACED_PHENOMENON_TABLES`, see [5.46.26]), so they're absent from
  the Sector Map, Galaxy Map and NAV. Done means all seven types are
  placed (schema migration plus generation in `phenomenonGen.py`/
  `SpaceSector.add_phenomenon`), returned by `phenomena_near_sector`,
  drawn by `sectormap.js` and link to `phenomenon.py`.

### System Map and shared rendering (`lib/systemmap.py`, `static/systemmap.js`, `static/bodyRendering.js`)

- [ ] **Star glow renders as an opaque shell, not a 3D glow.** The
  fresnel glow shell (`bodyRendering.makeGlowMaterial`, used by
  `systemmap.js` and `sectormap.js`, strengthened in [5.46.24]) should
  read as light fading out from the star's limb, but it shows as a solid
  opaque layer. Investigate the material's blending, transparency,
  `depthWrite` and alpha falloff. Verify with screenshots of both maps
  before and after (Playwright with the bundled Chromium works here).
- [ ] **Bring orbital paths back without drawing them over bodies.**
  Orbit lines disappeared after [5.46.24] split each scene into an
  orbits-only `<svg>` under the sphere `<canvas>` and a marker `<svg>`
  over it. They should be visible again while still being hidden behind
  the star/planet spheres. Check that `showScene` in `systemmap.js`
  shows the orbits layer and that the canvas isn't painting an opaque
  background over it.
- [ ] **Draw the Measure distance path and route it around obstacles.**
  "Measure distance" ([5.46.32]) reports a straight-line distance plus,
  when the line crosses the scene's central body, a tangent-and-arc
  detour around that one body (`computeMeasurement`,
  `routeAroundCircle`). It should also draw the path on the map, and the
  route should avoid every body it would pass through (planets, moons,
  either star of a binary), keep a safe distance from stars rather than
  just clearing the surface, and not thread between the two stars of a
  close binary. Done means the drawn path and reported distance agree
  and both respect those clearances.

### Web API

- [ ] **The API can't create a system inside an existing sector.**
  `POST /api/systems` only creates standalone systems (`sector_id =
  NULL`, see `docs/api.md`). Attaching one to a sector needs the sector's
  placement and Hill-sphere separation logic (`SpaceSector.add_system`),
  which was left out of the write API to keep the admin-auth change
  small.
- [ ] **The API can't edit a system's generated content.** `PATCH
  /api/systems/<id>` only renames. Changing stars/planets/moons/belts
  means `DELETE` then `POST` (regenerate). It may never need solving;
  kept here in case it does.

## Population and Politics

Exploratory ideas, not yet designed:

- [ ] Assign government ownership to star systems so that groups of
  systems form territories mapped in 3D space.
- [ ] Flag worlds with life for generated names of their dominant
  species.
- [ ] A database of spacefaring species.
- [ ] Model younger and older civilizations: what differs with a
  society's age and how to store and present it.
