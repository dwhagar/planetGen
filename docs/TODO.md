# planetGen Roadmap

## How to use this document

This file tracks **open work only**. Finished work is not listed here:
`CHANGELOG.md` and git history are the record. For how things work today,
read the code and the reference docs (`src/stellarObjects/schema.sql` and
`docs/database-schema.md` for the schema, `docs/api.md` and
`docs/html-interface.md` for the web interface,
`docs/design/galaxy-coordinate-system.md` for galaxy geometry). When an
item ships, delete it here, renumber the rest, and describe the change
in the PR's `changes/` note.

Each open item says what's wrong (the observed symptom), where to start
looking, and what "done" means, so it can be picked up cold.

## Background

The long-term goal is a fully populated galaxy (every sector, system,
planet, moon and belt) stored in MySQL and browsed through a web
interface served from Apache. The original roadmap phases are all done:
object-graph serialization, the relational schema and migrations, CLI
tools writing to the database, lazy galaxy-scale generation from a
density skeleton, and the Flask API (`src/html/api/`) with server-rendered
pages served by the same app (`src/html/web/`): Galaxy, Sector and System
maps, search, NAV, admin auth and wiki publishing.

## Open items

Numbered in the order to work them, from the smallest change to the
largest, with new bug fixes put first. The numbers run across groups;
renumber when items are added or finished.

### Plan: what to do first

- **Fix the asteroid belt overlap (1)** first: it's a generator bug
   that turns CI red at random.
- **Extend the cache (2)**, then do the System Map route (3).
- **Do the Galaxy Map rework (4) last** of the active work: it's the
   largest change and builds on the cache.

### Generator bugs (`src/stellarObjects/systemData.py`)

1. [ ] **An asteroid belt can overlap the next planet out.** Seen
   intermittently in CI (Python 3.9, run 36744097318, while testing PR
   #105, which didn't touch the generator):
   `src/tests/test_systems.py::test_each_tristate_flag_forced[K1III-False-COMETS]`
   failed in `assert_no_orbital_overlap` with "objects 12 (a) and 13 (t)
   overlap: gap=-1.6548277320167166, required>=0.05". The belt's
   `upper_limit` reached 1.65 AU past the next planet's orbit, where
   `MIN_ASTEROID_BELT_SEPARATION` requires at least 0.05 AU of clearance.
   It was a random K1III (giant) system with comets forced off; later runs
   passed. Start with `StarSystem.validate_system`'s correction pass
   (around the `MIN_ASTEROID_BELT_SEPARATION` case): it either doesn't
   account for a belt's width when pushing the following body out, or
   stops before it converges for a giant star's wide zones. Reproduce by
   looping that test's generation with fixed seeds until it fails. Done
   means a failing seed is found, it passes after the fix, and a
   regression test is pinned to it.

### Performance

2. [ ] **Add a cache so pages don't hit the database on every request.**
   Every page calls the Flask API through `html/lib/apiclient.py`,
   and every API route queries MySQL fresh, including results that rarely
   change (`/api/galaxy/sectors`, `/api/galaxy/shape`, sector and system
   detail). The 3D Galaxy Map already has one: its cube tiles are cached
   on disk by the web layer (`html/lib/tilecache.py`) and in the browser,
   and an edit refreshes only the tiles it touched (`GET /api/galaxy/
   changes`, from the v27 `modified_at` columns). Extend that pattern to
   the other pages, invalidating each page from its own rows'
   `modified_at`.

### System Map (`lib/systemmap.py`, `static/systemmap.js`)

3. [ ] **Draw the Measure distance path and route it around obstacles.**
   "Measure distance" ([5.46.32]) reports a straight-line distance plus,
   when the line crosses the scene's central body, a tangent-and-arc
   detour around that one body (`computeMeasurement`,
   `routeAroundCircle`). It should also draw the path on the map, and the
   route should avoid every body it would pass through (planets, moons,
   either star of a binary), keep a safe distance from stars rather than
   just clearing the surface, and not thread between the two stars of a
   close binary. Done means the drawn path and reported distance agree
   and both respect those clearances.

### Galaxy Map (`src/html/web/galaxy_views.py`, `lib/galaxymap3d.py`, `static/galaxymap3d.js`, `queryDb.galaxy_tiles`)

4. [ ] **Rework the Galaxy Map; it isn't useful in its current form.**
   Investigate a representation driven by the real galaxy/sector geometry
   (`stellarObjects/galaxyGeometry.py`: the cylindrical ring/layer/slot
   grid) instead of per-sector sprites plus an illustrative density
   cloud. Idea to evaluate first: group neighboring sectors and draw them
   as arc segments of rings (a band of rings and layers by an angular
   range), shaded by density, so the map shows the galaxy's
   structure at every zoom without enumerating individual sectors. Done
   means a written comparison of options, then an implementation that
   loads quickly at full-galaxy zoom. (Today's density cloud also renders
   black: its `InstancedMesh` material sets `vertexColors: true` with no
   per-vertex `color` attribute, which zeroes the per-instance colors.
   Turning that off alone saturates it to white, so it needs retuning or
   replacing as part of this.)

### Web API (`src/html/api/routes.py`)

Low priority; nobody is waiting on these.

5. [ ] **The API can't create a system inside an existing sector.**
    `POST /api/systems` only creates standalone systems (`sector_id =
    NULL`, see `docs/api.md`). Attaching one to a sector needs the sector's
    placement and Hill-sphere separation logic (`SpaceSector.add_system`),
    which was left out of the write API to keep the admin-auth change
    small.

6. [ ] **The API can't edit a system's generated content.** `PATCH
    /api/systems/<id>` only renames. Changing stars/planets/moons/belts
    means `DELETE` then `POST` (regenerate). It may never need solving;
    kept here in case it does.

## Population and Politics

Exploratory ideas, not yet designed. Each needs a design pass before it
can be ordered against the work above.

7. [ ] Assign government ownership to star systems so that groups of
    systems form territories mapped in 3D space.
8. [ ] Flag worlds with life for generated names of their dominant
    species.
9. [ ] A database of spacefaring species.
10. [ ] Model younger and older civilizations: what differs with a
    society's age and how to store and present it.
