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
density skeleton, and the Flask API (`src/html/api/`) behind a thin CGI
frontend (`src/html/`) with Galaxy, Sector and System maps, search, NAV,
admin auth and wiki publishing.

## Open items

Numbered in the order to work them: bug fixes first (the production
crash before anything else), then everything else from the smallest change
to the largest. The numbers run across groups; renumber when items are
added or finished.

### Plan: what to do first

- **Stop the crash (1).** Production stability comes before everything
   else, and its root cause shapes both the cache (3) and the Galaxy Map
   rework (6). The skeleton drawing fix (2) rides along with it.
- **Design the cache (3) once the crash cause is known**, then do the
   System Map route (4) and the phenomena schema work (5).
- **Do the Galaxy Map rework (6) last** of the active work: it's the
   largest change and builds on 1, 2 and 3.

### Bug fixes

Things that are broken today, most urgent first, then smallest to largest.

1. [ ] **Apache was OOM-killed on the production server (2026-09-24
   01:13 UTC).** `systemctl status apache2` showed `Result: oom-kill`
   after ~56 CPU-minutes. Just before it, `planetgen.error.log` shows
   `/galaxy_view.py` CGI requests timing out on their call to
   `/api/galaxy/view` (`apiclient._request`, 30 s `_TIMEOUT_SECONDS`) and
   then `Truncated or oversized response headers received from daemon
   process 'planetgen-api'` for every in-flight request, i.e. the WSGI
   daemon died. The access log shows the Galaxy Map's opening view asking
   for `/api/galaxy/view?cx=0&cy=0&cz=0&radius_pc=15082.8` (a 1.7 MB
   response) and `/api/galaxy/sectors` (508 KB) on each page load, and
   bursts of large-radius `/api/galaxy/view` calls returning 500/504.
   Suspect the Galaxy Map's view queries (`queryDb.galaxy_view`,
   `stellarObjects/galaxyViewport.py`) holding large result sets per
   request across several daemon threads. Being investigated in its own
   project thread; done means the cause is known and the map can't take
   the server down.

2. [ ] **The unfilled-sector skeleton draws as a sphere instead of a
   spiral galaxy.** The data is right: `tests/test_skeleton_shape.py`
   samples real slot addresses and confirms the qualifying ones sit on
   their shell radii in galaxy-frame parsecs and form a thin disk (out to
   ~14 kpc, 95% within ~0.9 kpc of the plane, thinning outward). The ball
   comes from the drawing: the "planned" tier
   (`galaxyViewport.planned_slots_in_view`) is every qualifying slot
   within `PLANNED_RADIUS_CAP_PC` (200 pc) of the view center, so it is
   always a sphere around wherever the camera is looking, and a solid one
   near the core, where every slot qualifies. Being fixed together with
   item 1, whose cube-tile fetching replaces that query; done means the
   planned tier shows qualifying slots across the whole view (or a
   per-tile sample of them at wide zoom) instead of one ball.

### Performance

3. [ ] **Add a cache so pages don't hit the database on every request.**
   Every CGI page calls the Flask API through `html/lib/apiclient.py`,
   and every API route queries MySQL fresh, including results that rarely
   change (`/api/galaxy/sectors`, `/api/galaxy/shape`, sector and system
   detail, the Galaxy Map's view payloads). Design a cache layer (where it
   lives: API process memory, a shared store, or HTTP caching headers;
   what the keys are; and how writes, generation and the orbit updater
   invalidate it). Wait for item 1's root cause first: repeated large
   galaxy queries are part of that picture.

### System Map (`lib/systemmap.py`, `static/systemmap.js`)

4. [ ] **Draw the Measure distance path and route it around obstacles.**
   "Measure distance" ([5.46.32]) reports a straight-line distance plus,
   when the line crosses the scene's central body, a tangent-and-arc
   detour around that one body (`computeMeasurement`,
   `routeAroundCircle`). It should also draw the path on the map, and the
   route should avoid every body it would pass through (planets, moons,
   either star of a binary), keep a safe distance from stars rather than
   just clearing the surface, and not thread between the two stars of a
   close binary. Done means the drawn path and reported distance agree
   and both respect those clearances.

### Sector Map and phenomena (`lib/starmap.py`, `static/sectormap.js`, `queryDb.py`, schema)

5. [ ] **Show every kind of stellar phenomenon on the Sector Map, all
   clickable.** Only nebulae, asteroid fields, black holes and neutron
   stars appear (`queryDb._PHENOMENON_TABLES`). Supernova remnants, rogue
   planets and interstellar comets have no galaxy-frame position columns
   (`_UNPLACED_PHENOMENON_TABLES`, see [5.46.26]), so they're absent from
   the Sector Map, Galaxy Map and NAV. Done means all seven types are
   placed (schema migration plus generation in `phenomenonGen.py`/
   `SpaceSector.add_phenomenon`), returned by `phenomena_near_sector`,
   drawn by `sectormap.js` and link to `phenomenon.py`.

### Galaxy Map (`src/html/galaxy.py`, `lib/galaxymap3d.py`, `static/galaxymap3d.js`, `queryDb.galaxy_view`)

6. [ ] **Rework the Galaxy Map; it isn't useful in its current form.**
   Investigate a representation driven by the real galaxy/sector geometry
   (`stellarObjects/galaxyGeometry.py`: concentric Fibonacci-sphere shells
   of Voronoi sector slots) instead of per-sector sprites plus an
   illustrative density cloud. Idea to evaluate first: group neighboring
   sectors and draw them as arc segments of shells (a band of shells by an
   angular range), shaded by density, so the map shows the galaxy's
   structure at every zoom without enumerating individual sectors. Done
   means a written comparison of options, then an implementation that
   loads quickly at full-galaxy zoom. (Today's density cloud also renders
   black: its `InstancedMesh` material sets `vertexColors: true` with no
   per-vertex `color` attribute, which zeroes the per-instance colors.
   Turning that off alone saturates it to white, so it needs retuning or
   replacing as part of this.)

### Web API (`src/html/api/routes.py`)

Low priority; nobody is waiting on these.

7. [ ] **The API can't create a system inside an existing sector.**
    `POST /api/systems` only creates standalone systems (`sector_id =
    NULL`, see `docs/api.md`). Attaching one to a sector needs the sector's
    placement and Hill-sphere separation logic (`SpaceSector.add_system`),
    which was left out of the write API to keep the admin-auth change
    small.

8. [ ] **The API can't edit a system's generated content.** `PATCH
    /api/systems/<id>` only renames. Changing stars/planets/moons/belts
    means `DELETE` then `POST` (regenerate). It may never need solving;
    kept here in case it does.

## Population and Politics

Exploratory ideas, not yet designed. Each needs a design pass before it
can be ordered against the work above.

9. [ ] Assign government ownership to star systems so that groups of
    systems form territories mapped in 3D space.
10. [ ] Flag worlds with life for generated names of their dominant
    species.
11. [ ] A database of spacefaring species.
12. [ ] Model younger and older civilizations: what differs with a
    society's age and how to store and present it.
