# API compatibility

<!-- Built from docs/reference/api-versions.json by scripts/build_api_compatibility.py. Do not edit by hand. -->

The API has a version of its own, one integer that is not the release number (`src/planetgen/api/version.py`, API.22). Every `/api` answer carries it in the `X-PlanetGen-API-Version` header, `GET /api/health` reports it as `api_version`, and the `User-Agent` of everything the software sends names it (`PlanetGen-API/<n>`). It goes up by one, by hand, in the pull request that makes a breaking change to an endpoint: a removed or renamed route or field, a changed meaning or type, a new required parameter. Additive changes (a new route, a new optional field or parameter) leave it alone. That pull request also adds the row below to `docs/reference/api-versions.json`, a line to the change list in [api.md](api.md#api-version-and-object-ids), the route snapshot `src/tests/fixtures/api_routes.json` and `api-version: changed` to its release note; the stamp fills in the release. `src/tests/test_api_compatibility.py` fails when any of them is missing.

## Versions

| API version | First release | Items | Change | Routes and payload formats changed |
|---|---|---|---|---|
| **12** | 8.0.2993 | GEN.1108 | A planet or moon gives its climate in seven bands of latitude, not in a few named places. | GET /api/systems/<id>: the regions of a planet or moon's equipment is {bands, mean_swing_k}: bands is the seven bands from the north pole to the south pole (the day centre to the night centre for a world that turns once an orbit), each {name, mean_k, swing_k, day_hours, night_hours, longest_day_days, longest_night_days, state}; polar, north, equator, south, day_side, terminator, night_side and the swings by those names are gone<br>GET /api/objects/<ref>/sheet: a planet or moon has bands, the table of its seven bands, and no temp_polar_k, temp_north_k, temp_equator_k, temp_south_k, temp_swing_k, temp_day_side_k, temp_terminator_k or temp_night_side_k rows |
| **11** | 8.0.2617 | API.33 | The sizes of a backup archive and of its parts say they are in bytes. | GET /api/admin/backups and GET /api/admin/backups/<name>: an archive and each of its parts have size_bytes where they had size |
| **10** | 8.0.2608 | API.37 | A planet or moon names its life by class, not by chemistry. | GET /api/systems/<id>: planets and moons have no life_chemical; they have life_class (the class number), life_class_variant, life_class_best, life_class_best_score, life_animal_grade and life_complexity (the columns), and life_class_id (such as 1a) and life_class_scores (the score of every class the body supports)<br>the species of GET /api/planets/<id>/species and GET /api/planets/<id>/habitable-species have life_class (the class id) in place of life_chemical |
| **9** | 8.0.2271 | MAP.287 | A star, planet or moon of a system scene carries a look record, not a colour. | GET /api/systems/<id>/scene: stars, planets and moons have look (the body look record, which holds the color and what the painters need) and no color; a comet keeps color |
| **8** | 8.0.2236 | MAP.262 | The Galaxy Map reads its data by sector-block node; the cube tiles are gone. | GET /api/galaxy/tiles removed; GET /api/galaxy/nodes is the map's data<br>GET /api/galaxy/changes answers nodes, no longer tiles |
| **7** | 8.0.2225 | API.26, API.31 | Every sector address is shell, band and wedge, where it was ring, layer and slot. | request parameters shell, band and wedge; answers shell_index, band_index and wedge_index (were ring, layer, slot; ring_index, layer_index, ring_slot_index)<br>no route that takes the old names is kept |
| **6** | 8.0.2215 | PERF.99 | The work queue endpoints answer runs held in Redis, not job trees. | GET /api/admin/work and GET /api/admin/work/<id> return runs with their stage, units and time left; <id> is a run id<br>POST /api/admin/work/<id>/control takes cancel or clear (pause and resume of a node are gone)<br>POST /api/admin/work/<id>/delete and clear-finished delete runs<br>GET /api/admin/work/<id>/log is new<br>POST /api/admin/work/lease/clear removed |
| **5** | 8.0.1783 | DB.51 | A system's description is a document built from stored data, not a table of paragraphs. | GET /api/systems/<uid>/sections removed |
| **4** | 8.0.1589 | ADM.50 | The galaxy settings files are gone. | GET /api/admin/galaxy-settings removed<br>GET /api/admin/galaxy-settings/<name> removed |
| **3** | 8.0.912 | API.23 | Stars, planets, moons, asteroid belts, comets and facilities are named by their object ID too. | fields id, star_id, planet_id, moon_id, belt_id, comet_id, host_id and homeworld_planet_id; <kind>:<ID> references<br>URLs /api/planets/<id>, /api/moons/<id>/class, /api/facilities/<id><br>POST /api/facilities and GET /api/facilities/orbit take the host's ID as host_id |
| **2** | 8.0.912 | API.23 | Sectors, systems and phenomena are named by their object ID, not by their database row number. | id, sector_id and ref fields of sector, system and phenomenon answers carry the printed object ID (0008000000-4000000-000)<br>URLs /api/systems/<id>, /api/sectors/<id>, /api/phenomena/<type>/<id>, /api/nebulae/<id>/...<br>from, to and course references such as system:0008000000-4000000-000<br>a row number in any of these is refused (404) |
| **1** | 8.0.866 | API.22 | The first numbered version: every /api answer carries it in the X-PlanetGen-API-Version header and /api/health reports it. | none (the first numbered version) |

## Clients and the servers they talk to

| Client | Where | Speaks | Talks to |
|---|---|---|---|
| The web pages | `static/apiversion.js` (adds the version header to each call a page makes), `web/lib/apiclient.py` | the API version of the release they ship in | the server they were served by: pages and server are one release, so they always agree |
| Remote generate (`generate.py --remote`, API.3) | not built yet | the API version of the release it is run from | a server of the same API version (API.5 widens this to the current major and, for one release cycle after a major bump, the previous one, for upload payloads only) |
| The command line (`planetgen ...`) | `src/planetgen/cli/` | none: it opens the database directly | no server involved |

## Which releases speak which version

Every client that exists today ships with the server, so the release a page or a command comes from decides the version it speaks:

| API version | Releases that speak it |
|---|---|
| **12** | 8.0.2993 and later |
| **11** | 8.0.2617 up to the release before 8.0.2993 |
| **10** | 8.0.2608 up to the release before 8.0.2617 |
| **9** | 8.0.2271 up to the release before 8.0.2608 |
| **8** | 8.0.2236 up to the release before 8.0.2271 |
| **7** | 8.0.2225 up to the release before 8.0.2236 |
| **6** | 8.0.2215 up to the release before 8.0.2225 |
| **5** | 8.0.1783 up to the release before 8.0.2215 |
| **4** | 8.0.1589 up to the release before 8.0.1783 |
| **3** | 8.0.912 up to the release before 8.0.1589 |
| **2** | none alone (replaced by version 3 in 8.0.912) |
| **1** | 8.0.866 up to the release before 8.0.912 |

A client of one API version used against a server of another is not supported until API.5 adds the check and the converters for older upload payloads (`GET /api/version`, the `PlanetGen-API-Version` header both ways, a refusal `400 api_version_unsupported`).
