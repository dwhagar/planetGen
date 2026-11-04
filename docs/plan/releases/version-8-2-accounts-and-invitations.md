# Version 8.2: User accounts with email invitations

Rebuilt on 2026-10-16 as a version map (Boss, 1:07 PM EDT: one version per headline, phases inside a version, blockers one phase earlier).
[docs/TODO.md](../../TODO.md) is the master file: it holds each item's full text, and its "Plan: versions" section indexes every version.
This file gives the version's goal, its phases and each phase's build threads with every item's prerequisites in order. A phase holds independent threads: an item's blockers sit in the same thread (earlier in the list) or in an earlier phase. Where the two disagree, TODO.md wins; when an item ships, it leaves TODO.md and its row here is deleted in the same PR (the progress ledger keeps its version). Research notes and judgment calls are in [notes.md](../notes.md); the headlines and the ideas not yet filed are in [release-roadmap.md](../release-roadmap.md).

Version: 8.2

Status: paused until 8.1 is complete (Boss, 7:44 PM EDT 2026-10-22: the five cores work 8.1 only; nothing in this version is started or owned until then, except the gates already running for UX.140, UX.128 and NAV.17).

## Goal

Two headlines: the continuous zoom (the Galaxy Map rebuilt so one camera drives every scene, MAP.334 to MAP.338, Boss 6:28 PM EDT 2026-10-28) and user accounts. User accounts, with sign-up by an email invitation that only an admin can send (USR.1): roles, SMTP settings, invitations, password emails and owner transfer first, then what a signed-in user gets: bookmarks in the account, API keys, the one-off system page and a home field.

## Phase 1: Accounts, invitations and bookmarks

Roles, the Owner tier of the settings page, SMTP, invitations, password emails, owner transfer and the bookmark manager.

### Accounts and sign-up

| ID | Item | Needs | Note |
|---|---|---|---|
| USR.26 | Roles in sessions, API keys and route guards; the Owner rules; the accounts CLI |  | Split from USR.2. |
| USR.3 | SMTP settings in the admin config | USR.26 |  |
| USR.28 | Account tokens: one issue-and-redeem library for invite, first-password, reset, email-change and owner-transfer links |  | Foundation found in the review. |
| USR.5 | Email loop for setting and resetting passwords | USR.3, USR.28 |  |
| USR.4 | Invite-only sign-up by unique link | USR.26, USR.3, USR.28, USR.5 |  |
| USR.6 | Owner transfer | USR.5, USR.28 |  |

### Bookmarks

| ID | Item | Needs | Note |
|---|---|---|---|
| UX.46 | Wiping the galaxy wipes the bookmarks |  |  |
| UX.179 | Sitemap entries and an Open Graph image for object pages (optional, can wait) |  |  |
| UX.206 | Link-preview (SEO and social) images for pages with a rendered image or map view are a cached snapshot of that view | UX.179 |  |
| VIEW.40 | System sky S22 (optional): cross-fade between the Galaxy Map's stars and the sky at the system, after the continuous zoom | VIEW.35 |  |
| VIEW.39 | System sky S21 (optional): aberration of the observer's motion, default off | VIEW.33 |  |
| VIEW.38 | System sky S20: design note, help text, citations and GEN.725 rows for the system sky | VIEW.35 |  |
| VIEW.37 | System sky S19: picking and hover text for the sky ('seen as it was N years ago'), with links for charted stars | VIEW.35 |  |
| VIEW.36 | System sky S18: sprites for nebulae and remnants, and luminous black holes, neutron stars and quasars as points | VIEW.28, VIEW.33 |  |
| VIEW.35 | System sky S17: wire the sky into systempage3d.js and createSystemView: on and off in preferences, a filling-in note and poll, a fallback without WebGL | VIEW.31, VIEW.33, VIEW.34 |  |
| VIEW.34 | System sky S16: the band layer: a sphere shader sampling the glow image with dither and grain | VIEW.31, VIEW.32 |  |
| VIEW.33 | System sky S15: skysphere.js: one points object centred on the camera, the vertex shader with the light-time passes driven by the orbit clock | VIEW.29, VIEW.32 |  |
| VIEW.32 | System sky S14: skylight.js: magnitude to opacity, size and colour, the GLSL twins, with a node test |  |  |
| VIEW.31 | System sky S13: GET /system/<uid>/sky: builds and caches the packed list and reports the tier state and counts, and takes ?limit= | VIEW.27, VIEW.29 |  |
| VIEW.30 | System sky S12: calibration test against the HYG catalogue on a developer machine: a generated sky has about 9,000 stars at 6.5 | VIEW.23, VIEW.25 |  |
| VIEW.29 | System sky S11: the packed sky list format: a Python writer, a JavaScript reader and one test that holds both | VIEW.23 |  |
| VIEW.28 | System sky S10: sprite data for nebulae and remnants: angular size, surface brightness, extinction | VIEW.21 |  |
| VIEW.27 | System sky S9: the glow band of the unresolved stars, with galactic dust only, ray-marched and cached per 50 pc cell | VIEW.20 |  |
| VIEW.26 | System sky S7b: a one-off command for sectors charted earlier: build the sky tiers on demand (Generate console and CLI) | VIEW.25 |  |
| VIEW.25 | System sky S7 (lock): charting builds the sky tiers: when a system is charted, the cells its sky needs are built with the existing backfill | VIEW.19, VIEW.24 |  |
| VIEW.24 | System sky S6: the sky tier table: the floor against distance from the visibility physics, registered in GEN.725 | VIEW.19, VIEW.20 |  |
| VIEW.23 | System sky S5: the sky catalogue galaxy/sky.py from stored stars: sky_from(observer, limit), excluding the observer system, with a tiers-reached report | VIEW.20, VIEW.21, VIEW.22 |  |
| VIEW.22 | System sky S4: sphere query by luminosity for charted stars, bright_stars, scattered phenomena and cluster members, on db/region.py |  |  |
| VIEW.21 | System sky S3: sky vectors on arrays with light-travel time: advance to the sky's time, retarded position, age at the retarded time, remnant visibility |  |  |
| VIEW.20 | System sky S2: photometry for the sky on arrays: apparent magnitude, the galactic dust law and the cool-star rule |  |  |
| VIEW.19 | System sky S1: measure the backfill cost of a 6.5 and a 7.5 magnitude sky on a real galaxy, and settle the luminosity-versus-mass ledger question |  |  |
| UX.47 | A bookmark manager: list, rename, sort, group and delete | UX.46 |  |
| UX.48 | Charted regions as bookmarks that frame and outline the region | UX.47 |  |
| UX.45 | Bookmark management | UX.46, UX.47, UX.48 | Parent; per browser until accounts (USR.7) move them. |
| UX.92 | A Bookmark button for planets, moons and belts | UX.45 |  |

### Robust tile caching

| ID | Item | Needs | Note |
|---|---|---|---|
| PERF.134 | The admin page shows the tile cache and can clear it: tiles per level, bytes, built and visited counts, the build queue, and clear and build-now actions |  | Boss, 2026-10-17, 4:31 PM EDT (tile caching, 8.2). |
| OPS.64 | The version code carries the database schema, generator epoch and API versions as well |  |  |
| PERF.232 | test_galaxy_seed::test_a_new_seed_replaces_the_old_only_before_any_sector takes over 10 minutes in its sectors step: cold per-process evolved-track evaluation per cell (bug) |  |  |
| PERF.236 | The prevalence draws with a forced habitable_world are too slow to finish: 600 systems did not finish in 5,000 s (about 0.9 s a system at baseline, several times that forced), and asteroid_belt +30 did not finish in 6,600 s (bug) |  |  |
| PERF.233 | One command builds the entire Galaxy Map tile cache from scratch, for the nightly maintenance run after the orbital update |  |  |
| PERF.237 | The warm-up's galaxy-nodes chunks took 37 s each right after an update although the nodes were already stored, against 4 s a chunk before it: find what a chunk does when nothing is built (bug) |  |  |
| PERF.235 | A dropped terminal ends the warm-up's progress view: add planetgen warm --follow to watch a running warm-up from the run record, and survive a hangup |  |  |
| PERF.234 | The warm-up's overall time left extrapolates the current step's speed over all remaining steps, so it reads hours when most steps are quick (bug) |  |  |
| DOC.103 | American spelling of stored values and API keys: cancelled, vapour, centre, flavour, ionising, oxidised, disc, link-neighbours |  |  |
| GEN.1247 | Outer asteroid belts store about 4.5e25 kg (7 Earth masses) across 11 to 74 AU at 6 to 9 Gy, far above a Kuiper-belt-like 0.03 Earth masses, while inner belts hold 1e20 to 1e24 kg: check the belt mass law at large distance and old age (bug) |  |  |
| TEST.562 | test_a_member_in_a_decayed_group_is_built_from_the_hierarchy_and_marks_no_group fails on main: 1 star where 2 are expected (bug) |  |  |
| GEN.1246 | The thin disc infall time scale follows the disc scale length for every galaxy mass, so downsizing is weak (young share 1.4 times higher at 1e10 than at 5e10 solar masses where it was 3): give the time scale its mass dependence from a source |  |  |
| GEN.1244 | Regular moons of giant planets form 46 to 130 giant radii out where the regular moons of the Solar System giants sit at 6 to 26: derive the regular-moon zone of a giant from its circumplanetary disc (bug) |  |  |
| GEN.1241 | Wire the system body view: the ledger and barycentre function read the typed tables for the system page |  |  |
| GEN.1240 | The greenhouse table is held at 1000 bar, so steam airs of 1e4 to 1e5 bar show 700 to 1300 K, too cool (bug) |  |  |
| GEN.1238 | A Venus twin (Venus mass at 0.72 AU, G2V) reaches 3,100 bar and 1,279 K against Venus's 92 bar and 737 K, with a tenth of the CO2 mass the pressure implies (bug) |  |  |
| GEN.1236 | Two class N planets sit at the 2000 K surface_air cap but keep the uncapped greenhouse factor, so stored temperature and greenhouse disagree (bug) |  |  |
| PERF.129 | Robust tile caching: the first zoom levels built ahead, tiles stored compressed, built proactively around every visit, kept until the data changes | PERF.134 | Boss, 2026-10-17, 4:31 PM EDT (tile caching, 8.2). |

## Phase 2: What a signed-in user gets

Account bookmarks, user API keys, the one-off system page, the home field and the umbrella.

Release gate (standing, Boss, 2026-10-21, 5:10 PM EDT): this version is not stamped until SEC.33 (the full security code audit) and TEST.261 (the full bug hunt audit) each have a dated report for it and no finding they mark blocking is open. Both are standing items: they stay open and are re-run for every major or minor version change (required whenever the major or minor version changes (8.1, 8.2 and on to 9.0 and 10.0), never for build numbers).

### User features

| ID | Item | Needs | Note |
|---|---|---|---|
| API.6 | User-level API keys, owned by the account that created them, that can read but not upload | USR.26 | Boss 01:31Z: every key belongs to the account that made it, so it follows USR.2. The API call logging he asked for is filed separately. |
| USR.7 | A user-level interface with bookmarks | USR.26 | Bookmarks of any object use NAV.7 references. |
| USR.27 | Admin Users page: list accounts, promote, demote, disable and delete | USR.26 | Split from USR.2. |
| USR.2 | Accounts with roles: user, moderator, admin and Owner | USR.26, USR.27 | Control schema; no blockers, placed late by priority. |
| USR.1 | User accounts | USR.2, USR.3, USR.4, USR.5, USR.6, USR.7 | Parent; closes with its subitems. |
| USR.8 | Every signed-in user can generate a one-off system | USR.26 | Boss 05:16Z. Today /admin/generate/system is admin-only and only admin accounts exist, so no earlier step is needed. Open question: per-user limit (default 30 an hour, admins unlimited). |
| USR.9 | `seo.privacy_note` text, "download my data" and "delete my account" on the account page | USR.7, USR.26 | Research suggestion; not requested. |
| USR.11 | A "home" field on user accounts: where the Galaxy Map zooms to when it opens | MAP.189 |  |

### The continuous zoom (headline, Boss 6:28 PM EDT 2026-10-28)

Boss: the zoom chain is a seminal feature of 8.2. Z0 (MAP.333, the quick fixes) is built in 8.1 and closes MAP.257; the rebuild below follows it. Owner: Core 2 for the map items, Core 5 for the preferences screen. Design: /mnt/project-files/design-drafts/galaxy-map-continuous-zoom-redesign-2026-10-28.md.

| ID | Item | Needs | Note |
|---|---|---|---|

### Halo mass research

| ID | Item | Needs | Note |
|---|---|---|---|

### Viewers and object pages

Boss (6:10 and 6:12 PM EDT, 2026-10-21): the viewing windows fit the browser window, and selecting anything opens the same info box with a View button to a dedicated page for every kind of item, with links to where it is.

| ID | Item | Needs | Note |
|---|---|---|---|
| PERF.204 | An admin's visit to an uncharted cell queues a one-sector generation run by itself, and the page says queued (the later option 2 of PERF.149 rule 2) |  |  |
| UX.143 | A dedicated page for every kind of item: star, planet, moon, belt, comet, each phenomenon, nebula, cluster, system, sector and region, reached from its Object ID |  | Boss's 8.2 design notes, 6:10 and 6:12 PM EDT 2026-10-21. |
| UX.144 | Every item page says where the item is: links up the chain (galaxy, sector, system, parent star or planet) and the nebula or other region that contains it | UX.143 | Boss's 8.2 design notes, 6:10 and 6:12 PM EDT 2026-10-21. |

### Search by physical values

| ID | Item | Needs | Note |
|---|---|---|---|

### Generation records

| ID | Item | Needs | Note |
|---|---|---|---|
| PERF.197 | Carry out the object timing decision: show the records with a test, or delete the recording, the tables, the migration and the fingerprint entry (G5, step 2) |  | Unused stored values audit (report of 10:33 PM EDT 2026-10-24), Boss 10:16 PM EDT: file each; chain PERF.195, Boss 10:35 PM EDT. |
| PERF.195 | Generation timing records are shown or no longer recorded: object_timing_stats and object_timing_samples (G5) | PERF.197 | Unused stored values audit (report of 10:33 PM EDT 2026-10-24), Boss 10:16 PM EDT: file each. |

## Phase 3: Settlements foundation

The design, storage and draw streams for nation-states, cities and ports on a planet with a technological civilization (Boss, 9:27 PM EDT 2026-10-24: earliest foundation in 8.2, nothing in 8.1). The generation is in 8.4, the editors in 8.6, the dynamic version in 9.2.

### Settlements and space-faring governments

| ID | Item | Needs | Note |
|---|---|---|---|
| POP.16 | Design: settlements and nation-states on a planet with a technological civilization (the tech gate, nation count and dominance laws, the city size law, the city to nation draw, space ports), with the sources and stand-in tiers |  | Boss, 9:27 PM EDT 2026-10-24: settlements tree; fixed until 9.2. |
| DB.103 | Storage for nation-states and cities: tables, Object IDs, column units and the migration | POP.16 | Boss, 9:27 PM EDT 2026-10-24: settlements tree; fixed until 9.2. |
| GEN.774 | Settlement draw streams: the named draw.py streams for nation count, dominance, city sizes, city names, city positions and the city to nation assignment, and their stand-in registry rows | POP.16 | Boss, 9:27 PM EDT 2026-10-24: settlements tree; fixed until 9.2. |
| POP.25 | Option: how many space-faring governments the galaxy has (a number, or natural), on the Generate page, the CLI, the API and the stored plan | POP.16 | Boss, 9:30 PM EDT 2026-10-24: space-faring governments, built back from 9.1, landing in 8.2. |
| DB.104 | Storage for a government's territory: homeworld system, growth state and a territory history, ready for 9.2 to make it dynamic | POP.16, DB.103 | Boss, 9:30 PM EDT 2026-10-24: space-faring governments, built back from 9.1, landing in 8.2. |
| POP.26 | Generate the requested number of space-faring governments: force civilizations where chance gives too few, hold back extra ones where it gives too many, homeworlds spread apart | POP.25, DB.104, GEN.774 | Boss, 9:30 PM EDT 2026-10-24: space-faring governments, built back from 9.1, landing in 8.2. |
| POP.27 | Territories of several worlds and sectors: every government holds a minimum, its homeworld is always inside, and the check is a test | POP.26 | Boss, 9:30 PM EDT 2026-10-24: space-faring governments, built back from 9.1, landing in 8.2. |
| POP.28 | Show the governments: the option and result on the Generate page and galaxy page, a government page with its worlds, systems and sectors, and the territories overlay legend | POP.27 | Boss, 9:30 PM EDT 2026-10-24: space-faring governments, built back from 9.1, landing in 8.2. |
| POP.29 | Design: the world population simulation and the trade routes between spaceports (partners, hop runs, the global list, map filters, nearest-route query), with sources and stand-in tiers | POP.16 | Boss, 9:38 PM EDT 2026-10-24: trade routes and world populations, done by 8.9, dynamic in 9.1. |
| DB.105 | Storage for world population and the global trade route list: tables, Object IDs, column units and the migration | POP.29, DB.103 | Boss, 9:38 PM EDT 2026-10-24: trade routes and world populations, done by 8.9, dynamic in 9.1. |

### Resource locations

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.798 | Design: resource locations (what counts as a valuable resource, how a belt, field, giant or icy body is scored, the sources and the stand-in tiers) |  | Boss, 10:29 PM EDT 2026-10-24: resource locations (valuable belts and fields as points of interest), foundations 8.2 to 9.0. |
| DB.106 | Storage for resource locations: grade, mass estimate and score per resource for belts, fields, giants and icy bodies, with column units and the migration | GEN.798 | Boss, 10:29 PM EDT 2026-10-24: resource locations (valuable belts and fields as points of interest), foundations 8.2 to 9.0. |
