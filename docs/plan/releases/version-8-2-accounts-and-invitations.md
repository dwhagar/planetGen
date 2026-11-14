# Version 8.2: User accounts, bookmarks and API keys, with the foundations for later versions

Rebuilt on 2026-10-16 as a version map (Boss, 1:07 PM EDT: one version per headline, phases inside a version, blockers one phase earlier).
[docs/TODO.md](../../TODO.md) is the master file: it holds each item's full text, and its "Plan: versions" section indexes every version.
This file gives the version's goal, its phases and each phase's build threads with every item's prerequisites in order. A phase holds independent threads: an item's blockers sit in the same thread (earlier in the list) or in an earlier phase. Where the two disagree, TODO.md wins; when an item ships, it leaves TODO.md and its row here is deleted in the same PR (the progress ledger keeps its version). Research notes and judgment calls are in [notes.md](../notes.md); the headlines and the ideas not yet filed are in [release-roadmap.md](../release-roadmap.md).

Version: 8.2

Status: paused until 8.1 is complete (Boss, 7:44 PM EDT 2026-10-22: the five cores work 8.1 only; nothing in this version is started or owned until then, except the gates already running for UX.140, UX.128 and NAV.17).

## Goal

Boss (5:36 PM EST on 2026-11-05): user accounts, bookmarks and API keys are the main features of 8.2; everything else in it is the foundations for later versions (the continuous zoom, tile caching, the settlements design and storage, the gravity grid). The description that follows is the detail of those features.

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
| VIEW.40 | System sky S22 (optional): cross-fade between the Galaxy Map's stars and the sky at the system, after the continuous zoom |  |  |
| VIEW.39 | System sky S21 (optional): aberration of the observer's motion, default off |  |  |
| VIEW.26 | System sky S7b: a one-off command for sectors charted earlier: build the sky tiers on demand (Generate console and CLI) |  |  |
| UX.47 | A bookmark manager: list, rename, sort, group and delete | UX.46 |  |
| UX.48 | Charted regions as bookmarks that frame and outline the region | UX.47 |  |
| UX.45 | Bookmark management | UX.46, UX.47, UX.48 | Parent; per browser until accounts (USR.7) move them. |
| UX.92 | A Bookmark button for planets, moons and belts | UX.45 |  |

### Robust tile caching

| ID | Item | Needs | Note |
|---|---|---|---|
| PERF.134 | The admin page shows the tile cache and can clear it: tiles per level, bytes, built and visited counts, the build queue, and clear and build-now actions |  | Boss, 2026-10-17, 4:31 PM EDT (tile caching, 8.2). |
| OPS.64 | The version code carries the database schema, generator epoch and API versions as well |  |  |

## Phase 2: What a signed-in user gets

Account bookmarks, user API keys, the one-off system page, the home field and the umbrella.

Release gate (standing, Boss, 2026-10-21, 5:10 PM EDT): this version is not stamped until SEC.33 (the full security code audit) and TEST.261 (the full bug hunt audit) each have a dated report for it and no finding they mark blocking is open. Both are standing items: they stay open and are re-run for every major or minor version change (required whenever the major or minor version changes (8.1, 8.2 and on to 9.0 and 10.0), never for build numbers).

### User features

| ID | Item | Needs | Note |
|---|---|---|---|
| API.6 | User-level API keys, owned by the account that created them, that can read but not upload | USR.26 | Boss 01:31Z: every key belongs to the account that made it, so it follows USR.2. The API call logging he asked for is filed separately. |
| USR.7 | A user-level interface with bookmarks | USR.26 | Bookmarks of any object use NAV.7 references. |
| UX.265 | Look and theme preference: a look (cockpit, chart, phosphor) times a theme (system, light, dark), stored on the account and in the browser, chosen in Preferences and the gear menu | USR.7 |  |
| API.36 | Optional: a region_names search table and a "Go to" box, refilled when the naming key changes (regions chain, step 6) | API.35 | Boss, 2:21 PM EST on 2026-11-07: slated for 8.2. |
| ADM.139 | In the admin console, an option to look at and reset manually edited items that were generated |  | Boss, 2:21 PM EST on 2026-11-07, while answering GEN.229. |
| MAP.394 | Map info panel: selecting an object on the Galaxy, Sector or System Map opens its Sensor view in the info panel |  |  |
| MAP.395 | Sensor pictures for the kinds with none: belt as a grain field, comet nucleus with coma and tail, rogue planet as a cold world with internal heat |  |  |
| MAP.396 | Sensor pictures for systems and for sectors, regions, galaxies and clusters: a top-down orrery plate and a star-field and gas patch from stored counts |  |  |
| TEST.629 | Browser test: every kind in the property sheet shows a Sensor view on its page and in the map info panel, with the readout line and the estimate label where an input is not stored | MAP.394, MAP.395, MAP.396 |  |
| DOC.107 | Docs: the Sensor view (what it shows for each kind, the estimate labels, the API route) in the help and the design notes | MAP.394, MAP.395, MAP.396 |  |
| UX.273 | Sensor (Rendered) view for every item: a rendered picture of the selected item beside the chart views | MAP.394, MAP.395, MAP.396, TEST.629, DOC.107 |  |
| UX.286 | Sensor view: a dotted line through the middle of a planet, moon or star along its axis of rotation, read from the stored spin axis and axial tilt |  |  |
| UX.287 | Sensor view: magnetic field lines drawn for a planet, moon, star or neutron star from its stored field strength and dipole orientation, when it has a field |  |  |
| UX.130 | Console design brief: colours, type, panel shapes and iconography, approved by Boss |  | Split from UX.43. |
| UX.131 | Console tokens, type and base template in both themes | UX.130, UX.265 | Split from UX.43. |
| USR.27 | Admin Users page: list accounts, promote, demote, disable and delete | USR.26 | Split from USR.2. |
| USR.2 | Accounts with roles: user, moderator, admin and Owner | USR.26, USR.27 | Control schema; no blockers, placed late by priority. |
| USR.1 | User accounts | USR.2, USR.3, USR.4, USR.5, USR.6, USR.7 | Parent; closes with its subitems. |
| USR.8 | Every signed-in user can generate a one-off system | USR.26 | Boss 05:16Z. Today /admin/generate/system is admin-only and only admin accounts exist, so no earlier step is needed. Open question: per-user limit (default 30 an hour, admins unlimited). |
| USR.11 | A "home" field on user accounts: where the Galaxy Map zooms to when it opens | MAP.189 |  |
| USR.12 | Permission catalogue and one `can()` check: user, moderator, admin and Owner as data | USR.2 | Roles from 8.2 (USR.2 already holds the moderator value) |
| USR.14 | Move every route, page and command from "is admin" checks to permission checks | USR.12 |  |

### The continuous zoom (headline, Boss 6:28 PM EDT 2026-10-28)

Boss: the zoom chain is a seminal feature of 8.2. Z0 (MAP.333, the quick fixes) is built in 8.1 and closes MAP.257; the rebuild below follows it. Owner: Lane 12 for the map items, Lane 15 for the preferences screen. Design: /mnt/project-files/design-drafts/galaxy-map-continuous-zoom-redesign-2026-10-28.md.

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.141 | Context around the selection: faint neighbours, and the sectors above and below |  | GitHub issues [#718](https://github.com/dwhagar/planetGen/issues/718) and [#716](https://github.com/dwhagar/planetGen/issues/716), one context view. The #716 bug label was overruled. |
| MAP.156 | View one layer or a range of layers top-down from the galaxy view, as a secondary option |  | Boss 2026-10-09 23:17Z; secondary option, not a main control. |

### Halo mass research

| ID | Item | Needs | Note |
|---|---|---|---|

### Viewers and object pages

Boss (6:10 and 6:12 PM EDT, 2026-10-21): the viewing windows fit the browser window, and selecting anything opens the same info box with a View button to a dedicated page for every kind of item, with links to where it is.

| ID | Item | Needs | Note |
|---|---|---|---|

### Search by physical values

| ID | Item | Needs | Note |
|---|---|---|---|

### Generation records

| ID | Item | Needs | Note |
|---|---|---|---|

### Gravity grid foundations (Boss, 5:14 PM EST on 2026-11-05)

| ID | Item | Needs | Note |
|---|---|---|---|

### Sky foundations for the view from a planet (Boss (5:44 PM EST on 2026-11-05): the view's design is approved, it lands in 9.0, so its foundations are built back to 8.2)

| ID | Item | Needs | Note |
|---|---|---|---|
| VIEW.4 | Constellation names in the name generator |  | Constellation names from the codec under the naming key, not sliced word lists. Floats: no blockers, can run any time. |
| VIEW.6 | Settle the handedness of the generated galaxy before mapping real sky coordinates |  | Research: needed before real sky coordinates. |
| VIEW.8 | A dust model for the sky: the galactic extinction law and `nebulae.extinction_av` |  | Research. |
| VIEW.11 | Photometry library: apparent magnitude, bolometric correction, blackbody colour and the extinction law |  | Split from VIEW.2. Extends GEN.416 (the blackbody colour, 9.0). |
| VIEW.12 | Sky catalogue: the stars an observer sees, with direction, brightness and light-time | VIEW.11 | Split from VIEW.2. |
| VIEW.13 | Galactic band: the glow of the unresolved stars and dust along each line of sight | VIEW.11, VIEW.8 | Split from VIEW.2. |

### Slated, not numbered: the rare and fictional materials (Boss, 5:26 PM EST on 2026-11-07)

Boss deferred the fictional materials and crystals (Kyber, Solar, Aurora, Tritanium, Dilithium, Pergium, Topaline, Kelbonite and Neutronium; Duranium, Trellium and Protomatter are out of scope) until 8.2 is ready to start, to be revisited after 8.1 is done. Nothing is numbered or scheduled. The plan is /mnt/project-files/research/materials-8-2-plan-2026-11-07.md (provisional parts MAT-1 to MAT-12, with 7 and 10 removed, and five open questions for Boss in its section 9). When 8.2 opens, the TODO thread assigns item numbers and places the parts as their own generator-side phase. They read the real rows the 8.1 real-material set builds (GEN.1355 to GEN.1383), and the Solar Crystal's carbon-13 and carbon-14 choices read the isotope record GEN.1377.

### Unified chain rebuild (moved from 8.1 by Boss, 9:19 PM EST on 2026-11-11)

The full recursion and the code and storage rebuild of the unified orbited-body chain (GEN.550), after the 8.1 fixes are in. Cross-check: reports/unified-chain-crossref-2026-11-12/.

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.565 | The mass object core and its three family classes: StellarObject, PlanetaryObject and SmallBody, with state and parent as data, replace the nine kind classes |  |  |
| DB.74 | One generic insert and load driven by the column list, one registry of kind tables, and the shared core table with three profile tables | GEN.565 |  |
| GEN.570 | Field declaration: every unit-valued field of the mass object classes declares its kind and native unit, and the one field table serves the serializer | GEN.565 |  |
| GEN.575 | Sweep the mass object classes: the 35 unit-suffixed attributes and every quantity held under two units get one declared native unit | GEN.565, GEN.570 |  |
| GEN.584 | The core's kinematic record: velocity, spin axis and rate, mass and epoch in the point, and the derived columns dropped | GEN.578, GEN.565 |  |
| GEN.581 | primary_first(members): before a system is validated the heaviest branch becomes the primary and the others are recalculated, one function for the tree, the fill, moons and free pairs | GEN.315 |  |
| GEN.578 | The point class gets the local frame and the barycentre system frame, and the system anchor is no longer the nearest star | GEN.576 |  |
| GEN.577 | Stored planet and moon states become the state of the barycentre of the body and its own members | GEN.576 |  |
| GEN.576 | system(host): the one function that groups bodies by their parent link and returns the barycentre sums, the central object's offset and the valid-system check | GEN.565, GEN.315 |  |
| GEN.587 | Replace the hand-written relation sites with the relations module, by file group, bit-identical on the golden set |  |  |
| GEN.571 | Replace the ad-hoc unit conversions with the registry, by file group, bit-identical on the golden set |  |  |
| GEN.555 | The pair builder for any masses: planet pairs, brown dwarf pairs and a planet with its moons share one stability test | GEN.315, GEN.331, GEN.332, GEN.581 |  |
| GEN.550 | Umbrella (placeholder): one recursive generation process for anything that can have something else orbiting it, from a planet with moons up to a group of star-like objects that each have a system |  |  |
| GEN.559 | Unbound comets: a size law and a cut, composition from the comet table, derived activity and velocity |  |  |
| GEN.556 | `populate(host)`: planets, belts, comets and moons come from the class tables by one recursive routine, and the system pass is absorbed | GEN.315, GEN.555, GEN.565, GEN.581 |  |
| GEN.349 | API output and stored-column readers take Quantity; the grep test fails on a hand-written conversion factor; the old helpers are removed |  | Split from GEN.300. Moved from 8.3. |
| TEST.240 | Round-trip and registry tests for the mass object core: every kind's stored row round-trips, and no kind table is named in a list the registry should own | GEN.565, DB.74 |  |
| TEST.239 | Tests for the unified chain: determinism by host, representation counts per kind and level, rogue and comet rates, derived remnant shares | GEN.554, GEN.558, GEN.559 |  |
| DOC.71 | Docs for the System rule: database-schema.md, the system view in docs/design and the header of system_tree.py | GEN.576, DB.78, GEN.578 |  |
| MAP.253 | Map and scene read the parent chain: systemscene, body_positions, the browser twin orbitpositions.js and the systemmap branches drop the three origin conventions | GEN.576, DB.78, GEN.577 |  |
| DB.78 | Central-object offsets replace the wobble: the reflex_offset columns, calculate_reflex_offset, the SQL mirror and the nine binary position columns go | GEN.576 |  |

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
| POP.26 | Generate toward the requested number of space-faring governments by tilting the odds that a species becomes space-faring (nothing forced, nothing cancelled), homeworlds spread apart, and report what came out | POP.25, DB.104, GEN.774 | Boss, 9:30 PM EDT 2026-10-24: space-faring governments, built back from 9.1, landing in 8.2. |
| POP.27 | Territories of several worlds and sectors: every government's homeworld is inside its territory, a government that really reaches only its own system holds one system and the page says why, and the check is a test | POP.26 | Boss, 9:30 PM EDT 2026-10-24: space-faring governments, built back from 9.1, landing in 8.2. |
| POP.28 | Show the governments: the option and result on the Generate page and galaxy page, a government page with its worlds, systems and sectors, and the territories overlay legend | POP.27 | Boss, 9:30 PM EDT 2026-10-24: space-faring governments, built back from 9.1, landing in 8.2. |
| POP.29 | Design: the world population simulation and the trade routes between spaceports (partners, hop runs, the global list, map filters, nearest-route query), with sources and stand-in tiers | POP.16 | Boss, 9:38 PM EDT 2026-10-24: trade routes and world populations, done by 8.9, dynamic in 9.1. |
| DB.105 | Storage for world population and the global trade route list: tables, Object IDs, column units and the migration | POP.29, DB.103 | Boss, 9:38 PM EDT 2026-10-24: trade routes and world populations, done by 8.9, dynamic in 9.1. |

### Resource locations

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.798 | Design: resource locations (what counts as a valuable resource, how a belt, field, giant or icy body is scored, the sources and the stand-in tiers) |  | Boss, 10:29 PM EDT 2026-10-24: resource locations (valuable belts and fields as points of interest), foundations 8.2 to 9.0. |
| DB.106 | Storage for resource locations: grade, mass estimate and score per resource for belts, fields, giants and icy bodies, with column units and the migration | GEN.798 | Boss, 10:29 PM EDT 2026-10-24: resource locations (valuable belts and fields as points of interest), foundations 8.2 to 9.0. |

### Fictional materials (Boss, 2026-11-08)

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.1408 | Fictional materials part 1: the fiction table (the four fictional materials only: Kyber crystal, Solar crystal, Aurora crystal and Tritanium; each with its kind, default inclusion and default name) |  | Boss, 10:28 AM EST on 2026-11-08. |
| GEN.1409 | Fictional materials part 2: name rules (a custom name per material, allowed characters, length, uniqueness, and the fallback to the default name) | GEN.1408 | Boss, 10:28 AM EST on 2026-11-08. |
| GEN.1410 | Fictional materials part 3: the settings object and its config default (per-material include flag and custom name) | GEN.1409 | Boss, 10:28 AM EST on 2026-11-08. |
| DB.116 | Fictional materials part 4: the galaxy_fiction table and the plan write (the chosen settings are stored with the galaxy at generation time; migration) | GEN.1410 | Boss, 10:28 AM EST on 2026-11-08. |
| ADM.143 | Fictional materials part 5: the Generate page fields, the CLI flags and the docs table for the per-material setting and custom names | GEN.1410 | Boss, 10:28 AM EST on 2026-11-08. |
| API.42 | Fictional materials part 6: read routes for a galaxy's fiction settings and the galaxy page line | DB.116 | Boss, 10:28 AM EST on 2026-11-08. |
| GEN.1411 | Fictional materials part 7: the biomass index (the gate for the three life-gated crystals) |  | Boss, 10:28 AM EST on 2026-11-08. |
| GEN.1412 | Fictional materials part 8: the material label resolver (shows the custom name, or the default, or nothing when the material is off) | GEN.1410 | Boss, 10:28 AM EST on 2026-11-08. |
| GEN.1413 | Fictional materials part 9: the crystal occurrence stage (draws where Kyber, Solar and Aurora crystal occur; takes the epoch) | GEN.1411, GEN.1412 | Boss, 10:28 AM EST on 2026-11-08. |
| DB.117 | Fictional materials part 10: the per-body rare-material table (migration) | GEN.1413 | Boss, 10:28 AM EST on 2026-11-08. |
| GEN.1414 | Fictional materials part 11: the crystal rows and their GEN.725 registry rows | GEN.1408 | Boss, 10:28 AM EST on 2026-11-08. |
| GEN.1415 | Fictional materials part 12: Tritanium rows and the reduced-body ore law | GEN.1408 | Boss, 10:28 AM EST on 2026-11-08. |
| UX.293 | Fictional materials part 15: pages, body JSON and the search filter show the included materials with their chosen names | GEN.1412, DB.117 | Boss, 10:28 AM EST on 2026-11-08. |
| ADM.144 | Fictional materials part 17: edit custom names after generation | DB.116, GEN.1412, ADM.145 | Boss, 10:28 AM EST on 2026-11-08. |
| ADM.145 | Fictional materials admin page A: list the four materials with default state and name from config.json and the galaxy's stored state | GEN.1410, DB.116, API.42 | Boss, 10:41 AM EST on 2026-11-08. |
| ADM.146 | Fictional materials admin page B: set the config.json defaults (on or off and default name per material) through the settings writer, with backup and audit | ADM.145, GEN.1410 | Boss, 10:41 AM EST on 2026-11-08. |
| ADM.147 | Fictional materials admin page C: show how many bodies hold each material | ADM.145, DB.117 | Boss, 10:41 AM EST on 2026-11-08. |
| TEST.699 | Fictional materials part 19: the toggle invariance and off-equals-absent suite (a galaxy with a material off is identical to one made without it) | GEN.1413 | Boss, 10:28 AM EST on 2026-11-08. |
| DOC.110 | Fictional materials part 20: reference docs and help for the setting | ADM.143 | Boss, 10:28 AM EST on 2026-11-08. |

### Flavor text settings and signals (Boss, 2026-11-08)

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.1435 | Ancient works setting: tag the derelict, wreck, megastructure and signal lines (57 of 82) and a per-galaxy on/off choice at galaxy planning, on the same screen and storage as the fictional materials setting | GEN.1430 | Boss, 2026-11-08. |
| GEN.1436 | Natural signal resolution: sensor readings for unknown-origin signals that resolve to a cause (pulsars, giants with volcanic moons, masers) | GEN.1430 | Boss, 2026-11-08. |
