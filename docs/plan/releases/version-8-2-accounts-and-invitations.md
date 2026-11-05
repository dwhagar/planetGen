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
| VIEW.30 | System sky S12: calibration test against the HYG catalogue on a developer machine: a generated sky has about 9,000 stars at 6.5 | VIEW.25 |  |
| VIEW.26 | System sky S7b: a one-off command for sectors charted earlier: build the sky tiers on demand (Generate console and CLI) | VIEW.25 |  |
| VIEW.25 | System sky S7 (lock): charting builds the sky tiers: when a system is charted, the cells its sky needs are built with the existing backfill | VIEW.19 |  |
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
| PERF.129 | Robust tile caching: the first zoom levels built ahead, tiles stored compressed, built proactively around every visit, kept until the data changes | PERF.134 | Boss, 2026-10-17, 4:31 PM EDT (tile caching, 8.2). |

## Phase 2: What a signed-in user gets

Account bookmarks, user API keys, the one-off system page, the home field and the umbrella.

Release gate (standing, Boss, 2026-10-21, 5:10 PM EDT): this version is not stamped until SEC.33 (the full security code audit) and TEST.261 (the full bug hunt audit) each have a dated report for it and no finding they mark blocking is open. Both are standing items: they stay open and are re-run for every major or minor version change (required whenever the major or minor version changes (8.1, 8.2 and on to 9.0 and 10.0), never for build numbers).

### User features

| ID | Item | Needs | Note |
|---|---|---|---|
| API.6 | User-level API keys, owned by the account that created them, that can read but not upload | USR.26 | Boss 01:31Z: every key belongs to the account that made it, so it follows USR.2. The API call logging he asked for is filed separately. |
| USR.7 | A user-level interface with bookmarks | USR.26 | Bookmarks of any object use NAV.7 references. |
| UX.265 | Look and theme preference: a look (cockpit, chart, phosphor) times a theme (system, light, dark), stored on the account and in the browser, chosen in Preferences and the gear menu | USR.7 |  |
| UX.277 | Map layer registry fields for the View menu: heading, keywords, common flag and a settings hook, every existing layer given its heading, and a guard test that fails for a layer with none |  |  |
| UX.278 | Shared View menu builder for the Galaxy, Sector and System Maps: DOM calls only, headings with counts, collapse state kept in the visitor preferences, and a mark where a setting differs from its default | UX.277 |  |
| UX.279 | Search filter for the View menu: a function with tests that filters entries by name and keywords and opens the matching groups, behind no visible control yet | UX.277 |  |
| UX.275 | Sensor registry and API: one registry of kind, look builder and painter (web/lib/sensor.py, static/sensorview.js) and GET /api/objects/<ref>/sensor, with a test that fails for any kind in describe/sheet.py that has no entry |  |  |
| UX.276 | Object page: the rendered example becomes the Sensor view, the first panel for every kind, with a Sensors show readout line and estimate labels | UX.275 |  |
| MAP.394 | Map info panel: selecting an object on the Galaxy, Sector or System Map opens its Sensor view in the info panel | UX.275, UX.276 |  |
| MAP.395 | Sensor pictures for the kinds with none: belt as a grain field, comet nucleus with coma and tail, rogue planet as a cold world with internal heat | UX.275 |  |
| MAP.396 | Sensor pictures for systems and for sectors, regions, galaxies and clusters: a top-down orrery plate and a star-field and gas patch from stored counts | UX.275 |  |
| TEST.629 | Browser test: every kind in the property sheet shows a Sensor view on its page and in the map info panel, with the readout line and the estimate label where an input is not stored | UX.276, MAP.394, MAP.395, MAP.396 |  |
| DOC.107 | Docs: the Sensor view (what it shows for each kind, the estimate labels, the API route) in the help and the design notes | UX.276, MAP.394, MAP.395, MAP.396 |  |
| UX.273 | Sensor (Rendered) view for every item: a rendered picture of the selected item beside the chart views | UX.275, UX.276, MAP.394, MAP.395, MAP.396, TEST.629, DOC.107 |  |
| TEST.628 | A docs test fails when an item's Version line disagrees with the version file its row sits in |  |  |
| UX.130 | Console design brief: colours, type, panel shapes and iconography, approved by Boss |  | Split from UX.43. |
| UX.131 | Console tokens, type and base template in both themes | UX.130, UX.265 | Split from UX.43. |
| USR.27 | Admin Users page: list accounts, promote, demote, disable and delete | USR.26 | Split from USR.2. |
| USR.2 | Accounts with roles: user, moderator, admin and Owner | USR.26, USR.27 | Control schema; no blockers, placed late by priority. |
| USR.1 | User accounts | USR.2, USR.3, USR.4, USR.5, USR.6, USR.7 | Parent; closes with its subitems. |
| USR.8 | Every signed-in user can generate a one-off system | USR.26 | Boss 05:16Z. Today /admin/generate/system is admin-only and only admin accounts exist, so no earlier step is needed. Open question: per-user limit (default 30 an hour, admins unlimited). |
| USR.9 | `seo.privacy_note` text, "download my data" and "delete my account" on the account page | USR.7, USR.26 | Research suggestion; not requested. |
| USR.11 | A "home" field on user accounts: where the Galaxy Map zooms to when it opens | MAP.189 |  |
| USR.12 | Permission catalogue and one `can()` check: user, moderator, admin and Owner as data | USR.2 | Roles from 8.2 (USR.2 already holds the moderator value) |
| USR.13 | Moderator role: assign and remove it on the admin accounts page, with the audit entry | USR.12, USR.27 |  |
| USR.14 | Move every route, page and command from "is admin" checks to permission checks | USR.12 |  |

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

### Search by physical values

| ID | Item | Needs | Note |
|---|---|---|---|

### Generation records

| ID | Item | Needs | Note |
|---|---|---|---|

### Gravity grid foundations (Boss, 5:14 PM EST on 2026-11-05)

| ID | Item | Needs | Note |
|---|---|---|---|
| PERF.72 | Research the cost of a sector's gravity grid and where to cut between exact and aggregated sources |  | Research; the study is done (reports/perf72-gravity-grid-cost-2026-11-05.md), the note in gravity-map.md is still to be extended. |

### Sky foundations for the view from a planet (Boss (5:44 PM EST on 2026-11-05): the view's design is approved, it lands in 9.0, so its foundations are built back to 8.2)

| ID | Item | Needs | Note |
|---|---|---|---|
| VIEW.4 | Constellation names in the name generator |  | Constellation names from the codec under the naming key, not sliced word lists. Floats: no blockers, can run any time. |
| VIEW.6 | Settle the handedness of the generated galaxy before mapping real sky coordinates |  | Research: needed before real sky coordinates. |
| VIEW.1 | View from a planet |  | Research session with Boss. |
| VIEW.7 | Sky tiers and the backfill cost measurement | VIEW.1 | Research: the sky is about 8% to 10% complete today. |
| VIEW.8 | A dust model for the sky: the galactic extinction law and `nebulae.extinction_av` | VIEW.1 | Research. |
| VIEW.11 | Photometry library: apparent magnitude, bolometric correction, blackbody colour and the extinction law | VIEW.1 | Split from VIEW.2. Extends GEN.416 (the blackbody colour, 9.0). |
| VIEW.12 | Sky catalogue: the stars an observer sees, with direction, brightness and light-time | VIEW.11, VIEW.7 | Split from VIEW.2. |
| VIEW.13 | Galactic band: the glow of the unresolved stars and dust along each line of sight | VIEW.11, VIEW.8 | Split from VIEW.2. |
| VIEW.43 | Dust model built from the galaxy's own gas: planetgen/physics/dust_model.py with the column and layer-thickness tables, and the stand-ins DU1 to DU5 |  |  |
| VIEW.44 | The sky reads the gas-built dust model: galaxy/sky.py, galaxy/sky_band.py and generation/sky_tiers.py take the two tables through the dust dict, and the table hash is in their cache keys | VIEW.43 |  |
| VIEW.45 | Re-baseline the HYG calibration test (VIEW.30) for the new dust, in dust-free and default-dust modes, in the small and the Milky Way sized galaxy | VIEW.44, VIEW.30 |  |
| GEN.1290 | Gas and star metallicity depend on the galaxy's stellar mass (the mass-metallicity relation), so small galaxies are metal-poor and dustier-poor |  |  |
| GEN.1291 | The gas metallicity gradient scales with the galaxy: 0.088 dex per disk scale length in place of 0.034 dex per kpc (row NBC1) |  |  |

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
