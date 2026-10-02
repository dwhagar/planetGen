# Phase 3: Interface, API and work queue

This is one phase of the plan drawn up on 2026-10-01 from Boss's list
of that evening and his research notes. [docs/TODO.md](../TODO.md) is
the master file: it holds each item's full text (what's wrong, where to
look, what "done" means), and its "Plan: phases" section indexes every
phase. This file holds what the phase needs beyond that: its goal, the
order and dependencies, how to split it into build threads, and the
research notes that apply. Where the two disagree, TODO.md wins; when an
item ships, it leaves TODO.md and its row here is deleted in the same PR.

## Goal

Modernize the pages (actions behind menus, icons, structured planet
details instead of the Markdown render, simpler rows, real phenomena
filters, units everywhere), build remote generation through the API,
and move heavy work onto the work queue with smarter backfill.

## Items

### Page modernization

| ID | Item | Parent |
|---|---|---|
| UX.28 (new) | Investigate icons instead of words on buttons |  |
| UX.26 (new) | Edit and admin actions as a button that opens a menu (bug) |  |
| UX.31 (new) | Editing a star system: an edit button with a quick menu, not a long panel (bug) | UX.26 |
| UX.27 (new) | System page: the system and navigation buttons on one row that doesn't overlap (bug) |  |
| UX.25 (new) | Rogue planets: octant and a small map symbol beside each name (bug) |  |
| UX.24 (new) | Sector contents: rogue planets after systems and phenomena, and expanded rows the full table width (bug) |  |
| UX.29 (new) | Every comet in a system shows its type as a link (bug) |  |
| UX.32 (new) | Planet rows show the class only, without the type and moon labels |  |
| UX.30 (new) | Planet information without the Markdown render |  |
| UX.33 (new) | Filter phenomena by their classes and types (bug) |  |
| UX.2 | Menus sized to what they hold (bug) |  |
| UX.3 | Warn every visitor while a background job changes the galaxy |  |
| ADM.14 | Line up the Generate page's text boxes, not their headings (bug) |  |

UX.28 (which buttons become icons, and the icon set) goes first, since
UX.25, UX.26, UX.27 and UX.31 use its icons. UX.26 (actions behind a
menu button) before UX.31 (the system page's case). UX.24, UX.25 and
UX.29 are small page fixes. UX.30 (planet details as structured HTML)
is the biggest; build it after UX.22's ladders so it uses them, and with
UX.32 (rows show the class only), which changes the same rows. UX.33
needs new API filters and `queryDb` columns. UX.2, UX.3 and ADM.14 are
independent.

### Units

| ID | Item | Parent |
|---|---|---|
| UX.22 | Meaningful units for every measurement |  |
| UX.23 | A shared unit-ladder module | UX.22 |

UX.23 (the shared ladder module) first, then one quantity family per
PR.

### Remote generation through the API

| ID | Item | Parent |
|---|---|---|
| API.4 | API compatibility data in the docs |  |
| API.5 | API version and compatibility checking |  |
| API.6 | User-level API keys, owned by the account that created them, that can read but not upload |  |
| API.9 | Key scopes | API.3 |
| API.7 | Investigate and plan upload limits |  |
| API.3 | Remote generate: generate on a local machine, upload through the API |  |
| API.10 | Reservations: claimed sectors and id blocks per run | API.3 |
| API.11 | Staging tables | API.3 |
| API.12 | The download: seed, skeleton and name state | API.3 |
| API.13 | Generation without a database | API.3 |
| API.14 | Upload routes, compressed, in batches | API.3 |
| API.8 | Verify uploaded data before it is finalized |  |
| ADM.13 | Incomplete uploads page |  |

API.4 (compatibility data) and API.5 (version checking) first, then
key scopes (API.6, API.9; API.6's link from user-level keys to user
accounts builds after USR.2, so that part waits for phase 4's accounts)
and upload limits (API.7), then API.3's
pieces (reservations, staging tables, the download, generation without
a database, compressed batch uploads), verification (API.8) and the
incomplete uploads page (ADM.13).

### Work queue and backfill

| ID | Item | Parent |
|---|---|---|
| PERF.19 | Everything the API or web site starts runs on the work queue (investigate) |  |
| PERF.18 | Run the GEN.30 bright-star backfill in parallel on the work queue |  |
| PERF.20 | Short-term caching through the work queue and API (needs planning) |  |
| ADM.15 | Change the worker count from the Queue page, with a "Ludicrous Speed" mode |  |
| GEN.40 | Weed out sectors by star density before the bright-star backfill |  |
| GEN.41 | Investigate: how much backfill work a density pre-pass would save | GEN.40 |
| GEN.42 | A pass that drops sectors from a region by probability | GEN.40 |
| GEN.43 | Don't over-filter: keep bright stars in odd places | GEN.40 |
| PERF.1 | Generation at scale |  |

PERF.19 (audit what the API and site start, inline or queued) first;
then PERF.18 (the backfill in parallel on the queue), PERF.20 (short-term
caching, planned with Boss) and ADM.15 (worker count and Ludicrous
Speed). GEN.40's investigation (GEN.41) decides whether GEN.42 and
GEN.43 are built; all of them use phase 1's GEN.44 levels to skip
finished sectors, and PERF.11's stored densities (both in phase 1's
per-sector stats table). Build PERF.18 and
GEN.42 in one thread: both change `backfill_bright_stars_around`.

### Reproducible galaxies, finished

| ID | Item | Parent |
|---|---|---|
| ADM.17 (new) | The Generate page shows the galaxy's seed and version | GEN.55 |
| API.16 (new) | The API reports the galaxy's seed, version and run history | GEN.55 |
| API.17 (new) | Remote generation reproduces what the server would make | GEN.55 |
| GEN.59 (new) | Edits and time evolution recorded as layers on top of the seed | GEN.55 |

Phase 1 builds the seeds, storage, fingerprint and golden test (DB.6,
OPS.10, OPS.11, GEN.56 to GEN.58, DB.7, TEST.77). Here ADM.17 and API.16
show the seed and version (API.12's download uses API.16's fields),
API.17 makes remote generation reproduce the server's (after API.12 and
API.13; API.8 re-runs a sample to verify uploads), and GEN.59 records
edits and time evolution as layers. OPS.12 (`generate.py reproduce`,
the end goal) follows in phase 4. GEN.42 and PERF.18 in
this phase must use the derived seeds and keep TEST.77 green.

## Research notes

Boss's research notes (kept in the project's shared files under `todo-tasks/research/`) proposed fixes and numbers. They were checked against the code on 2026-10-01; where they were wrong about the code, the correction is given. Their numbers are starting points to tune, not requirements.

- **Structured planet details (UX.30).** The notes sketch a planet card:
  a header with class, sector and habitability; columns for physical
  properties, an atmosphere profile with a bar per gas, and satellites;
  and an actions row. They suggest HTML custom elements (Web
  Components). The site renders on the server with Jinja today, so
  server-built HTML with the same layout fits better unless a part
  needs to update in the browser.
- **Action menus (UX.26, UX.31).** One "Actions" button per row opening
  a floating menu that doesn't move the page, with entries such as Edit,
  Regenerate, Backfill sector and Delete; editing opens a focused dialog
  with Save and Cancel.
- **One-row header (UX.27).** Details and actions on one row; when it
  is too narrow, the navigate buttons fold into one "Navigation" menu
  with "From here" and "To here". The notes say below 768 px; the site's
  UX rules use container queries and the size classes in TODO.md's UX
  section, so use those.
- **Icons (UX.28).** Starting list: edit (pencil), delete (trash), show
  on map (crosshair), filter (funnel), navigate (route or compass),
  actions menu (ellipsis or gear), each with `title` and `aria-label`.
- **Rows (UX.32).** One "Class" chip (for example "Class M") and a moon
  count badge that opens the moons, instead of the type chip and moon
  label.
- **Sector table (UX.24).** Order: star systems, then non-stellar
  phenomena (nebulae, belts, comets, anomalies), then rogue planets; a
  detail row is `<tr><td colspan=...>` across every column.
- **Rogue planet octants (UX.25).** The notes number the octants 1 to 8
  by the signs of (dx, dy, dz) from the sector center: 1 (+,+,+),
  2 (-,+,+), 3 (-,-,+), 4 (+,-,+), 5 (+,+,-), 6 (-,+,-), 7 (-,-,-),
  8 (+,-,-). Check it against the octant the site already shows before
  changing anything.
- **Phenomena filters (UX.33).** Filter by kind, then by class: rogue
  planets by class (and terrestrial or gas giant, mass range, octant),
  comets by orbit type (short period, long period, hyperbolic or
  interstellar, sungrazer), and the other kinds by class. Correction:
  the notes' single `phenomena` table with `category` and `class`
  columns doesn't exist; each kind has its own table
  (`rogue_planets.planet_class`, and so on), so the filter is per kind in
  `queryDb`.
- **Comet links (UX.29).** Every comet resolves to a class and links to
  its reference page (period, eccentricity, perihelion). Parabolic
  comets have no period class today.

## Build threads

Each thread is briefed with its exact item IDs and takes no others.

1. Icons and actions: UX.28, UX.26, UX.31, UX.27, UX.25.
2. Sector and system page fixes: UX.24, UX.29, UX.32, then UX.30 (after
   thread 3's ladders).
3. Units: UX.23, then UX.22 by quantity.
4. Phenomena filters: UX.33.
5. API remote generation: API.4 onward, in the order above.
6. Work queue: PERF.19, PERF.20, ADM.15.
7. Backfill: GEN.41, then PERF.18 with GEN.42 and GEN.43.
8. Independent: UX.2, UX.3, ADM.14.
9. Reproducible galaxies: API.16, ADM.17, GEN.59, then API.17 (after
   thread 5's API.12 and API.13).

## Open questions for Boss

- UX.30: does the system overview stay prose?
- UX.32: does "Habitable" stay as a chip?
- UX.3, API.7, PERF.20: the open questions in their TODO.md entries.
