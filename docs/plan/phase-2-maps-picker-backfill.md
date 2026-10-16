# Phase 2: Classes, orbits, courses and recipes

Rebuilt on 2026-10-07 from Boss's lists of 2026-10-03 and 2026-10-07
and the dependency tree (every open item, its prerequisites and the
files it shares). [docs/TODO.md](../TODO.md) is the master file: it holds each
item's full text, and its "Plan: phases" section indexes every phase.
This file gives the phase's goal, its build threads with each item's
prerequisites in order, and its open questions. Where the two disagree,
TODO.md wins; when an item ships, it leaves TODO.md and its row here is
deleted in the same PR. Research notes, the files several items share and
the judgment calls behind the placement are in [notes.md](notes.md).

## Goal

The planet class refactor around the habitability index (with GEN.33, GEN.28, GEN.27 and GEN.29), nebula conditions on planets, n-body orbital updates and rogue collisions, editable trajectories, courses and waypoints, generate-by-recipe in the API, the backfill density pass, daily maintenance, the pilot-style visual design, light-travel positions, the asteroid-field and anomaly plans, and the API pieces remote generation needs first.

## Threads

Each thread is briefed with its exact item IDs and takes no others. Items
run top to bottom inside a thread; "Needs" lists what must land first
(from this phase or an earlier one).

### Classes

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.33 | One class per PR, each with its tests |  | Moved to phase 2 under GEN.90, after the habitability score. One class per PR (R and S first). Built on fixed physics so new classes aren't tuned to wrong masses, moons or zones. |
| GEN.28 | Seven new planet classes in the letter gaps (R, S, U, W, X, Y, Z) | GEN.33 | Under GEN.90; class Z is Boss's Earth-size world that never had life (2026-10-03). Class S landed with GEN.38 (PR #415); the other six classes here. Closes with GEN.33's PRs; PLANET_CLASSES in program_constants.py. |
| GEN.29 | Sweep every planet class for sense once the new ones are in (bug) | GEN.28, GEN.27, GEN.91, GEN.92 | Under GEN.90: the refactor is the sweep. Bug, but by definition a sweep after the new classes; it can't go earlier. Includes rocky rogues of 10-16 Earth masses (up to 17,600 km) that get S as nearest fit (PR #415). |
| GEN.145 | Class S atmosphere rule: S keeps air unless the shoreline ratio is over 30 | GEN.91 | Research follow-up to class S (built). |
| GEN.27 | Class P (glaciated world) only in the habitable zone, and fitting there | GEN.33 | Under GEN.90. Same reconcile/zone code as phase 0's physics fixes. |
| GEN.91 | Classes like S and V in the hot and cold zones | GEN.33 |  |
| GEN.146 | Teff-dependent habitable zone from the Kopparapu table | GEN.91 | Research: in GEN.91's dependency chain. |
| GEN.92 | Life and its highest stage follow the habitability score | GEN.28 |  |
| GEN.90 | Refactor the planet classes around the habitability index | GEN.33, GEN.28, GEN.27, GEN.29, GEN.91, GEN.92 | Parent; takes GEN.33, GEN.28, GEN.27 and GEN.29 as its subitems. |

### Galaxy Map

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.75 | The mini map as a second engine view |  | Locked second camera on MAP.64. |
| MAP.59 | Make it plain that a zoomed-in slab is a slab, not a wedge | MAP.75 |  |
| MAP.132 | Overlay markers for black holes, nebulae and habitable worlds |  | Later (Boss 2026-10-08 01:59Z colors). |

### Courses

| ID | Item | Needs | Note |
|---|---|---|---|
| NAV.21 | Fit the view to the whole course |  | Same fit as MAP.53, exempt from any zoom limit (MAP.58 was superseded by MAP.125). |
| NAV.17 | A saved course record with both forms |  |  |
| NAV.18 | Save, list, open, rename and delete, per browser | NAV.17 | Sibling of bookmarks.js. |
| NAV.4 | Save a course | NAV.17, NAV.18 | Per browser now; NAV.19 moves it into accounts in phase 3+. |
| NAV.36 | Unknown-space jumps drawn red and glowing |  | Route strip, navmap.py and the Galaxy Map course; legend; reduced motion; both themes. |
| NAV.39 | Saved courses remember their unknown-space jumps and check them again | NAV.17 | With NAV.17. |
| NAV.49 | Waypoints: pick objects in Star select mode and plot a course through them, kept on the map until cleared | MAP.122, NAV.17 | Merges "Plotted courses should appear on the galactic map and stay until cleared". |
| NAV.22 | Courses inside a sector and a system |  |  |
| NAV.23 | Open a saved course on the map | NAV.4, NAV.21 |  |
| NAV.5 | Show a course on the Galaxy Map | NAV.21, NAV.22, NAV.23 | Courses stay drawn until cleared (NAV.49). Parent; most of it exists (MAP.27). |
| NAV.25 | Find the obstacles along a path |  | Corridor query from NAV.10; sectors along the line from NAV.38. |
| NAV.26 | Bend the path around keep-out spheres | NAV.25 |  |
| NAV.56 | Census of overlapping keep-out spheres in a generated galaxy | NAV.51 | Research: decides the overlap policy. |
| NAV.55 | A tuning block for the keep-out knobs |  | Research: keep-out knobs. |
| NAV.27 | Moving bodies inside a system | NAV.26 |  |
| NAV.28 | Show and save the adjusted course | NAV.26, NAV.4 |  |
| NAV.51 | Courses route around asteroid fields | NAV.25, NAV.26 | Boss 2026-10-09 01:02Z. Parked with the NAV chain (2026-10-09). |
| NAV.54 | Keep-out radii for asteroid fields, supermassive holes, moons and nebulae (NAV.24 built) | NAV.51 | Research: NAV.24 (built) follow-ups. |
| NAV.6 | Courses that steer clear of gravity wells | NAV.25, NAV.26, NAV.27, NAV.28, NAV.51 | Parent; closes with its subitems. |
| NAV.59 | Headline of 8.5: courses that bend around gravity wells, with waypoints, editable and saved courses | NAV.6, NAV.49, NAV.58, NAV.4, NAV.23 | 8.5 headline (Boss 2026-10-10 21:05Z). Parent. |
| NAV.58 | Edit a course: add, remove, reorder and drag its waypoints, then save the edit | NAV.49, NAV.4 | 8.5 headline (Boss 2026-10-10 21:05Z). |

### Pages

| ID | Item | Needs | Note |
|---|---|---|---|
| UX.32 | Planet rows show the class only, without the type and moon labels |  | Same rows as UX.105 (absorbed UX.30); one thread. |
| UX.113 | The sector page shows the profile of any address in the reach, and heavy cells and far neighbourhoods warn before they run | GEN.235, API.27 | Whole-sphere build train, Boss 2026-10-15. Step 8. Interface. |
| UX.115 | Interface update: every page and control works with the spherical grid, 96-bit IDs, the new phenomena and the region queries (stage after the sphere and generation stages) | GEN.235, GEN.237, GEN.240, NAV.61, API.27, UX.113 | Interface; queued after sphere stages (Boss, 12:52 AM EDT) |
| UX.116 | Sliders for every scatter and fill limit: one limit_slider with a number box, in three new segments under the scatter-limits row of the Generate page | GEN.281, GEN.283 | Core 3; sliders for the scatter and fill limits |
| UX.114 | Scatter-row pages list the stored members and the map shows a system mass |  | Multi-star systems at scatter time, Boss 2026-10-15. Step 8. Interface. |
| DOC.4 | Correct the stale statements the research found in docs, docstrings and comments |  | Docs and comments only; no behaviour change. |
| DOC.48 | Docs: the reach, the halo, the profile rules and the Wolf-Rayet and nebula-contents rules in the design notes | GEN.234, GEN.235, GEN.237, GEN.241, GEN.243 | Whole-sphere build train, Boss 2026-10-15. Step 14. Quality. |
| DOC.50 | Docs: the schema, the multi-star note and the catalogue say where the members live | GEN.262, GEN.264, GEN.266, TEST.164 | Multi-star systems at scatter time, Boss 2026-10-15. Step 10. Quality. |
| DOC.51 | Docs: the cell key in the coordinate system and schema docs, the region note moves to built, the API states the cap and the ladder | GEN.276, NAV.60 | Spherical region query, stage C. Step 8. Quality. |
| DOC.49 | Docs: point the rate, anomaly, multi-star and nebula notes at the phenomena catalogue | GEN.251, GEN.130 | Phenomena catalogue, Boss 2026-10-15. Step 9. Quality. |
| UX.43 | A visual design built like a pilot's starmap and navigation console |  |  |
| UX.82 | Theme checks after PR #800: SVG currentColor, two Shoelace contrast failures, alpha in --bg-subtle |  | Research follow-up to UX.76 (built). |

### Backfill

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.42 | A pass that drops sectors from a region by probability |  | Same function as PERF.18 (backfill_bright_stars_around); one thread. |
| PERF.35 | An interval or chunk ledger for untouched sectors once block-first backfill lands | GEN.42 | Research: ledger size at full galaxy. |
| PERF.18 | Run the GEN.30 bright-star backfill in parallel on the work queue |  | Backfill blocks become RQ jobs. Same stars as the one-process backfill for one seed needs GEN.39. |
| PERF.48 | Low priority: a numeric-only INSERT formatter or C driver for bright_stars and phenomenon_scatter |  | Generation performance study. |
| PERF.36 | Memory and request guard: never list more than about 50,000 candidate cells, and refuse huge enumerations in a web request |  | Research: protects the Generate page estimates. |
| GEN.40 | Weed out sectors by star density before the bright-star backfill | GEN.42 | Parent; closes with its subitems. |

### Queue

| ID | Item | Needs | Note |
|---|---|---|---|
| PERF.20 | Short-term caching through the work queue and API (needs planning) |  | Plans on RQ results and the new caches. Plan with Boss. |

### API

| ID | Item | Needs | Note |
|---|---|---|---|
| API.5 | API version and compatibility checking | API.4 |  |
| API.10 | Reservations: claimed sectors and id blocks per run |  | Reserved id blocks build on DB.3's id-block fix (PR #347). |
| API.11 | Staging tables | API.10 | Galaxy schema migration (staging); after NAV.10 in the writer queue. |
| API.12 | The download: seed, skeleton and name state | API.5 | Downloads the naming key, not name registries. Downloads the seed (what a seed means is GEN.39) and the name state (rules from GEN.46, done in PR #370). |

### 3D system

| ID | Item | Needs | Note |
|---|---|---|---|

### Reproducible galaxies

| ID | Item | Needs | Note |
|---|---|---|---|
| OPS.15 | Each update says whether it changes generated output |  | Needs the fingerprint, so phase 2. |

### Daily maintenance

| ID | Item | Needs | Note |
|---|---|---|---|
| OPS.16 | A daily maintenance script for Linux, macOS and Windows |  | scripts/maintenance.sh and .ps1: positional update, delta merge, rotation; lock; optional OPS.15 check. |
| OPS.30 | A lock helper for the maintenance run |  | Research: used by OPS.16 (ADM.20's admin merge was dropped). |
| OPS.17 | Install and update set up the daily maintenance schedule | OPS.16 | Same scripts as OPS.7/OPS.8/OPS.13 (install/update, deploy-common), after them. |

### Picker

| ID | Item | Needs | Note |
|---|---|---|---|
| NAV.45 | "What's within N pc" from the Galaxy Map and Sector Map |  | Boss 04:39Z. Map action opening NAV.44; drawn sphere optional. |
| NAV.61 | Nearest-first by a radius ladder, and a dense-core guard with an early-stop count | NAV.60 | Spherical region query, stage C. Step 4. Core. |
| NAV.60 | What is nearby reads by the cell key: one statement for every kind, including objects with no sector and unbuilt scatter rows | GEN.276 | Spherical region query, stage C. Step 3. Core. |
| GEN.276 | planetgen/galaxy/cellkey.py: the cell key function and the sphere-to-key-runs cover, numpy only |  | Spherical region query, stage C. Step 2. Core. |

### Bookmarks

| ID | Item | Needs | Note |
|---|---|---|---|
| UX.48 | Charted regions as bookmarks that frame and outline the region | UX.47 |  |
| UX.92 | A Bookmark button for planets, moons and belts | UX.45 |  |
| UX.45 | Bookmark management | UX.46, UX.47, UX.48 | Parent; per browser until accounts (USR.7) move them. |

### Orbital updates

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.115 | The galaxy's own gravity: a smooth disk, bulge and halo potential |  | Model from Boss's "Computational Astrodynamics.md" (2026-10-07): bulge, disk and halo potential. |
| GEN.243 | Wolf-Rayet stars and ring nebulae: a W state for the last 0.3 Myr of very massive stars, and nebula class X | GEN.235, GEN.158, GEN.241 | Whole-sphere build train, Boss 2026-10-15. Step 13. Core. |
| GEN.242 | A pre-main-sequence state for young stars, so protostars can exist as real stars | GEN.235 | Whole-sphere build train, Boss 2026-10-15. Step 11. Core. |
| GEN.241 | Every nebula gets the stars it needs as real systems made with the nebula, inside its bounds | GEN.235, GEN.236, GEN.242, GEN.251 | Whole-sphere build train, Boss 2026-10-15. Step 12. Core. Folds in MAP.182. |
| API.27 | A read-only profile lookup for any address in the reach, so an uncharted sector can be previewed | GEN.235 | Whole-sphere build train, Boss 2026-10-15. Step 8. Interface. |
| GEN.240 | The three scatters fill the whole sphere with one algorithm: band groups, the phenomena pass on the object-first sampler, and a component pass | GEN.234, GEN.235, GEN.237 | Whole-sphere build train, Boss 2026-10-15. Step 6. Core. |
| GEN.239 | High-velocity clouds and the Magellanic Stream as placed gas clouds (later) | GEN.237 | Whole-sphere build train, Boss 2026-10-15. Step 10. Core, later. |
| GEN.259 | Intermediate-mass black hole share gives about 100,000 in the Galaxy against a handful of candidates (bug) |  | Phenomena catalogue, Boss 2026-10-15. Question for Boss. Core. |
| GEN.258 | Protostars embedded in giant molecular clouds | GEN.241, GEN.242 | Phenomena catalogue, Boss 2026-10-15. Step 8, later. Core. |
| GEN.257 | A Fermi-bubble overlay on the Galaxy Map |  | Phenomena catalogue, Boss 2026-10-15. Step 7d, later. Interface. |
| GEN.256 | A radiation hazard score derived from distance to hot stars, Wolf-Rayet stars, remnants and magnetars | GEN.243, GEN.153 | Phenomena catalogue, Boss 2026-10-15. Step 7c, later. Core. |
| GEN.255 | Drawn jets: Herbig-Haro flows, microquasar jets and pulsar jets along their spin axes | GEN.242 | Phenomena catalogue, Boss 2026-10-15. Step 7b, later. Core and Interface. |
| GEN.254 | Superbubble field beside the nebula field (the Local Bubble cavity) | GEN.235 | Phenomena catalogue, Boss 2026-10-15. Step 7a, later. Core. |
| GEN.253 | White-dwarf pollution flag and cool white dwarf presets | GEN.130 | Phenomena catalogue, Boss 2026-10-15. Step 6, later. Core. |
| GEN.252 | Brown dwarfs as stars: an L, T and Y state for companions and primaries, with the brown-dwarf desert | GEN.251 | Phenomena catalogue, Boss 2026-10-15. Step 5, later. Core. |
| GEN.249 | Be stars, carbon stars and Mira variables as star attributes (later) | GEN.235 | Whole-sphere build train, Boss 2026-10-15. Later group, section 5.12. Core, later. |
| GEN.248 | Tidal streams as placed components of the halo (later) | GEN.237, GEN.234 | Whole-sphere build train, Boss 2026-10-15. Later group, section 5.12. Core, later. |
| GEN.247 | Accreting white dwarf binaries: novae, cataclysmic and symbiotic variables, Type Ia progenitors (later) | GEN.129, GEN.130 | Whole-sphere build train, Boss 2026-10-15. Later group, section 5.12. Core, later. |
| GEN.246 | Luminous blue variables and hypergiants (later) | GEN.243 | Whole-sphere build train, Boss 2026-10-15. Later group, section 5.12. Core, later. |
| GEN.245 | Debris disks around A and F stars (later) | GEN.242 | Whole-sphere build train, Boss 2026-10-15. Later group, section 5.12. Core, later. |
| GEN.244 | Young open clusters and OB associations as components (later) | GEN.237, GEN.241 | Whole-sphere build train, Boss 2026-10-15. Later group, section 5.12. Core, later. |
| GEN.238 | Wanderers: hypervelocity stars reach out to the grid, and kicked neutron stars and black holes are an ejecta population | GEN.235 | Whole-sphere build train, Boss 2026-10-15. Step 6. Core. |
| GEN.237 | Satellite and cluster components: one table for the Magellanic Clouds, Sagittarius, dwarf galaxies and globular clusters | GEN.235, GEN.159, GEN.160 | Whole-sphere build train, Boss 2026-10-15. Step 5. Core. |
| GEN.236 | Kind factors in the profile: the rate per star of each rare kind times a factor from the cell, and the planet cull by metallicity | GEN.235, GEN.158 | Whole-sphere build train, Boss 2026-10-15. Steps 2 and 4. Core. |
| GEN.235 | One rule set for every cell: a sector profile computed from the address, and any asked sector in the reach is generated from it and saved | GEN.234 | Whole-sphere build train, Boss 2026-10-15. Step 2. Core. |
| GEN.234 | A measured stellar halo replaces the flat density floor, with the halo, cluster and satellite populations |  | Whole-sphere build train, Boss 2026-10-15. Step 3. Core. |
| GEN.251 | Sector fill stages on their own random streams, and a system pass so every object that can have a system rolls for one | GEN.235 | Phenomena catalogue, Boss 2026-10-15. Step 2. Core. |
| GEN.267 | A MAX_STARS setting and a --stars N directive, with a prevalence knob for five or more stars |  | Multi-star systems at scatter time, Boss 2026-10-15. Step 11. Core. |
| GEN.275 | Substellar and circumbinary types: brown-dwarf and planetary-mass companions, circumbinary disks and planets |  | Multi-star catalogue, Boss 2026-10-15 6:08 PM EDT. Step 19. Core. |
| GEN.274 | Cluster-hosted compact binaries and runaway bookkeeping: `runaway_origin` and the compact companion left by a supernova |  | Multi-star catalogue, Boss 2026-10-15 6:08 PM EDT. Step 18. Core. |
| GEN.273 | Young groups and trapezia: the `group` node, offset and velocity columns, a decay age and the remnant draw |  | Multi-star catalogue, Boss 2026-10-15 6:08 PM EDT. Step 17. Core. |
| GEN.272 | Accreting compact binaries: X-ray binaries, accretion luminosity, microquasars, spider pulsars and compact triples | GEN.269, GEN.266 | Multi-star catalogue, Boss 2026-10-15 6:08 PM EDT. Step 16. Core. |
| GEN.271 | White-dwarf binary types: close double white dwarfs, post-common-envelope pairs, cataclysmic variables, AM CVn, symbiotic and barium stars | GEN.269 | Multi-star catalogue, Boss 2026-10-15 6:08 PM EDT. Step 15. Core. |
| GEN.270 | Close-binary types: contact, Algol, spectroscopic and eclipsing flags, tertiary shares, merger products and Be stars | GEN.269 | Multi-star catalogue, Boss 2026-10-15 6:08 PM EDT. Step 14. Core. |
| GEN.269 | Interaction states in the star model: a Roche classifier, `interaction` on every pair, member state and flags, rejuvenated age and the new kinds |  | Multi-star catalogue, Boss 2026-10-15 6:08 PM EDT. Step 13. Core. |
| GEN.266 | Compact-binary orbit rules: a black hole or neutron star with a star, tight compact pairs and tertiaries |  | Multi-star systems at scatter time, Boss 2026-10-15. Step 7. Core. |
| GEN.265 | Nightly motion of a system: the barycentre moves with the system mass, members are computed from their stored phase on read |  | Multi-star systems at scatter time, Boss 2026-10-15. Step 6. Core. |
| PERF.103 | A per-night cost ceiling on the orbital update: warn with the backlog size and move the rest the next night |  | Orbital update, Boss 2026-10-15 6:26 PM EDT. Small. Core. |
| PERF.104 | Cache warm-up for the whole site: build every list, map, sector, system and search index cache after an update, then rebuild only what a change invalidates |  |  |
| GEN.264 | The fill builds a system from its stored members; the system pass keeps only planets, moons and belts | GEN.262, GEN.251 | Multi-star systems at scatter time, Boss 2026-10-15. Step 5. Core. |
| GEN.263 | Promote every bound compact pair to the scatter, whatever its mass | GEN.262 | Multi-star systems at scatter time, Boss 2026-10-15. Step 4. Core. |
| GEN.262 | Companions in the scatter passes: the mass, luminosity and phenomenon passes draw members, with a rate correction so pairs are not double counted |  | Multi-star systems at scatter time, Boss 2026-10-15. Step 3. Core. |
| GEN.144 | One distance for the Sun from the galactic centre across the constants, the density model and the design docs |  | Research: three Sun distances today. |
| GEN.109 | N-body influence from every object inside the largest nearby Hill sphere plus the galactic gradient, with a Hill-radius warning | GEN.115 |  |
| GEN.291 | A bound planet becomes a rogue planet when an encounter unbinds it: the unbind and capture test and the event log | GEN.289, GEN.290 |  |
| GEN.290 | An encounter search that cannot miss a pass at large steps (100 years per second) | GEN.109 |  |
| GEN.289 | Planets and moons stored as a position and velocity relative to their star, so perturbations have somewhere to land | GEN.109 |  |
| GEN.140 | Orbital math guards the edge-case table adds (GEN.108 built) |  | Research follow-up to GEN.108 (built). |
| GEN.139 | Orbit-update thresholds: per-object epoch, path-length rule and what the 0.01 mpc applies to (GEN.106 built) |  | Research follow-up to GEN.106 (built). |
| ADM.36 | Change an object's trajectory vector | GEN.109 |  |
| GEN.110 | Rogue planet collisions: asteroid fields, merged giants and new stars | GEN.109, GEN.142 |  |
| GEN.143 | A collision_events table, an admin report and a test that runs the whole collision path | GEN.110 | Research: makes the collision code testable. |
| GEN.197 | Object IDs on ejection, merger and split events | GEN.143 |  |
| GEN.142 | Peculiar velocity for rogue planets and asteroid fields |  | Research: prerequisite of GEN.110. |
| GEN.105 | Orbital updates | GEN.115, GEN.109, GEN.110 | Parent of the orbital update work. |
| GEN.293 | Long-step rules for hypervelocity stars that leave the sector grid and for parabolic comets |  |  |
| GEN.292 | Run the orbital update for a stated span in one go from the command line (elapsed years) |  |  |
| GEN.198 | The gravity field evaluator: pull, potential and tidal gradient at any point from the point masses and the galaxy potential | GEN.115, GEN.109 | Gravity map foundation (Boss 2026-10-10). |
| TEST.131 | Tests for the gravity field evaluator: two bodies, Lagrange points, softening and the rotation curve | GEN.198 |  |
| TEST.159 | Whole-sphere tests: required stars per nebula kind, Wolf-Rayet counts, tiling at the reach, profile and empty-sector determinism, and the pre-placed share per population | GEN.235, GEN.236, GEN.241, GEN.243 | Whole-sphere build train, Boss 2026-10-15. Step 15. Quality. |
| TEST.162 | Multiplicity audit: companion states, row bookkeeping, stability, determinism, no planets in a scatter, fill reuse, motion and the census anchors | GEN.262, GEN.264 | Multi-star systems at scatter time, Boss 2026-10-15. Step 9. Quality. |
| TEST.164 | Catalogue type audit: type counts, region gates, the N tail, interaction consistency and group decay | GEN.269, GEN.270, GEN.271, GEN.272, GEN.273, GEN.274, GEN.275 | Multi-star catalogue, Boss 2026-10-15 6:08 PM EDT. Step 20. Quality. |
| TEST.166 | The cell key follows an object across a cube boundary, and the point columns and the key round-trip bit for bit |  | Spherical region query, stage C. Step 7. Quality. |
| TEST.165 | Re-measure the region query on a real generated galaxy, warm and cold | GEN.276 | Spherical region query, stage C. Step 6. Quality. |
| TEST.168 | Control sweep: a browser test drives every control on every page against a spherical galaxy with 96-bit IDs | UX.115 | Quality; queued after UX.115 |
| TEST.169 | Admin page browser tests: every admin page and task works end to end after the streamline | ADM.56 | Quality; queued admin chain |
| TEST.161 | Fill-stage and system-pass tests: stream independence, same seed same system, and the probability table | GEN.251 | Phenomena catalogue, Boss 2026-10-15. Step 2 tests. Quality. |
| TEST.160 | Representation audit: seeded sectors in each region class show every kind at its expected frequency, with systems and required contents checked | GEN.251 | Phenomena catalogue, Boss 2026-10-15. Step 4. Quality. |
| PERF.72 | Research the cost of a sector's gravity grid and where to cut between exact and aggregated sources | GEN.198 | Research. |
| GEN.199 | Per-zone gravity grid for a sector: stored in a cache and rebuilt when the orbits update | GEN.198, PERF.72 |  |
| MAP.183 | Optional, default not built for the first release: near-view P and Q nebulae from the seeded field below about 250 pc |  | Phase 2. Optional, default not built. |
| MAP.181 | Dust colour for the dark nebula family on both themes |  | Phase 2. Foundations lane 2. |
| MAP.180 | Picking, hover and fly-to for nebula regions and cover | MAP.176 | Phase 2. Foundations lane 2. |
| MAP.179 | Scattered nebulae as markers from level 8 | GEN.202, MAP.173 | Phase 2. Foundations lane 2. |
| MAP.178 | Field nebula clouds drawn one by one from the nebula table | MAP.173 | Phase 2. Foundations lane 2. |
| MAP.177 | Color by "Nebula cover" on the Galaxy Map | MAP.174 | Phase 2. Foundations lane 2. |
| MAP.176 | Nebula regions on the Galaxy Map: grouped soft sprites in the dominant family colour | MAP.173 | Phase 2. Foundations lane 2. |
| MAP.207 | The Galaxy Map draws satellites, clusters and wanderers as points beyond the disk outline | GEN.237, GEN.238 | Whole-sphere build train, Boss 2026-10-15. Step 9. Interface. |
| MAP.211 | Hypervelocity stars show no accessible information: show their type, and a button to the phenomenon or the star system (bug) |  |  |
| MAP.210 | Nothing on the map can be selected, even in sector view (bug) |  |  |
| MAP.209 | Opening a sector from the Galaxy Map flashes one view, reloads into another, and leaves the details panel empty (bug) |  |  |
| MAP.208 | Stars and all phenomena on the three maps draw as points of light until near sector level, phenomena as bright as a star of similar mass |  |  |
| POP.13 | A species editor: create and edit species and their environment limits |  |  |
| GEN.141 | Faster Kepler solver (Mikkola or Markley) with brentq as fallback |  | Research: optional speed-up. |

### Nebula planets

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.95 | Nebula conditions applied when planets and surfaces are generated |  |  |
| GEN.151 | Supernova remnant sizes from the density-dependent Sedov-Taylor law (GEN.10 follow-up) |  | Research follow-up to GEN.10 (built). |
| GEN.93 | Nebula conditions in planet generation | GEN.95 | Parent. |
| GEN.149 | Planetary-nebula central stars: 0.5 to 0.7 Msun, 1e2 to 1e4 Lsun, up to 2e5 K |  | Research: sub-item of GEN.93. |

### Asteroid fields

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.112 | Plan asteroid fields and belts as object systems for rendering |  | Plan only. |

### Anomalies

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.154 | Show the Einstein radius on compact-object pages |  | Research: Tier 1 of GEN.113. |
| GEN.153 | Magnetar subtype of neutron star, and an age-dependent pulsar fraction |  | Research: Tier 1 of GEN.113. |
| GEN.280 | Galactic core density: nuclear star cluster and nuclear disc terms in the density model, replacing the Gaussian core of ADM.49 | GEN.235 | Core 3; step 1 of the galactic core note; needs a Core 1 migration revision |
| GEN.281 | Galactic core cell star ceiling and fill luminosity cap: a dense cell fills only stars at or above the cap, and the unresolved count and mass are computed from the address | GEN.280, GEN.276 | Core 3; step 2 of the galactic core note; shares one Core 1 migration with GEN.283 |
| GEN.282 | Galactic core gas and young terms: the central molecular zone, the nebula rate and young clusters (Arches, Quintuplet, the central parsec) | GEN.237, GEN.244 | Core 3; step 3 of the galactic core note |
| GEN.283 | Scatter limits by place: a core luminosity floor and mass limit, a halo luminosity floor (default 500, lowest 50), and the baseline everywhere else | GEN.235, GEN.158 | Core 3; step 3b of the galactic core note; shares one Core 1 migration with GEN.281 |
| GEN.284 | Planet survival at the galactic core: a metallicity factor and a survival formula instead of the hard cut | GEN.280, GEN.162 | Core 3; step 4 of the galactic core note |
| GEN.285 | Nuclear mass components in the galactic potential for the orbital update | GEN.280 | Core 3; step 5 of the galactic core note |
| GEN.286 | Horizontal-branch and AGB phases in the stellar evolution model |  | Core 3; HB and AGB phases (Boss, 1:12 AM EDT); plan in note 5.3; no migration; interim [Fe/H] defaults until GEN.158 |
| GEN.287 | The stellar evolution model spreads giant-branch stars evenly in log luminosity, which probably overstates the bright giants by about a factor of two against real halo counts (bug) |  | Core 3; folded into GEN.286 unless Core 3 says otherwise |
| GEN.278 | Magnetar scatter class: magnetars drawn as scatter rows of young neutron stars | GEN.250, GEN.235, GEN.153 | Owner Core 3; Core 3 note default stands. |

### Maps

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.121 | Every map shows and steps to its neighbouring regions, on one map engine |  | Folds MAP.127, the Galaxy Map's neighbouring blocks and slabs (Boss 2026-10-07 16:27Z: one map engine for all maps). |
| MAP.156 | View one layer or a range of layers top-down from the galaxy view, as a secondary option |  | Boss 2026-10-09 23:17Z; secondary option, not a main control. |
| MAP.195 | Galaxy Map finer tiles leave out the stars a coarser tile already sent |  |  |
| UX.102 | One logging standard for every output the user or an admin reads: the same timestamps, format, dividers and step markers everywhere |  |  |
| GEN.212 | Galactic move in the orbital update as NumPy batches (3 to 4 times on the whole update) |  |  |
| PERF.86 | Redo the work queue: drop the task and subtask logic and use Redis simply |  |  |
| ADM.52 | A star system editor: plan a system as a 2D diagram or from the command line |  |  |
| GEN.210 | The galactic core is one round disc sector in the middle of the grid, 4 pc across and 4 pc tall |  |  |
| UX.101 | PHI-6 symbols in the planet and moon lists, coloured by rating | UX.100 |  |
| UX.100 | Planet and moon lists show habitability as symbols, not notes |  |  |
| MAP.193 | Show every object's axis of rotation, direction of rotation and axial tilt |  |  |
| MAP.192 | See the objects whose positions changed in the last orbital update |  |  |
| MAP.191 | Galaxy Map default view: no density shading, only nebulae shade the stars behind them |  |  |
| MAP.190 | Several separate regions of created sectors: which one the Galaxy Map opens on | MAP.189 |  |
| GEN.189 | Gamma-ray burst and AGN ozone loss for the ozone_loss_flag (low priority) |  | Foundations lane 2 report 06:20Z. |
| DOC.16 | The API help page for visitors: what the API is, how to get a key, and where the reference lives | DOC.6 | Boss 23:53Z; Phase 2. |
| DOC.15 | A glossary and units help page: coordinates, scales, sector paths, object IDs, time and the in-universe wording | DOC.6 | Boss 23:53Z; Phase 2. |
| DOC.14 | The admin help pages: the queue, the stats page, settings, lockouts and the naming key | DOC.6 | Boss 23:53Z; Phase 2. |
| DOC.13 | The Generate page help pages: layer specs, spans, radial fills, directives, one-off systems and jobs | DOC.6 | Boss 23:53Z; Phase 2. |
| DOC.12 | The account help pages: signing in, two-factor, the account page, API keys and bookmarks | DOC.6 | Boss 23:53Z; Phase 2. |
| DOC.11 | The reference browser help pages: species, polities, object classes and phenomena | DOC.6 | Boss 23:53Z; Phase 2. |
| DOC.10 | The search and navigation help pages: search, nearby, the nav page and routes | DOC.6 | Boss 23:53Z; Phase 2. |
| DOC.9 | The star system help pages: the system list, a system page, the system map and the planets, moons and belts shown | DOC.6 | Boss 23:53Z; Phase 2. |
| DOC.8 | The sector help pages: the sector list, a sector page, the sector map and the sector scene | DOC.6 | Boss 23:53Z; Phase 2. |
| DOC.7 | The Galaxy Map help page: layers, zoom, fly-through, Color by, select modes, bookmarks and the locate box | DOC.6 | Boss 23:53Z; Phase 2. |
| DOC.6 | A static help section in the web interface: page template, index, per-page help links and a coverage test |  | Boss 23:53Z; the framework the per-feature help pages use. |
| MAP.160 | Quantise the Galaxy Map GPU buffers (deferred) |  | Deferred; MAP.147 recommendation step 3. |

### Recipes

| ID | Item | Needs | Note |
|---|---|---|---|
| API.18 | Generate by recipe: JSON for sectors, systems, planets, moons and phenomena |  | Validation through the Pydantic models; 400 for nonsense, 422 for validation failures. |

### View

| ID | Item | Needs | Note |
|---|---|---|---|

### Features from the GitHub issues

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.139 | The Galaxy View uses its spare space: an info box with a Details link, and menu items |  | GitHub issues [#758](https://github.com/dwhagar/planetGen/issues/758) and [#715](https://github.com/dwhagar/planetGen/issues/715), one layout change. |
| MAP.140 | Double-click on a selected object goes there and opens its information |  | GitHub issue [#714](https://github.com/dwhagar/planetGen/issues/714). |
| MAP.141 | Context around the selection: faint neighbours, and the sectors above and below |  | GitHub issues [#718](https://github.com/dwhagar/planetGen/issues/718) and [#716](https://github.com/dwhagar/planetGen/issues/716), one context view. The #716 bug label was overruled. |
| MAP.142 | Nebulae have fuzzy, fading boundaries |  | GitHub issue [#713](https://github.com/dwhagar/planetGen/issues/713). |
| MAP.143 | Color sectors by their number of habitable locations |  | GitHub issue [#717](https://github.com/dwhagar/planetGen/issues/717). |
| GEN.129 | Multi-star systems of up to twelve stars |  | GitHub issue [#777](https://github.com/dwhagar/planetGen/issues/777). Build after the GEN.128 design. |
| GEN.130 | Compact-remnant systems at sector fill: companions and pulsar planets for neutron stars and black holes | GEN.251, GEN.266 | Phenomena catalogue, Boss 2026-10-15. Step 3, rescoped. Core. |
| DB.18 | Migration helpers for slow DDL: online indexes, instant columns and batched updates |  | Research: slow DDL for the big tables. |
| ADM.43 | A full configuration page under Admin |  | GitHub issue [#515](https://github.com/dwhagar/planetGen/issues/515). |
| ADM.44 | Web, Open Graph and SEO settings | ADM.43 | GitHub issue [#743](https://github.com/dwhagar/planetGen/issues/743). |
| GEN.134 | Tune the star populations to the observed star-formation profile by galactic radius |  | From GEN.133's four unbuilt proposals (star-types-by-galactic-radius.md); Boss to confirm which he wants. |

## Open questions for Boss

- NAV.24: A keep-out radius for every kind of object, see its entry in TODO.md.
- UX.32: Planet rows show the class only, without the type and moon labels, see its entry in TODO.md.
