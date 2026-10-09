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
| GEN.33 | One class per PR, each with its tests | GEN.89 | Moved to phase 2 under GEN.90, after the habitability score. One class per PR (R and S first). Built on fixed physics so new classes aren't tuned to wrong masses, moons or zones. |
| GEN.28 | Seven new planet classes in the letter gaps (R, S, U, W, X, Y, Z) | GEN.33 | Under GEN.90; class Z is Boss's Earth-size world that never had life (2026-10-03). Class S landed with GEN.38 (PR #415); the other six classes here. Closes with GEN.33's PRs; PLANET_CLASSES in program_constants.py. |
| GEN.145 | Class S atmosphere rule: S keeps air unless the shoreline ratio is over 30 | GEN.91 | Research follow-up to class S (built). |
| GEN.27 | Class P (glaciated world) only in the habitable zone, and fitting there | GEN.33 | Under GEN.90. Same reconcile/zone code as phase 0's physics fixes. |
| GEN.91 | Classes like S and V in the hot and cold zones | GEN.33 |  |
| GEN.146 | Teff-dependent habitable zone from the Kopparapu table | GEN.91 | Research: in GEN.91's dependency chain. |
| GEN.92 | Life and its highest stage follow the habitability score | GEN.89, GEN.28 |  |
| GEN.29 | Sweep every planet class for sense once the new ones are in (bug) | GEN.28, GEN.27, GEN.91, GEN.92 | Under GEN.90: the refactor is the sweep. Bug, but by definition a sweep after the new classes; it can't go earlier. Includes rocky rogues of 10-16 Earth masses (up to 17,600 km) that get S as nearest fit (PR #415). |
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
| NAV.36 | Unknown-space jumps drawn red and glowing | UX.35 | Route strip, navmap.py and the Galaxy Map course; legend; reduced motion; both themes. |
| NAV.39 | Saved courses remember their unknown-space jumps and check them again | NAV.17 | With NAV.17. |
| NAV.49 | Waypoints: pick objects in Star select mode and plot a course through them, kept on the map until cleared | MAP.122, NAV.17 | Merges "Plotted courses should appear on the galactic map and stay until cleared". |

### Pages

| ID | Item | Needs | Note |
|---|---|---|---|
| UX.32 | Planet rows show the class only, without the type and moon labels |  | Same rows as UX.30; one thread. |
| UX.30 | Planet information without the Markdown render | UX.22, UX.32, UX.23 | Uses the unit ladders; shows the composition rows DB.2 now reads (PR #347). |
| DOC.4 | Correct the stale statements the research found in docs, docstrings and comments |  | Docs and comments only; no behaviour change. |
| UX.43 | A visual design built like a pilot's starmap and navigation console | UX.42 |  |
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
| API.10 | Reservations: claimed sectors and id blocks per run | API.9 | Reserved id blocks build on DB.3's id-block fix (PR #347). |
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
| OPS.34 | Windows Redis in WSL: fix the keep-alive advice and add a Start-RedisInWsl remedy |  | Research: corrects OPS.21 and OPS.27 (built). |
| OPS.31 | Lint every example plist, XML and service file in CI |  | Research: found with the macOS plist bug. |

### Picker

| ID | Item | Needs | Note |
|---|---|---|---|
| NAV.45 | "What's within N pc" from the Galaxy Map and Sector Map |  | Boss 04:39Z. Map action opening NAV.44; drawn sphere optional. |

### Bookmarks

| ID | Item | Needs | Note |
|---|---|---|---|
| UX.48 | Charted regions as bookmarks that frame and outline the region | UX.47 |  |
| UX.45 | Bookmark management | UX.46, UX.47, UX.48 | Parent; per browser until accounts (USR.7) move them. |

### Orbital updates

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.115 | The galaxy's own gravity: a smooth disk, bulge and halo potential |  | Model from Boss's "Computational Astrodynamics.md" (2026-10-07): bulge, disk and halo potential. |
| GEN.144 | One distance for the Sun from the galactic centre across the constants, the density model and the design docs |  | Research: three Sun distances today. |
| GEN.109 | N-body influence from every object inside the largest nearby Hill sphere plus the galactic gradient, with a Hill-radius warning | GEN.115 |  |
| GEN.140 | Orbital math guards the edge-case table adds (GEN.108 built) |  | Research follow-up to GEN.108 (built). |
| GEN.139 | Orbit-update thresholds: per-object epoch, path-length rule and what the 0.01 mpc applies to (GEN.106 built) |  | Research follow-up to GEN.106 (built). |
| ADM.36 | Change an object's trajectory vector | GEN.109 |  |
| GEN.110 | Rogue planet collisions: asteroid fields, merged giants and new stars | GEN.109, GEN.142 |  |
| GEN.143 | A collision_events table, an admin report and a test that runs the whole collision path | GEN.110 | Research: makes the collision code testable. |
| GEN.142 | Peculiar velocity for rogue planets and asteroid fields |  | Research: prerequisite of GEN.110. |
| GEN.105 | Orbital updates | GEN.115, GEN.109, GEN.110 | Parent of the orbital update work. |
| GEN.141 | Faster Kepler solver (Mikkola or Markley) with brentq as fallback |  | Research: optional speed-up. |

### Nebula planets

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.95 | Nebula conditions applied when planets and surfaces are generated | GEN.89 |  |
| GEN.152 | Nebula cloud field is 10 to 40 times too full; lower it to the observed filling (GEN.47 rate check) |  | Research: reopens GEN.47 as a rate check. |
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
| GEN.113 | Analyze the anomaly docs: which anomalies to add and how |  | Boss: "probably Phase 2 and 3 just to analyze those". |
| GEN.154 | Show the Einstein radius on compact-object pages |  | Research: Tier 1 of GEN.113. |
| GEN.153 | Magnetar subtype of neutron star, and an age-dependent pulsar fraction |  | Research: Tier 1 of GEN.113. |

### Maps

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.121 | Every map shows and steps to its neighbouring regions, on one map engine |  | Folds MAP.127, the Galaxy Map's neighbouring blocks and slabs (Boss 2026-10-07 16:27Z: one map engine for all maps). |
| MAP.146 | Zoom drill-down centred on the clicked point, not on fixed wedges, blocks and slabs | ADM.29, MAP.122 | Boss 2026-10-09 22:24Z. Replaces the fixed ladder in galaxy/drill.py and galaxyprisms.js. |
| MAP.147 | The Galaxy Map wire format: measure what the browser downloads and compare smaller options |  | Boss 2026-10-09 22:41Z. Starts with the investigation (Research Lane 3). Decide with MAP.146's tile keys. |

### Recipes

| ID | Item | Needs | Note |
|---|---|---|---|
| API.18 | Generate by recipe: JSON for sectors, systems, planets, moons and phenomena | GEN.96, API.9 | Validation through the Pydantic models; 400 for nonsense, 422 for validation failures. |

### View

| ID | Item | Needs | Note |
|---|---|---|---|

### Features from the GitHub issues

| ID | Item | Needs | Note |
|---|---|---|---|
| PERF.33 | Progress bars and ETAs from measured performance | PERF.32 | GitHub issue [#661](https://github.com/dwhagar/planetGen/issues/661) (the use half). |
| MAP.139 | The Galaxy View uses its spare space: an info box with a Details link, and menu items |  | GitHub issues [#758](https://github.com/dwhagar/planetGen/issues/758) and [#715](https://github.com/dwhagar/planetGen/issues/715), one layout change. |
| MAP.140 | Double-click on a selected object goes there and opens its information |  | GitHub issue [#714](https://github.com/dwhagar/planetGen/issues/714). |
| MAP.141 | Context around the selection: faint neighbours, and the sectors above and below |  | GitHub issues [#718](https://github.com/dwhagar/planetGen/issues/718) and [#716](https://github.com/dwhagar/planetGen/issues/716), one context view. The #716 bug label was overruled. |
| MAP.142 | Nebulae have fuzzy, fading boundaries |  | GitHub issue [#713](https://github.com/dwhagar/planetGen/issues/713). |
| MAP.143 | Color sectors by their number of habitable locations | GEN.89 | GitHub issue [#717](https://github.com/dwhagar/planetGen/issues/717). |
| GEN.129 | Multi-star systems of up to seven stars | GEN.128 | GitHub issue [#777](https://github.com/dwhagar/planetGen/issues/777). Build after the GEN.128 design. |
| GEN.130 | Exotic star systems: a black hole, neutron star or similar at the center | GEN.128 | GitHub issue [#778](https://github.com/dwhagar/planetGen/issues/778). Build after the GEN.128 design. |
| GEN.165 | A comet's time-0 scene position disagrees with its stored position in a rare random system (bug) |  | Bugfix lane. |
| DB.15 | A migration progress bar with the time remaining | PERF.32 | GitHub issue [#727](https://github.com/dwhagar/planetGen/issues/727). |
| DB.18 | Migration helpers for slow DDL: online indexes, instant columns and batched updates | DB.15 | Research: slow DDL for the big tables. |
| ADM.43 | A full configuration page under Admin | ADM.42 | GitHub issue [#515](https://github.com/dwhagar/planetGen/issues/515). |
| ADM.44 | Web, Open Graph and SEO settings | ADM.42, ADM.43 | GitHub issue [#743](https://github.com/dwhagar/planetGen/issues/743). |
| GEN.134 | Tune the star populations to the observed star-formation profile by galactic radius |  | From GEN.133's four unbuilt proposals (star-types-by-galactic-radius.md); Boss to confirm which he wants. |

## Open questions for Boss

- NAV.24: A keep-out radius for every kind of object, see its entry in TODO.md.
- UX.32: Planet rows show the class only, without the type and moon labels, see its entry in TODO.md.
- UX.30: Planet information without the Markdown render, see its entry in TODO.md.
