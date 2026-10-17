# Version 8.5: Courses that bend around gravity wells

Rebuilt on 2026-10-16 as a version map (Boss, 1:07 PM EDT: one version per headline, phases inside a version, blockers one phase earlier).
[docs/TODO.md](../TODO.md) is the master file: it holds each item's full text, and its "Plan: versions" section indexes every version.
This file gives the version's goal, its phases and each phase's build threads with every item's prerequisites in order. A phase holds independent threads: an item's blockers sit in the same thread (earlier in the list) or in an earlier phase. Where the two disagree, TODO.md wins; when an item ships, it leaves TODO.md and its row here is deleted in the same PR (the progress ledger keeps its version). Research notes and judgment calls are in [notes.md](notes.md); the headlines and the ideas not yet filed are in [release-roadmap.md](release-roadmap.md).

Version: 8.5

## Goal

Courses that bend around gravity wells, with waypoints, editable and saved courses (NAV.59): saved courses first, then the maps and the keep-out steering, then waypoints and editing.

## Phase 1: Saved courses

The saved course record and its list, open and rename.

### Saved courses

| ID | Item | Needs | Note |
|---|---|---|---|
| NAV.17 | A saved course record with both forms |  |  |
| NAV.18 | Saved courses list: open, rename, delete and the Direct/Route switch, over the saved-list store | NAV.17, UX.129 | Sibling of bookmarks.js. |
| NAV.4 | Save a course | NAV.17, NAV.18 | Per browser now; NAV.19 moves it into accounts in phase 3+. |
| NAV.19 | Saved courses in the account (after USR.7) | USR.7, NAV.18 |  |
| NAV.21 | Fit the view to the whole course |  | Same fit as MAP.53, exempt from any zoom limit (MAP.58 was superseded by MAP.125). |
| NAV.23 | Open a saved course on the map | NAV.4, NAV.21 |  |
| NAV.39 | Saved courses remember their unknown-space jumps and check them again | NAV.17 | With NAV.17. |

### Halo and intergalactic switch

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.537 | Design the plan configuration: the boundaries of the sector plan, how far scatters go, which sector addresses are valid, and the two switches halo and satellite galaxies (both default off) |  | Boss, 2026-10-17, 4:48 PM EDT (two switches, plan configuration). |
| GEN.538 | The plan configuration is the single source of the plan's boundaries, scatter reach and valid sector addresses: scatter and sector generation read it and never decide them | GEN.537, GEN.301, PERF.121 | Boss, 2026-10-17, 4:48 PM EDT (two switches, plan configuration). |
| GEN.539 | The halo command: its own command that switches the halo on or off in the plan configuration (default off) and sets how far the halo reaches | GEN.538, PERF.121 | Boss, 2026-10-17, 4:48 PM EDT (two switches, plan configuration). |
| GEN.540 | The satellite galaxies command: its own command that switches the satellite galaxies on or off in the plan configuration (default off), separate from the halo | GEN.538, PERF.121 | Boss, 2026-10-17, 4:48 PM EDT (two switches, plan configuration). |
| MAP.248 | The maps and pages follow the plan configuration: with the halo or the satellite galaxies off, nothing is drawn, picked or offered beyond the configured boundaries | GEN.538, GEN.539, GEN.540, MAP.243 | Boss, 2026-10-17, 4:48 PM EDT (two switches, plan configuration). |
| UX.137 | The Generate page, the settings and the command line set the halo and the satellite galaxies switches, each through its own command, and say what they are | GEN.539, GEN.540 | Boss, 2026-10-17, 4:48 PM EDT (two switches, plan configuration). |
| TEST.226 | Tests for the plan configuration and its two switches: the interior is identical with each combination, off generates nothing beyond the boundaries, addresses outside them are refused, the 8.1 sphere is the all-on case | GEN.538, GEN.539, GEN.540 | Boss, 2026-10-17, 4:48 PM EDT (two switches, plan configuration). |

## Phase 2: Courses on the maps and around wells

Showing courses on the maps and bending them around keep-out spheres.

### Courses on the maps

| ID | Item | Needs | Note |
|---|---|---|---|
| NAV.22 | Courses inside a sector and a system |  |  |
| NAV.5 | Show a course on the Galaxy Map | NAV.21, NAV.22, NAV.23 | Courses stay drawn until cleared (NAV.49). Parent; most of it exists (MAP.27). |
| NAV.36 | Unknown-space jumps drawn red and glowing |  | Route strip, navmap.py and the Galaxy Map course; legend; reduced motion; both themes. |

### Steering clear of gravity wells

| ID | Item | Needs | Note |
|---|---|---|---|
| NAV.54 | One keep-out table for every kind of object: fields, supermassive holes and quasars, moons, stars in their own system, nebulae and remnants |  | Research: NAV.24 (built) follow-ups. |
| NAV.62 | Gravity of the stars and objects a cell leaves out: work out how a course and the orbital update account for the unresolved mass, derived from the density model | NAV.54 | Boss 5:01 PM EDT: gravity of left-out mass. |
| NAV.25 | Find the obstacles along a path | NAV.54 | Corridor query from NAV.10; sectors along the line from NAV.38. |
| NAV.56 | Census of overlapping keep-out spheres in a generated galaxy | NAV.54 | Research: decides the overlap policy. |
| NAV.26 | Bend the path around keep-out spheres | NAV.25, NAV.56, NAV.62, NAV.54 |  |
| NAV.27 | Moving bodies inside a system | NAV.26, GEN.289 |  |
| NAV.28 | Show and save the adjusted course | NAV.26, NAV.4 |  |
| NAV.51 | Courses route around asteroid fields | NAV.25, NAV.26, NAV.54 | Boss 2026-10-09 01:02Z. Parked with the NAV chain (2026-10-09). |
| NAV.6 | Courses that steer clear of gravity wells | NAV.25, NAV.26, NAV.27, NAV.28, NAV.51 | Parent; closes with its subitems. |

## Phase 3: Waypoints and the headline

Waypoints, editing a course and the umbrella.

### Waypoints and editing

| ID | Item | Needs | Note |
|---|---|---|---|
| NAV.49 | Waypoints: pick objects in Star select mode and plot a course through them, kept on the map until cleared | MAP.122, NAV.17 | Merges "Plotted courses should appear on the galactic map and stay until cleared". |
| NAV.58 | Edit a course: add, remove, reorder and drag its waypoints, then save the edit | NAV.49, NAV.4 | 8.5 headline (Boss 2026-10-10 21:05Z). |
| NAV.59 | Headline of 8.5: courses that bend around gravity wells, with waypoints, editable and saved courses | NAV.6, NAV.49, NAV.58, NAV.4, NAV.23 | 8.5 headline (Boss 2026-10-10 21:05Z). Parent. |
