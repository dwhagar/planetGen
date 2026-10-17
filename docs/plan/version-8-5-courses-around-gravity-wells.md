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
| GEN.537 | Design the halo and intergalactic switch: what it covers, where it is stored, what off means and what turns on later |  | Boss, 2026-10-17, 4:31 PM EDT (halo switch, 8.5). |
| GEN.538 | The halo and intergalactic switch: a stored galaxy setting, default off, that decides whether the halo and intergalactic commands are in the step list | GEN.537, GEN.301, PERF.121 | Boss, 2026-10-17, 4:31 PM EDT (halo switch, 8.5). |
| MAP.248 | The maps and pages follow the halo and intergalactic switch: nothing is drawn, picked or offered beyond the outline when it is off | GEN.538, MAP.243 | Boss, 2026-10-17, 4:31 PM EDT (halo switch, 8.5). |
| UX.137 | The Generate page, the settings and the command line set the halo and intergalactic switch, and say what it is | GEN.538 | Boss, 2026-10-17, 4:31 PM EDT (halo switch, 8.5). |
| TEST.226 | Tests for the halo and intergalactic switch: the interior is identical with it on and off, off generates nothing beyond the outline, on gives the 8.1 sphere | GEN.538 | Boss, 2026-10-17, 4:31 PM EDT (halo switch, 8.5). |

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
