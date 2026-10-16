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
| NAV.18 | Save, list, open, rename and delete, per browser | NAV.17 | Sibling of bookmarks.js. |
| NAV.4 | Save a course | NAV.17, NAV.18 | Per browser now; NAV.19 moves it into accounts in phase 3+. |
| NAV.19 | Saved courses in the account (after USR.7) | USR.7, NAV.18 |  |
| NAV.21 | Fit the view to the whole course |  | Same fit as MAP.53, exempt from any zoom limit (MAP.58 was superseded by MAP.125). |
| NAV.23 | Open a saved course on the map | NAV.4, NAV.21 |  |
| NAV.39 | Saved courses remember their unknown-space jumps and check them again | NAV.17 | With NAV.17. |

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
| NAV.25 | Find the obstacles along a path |  | Corridor query from NAV.10; sectors along the line from NAV.38. |
| NAV.26 | Bend the path around keep-out spheres | NAV.25 |  |
| NAV.27 | Moving bodies inside a system | NAV.26 |  |
| NAV.28 | Show and save the adjusted course | NAV.26, NAV.4 |  |
| NAV.51 | Courses route around asteroid fields | NAV.25, NAV.26 | Boss 2026-10-09 01:02Z. Parked with the NAV chain (2026-10-09). |
| NAV.6 | Courses that steer clear of gravity wells | NAV.25, NAV.26, NAV.27, NAV.28, NAV.51 | Parent; closes with its subitems. |
| NAV.54 | Keep-out radii for asteroid fields, supermassive holes, moons and nebulae (NAV.24 built) | NAV.51 | Research: NAV.24 (built) follow-ups. |
| NAV.55 | A tuning block for the keep-out knobs |  | Research: keep-out knobs. |
| NAV.56 | Census of overlapping keep-out spheres in a generated galaxy | NAV.51 | Research: decides the overlap policy. |

## Phase 3: Waypoints and the headline

Waypoints, editing a course and the umbrella.

### Waypoints and editing

| ID | Item | Needs | Note |
|---|---|---|---|
| NAV.49 | Waypoints: pick objects in Star select mode and plot a course through them, kept on the map until cleared | MAP.122, NAV.17 | Merges "Plotted courses should appear on the galactic map and stay until cleared". |
| NAV.58 | Edit a course: add, remove, reorder and drag its waypoints, then save the edit | NAV.49, NAV.4 | 8.5 headline (Boss 2026-10-10 21:05Z). |
| NAV.59 | Headline of 8.5: courses that bend around gravity wells, with waypoints, editable and saved courses | NAV.6, NAV.49, NAV.58, NAV.4, NAV.23 | 8.5 headline (Boss 2026-10-10 21:05Z). Parent. |
