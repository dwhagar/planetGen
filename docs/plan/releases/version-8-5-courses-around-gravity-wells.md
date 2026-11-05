# Version 8.5: Courses that bend around gravity wells

Rebuilt on 2026-10-16 as a version map (Boss, 1:07 PM EDT: one version per headline, phases inside a version, blockers one phase earlier).
[docs/TODO.md](../../TODO.md) is the master file: it holds each item's full text, and its "Plan: versions" section indexes every version.
This file gives the version's goal, its phases and each phase's build threads with every item's prerequisites in order. A phase holds independent threads: an item's blockers sit in the same thread (earlier in the list) or in an earlier phase. Where the two disagree, TODO.md wins; when an item ships, it leaves TODO.md and its row here is deleted in the same PR (the progress ledger keeps its version). Research notes and judgment calls are in [notes.md](../notes.md); the headlines and the ideas not yet filed are in [release-roadmap.md](../release-roadmap.md).

Version: 8.5

Status: paused until 8.1 is complete (Boss, 7:44 PM EDT 2026-10-22: the five cores work 8.1 only; nothing in this version is started or owned until then, except the gates already running for UX.140, UX.128 and NAV.17).

## Goal

Courses that bend around gravity wells, with waypoints, editable and saved courses (NAV.59): saved courses first, then the maps and the keep-out steering, then waypoints and editing.

## Phase 1: Saved courses

The saved course record and its list, open and rename.

### Saved courses

| ID | Item | Needs | Note |
|---|---|---|---|
| NAV.4 | Save a course |  | Per browser now; NAV.19 moves it into accounts in phase 3+. |
| NAV.19 | Saved courses in the account (after USR.7) | USR.7 |  |
| DB.111 | Storage for vessels and stations: the vessel table (Object ID, owner, name, kind, speed table entry, state, anchor, course id, departure and arrival time, finished legs) and the owner and free-point columns on facilities if it lacks them | NAV.4 |  |

### Halo and intergalactic switch

| ID | Item | Needs | Note |
|---|---|---|---|
| GEN.537 | Design the plan configuration: the boundaries of the sector plan, how far scatters go, which sector addresses are valid, and the two switches halo and satellite galaxies (both default off) |  | Boss, 2026-10-17, 4:48 PM EDT (two switches, plan configuration). |
| GEN.538 | The plan configuration is the single source of the plan's boundaries, scatter reach and valid sector addresses: scatter and sector generation read it and never decide them | GEN.537, GEN.301 | Boss, 2026-10-17, 4:48 PM EDT (two switches, plan configuration). |
| GEN.539 | The halo command: its own command that switches the halo on or off in the plan configuration (default off) and sets how far the halo reaches | GEN.538 | Boss, 2026-10-17, 4:48 PM EDT (two switches, plan configuration). |
| GEN.540 | The satellite galaxies command: its own command that switches the satellite galaxies on or off in the plan configuration (default off), separate from the halo | GEN.538 | Boss, 2026-10-17, 4:48 PM EDT (two switches, plan configuration). |
| GEN.620 | The region removal command: its own admin command that deletes the halo or the satellite galaxies from a made galaxy and makes their cells invalid addresses (region argument halo or satellites, needs the switch already off) | GEN.538, GEN.539, GEN.540 |  |
| MAP.248 | The maps and pages follow the plan configuration: with the halo or the satellite galaxies off, nothing is drawn, picked or offered beyond the configured boundaries | GEN.538, GEN.539, GEN.540 | Boss, 2026-10-17, 4:48 PM EDT (two switches, plan configuration). |
| UX.137 | The Generate page, the settings and the command line set the halo and the satellite galaxies switches, each through its own command, and say what they are | GEN.539, GEN.540 | Boss, 2026-10-17, 4:48 PM EDT (two switches, plan configuration). |
| TEST.226 | Tests for the plan configuration and its two switches: the interior is identical with each combination, off generates nothing beyond the boundaries, addresses outside them are refused, the 8.1 sphere is the all-on case | GEN.538, GEN.539, GEN.540 | Boss, 2026-10-17, 4:48 PM EDT (two switches, plan configuration). |

## Phase 2: Courses on the maps and around wells

Showing courses on the maps and bending them around keep-out spheres.

### Courses on the maps

| ID | Item | Needs | Note |
|---|---|---|---|
| NAV.70 | Schematic route map part 6 (optional): a distance-keeping layout (layout L2) for the neighbourhood, if the ribbon is judged too plain |  |  |
| NAV.77 | vessel_position(vessel, t): the position of a vessel at any galaxy time from its stored course, departure time and speed, agreeing with the course planner's travel times to the second | DB.111, NAV.4, GEN.105 |  |
| MAP.389 | Borders on the maps: hard lines on the Sector and System Maps where two owners meet, a soft haze on the Galaxy Map, and a reserved hatch style for contested systems, drawn from the stored owners and cached by the territory epoch | POP.28 |  |

### Steering clear of gravity wells

| ID | Item | Needs | Note |
|---|---|---|---|
| NAV.56 | Census of overlapping keep-out spheres in a generated galaxy |  | Research: decides the overlap policy. |
| NAV.6 | Courses that steer clear of gravity wells |  | Parent; closes with its subitems. |

### Gravity layers on the Sector and System Maps (Boss, 5:14 PM EST on 2026-11-05)

| ID | Item | Needs | Note |
|---|---|---|---|
| MAP.168 | The Sector Map gravity layer: coloured zones for pull, well depth and tidal strength | API.24 |  |
| MAP.169 | The System Map gravity layer: the orbital plane as a heat map with Lagrange points and Hill spheres | MAP.168 |  |
| UX.94 | Gravity map wording: the legend, the mode names and the Map help text | MAP.168 |  |

### Console look: icons and the maps' chrome

| ID | Item | Needs | Note |
|---|---|---|---|
| UX.269 | Console look on the icon sprite and the maps' chrome (Galaxy, Sector and System Map controls), by tokens only | UX.268 |  |
| UX.274 | Streamline the View menu further: fewer entries, grouped by task, the same on the Galaxy, Sector and System Maps |  |  |

### Edit everything foundations: aspects and the recipe schema (spread back from 8.6)

| ID | Item | Needs | Note |
|---|---|---|---|
| ADM.62 | Aspect catalogue for systems, stars, planets, moons and belts: every stored field and every generation input, with how they depend on each other |  |  |
| ADM.63 | Aspect catalogue for phenomena: black holes, neutron stars, nebulae, remnants, rogue planets, comets, asteroid fields and quasars | ADM.62 |  |
| ADM.64 | Keep-or-regenerate resolution: from the admin's marks to the exact set of aspects redrawn | ADM.62 |  |
| API.28 | Recipe schema and runner: a recipe fixes, ranges or randomises any aspect, and a library turns it into objects | ADM.62, ADM.64 | Shared with 8.9 (API.18, API.19) |

## Phase 3: Waypoints and the headline

Waypoints, editing a course and the umbrella.

Release gate (standing, Boss, 2026-10-21, 5:10 PM EDT): this version is not stamped until SEC.33 (the full security code audit) and TEST.261 (the full bug hunt audit) each have a dated report for it and no finding they mark blocking is open. Both are standing items: they stay open and are re-run for every major or minor version change (required whenever the major or minor version changes (8.1, 8.2 and on to 9.0 and 10.0), never for build numbers).

### Waypoints and editing

| ID | Item | Needs | Note |
|---|---|---|---|
| NAV.59 | Headline of 8.5: courses that bend around gravity wells, with waypoints, editable and saved courses | NAV.6, NAV.4 | 8.5 headline (Boss 2026-10-10 21:05Z). Parent. |

## Phase 4: Trade routes

Spaceport trade routes on the navigation graph, the global list and the map tools (Boss, 9:38 PM EDT 2026-10-24).

Release gate (standing, Boss, 2026-10-21, 5:10 PM EDT): this version is not stamped until SEC.33 (the full security code audit) and TEST.261 (the full bug hunt audit) each have a dated report for it and no finding they mark blocking is open. Both are standing items: they stay open and are re-run for every major or minor version change (required whenever the major or minor version changes (8.1, 8.2 and on to 9.0 and 10.0), never for build numbers).

### Trade routes

| ID | Item | Needs | Note |
|---|---|---|---|
| POP.32 | Trade partners: which spaceports trade with which, by the gravity model on population and distance, with a cap per port and hubs | POP.21, POP.30, DB.105 | Boss, 9:38 PM EDT 2026-10-24: trade routes and world populations, done by 8.9, dynamic in 9.1. |
| POP.33 | Generate the trade routes: a hop-based run on the navigation graph between each pair of spaceports, saved in the one global route list | POP.32 | Boss, 9:38 PM EDT 2026-10-24: trade routes and world populations, done by 8.9, dynamic in 9.1. |
| POP.34 | The global trade route list: a page and an API every account can read, with filters by polity, route length, hop count and spaceport | POP.33 | Boss, 9:38 PM EDT 2026-10-24: trade routes and world populations, done by 8.9, dynamic in 9.1. |
| NAV.71 | Schematic route map part 7: many routes as lines through hubs (trade routes), with the octilinear look (layout L3) | POP.33, POP.34, DB.105 |  |
| MAP.309 | Star map filter for trade routes: draw the routes as paths on the Galaxy Map, filtered by polity, length or hop count | POP.34 | Boss, 9:38 PM EDT 2026-10-24: trade routes and world populations, done by 8.9, dynamic in 9.1. |
| MAP.310 | Routes away from a selected object: select a system, spaceport or world on the star map and show only the routes leaving it | MAP.309 | Boss, 9:38 PM EDT 2026-10-24: trade routes and world populations, done by 8.9, dynamic in 9.1. |
| NAV.63 | Find the nearest trade route: by spaceport, and by where someone could put themselves in the path of a route | POP.33, MAP.309 | Boss, 9:38 PM EDT 2026-10-24: trade routes and world populations, done by 8.9, dynamic in 9.1. |
| TEST.382 | First live test pass of the world population and the trade routes on a generated galaxy, before the 8.5 stamp | POP.33, MAP.310, NAV.63 | Boss, 9:38 PM EDT 2026-10-24: trade routes and world populations, done by 8.9, dynamic in 9.1. |
| NAV.64 | Find the nearest resource location from a position or along a course, by resource kind and minimum score | MAP.312, NAV.63 | Boss, 10:29 PM EDT 2026-10-24: resource locations (valuable belts and fields as points of interest), foundations 8.2 to 9.0. |
