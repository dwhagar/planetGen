# Version 9.3: Ships and stations

Rebuilt on 2026-10-16 as a version map (Boss, 1:07 PM EDT: one version per headline, phases inside a version, blockers one phase earlier).
[docs/TODO.md](../TODO.md) is the master file: it holds each item's full text, and its "Plan: versions" section indexes every version.
This file gives the version's goal, its phases and each phase's build threads with every item's prerequisites in order. A phase holds independent threads: an item's blockers sit in the same thread (earlier in the list) or in an earlier phase. Where the two disagree, TODO.md wins; when an item ships, it leaves TODO.md and its row here is deleted in the same PR (the progress ledger keeps its version). Research notes and judgment calls are in [notes.md](notes.md); the headlines and the ideas not yet filed are in [release-roadmap.md](release-roadmap.md).

Version: 9.3

## Goal

Users place their own ships and stations and follow them across the map (USR.10, POP.11). Comes after politics and territories (Boss, 2026-10-16).

## Phase 1: Ships and stations

Users add their own objects and the scatter of vessels and stations for the races.

Release gate (standing, Boss, 2026-10-21, 5:10 PM EDT): this version is not stamped until SEC.33 (the full security code audit) and TEST.261 (the full bug hunt audit) each have a dated report for it and no finding they mark blocking is open. Both are standing items: they stay open and are re-run for every point release.

### Ships and stations

| ID | Item | Needs | Note |
|---|---|---|---|
| DB.63 | Table for ships and stations, user-owned or generated (migration only) |  | Split from USR.10. Core 1 owns migrations. |
| USR.29 | Vessel and station objects: add, move along a course and list, with ownership | DB.63, NAV.4, GEN.105 | Split from USR.10. |
| POP.11 | Scatter space vessels and stations for the races in their territories | DB.63, USR.29 |  |
| USR.30 | Add a ship or station from a page and see it and its track on the maps | USR.29 | Split from USR.10. |
| USR.10 | Users add their own objects (ships and stations) and track their movement on the map | DB.63, USR.29, USR.30 |  |
| MAP.221 | Station and vessel icons on the Galaxy Map | MAP.219, USR.29, POP.11 | Split from MAP.132. |
| MAP.132 | Overlay markers for black holes, nebulae and habitable worlds | MAP.219, MAP.220, MAP.221 | Later (Boss 2026-10-08 01:59Z colors). |
