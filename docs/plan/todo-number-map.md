# TODO number map: old numbers and tree IDs to category IDs

## Next free IDs

A new item, subitem or bug takes its category's next free ID, and the
same PR moves that row up by one. Each category's counter is its next free number minus one (the highest
ID ever issued in it, which is also its item count). The table only
allocates IDs; since 2026-10-10 it no longer feeds the version (BUILD
counts stamped releases, see `changes/README.md`).
`scripts/bump_version.py --check` fails if `docs/TODO.md` uses an ID at
or past a category's next free one, so a stale row is caught before a
release is stamped.

| Category | Next free ID |
|---|---|
| UX | UX.257 |
| MAP | MAP.364 |
| NAV | NAV.65 |
| GEN | GEN.1236 |
| PERF | PERF.231 |
| DB | DB.110 |
| API | API.38 |
| ADM | ADM.138 |
| SEC | SEC.37 |
| TEST | TEST.553 |
| USR | USR.32 |
| OPS | OPS.94 |
| DOC | DOC.101 |
| VIEW | VIEW.19 |
| POP | POP.46 |

## Why this exists

Until 2026-10-01, `docs/TODO.md` numbered its items with running numbers
("1." to "111."). Those numbers are now replaced by permanent category IDs
such as `UX.1` or `MAP.16` (see [TODO.md](../TODO.md) for the scheme):
the category and a plain running count within it, like the database
schema version.

For a short time on 2026-10-01 (PR #189, merged about 06:05Z, to the
flat-ID PR) the IDs were a tree with dots for subitems, such as
`MAP.2.1.1`. Boss then asked for flat IDs ("CATEGORY.NUMBER (sequential,
like the DB version, no points on the category items, just as item coutn
on each, simple)"). Every dotted ID was replaced; "Tree IDs to flat IDs"
below maps them, for the commits, PR titles (#189, #190) and the old
tags that cite them. Everything else in this document uses the flat IDs.

The old numbers live on in places that can't be edited:
[CHANGELOG.md](../../CHANGELOG.md) entries ("TODO items 40 and 41"),
commit messages and pull request titles on GitHub ("(TODO 70, 71)"), and
older versions of the docs and code comments. This document maps every
old number ever used in `docs/TODO.md`, open or finished, to its new ID.

One old number can mean several different items. Before about 21:00 UTC
on 2026-09-30 the file was renumbered many times ("renumber the rest"),
and several branches carried different numberings at the same time.
After that, finished items were deleted without renumbering, but a few
numbers were still reused or collided (see the notes). So an old number
only means something together with the date it was written.

This map was built from the full history of `docs/TODO.md` on `main`
and every branch (126 commits, 2026-09-23 to 2026-10-01, up to release
7.40.0 and PR #160, then updated for the 30 commits that touched it up to
release 7.58.2 and PR #180), from CHANGELOG.md, from commit and PR titles, and
from the docs and code comments that cite numbers.

## How to use it

1. Find the date (UTC) of the thing that cites the number: the commit
   or the PR (GitHub shows both). Every known CHANGELOG citation is
   already resolved under "Where old numbers are cited".
2. Look up the number in the main table and pick the row whose date
   range covers that date.
3. If two rows overlap (parallel branches), pick the one whose title
   matches what the citing text describes. The notes list every known
   overlap.

Date ranges run from the first commit that gave the number that meaning
to the last commit (on any branch) that still had it. A number kept its
meaning until the next commit that changed it, so a citation a little
after the end of a range still belongs to it. A range that ends at
2026-10-01 05:29Z (the 7.58.2 release stamp, the last change before this
update) means the number still had that meaning in TODO.md at 7.58.2.
For an item finished after 7.40.0, the range ends when the PR that
deleted it merged to `main` (its branch had deleted it a little
earlier). All times are UTC.

IDs for open items are the ones in the renumbered TODO.md. Finished items
got IDs too, in the same categories, so a citation of a finished item
still resolves; those IDs are permanent and never reused. Statuses give
the release in CHANGELOG.md and the PR that shipped the item. Every
release note that was pending at 7.40.0 has since been released (7.40.1
to 7.42.0), so no status says "next release" any more.

## The numberings, by date

| From (UTC) | Numbering |
|---|---|
| 2026-09-23 and before | No numbers. "Phase 0" to "Phase 5" in older docs are roadmap phases, not items, and are all complete. |
| 2026-09-24 01:32Z (PR #71) | First numbering, 1 to 15. |
| 2026-09-24 01:57Z to 02:27Z | The PR #73 branch renumbered its copy (1 OOM, 2 skeleton, 3 cache ...) while `main` kept the first numbering and only dropped item 1 at 02:18Z (PR #72). See note 1. |
| 2026-09-24 02:27Z (PR #73) | Renumbered 1 to 10. |
| 2026-09-24 02:53Z (branch; `main` at 05:51Z, PR #81) | Renumbered 1 to 9. Unchanged until 2026-09-30. |
| 2026-09-30 16:44Z (branch; `main` at 18:09Z, PR #106) | The Galaxy Map plan inserted as 3 to 12; 1 to 18. |
| 2026-09-30 16:49Z and 16:51Z (same branch) | A new 8 and a new 14 inserted; 1 to 20. |
| 2026-09-30 18:14Z to 18:39Z (branch only) | Boss's notes inserted; 1 to 37. See note 9. |
| 2026-09-30 18:39Z (`main` at 18:41Z, PR #109) | 1 to 42. Items 1 to 38 keep these numbers from here on. |
| 2026-09-30 19:02Z (`main` at 19:16Z, PR #110) | Security findings 39 to 52, generation bugs 53 to 58, Population 59 to 62. See note 7. |
| 2026-09-30 20:01Z (branch; `main` at 20:04Z, PR #111) | Web pages 59 to 62, Population 63 to 66. |
| 2026-09-30 20:07Z (`main` at 20:23Z, PR #112) | Security findings fixed: 39 Hardening, 40 to 45 bugs, 46 to 49 web pages, 50 to 53 Population. |
| 2026-09-30 20:43Z (`main` at 20:44Z, PR #114) | 50 installers, Population 51 to 54 (final). |
| 2026-09-30 20:48Z (`main` at 21:04Z, PR #115) | Windows jobs written as 54 on its branch, 55 after the merge. |
| 2026-09-30 about 21:00Z onward | No more renumbering: finished items deleted, new items appended (55 to 92). Exceptions: notes 2, 5 and 6. |
| 2026-10-01 04:50Z to 05:27Z (branch; `main` at 05:27Z, PR #175) | New items appended: 93 to 95, the Galaxy Map bugs 96 to 99 (added in merge commit 8767450), 100, 101, 102 and 103. No renumbering. See note 13. |
| 2026-10-01 05:40Z to 05:44Z (branch; `main` at 05:46Z, PR #185) | New items appended: 104 to 108. No renumbering. |
| 2026-10-01 05:50Z to 05:56Z (branch; `main` at 05:56Z, PR #187) | New items appended: 109 to 111. No renumbering. The category IDs replaced the numbers in the docs-refresh PR, after `main` at 05:56Z (PR #187). |

## Main table: old number to new ID

Sorted by old number, then date.

| Old | Date range (UTC) | New ID | Title | Status |
|---|---|---|---|---|
| 1 | 2026-09-24 01:32Z to 02:02Z | OPS.2 | Apache OOM-killed on the production server | done in 5.47.0, PR #72 |
| 1 | 2026-09-24 02:25Z to 2026-09-30 18:09Z | PERF.2 | Cache so pages don't hit the database every request | done in 7.56.0, PR #178 |
| 1 | 2026-09-30 18:14Z to 22:15Z | UX.6 | Every distance in its most meaningful unit | done in 7.10.0, PR #122 |
| 2 | 2026-09-24 01:32Z to 02:18Z | MAP.7 | Remove the large sphere marker for a star in a sector | done in 5.47.1, PR #73 |
| 2 | 2026-09-24 01:57Z to 02:02Z | MAP.6 | Unfilled-sector skeleton draws as a sphere | done in 5.47.0, PR #72 (test added in 5.47.1, PR #73) |
| 2 | 2026-09-24 02:25Z to 2026-09-30 18:09Z | MAP.10 | Draw the Measure distance path around obstacles | done in 7.22.0, PR #135 |
| 2 | 2026-09-30 18:14Z to 22:36Z | UX.7 | Planet list: type chip, habitable-moon chip, belt distances | done in 7.12.0, PR #126 |
| 3 | 2026-09-24 01:32Z to 02:18Z | MAP.8 | System Map orbital paths back, not over bodies | done in 5.47.1, PR #73 |
| 3 | 2026-09-24 01:57Z to 02:02Z | PERF.2 | Cache so pages don't hit the database every request | done in 7.56.0, PR #178 |
| 3 | 2026-09-24 02:25Z to 05:38Z | MAP.11 | Every kind of phenomenon on the Sector Map, clickable | done in 5.51.0, PR #81 |
| 3 | 2026-09-24 02:53Z to 2026-09-30 17:58Z | MAP.5 | Rework the Galaxy Map | replaced on 2026-09-30 by MAP.31 to MAP.42 and MAP.1 |
| 3 | 2026-09-30 16:44Z to 18:09Z | MAP.31 | Spiral arms stand out in the density shading | done in 7.9.0, PR #120 |
| 3 | 2026-09-30 18:14Z to 22:36Z | UX.8 | A less dense top bar (gear menu) | done in 7.11.0, PR #126 |
| 4 | 2026-09-24 01:32Z to 02:18Z | MAP.9 | Star glow renders as an opaque shell | done in 5.47.1, PR #73 |
| 4 | 2026-09-24 01:57Z to 02:02Z | MAP.10 | Draw the Measure distance path around obstacles | done in 7.22.0, PR #135 |
| 4 | 2026-09-24 02:25Z to 05:38Z | MAP.5 | Rework the Galaxy Map | replaced on 2026-09-30 by MAP.31 to MAP.42 and MAP.1 |
| 4 | 2026-09-24 02:53Z to 2026-09-30 17:58Z | API.1 | Create a system inside an existing sector | done in 7.43.0, PR #170 |
| 4 | 2026-09-30 16:44Z to 18:09Z | MAP.32 | Scale readout in sectors, pc and ly | done in 7.9.0, PR #120 |
| 4 | 2026-09-30 18:14Z to 2026-10-01 00:34Z | UX.9 | Tag search: collapsible groups and phenomena | done in 7.22.1, PR #135, and 7.29.0, PR #143 |
| 5 | 2026-09-24 01:32Z to 02:18Z | MAP.6 | Unfilled-sector skeleton draws as a sphere | done in 5.47.0, PR #72 (test added in 5.47.1, PR #73) |
| 5 | 2026-09-24 01:57Z to 02:02Z | MAP.11 | Every kind of phenomenon on the Sector Map, clickable | done in 5.51.0, PR #81 |
| 5 | 2026-09-24 02:25Z to 05:38Z | API.1 | Create a system inside an existing sector | done in 7.43.0, PR #170 |
| 5 | 2026-09-24 02:53Z to 2026-09-30 17:58Z | API.2 | Edit a system's generated content | done in 7.43.0, PR #170 |
| 5 | 2026-09-30 16:44Z to 18:09Z | MAP.33 | Hybrid master-wedge slot rule (schema v35) | done in 7.13.0, PR #129 |
| 5 | 2026-09-30 18:14Z to 23:21Z | GEN.1 | Real-world rates for interstellar objects | done in 7.18.0, PR #136 |
| 6 | 2026-09-24 01:32Z to 02:18Z | PERF.2 | Cache so pages don't hit the database every request | done in 7.56.0, PR #178 |
| 6 | 2026-09-24 01:57Z to 02:02Z | MAP.5 | Rework the Galaxy Map | replaced on 2026-09-30 by MAP.31 to MAP.42 and MAP.1 |
| 6 | 2026-09-24 02:25Z to 05:38Z | API.2 | Edit a system's generated content | done in 7.43.0, PR #170 |
| 6 | 2026-09-24 02:53Z to 2026-09-30 17:58Z | POP.1 | Government ownership of systems | done in 7.49.0, PR #169 (optional since 7.58.2, PR #180) |
| 6 | 2026-09-30 16:44Z to 18:09Z | MAP.34 | Mega-blocks sized from the pixel scale | done in 7.17.0, PR #134 |
| 6 | 2026-09-30 18:14Z to 23:21Z | GEN.2 | Rogue planet mass bins | done in 7.18.0, PR #136 |
| 7 | 2026-09-24 01:32Z to 02:18Z | MAP.10 | Draw the Measure distance path around obstacles | done in 7.22.0, PR #135 |
| 7 | 2026-09-24 01:57Z to 02:02Z | API.1 | Create a system inside an existing sector | done in 7.43.0, PR #170 |
| 7 | 2026-09-24 02:25Z to 05:38Z | POP.1 | Government ownership of systems | done in 7.49.0, PR #169 (optional since 7.58.2, PR #180) |
| 7 | 2026-09-24 02:53Z to 2026-09-30 17:58Z | POP.2 | Names for dominant species on living worlds | done in 7.49.0, PR #169 |
| 7 | 2026-09-30 16:44Z to 18:09Z | MAP.35 | Continuous blocks and a Slice control | done in 7.9.0, PR #120 |
| 7 | 2026-09-30 18:14Z to 22:54Z | GEN.3 | A supermassive black hole in every galaxy | done in 7.15.0, PR #132 |
| 8 | 2026-09-24 01:32Z to 02:18Z | MAP.11 | Every kind of phenomenon on the Sector Map, clickable | done in 5.51.0, PR #81 |
| 8 | 2026-09-24 01:57Z to 02:02Z | API.2 | Edit a system's generated content | done in 7.43.0, PR #170 |
| 8 | 2026-09-24 02:25Z to 05:38Z | POP.2 | Names for dominant species on living worlds | done in 7.49.0, PR #169 |
| 8 | 2026-09-24 02:53Z to 2026-09-30 17:58Z | POP.3 | Database of spacefaring species | done in 7.49.0, PR #169 |
| 8 | 2026-09-30 16:44Z | MAP.38 | Block info on click | done in 7.25.0, PR #142 |
| 8 | 2026-09-30 16:49Z to 18:09Z | MAP.36 | One solid of blocks for filled and unfilled sectors | done in 7.25.0, PR #142 |
| 8 | 2026-09-30 18:14Z to 2026-10-01 05:05Z | PERF.2 | Cache so pages don't hit the database every request | done in 7.56.0, PR #178 |
| 9 | 2026-09-24 01:32Z to 02:18Z | MAP.5 | Rework the Galaxy Map | replaced on 2026-09-30 by MAP.31 to MAP.42 and MAP.1 |
| 9 | 2026-09-24 01:57Z to 02:02Z | POP.1 | Government ownership of systems | done in 7.49.0, PR #169 (optional since 7.58.2, PR #180) |
| 9 | 2026-09-24 02:25Z to 05:38Z | POP.3 | Database of spacefaring species | done in 7.49.0, PR #169 |
| 9 | 2026-09-24 02:53Z to 2026-09-30 17:58Z | POP.4 | Younger and older civilizations | done in 7.49.0, PR #169 |
| 9 | 2026-09-30 16:44Z | MAP.39 | Smooth zooming | done in 7.32.0, PR #147 |
| 9 | 2026-09-30 16:49Z to 18:09Z | MAP.38 | Block info on click | done in 7.25.0, PR #142 |
| 9 | 2026-09-30 18:14Z to 23:48Z | MAP.10 | Draw the Measure distance path around obstacles | done in 7.22.0, PR #135 |
| 10 | 2026-09-24 01:32Z to 02:18Z | API.1 | Create a system inside an existing sector | done in 7.43.0, PR #170 |
| 10 | 2026-09-24 01:57Z to 02:02Z | POP.2 | Names for dominant species on living worlds | done in 7.49.0, PR #169 |
| 10 | 2026-09-24 02:25Z to 05:38Z | POP.4 | Younger and older civilizations | done in 7.49.0, PR #169 |
| 10 | 2026-09-30 16:44Z | MAP.40 | Keep three.js; record why | done, PR #153 (docs only, shipped with 7.36.0) |
| 10 | 2026-09-30 16:49Z to 18:09Z | MAP.39 | Smooth zooming | done in 7.32.0, PR #147 |
| 10 | 2026-09-30 18:14Z to 22:03Z | MAP.31 | Spiral arms stand out in the density shading | done in 7.9.0, PR #120 |
| 11 | 2026-09-24 01:32Z to 02:18Z | API.2 | Edit a system's generated content | done in 7.43.0, PR #170 |
| 11 | 2026-09-24 01:57Z to 02:02Z | POP.3 | Database of spacefaring species | done in 7.49.0, PR #169 |
| 11 | 2026-09-30 16:44Z | MAP.1 | Galaxy Map follow-ups (edge cases) | done in 7.42.1, PR #168 |
| 11 | 2026-09-30 16:49Z to 18:09Z | MAP.40 | Keep three.js; record why | done, PR #153 (docs only, shipped with 7.36.0) |
| 11 | 2026-09-30 18:14Z to 22:03Z | MAP.32 | Scale readout in sectors, pc and ly | done in 7.9.0, PR #120 |
| 12 | 2026-09-24 01:32Z to 02:18Z | POP.1 | Government ownership of systems | done in 7.49.0, PR #169 (optional since 7.58.2, PR #180) |
| 12 | 2026-09-24 01:57Z to 02:02Z | POP.4 | Younger and older civilizations | done in 7.49.0, PR #169 |
| 12 | 2026-09-30 16:44Z | MAP.41 | Remove the server's leftover density sampling | done in 7.9.0, PR #120 |
| 12 | 2026-09-30 16:49Z to 18:09Z | MAP.1 | Galaxy Map follow-ups (edge cases) | done in 7.42.1, PR #168 |
| 12 | 2026-09-30 18:14Z to 22:42Z | MAP.33 | Hybrid master-wedge slot rule (schema v35) | done in 7.13.0, PR #129 |
| 13 | 2026-09-24 01:32Z to 02:18Z | POP.2 | Names for dominant species on living worlds | done in 7.49.0, PR #169 |
| 13 | 2026-09-30 16:44Z | API.1 | Create a system inside an existing sector | done in 7.43.0, PR #170 |
| 13 | 2026-09-30 16:49Z to 18:09Z | MAP.41 | Remove the server's leftover density sampling | done in 7.9.0, PR #120 |
| 13 | 2026-09-30 18:14Z to 23:14Z | MAP.34 | Mega-blocks sized from the pixel scale | done in 7.17.0, PR #134 |
| 14 | 2026-09-24 01:32Z to 02:18Z | POP.3 | Database of spacefaring species | done in 7.49.0, PR #169 |
| 14 | 2026-09-30 16:44Z | API.2 | Edit a system's generated content | done in 7.43.0, PR #170 |
| 14 | 2026-09-30 16:49Z | API.1 | Create a system inside an existing sector | done in 7.43.0, PR #170 |
| 14 | 2026-09-30 16:51Z to 18:09Z | UX.10 | Timestamps in the viewer's own time zone | done in 7.21.0, PR #135 |
| 14 | 2026-09-30 18:14Z to 22:03Z | MAP.35 | Continuous blocks and a Slice control | done in 7.9.0, PR #120 |
| 15 | 2026-09-24 01:32Z to 02:18Z | POP.4 | Younger and older civilizations | done in 7.49.0, PR #169 |
| 15 | 2026-09-30 16:44Z | POP.1 | Government ownership of systems | done in 7.49.0, PR #169 (optional since 7.58.2, PR #180) |
| 15 | 2026-09-30 16:49Z | API.2 | Edit a system's generated content | done in 7.43.0, PR #170 |
| 15 | 2026-09-30 16:51Z to 18:09Z | API.1 | Create a system inside an existing sector | done in 7.43.0, PR #170 |
| 15 | 2026-09-30 18:14Z to 2026-10-01 00:35Z | MAP.36 | One solid of blocks for filled and unfilled sectors | done in 7.25.0, PR #142 |
| 16 | 2026-09-30 16:44Z | POP.2 | Names for dominant species on living worlds | done in 7.49.0, PR #169 |
| 16 | 2026-09-30 16:49Z | POP.1 | Government ownership of systems | done in 7.49.0, PR #169 (optional since 7.58.2, PR #180) |
| 16 | 2026-09-30 16:51Z to 18:09Z | API.2 | Edit a system's generated content | done in 7.43.0, PR #170 |
| 16 | 2026-09-30 18:14Z to 2026-10-01 00:35Z | MAP.38 | Block info on click | done in 7.25.0, PR #142 |
| 17 | 2026-09-30 16:44Z | POP.3 | Database of spacefaring species | done in 7.49.0, PR #169 |
| 17 | 2026-09-30 16:49Z | POP.2 | Names for dominant species on living worlds | done in 7.49.0, PR #169 |
| 17 | 2026-09-30 16:51Z to 18:09Z | POP.1 | Government ownership of systems | done in 7.49.0, PR #169 (optional since 7.58.2, PR #180) |
| 17 | 2026-09-30 18:14Z to 2026-10-01 01:44Z | MAP.39 | Smooth zooming | done in 7.32.0, PR #147 |
| 18 | 2026-09-30 16:44Z | POP.4 | Younger and older civilizations | done in 7.49.0, PR #169 |
| 18 | 2026-09-30 16:49Z | POP.3 | Database of spacefaring species | done in 7.49.0, PR #169 |
| 18 | 2026-09-30 16:51Z to 18:09Z | POP.2 | Names for dominant species on living worlds | done in 7.49.0, PR #169 |
| 18 | 2026-09-30 18:14Z to 2026-10-01 02:31Z | MAP.40 | Keep three.js; record why | done, PR #153 (docs only, shipped with 7.36.0) |
| 19 | 2026-09-30 16:49Z | POP.4 | Younger and older civilizations | done in 7.49.0, PR #169 |
| 19 | 2026-09-30 16:51Z to 18:09Z | POP.3 | Database of spacefaring species | done in 7.49.0, PR #169 |
| 19 | 2026-09-30 18:14Z to 2026-10-01 04:16Z | MAP.1 | Galaxy Map follow-ups (edge cases) | done in 7.42.1, PR #168 |
| 20 | 2026-09-30 16:51Z to 18:09Z | POP.4 | Younger and older civilizations | done in 7.49.0, PR #169 |
| 20 | 2026-09-30 18:14Z to 22:03Z | MAP.41 | Remove the server's leftover density sampling | done in 7.9.0, PR #120 |
| 21 | 2026-09-30 18:14Z to 22:03Z | MAP.42 | Wedge lines from the center | done in 7.9.0, PR #120 |
| 22 | 2026-09-30 18:14Z to 23:48Z | UX.10 | Timestamps in the viewer's own time zone | done in 7.21.0, PR #135 |
| 23 | 2026-09-30 18:14Z to 2026-10-01 00:34Z | ADM.2 | Generate a column and a shell | done in 7.27.0, PR #143 |
| 24 | 2026-09-30 18:14Z to 2026-10-01 01:58Z | ADM.3 | Generate buttons on unfilled sectors (Sector Map, Galaxy Map) | done in 7.27.0 (PR #143) and 7.34.0 (PR #151) |
| 25 | 2026-09-30 18:14Z to 2026-10-01 05:05Z | UX.17 | A view that suits each phenomenon | done in 7.57.0, PR #178 |
| 26 | 2026-09-30 18:14Z to 2026-10-01 04:37Z | UX.18 | Show phenomena's octant and nearest systems | done: storage in 7.33.0 (PR #149), pages in 7.48.0 (PR #167) |
| 27 | 2026-09-30 18:14Z | GEN.6 | Finish the correlative update (galactic motion) | done in 7.37.0, PR #157 |
| 27 | 2026-09-30 18:39Z to 2026-10-01 03:53Z | GEN.10 | Generate nebulae and remnants with their stars, and map them | done in 7.30.0 (PR #144), 7.36.0 (PR #153) and 7.41.1 (PR #148) |
| 28 | 2026-09-30 18:14Z | NAV.1 | Courses in "bearing mark mark" on nested frames | done in 7.14.0, PR #130 (see note 4) |
| 28 | 2026-09-30 18:39Z to 23:32Z | GEN.11 | Class nebulae and remnants A-W | done in 7.19.0, PR #138 |
| 29 | 2026-09-30 18:14Z | NAV.2 | Warp and fold speeds | done in 7.8.0, PR #121 |
| 29 | 2026-09-30 18:39Z to 2026-10-01 00:35Z | GEN.12 | Record what sits inside a nebula | done in 7.24.0, PR #140 |
| 30 | 2026-09-30 18:14Z | DB.1 | Starbases, colonies and outposts in the database | done in 7.35.0, PR #152 |
| 30 | 2026-09-30 18:39Z to 2026-10-01 01:34Z | GEN.13 | Names that follow one standard | done in 7.31.0, PR #145 |
| 31 | 2026-09-30 18:14Z | UX.5 | Place facilities from the web interface | done in 7.47.0, PR #167 |
| 31 | 2026-09-30 18:39Z to 23:32Z | GEN.14 | Class asteroid fields | done in 7.19.0, PR #138 |
| 32 | 2026-09-30 18:14Z | API.1 | Create a system inside an existing sector | done in 7.43.0, PR #170 |
| 32 | 2026-09-30 18:39Z to 2026-10-01 02:35Z | GEN.6 | Finish the correlative update (galactic motion) | done in 7.37.0, PR #157 |
| 33 | 2026-09-30 18:14Z | API.2 | Edit a system's generated content | done in 7.43.0, PR #170 |
| 33 | 2026-09-30 18:39Z to 2026-10-01 02:57Z | NAV.1 | Courses in "bearing mark mark" on nested frames | done in 7.14.0, PR #130 (see note 4) |
| 34 | 2026-09-30 18:14Z | POP.1 | Government ownership of systems | done in 7.49.0, PR #169 (optional since 7.58.2, PR #180) |
| 34 | 2026-09-30 18:39Z to 21:54Z | NAV.2 | Warp and fold speeds | done in 7.8.0, PR #121 |
| 35 | 2026-09-30 18:14Z | POP.2 | Names for dominant species on living worlds | done in 7.49.0, PR #169 |
| 35 | 2026-09-30 18:39Z to 2026-10-01 02:24Z | DB.1 | Starbases, colonies and outposts in the database | done in 7.35.0, PR #152 |
| 36 | 2026-09-30 18:14Z | POP.3 | Database of spacefaring species | done in 7.49.0, PR #169 |
| 36 | 2026-09-30 18:39Z to 2026-10-01 04:37Z | UX.5 | Place facilities from the web interface | done in 7.47.0, PR #167 |
| 37 | 2026-09-30 18:14Z | POP.4 | Younger and older civilizations | done in 7.49.0, PR #169 |
| 37 | 2026-09-30 18:39Z to 2026-10-01 04:24Z | API.1 | Create a system inside an existing sector | done in 7.43.0, PR #170 |
| 38 | 2026-09-30 18:39Z to 2026-10-01 04:24Z | API.2 | Edit a system's generated content | done in 7.43.0, PR #170 |
| 39 | 2026-09-30 18:39Z to 18:41Z | POP.1 | Government ownership of systems | done in 7.49.0, PR #169 (optional since 7.58.2, PR #180) |
| 39 | 2026-09-30 19:02Z to 20:27Z | SEC.3 | Seeded admin/password login claimable | done in 7.5.0, PR #112 |
| 39 | 2026-09-30 20:07Z to 22:53Z | SEC.16 | Hardening (HSTS, lock file, login backoff, input bounds) | done: HSTS in 7.5.0 (PR #112); rest as SEC.17 to SEC.19 |
| 39b | 2026-09-30 21:37Z to 22:27Z | SEC.18 | Upper bounds on admin generation inputs (cited as "39b") | done in 7.10.2, PR #118 (see note 8) |
| 39c | 2026-09-30 21:39Z to 21:54Z | SEC.19 | Hashed lock file for pip dependencies (cited as "39c") | done in 7.7.0, PR #119 (see note 8) |
| 39a | 2026-09-30 22:33Z to 22:51Z | SEC.17 | Per-username login backoff (cited as "39a") | done in 7.14.1, PR #131 (see note 8) |
| 40 | 2026-09-30 18:39Z to 18:41Z | POP.2 | Names for dominant species on living worlds | done in 7.49.0, PR #169 |
| 40 | 2026-09-30 19:02Z to 20:27Z | SEC.4 | Web-user compromise can become root via the installer | done in 7.5.0, PR #112 |
| 40 | 2026-09-30 20:07Z to 23:25Z | GEN.15 | Moons orbit outside their planet's Hill sphere (bug) | done in 7.18.1, PR #137 |
| 41 | 2026-09-30 18:39Z to 18:41Z | POP.3 | Database of spacefaring species | done in 7.49.0, PR #169 |
| 41 | 2026-09-30 19:02Z to 20:27Z | SEC.5 | HTML pages have no rate limit | done in 7.5.0, PR #112 |
| 41 | 2026-09-30 20:07Z to 23:25Z | GEN.16 | Moons can orbit inside their planet (bug) | done in 7.18.1, PR #137 |
| 42 | 2026-09-30 18:39Z to 18:41Z | POP.4 | Younger and older civilizations | done in 7.49.0, PR #169 |
| 42 | 2026-09-30 19:02Z to 20:27Z | SEC.6 | /tmp fallback for jobs and tile cache can be hijacked | done in 7.5.0, PR #112 |
| 42 | 2026-09-30 20:07Z to 22:43Z | GEN.17 | A close binary's planets inside the binary (bug) | done in 7.13.1, PR #128 |
| 43 | 2026-09-30 19:02Z to 20:27Z | SEC.7 | Debug log is mode 0666 | done in 7.5.0, PR #112 |
| 43 | 2026-09-30 20:07Z to 23:09Z | GEN.18 | A planet's Hill sphere overlaps the belt inside it (bug) | done in 7.16.1, PR #133 |
| 44 | 2026-09-30 19:02Z to 20:27Z | SEC.8 | Login time reveals which usernames exist | done in 7.5.0, PR #112 |
| 44 | 2026-09-30 20:07Z to 22:08Z | GEN.19 | A binary's secondary outweighs its primary (bug) | done in 7.9.1, PR #123 |
| 45 | 2026-09-30 19:02Z to 20:27Z | SEC.9 | Changing credentials leaves other sessions logged in | done in 7.5.0, PR #112 |
| 45 | 2026-09-30 20:07Z to 21:53Z | GEN.20 | Sector growth ignores black holes and neutron stars (bug) | done in 7.6.1, PR #116 |
| 46 | 2026-09-30 19:02Z to 20:27Z | SEC.10 | Sector wiki link accepts any scheme | done in 7.5.0, PR #112 |
| 46 | 2026-09-30 20:07Z | POP.1 | Government ownership of systems | done in 7.49.0, PR #169 (optional since 7.58.2, PR #180) |
| 46 | 2026-09-30 20:08Z to 23:48Z | UX.11 | Paginated list of every system | done in 7.20.0, PR #135 |
| 47 | 2026-09-30 19:02Z to 20:27Z | SEC.11 | /api/databases lists the control schema | done in 7.5.0, PR #112 |
| 47 | 2026-09-30 20:07Z | POP.2 | Names for dominant species on living worlds | done in 7.49.0, PR #169 |
| 47 | 2026-09-30 20:08Z to 2026-10-01 00:34Z | MAP.12 | Arc-segment wireframe on the Sector Map | done in 7.29.1, PR #143 |
| 48 | 2026-09-30 19:02Z to 20:27Z | SEC.12 | Public endpoints return raw database errors | done in 7.5.0, PR #112 |
| 48 | 2026-09-30 20:07Z | POP.3 | Database of spacefaring species | done in 7.49.0, PR #169 |
| 48 | 2026-09-30 20:08Z to 22:36Z | UX.12 | System page: one ordered list of everything in orbit | done in 7.12.0, PR #126 |
| 49 | 2026-09-30 19:02Z to 20:27Z | SEC.13 | CSRF token not tied to the login session | done in 7.5.0, PR #112 |
| 49 | 2026-09-30 20:07Z | POP.4 | Younger and older civilizations | done in 7.49.0, PR #169 |
| 49 | 2026-09-30 20:08Z to 2026-10-01 04:37Z | MAP.4 | System Map names never overlap | done in 7.21.1, PR #135 (see note 3) |
| 50 | 2026-09-30 19:02Z to 20:27Z | SEC.14 | New API key's value rides in the flash cookie | done in 7.5.0, PR #112 |
| 50 | 2026-09-30 20:08Z to 20:48Z | POP.1 | Government ownership of systems | done in 7.49.0, PR #169 (optional since 7.58.2, PR #180) |
| 50 | 2026-09-30 20:43Z to 2026-10-01 02:57Z | OPS.3 | PowerShell installers and macOS-safe bash scripts | done in 7.16.0, PR #125 (see note 4) |
| 51 | 2026-09-30 19:02Z to 20:27Z | SEC.15 | Nothing sets config.json's permissions | done in 7.5.0, PR #112 |
| 51 | 2026-09-30 20:08Z to 20:48Z | POP.2 | Names for dominant species on living worlds | done in 7.49.0, PR #169 |
| 51 | 2026-09-30 20:43Z to 2026-10-01 04:37Z | POP.1 | Government ownership of systems | done in 7.49.0, PR #169 (optional since 7.58.2, PR #180) |
| 52 | 2026-09-30 19:02Z to 20:27Z | SEC.16 | Hardening (HSTS, lock file, login backoff, input bounds) | done: HSTS in 7.5.0 (PR #112); rest as SEC.17 to SEC.19 |
| 52 | 2026-09-30 20:08Z to 20:48Z | POP.3 | Database of spacefaring species | done in 7.49.0, PR #169 |
| 52 | 2026-09-30 20:43Z to 2026-10-01 04:37Z | POP.2 | Names for dominant species on living worlds | done in 7.49.0, PR #169 |
| 53 | 2026-09-30 19:02Z to 20:27Z | GEN.15 | Moons orbit outside their planet's Hill sphere (bug) | done in 7.18.1, PR #137 |
| 53 | 2026-09-30 20:08Z to 20:48Z | POP.4 | Younger and older civilizations | done in 7.49.0, PR #169 |
| 53 | 2026-09-30 20:43Z to 2026-10-01 04:37Z | POP.3 | Database of spacefaring species | done in 7.49.0, PR #169 |
| 54 | 2026-09-30 19:02Z to 20:27Z | GEN.16 | Moons can orbit inside their planet (bug) | done in 7.18.1, PR #137 |
| 54 | 2026-09-30 20:43Z to 2026-10-01 04:37Z | POP.4 | Younger and older civilizations | done in 7.49.0, PR #169 |
| 54 | 2026-09-30 20:48Z | OPS.4 | Generate page jobs on native Windows | done in 7.9.2, PR #124 |
| 55 | 2026-09-30 19:02Z to 20:27Z | GEN.17 | A close binary's planets inside the binary (bug) | done in 7.13.1, PR #128 |
| 55 | 2026-09-30 20:48Z to 22:12Z | OPS.4 | Generate page jobs on native Windows | done in 7.9.2, PR #124 |
| 55 | 2026-09-30 23:48Z to 2026-10-01 02:41Z | GEN.21 | Star population model: population ages at sector fill | done in 7.23.0 (PR #141) and 7.38.0 (PR #159) |
| 55 | 2026-10-01 00:07Z to 00:34Z | MAP.13 | Sector Map star dots sized to giants and white dwarfs | done in 7.29.2, PR #143 |
| 56 | 2026-09-30 19:02Z to 20:27Z | GEN.18 | A planet's Hill sphere overlaps the belt inside it (bug) | done in 7.16.1, PR #133 |
| 56 | 2026-09-30 23:48Z to 2026-10-01 02:41Z | GEN.22 | Pre-place bright stars at plan time | done in 7.38.0, PR #159 (uncertain, see note 2) |
| 56 | 2026-10-01 01:19Z to 04:37Z | UX.1 | Class reference pages | done in 7.46.0, PR #167 |
| 57 | 2026-09-30 19:02Z to 20:27Z | GEN.19 | A binary's secondary outweighs its primary (bug) | done in 7.9.1, PR #123 |
| 57 | 2026-10-01 01:15Z to 05:29Z | ADM.5 | Central validate module | done, PR #235 |
| 58 | 2026-09-30 19:02Z to 20:27Z | GEN.20 | Sector growth ignores black holes and neutron stars (bug) | done in 7.6.1, PR #116 |
| 58 | 2026-10-01 01:15Z to 05:29Z | ADM.6 | Override a planet's or moon's class | done, PR #260 |
| 59 | 2026-09-30 19:02Z to 19:17Z | POP.1 | Government ownership of systems | done in 7.49.0, PR #169 (optional since 7.58.2, PR #180) |
| 59 | 2026-09-30 20:01Z to 20:27Z | UX.11 | Paginated list of every system | done in 7.20.0, PR #135 |
| 59 | 2026-10-01 01:15Z to 05:29Z | ADM.7 | Override a star | done, PR #260 |
| 60 | 2026-09-30 19:02Z to 19:17Z | POP.2 | Names for dominant species on living worlds | done in 7.49.0, PR #169 |
| 60 | 2026-09-30 20:01Z to 20:27Z | MAP.12 | Arc-segment wireframe on the Sector Map | done in 7.29.1, PR #143 |
| 60 | 2026-10-01 01:15Z to 05:29Z | ADM.8 | Delete and regenerate buttons, sector down | done, PR #244 |
| 61 | 2026-09-30 19:02Z to 19:17Z | POP.3 | Database of spacefaring species | done in 7.49.0, PR #169 |
| 61 | 2026-09-30 20:01Z to 20:27Z | UX.12 | System page: one ordered list of everything in orbit | done in 7.12.0, PR #126 |
| 61 | 2026-10-01 01:15Z to 05:29Z | SEC.1 | Lock out an IP after failed logins | done, PR #220 |
| 62 | 2026-09-30 19:02Z to 19:17Z | POP.4 | Younger and older civilizations | done in 7.49.0, PR #169 |
| 62 | 2026-09-30 20:01Z to 20:27Z | MAP.4 | System Map names never overlap | done in 7.21.1, PR #135 (see note 3) |
| 62 | 2026-10-01 01:44Z to 05:29Z | UX.2 | Menus sized to what they hold | done, PR #528 |
| 63 | 2026-09-30 20:01Z to 20:27Z | POP.1 | Government ownership of systems | done in 7.49.0, PR #169 (optional since 7.58.2, PR #180) |
| 63 | 2026-10-01 01:44Z to 05:05Z | MAP.3 | A bigger Galaxy Map with controls underneath | done in 7.55.0, PR #178 |
| 64 | 2026-09-30 20:01Z to 20:27Z | POP.2 | Names for dominant species on living worlds | done in 7.49.0, PR #169 |
| 64 | 2026-10-01 02:13Z to 05:29Z | USR.2 | Roles: user, admin, Owner | open |
| 64 | 2026-10-01 02:24Z to 02:26Z | MAP.28 | Nested ladder geometry (drill-down section 3) | done in 7.41.2, PR #160 |
| 65 | 2026-09-30 20:01Z to 20:27Z | POP.3 | Database of spacefaring species | done in 7.49.0, PR #169 |
| 65 | 2026-10-01 02:13Z to 05:29Z | USR.3 | SMTP settings | open |
| 65 | 2026-10-01 02:24Z to 02:26Z | MAP.29 | Stage contents API (drill-down section 7) | done in 7.41.3, PR #160 |
| 66 | 2026-09-30 20:01Z to 20:27Z | POP.4 | Younger and older civilizations | done in 7.49.0, PR #169 |
| 66 | 2026-10-01 02:13Z to 05:29Z | USR.4 | Invite-only sign-up | open |
| 66 | 2026-10-01 02:24Z to 02:26Z | MAP.16 | Drill-down stages | done in 7.44.0, PR #171 (bug MAP.17 fixed, PR #208) |
| 67 | 2026-10-01 02:13Z to 05:29Z | USR.5 | Email loop for passwords | open |
| 67 | 2026-10-01 02:24Z to 02:26Z | MAP.20 | Generate from the sector level | open (map buttons and radius dialog done in 7.53.0, PR #177) |
| 68 | 2026-10-01 02:13Z to 05:29Z | USR.6 | Owner transfer | open |
| 68 | 2026-10-01 02:24Z to 02:26Z | MAP.21 | Sector Map pick mode and Nav links | done in 7.58.0, PR #178 |
| 69 | 2026-10-01 02:13Z to 05:29Z | USR.7 | User-level interface with bookmarks | open |
| 69 | 2026-10-01 02:24Z to 02:26Z | MAP.22 | NAV page picks on the map | done, Bookmarks select on NAV, PR #234 |
| 70 | 2026-10-01 02:24Z to 02:26Z | MAP.23 | Bookmarks | done, per-browser bookmarks (static/bookmarks.js), PR #234 |
| 70 | 2026-10-01 02:27Z to 03:54Z | MAP.28 | Nested ladder geometry (drill-down section 3) | done in 7.41.2, PR #160 |
| 71 | 2026-10-01 02:24Z to 02:26Z | MAP.24 | Address bar | done in 7.50.0, PR #172 |
| 71 | 2026-10-01 02:27Z to 03:54Z | MAP.29 | Stage contents API (drill-down section 7) | done in 7.41.3, PR #160 |
| 72 | 2026-10-01 02:24Z to 02:26Z | MAP.25 | "Show on Galaxy Map" links | open |
| 72 | 2026-10-01 02:27Z to 04:33Z | MAP.16 | Drill-down stages | done in 7.44.0, PR #171 (bug MAP.17 fixed, PR #208) |
| 73 | 2026-10-01 02:24Z to 02:26Z | MAP.27 | NAV course on the Galaxy Map | done in 7.52.0, PR #176 |
| 73 | 2026-10-01 02:27Z to 05:29Z | MAP.20 | Generate from the sector level | done: map buttons and radius dialog in 7.53.0 (PR #177), block and layer generate after 7.58.2 (PR #182, #183) |
| 74 | 2026-10-01 02:27Z to 05:05Z | MAP.21 | Sector Map pick mode and Nav links | done in 7.58.0, PR #178 |
| 75 | 2026-10-01 02:27Z to 05:29Z | MAP.22 | NAV page picks on the map | done, Bookmarks select on NAV, PR #234 |
| 76 | 2026-10-01 02:27Z to 05:29Z | MAP.23 | Bookmarks | done, per-browser bookmarks (static/bookmarks.js), PR #234 |
| 77 | 2026-10-01 02:27Z to 04:41Z | MAP.24 | Address bar | done in 7.50.0, PR #172 |
| 78 | 2026-10-01 02:27Z to 05:29Z | MAP.25 | "Show on Galaxy Map" links | done after 7.58.2, PR #188 (bug MAP.26 fixed, PR #208) |
| 79 | 2026-10-01 02:27Z to 04:51Z | MAP.27 | NAV course on the Galaxy Map | done in 7.52.0, PR #176 |
| 80 | 2026-10-01 02:41Z to 05:29Z | DOC.1 | Number TODO items by category (this renumbering) | done in the docs-refresh PR (version-scheme questions moved to OPS.1) |
| 81 | 2026-10-01 02:51Z to 05:29Z | DOC.2 | Architecture document | done in the docs-refresh PR |
| 82 | 2026-10-01 02:51Z to 05:29Z | DOC.3 | Design documents current, with reasons | done in the docs-refresh PR |
| 83 | 2026-10-01 02:55Z to 05:29Z | VIEW.2 | Starmap seen from a planet | open |
| 84 | 2026-10-01 02:55Z to 05:29Z | VIEW.3 | Render the view as a PNG with constellations | open |
| 85 | 2026-10-01 02:55Z to 05:29Z | VIEW.4 | Constellation names in the name generator | open |
| 86 | 2026-10-01 03:15Z to 05:29Z | PERF.3 | Estimate size and time before bulk generation | done, PR #238 (stats in control schema v6) |
| 87 | 2026-10-01 03:26Z to 05:29Z | UX.3 | Warn visitors while a background job changes the galaxy | done, PR #2086 |
| 88 | 2026-10-01 03:26Z to 05:29Z | PERF.4 | Second progress bar for slow plan layers | done, PR #258 |
| 89 | 2026-10-01 03:36Z to 05:29Z | PERF.5 | Scatter bright stars in stages | done, PR #229 |
| 90 | 2026-10-01 03:46Z to 05:29Z | PERF.6 | Rate-limit SQL calls, do more per call | done: investigation, then PR #222, #223 and #225 (PERF.8 caps the writers) |
| 91 | 2026-10-01 03:46Z to 05:29Z | PERF.7 | Parallelize sector and system generation | done, PR #225 and #227 |
| 92 | 2026-10-01 03:46Z to 05:29Z | PERF.8 | Parallel background work queue in the API | done, PR #225 and #227 |
| 93 | 2026-10-01 04:50Z to 05:29Z | PERF.9 | Weight the bright-star ETA by the shape of the galaxy | done, PR #258 |
| 94 | 2026-10-01 04:58Z to 05:29Z | PERF.10 | Record generation speed across a log scale of densities | done, PR #238 |
| 95 | 2026-10-01 04:58Z to 05:29Z | PERF.11 | Store each sector's expected and actual density | open |
| 96 | 2026-10-01 05:15Z to 05:29Z | MAP.43 | Wedge lines run past the galaxy's edge (bug) | done, PR #201 |
| 97 | 2026-10-01 05:15Z to 05:29Z | MAP.47 | Bright stars vanish when zoomed out (bug) | done, PR #201 |
| 98 | 2026-10-01 05:15Z to 05:29Z | MAP.48 | Stars take a while to appear after a zoom (bug) | done, PR #201 |
| 99 | 2026-10-01 05:15Z to 05:29Z | MAP.37 | Generated systems are hard to find on the map (bug) | done, PR #201 |
| 100 | 2026-10-01 05:22Z to 05:29Z | MAP.17 | No free camera: drill down from a top-down view by wedge, slice and block (bug) | done, PR #208 |
| 101 | 2026-10-01 05:24Z to 05:29Z | MAP.26 | "Show on Galaxy Map" opens at the sector; map Back and Forward (bug) | done, PR #208 |
| 102 | 2026-10-01 05:27Z to 05:29Z | UX.13 | One meaningful-unit ladder for speeds | done, PR #234 |
| 103 | 2026-10-01 05:27Z to 05:29Z | UX.14 | One meaningful-unit ladder for time periods | done, PR #234 |
| 104 | 2026-10-01 05:40Z to 05:54Z | MAP.45 | Rogue planets (and maybe other objects) drawn outside the sector's wireframe (bug) | done, PR #200 |
| 105 | 2026-10-01 05:43Z to 05:54Z | UX.15 | Put an object's data beside its 3D render when there's room (bug) | done, PR #195 |
| 106 | 2026-10-01 05:43Z to 05:54Z | GEN.8 | Give rogue planets a planet class, with a rogue flag in the class constants | done, PR #245 |
| 107 | 2026-10-01 05:44Z to 05:54Z | MAP.46 | Rogue planets are hard to find on the Sector Map (bug) | done, PR #200 |
| 108 | 2026-10-01 05:44Z to 05:54Z | MAP.15 | Stars and glowing phenomena as points of light on the Sector Map | done, PR #234 |
| 109 | 2026-10-01 05:50Z to 05:56Z | MAP.30 | Slab list to the left of the map, and a 3:4 map | done (slab slider to the right of a 4:3 map), PR #234 |
| 110 | 2026-10-01 05:50Z to 05:56Z | MAP.44 | Wedge lines and ring circles run far past a zoomed-in block (bug) | done, PR #208 |
| 111 | 2026-10-01 05:50Z to 05:56Z | MAP.18 | The block under the pointer is too hard to see from above (bug) | done, PR #208 |

## Reverse lookup: new ID to old numbers

Parents marked "new parent" had no old number of their own.

| New ID | Title | Old numbers (date range, UTC) | Status |
|---|---|---|---|
| ADM.1 | Admin editing: overrides, delete, regenerate (new parent) | none | done (all subitems: PRs #235, #244, #260) |
| ADM.2 | Generate a column and a shell | 23 (2026-09-30 18:14Z to 2026-10-01 00:34Z) | done in 7.27.0, PR #143 |
| ADM.3 | Generate buttons on unfilled sectors (Sector Map, Galaxy Map) | 24 (2026-09-30 18:14Z to 2026-10-01 01:58Z) | done in 7.27.0 (PR #143) and 7.34.0 (PR #151) |
| ADM.4 | Collapsible Generate page sections; pick the center sector | none | done, PR #279 and #301 |
| ADM.5 | Central validate module | 57 (2026-10-01 01:15Z to 05:29Z) | done, PR #235 |
| ADM.6 | Override a planet's or moon's class | 58 (2026-10-01 01:15Z to 05:29Z) | done, PR #260 |
| ADM.7 | Override a star | 59 (2026-10-01 01:15Z to 05:29Z) | done, PR #260 |
| ADM.8 | Delete and regenerate buttons, sector down | 60 (2026-10-01 01:15Z to 05:29Z) | done, PR #244 |
| ADM.9 | "Place a facility": host by placement, log orbit slider, moving belt facilities (bug) | none | done, PR #211 |
| ADM.10 | Admin page to view and manage the work queue | none | done, PR #294 |
| ADM.11 | Jobs keep running after the browser closes | none | done, PR #302 |
| ADM.12 | Jobs as a tree, with timing for every node | none | done, PR #285 |
| ADM.13 | Incomplete uploads page | none | open |
| ADM.14 | Line up the Generate page's text boxes, not their headings (bug) | none | done, PR #544 |
| ADM.15 | Change the worker count from the Queue page, with a "Ludicrous Speed" mode | none | done, PR #1506 |
| ADM.16 | Prevalence controls on the Generate page | none | done, PR #506 |
| ADM.17 | The Generate page shows the galaxy's seed and version | none | dropped (Boss, 2026-10-09 20:42Z) |
| ADM.18 | The galaxy's creation settings saved as a JSON file, downloadable from the Admin dashboard | none | done, PR #816 |
| ADM.19 | The Admin dashboard lists the 18 settings backups for download | none | dropped (Boss, 2026-10-09 20:42Z) |
| ADM.20 | A "merge now" button on the Admin dashboard (low priority) | none | dropped (Boss, 2026-10-09 20:42Z) |
| ADM.21 | Input validation on Pydantic models | none | done, PR #748 |
| ADM.22 | Job logs streamed over SSE into Xterm.js, with native progress bars | none | done, PR #551 |
| ADM.23 | Log output wraps with hard line breaks (bug) | none | done, PR #508 |
| ADM.24 | A failed action's log closes before it can be read (bug) | none | done, PR #560 |
| ADM.25 | Error tracebacks don't reach the console and the web log window (bug) | none | done, PR #560 |
| ADM.26 | The bright-star backfill shows no progress bar on the web (bug) | none | done, PR #560 |
| ADM.27 | Changing a planet's class doesn't regenerate its surface conditions (bug) | none | done, PR #448 |
| ADM.28 | A simpler Generate page: layer specs, a Customize window and plain controls | none | done, PR #991 |
| ADM.29 | Fill a span of layers, rings or columns | none | done, PR #926 |
| ADM.30 | Radial generation: a cylinder of N sectors around a point | none | done, PR #934 |
| ADM.31 | Every generate action offers to show what it made on the Galaxy Map | none | done, PR #976 |
| ADM.32 | Add a star system to a sector: at the emptiest spot, at given coordinates, or at random outside every Hill sphere | none | open |
| ADM.33 | The owner can override "no room" warnings and generate anyway | none | done, PR #476 |
| ADM.34 | One admin menu per screen, holding only that screen's actions | none | done, PR #562 |
| ADM.35 | Full control from every screen: edit anything, regenerate with every input, backfill or erase what is in view | none | open |
| ADM.36 | Change an object's trajectory vector | none | open |
| ADM.37 | The Generate page's prevalence fields show "0% change" instead of each feature's real share (bug) | none | done, PR #514 |
| ADM.38 | Worker exceptions fail on Python 3.9 because _picklable_error calls Exception.add_note, which needs 3.11 (bug) | none | done, PR #627 |
| ADM.39 | Running queue jobs show an ETA (bug) | none | done, PR #800 |
| ADM.40 | The Generate page stops reporting a lost connection (bug) | none | done, PR #800 |
| ADM.41 | Every job is easy to find and cancel in the web queue, and past jobs are paginated and clean (bug) | none | done, PR #800 |
| ADM.42 | One settings model describes every config.json option | none | done, PR #912 |
| ADM.43 | A full configuration page under Admin | none | done, PR #1921 |
| ADM.44 | Web, Open Graph and SEO settings | none | done, PR #1921 |
| ADM.45 | Prevalence fields take the override share directly and must total 100% | none | done, PR #991 |
| ADM.46 | Generate page: a progress line and per-layer counts instead of one line per sector | none | dropped |
| ADM.47 | Generating a neighbourhood from the Generate page shows no per-sector stats and looks slow or silent (bug) | none | done, PR #846 |
| ADM.48 | Two test_api_auth_sweep tests fail: /admin/stats/galaxy-settings/<name> answers 302, not 403, to an unauthorised caller (bug) | none | done, PR #863 |
| ADM.49 | Galaxy shape density settings: the user changes the density range of the spiral arms, the inter-arm space, the core and the bulge | none | done, PR #1116 |
| ADM.50 | Remove the galaxy settings files and everything around them; keep only an internal record of the current galaxy | none | done, PR #1896 |
| ADM.51 | Regenerating with a class the object cannot support asks to confirm, then does it (bug) | none | done, PR #1374 |
| ADM.52 | A star system editor: plan a system as a 2D diagram or from the command line | none | open |
| ADM.53 | Admin login returns an error 500 when another job holds the row it updates, and slow pages time out (bug) | none | done, PR #1318 |
| ADM.54 | A Redis health check on the admin page: answering, append-only file on, memory policy, memory used, last save, worker count (optional) | none | done, PR #1401 |
| ADM.55 | Admin usability audit: walk every admin page and task, rank what is hard to use, and propose one layout | none | done, PR #2498 |
| ADM.56 | Streamline and fix the admin pages from the audit: one admin shell, one navigation and one pattern for tables, forms, confirmations and job logs | none | done, PR #2498 |
| ADM.57 | Headline of 8.6: edit everything, the admin can change anything and the world adjusts after the admin approves a generation plan | none | open |
| ADM.58 | Change plan engine: operations, the effect diff and the plan states, as a library | none | done, PR #3677 |
| ADM.59 | Change plan store: save, load and list plans, and mark them stale | none | open |
| ADM.60 | Apply, undo and audit a plan: one job, one activity-log entry per operation | none | open |
| ADM.61 | Edits survive repair and regeneration: replay the journal after a sector is rebuilt from its seed or from parity | none | open |
| ADM.62 | Aspect catalogue for systems, stars, planets, moons and belts: every stored field and every generation input, with how they depend on each other | none | open |
| ADM.63 | Aspect catalogue for phenomena: black holes, neutron stars, nebulae, remnants, rogue planets, comets, asteroid fields and quasars | none | open |
| ADM.64 | Keep-or-regenerate resolution: from the admin's marks to the exact set of aspects redrawn | none | open |
| ADM.65 | Locks: keep marks that stay on an object or one of its aspects, shown on the pages and honoured by every plan | none | open |
| ADM.66 | Plan approval screen: every operation with its before and after, approve, reject one, or cancel | none | open |
| ADM.67 | Plan diagram: the system before and after as the 2D diagram, and the sector with moved stars marked | none | open |
| ADM.68 | Plan history: every plan with who asked, approved, applied and undid it, and an undo button | none | open |
| ADM.69 | Keep-or-regenerate matrix: the one screen component every regen screen uses | none | open |
| ADM.70 | Specification form: "regenerate to this specification" built from the aspect catalogue and the system generator options | none | open |
| ADM.71 | Move a star system or phenomenon within its sector or to another sector, as a plan | none | open |
| ADM.72 | Edit a body's orbit: distance, eccentricity, inclination and angles, with its moons following, as a plan | none | open |
| ADM.73 | Density overrides for any region: store them, show them, and apply them to later fills | none | merged into MAP.218 (Boss, 2026-10-16, 11:58 PM EDT: consolidation sweep) |
| ADM.74 | Adjust the world to a density change: a plan that adds or removes systems so a region matches its new density | none | open |
| ADM.75 | Rearrange a system around a change to a planet, moon or belt: the plan that fewest changes make valid again | none | open |
| ADM.76 | Rearrange a system around a change to a star: class, mass, luminosity or companions | none | open |
| ADM.77 | Regenerate a planet, moon or belt to a specification, with a keep set, rearranging the system if it must | none | open |
| ADM.78 | Regenerate a star or a whole system to a specification, with a keep set | none | open |
| ADM.79 | Regenerate a phenomenon to a specification, with a keep set | none | open |
| ADM.80 | Regenerate a sector or a region keeping chosen aspects: locked objects, positions, names | none | open |
| ADM.81 | Edit any stored field of any object as a plan: one generic editor driven by the aspect catalogue | none | open |
| ADM.82 | Unit finalizer: validate a system or sector, write it in one transaction, correct only a clashing name and derived columns, or reject with reasons | none | open |
| ADM.83 | Settings writer: validate, write and back up the user's config.json, with an audit row per change | none | done, PR #1857 |
| ADM.84 | Settings page under Admin: every option with its help, edit, validate and save, with a diff before saving | none | done, PR #1875 |
| ADM.85 | Settings backups, restore, import and export on the settings page | none | done, PR #1905 |
| ADM.86 | Owner-tier settings: Owner-only options behind a recent-password prompt | none | done, PR #3083 |
| ADM.87 | Backup archive: one bz2 tarball per backup holding the settings, the control database and the parity data as parts, with retention and restore of the whole or one part | none | done, PR #1857 |
| ADM.88 | Control database backup part: a dump of the planetgen_control schema (admin logins, API keys, run records) in the backup archive | none | open |
| ADM.89 | Control database restore: load a backed-up planetgen_control into a scratch schema, check it, then replace the live schema, from the command and from the settings page | none | open |
| ADM.90 | Database health and repair page: run the check, see damaged and unprotected sectors, repair from the live parity or from an archive, restore parts, as a queue job | none | open |
| ADM.91 | Backups and repair, target 8.9: one bz2 archive of the settings, the control database and the parity data, and the ways to use it (restore the control database, restore parity, repair the galaxy, the repair page) | none | open |
| ADM.92 | Admin shell: one layout, one navigation and the shared table, form, confirmation and job-log components | none | done, PR #2494 |
| ADM.93 | Move every existing admin page onto the shell and fix the defects the audit ranks above its line | none | done, PR #2494 |
| ADM.94 | Density override store and the fill read: set, undo and combine overrides for any region, and have later fills use them | none | open |
| ADM.95 | A remote run is a job: it shows in the job tree and the job page with its progress and report | none | open |
| ADM.96 | One configuration file: the user's own config.json holds every setting, and the settings.json overlay goes | none | done, PR #1857 |
| ADM.97 | The settings page can show a number as a slider: the cache size is the first one | none | done, PR #2456 |
| ADM.98 | Admin control suite research: how to pack many options into one quiet interface (admin pages, the Generate suite and the settings), a report and decision cards before anything is built | none | done, PR #2411 |
| ADM.99 | Generate suite redesign: the Generate page rebuilt inside the admin control suite so a long list of options is easy to use | none | done, PR #2498 |
| ADM.100 | Settings page redesign: plain-language groups and descriptions, and settings that cannot be changed hidden by a toggle that is on by default | none | done, PR #2456 |
| ADM.101 | Settings: words for every option: label, one-line summary, group, control and advanced flag in the settings model | none | done, PR #2411 |
| ADM.102 | Settings: view classes: editable, owner, server and overridden, with the reason an option is fixed | none | done, PR #2456 |
| ADM.103 | Settings: the new page: grouped rows, a left rail, search, a switch that shows settings that cannot be changed (off by default), a sticky save bar, per-row reset and Details | none | done, PR #2456 |
| ADM.104 | Settings: controls: switch, slider with a number box, duration, rate limit, chips and colour picker | none | done, PR #2456 |
| ADM.105 | Settings: the review page in plain language: 'Cache size: 1,024 MB to 2,048 MB', the restart list and the overridden list | none | done, PR #2456 |
| ADM.106 | Settings: the server-files disclosure and the overridden-option badge | none | done, PR #2456 |
| ADM.107 | Generate: the frame: status strip, pinned current-job panel (idle shows the last result), History tab and the Work Queue link | none | done, PR #2494 |
| ADM.108 | Generate: the task flow engine: a Start screen with six cards, a step rail, the draft carried in the URL, a Review step with the estimate, summary and command line, and Run | none | done, PR #2494 |
| ADM.109 | Generate: the Add sectors task: two choices, the four-tab place picker and how-big with a live count | none | done, PR #2494 |
| ADM.110 | Generate: the Start a new galaxy task: seed-decides shape fields, the first neighbourhood and a typed confirmation | none | done, PR #2498 |
| ADM.111 | Generate: the Scatter task: rebuild, redo some, add a dimmer step, with one Advanced block of scatter limits | none | done, PR #2498 |
| ADM.112 | Generate: prevalence, star mix and overrides as one shared Advanced block, rendered once | none | done, PR #2498 |
| ADM.113 | Generate: the one-off system on the frame, with tri-state rows as segmented controls | none | done, PR #2494 |
| ADM.114 | Generate: the Check and Wipe tasks and the History tab | none | done, PR #2498 |
| ADM.115 | Generate: Run again from History | none | done, PR #2494 |
| ADM.116 | Generate: remove the old Generate page pieces and update docs/reference/html-interface.md and the architecture table | none | done, PR #2498 |
| ADM.117 | Generate: browser tests for every task end to end, in both themes, without script and at phone width | none | done, PR #2498 |
| ADM.118 | Admin registry: pages, tasks and settings declare label, summary, keywords, role, danger level and URL, and the rail, hub and search read it | none | done, PR #2494 |
| ADM.119 | Job strip on every admin page: a thin strip that shows a running job, Cancel and the log | none | done, PR #2494 |
| ADM.120 | Command palette and /admin/search: server search over the registry, Ctrl+K and /, results for pages, tasks, settings, jobs and sectors | none | done, PR #2494 |
| ADM.121 | Typed commands that open a task's Review step pre-filled and never run by themselves | none | done, PR #2494 |
| ADM.122 | Admin keyboard set and a help overlay | none | done, PR #2494 |
| ADM.123 | Queue one next job: a single next slot that starts when the running job ends, with Cancel, shown in the job strip and History | none | done, PR #2494 |
| ADM.124 | Generate: the shape dials panel, built from the lever table: five dials of sliders with Auto, the pinned-versus-drawn read-out, any-two-of-three for mass, size and density, and the Shape dropdown | none | done, PR #2494 |
| ADM.125 | Generate: named recipes saved per admin in the control database | none | done, PR #2646 |
| ADM.126 | Generate: a live 'how big' count while typing a command or filling the form | none | done, PR #2677 |
| ADM.127 | Generate form for a new galaxy: tighter star brightness section, sliders above their boxes, a box that invokes no scatters, one degrees on arm pitch, Edge as one typed value with units (bug) | none | done, PR #2608 |
| ADM.128 | The admin one-off system page and the API schema offer the star count and the cap (stars and max_stars) | none | open |
| ADM.129 | Admin page buttons that start the habitability re-score and the map luminosity fill through the queue | none | done, PR #3181 |
| ADM.130 | Aspects for civilizations: species, nations, cities, ports and stations in the aspect catalogue, with keep marks and locks | none | open |
| ADM.131 | Editable civilizations: edit a species, rename or split a nation, move a city, transfer a city between nations, set the dominant nation and place or remove a port, as approved plans | none | open |
| ADM.132 | After an admin edit of a planet or moon, its habitability, Suitability and life scores are recomputed, and a test proves it (check) (bug) | none | done, PR #3268 |
| ADM.133 | Admin edits of resource locations: change a belt, field or body's resource grades and scores as an aspect of the edit-everything plan, with a lock so regeneration keeps them | none | open |
| ADM.134 | /admin/queue/<run_id>/log answers 302 to a visitor, but the auth sweep expects 403 for a json_only view: test_admin_pages_send_an_admin_on_seeded_credentials_to_account and test_admin_pages_send_anonymous_visitors_to_login fail on clean main (bug) | none | done, already fixed on main (admin_queue_log aborts 403 for a visitor; Core 5 confirmed, no change) |
| ADM.135 | The data check in front of a Generate job no longer fails because another check holds the check place, and a refusal names the run that holds it (bug) | none | done, PR #4675 |
| ADM.136 | A New galaxy job stalls waiting on a generation run started by another process, and the holder is not identified (bug) | none | open |
| ADM.137 | Plan operations crash: admin/plans.py calls stabilize_star_system with a pinned keyword it no longer accepts, so 14 test_plans tests fail (bug) | none | open |
| API.1 | Create a system inside an existing sector | 10 (2026-09-24 01:32Z to 02:18Z); 7 (2026-09-24 01:57Z to 02:02Z); 5 (2026-09-24 02:25Z to 05:38Z); 4 (2026-09-24 02:53Z to 2026-09-30 17:58Z); 13 (2026-09-30 16:44Z); 14 (2026-09-30 16:49Z); 15 (2026-09-30 16:51Z to 18:09Z); 32 (2026-09-30 18:14Z); 37 (2026-09-30 18:39Z to 2026-10-01 04:24Z) | done in 7.43.0, PR #170 |
| API.2 | Edit a system's generated content | 11 (2026-09-24 01:32Z to 02:18Z); 8 (2026-09-24 01:57Z to 02:02Z); 6 (2026-09-24 02:25Z to 05:38Z); 5 (2026-09-24 02:53Z to 2026-09-30 17:58Z); 14 (2026-09-30 16:44Z); 15 (2026-09-30 16:49Z); 16 (2026-09-30 16:51Z to 18:09Z); 33 (2026-09-30 18:14Z); 38 (2026-09-30 18:39Z to 2026-10-01 04:24Z) | done in 7.43.0, PR #170 |
| API.3 | Remote generate: generate on a local machine, upload through the API | none | open |
| API.4 | API compatibility data in the docs | none | done, PR #3543 |
| API.5 | API version and compatibility checking | none | open |
| API.6 | User-level API keys, owned by the account that created them, that can read but not upload | none | open |
| API.7 | Investigate and plan upload limits | none | done, PR #1849 |
| API.8 | Verify uploaded data before it is finalized | none | open |
| API.9 | Key scopes | none | done, PR #1033 |
| API.10 | Reservations: claimed sectors and id blocks per run | none | open |
| API.11 | Staging tables | none | open |
| API.12 | The download: seed, skeleton and name state | none | open |
| API.13 | Generation without a database | none | open |
| API.14 | Upload routes, compressed, in batches | none | open |
| API.15 | Log every API call with its user, how it came in, and its HTTP response code | none | done, PR #3136 |
| API.16 | The API reports the galaxy's seed, version and run history | none | dropped (Boss, 2026-10-09 20:42Z) |
| API.17 | Remote generation reproduces what the server would make | none | open |
| API.18 | Generate by recipe: JSON for sectors, systems, planets, moons and phenomena | none | open |
| API.19 | Galaxy-scale recipes: build a whole galaxy, piece by piece, from JSON | none | open |
| API.20 | `require_json_body` returns 500 for a deeply nested JSON body (bug) | none | done, PR #865 |
| API.21 | Flask-Limiter puts `Retry-After` on successful responses (bug) | none | done, PR #865 |
| API.22 | An API version number: one sequential integer, shown in admin and in the status response | none | done, PR #1119 |
| API.23 | The object ID as the public reference: pages, URLs, the API, wiki links and objectref use it in place of row ids | none | done, PR #1133 |
| API.24 | GET /api/sectors/<id>/gravity: the per-zone gravity grid as JSON | none | open |
| API.25 | A structured description endpoint, and Markdown and wikitext on request only | none | done, PR #1527 |
| API.26 | Bump API_VERSION for the shell, band and wedge sector address | none | done, PR #3047 |
| API.27 | A read-only profile lookup for any address in the reach, so an uncharted sector can be previewed | none | done, PR #1919 |
| API.28 | Recipe schema and runner: a recipe fixes, ranges or randomises any aspect, and a library turns it into objects | none | open |
| API.29 | Upload spool and unit states: received batches on disk, units assembled at seal, finalize as a queue job | none | open |
| API.30 | generate.py --remote: the client driver for sector and galaxy runs | none | open |
| API.31 | Every outgoing web request sends a standard User-Agent with the software, its version and the API version | none | done, PR #3047 |
| API.32 | test_fuzz_web_routes fails on a machine that has a real database named planetgen_control: every galaxy route answers 500 for db=planetgen_control (bug) | none | done, PR #2659 |
| API.33 | API field descriptions, the help pages and the CLI output: every number states its unit (scan) (bug) | none | open |
| API.34 | API and remote uploads carry resource locations: the resource fields are in the API, and an uploaded galaxy is checked for them | none | open |
| API.35 | Region names on the wire: a name on charted nodes and the stage and node endpoints, galaxy-map-wire-format.md (regions chain, step 2) | none | done, PR #3396 |
| API.36 | Optional: a region_names search table and a "Go to" box, refilled when the naming key changes (regions chain, step 6) | none | open |
| API.37 | API fields for the life class: life_class, variant, best, best score, per-class scores and the animal grade, with an API version bump | none | done, PR #3871 |
| DB.1 | Starbases, colonies and outposts in the database | 30 (2026-09-30 18:14Z); 35 (2026-09-30 18:39Z to 2026-10-01 02:24Z) | done in 7.35.0, PR #152 |
| DB.2 | Asteroid field and comet composition rows are written but never read (bug) | none | done, PR #347 |
| DB.3 | resetDb while another process holds id blocks can duplicate primary keys (bug) | none | done, PR #347 |
| DB.4 | A database with an emptied schema_migrations table is treated as current (bug) | none | done, PR #342 |
| DB.5 | Several first connections to an empty database race to create the schema (bug) | none | done, PR #342 |
| DB.6 | Store the galaxy's 128-bit seed, the version that made it, and every generation run | none | done, PR #387 |
| DB.7 | The version that generated each sector, and a warning for mixed-version galaxies | none | done, PR #813 |
| DB.8 | Check a galaxy database and say whether it is damaged | none | done, PR #893 |
| DB.9 | Repair a damaged galaxy database from a parity file | none | done, PR #3460 |
| DB.10 | Repair reads the newest settings JSON and the pending deltas | none | dropped (Boss, 2026-10-09 20:42Z) |
| DB.11 | The database layer and migrations on SQLAlchemy and Alembic | none | done, PR #751 |
| DB.12 | MariaDB reads a stored -0.0 back as +0.0 (bug) | none | done, PR #442 |
| DB.13 | Every stored value in its own column, not in JSON blocks, and indexed for search | none | done, PR #766 |
| DB.14 | Sector stats keep the raw star statistics (mean age, summed luminosity, star count) and no baked color | none | done, PR #556 |
| DB.15 | A migration progress bar with the time remaining | none | done, PR #924 |
| DB.16 | Store the generator epoch and run id on each sector instead of four version text columns | none | done, PR #2522 |
| DB.17 | Repair by regenerating a damaged sector from its seed when parity cannot rebuild it | none | open |
| DB.18 | Migration helpers for slow DDL: online indexes, instant columns and batched updates | none | done, PR #1901 |
| DB.19 | Compact or derive phenomenon rows (needed only if the mass cut is lowered to 10 solar masses or less) | none | closed, not needed at the 20 solar mass cut (PR #866); reopen if the cut is lowered to 10 solar masses or less |
| DB.20 | Object IDs in the schema: uid becomes BINARY(10), unique on its own, plus an id_counters table | none | done, PR #1055 |
| DB.21 | A deep pass for the database check: validate every star system, with the estimated time shown first | none | done, PR #1082 |
| DB.22 | Sectors saved by hand without a grid address have no stored object ID until they are saved again (bug) | none | done, PR #1127 |
| DB.23 | Store a "map luminosity" for phenomena that are faint but massive (Alembic migration) | none | done, PR #1365 |
| DB.24 | The nebula_field table: one row per field nebula cloud, built at plan time (Alembic migration) | none | done, PR #1317 |
| DB.25 | A per-level aggregate table for the 3-ary block ladder so a Galaxy Map window needs one query (Alembic migration) | none | done, PR #1376 |
| DB.26 | One migration for the habitability columns: flux variation and the air, water and suitability ratings | none | done, PR #1423 |
| DB.27 | One migration for the species tech-level columns: six indices, tech level and tech source | none | done, PR #1446 |
| DB.28 | One migration for the facility type tables and the new facility columns | none | done, PR #1476 |
| DB.29 | update.sh and every database migration refuse to start while jobs are active, unless --yes is passed or the user confirms cancelling them | none | done, PR #3645 |
| DB.30 | One migration for the per-species environment limits | none | done, PR #2520 |
| DB.31 | Store structure, not prose, for life and flavour | none | done, PR #2206 |
| DB.32 | Schema: rename ring, layer and slot to shell, band and wedge, and convert the existing galaxy in place without moving anything | none | done, PR #1607 |
| DB.33 | MySQL 8.0.46 reflects the foreign keys fk_facilities_type and fk_facility_field_values_field as RESTRICT, so check-db calls a healthy database damaged (bug) | none | done, PR #2436 |
| DB.34 | Cell key and covering index on every table that holds a position, in the shared 0093 migration | none | done, PR #1757 |
| DB.35 | Scatter rows keep their address columns and gain the point columns and the cell key | none | done, PR #1757 |
| DB.36 | Tables for change plans, their operations and the edit journal (migration only) | none | done, PR #2523 |
| DB.37 | Lock marks table (migration only) | none | done, PR #2524 |
| DB.38 | Region density overrides table (migration only) | none | done, PR #2525 |
| DB.39 | Edit request support tables (migration only): request comments and notification records | none | done, PR #2527 |
| DB.40 | Parity data backup part: the parity file of the DB repair copied into the backup archive | none | open |
| DB.41 | Parity data restore: put a backed-up parity file back and check it against the database | none | open |
| DB.42 | Repair the galaxy from an archive: `repair-db --from-archive` uses the parity part of a backup when the live parity file is missing or damaged | none | open |
| DB.43 | Reed-Solomon codec in the standard library: planetgen/util/rscodec.py | none | done, PR #2423 |
| DB.44 | Migration only: the sector storage checksum and the parity dirty flag | none | done, PR #2420 |
| DB.45 | Storage checksum and the dirty flag: the export, the checksum and every writer | none | done, PR #3414 |
| DB.46 | The parity file: groups, records, a delta update on every save and the stale-record rule | none | done, PR #3451 |
| DB.47 | planetgen repair-db: rebuild damaged sectors from the parity file | none | done, PR #3456 |
| DB.48 | repair-db: rebuild the table when a clustered-index page is unreadable | none | done, PR #3460 |
| DB.49 | Migration only: the leaf digest column and the shell-band tree table | none | done, PR #2521 |
| DB.50 | Migration only: structure columns for life and flavour | none | done, PR #2206 |
| DB.51 | Generation stores structure, not prose; the readers and the string searches go | none | done, PR #2206 |
| DB.52 | Migration only: drop the prose columns and tables | none | done, PR #2206 |
| DB.53 | Built-system star storage: the orbit_nodes table, stars.parent_node_id, planets.host_node_id, the catalogue columns and the disk descriptor (one Core 1 revision) | none | done, PR #1952 |
| DB.54 | The three scatter-limit settings as galaxy_shape columns: core luminosity floor, core mass limit and halo luminosity floor | none | done, PR #1887 |
| DB.55 | One additive migration for the 8.3 schema: point columns for every kind, heaviest object per sector, asteroid-field provenance | none | done, PR #1952 |
| DB.56 | The orbital_events table (migration only): collision outcomes, mergers, splits, ejections, captures, unbinds and warnings | none | done, PR #2505 |
| DB.57 | The update records what it moved: a replace-each-run table of previous and new positions | none | open |
| DB.58 | Control-schema revision for accounts: role, email and status, the single-Owner index, account_tokens and invite_uses | none | done, PR #2536 |
| DB.59 | Nebula age and ionizing rate, and remnant age, stored | none | done, already built (v87, DB.26 and GEN.150), found by Core 2 2026-10-16 6:05 PM EDT |
| DB.60 | Migration: record the density multiplier used on each filled sector (migration only) | none | done, PR #2526 |
| DB.61 | Migration: unnamed clusters and a star-to-cluster parent link (migration only) | none | done, PR #2528 |
| DB.62 | Upload bookkeeping tables in the control database (migration only) | none | done, PR #2537 |
| DB.63 | Table for ships and stations, user-owned or generated (migration only) | none | done, PR #2530 |
| DB.64 | Table for stored constellations (migration only) | none | done, PR #2529 |
| DB.65 | Clamp the cell key index at 2^21-1 in `cellkey.py` (`cell_index`, `cell_keys`) and in the SQL expression the migration generates | none | done, PR #2025 |
| DB.66 | The shape row carries the galaxy family and one column per density component, as a migration only | none | done, PR #2439 |
| DB.67 | Planets and moons store the albedo, the greenhouse multiplier and the molar density they were built with | none | done, PR #3474 |
| DB.68 | Planets and moons store the ozone column and the haze, dust and cloud optical depths | none | done, PR #2856 |
| DB.69 | Planets and moons store the composition fractions | none | done, PR #2751 |
| DB.70 | Planets, moons and small bodies store the surface and small-body composition shares | none | done, PR #3920 |
| DB.71 | Planets and moons store the partial pressures of the added gases, not only the ten stored today | none | done, PR #3110 |
| DB.72 | The scatter tables hold band_index as SMALLINT, which overflows from shell 20,861 outward, so the whole-sphere scatter stops at 97 percent (bug) | none | done, PR #2060 |
| DB.73 | Storage for one orbited-body model: a shared parent link for every kind, free parents, the orbit node table filled, moons storable under a rogue | none | done, PR #2863 |
| DB.74 | One generic insert and load driven by the column list, one registry of kind tables, and the shared core table with three profile tables | none | open |
| DB.75 | Object IDs stay stable: a refill that changes a system's planet or moon count keeps every surviving body's ID, and a body that changes parent keeps its ID | none | open |
| DB.76 | Column registry: every numeric column declares its table, kind and stored unit, loaded and stored through one list | none | done, PR #2487 |
| DB.77 | Precision table: each length column declares its reference frame and the largest distance at which a 10 m body keeps a metre | none | done, PR #2534 |
| DB.78 | Central-object offsets replace the wobble: the reflex_offset columns, calculate_reflex_offset, the SQL mirror and the nine binary position columns go | none | open |
| DB.79 | The object index table (migration only): one row per object of every kind the site lists, with keyset order keys and a name index | none | done, PR #2361 |
| DB.80 | `VALUES ROW(...)` key-run lists are rejected by MariaDB 10.11 although INSTALL.md supports MariaDB 10.4+ and MySQL 8.0.16+ (bug) | none | done, PR #2406 |
| DB.81 | Shape dials: the plan stores which sliders were pinned and at what value, so the galaxy page can say 'mass pinned at X, the rest drawn' | none | done, PR #2403 |
| DB.82 | The object index builder uses VALUES() in ON DUPLICATE KEY UPDATE, which test_no_new_values_function_in_on_duplicate_key_update rejects (bug) | none | done, PR #2390 |
| DB.83 | `galaxy_block_stats` has no maintainer: nothing keeps the block totals current after a scatter, a fill or an edit (bug) | none | done, PR #3476 |
| DB.84 | Every schema change on `sectors` rebuilds the whole table because of its FULLTEXT index: decide how the later DB items pay for it | none | done, PR #2984 |
| DB.85 | Globular clusters and superbubbles stored as containers: which container holds a point, and what a container holds | none | done, PR #2582 |
| DB.86 | Six columns in the column registry are declared 'unspecified': atm_density, atm_molar_density, formed_at, parabolic_mean_anomaly and two more | none | done, PR #2801 |
| DB.87 | The precision table lists the vessels table but the documented floor table does not: test_precision fails on main (bug) | none | done, PR #2541 |
| DB.88 | The vessels table has an Object ID but is not in db/kind_tables.py: test_kind_tables fails on main (bug) | none | done, PR #2545 |
| DB.89 | Giants store their interior layers: the layer boundaries, the helium share per layer, the core flag, core mass and core composition | none | done, PR #3089 |
| DB.90 | Planets and moons store the derived class and a reclassification step keeps it true after an edit or a rebuild | none | open |
| DB.91 | Planets and moons store the four life scores and drop phi_bio, phi_cpx and phi_tech | none | done, PR #2930 |
| DB.92 | Planets and moons store the depth to liquid water and the thickness of the water-bearing layer (aquifers) | none | done, PR #3221 |
| DB.93 | Planets and moons store the equipment levels and regional temperatures and drop equipment_tier; the API serves one equipment object | none | done, PR #2766 |
| DB.94 | The object index, resolve_object and the search panels read the moons of free planets | none | done, PR #2881 |
| DB.95 | The sector_objects view includes the moons of free planets | none | dropped (lane B finding, 2:53 AM EDT 2026-10-24: nothing reads the sector_objects view, schema.sql is applied before migrations so a view naming moons.rogue_planet_id would break opening an older database, and DB.94 lists free moons by sector) |
| DB.96 | Store the atmosphere colour and the cloud cover on planets and moons | none | done, PR #3393 |
| DB.97 | Remove store.detect_schema_version and _VERSION_MARKERS: a database with no version row is refused with a clear message | none | done, PR #2905 |
| DB.98 | Re-baseline the Alembic migration history at the schema the server is on: delete the 72 old-schema fixtures and the per-revision migration tests, keep one upgrade-from-the-previous-release-to-head test | none | done, PR #2931 |
| DB.99 | Store the ground colour and the ice fraction on planets and moons | none | open |
| DB.100 | Store all the life and evolutionary data of a world in the database and build the life paragraph from the stored results | none | done, PR #3184 |
| DB.101 | Storage and units for every physics rectification output: a column, a column_units entry and a migration shipped with its item | none | open |
| DB.102 | The version key stored with the galaxy seed carries the API version, the database schema and the generator epoch after the release version, with the zero padding shortened to make room | none | done, PR #3323 |
| DB.103 | Storage for nation-states and cities: tables, Object IDs, column units and the migration | none | open |
| DB.104 | Storage for a government's territory: homeworld system, growth state and a territory history, ready for 9.2 to make it dynamic | none | open |
| DB.105 | Storage for world population and the global trade route list: tables, Object IDs, column units and the migration | none | open |
| DB.106 | Storage for resource locations: grade, mass estimate and score per resource for belts, fields, giants and icy bodies, with column units and the migration | none | open |
| DB.107 | Region stats: the region_stats table, its save and delete hook and the charted percent (regions chain, step 1b) | none | done, PR #3383 |
| DB.108 | Total mass per sector and per region: mass_sum on the sector block stats row and on region_stats, kept by the same save and delete hook, with a backfill command (regions chain, step 1e) | none | done, PR #3465 |
| DB.109 | Optional parity follow-ups: refresh the parity file after a generation run, and keep a before-image on the edit paths so an edit is a parity delta, not a group recompute | none | open |
| DOC.1 | Number TODO items by category (this renumbering) | 80 (2026-10-01 02:41Z to 05:29Z) | done in the docs-refresh PR (version-scheme questions moved to OPS.1) |
| DOC.2 | Architecture document | 81 (2026-10-01 02:51Z to 05:29Z) | done in the docs-refresh PR |
| DOC.3 | Design documents current, with reasons | 82 (2026-10-01 02:51Z to 05:29Z) | done in the docs-refresh PR |
| DOC.4 | Correct the stale statements the research found in docs, docstrings and comments | none | done, PR #1849 |
| DOC.5 | Rewrite the object ID docs: object-ids.md, database-schema.md and api.md | none | done, PR #1075 |
| DOC.6 | A static help section in the web interface: page template, index, per-page help links and a coverage test | none | done, PR #4651 |
| DOC.7 | The Galaxy Map help page: layers, zoom, fly-through, Color by, select modes, bookmarks and the locate box | none | open |
| DOC.8 | The sector help pages: the sector list, a sector page, the sector map and the sector scene | none | open |
| DOC.9 | The star system help pages: the system list, a system page, the system map and the planets, moons and belts shown | none | open |
| DOC.10 | The search and navigation help pages: search, nearby, the nav page and routes | none | open |
| DOC.11 | The reference browser help pages: species, polities, object classes and phenomena | none | open |
| DOC.12 | The account help pages: signing in, two-factor, the account page, API keys and bookmarks | none | open |
| DOC.13 | The Generate page help pages: layer specs, spans, radial fills, directives, one-off systems and jobs | none | open |
| DOC.14 | The admin help pages: the queue, the stats page, settings, lockouts and the naming key | none | open |
| DOC.15 | A glossary and units help page: coordinates, scales, sector paths, object IDs, time and the in-universe wording | none | open |
| DOC.16 | The API help page for visitors: what the API is, how to get a key, and where the reference lives | none | open |
| DOC.17 | Describe the object ID in api.md | none | done, PR #1133 |
| DOC.18 | Last pass: re-run the documentation sweep when the other DOC items are done and Phase 1 is complete | none | open |
| DOC.19 | After UX.102: document the log output standard for the CLI, the Generate page and admin logs | none | done, PR #4453 |
| DOC.20 | After PERF.86: rewrite the queue and scatter notes and the queue docstrings | none | done, PR #3661 |
| DOC.21 | Rewrite the package docstrings that are wrong or missing, and the queryDb mentions | none | done, PR #1698 |
| DOC.22 | Update the orbital, slow-reads, login, ops, test-suite and db-check notes to their built state | none | done, PR #1599 |
| DOC.23 | Update the nebula, object ID and map-data-layer notes to what merged | none | done, PR #1599 |
| DOC.24 | Bring Python 3.12, numba and astropy into the research notes and three source comments | none | done, PR #2086 |
| DOC.25 | Refresh docs/design/architecture/design-decisions.md: mark the reversed and built decisions | none | done, PR #1594 |
| DOC.26 | Give every design note a status line and refresh the stale ones | none | done, PR #1596 |
| DOC.27 | Update the architecture.md flows and add a design-notes index at docs/design/README.md | none | done, PR #1628 |
| DOC.28 | Rewrite the code map in docs/design/architecture/architecture.md for the planetgen packages | none | done, PR #1628 |
| DOC.29 | A docs link check in the test suite, and a note that the web interface has no help or about page yet | none | done, PR #1581 |
| DOC.30 | Fix the system-file-format.md class table, the moon blacklist and the example-systems.md module names | none | done, PR #1586 |
| DOC.31 | Fix docs/reference/html-interface.md: routes, header sections, search parameters, error pages and rate limits | none | done, PR #1625 |
| DOC.32 | Update database-schema.md: control schema, missing tables and columns, old module names | none | done, PR #1625 |
| DOC.33 | Rewrite the database-schema.md Versioning, Persistence and Planned changes sections for Alembic and v83 | none | done, PR #1625 |
| DOC.34 | Bring docs/reference/api.md up to date: endpoints, object references, errors and the Not done yet list | none | done, PR #1623 |
| DOC.35 | Document the helper command-line tools: orbits, warm_map, worker, query, migrate, reset, lockouts, dedupe, config | none | done, PR #1588 |
| DOC.36 | Sync docs/reference/cli.md with the 11 subcommands, the options and the real bright-star defaults | none | done, PR #1616 |
| DOC.37 | Rewrite the install docs for the venv install, in or right after OPS.45 | none | done, PR #1594 |
| DOC.38 | Refresh docs/reference/config.md prose after PERF.64 and PERF.24: page_cache, stored_counts, PLANETGEN_WORKERS | none | done, PR #1586 |
| DOC.39 | Fix the README facts, the documentation table and the headline rule (Boss decision) | none | done, PR #1590 |
| DOC.40 | Bring docs/reference/testing.md and docs/reference/ci-runners.md in line with ci.yml and scripts/test.py | none | done, PR #1581 |
| DOC.41 | Fix the activity-log file name in the fail2ban guide and the shipped jail | none | done, PR #1584 |
| DOC.42 | Describe what update.sh and install.sh really do today in INSTALL.md and the Apache and macOS guides | none | done, PR #1584 |
| DOC.43 | Fix the Python 3.12 and numba leftovers in the install and CI docs and the 3.10 check in the installer | none | done, PR #1584 |
| DOC.44 | Add Redis to the install and deploy docs: it is required, and no guide lists it | none | done, PR #1584 |
| DOC.45 | Reference docs follow the description layer | none | done, PR #2240 |
| DOC.46 | Docs for the spherical sector grid | none | done, PR #1627 |
| DOC.47 | Rebuild docs/design/architecture/architecture.html from the new architecture.md | none | done, PR #2037 |
| DOC.48 | Docs: the reach, the halo, the profile rules and the Wolf-Rayet and nebula-contents rules in the design notes | none | done, PR #3537 |
| DOC.49 | Docs: point the rate, anomaly, multi-star and nebula notes at the phenomena catalogue | none | merged into DOC.50 (Boss, 2026-10-16, 11:58 PM EDT: consolidation sweep) |
| DOC.50 | Docs: the schema, the multi-star note and the catalogue say where the members live | none | open |
| DOC.51 | Docs: the cell key in the coordinate system and schema docs, the region note moves to built, the API states the cap and the ladder | none | open |
| DOC.52 | A docs index page: links to every HTML page under docs/, each with a description of what it contains | none | done, PR #2154 |
| DOC.53 | The tuning module header names the rule for a new signed constant: the fuzz tests accept only a dict or a positive "fall" | none | done, PR #2070 |
| DOC.54 | Admin help: how to edit, regenerate, lock, approve and undo | none | open |
| DOC.55 | Help: roles, edit requests and how an approval works | none | open |
| DOC.56 | Backups, restore and repair: the deployment guide and the help page | none | open |
| DOC.57 | Help section framework: page template, index, '?' links, help search and the link check | none | done, PR #4651 |
| DOC.58 | Help coverage test: every user-facing route has a help page and every help page has a route | none | open |
| DOC.59 | Docs: the galaxy types note moves to built and the density, potential, schema, API and Generate page notes say the shape holds a family and components | none | open |
| DOC.60 | Settings note: read sections A2 and A6 of the settings and accounts design against the built one-file settings and fix what is stale | none | done, PR #2070 |
| DOC.61 | Docs: how the bubbles shape the stars, the setting values and what the drawn layer now is | none | open |
| DOC.62 | Docs: the generation cost figures and the redo scope are stale (the phenomena cut, rows, bytes and time; `--redo-scatters` leaves satellites and clusters alone) | none | done, PR #2070 |
| DOC.63 | Docs: the atmosphere data model, the gas table and how the values agree | none | done, PR #3907 |
| DOC.64 | Docs: the atmosphere colour, its method and its sources | none | done, PR #2808 |
| DOC.65 | Docs: bulk composition, the materials table and how density, mass and gravity follow | none | done, PR #3913 |
| DOC.66 | Docs: the composition model of stars, bodies and atmospheres, the radiation it derives and the sources of every number | none | open |
| DOC.67 | An architecture map Boss can drill down: project areas, modules and functions with the call flow between them, generated from the code | none | done, PR #2154 |
| DOC.68 | Generate docs/generated/execution-plan.html from scripts/build_todo_docs.py like the other plan pages | none | done, PR #2157 |
| DOC.69 | Global documents are regenerated only by the TODO thread: a stale generated page never fails a lane's gate and a lane's PR never touches one | none | done, PR #2172 |
| DOC.70 | Docs for the unit registry: units-and-number-formatting.md and database-schema.md updated, the unit rules in the registry's file header | none | done, PR #2800 |
| DOC.71 | Docs for the System rule: database-schema.md, the system view in docs/design and the header of system_tree.py | none | open |
| DOC.72 | The README lists every external file the program needs, with its URL and how to install it, and a check keeps the list complete | none | done, PR #3964 |
| DOC.73 | INSTALL.md recommends a larger innodb_redo_log_capacity for a machine that runs the full test gate | none | done, PR #3965 |
| DOC.74 | Design part 01, section 7.4: the CO2 row and its recommended CO2-2 sublimation fit are wrong in the cold; the table now says to use CO2-1 below 195.3 K (bug) | none | done, PR #2820 |
| DOC.75 | Help and docs text for the whole habitability system: the score, the symbols, the lists, the sector colouring and the map markers | none | done, PR #3052 |
| DOC.76 | Docs for the condition-based classes: the rule table, the planet class pages and the help | none | open |
| DOC.77 | Docs for the four life scores: the help text, docs/design/habitability/habitability-index.md and the database schema document | none | done, PR #3126 |
| DOC.78 | One end-to-end draw-order table from the star to the PHI-6 for the contracts file: blocks, streams and the reproducibility rule | none | done, PR #2813 |
| DOC.79 | The compressed-water-layer design document has inconsistent x=0.2 reference rows in part C | none | done, PR #3262 |
| DOC.80 | Progress docs rule: nothing after a progress block may run longer than the bar threshold without its own Step | none | done, PR #3125 |
| DOC.81 | Help page: resource locations, what the scores mean and how to find them on the map | none | open |
| DOC.82 | Region names docs and the audit check "region names stay unique among neighbours" (regions chain, step 7) | none | done, PR #3400 |
| DOC.83 | Plain-words style guide for every sentence the site writes about a planet or a star, with a test that enforces it | none | done, PR #4478 |
| DOC.84 | Help page suite for life: the four kinds of life, every life class, and the scales we use (abundance, trend, security, suitability) | none | done, PR #4506 |
| DOC.85 | Refresh the class table in docs/reference/system-file-format.md: the F, G, B and N meanings and the zone and habitable columns are from the old zone-drawn classes | none | done, PR #4499 |
| DOC.86 | Oxford comma everywhere: check lists of three or more for the serial comma, and fix every template, help page and generated sentence | none | open |
| DOC.87 | References R1: link helper (reference_link(kind, id)), the anchor and URL scheme, the banned-words rule for out-of-universe words, and the crawl test | none | done, PR #4602 |
| DOC.88 | References R2: the References shell and generator (Help menu entry, /references with the Out of universe banner, scripts/build_reference_pages.py with --check, built by update.sh) | none | done, PR #4624 |
| DOC.89 | References R3: Sources pages (papers, datasets, archives, licenses, credits), each paper linking to where it is used | none | open |
| DOC.90 | References R4: stand-in registry pages (every GEN.725 row, one page per stand-in with how it was made and decided, tier and level legends, filters, and the fitted-to-a-result page for L2 to L4 with tuning targets) | none | open |
| DOC.91 | References R5: method pages by area, what we have and lack, and the computed scorecard; /help/properties moves here with no redirect | none | open |
| DOC.92 | References R6: move the speculative chemistry out of the life help and the air chart into /references/speculative, with per-class links (including the multi-class section for GEN.1209 and GEN.1213) | none | open |
| DOC.93 | References R7: strip model hints from the leftover pages (maps, sector estimates text, API docs strings) and add links at the important places | none | open |
| DOC.94 | References R8: test lock: the page crawl proves no out-of-universe word outside References and every link resolves, in the gate for changed pages | none | open |
| DOC.95 | References R3b (very important, first, ahead of the page build): tie every unmatched citation to its real paper, or label it not a paper, in batches by heaviest first | none | done, PR #4603 |
| DOC.96 | Remove every reference to a paper that did not resolve, and rewrite the text so it reads naturally without it | none | done, PR #4638 |
| DOC.97 | Citation links in Help > References: every citation on a References page is a link to its entry on the reference list and to its document information | none | open |
| DOC.98 | Citations live only in the central registry and references.md: move them out of every other page, design document and code comment | none | open |
| DOC.99 | Test lock: a citation written outside the registry and references.md fails the gate | none | open |
| DOC.100 | A stable citation ID for every paper in the registry, shown in the registry, docs/REFERENCES.md and the Papers page | none | done, PR #4648 |
| GEN.1 | Real-world rates for interstellar objects | 5 (2026-09-30 18:14Z to 23:21Z) | done in 7.18.0, PR #136 |
| GEN.2 | Rogue planet mass bins | 6 (2026-09-30 18:14Z to 23:21Z) | done in 7.18.0, PR #136 |
| GEN.3 | A supermassive black hole in every galaxy | 7 (2026-09-30 18:14Z to 22:54Z) | done in 7.15.0, PR #132 |
| GEN.4 | Nebulae, remnants and asteroid fields (new parent) | none | done |
| GEN.5 | Known generation bugs of 2026-09-30 (new parent) | none | done |
| GEN.6 | Finish the correlative update (galactic motion) | 27 (2026-09-30 18:14Z); 32 (2026-09-30 18:39Z to 2026-10-01 02:35Z) | done in 7.37.0, PR #157 |
| GEN.7 | Star population and bright stars (new parent) | none | done |
| GEN.8 | Give rogue planets a planet class, with a rogue flag in the class constants | 106 (2026-10-01 05:43Z to 05:54Z) | done, PR #245 |
| GEN.9 | Plan for more than one galaxy in the database | none | done, PR #1849 |
| GEN.10 | Generate nebulae and remnants with their stars, and map them | 27 (2026-09-30 18:39Z to 2026-10-01 03:53Z) | done in 7.30.0 (PR #144), 7.36.0 (PR #153) and 7.41.1 (PR #148) |
| GEN.11 | Class nebulae and remnants A-W | 28 (2026-09-30 18:39Z to 23:32Z) | done in 7.19.0, PR #138 |
| GEN.12 | Record what sits inside a nebula | 29 (2026-09-30 18:39Z to 2026-10-01 00:35Z) | done in 7.24.0, PR #140 |
| GEN.13 | Names that follow one standard | 30 (2026-09-30 18:39Z to 2026-10-01 01:34Z) | done in 7.31.0, PR #145 |
| GEN.14 | Class asteroid fields | 31 (2026-09-30 18:39Z to 23:32Z) | done in 7.19.0, PR #138 |
| GEN.15 | Moons orbit outside their planet's Hill sphere (bug) | 53 (2026-09-30 19:02Z to 20:27Z); 40 (2026-09-30 20:07Z to 23:25Z) | done in 7.18.1, PR #137 |
| GEN.16 | Moons can orbit inside their planet (bug) | 54 (2026-09-30 19:02Z to 20:27Z); 41 (2026-09-30 20:07Z to 23:25Z) | done in 7.18.1, PR #137 |
| GEN.17 | A close binary's planets inside the binary (bug) | 55 (2026-09-30 19:02Z to 20:27Z); 42 (2026-09-30 20:07Z to 22:43Z) | done in 7.13.1, PR #128 |
| GEN.18 | A planet's Hill sphere overlaps the belt inside it (bug) | 56 (2026-09-30 19:02Z to 20:27Z); 43 (2026-09-30 20:07Z to 23:09Z) | done in 7.16.1, PR #133 |
| GEN.19 | A binary's secondary outweighs its primary (bug) | 57 (2026-09-30 19:02Z to 20:27Z); 44 (2026-09-30 20:07Z to 22:08Z) | done in 7.9.1, PR #123 |
| GEN.20 | Sector growth ignores black holes and neutron stars (bug) | 58 (2026-09-30 19:02Z to 20:27Z); 45 (2026-09-30 20:07Z to 21:53Z) | done in 7.6.1, PR #116 |
| GEN.21 | Star population model: population ages at sector fill | 55 (2026-09-30 23:48Z to 2026-10-01 02:41Z) | done in 7.23.0 (PR #141) and 7.38.0 (PR #159) |
| GEN.22 | Pre-place bright stars at plan time | 56 (2026-09-30 23:48Z to 2026-10-01 02:41Z) | done in 7.38.0, PR #159 (uncertain, see note 2) |
| GEN.23 | Generate a smaller sphere, then backfill bright stars around it per sector block | none | done (schema v49, PR #226) |
| GEN.24 | Generate the galactic core on layer 0 | none | dropped (Boss, 2026-10-16, 3:00 PM EDT; superseded: layers are gone and the core is shell 0, filled by the shell mode) |
| GEN.25 | A moon reclassified after its planet moves can be too large for its planet (bug) | none | done, PR #350 |
| GEN.26 | Rogue planet surface conditions | none | done, PR #263 (schema v48; design docs/design/systems/rogue-planet-surface.md) |
| GEN.27 | Class P (glaciated world) only in the habitable zone, and fitting there | none | dropped (classes are derived from the finished body (cut-over epoch 121); Core 4 classifier re-read 2026-10-30, reports/core4-classifier-reread-2026-10-30.md: the clause is a fixed letter rule the derived classes make moot) |
| GEN.28 | Seven new planet classes in the letter gaps (R, S, U, W, X, Y, Z) | none | dropped (revoked by Boss's full-build sweep 8:16 PM EDT 10-29: new class letters R, U, W, X, Y and Z conflict with conditions-based classes (10-23) and 'all 20 letters stay' (10-28); the part that stands is GEN.963) |
| GEN.29 | Sweep every planet class for sense once the new ones are in (bug) | none | done (sweep report reports/core4-class-sweep-2026-10-30.md, items GEN.1131, GEN.1132, GEN.1133) |
| GEN.30 | Bright-star thresholds: 1000 L_sun galaxy-wide, tiered backfill around generated sectors | none | done, PRs #295, #296 |
| GEN.31 | A point just under layer 0's top face lands in layer 1 (bug) | none | done, PR #353 |
| GEN.32 | Re-running an interrupted bright-star band draws it twice (bug) | none | done, PR #371 |
| GEN.33 | One class per PR, each with its tests | none | dropped (revoked by Boss's full-build sweep 8:16 PM EDT 10-29: new class letters R, U, W, X, Y and Z conflict with conditions-based classes (10-23) and 'all 20 letters stay' (10-28); the part that stands is GEN.963) |
| GEN.34 | Gas and ice giants come out too light, so there are no super-Jupiters (bug) | none | done, PR #350 |
| GEN.35 | Rocky planets only ever get Class D moons (bug) | none | done, PR #350 |
| GEN.36 | Moon regeneration can produce gas-giant or blacklisted moon classes (bug) | none | done, PR #350 |
| GEN.37 | 97% of planets land in the cold zone (bug) | none | done, PR #350 |
| GEN.38 | Rocky rogue planets over 10,000 km are still classed C (bug) | none | done, PR #415 |
| GEN.39 | The same seed can't reproduce the same galaxy (bug) | none | done, PR #381 |
| GEN.40 | Weed out sectors by star density before the bright-star backfill | none | open |
| GEN.41 | Investigate: how much backfill work a density pre-pass would save | none | done, PR #812 |
| GEN.42 | A block-first draw for the backfill (umbrella) | none | open |
| GEN.43 | Don't over-filter: keep bright stars in odd places | none | done, PR #812 |
| GEN.44 | Store each sector's backfill level so finished sectors drop out of any backfill | none | done, PR #425 |
| GEN.45 | Check the rogue planet mix of terrestrial and gas giants (bug) | none | done, PR #350 (mass mix uses dN/dM ∝ M^-0.65 read per log mass; per unit mass gave 87% gas giants) |
| GEN.46 | Star system names of at most two words (bug) | none | done, PR #370 |
| GEN.47 | Nebulae almost never appear (bug) | none | done, PR #419 |
| GEN.48 | Forcing options are impractical for whole sectors; replace them with prevalence controls (bug) | none | done, PR #502 |
| GEN.49 | `+habitable_world` silently fails on hot stars (bug) | none | done, PR #373 |
| GEN.50 | `-planets +asteroid_belt` still makes an asteroid belt (bug) | none | done, PR #373 |
| GEN.51 | Forcing options only for single-system generation | none | done, PR #398 |
| GEN.52 | Prevalence controls for sector and galaxy runs | none | done, PR #502 |
| GEN.53 | The two stars of a binary don't share one age (bug) | none | done, PR #367 |
| GEN.54 | A `--star-type` secondary gets a mass that doesn't fit its type (bug) | none | done, PR #367 |
| GEN.55 | Same seed, same data: a sector's contents depend only on the seed, the version and its address (internal) | none | done, PR #3152 |
| GEN.56 | Every random draw in generation comes from the derived seeds | none | done, PR #791 |
| GEN.57 | A sector's contents depend only on the seed, the version and its address | none | done, PR #3043 |
| GEN.58 | A fingerprint of a galaxy's generated content | none | done, PR #809 |
| GEN.59 | Admin changes stored as a net difference from the generated galaxy | none | dropped (Boss, 2026-10-09 20:42Z) |
| GEN.60 | Rogue gas giants get a Jupiter-sized radius at every mass (bug) | none | done, PR #415 |
| GEN.61 | The daily merge folds pending admin changes into a new JSON file | none | dropped (Boss, 2026-10-09 20:42Z) |
| GEN.62 | Binary stars: two-word names sharing the first word, the companion's word drawn from "small" and "child" sounds, planets named for one word (bug) | none | done, PR #393 |
| GEN.63 | Planet names are unique within a sector | none | dropped: names come from IDs (GEN.67, 2026-10-07) |
| GEN.64 | A packed position ID as the name of every interstellar object and bright-sweep system | none | done, PR #406 |
| GEN.65 | A generation run fails from the web UI but not from the CLI (bug) | none | done, PR #476 |
| GEN.66 | Physics on scipy, and astropy constants and units | none | done, PR #767 |
| GEN.67 | Names from IDs for objects that have no star-derived name | none | done, PR #740 |
| GEN.68 | Research: the cheapest unique IDs for every object, unfilled sectors included | none | done, PR #623 |
| GEN.69 | A unique ID for every object, star systems and unfilled sectors included | none | done, PR #623 |
| GEN.70 | A naming key in the control database, made at galaxy creation and changeable by admin | none | done, PR #731 |
| GEN.71 | Name interstellar objects, phenomena and constellations from the codec | none | done, PR #738 |
| GEN.72 | A backfilled bright star should get a name only when its sector is generated (bug) | none | done, PR #657 |
| GEN.73 | Nebulae don't get unique names (bug) | none | folded into GEN.71 |
| GEN.74 | One point-in-space object that keeps every coordinate system in step, used by every object | none | done, PR #665 |
| GEN.75 | A nebula shape from metaballs and warped noise, as a mesh | none | done, PR #580 |
| GEN.76 | A sector with no qualifying stars never generates or is marked generated (bug) | none | done, PR #442 |
| GEN.77 | Neighborhood generation fails when its first sector is below the star threshold (bug) | none | done, PR #461 |
| GEN.78 | Some regions have a star probability of zero (bug) | none | done, PR #461 |
| GEN.79 | Bright stars only land between layers -121 and 121, so the bulge never shows (bug) | none | done, PR #461 |
| GEN.80 | Population and species generation runs on worlds without a technological civilization (bug) | none | done, PR #448 |
| GEN.81 | The console refuses runs instead of warning and doing what was asked (bug) | none | done, PR #467 |
| GEN.82 | Black holes show a Hawking temperature and luminosity of zero (bug) | none | done, PR #448 |
| GEN.83 | A planetary habitability index (PHI) | none | done, PR #1025 |
| GEN.84 | Habitability design: one score structure and reconciled thresholds | none | done, PR #803 |
| GEN.85 | Atmosphere species, partial pressures and mantle redox for every planet | none | done, PR #848 |
| GEN.86 | Stellar activity (XUV, flares) and planetary magnetic fields | none | done, PR #908 |
| GEN.87 | Surface radiation dose | none | done, PR #985 |
| GEN.88 | Hydrosphere and ocean chemistry | none | done, PR #963 |
| GEN.89 | The habitability score for every planet and moon | none | done, PR #1025 |
| GEN.90 | Refactor the planet classes around the habitability index | none | dropped (classes are derived from the finished body (cut-over epoch 121); Core 4 classifier re-read 2026-10-30, reports/core4-classifier-reread-2026-10-30.md: the clause is a fixed letter rule the derived classes make moot) |
| GEN.91 | Classes like S and V in the hot and cold zones | none | open |
| GEN.92 | Life and its highest stage follow the habitability score | none | done, PR #2833 |
| GEN.93 | Nebula conditions in planet generation | none | open |
| GEN.94 | Feasibility study: can planets form in each nebula class, and what changes | none | done, PR #823 |
| GEN.95 | Nebula conditions applied when planets and surfaces are generated | none | open |
| GEN.96 | Generation directives for a sector (an override button) | none | done, PR #915 |
| GEN.97 | Generate N random neighborhoods | none | done, PR #950 |
| GEN.98 | Bright-star backfill from the farthest generated boundary outward | none | done, PR #795 |
| GEN.99 | Nebula volume backfill with the star types the nebula needs | none | dropped (Boss, 2026-10-16, 2:51 PM EDT: the required stars are made with the nebula, GEN.241) |
| GEN.100 | Phenomena placed galaxy-wide first and kept when sectors fill | none | done, PR #793 |
| GEN.101 | Fill order: nearest sectors first along a pruned Hilbert octree curve | none | done, PR #1807 |
| GEN.102 | Investigate filling all near-zero-density void space at once | none | done, decided no pre-fill (Boss, 2026-10-16, 5:55 PM EDT) |
| GEN.103 | Research where each star type and phenomenon belongs in the galaxy's structure | none | done, PR #1953 |
| GEN.104 | A spin vector and a realistic axial tilt for every rotating object | none | done, PR #826 |
| GEN.105 | Orbital updates | none | open |
| GEN.106 | Movement thresholds and a next-update-due column | none | done, PR #802 |
| GEN.107 | The update reports how many objects moved, changed sector, or entered or left a nebula | none | done, PR #817 |
| GEN.108 | Orbital math limits: where each method breaks down and what happens there | none | done, PR #815 |
| GEN.109 | N-body influence from every object inside the largest nearby Hill sphere plus the galactic gradient, with a Hill-radius warning | none | open |
| GEN.110 | Rogue planet collisions: asteroid fields, merged giants and new stars | none | open |
| GEN.111 | Email the admin when two objects are inside each other's Hill radius | none | open |
| GEN.112 | Plan asteroid fields and belts as object systems for rendering | none | done, PR #1849 |
| GEN.113 | Analyze the anomaly docs: which anomalies to add and how | none | done, PR #1721 |
| GEN.114 | Add the chosen anomalies to the starmap | none | split into GEN.254, GEN.255, GEN.256 and GEN.257 (2026-10-15) |
| GEN.115 | The galaxy's own gravity: a smooth disk, bulge and halo potential | none | open |
| GEN.116 | Error when generating a neighbourhood centred near the galaxy edge, awaiting Boss's error text (bug) | none | closed, not reproduced (Boss 2026-10-08) |
| GEN.117 | The Galaxy Map shows bright stars only in a thin band on the galactic plane (bug) | none | done, PR #600 |
| GEN.118 | The galaxy bulge is about 40 times too light, so edge-on views show no bulge (bug) | none | done, PR #491 |
| GEN.119 | The galaxy density model has no thick disk (bug) | none | done, PR #491 |
| GEN.120 | Integrate gatedPhonemeCodec.py into the naming package and remove it from the repo root | none | done, PR #554 |
| GEN.121 | A velocity on every object, filled at generation and stored with an epoch | none | done, PR #735 |
| GEN.122 | Orbital elements for planets, moons and comets, kept in step with the state vector | none | done, PR #741 |
| GEN.123 | The projected path of a body through a sector, saved as a spline | none | done, PR #762 |
| GEN.124 | Every object knows its sector address (ring, layer, slot), recalculated whenever its position changes | none | done, PR #729 |
| GEN.125 | Stand-alone facilities store a velocity | none | done, PR #782 |
| GEN.126 | Run an orbital update as the last step of a generation run | none | done, PR #771 |
| GEN.127 | A sector generated around a backfilled bright star gives that star a planetary system (bug) | none | done, PR #805 |
| GEN.128 | Design: multi-star hierarchies and compact-object primaries | none | done, PR #1186 |
| GEN.129 | Multi-star systems of up to twelve stars | none | open |
| GEN.130 | Compact-remnant systems at sector fill: companions and pulsar planets for neutron stars and black holes | none | done, PR #2307 |
| GEN.131 | Bright-star scatter logs how many stars it added to each layer, by type | none | done, PR #805 |
| GEN.132 | Per-sector regional rates for neutron stars and black holes, extended system-wide where the research supports it | none | done, PR #807 |
| GEN.133 | Analysis of every star type's rate against its distance from the galactic core | none | done, PR #807 |
| GEN.134 | Tune the star populations to the observed star-formation profile by galactic radius | none | done, PR #1972 |
| GEN.135 | Deterministic math helpers for every stored float, with a lint and frozen constants | none | done, PR #1927 |
| GEN.136 | Fingerprint tree: a stored leaf digest per sector and a shell-band tree of digests | none | done, PR #2617 |
| GEN.137 | Placed objects move along a galactic orbit but never by velocity times time, and phenomenon_scatter has no plan time (bug) | none | done, PR #843 |
| GEN.138 | Moon `hill_radius_km` uses the star's mass, so moon spacing and the orbit slider are wrong (bug) | none | done, PR #865 |
| GEN.139 | Orbit-update thresholds: per-object epoch, path-length rule and what the 0.01 mpc applies to (GEN.106 built) | none | open |
| GEN.140 | Orbital math guards the edge-case table adds (GEN.108 built) | none | done, PR #2473 |
| GEN.141 | Faster Kepler solver (Mikkola or Markley) with brentq as fallback | none | merged into GEN.212 (Boss, 2026-10-16, 11:58 PM EDT: consolidation sweep) |
| GEN.142 | Peculiar velocity for rogue planets and asteroid fields | none | superseded by GEN.557 (folded into the wanderer body draw), 2026-10-21 |
| GEN.143 | The orbital-event admin report, the last-collisions view and a test that runs the whole collision path | none | open |
| GEN.144 | One distance for the Sun from the galactic centre across the constants, the density model and the design docs | none | done, PR #2473 |
| GEN.145 | Class S atmosphere rule: S keeps air unless the shoreline ratio is over 30 | none | dropped (classes are derived from the finished body (cut-over epoch 121); Core 4 classifier re-read 2026-10-30, reports/core4-classifier-reread-2026-10-30.md: the clause is a fixed letter rule the derived classes make moot) |
| GEN.146 | Teff-dependent habitable zone from the Kopparapu table | none | done, PR #2685 |
| GEN.147 | Classes N and Q carry a life chemical and an uncapped life timeline though they are lifeless (bug) | none | done, PR #865 |
| GEN.148 | Habitability index follow-ups from the research (GEN.84 built) | none | done, PR #1457 |
| GEN.149 | Planetary-nebula central stars: 0.5 to 0.7 Msun, 1e2 to 1e4 Lsun, up to 2e5 K | none | done, PR #1760 |
| GEN.150 | H II region radius and density from the ionizing photon rate, and IMF-based nebula hosts | none | done, PR #1442 |
| GEN.151 | Supernova remnant sizes from the density-dependent Sedov-Taylor law (GEN.10 follow-up) | none | done, PR #1795 |
| GEN.152 | Nebula cloud field is 10 to 40 times too full; lower it to the observed filling (GEN.47 rate check) | none | done, PR #1239 |
| GEN.153 | Magnetar subtype of neutron star, and an age-dependent pulsar fraction | none | done, PR #1792 |
| GEN.154 | Show the Einstein radius on compact-object pages | none | done, PR #1799 |
| GEN.155 | A nuclear-cluster object for the Sgr A* sector (optional) | none | dropped (not built; decision note, PR #1716) |
| GEN.156 | Pin astropy to CODATA 2018 and IAU 2015 constants before its first import (GEN.66 follow-up) | none | done, PR #1310 |
| GEN.157 | The `neighbor_galaxies` table and a verified data file of about 25 real galaxies | none | open |
| GEN.158 | Add a metallicity value to stars | none | done, PR #1789 |
| GEN.159 | Globular clusters: cluster table, King tables and the Milky Way catalogue | none | done, PR #1810 |
| GEN.160 | Cluster density in the sector gate, with a "cluster" population | none | done, PR #1817 |
| GEN.161 | Bright-first fill for cluster sectors | none | open |
| GEN.162 | Planet cull and blue stragglers in clusters | none | open |
| GEN.163 | Type-B pulsar planets in globular clusters (GEN.130 follow-on) | none | open |
| GEN.164 | Synthetic globular-cluster systems for generated galaxies | none | open |
| GEN.165 | A comet's time-0 scene position disagrees with its stored position in a rare random system (bug) | none | done, PR #948 |
| GEN.166 | A mass_range argument on NeutronStar and BlackHole, and intermediate-mass black holes as their own kind | none | done, PR #866 |
| GEN.167 | A lowest-mass option for the phenomenon scatter: --phenomenon-min-mass, default 20 solar masses | none | done, PR #866 |
| GEN.168 | The sector fill draws the phenomena below the scatter cut | none | done, PR #866 |
| GEN.169 | Decide the phenomenon scatter rates: regional factors and the 0.1% intermediate-mass black holes | none | done, PR #1186 |
| GEN.170 | Object ID layout: an 80-bit ID of birth sector, serial and body number, with pack, unpack, format and parse functions | none | done, PR #1017 |
| GEN.171 | The sector fill gives object IDs by generation rank | none | done, PR #1055 |
| GEN.172 | Run-time births get object IDs from the counters | none | done, PR #1075 |
| GEN.173 | Deleting a body and then adding one fails with IntegrityError 1062 on uq_planets_uid (bug) | none | done, PR #987 |
| GEN.174 | Bodies an admin adds are saved with a NULL uid (bug) | none | done, PR #987 |
| GEN.175 | Regenerating a phenomenon sets its uid to NULL (bug) | none | done, PR #960 |
| GEN.176 | A nebula or remnant is born in the sector holding the centre of the space it occupies | none | done, PR #1075 |
| GEN.177 | Planetary magnetic fields: a stagnant-lid factor | none | done, PR #1410 |
| GEN.178 | Magnetic fields: the induced field of an ocean moon | none | done, PR #1447 |
| GEN.179 | Store each sector's generation directive and attempt record with the sector | none | open |
| GEN.180 | Directives: a forced fill (met_forced) after K failed draws, and refusal of impossible requests up front | none | done, PR #3188 |
| GEN.181 | Directives: refuse impossible requests up front from the compound-Poisson tables | none | merged into GEN.180 (Boss, 2026-10-16, 11:58 PM EDT: consolidation sweep) |
| GEN.182 | Some comets' orbits do not bring them back: they are ejected into space (bug) | none | done, PR #938 |
| GEN.183 | A mass cut the user sets: preset values on a slider from 8 to 20 solar masses | none | done, PR #969 |
| GEN.184 | A luminosity floor the user sets: presets from 2500 to 4 million solar luminosities, default 3000 | none | done, PR #965 |
| GEN.185 | The scatter in five passes: mass-limit objects, mass-limit stars, brightest-star sector marks, the luminosity pass, then other phenomena | none | done, PR #953 |
| GEN.186 | Random neighborhoods: an option to keep away from filled space | none | done, PR #1088 |
| GEN.187 | Bright-star back scatter by mass: rings of 1, 2, 5 and 8 solar masses around filled space | none | done, PR #1059 |
| GEN.188 | Mass limit default is 8, and the mass slider and luminosity dropdown sit side by side on Generate and New galaxy | none | done, PR #1007 |
| GEN.189 | Gamma-ray burst and AGN ozone loss for the ozone_loss_flag (low priority) | none | done, PR #1800 |
| GEN.190 | The Phenomena table is empty after a web-generated galaxy: the phenomena scatter pass never runs (bug) | none | done, PR #1008 |
| GEN.191 | New galaxy ignores the mass limit slider: the plan step does not store the limit, so the scatter uses the default whatever the form says (bug) | none | done, PR #1008 |
| GEN.192 | Phenomena scatter log is thin and the Phenomena table page needs checking after a run (bug) | none | done, PR #1013 |
| GEN.193 | Phenomena table stays empty after the scatter: scattered unbuilt phenomena are not listed (bug) | none | done, PR #1019 |
| GEN.194 | Default mass limit 14 solar masses and default luminosity floor 9,000 solar luminosities | none | done, PR #1051 |
| GEN.195 | A separate mass limit for neutron stars and black holes, and the central black hole or quasar always created | none | done, PR #1081 |
| GEN.196 | One "Redo scatters" box on Generate: choose which scatters to redo, with new settings for each | none | done, PR #1093 |
| GEN.197 | Orbital event record: the `record_event` writer and the object-ID rules for ejection, capture, merger and split | none | open |
| GEN.198 | The gravity field evaluator: pull, potential and tidal gradient at any point from the point masses and the galaxy potential | none | open |
| GEN.199 | Per-zone gravity grid for a sector: stored in a cache and rebuilt when the orbits update | none | open |
| GEN.200 | One shared function for the "map luminosity" of a mass: the luminosity a main-sequence star of that mass would have | none | done, PR #1252 |
| GEN.201 | Fill the map luminosity in the scatter and in existing galaxies | none | done, PR #1365 |
| GEN.202 | Scattered nebulae (planetary nebulae and remnants) join the nebula index, deduplicated by object id | none | open |
| GEN.203 | Supernova-remnant and planetary-nebula scatter densities set to the catalogued counts | none | done, PR #1239 |
| GEN.204 | Nebula host chance per O and B star lowered to the catalogued one-per-few-O-stars | none | done, PR #1239 |
| GEN.205 | `planetgen plan --redo-nebula-field`: apply the new nebula rate to an existing galaxy | none | done, PR #1325 |
| GEN.206 | A planet's description claims life where the habitability index says none can live (bug) | none | done, PR #1278 |
| GEN.207 | PHI-6: split chemistry into air and water and add a Suitability dimension to the habitability index | none | done, PR #1429 |
| GEN.208 | The habitability description is one paragraph that says where the planet falls in each category and what protection it needs | none | done, PR #1433 |
| GEN.209 | Globular clusters are parent objects for their stars and are named when their sector is filled | none | open |
| GEN.210 | The galactic core is one round disc sector in the middle of the grid, 4 pc across and 4 pc tall | none | dropped (Boss, 2026-10-16, 3:00 PM EDT; superseded: the spherical grid replaced the cylinder grid and shell 0 is the centre ball) |
| GEN.211 | Constants table instead of astropy: faster start, 29 MB less memory per process | none | done, PR #1310 |
| GEN.212 | Galactic move in the orbital update as NumPy batches (3 to 4 times on the whole update) | none | open |
| GEN.213 | Mercury-like bodies never keep their magnetic field (bug) | none | done, PR #1412 |
| GEN.214 | Planet descriptions still say PHI-4 where they should say PHI-6 (bug) | none | done, PR #1472 |
| GEN.215 | Air shows Green with 0 Air when the air pressure is 0 (bug) | none | done, PR #1478 |
| GEN.216 | Habitability text is too wordy and the Suitability sentence is oddly worded (bug) | none | done, PR #1482 |
| GEN.217 | Refresh the stored habitability numbers on existing planets and moons | none | done, PR #1492 |
| GEN.218 | One habitability describer, in the short wording | none | done, PR #2206 |
| GEN.219 | The habitability levels page and tooltips are generated from the limits | none | done, PR #2206 |
| GEN.220 | "Habitable" means the index | none | done, PR #2206 |
| GEN.221 | Short plain prose for the fixed text | none | done, PR #2206 |
| GEN.222 | Generation text output for nebulae does not list the kinds of nebulae and the layers in each slab (bug) | none | done, PR #1529 |
| GEN.223 | The nightly orbital update moves every object in the galaxy, and each object is stored as a point in space in its own columns (bug) | none | done, PR #1757 |
| GEN.224 | Spherical grid module: shells, bands, wedges, caps and the core ball | none | done, PR #1549 |
| GEN.225 | Outline and bounds on shells | none | done, PR #1551 |
| GEN.226 | Addresses on shells: designation, sector ID, seeds and object ID | none | done, PR #1607 |
| GEN.227 | Generation paths on shells | none | done, PR #1607 |
| GEN.228 | Wide binaries must be stable: the two stars and the far orbit of each keep out of each other (bug) | none | done, PR #1610 |
| GEN.229 | Repair or re-roll the wide binaries that already break the stability rule | none | open |
| GEN.230 | Sector designations and object-ID birth fields cannot name bands at or beyond 2048 (bug) | none | done, PR #1757 |
| GEN.231 | Generate a whole block on the shell ladder: the generate --block mode and the block buttons come back | none | done, PR #3567 |
| GEN.232 | Cap the sector-path integration for fast bodies so one hypervelocity star cannot run for minutes (bug) | none | done, PR #1660 |
| GEN.233 | Grid reach: every cell out to 262 kpc has defined contents, and a plan-time reach setting widens the object ID to 96 bits | none | done, PR #1757 |
| GEN.234 | A measured stellar halo replaces the flat density floor, with the halo, cluster and satellite populations | none | done, PR #1763 |
| GEN.235 | One rule set for every cell: a sector profile computed from the address, and any asked sector in the reach is generated from it and saved | none | done, PR #1784 |
| GEN.236 | Kind factors in the profile: the rate per star of each rare kind times a factor from the cell, and the planet cull by metallicity | none | done, PR #1824 |
| GEN.237 | Satellite and cluster components: one table for the Magellanic Clouds, Sagittarius, dwarf galaxies and globular clusters | none | done, PR #1851 |
| GEN.238 | Wanderers: hypervelocity stars reach out to the grid, and kicked neutron stars and black holes are an ejecta population | none | done, PR #2132 |
| GEN.239 | High-velocity clouds and the Magellanic Stream as placed gas clouds (later) | none | open |
| GEN.240 | The three scatters fill the whole sphere with one algorithm: band groups, the phenomena pass on the object-first sampler, and a component pass | none | done, PR #1994 |
| GEN.241 | Every nebula gets the stars it needs as real systems made with the nebula, inside its bounds | none | done, PR #2144 |
| GEN.242 | A pre-main-sequence state for young stars, so protostars can exist as real stars | none | done, PR #2009 |
| GEN.243 | Wolf-Rayet stars and ring nebulae: a W state for the last 0.3 Myr of very massive stars, and nebula class X | none | done, PR #2178 |
| GEN.244 | Young open clusters and OB associations as components (later) | none | done, PR #2253 |
| GEN.245 | Debris disks around A and F stars (later) | none | open |
| GEN.246 | Luminous blue variables and hypergiants (later) | none | open |
| GEN.247 | Accreting white dwarf binaries: novae, cataclysmic and symbiotic variables, Type Ia progenitors (later) | none | open |
| GEN.248 | Tidal streams as placed components of the halo (later) | none | open |
| GEN.249 | Be stars, carbon stars and Mira variables as star attributes (later) | none | open |
| GEN.250 | A phenomena registry: one table for every kind of object, with a completeness test | none | done, PR #1700 |
| GEN.251 | Sector fill stages on their own random streams, and a system pass so every object that can have a system rolls for one | none | done, PR #1815 |
| GEN.252 | Brown dwarfs as stars: an L, T and Y state for companions and primaries, with the brown-dwarf desert | none | done, PR #2164 |
| GEN.253 | White-dwarf pollution flag and cool white dwarf presets | none | done, PR #2681 |
| GEN.254 | Superbubble field beside the nebula field (the Local Bubble cavity) | none | done, PR #1814 |
| GEN.255 | Drawn jets: Herbig-Haro flows, microquasar jets and pulsar jets along their spin axes | none | open |
| GEN.256 | A radiation hazard score derived from distance to hot stars, Wolf-Rayet stars, remnants and magnetars | none | open |
| GEN.257 | A Fermi-bubble overlay on the Galaxy Map | none | open |
| GEN.258 | Protostars embedded in giant molecular clouds | none | done, PR #2674 |
| GEN.259 | Intermediate-mass black hole share gives about 100,000 in the Galaxy against a handful of candidates (bug) | none | done, PR #2112 |
| GEN.260 | Scatter member columns and tables: a scatter row is a system with a barycentre, a total mass and a stored tree of members | none | done, PR #1757 |
| GEN.261 | The member-tree generator: stars, brown dwarfs, remnants and protostars in stable hierarchies up to twelve members | none | done, PR #1705 |
| GEN.262 | Companions in the scatter passes: the mass, luminosity and phenomenon passes draw members, with a rate correction so pairs are not double counted | none | done, PR #1793 |
| GEN.263 | Promote every bound compact pair to the scatter, whatever its mass | none | open |
| GEN.264 | The fill builds a system from its stored members; the system pass keeps only planets, moons and belts | none | open |
| GEN.265 | Nightly motion of a system: the barycentre moves with the system mass, members are computed from their stored phase on read | none | done, PR #1766 |
| GEN.266 | Compact-binary orbit rules: a black hole or neutron star with a star, tight compact pairs and tertiaries | none | done, PR #1771 |
| GEN.267 | A MAX_STARS setting and a --stars N directive, with a prevalence knob for five or more stars | none | done, PR #2704 |
| GEN.268 | System-type registry: a SYSTEM_TYPE_RATES dict with one entry per catalogue type, and a completeness test | none | done, PR #1702 |
| GEN.269 | Interaction states in the star model: a Roche classifier, `interaction` on every pair, member state and flags, rejuvenated age and the new kinds | none | done, PR #2481 |
| GEN.270 | Close-binary types: contact, Algol, spectroscopic and eclipsing flags, tertiary shares, merger products and Be stars | none | done, PR #2879 |
| GEN.271 | White-dwarf binary types: close double white dwarfs, post-common-envelope pairs, cataclysmic variables, AM CVn, symbiotic and barium stars | none | done, PR #2928 |
| GEN.272 | Accreting compact binaries: X-ray binaries, accretion luminosity, microquasars, spider pulsars and compact triples | none | done, PR #3055 |
| GEN.273 | Young groups and trapezia: the `group` node, offset and velocity columns, a decay age and the remnant draw | none | open |
| GEN.274 | Cluster-hosted compact binaries and runaway bookkeeping: `runaway_origin` and the compact companion left by a supernova | none | open |
| GEN.275 | Substellar and circumbinary types: brown-dwarf and planetary-mass companions, circumbinary disks and planets | none | open |
| GEN.276 | The sphere-to-key-runs cover: `sphere_runs` over the existing cell key (at most 256 merged runs) | none | done, PR #1966 |
| GEN.277 | Design note: which of the later phenomena are scattered, backfilled or stay sector-fill-only | none | done, PR #1721 |
| GEN.278 | Magnetar scatter class: magnetars drawn as scatter rows of young neutron stars | none | done, PR #1805 |
| GEN.279 | A realistic galactic core from the density and generation model (replaces the dropped one-off nuclear-cluster object GEN.155) | none | done, PR #1728 |
| GEN.280 | Galactic core density: nuclear star cluster and nuclear disc terms in the density model, replacing the Gaussian core of ADM.49 | none | done, PR #1855 |
| GEN.281 | Galactic core cell star ceiling and fill luminosity cap: a dense cell fills only stars at or above the cap, and the unresolved count and mass are computed from the address | none | done, PR #1820 |
| GEN.282 | Galactic core gas and young terms: the central molecular zone, the nebula rate and young clusters (Arches, Quintuplet, the central parsec) | none | done, PR #2255 |
| GEN.283 | Scatter limits by place: a core luminosity floor and mass limit, a halo luminosity floor (default 500, lowest 50), and the baseline everywhere else | none | done, PR #1904 |
| GEN.284 | Planet survival at the galactic core: a metallicity factor and a survival formula instead of the hard cut | none | done, PR #2640 |
| GEN.285 | Nuclear mass components in the galactic potential for the orbital update | none | open |
| GEN.286 | Horizontal-branch and AGB phases in the stellar evolution model | none | done, PR #1932 |
| GEN.287 | The stellar evolution model spreads giant-branch stars evenly in log luminosity, which probably overstates the bright giants by about a factor of two against real halo counts (bug) | none | done, PR #1932 |
| GEN.288 | Rename the PHI-6 equipment tier "Sealed suit" to "Environment suit" everywhere | none | done, PR #1779 |
| GEN.289 | Planets and moons kept as the universal point-in-space object, with point mass and movement vector stored, so perturbations have somewhere to land | none | done, PR #2108 |
| GEN.290 | An encounter search that cannot miss a pass at large steps (100 years per second) | none | open |
| GEN.291 | A bound planet becomes a rogue planet when an encounter unbinds it: the unbind and capture test and the event log | none | open |
| GEN.292 | Run the orbital update for a stated span in one go from the command line (elapsed years) | none | open |
| GEN.293 | Long-step rules for hypervelocity stars that leave the sector grid and for parabolic comets | none | open |
| GEN.294 | Every object in the galaxy is held as one SpatialPosition3D with its position, point mass and movement vector stored as columns | none | done, PR #1952 |
| GEN.295 | Nebulae as point groups whose shape follows the orbital update (design and research) | none | open |
| GEN.296 | The uncharted-star information panel shows the star's metallicity [Fe/H], which is already stored on bright_stars | none | open |
| GEN.297 | Calibrate the cluster generator against the Harris file: its distributions guide the laws, and no real cluster is placed in any galaxy | none | open |
| GEN.298 | Remnants read the superbubble: the remnant size at build and the remnant factor use the bubble a supernova sits in | none | done, PR #2132 |
| GEN.299 | The Quantity class (the unit class): one value, every unit version always available, all conversions and the shared ladders | none | done, PR #2490 |
| GEN.300 | Every quantity in the code moves onto the Quantity class, and hand conversions are removed | none | open |
| GEN.301 | Headline of 8.1: the new generate scatters, accurate galactic generation and scattering including the halo | none | open |
| GEN.302 | A region object: a sphere, box, shell, band, wedge or set of sectors, and the systems and phenomena inside it | none | open |
| GEN.303 | Per-cell system cap as a configuration option: `generation.cell_cap`, default 5,000, lowest 500, highest 50,000 | none | open |
| GEN.304 | The canonical fingerprint encoder: floats to 9 significant digits, SHA-256, structural and float digests | none | done, PR #2578 |
| GEN.305 | Store a leaf digest at save time and build the tree: roots and region digests | none | done, PR #2617 |
| GEN.306 | Population seeds keyed on the object, not on row ids | none | done, PR #1942 |
| GEN.307 | The object-ID collision number keyed on the lower address | none | done, PR #2595 |
| GEN.308 | A block-first sampler for the backfill: a Poisson candidate count from the block's density bound, then thinning | none | done, PR #2070 |
| GEN.309 | The backfill picks the draw path per block and gives the same stars wherever it lands | none | done, PR #2934 |
| GEN.310 | Fill a whole block of the shell ladder: planetgen generate --block and the queue job | none | done, PR #3159 |
| GEN.311 | Built systems store and load their orbit tree: insert_star_system, load_star_system and serialization read and write orbit_nodes | none | done, PR #2481 |
| GEN.312 | Naming for trees: path codes (Aa, Ab, Ba) with branch words at wide splits, and planet stems such as Pikkita Ba I | none | open |
| GEN.313 | Orbital update of built systems: advance_orbital_phases for nodes and the GEN.109 hooks | none | open |
| GEN.314 | Planet zones in a hierarchy and flux sums: S-type and P-type windows per node (Holman-Wiegert) and the flux a planet receives from every member | none | done, PR #2125 |
| GEN.315 | A fill builds a system from its stored members and nodes; the system pass keeps only planets, moons, belts and disks | none | open |
| GEN.316 | Stages D and G draw a member tree for every object the fill creates, and the two-star classes are retired | none | open |
| GEN.317 | Pulsar planets and the radiation environment of neutron stars | none | done, PR #2265 |
| GEN.318 | White-dwarf and planetary-nebula-star rows of the system pass: engulfment of inner planets and the planetary nebula's circumbinary disk | none | done, PR #2307 |
| GEN.319 | Nebula rates read the profile: young and old shares, the cell's gas and the Type Ia share as profile fields | none | done, PR #2003 |
| GEN.320 | build_nebula_with_contents: every hosted nebula gets its stars as real systems inside its bounds (classes C, G, F, planetary, remnants, P and Q) | none | done, PR #2014 |
| GEN.321 | Emission regions D and E: the ionizing members of a young cluster in a Plummer sphere, cut by sector cube | none | done, PR #2144 |
| GEN.322 | Band groups and the halo term: the bright-star and mass passes cover the 262 kpc reach | none | done, PR #1991 |
| GEN.323 | The phenomena pass runs on the object-first sampler | none | done, PR #1992 |
| GEN.324 | A component pass: clusters, the Clouds and Sagittarius place their stars, and the fill adds the component's density to lambda | none | done, PR #1993 |
| GEN.325 | Planet cull in dense components by the survival formula and the metallicity factor | none | done, PR #2624 |
| GEN.326 | Blue stragglers in clusters: N_BSS scales as the core mass to the 0.4 | none | open |
| GEN.327 | The group node and the remnant draw: offset and velocity columns, a decay age, and the lightest member ejected | none | done, PR #3217 |
| GEN.328 | Groups at fill in star-forming cores and OB regions | none | done, PR #3920 |
| GEN.329 | Runaway bookkeeping: runaway_origin and the compact companion left by a supernova | none | done, PR #2634 |
| GEN.330 | Cluster-hosted compact binaries: the cluster rates replace the field fractions | none | open |
| GEN.331 | Substellar companions: brown-dwarf close (the desert) and wide companions, brown-dwarf pairs and planetary-mass companions | none | open |
| GEN.332 | Circumbinary and S-type planets, circumbinary disks and rogue planet pairs at fill | none | open |
| GEN.333 | Scatter limits by place: the tier chosen per band by the dominant lambda term, applied in the passes, the backfill and the fill | none | done, PR #1904 |
| GEN.334 | Pair closest-approach library: sweep, pair windows, window cap and chunking | none | done, PR #2473 |
| GEN.335 | One definition of the stored point columns, and a test that lists every object table | none | done, PR #1883 |
| GEN.336 | Compact and free objects hold a velocity vector and a mass: stars, black holes, neutron stars, rogue planets, interstellar comets, asteroid fields, nebulae, remnants | none | done, PR #1952 |
| GEN.337 | Epoch, vector and mass for the rest: belts, quasars, comets, facilities, scatter points and members | none | done, PR #1952 |
| GEN.338 | The galactic potential as one pure module: Hernquist bulge, Miyamoto-Nagai disk and NFW halo, with settings and the rotation-curve tests | none | done, PR #4526 |
| GEN.339 | Generated and stored galactic velocities follow the potential: the old curve is replaced, dispersions are seeded, stored velocities are reseeded | none | open |
| GEN.340 | One star stepper: kick-drift-kick in the potential, a step per object, on arrays | none | open |
| GEN.341 | Influence-set search: every object inside the largest nearby Hill sphere, found by sector bounds, with cover lists and a brute-force check | none | done, PR #2473 |
| GEN.342 | The acceleration kernel: the pull of the influence set plus the galactic gradient, as one function, and the velocity update from it | none | done, PR #4530 |
| GEN.343 | The perturbed pass for a system: host, siblings, moons and outside perturbers integrated through an encounter | none | open |
| GEN.344 | Hill-radius warnings by the gravitational criterion, to the debug and activity logs | none | open |
| GEN.345 | Collision outcomes as a pure function: regime table, giant merger, fusion check, conservation | none | done, PR #2473 |
| GEN.346 | Applying a collision to the galaxy: new asteroid field, merged giant, merger-remnant star, retired bodies, events | none | open |
| GEN.347 | The last-collisions admin view and the test that runs the whole collision path | none | open |
| GEN.348 | Describer, CLI query lines and class page ranges take Quantity | none | done, PR #2539 |
| GEN.349 | API output and stored-column readers take Quantity; the grep test fails on a hand-written conversion factor; the old helpers are removed | none | open |
| GEN.350 | Class R, sub-Neptune, with the envelope-survival rule and the Hycean flag (class T stays: it is the sub-Neptune row of the classifier) | none | dropped (revoked by Boss's full-build sweep 8:16 PM EDT 10-29: new class letters R, U, W, X, Y and Z conflict with conditions-based classes (10-23) and 'all 20 letters stay' (10-28); the part that stands is GEN.963) |
| GEN.351 | Class U, icy world | none | dropped (revoked by Boss's full-build sweep 8:16 PM EDT 10-29: new class letters R, U, W, X, Y and Z conflict with conditions-based classes (10-23) and 'all 20 letters stay' (10-28); the part that stands is GEN.963) |
| GEN.352 | Classes W and X, subsurface ocean bodies | none | dropped (revoked by Boss's full-build sweep 8:16 PM EDT 10-29: new class letters R, U, W, X, Y and Z conflict with conditions-based classes (10-23) and 'all 20 letters stay' (10-28); the part that stands is GEN.963) |
| GEN.353 | Class Y, Titan-like world | none | dropped (revoked by Boss's full-build sweep 8:16 PM EDT 10-29: new class letters R, U, W, X, Y and Z conflict with conditions-based classes (10-23) and 'all 20 letters stay' (10-28); the part that stands is GEN.963) |
| GEN.354 | Class Z, lifeless temperate world | none | dropped (revoked by Boss's full-build sweep 8:16 PM EDT 10-29: new class letters R, U, W, X, Y and Z conflict with conditions-based classes (10-23) and 'all 20 letters stay' (10-28); the part that stands is GEN.963) |
| GEN.355 | Atmosphere retention module with its test table | none | done, PR #4132 |
| GEN.356 | Class N re-scoped as the hot, lifeless world; its life chemical and hidden timeline removed | none | dropped (classes are derived from the finished body (cut-over epoch 121); Core 4 classifier re-read 2026-10-30, reports/core4-classifier-reread-2026-10-30.md: the clause is a fixed letter rule the derived classes make moot) |
| GEN.357 | Generation knows the cloud it sits in: surrounding_cloud extended, containers found early, StarSystem.apply_cloud | none | done, PR #3367 |
| GEN.358 | GEN.94's rules and the four exotic-environment additions reach planet formation | none | open |
| GEN.359 | Sector fill links cluster stars to their cluster and names the cluster at its first fill | none | open |
| GEN.360 | Asteroid fields and belts render from their uid: one seeded generator in a shared asteroidfield.js for dust, rocks and heroes | none | open |
| GEN.361 | The asteroid field page has a 3D view | none | open |
| GEN.362 | Rocks are instanced meshes in chunks with three levels of detail, within a frame budget | none | open |
| GEN.363 | Belt bodies move on the GPU: the vertex shader solves Kepler from per-instance elements | none | open |
| GEN.364 | Adaptive per-cell system cap: above the cap the luminosity floor rises per cell until what is resolved fits, brightest stars first, with no step at the ceiling | none | done, PR #1868 |
| GEN.365 | Hot white dwarf and planetary-nebula central star radius disagrees with the radius that its luminosity and temperature imply (bug) | none | done, PR #4385 |
| GEN.366 | Check the reference satellite table against McConnachie 2012: it only tests that the derived population looks plausible | none | done, PR #4357 |
| GEN.367 | The Magellanic Clouds' gas is not in `kind_factors` and `nebula_field`: the Clouds hold stars but no nebulae | none | open |
| GEN.368 | Plan the appearance rates of everything the Milky Way supplies as fixed content: the count, place and size laws, with sources, so a random galaxy draws them from its seed | none | done, PR #1960 |
| GEN.369 | Satellites drawn from the galaxy seed: count, distance, mass, size, profile, [Fe/H] and age are calculated from stated laws, with no Milky Way option | none | done, PR #2638 |
| GEN.370 | Hill spheres, keep-out spheres and the default tidal distance use the Milky Way's mass and Sol's distance for every galaxy: take the galaxy's own enclosed mass and the cell's own distance (bug) | none | done, PR #2460 |
| GEN.371 | Feasibility study: generate the satellite galaxies through the scatter, making only the region each one occupies | none | open |
| GEN.372 | Numeric precision audit at galactic scale: every float coordinate and sum, how far it gets from its origin, where the first error shows, and what needs a wider or split representation | none | done, PR #1959 |
| GEN.373 | Massive multiples sit in young regions: set `prevalence` from the cell's young share, with a derived cap on the rise | none | dropped (Boss, 2026-10-17, 10:07 AM EDT, his card answer as Core 4 reported it) |
| GEN.374 | Orbit phase formulas reduce first, and the JavaScript Kepler solver uses tolerance 1e-12: Python, SQL and JavaScript agree | none | done, PR #2473 |
| GEN.375 | Field kernel subtracts a local origin before its loop in `physics/sector_path.py` | none | open |
| GEN.376 | Region test uses the haversine form for angles under 1e-3 rad in `galaxy/region.py` | none | open |
| GEN.377 | Scatter limits by place for the neutron star and black hole scatter below the cut: the core and halo tiers apply to it, not one global mass cut | none | open |
| GEN.378 | The galaxy type families: the checked design with its scaling relation table and the draw order | none | done, PR #2368 |
| GEN.379 | The density model as a list of components, each with an amplitude that may be 0, and a mass for each component | none | done, PR #2531 |
| GEN.380 | The galaxy type draw: stellar mass from a Schechter function, then the family from the mass-dependent fractions, then the halo mass | none | open |
| GEN.381 | A Sérsic spheroid as a density component: index, effective radius, axis ratios and orientation, with a mass | none | open |
| GEN.382 | The spiral family: scale length and height, bulge-to-total, bar, arm count and pitch angle from mass and Hubble stage | none | open |
| GEN.383 | The elliptical family: Sérsic index and size from mass, intrinsic axis ratios, no disk, arms or bar, an old population | none | open |
| GEN.384 | The lenticular family: a smooth disk and a large bulge, no arms, old stars and almost no gas | none | open |
| GEN.385 | The irregular and dwarf families: a clumpy exponential or Sérsic body with young stars, and an old dwarf spheroid | none | open |
| GEN.386 | Population mix, star-formation history and gas fraction by family: the young share that the nebula and profile rules read | none | open |
| GEN.387 | Globular clusters, the nuclear cluster and the central black hole follow the family's spheroid and halo | none | open |
| GEN.388 | The home position by family: one rule of the shape, with the spiral rule as the disk case | none | dropped (Boss, 2026-10-16, 9:43 PM EDT: the galaxy stores no home position; the map starts where sectors are generated) |
| GEN.389 | Galaxy creation draws the family from the seed and stores it, with the family in the log, the API and the version check | none | open |
| GEN.390 | Headline of 9.4: generate a galaxy of any known type, with spiral, lenticular, elliptical, irregular and dwarf families from one set of rules | none | open |
| GEN.391 | Superbubbles shape the stars: the design of the bubble field as a density modulation, with its bounds and touch points checked | none | open |
| GEN.392 | A bubble field handle on the galaxy shape: seed-carried bubbles with a box cut, a certified maximum factor and a bubble-free density read | none | open |
| GEN.393 | The density terms apply the bubble factors at one call site per form, with a golden test that no bubble means identical density | none | open |
| GEN.394 | The scatter bounds and the phenomena pass respect the bubble factor: a certified bound in band_majorants and per-position thinning | none | open |
| GEN.395 | Sector build follows the bubbles: the fill count from the cell mean, star positions weighted by the factor, shown in the profile and the region query | none | open |
| GEN.396 | Revisit the Yerkes class III luminosity range: 50 to 1000 solar luminosities sits above many dim red giants | none | open |
| GEN.397 | `planetgen plan` stage 2 prints the progress weight as the number of phenomena to place (bug) | none | done, PR #2124 |
| GEN.398 | The massive-star and bright-star passes place more stars than their log lines expect: 70 and 12 percent over (bug) | none | done, PR #3162 |
| GEN.399 | Research the halo mass function and compare it with the halo mass approximation GEN.368 uses | none | done, report gen399-halo-mass-function-2026-10-31.md |
| GEN.400 | The galaxy's mass chain: stellar mass from the Schechter function, the halo mass by the approximation, the virial radius and the concentration | none | done, PR #2011 |
| GEN.401 | The shape of the stars from the stellar mass: disk length and height, thick disk, bulge, bar and arms replace the generate defaults | none | done, PR #2051 |
| GEN.402 | The central black hole mass and the nuclear cluster follow the bulge: the M-bulge and M-sigma relations replace the fixed ranges | none | done, PR #2224 |
| GEN.403 | The globular cluster count follows the halo mass: the system mass from the halo mass replaces the fixed count of 157 (bug) | none | done, PR #2638 |
| GEN.404 | The Sun radii are removed: every use of a home or Sun radius becomes a calculation at the object's own position | none | done, PR #2270 |
| GEN.405 | Metallicity and rates by position from the galaxy's own model: the mass-metallicity relation and the local rates replace the fixed gradient and the quoted-at-the-Sun values | none | open |
| GEN.406 | Research the population ages and velocity dispersions by component, the thick disk scaling and the stellar halo fraction so they become laws | none | done, PR #2022 |
| GEN.407 | Correct the horizontal-branch and red giant tables to the checked sources: the RGB time law, the evolutionary flux, the HB colour row, the extreme HB luminosity and the tip | none | done, PR #1985 |
| GEN.408 | Research: verify the atmosphere sources and say what is missing | none | done, PR #1951 |
| GEN.409 | The atmosphere colour design: the checked plan for the gas table, the mixture, the stored values, ozone, aerosols and the colour model | none | done, PR #2063 |
| GEN.410 | One table of the gases and compounds that can be in an atmosphere, with every compound the descriptions mention | none | done, PR #2508 |
| GEN.411 | The atmosphere's molar density comes from its per-gas partial pressures, so mixture, pressure and density always agree | none | done, PR #3311 |
| GEN.412 | Generation records the albedo and the greenhouse multiplier it drew, and the temperature is exactly what they give | none | done, PR #2812 |
| GEN.413 | Ozone in the mixture: the column and partial pressure from the oxygen, the star's UV and the Chapman balance | none | done, PR #2761 |
| GEN.414 | Haze, dust and cloud: derived or drawn optical depths for every planet and moon | none | done, PR #2856 |
| GEN.415 | Every compound the planet descriptions mention comes from the gas table, and the atmosphere text is written from the stored mixture | none | done, PR #3472 |
| GEN.416 | Blackbody colour of a star: Planck to the CIE curves to sRGB, as the first part of the photometry library | none | done, PR #2566 |
| GEN.417 | Rayleigh optical depth of the mixture: per-gas cross sections from the table and the air column from pressure and gravity | none | done, PR #2724 |
| GEN.418 | Atmosphere colour of a solid world: the absolute colour of the air with ozone, haze, dust, cloud and the ground's bounce | none | done, PR #2768 |
| GEN.419 | Atmosphere colour of gas and ice giants, read at a reference level in the cloud tops | none | done, PR #2773 |
| GEN.420 | The descriptions show the atmosphere colour: a swatch and a word in the planet and moon facts, in every output | none | done, PR #2791 |
| GEN.421 | Headline of 9.0: what every body and atmosphere is made of, and the atmosphere colour that follows from it | none | open |
| GEN.422 | Stop the recalculation redraw and order the steps so pressure, mixture and temperature agree after the cold-side freeze-out | none | done, PR #3311 |
| GEN.423 | Atmosphere density, scale height and the greenhouse multiplier follow from the mixture | none | done, PR #3474 |
| GEN.424 | The balance gas is a named gas with its own molar mass, and habitability reads and writes the mean molar mass of the stored mixture | none | done, PR #3110 |
| GEN.425 | Research: verify the mass-radius and equation-of-state sources and say what is missing | none | done, PR #1957 |
| GEN.426 | The bulk composition design: the fractions, the equation of state, the draws and the derived density, mass and gravity | none | done, PR #2063 |
| GEN.427 | One table of the bulk materials with their equation-of-state constants and the source of each | none | done, PR #2580 |
| GEN.428 | Planets and moons draw their bulk composition fractions instead of carrying a fixed text | none | done, PR #2751 |
| GEN.429 | Bulk density, mass, radius and gravity are derived from the composition fractions with the mass-radius relation | none | done, PR #2751 |
| GEN.430 | The planet and moon descriptions state the composition from the stored fractions | none | done, PR #3378 |
| GEN.431 | Research: what is in the atmospheres of every kind of world, and the numbers each compound needs | none | done, PR #1970 |
| GEN.432 | Research: what rocky, icy and giant bodies and their small relatives are made of, and the numbers each material needs | none | done, PR #1972 |
| GEN.433 | Research: what stars and stellar remnants are made of, and the numbers each needs | none | done, PR #1972 |
| GEN.434 | Research: radiation through atmospheres and at the ground, and the numbers it needs | none | done, PR #1972 |
| GEN.435 | The stellar composition design: the abundance relations, what is stored and what is derived, and the remnants | none | done, PR #2063 |
| GEN.436 | The radiation design: atmospheric transmission and greenhouse by band, shielding by column, and the ground's own dose | none | done, PR #2063 |
| GEN.437 | The composition model: how the element, compound, material and body tables link, and what is stored against derived | none | done, PR #2063 |
| GEN.438 | One table of the elements and their isotopes: atomic data and the radioactive isotopes' half-lives and decay energies | none | done, PR #2184 |
| GEN.439 | The stellar composition relations: abundances by population, age and [Fe/H], and the remnants' layers | none | done, PR #2208 |
| GEN.440 | Star composition is derived from the stored [Fe/H], age and population, for every star, remnant and brown dwarf | none | done, PR #2208 |
| GEN.441 | Surface and crust composition: the minerals, ices and ocean on the ground, from the bulk fractions, the temperature and the air | none | done, PR #3335 |
| GEN.442 | Asteroids, comets, rings and dwarf planets draw their composition from the formation distance and the snow line | none | done, PR #2956 |
| GEN.443 | Atmospheric radiation: transmission by band and the greenhouse effect from the stored mixture, the aerosols and the cloud | none | done, PR #2910 |
| GEN.444 | Surface radiation: the ground-level dose from the air's shielding and the ground's own radiogenic content | none | open |
| GEN.445 | Habitability, toxicity and the UV index read the derived radiation instead of the proxies | none | done, PR #FOLDED |
| GEN.446 | The descriptions state the composition of stars, surfaces and small bodies, and the radiation at the surface | none | open |
| GEN.447 | The wide-pair encounter limit uses the relative-velocity dispersion and the mass-weighted perturber mass | none | done, PR #2100, #2101, #2102 |
| GEN.448 | Massive stars get binaries between 0.25 and 50 AU, so the metallicity trend reaches the separations where it was measured | none | done, PR #2100, #2101, #2102 |
| GEN.449 | The companion table for mid and early B stars is raised to the observed frequencies | none | done, PR #2100, #2101, #2102 |
| GEN.450 | The star-formation profile's radial slope: reconcile the Schmidt index with the observed gradient, with the choice named | none | done, PR #2776 |
| GEN.451 | The mantle redox mean follows the meteorite building blocks: a young rocky planet starts reduced, not at the modern Earth's 3.5 above the iron-wuestite buffer | none | done, PR #2820 |
| GEN.452 | The condensation table is replaced by the checked saturation curves: the hydrogen sulphide row is a copy of the sulphur dioxide row, and carbon dioxide, sulphur dioxide and the cold gases are 2 to 20 K off (bug) | none | done, PR #3041 |
| GEN.453 | The radiation constants are corrected: the dipole cutoff scales with moment over radius squared, the dose curve against the column and the field factor follow the measured and simulated tables (bug) | none | done, PR #2726 |
| GEN.454 | Giant planets, ice giants and brown dwarfs draw their atmosphere mixture from the equilibrium chemistry sequence: by temperature, pressure, metallicity and quench | none | done, PR #2820 |
| GEN.455 | Airless bodies and comets draw an exosphere or coma from their surface composition, temperature and the star's wind: the sodium and potassium tails and the production rates | none | open |
| GEN.456 | Lava worlds and hot rocky planets draw a vapour atmosphere over the magma from the surface composition and temperature: sodium, potassium, silicon monoxide, magnesium, iron and the oxides | none | done, PR #2820 |
| GEN.457 | The small bodies and moons composition design: the taxon and comet draws, compaction, the moon origin types and the ice fraction | none | done, PR #2374 |
| GEN.458 | Each system draws its snow line, 0.2 dex around the 170 K blackbody position, and stores it once; the fixed 2.7 AU constant goes | none | done, PR #2655 |
| GEN.459 | One radiogenic module: heat from the isotope-wise decay of potassium, thorium, uranium and aluminium-26, drawn from the host star's abundances, replaces the bulk-silicate values that are 37 to 50 percent too high and the separate ground-dose draw (bug) | none | done, PR #2868 |
| GEN.460 | Giant planets and brown dwarfs get their radius and heavy-element content from the checked laws: the heavy-element mass law, the Fortney and Baraffe tables with inflation, and the burning limits shifted by metallicity and helium (bug) | none | done, PR #2912 |
| GEN.461 | The hydrosphere reads the stored water mass fraction: the compressed water-layer depth, the base pressure, the high-pressure ice and the liquidus are recomputed, and the class label follows from the water share (bug) | none | done, PR #2993 |
| GEN.462 | Ground albedo by band: the area-weighted colour of the ground, ice and ocean cover with space weathering, handed to the atmosphere colour once | none | done, PR #3920 |
| GEN.463 | The neutron star's radius is derived from sourced mass-radius relations instead of the fixed range of 10 to 13 km (bug) | none | done, PR #2112 |
| GEN.464 | The disk solids law reads the star's metal mass fraction Z over Z0, not 10 to the power of [Fe/H], which is off by 1.7 to 3 times in metal-poor alpha-enhanced stars | none | done, PR #3493 |
| GEN.465 | A second extraction pass over the sources the body and atmosphere research could not reach and that are reachable now | none | done, PR #2645 |
| GEN.466 | Moons draw their origin type, ice fraction and composition from the parent, the circumplanetary ice line and the age | none | done, PR #3282 |
| GEN.467 | The main-sequence lifetime is derived from sourced stellar relations instead of the fixed 10 Gyr for one solar mass and a power law of -2.5 (bug) | none | done, PR #2112 |
| GEN.468 | Settle the water share beyond the snow line: the Marboeuf-derived 0.14 to 0.39 of the solids against the convention of about 0.5 | none | done, PR #2644 |
| GEN.469 | The disk's star-formation history from the gas law: the local B-star share is 0.2 percent against the observed 0.04 percent | none | done, PR #2139 |
| GEN.470 | One saturation-pressure function p_sat(species, T) for the gases and ices: the Fray and Schmitt polynomials, water by Feistel and Wagner, and its tests | none | done, PR #2590 |
| GEN.471 | Clathrate hydrates and eutectic brines for icy crusts and ocean moons: the dissociation pressures and the liquid range below 273 K | none | open |
| GEN.472 | One gated possibility-table format and loader: every composition category is a table of possibilities with a class gate and a weight, drawn by the normalised weights of the rows that pass | none | done, PR #2061 |
| GEN.473 | Replace the asteroid taxon shares and the zone kernel, measured in the Solar System belt, by a derived or general distribution | none | open |
| GEN.474 | Replace the Earth rock-type library, measured on Earth, by ground types derived from the bulk composition and the melting of the mantle | none | open |
| GEN.475 | Meteorite groups become outputs and not inputs: regions of composition space reachable from the drawn bulk, with the measured groups kept as calibration | none | open |
| GEN.476 | Aluminium-26 and iron-60 follow from the star's formation environment, not from the Solar System's own value and range | none | open |
| GEN.477 | The asteroid belt and field composition is drawn from the star, the distance and the belt's mass, not four names drawn uniformly from a list of 84 (bug) | none | open |
| GEN.478 | The atmosphere possibility tables: the 16 categories and 104 rows of part 09, loaded through the gated table loader, with every row that has no sourced weight flagged and drawn from its bracket | none | done, PR #2587 |
| GEN.479 | The class atmosphere mixes contradict the checked chemistry: classes C, D and S have no mix, free oxygen sits in seven classes with no life rule, class L is 60 percent argon, E and F hold ammonia and methane at 340 to 400 K, N reads reducing and is oxidised, and I, J and T share one mix in every zone (bug) | none | done, PR #3729 |
| GEN.480 | Source the occurrence weights of bare rock, hydrogen-helium envelope and water-rich worlds and of the envelope metal fraction, by radius, mass and insolation (AR1) | none | done, PR #3853 |
| GEN.481 | Source the volatile delivery against formation distance: the water mass fraction by a_form and the snow-line law over time (AR2) | none | done, PR #3220 |
| GEN.482 | Source the magma-ocean fO2 prior against mass and the oxidation law of the later mantle (AR3) | none | done, PR #3220 |
| GEN.483 | Source the cosmic-shoreline constants and the pressure caps by retention regime (AR4) | none | open |
| GEN.484 | Source the giant planets' C/O distribution and the fraction of carbon-rich hosts (AR5) | none | done, PR #3840 |
| GEN.485 | Source the hot-Jupiter thermochemistry above 2200 K, a Kzz distribution and the edges of the 1500 K and 2200 K regimes (AR6) | none | done, PR #3840 |
| GEN.486 | Source the lava-world vapour regimes and the fraction of lava worlds that keep volatiles (AR7) | none | done, PR #3840 |
| GEN.487 | Source the evolved and frost atmospheres: cold vapour pressures, the habitable-zone coefficients, the occurrence of Titan-, Pluto- and Venus-like bodies and the abiotic O2 frequency (AR8) | none | done, PR #3848 |
| GEN.488 | Source the comet coma classes and their volatile ratios (AR9) | none | done, PR #3840 |
| GEN.489 | Source the white-dwarf and pulsar planet options (AR10) | none | done, PR #3853 |
| GEN.490 | Source the trace-gas presence weights (AR11) | none | done, PR #3853 |
| GEN.491 | Source the cloud and haze windows and yields (AR12) | none | done, PR #3853 |
| GEN.492 | Source the ground radioisotope abundances by host and the potassium depletion by body class (AR13) | none | open |
| GEN.493 | Source the initial-rotation distribution and the saturation time against stellar mass for the host-star UV regime (AR14) | none | open |
| GEN.494 | Source the brown-dwarf type against temperature and the T/Y boundary near 500 to 700 K (AR15) | none | open |
| GEN.495 | The comet population table: species ranges and log-normal scatters replace the draw centred on one comet, with a carbon-chain-depleted flag | none | done, PR #2596 |
| GEN.496 | The photolysis lifetime table: species lifetimes at 1 AU for the quiet Sun, the active-Sun factors and the species-specific stellar-spectrum ratios | none | done, PR #2585 |
| GEN.497 | The cold-end vapour pressures: the laboratory-corrected N2, CO, CH4, O2 and NH3 below 50 K and the Span and Wagner carbon dioxide sublimation curve | none | done, PR #2671 |
| GEN.498 | One oxygen-fugacity scale for the generator: the reference fixed, delta-IW stored on it, and the buffers converted | none | done, PR #2583 |
| GEN.499 | The lava-world vapour species and mean molar mass from temperature and composition, replacing the heavy-oxide range of 40 to 60 g/mol | none | done, PR #2820 |
| GEN.500 | The licence audit: the licence column filled for every data source, and the permission list for NIST, JPL, Leiden, Bremen and Karkoschka | none | dropped (Boss, 2026-10-17, 2:35 PM EDT, typed in the project chat: every acquired paper and dataset is usable, so no licence audit) |
| GEN.501 | Source the activity factor of the NH3 photolysis rate from the Huebner and Mukherjee rate tables | none | done, PR #2824 |
| GEN.502 | The giant branch finishes the stellar tables check: the red giant time law as a function of core mass, the tip by metallicity and the branch duration against the observed horizontal-branch ratio | none | done, PR #2112 |
| GEN.503 | The shared tuning harness: a target table, a simulate-and-compare function and a loader that refuses a tunable row without a target | none | done, PR #2753 |
| GEN.504 | The Kepler transit-plus-recovery selection model | none | done, PR #2756 |
| GEN.505 | The mass-radius sample selection model | none | done, PR #2756 |
| GEN.506 | The JWST atmosphere target-selection model | none | done, PR #2756 |
| GEN.507 | The comet admission model by heliocentric distance | none | done, PR #2756 |
| GEN.508 | The white-dwarf pollution detection and settling model | none | done, PR #2756 |
| GEN.509 | Tuned draw builder B1: the radius valley and the envelope | none | done, PR #3025 |
| GEN.510 | Tuned draw builder B2: delivery and water worlds | none | open |
| GEN.511 | Tuned draw builder B3: retention and lava worlds | none | done, PR #2756 |
| GEN.512 | Tuned draw builder B4: the iron branch | none | done, PR #3250 |
| GEN.513 | Tuned draw builder B5: redox and white-dwarf debris | none | done, PR #2756 |
| GEN.514 | Tuned draw builder B6: radioisotopes | none | done, PR #2756 |
| GEN.515 | Tuned draw builder B7: ground and petrology | none | open |
| GEN.516 | Tuned draw builder B8: giant planets | none | done, PR #2756 |
| GEN.517 | Tuned draw builder B9: brown dwarfs | none | done, PR #2756 |
| GEN.518 | Tuned draw builder B10: the chemistry module | none | done, PR #2765 |
| GEN.519 | Tuned draw builder B11: aerosol amplitude | none | open |
| GEN.520 | Tuned draw builder B12: the comet sampler | none | done, PR #2760 |
| GEN.521 | Tuned draw builder B13: belts | none | open |
| GEN.522 | Tuned draw builder B14: moons and Kuiper belt objects | none | open |
| GEN.523 | Tuned draw builder B15: host-star activity | none | open |
| GEN.524 | Tuned draw builder B16: halogens and trace gases | none | open |
| GEN.525 | Tuned draw builder B17: the design-flag registry | none | done, PR #2756 |
| GEN.526 | Satellites with gas get young and intermediate populations derived from their age history and gas, so the Clouds hold massive stars | none | done, PR #2043 |
| GEN.527 | The second star of a wide pair in a saved runaway system loads back 0.4 m/s off the system's velocity: find the cause instead of loosening the test (bug) | none | done, PR #2075 |
| GEN.528 | The rest of GEN.319: a Magellanic cell shows every class, and the scatter's planetary nebulae are weighted by the old share | none | done, PR #2077 |
| GEN.529 | The pre-main-sequence track takes its Hayashi temperature, birthline radius and energy factor from a sourced table, and the state reaches stars still contracting after 10 Myr | none | done, PR #2112 |
| GEN.530 | A pre-placed star's hosted nebula shows from its scatter row on the maps and pages | none | done, PR #2388 |
| GEN.531 | Boss's research requests from the galaxy laws research (R1 to R15): the observed stellar halo, the dwarf range, the thick disk fits and the cluster system | none | open |
| GEN.532 | The galaxy shape's thick disk, scale heights and stellar halo come from the research laws instead of the tuning.py numbers (the GEN.401 follow-up) | none | done, PR #2186 |
| GEN.533 | The galaxy shape draws bar presence, axis ratios and angle from the seed, rebuilt each time and not stored in a column | none | done, PR #2186 |
| GEN.534 | A satellite's young stars reach the system-type tables: the young-star system types name `satellite_young` | none | open |
| GEN.535 | Check the Magellanic gas fraction and supernova-remnant count against the observed Clouds, and the dwarfs that still come out with gas (Leo I) | none | done, report gen535 (GEN.1191, GEN.1192, GEN.1193, TEST.508) |
| GEN.536 | Research the velocity dispersion of molecular gas (CO) in nearby spirals to set the young stars' scale height | none | done, report gen536 (GEN.1181, GEN.1182) |
| GEN.537 | Design the plan configuration: the boundaries of the sector plan, how far scatters go, which sector addresses are valid, and the two switches halo and satellite galaxies (both default off) | none | open |
| GEN.538 | The plan configuration is the single source of the plan's boundaries, scatter reach and valid sector addresses: scatter and sector generation read it and never decide them | none | open |
| GEN.539 | The halo command: its own command that switches the halo on or off in the plan configuration (default off) and sets how far the halo reaches | none | open |
| GEN.540 | The satellite galaxies command: its own command that switches the satellite galaxies on or off in the plan configuration (default off), separate from the halo | none | open |
| GEN.541 | The wide-pair perturber mass: black holes carry about twice the second moment of all stars; decide whether remnants stay in it, with the black-hole mass function checked | none | done, PR #2164 |
| GEN.542 | scripts/star_type_by_radius.py fails with a TypeError because the plan arguments now default the disk and arm options to None (bug) | none | done, PR #4309 |
| GEN.543 | Define "usable" for a P-type window so the design table's pair-node row (34% falling to 18%) can be reproduced | none | done, PR #2238 |
| GEN.544 | Revisit the disc's vertical profile: the sech-squared core has half the midplane density of the exponential fits its scale heights come from | none | done, PR #2323 |
| GEN.545 | The mean system mass (0.4 solar, recalled within a factor of 1.5) sets every galaxy's system count once GEN.404 lands: derive it from the generator's own mass distribution | none | done, PR #2434 |
| GEN.546 | Draw a disk star's age from the star formation history inside its window, not uniformly | none | done, PR #2268 |
| GEN.547 | The local B-star share is 0.12 percent against 0.04 observed: find what keeps it three times too high (strict expected failure) | none | done, PR #2268 |
| GEN.548 | Build-time downloads for the isochrone grids, the white-dwarf layers and the neutron-star equation of state, each with a checksum and a labelled stand-in until it lands | none | done, PR #2236 |
| GEN.549 | L dwarfs hotter than 1,500 K are 2 to 4 times too rare against the 20 pc census: use a cooling track that reaches them (strict expected failure) | none | done, PR #2247 |
| GEN.550 | Umbrella (placeholder): one recursive generation process for anything that can have something else orbiting it, from a planet with moons up to a group of star-like objects that each have a system | none | open |
| GEN.551 | Approve the unified orbited-body design note and its eight cards (D1 to D8) before any build of the chain | none | done, PR #2171 |
| GEN.552 | Host registry: the five fixed `system_model` values become pointers to one gated orbiter table per host class, with no change in behaviour | none | done, PR #2214 |
| GEN.553 | One star-like body draw: origin, birth, state and [Fe/H], with the initial-final relation, the remnant mass, radius and age rules in one place | none | open |
| GEN.554 | Remnant counts follow from the mass function and the star formation history, with a kick step, and Boss's research densities become test targets | none | open |
| GEN.555 | The pair builder for any masses: planet pairs, brown dwarf pairs and a planet with its moons share one stability test | none | open |
| GEN.556 | `populate(host)`: planets, belts, comets and moons come from the class tables by one recursive routine, and the system pass is absorbed | none | open |
| GEN.557 | Wanderers: a rogue planet becomes a planet with no parent, drawn by the one body chain with an origin table, a notional formation host, the population's age and motion | none | open |
| GEN.558 | Moons of wanderers: ejection survivors and isolated-disk moons are table rows, and the fixed 20% boolean is gone | none | open |
| GEN.559 | Unbound comets: a size law and a cut, composition from the comet table, derived activity and velocity | none | open |
| GEN.560 | Research list for the chain: ejection survival of moons, circumplanetary disks of isolated planetary-mass objects, failed-supernova and direct-collapse shares, the size law of interstellar objects | none | done, PR #2924 |
| GEN.561 | The gear and equipment tier follows the words: the kit named on a body is derived from the habitability sentence, not drawn on its own | none | done, PR #2240 |
| GEN.562 | One function for the white-dwarf initial-final relation and for the neutron star and black hole mass draw, with no change in what any seed generates | none | done, PR #2220 |
| GEN.563 | The black hole density follows the mass function and the verified failed-supernova share: one constant changes, the neutron star density stays | none | done, PR #2256 |
| GEN.564 | A refill job for a galaxy under an older generator epoch: clear the fill-level rows of every sector and rebuild them without a rescatter | none | folded into GEN.764 (Core 1, 7:32 PM EDT 2026-10-24) |
| GEN.565 | The mass object core and its three family classes: StellarObject, PlanetaryObject and SmallBody, with state and parent as data, replace the nine kind classes | none | open |
| GEN.566 | Research: the size-frequency law of small bodies, for belt members drawn on demand and the comet size cut, beside the interstellar size law of GEN.560 | none | done, PR #2800 |
| GEN.567 | Design note: the metallicity-dependent black hole density: the form of f(Z), the birth [Fe/H] of a remnant, and what the scatter reads | none | done, PR #2190 |
| GEN.568 | One registry of every unit the project uses | none | done, PR #2476 |
| GEN.569 | The Quantity value class: the number as given, every other unit answered on first read and kept, declared relations, checks | none | done, PR #2490 |
| GEN.570 | Field declaration: every unit-valued field of the mass object classes declares its kind and native unit, and the one field table serves the serializer | none | open |
| GEN.571 | Replace the ad-hoc unit conversions with the registry, by file group, bit-identical on the golden set | none | open |
| GEN.572 | The two rounded unit factors become exact: 63241.077 AU per light-year and 3.26156 light-years per parsec, as one epoch-changing change | none | done, PR #2243 |
| GEN.573 | The position class reads its units from the registry: length unit by key, no bare 1000 or copy of c, mass and mu a declared relation | none | done, PR #2490 |
| GEN.574 | Describe and API: qty() builds a Quantity and the {kind, si} shape becomes {kind, value, unit} | none | done, PR #2806 |
| GEN.575 | Sweep the mass object classes: the 35 unit-suffixed attributes and every quantity held under two units get one declared native unit | none | open |
| GEN.576 | system(host): the one function that groups bodies by their parent link and returns the barycentre sums, the central object's offset and the valid-system check | none | open |
| GEN.577 | Stored planet and moon states become the state of the barycentre of the body and its own members | none | open |
| GEN.578 | The point class gets the local frame and the barycentre system frame, and the system anchor is no longer the nearest star | none | open |
| GEN.579 | The hosting gate: can_host(body) in three parts and the gated table of the satellite share by size class, as labelled stand-in rows | none | done, PR #2800 |
| GEN.580 | Research: the binary and satellite fraction of small bodies by size, and the survival of sub-satellites of moons, to replace the stand-in cut of the hosting gate | none | done, PR #2800 |
| GEN.581 | primary_first(members): before a system is validated the heaviest branch becomes the primary and the others are recalculated, one function for the tree, the fill, moons and free pairs | none | open |
| GEN.582 | Origin and present-state split: fields marked origin or state in the Field table, the draw record held until validation finishes, and the zone and class re-run from the seed | none | open |
| GEN.583 | reparent: the state-preserving change of a body's parent, shared by re-primary at an event, an ejection, a capture and a merger | none | open |
| GEN.584 | The core's kinematic record: velocity, spin axis and rate, mass and epoch in the point, and the derived columns dropped | none | open |
| GEN.585 | Dimensions and unit systems in the registry: every unit carries its exponent vector, the derived kinds are listed, and a constant's number in any coherent system follows from its dimension | none | done, PR #2542 |
| GEN.586 | physics/relations.py: the named laws declared once as data, with a fast float kernel in SI and a Quantity wrapper generated from the same declaration, and a dimension test in the mathcheck gate | none | done, PR #2542 |
| GEN.587 | Replace the hand-written relation sites with the relations module, by file group, bit-identical on the golden set | none | open |
| GEN.588 | Optional: operator arithmetic on Quantity for the edges (mass divided by volume gives a density), unscheduled | none | open |
| GEN.589 | Research: the central molecular zone's upper cluster mass and its cloud-shock disruption, so the young clusters of the core match the few that are observed | none | done, research report (gen589-cmz-clusters-2026-10-31.md) |
| GEN.590 | Wire the downloaded isochrones, white-dwarf sequences and neutron-star tables into star composition (after Boss answers the Require / Stand-in card) | none | done, PR #2282 |
| GEN.591 | An elliptical comet's aphelion can exceed its star's Hill sphere in about 10% of runs: the comet bound test fails on main (bug) | none | done, PR #2306 |
| GEN.592 | A capped dense-core cell generates 0 systems when the star scatter was not drawn, while the resolution says 4 to 30 (observation, not traced) | none | done, PR #2327 |
| GEN.593 | Pulsar planets: the disc range stops short of PSR B1257+12's planets, and two of the three formation paths (a captured giant, an ablated companion) leave a pulsar with no planets | none | open |
| GEN.594 | Each AD dataset (AD1 to AD8) ships its best possible fallback, used only when a run is told explicitly to disregard missing files | none | open |
| GEN.595 | Replace the Solar System comet population shape (GEN.520 D6) with a derived distribution | none | open |
| GEN.596 | A companion below 0.01 solar masses makes multistar.draw_member_tree raise 'no cooling track' (bug) | none | done, PR #2428 |
| GEN.597 | A rocky planet of a wide binary gets flux_variation above 1 (the insolation test passes only on seeds 0 to 59) (bug) | none | done, PR #2399 |
| GEN.598 | CLUSTER_MEAN_SYSTEM_MASS_SOL (0.4, recalled) sets the system count of globular clusters and satellites: derive it as GEN.545 derives the galaxy's | none | done, PR #2478 |
| GEN.599 | A run with missing downloaded data files refuses unless it is told explicitly to disregard them, and then uses the best possible fallback and says so in its header | none | done, PR #2560 |
| GEN.600 | The updater and the downloader print the URL of every file they load | none | done, PR #2464 |
| GEN.601 | Shape dials: the lever table, one registry row per slider (name, flag, label, unit, minimum, maximum, step, Auto default, the law it moves, range source) and the parity test | none | done, PR #2395 |
| GEN.602 | Shape dials: the command-line flags built from the lever table, and the old absolute shape options retired | none | done, PR #2403 |
| GEN.603 | Shape dials: the Scale dial: mass override, compactness in sigma, density, and the any-two-of-three pin logic with the third shown | none | done, PR #2403 |
| GEN.604 | Shape dials: the Structure dial: flatness (scale length over scale height), bulge share and bulge size (scale radius over L) | none | done, PR #2403 |
| GEN.605 | Shape dials: the Edge dial: where the galaxy ends, with the resulting radius and the expected sector and system counts shown before a plan | none | done, PR #2403 |
| GEN.606 | Shape dials: the Halo dial: dark halo mass over stellar mass and the stellar halo share, with r200 shown and a size-halo warning | none | done, PR #2403 |
| GEN.607 | Shape dials: the Census dial: star richness (systems per solar mass, equivalently the mean system mass), advanced and hidden by default | none | done, PR #2443 |
| GEN.608 | Shape dials: the Shape selector, added to the table and the Generate form unfilled: key `spiral`, flag `--shape`, one option for now | none | done, PR #2395 |
| GEN.609 | Shape dials: check the range source of every slider (size-law mass span, flatness, bulge share, bulge size, edge, halo ratios) and record each as verified, derived, bracket or placeholder | none | done, PR #2403 |
| GEN.610 | The `--nuclear-amplitude` help text defines 1 as 'the observed' nuclear cluster and disc (2.5e7 and 1.5e9 solar masses), which are Milky Way values: replace them with a draw from the galaxy's own laws | none | done, PR #2403 |
| GEN.611 | Shape dials: the pattern dial: arm contrast (crest minus trough over crest plus trough, default 0.4), with arm count and pitch angle, replacing the two arm densities | none | done, PR #2403 |
| GEN.612 | Extend the halo, stellar-to-halo and size laws down to dwarf and ultra-faint masses (below 1e8 solar masses), so the shape dials' minimum mass is reachable | none | open |
| GEN.613 | Investigate the limits of the density algorithm: for which mass, size, flatness and bulge share the galaxy density model stays valid, and what it does outside them | none | open |
| GEN.614 | The nuclear disc is not softened to the grid: switch it off when its length (0.057 times the bulge radius) is under 2 pc, or soften it like the nuclear cluster | none | open |
| GEN.615 | Galaxies below the low-mass cut are one disc plus a bulge: no thick disc, bar, stellar halo or nuclear structures | none | open |
| GEN.616 | Average the galaxy density over each planned sector (or a 3 by 3 by 3 sample) instead of taking its centre point, for galaxies under 28 pc scale length | none | dropped (Boss, 2026-10-21 8:08 PM EDT, replaced by GEN.619) |
| GEN.617 | Show the density ratio D, the share of stars in planned sectors and the outline size under the dials and in the plan output | none | done, PR #2469 |
| GEN.618 | The stellar-to-halo inversion gives a 6e13 solar-mass halo for 2e11 solar masses of stars, the shape dials' mass ceiling: check the Moster inversion at the top end | none | open |
| GEN.619 | A variable sector edge for small galaxies: 4, 2 or 1 pc chosen by the thin disc scale length, with the one-star outline rule as a density threshold | none | open |
| GEN.620 | The region removal command: its own admin command that deletes the halo or the satellite galaxies from a made galaxy and makes their cells invalid addresses (region argument halo or satellites, needs the switch already off) | none | open |
| GEN.621 | The Scale dial's Auto default draws a galaxy about half the Milky Way's size, with an epoch bump | none | done, PR #2418 |
| GEN.622 | pulsar_radial_factor divides by zero beyond about 580 kpc galactic radius (bug) | none | done, PR #2739 |
| GEN.623 | Data tables carry no licence column or licence check; the licence information lives in DATA-LICENSE.md | none | done, PR #2596 |
| GEN.624 | Satellite galaxies scale with the host galaxy's mass: a small galaxy no longer pays for 11 satellites (decided by Boss; needs a generator epoch) (bug) | none | done, PR #2638 |
| GEN.625 | Tune the Auto galaxy defaults so the plan coverage warning never fires and the average default galaxy is at most half the Milky Way's size (decided by Boss; needs a generator epoch) (bug) | none | done, PR #2627 |
| GEN.626 | The plan's crest line and planned-share figure are wrong: print them from the built density, not a table and a fixed system mass (bug) | none | done, PR #2598 |
| GEN.627 | The outline takes satellites and clusters in by tidal reach, so the candidate sector count is 16 times too high and random starts land in empty cells (bug) | none | done, PR #2638 |
| GEN.628 | The nuclear disc is a fixed 61 times the nuclear cluster, 21 percent of a small galaxy's stars against about 3 percent for the Milky Way: draw it from a law (bug, needs a generator epoch) | none | done, PR #2620 |
| GEN.629 | The comet volatile ratio S=(CO+CO2)/H2O scatters by only 0.1 dex, too narrow against single comets: 67P's 7.8 percent is 3.6 sigma low, so the draw seldom gives its CO and CO2 together | none | done, PR #2744 |
| GEN.630 | GEN.316 part (a): the N-star StarSystem host and the stars.role 'member' migration, still building only the pairs, with readers and loaders updated | none | done, PR #2689 |
| GEN.631 | GEN.316 part (b), narrowed: stages D and G draw their members from draw_member_tree, trees of 3 or more stars build as member stars, and the old companion draw and the two-star probability table are deleted | none | done, PR #2915 |
| GEN.632 | Two code paths still hold their own vapour-pressure constants: rogue_surface.nitrogen_vapor_pressure_pa and the water curve in hydrosphere.py should read the one saturation function | none | done, PR #2938 |
| GEN.633 | The backfill after a run puts its dimmer bright stars only in sectors not yet generated, so a sector holds different systems depending on the order the sectors were generated in (bug) | none | done, PR #3043 |
| GEN.634 | The plan outline takes a whole shell in by its highest band, so the candidate sector count is 3.5 billion where about 1.0 billion expect a system (interval outline) | none | open |
| GEN.635 | The cluster system has no radial profile yet, and the dwarf end of the globular cluster count uses no occupation draw (the Carlsten fits are not built in GEN.403) | none | open |
| GEN.636 | StarSystem._orbit_ceiling_au no longer equals validation.orbit_ceiling_au after the planet survival limit, so test_system_methods_still_delegate fails on main (bug) | none | done, PR #2653 |
| GEN.637 | Re-run the GEN.625 1,000-seed all-Auto sweep on main after GEN.624 and GEN.403 (#2638) | none | done, sweep report gen637-all-auto-sweep-2026-10-28.md, no code |
| GEN.638 | Planets and belts can sit past the survival ceiling because the age is raised after they are placed, so test_wide_binary_system_generates_and_stays_within_bounds[G2V] fails in about 0.4 percent of draws (bug) | none | done, PR #2695 |
| GEN.639 | The planet-mass habitable zone limits and the pre-main-sequence runaway history left open by GEN.146 (design note 4.1) | none | done, PR #3020 |
| GEN.640 | GEN.316 part (c): delete BinaryStarProxy, WideBinaryPair and the writers of the binary_* columns, once circumbinary and S-type planets come from the tree | none | open |
| GEN.641 | Research: the interior of a giant planet, ice giant and brown dwarf: the layers (molecular, liquid and metallic hydrogen and helium, helium rain, a dilute core), whether there is a rocky or icy core, and what the core is made of | none | done, PR #2800 |
| GEN.642 | Giants, ice giants and brown dwarfs get a full interior composition: the molecular envelope, the liquid and metallic hydrogen and helium layers, the dilute zone, and a rocky, iron or icy core with its own composition, drawn and described | none | done, PR #2912 |
| GEN.643 | Habitable means PHI-6 overall Blue or Green: one rule, one function, read by the sector counts, the markers, the lists, the search and the descriptions | none | done, PR #2886 |
| GEN.644 | Proposal: tie the PHI-6 into the planet classification (class letters, a habitability tier, or a suffix) and the rule for where a class and a PHI-6 colour must agree | none | open |
| GEN.645 | Design: condition-based planet and moon classes: the rule table, the inputs, the legacy letters and what happens to GEN.28, GEN.90 and GEN.91 | none | done, PR #3080 |
| GEN.646 | Class designation leaves generation: planets and moons are drawn from their zone, mass, formation distance and composition recipes with no class letter | none | done, PR #3729 |
| GEN.647 | The class is assigned after the body is complete: classify(body) from the stored atmosphere, temperature, terrestrial composition and gas composition, for planets and moons | none | open |
| GEN.648 | The generation request for a habitable world (+habitable_world) steers by conditions, not by a class list: redraw water, life, temperature and atmosphere inside their ranges until the PHI-6 overall is Blue or Green | none | open |
| GEN.649 | The four life scores: Microbial, Flora, Fauna and Sentient, each 0 (none) to 9 (abundant), from the evolutionary stage and the surface conditions; human operability is removed | none | done, PR #2930 |
| GEN.650 | Design: the planet mass law once the class no longer supplies the mass: the isolation mass at the formation distance, the snow line and the host, replacing the class mass ranges of decision 28 (design note, prerequisite of GEN.646) | none | done, PR #2947 |
| GEN.651 | Wire the star's [Fe/H] into the star evolution: birth [Fe/H] is drawn before the state and passed to the lifetime and the giant branch, and the composition function is read by more than the white-dwarf flag | none | done, PR #2952 |
| GEN.652 | The two main-sequence tables disagree: evolve_star reads mass-only piecewise laws and star_makeup reads MIST; one becomes the source of the other (bug) | none | done, PR #3304 |
| GEN.653 | Suitability's food rung reads the Flora and Fauna scores: Flora is always good food, both are better, Fauna includes all animal life | none | done, PR #2930 |
| GEN.654 | Aquifers: a world with a dry, frozen or too hot surface can hold liquid water underground, within reach of the 273 K level; the depth feeds the Water and Suitability rungs and the Microbial score only | none | done, PR #3221 |
| GEN.655 | Survival equipment: one physics/equipment.py returns the list of needs (required or recommended) for the four slots, and replaces habitability.equipment and protection.requirements | none | done, PR #2766 |
| GEN.656 | Regional temperatures (Polar, Northern, Equatorial, Southern; day side, terminator and night side for locked worlds) and the day and night swing, with the regional sentences of the detail text | none | done, PR #2844 |
| GEN.657 | Toxic-species rows in the equipment limits table, wired to the stored gases: required gas filter grade by IDLH, with every recalled limit checked against NIOSH and ACGIH | none | done, PR #3702 |
| GEN.658 | Dust loading and dust toxicity feed the dust filter and seals level of the equipment needs | none | done, PR #2970 |
| GEN.659 | The interstellar debris density stand-in of 1e12 per pc^3 is 700 to 1,800 times too low: the measured value is 7e14 to 1.8e15 per pc^3 (bug) | none | done, PR #2801 |
| GEN.660 | The complex_organics row of the materials table has a density of 1.5 g/cm3 while Choukroun 2020, used by the comet draw, has about 2 g/cm3 (bug) | none | done, PR #4310 |
| GEN.661 | Bright-star counts still use the reference-composition fraction (stand-in, bug) | none | done, PR #3920 |
| GEN.662 | Cluster and satellite stars draw a halo [Fe/H] for their lifetimes (stand-in, bug) | none | done, PR #3920 |
| GEN.663 | A chemical-evolution birth-epoch model for thorium and uranium: row o4 draws the formation time from 3 to 13 Gyr | none | done, PR #2822 |
| GEN.664 | The thin-disc offset rows disagree by 0.12 dex: Botelho solar twins (+0.08 dex) against the Hypatia thin-disc statistic (-0.044) | none | done, PR #2815 |
| GEN.665 | Two tuning targets are weaker than they look: the white-dwarf redox sample is 7 stars and the mass-radius cut barely filters | none | done, PR #2815 |
| GEN.666 | The life text on a planet page contradicts itself: one paragraph says life would evolve at a normal pace, the next says hyper-accelerated (bug) | none | done, PR #2799 |
| GEN.667 | Ozone reads the oxygen partial pressure through the mixture accessor once the child table lands | none | done, PR #3373 |
| GEN.668 | Atmosphere loss events (GEN.189) do not yet reduce the stored ozone column | none | done, PR #3307 |
| GEN.669 | Residuals of the generator laws after GEN.450 and GEN.270: the star-formation profile shares and the Algol rate against their observed windows | none | open |
| GEN.670 | Research for body rendering: band, storm and Rossby scale laws, vegetation tint by host spectrum, ring colour by composition | none | open |
| GEN.671 | Where a hardened suit works against the radiation, the needs say "required hardened suit" and "recommended hardened habitat" | none | done, PR #2883 |
| GEN.672 | Two gas-equilibrium modules with two copies of the NASA thermo data are on main: merge them into one table and one solver | none | done, PR #3061 |
| GEN.673 | Every giant stores its photosphere pressure (about 0.1 bar on the day side) before draw_mixture, which the hot-Jupiter rows gate on | none | done, PR #3373 |
| GEN.674 | Atmosphere mixture recipe gaps: no perihelion gate on rows k, a near-zero row normalised to 1, rows j3 and j4 never drawn, no carbonate-silicate CO2 sink, no Pluto-like family, missing lava-vapour gases (bug) | none | done, PR #3072 |
| GEN.675 | Register the atmosphere mixture stand-ins in the gap plan registry: kind, evidence label, knobs, validity range and falsifier test for each | none | done, PR #3037 |
| GEN.676 | Check the NH3 photolysis activity factor (1.078) against the published Huebner photodissociation table and replace the Boss-supplied value | none | dropped (Boss, 11:53 PM EDT 2026-10-23: the gap cannot be filled, no source is reachable; the factor 1.078 stays as a labelled stand-in) |
| GEN.677 | The phenomena scatter places about 86,000 stellar black holes per shell pair and no ordinary neutron stars (only magnetars) (bug) | none | done, PR #2968 |
| GEN.678 | The bright-star units expected about 23,000 stars but placed 1,000 to 3,000, so the run placed 229,596 stars against 365,093 expected and took 20 minutes (bug) | none | done, PR #2997 |
| GEN.679 | classify(body) and planet_classes.json: the rule table as data and one pure function, with a total-partition test | none | done, PR #2781 |
| GEN.680 | Class shadow sweep: compare classify with the drawn class and fit the thresholds, and report the shares | none | done, PR #2781 |
| GEN.681 | The kind of a planet or moon comes from formation (mass against the snow line and the host) and replaces the class draw | none | open |
| GEN.682 | Radius, density and gravity come from the bulk composition; delete radius_range, size_mode and density_range | none | open |
| GEN.683 | Atmosphere recipe gates move from class to kind and regime; delete CLASS_MIXES, atm_density_range and atm_molar_density_range | none | done, PR #3311 |
| GEN.684 | Albedo and greenhouse come from the composition; delete the two class ranges | none | open |
| GEN.685 | Water comes from the composition and the origin; delete WATER_MASS_FRACTION[class] | none | open |
| GEN.686 | Age needs come from the processes; delete age_ranges | none | done, PR #3080 |
| GEN.687 | Rogue planets get their class from classify; delete choose_rogue_planet_class | none | open |
| GEN.688 | Generation requests name conditions, not class letters: slots and the admin class change take conditions, +habitable_world steers on the PHI-6, the directive counts is_habitable | none | open |
| GEN.689 | The class cut-over: classify at the end of generation, validation becomes total and derived, one generator epoch, a migration backfills the column | none | open |
| GEN.690 | Delete PLANET_CLASSES generation fields, PLANET_CLASS_PROBABILITIES, LIFE_PLANET_CLASSES, PLANET_CLASS_MAX_LIFE_STAGE and MOON_BLACKLIST; update the golden outputs | none | open |
| GEN.691 | The Suitability sentence says "no food readily available, but there is water and shelter available" while the Water row is Red and "undrinkable without heavy treatment" (bug) | none | done, PR #3162 |
| GEN.692 | The planet page Atmosphere table shows "Oxygen 100%" for a class B planet whose text says helium, sodium and oxygen, because it divides only the stored gases (bug) | none | done, PR #2799 |
| GEN.693 | Light-gas escape filter after outgassing: the mixture draw gives a small warm body hydrogen as its main gas | none | done, PR #2862 |
| GEN.694 | Oxygen caps in the atmosphere draw: biotic oxygen at most 35 percent, abiotic at most 95 percent with a nitrogen and argon floor, and no mixture shown as 100 percent | none | done, PR #2862 |
| GEN.695 | Class B keeps some hydrogen and helium in its atmosphere mix next to the sodium and oxygen of a molten world | none | dropped (Boss, 1:22 AM EDT 2026-10-24, card option 3: class B never has hydrogen or helium; a stripped core is at least 5 Earth masses) |
| GEN.696 | Research source: the leftover envelope mass fraction of a stripped giant core | none | dropped (Boss, 1:22 AM EDT 2026-10-24: class B never has hydrogen or helium) |
| GEN.697 | The generator records an origin quantity stripped_core for a small hot body that began as a giant | none | dropped (Boss, 1:22 AM EDT 2026-10-24: class B never has hydrogen or helium) |
| GEN.698 | Recipe rows for the stripped core: a hydrogen and helium envelope share for a small hot body, replacing the gating of rows e1 to e6 that stops at 1 Earth mass | none | dropped (Boss, 1:22 AM EDT 2026-10-24: class B never has hydrogen or helium; lane L's research note stays as reports/stripped-core-envelope-2026-10-24.md) |
| GEN.699 | Model formation (primordial) heat separately instead of dividing the radiogenic heat by a drawn Urey ratio of 0.1 to 0.9 | none | done, PR #3046 |
| GEN.700 | Tidal heating of moons counts in the energy term of the habitability score: keep the tidal flux on the moon and add it to the heat flow | none | done, PR #2989 |
| GEN.701 | Rogue planets get habitability scores: score columns on rogue_planets, a World built from the rogue surface regime, and the API fields | none | done, PR #3265 |
| GEN.702 | Some equipment items sit in the wrong habitability row: self-contained air under a Blue Air row on cold worlds, and a hardened pressure suit under Radiation while the Pressure row names no gear (bug) | none | done, PR #2920 |
| GEN.703 | The greenhouse optical depth (tau0) comes from the stored mixture and the band opacities, replacing the per-gas strengths | none | done, PR #2992 |
| GEN.704 | Moons of a moved planet keep the dose, solar particle exposure, UV and habitability of their old distance from the star (bug) | none | done, PR #2989 |
| GEN.705 | The metallicity ([Fe/H]) factor on close companions no longer reaches the fill since the tree draw replaced the old companion draw (bug) | none | done, PR #4385 |
| GEN.706 | Low-mass M dwarfs get fewer companions than the star-count table, forced pairs around M8 primaries come back unmet, and compact companions of living stars are cut from fill trees (bug) | none | done, PR #4385 |
| GEN.707 | Register the white dwarf binary stand-ins in the gap plan registry: kind, evidence label, knobs, validity range and falsifier test for each | none | open |
| GEN.708 | White dwarf merger products draw their 0.8 to 1.1 solar mass bump: the star model gives the white dwarf its mass before the member tree, so the tree only flags 10 percent as merger products | none | open |
| GEN.709 | Stored close-pair star positions disagree with the orbit_nodes tree after the multi-star work: a circular pair keeps a random arg_periapsis, and deeper trees are far off (bug) | none | done, PR #3017 |
| GEN.710 | A pre-placed row with a stored tree can have a member heavier than the row star at stars[0], so the nebula it lights and the star page disagree with the built sector (bug) | none | done, PR #3033 |
| GEN.711 | Register the planet mass law stand-ins in the gap plan registry: kind, evidence label, knobs, validity range and falsifier test for each | none | done, PR #3596 |
| GEN.712 | The materials row complex_organics has a density of 1.5 g/cm3 but Choukroun 2020 gives about 2 g/cm3 (bug) | none | done, PR #2995 |
| GEN.713 | Research the size law and the cut for long-period and interstellar comets (blocks GEN.559) | none | done, report reports/gen713-comet-size-law-2026-10-24.md (Lane K, about 2:30 PM EDT 2026-10-24 as the lane reported) |
| GEN.714 | Moon inclination is drawn from 0 to 25 degrees, so captured moons (retrograde, high inclination) are unreachable (bug) | none | done, PR #3084 |
| GEN.715 | Real infrared absorption tables (k-tables or line lists) for the greenhouse opacity, replacing the labelled band-strength stand-ins of GEN.703 | none | done, PR #3507 |
| GEN.716 | add_star_hosted_nebulae lights a nebula from stars[0], which for a pre-placed row with a stored tree is the structural primary and not the row star, so the built nebula can differ from the predicted one (bug) | none | done, PR #3033 |
| GEN.717 | The greenhouse opacity has no water-vapour vertical profile, so the 15 micron CO2 band is not masked by water and a CO2 doubling warms 0.3 to 0.6 K instead of about 1.1 K | none | done, PR #3507 |
| GEN.718 | Cold haze worlds like Titan come out about twice too warm in the spectral greenhouse opacity | none | done, PR #3507 |
| GEN.719 | A white dwarf primary makes draw_member_tree ignore n_stars=1, so a binary_system prevalence of -100 still yields a wide double white dwarf (bug) | none | done, PR #3162 |
| GEN.720 | star.age is adjusted after the planets are generated (adjust_age_for_planets), so any age-dependent planet value goes stale (bug) | none | done, PR #4313 |
| GEN.721 | A tectonic-mode switch (mobile lid or stagnant lid) for the planetary heat model | none | done, PR #3510 |
| GEN.722 | Radius valley follow-up to GEN.509: the mechanism choice, the stellar-mass slope (e-link) and the width stand-ins | none | done, built with GEN.822 (radius_valley.py, rows RV1 to RV4) |
| GEN.723 | Refresh the hot-Jupiter recipe rows h1 to h3: their edge points were computed with the deleted gas-equilibrium module and sit up to 16 K (0.6 percent) above the surviving solver | none | done, PR #3493 |
| GEN.724 | Register the GEN.272 stand-ins in the gap plan registry and correct the registry row millisecond-pulsar-companion (5 percent) to the 0.2 percent the code keeps | none | done, PR #3808 |
| GEN.725 | Stand-in values carry a source tier and a tweak level: two new columns in the gap plan registry and one document in docs listing every stand-in, where it came from, its tier and how far it was tuned for outcome | none | done, PR #3163 |
| GEN.726 | Surface pressure supply law: atmosphere mass from delivered volatiles, degassing, escape and sinks, then pressure from mass and gravity | none | done, PR #3111 |
| GEN.727 | Eccentricity and inclination for planets and moons, with the seasons, eccentric insolation and tidal heating that follow | none | done, PR #3380 |
| GEN.728 | Planet count from the observed multiplicity and occurrence distributions by star mass and metallicity | none | done, PR #4105 |
| GEN.729 | Planet spacing from the observed period-ratio and size-uniformity statistics, with stability checked | none | done, PR #4105 |
| GEN.730 | Formation-era dynamics: giant impacts, scattering and instability instead of one survival constant | none | done, PR #4105 |
| GEN.731 | Planet spin: day length from accretion angular momentum, tidal despinning with per-body Q and k2, and spin-orbit resonances | none | done, PR #3380 |
| GEN.732 | Planetary rings: existence, mass, extent and composition from the Roche limit and a disrupted moon or captured body | none | done, PR #3348 |
| GEN.733 | Hydrogen-helium equation-of-state tables for giant interiors: download job and use in the interior model | none | done, PR #4436 |
| GEN.734 | Open weight tables of the atmosphere and bulk recipes: planet-minus-star Fe/Mg scatter, core light elements, rogue-planet stripping, moon composition against distance | none | done, PR #3208 |
| GEN.735 | Interior heat beyond the stagnant lid: mobile-lid cooling, a separate core, a pressure-dependent liquidus and coupling to climate | none | done, PR #3510 |
| GEN.736 | Magnetic field from a dynamo: the core's heat flow and rotation (buoyancy-flux scaling), reading the interior heat model | none | done, PR #4436 |
| GEN.737 | Tidal heating from the orbit's real eccentricity for rocky bodies and planets, with the random 1 to 50 boost removed | none | done, PR #3380 |
| GEN.738 | Atmospheric escape: a thermosphere model, hydrodynamic escape and time-integrated XUV loss | none | done, PR #3260 |
| GEN.739 | Impact erosion and late delivery of atmospheres, with history | none | done, PR #3246 |
| GEN.740 | Life and the air: a coupling between the life stage and the gases (oxygen from a biosphere) | none | done, PR #3862 |
| GEN.741 | Post-main-sequence stellar tracks: subgiant, red giant, horizontal branch, asymptotic giant and supergiant from mass and age instead of random draws | none | dropped (revoked by the full-build sweep 8:16 PM EDT 10-29: superseded by GEN.935, which reads the MIST post-main-sequence rows; its synthetic-cluster test moved to GEN.935) |
| GEN.742 | White dwarf cooling sequences by mass and atmosphere, including crystallisation and spectral type against temperature | none | done, PR #3493 |
| GEN.743 | Neutron star equation of state: integrate the tables for radius and structure, and add helium-core and oxygen-neon white dwarfs | none | open |
| GEN.744 | Alpha-enhanced stellar tracks and an age-metallicity relation from the galaxy's own history | none | open |
| GEN.745 | Stellar spectra by temperature, gravity and activity (model atmospheres) instead of a blackbody with a seeded UV share | none | open |
| GEN.746 | Cosmic-ray dose below 100 g/cm2 of air from shower simulations instead of a straight stand-in between the Moon and Mars | none | done, PR #3920 |
| GEN.747 | Solar energetic particle events: the full event set and the Weibull parameters | none | done, PR #3920 |
| GEN.748 | Ozone: temperature-dependent Chappuis cross sections and a spread no wider than the paper's | none | done, PR #3920 |
| GEN.749 | Trapped radiation belts and the magnetosphere's shielding beyond a dipole size | none | done, PR #3410 |
| GEN.750 | Evolutionary timeline: a model of the time to each life stage from conditions (Boss ruled: model it from physics) | none | done, PR #3663 |
| GEN.751 | Add the lava-vapour gases Al, K2, KO, K2O and Na2O to the gas table, from the research parts, so the lava-world air has its full species list | none | done, PR #3222 |
| GEN.752 | Climate feedback on the surface temperature fixed point: ice-albedo, water vapour and cloud albedo | none | open |
| GEN.753 | Validation sweep for each physics area: compare generated populations with published statistics | none | open |
| GEN.754 | Wiring guard: the stand-in registry says what the generator reads, and a test fails for a table or module nothing reads | none | done, PR #3808 |
| GEN.755 | Check whether the downloaded MIST grid holds post-main-sequence phases; if not, add a setup fetch under GEN.741 | none | done, absorbed by GEN.935 (the MIST grid holds the post-main-sequence rows; datasets/isochrones.py keeps only the main sequence) |
| GEN.756 | Piece of GEN.411: the body-to-mixture-subject builder (physics/mixture_subject.py) | none | done, PR #3129 |
| GEN.757 | Piece of GEN.411: the ocean water-vapour cap and the sulfur uptake functions (physics/mixture_sinks.py) | none | done, PR #3129 |
| GEN.758 | Piece of GEN.411: mechanical test fixes for the deleted class ranges | none | done, PR #3311 |
| GEN.759 | Wire the r-process Th/U enrichment (rprocess_enrichment, GEN.663) into the radiogenic heat module | none | done, PR #3510 |
| GEN.760 | Wire the brown dwarf spectral-type map (brown_dwarf_types, GEN.517) into the brown dwarf draw | none | done, PR #3317 |
| GEN.761 | Wire the satellite hosting gate (hosting_gate, GEN.579) into the moon draw | none | done, PR #3258 |
| GEN.762 | Wire the clathrate and eutectic brine limits (clathrate_brine, GEN.471) into the icy body atmosphere and ocean draw | none | done, PR #3128 |
| GEN.763 | The Urey ratio can exceed 1 for small cold-start bodies (a 0.0004 Earth-mass moon gives 3.02) and the docs table covers only Earth-size bodies (bug) | none | done, PR #3510 |
| GEN.764 | Bring a galaxy to the current generator epoch without a full regeneration: record the epoch each stored result was made under and redo only what an epoch bump changed | none | done, PR #3439 |
| GEN.765 | A world with Flora 5, Fauna 5 and Sentient 5 scores Suitability Red with "no food readily available": the food rung ignores the life the scores say is there (bug) | none | done, PR #3203 |
| GEN.766 | Generated life text states what is detected, in the present tense: no "speculative", "could", "would" or "boasts" (bug) | none | done, PR #3203 |
| GEN.767 | Day and night time shares for every object that has a day or night, with the extra night when the parent planet blocks the sun | none | done, PR #3510 |
| GEN.768 | Hours a person can be outside in a day: daylight hours limited by heat, cold, pressure, air supply and radiation dose, computed and stored per planet and moon | none | done, PR #4321 |
| GEN.769 | Small lengths read in a sensible unit: an ice shell of 0.00109 km shows as 1.09 m, not in km (bug) | none | done, PR #3478 |
| GEN.770 | Water rows state their unit: share of the mass, surface under ocean, surface above sea level, water activity, ocean pH (bug) | none | done, PR #3478 |
| GEN.771 | Property sheet: the 40 ratio and 7 index rows state a unit or a scale (shares in percent or ppm, scores as "n of 9", ratios by what they divide) (bug) | none | done, PR #3478 |
| GEN.772 | Number-display guard: every displayed quantity goes through the unit registry formatter, and a test fails on a bare number in a page, an API document or a command output (bug) | none | done, PR #3478 |
| GEN.773 | A temperature difference is never converted as a temperature: the day and night swing of 9.07 K shows as 9.07 K (or 9.07 C degrees), not -264 C (bug) | none | done, PR #3478 |
| GEN.774 | Settlement draw streams: the named draw.py streams for nation count, dominance, city sizes, city names, city positions and the city to nation assignment, and their stand-in registry rows | none | open |
| GEN.775 | Replace the planet-minus-star Fe/Mg scatter and super-Mercury tail (registry row BC8) with a larger sample of rocky planets with precise masses and stellar Fe/Mg | none | done, PR #3846 |
| GEN.776 | Replace the rogue-planet volatile loss at ejection (registry row BC11) with a measured atmosphere or free-floating-planet volatile survey | none | closed, no measurement of a rocky free-floating planet's volatiles exists (Sources 1 searched) |
| GEN.777 | Replace the core light-element split (registry row BC18) with a second body that has seismic or core-density data | none | done, Mars check made; the wider core light-element range is filed as a wiring item |
| GEN.778 | Fill tidal_flux_w_m2 for planets and moons: tidal heating flux at the surface from eccentricity, spin and the body's interior (never filled in the sample) | none | done, PR #4323 |
| GEN.779 | Fill induced_moment_a_m2: the induced magnetic dipole moment of conducting oceans and cores in a varying field (never filled; every planet empty) | none | done, PR #3796 |
| GEN.780 | Fill porosity for rocky and icy bodies: the surface and crust porosity from composition and gravity (never filled) | none | done, PR #4327 |
| GEN.781 | Fill hp_ice_km: the thickness of the high-pressure ice layer under an ocean (never filled) | none | done, PR #4327 |
| GEN.782 | Fill ice_shell_km for every world with a surface ice shell, and state in the sheet rule which worlds have none | none | done, PR #3629 |
| GEN.783 | Fill shelter_column_g_cm2: the column mass of rock or ice that shields a shelter from radiation, for every world where it applies | none | done, PR #3629 |
| GEN.784 | Fill energy_flux_w_m2 and flux_variation for giant planets: internal and absorbed energy flux and its variation (never filled for classes I, J, T) | none | done, PR #4331 |
| GEN.785 | Fill star decay_age_yr and lethal_event_rate_per_gyr (never filled for any star) | none | done, PR #4105 |
| GEN.786 | Fill star state (never filled for any star): the evolutionary state of each star | none | done, PR #4105 |
| GEN.787 | Fill star disk_kind, disk_mass_solar, disk_inner_au and disk_outer_au for the stars that have a disk (never filled except for one G star) | none | done, PR #3984 |
| GEN.788 | Fill star flare_alpha, uv_c_share, log_lx_lbol and lifespan_gy for stars of classes A to F and the G and K stars where they are empty | none | done, PR #3824 |
| GEN.789 | Write the sheet rule for wide_binary_a_crit_km: not applicable to a single star, filled for a binary (348 of 522 stars empty) | none | done, PR #4317 |
| GEN.790 | Write the sheet rule for the comet columns mean_anomaly_deg, orbital_period_years and period_class (54 of 234 comets empty) | none | done, PR #4317 |
| GEN.791 | Generator switch that forces a planet or moon class so a sample can hold at least 30 of every class (random generation gives under 30 of classes O, P, T, V, E, K, L, N, F, G) | none | done, PR #3592 |
| GEN.792 | File the GEN.725 registry rows from the second-hand citations report: first-hand sources, second-hand fallbacks and the final open list | none | done, PR #3212 |
| GEN.793 | Settle two source flags in the registry: CFC-12 forcing coefficient 0.32 (IPCC TAR) against 0.33 (gap plan), and the NH3 1.078 factor against Irwin 2019 | none | done, PR #3212 |
| GEN.794 | Replace the NH3 photolysis activity factor 1.078 (GEN.501, registry stand-in 1.05 to 1.10) with the Huebner and Mukherjee 2015 quiet and active NH3 rates, or a computation from the Cheng 2006 cross section and measured solar spectra | none | done, PR #4069 |
| GEN.795 | The day and night temperature swing reads as a range: "±9.07 K" (or "±9.07 C degrees") or the full low to high temperatures, never a bare difference (bug) | none | done, PR #3487 |
| GEN.796 | Each regional temperature (polar, northern, equatorial, southern, mean and the rest) shows its actual range: the low and the high with the stated value at the middle | none | done, PR #3487 |
| GEN.797 | Register the GEN.327 young group stand-ins in the GEN.725 registry: member count weights, group radius, 2.5 km/s velocity dispersion, 300 AU remnant cap, with tier and tweak level | none | done, PR #3384 |
| GEN.798 | Design: resource locations (what counts as a valuable resource, how a belt, field, giant or icy body is scored, the sources and the stand-in tiers) | none | open |
| GEN.799 | Resource scan: compute the grade, mass and score of every resource for belts, fields, giants and icy bodies from their stored composition, as its own command and galaxy stage | none | open |
| GEN.800 | Wire the stored magnetosphere and stellar activity outputs into escape, dose and habitability (G2) | none | done, PR #4105 |
| GEN.801 | Wire the giant-planet interior outputs into the field, the heat and the description (G3) | none | open |
| GEN.802 | Wire the formation-distance water law into the atmosphere recipes b1 to b4 and reconcile it with the GEN.510 mass law (GEN.481b) | none | done, PR #3278 |
| GEN.803 | Wire the magma-ocean redox findings: the a1 width and its extension above 1 Earth mass, link a2 to a1, split a3 by building-block type (GEN.482b) | none | done, PR #3240 |
| GEN.804 | Atmosphere escape and retention read the stored XUV exposure, flare irradiation and magnetopause (G2, step 1) | none | open |
| GEN.805 | The surface dose and the habitability radiation component read the stored flare irradiation index and magnetopause (G2, step 2) | none | open |
| GEN.806 | Delete the duplicate XUV recomputation in physics/mixture_subject.py and read xuv_exposure_index (G2, step 3) | none | open |
| GEN.807 | Stellar activity (flare rate, XUV saturation) follows the star's rotation: rotation_period_hours, flare_alpha and xuv_saturated get readers (G2, step 4) | none | open |
| GEN.808 | The ozone-loss flag feeds the UV limit on life and the surface look (G2, step 5) | none | open |
| GEN.809 | A giant's dipole moment and class come from the metallic-hydrogen shell and the dilute core (G3, step 1) | none | done, PR #4436 |
| GEN.810 | A giant's internal heat reads the interior layers and the helium shell (G3, step 2) | none | done, PR #4436 |
| GEN.811 | The description of a giant says what its interior is (G3, step 3) | none | done, PR #3952 |
| GEN.812 | Greenhouse module: a water-vapour profile that dries with height, and pressure broadening checked against Venus and Earth, so an Earth twin at its own fixed point comes out near 288 K (finding 1) | none | done, PR #3551 |
| GEN.813 | Volatile inventory brackets recalibrated against the totals of the Earth, Venus, Mars and Titan, not only the ratios (finding 2) | none | done, PR #3254 |
| GEN.814 | A water-rich world (class O, 10 to 50 wt% water) gets a carbon sink: seafloor weathering and the high-pressure ice that seals the mantle (finding 3) | none | done, PR #3286 |
| GEN.815 | Sulfur sink: a sulfate and sulfide balance replaces the stand-in SULFUR_SURFACE_KEPT (finding 4) | none | done, PR #3373 |
| GEN.816 | A giant's temperature at 1 bar includes its internal heat, not only the greenhouse of the absorbed starlight (finding 5) | none | done, PR #3920 |
| GEN.817 | Rogue and starless bodies: their air goes through the same fixed point (finding 7) | none | folded into GEN.557 |
| GEN.818 | Aluminium in the melt model: an Al2O3 activity (with a tier) so a lava draw can make aluminium vapour (GEN.751 follow-up) | none | done, PR #3264 |
| GEN.819 | Wire impact_inventory (GEN.739) into the pressure supply law: late delivery and impact erosion enter the air of a rocky planet | none | open |
| GEN.820 | Residuals of the bulk iron branch table (GEN.512): the coreless branch, a slope other than 1, the 0.05 dex scatter stand-in and a larger tuning sample | none | open |
| GEN.821 | Ice-moon nitrogen inventory: ammonia ice as the nitrogen source for Titan-like moons, so the titan family derives its air instead of taking its pressure | none | open |
| GEN.822 | Wire hydro_based_mass_loss_g_s (GEN.738) into the sub-Neptune envelope rows (radius_valley, runaway_history): escape loss changes the stored radii | none | done, PR #3354 |
| GEN.823 | Wire the water-world carbon sink (GEN.814, waterworld_carbon) into the air solve and show the water-world temperatures sit inside the published range | none | done, PR #3504 |
| GEN.824 | The class O CO2 and H2O inventory bracket stops at 0.1 of the water mass, not 0.2 | none | open |
| GEN.825 | Region names: region_name() for the codec kind region, the anchor and sphere geometry of every 9-block, and a measurement of the squeezed blocks (regions chain, step 1) | none | done, PR #3375 |
| GEN.826 | Random-start and Generate "around a place" fill one anchor sphere, so a normal fill charts a whole region (regions chain, step 1c) | none | done, PR #3490 |
| GEN.827 | Region extent: the bounding sphere of a cell's generated sectors (never under the 12 pc fill sphere) and the measurement of two far-apart clumps (regions chain, step 1d) | none | done, PR #3375 |
| GEN.828 | Giant planet occurrence falls with host mass and rises with metallicity, and hot Jupiters are rare: giants per M dwarf are about 0.48 today against 0.03 to 0.06 observed, hot Jupiters 0.117 per M dwarf and 0.026 per FGK star against under 0.002 and 0.005 to 0.01 | none | done, PR #4105 |
| GEN.829 | A 2 to 20 Earth-mass population exists: sub-Neptunes with a hydrogen envelope, 0.5 to 1 planets per FGK star inside 1 AU against 0.195 today, and the mass function falls steadily above about 10 Earth masses | none | done, PR #4105 |
| GEN.830 | A class O (ocean) planet in zone e reaches a surface pressure of 9.05e9 Pa (about 90,000 bar) with the solved air: unphysical for an air above a liquid ocean (bug) | none | done, PR #3504 |
| GEN.831 | The Wolf-Rayet class W has no kind in stellar_composition, so a W star has no composition record (bug) | none | done, PR #3479 |
| GEN.832 | tests/test_gen_compact_orbit_laws.py fails 4 tests on main: a compact object with a star yields no periods (bug) | none | done, PR #3433 |
| GEN.833 | Cold icy moons get sparse ice cover because the hydrosphere's toy ocean_fraction feeds the ground draw (bug) | none | done, PR #3973 |
| GEN.834 | Non-water ice species (CO2, CH4, N2 and others) are drawn into the ground but not stored per species | none | done, PR #3973 |
| GEN.835 | The tidal heat flux (tidal_flux_w_m2) is wired into the ground draw's sulfur rule | none | open |
| GEN.836 | Moon ring hooks and rogue planets use the ring fields: has_ring for moon origin in moon_composition, and rogue planets draw rings | none | done, PR #3523 |
| GEN.837 | Forced class M and class P planets in the habitable zone carry far too much steam: about a third of 300 draws sit above 400 K because the wet volatile inventories give hundreds of bar (bug) | none | done, PR #3504 |
| GEN.838 | A separate core temperature: the interior heat model keeps one mantle temperature, so the core heat flow, inner-core growth and the dynamo cannot follow their own cooling | none | done, PR #4436 |
| GEN.839 | The hydrosphere ice shell and the surface composition volcanism gate still read the provisional-age heat, so they disagree with the final-age heat GEN.763 stores (bug) | none | done, PR #4313 |
| GEN.840 | The atmosphere recipe tables still gate each mixture row on class letters: change the gate to kind plus conditions and drop the letter from the mixture subject (the recipe-gate half of GEN.683) | none | done, PR #3729 |
| GEN.841 | thermosphere.exobase can loop for over 20 minutes in the planet atmosphere draw: a planet whose exobase never closes steps through a huge range (bug) | none | done, PR #3405 |
| GEN.842 | Follow-ups to elliptical orbits, tidal heating and planet spin: the approximations GEN.727, GEN.737 and GEN.731 left | none | open |
| GEN.843 | The stored magnetopause_rp is in Earth radii but is labelled planet radii: wrong for Jupiter and every non-Earth-sized world (bug) | none | done, PR #3783 |
| GEN.844 | A random planet's colour crashes generation: 'the spectrum has no light between 380 and 780 nm' (bug) | none | done, PR #3443 |
| GEN.845 | Redraw a stale scatter inside the bring-to-epoch run instead of refusing it (follow-up to GEN.764) | none | done, PR #3910 |
| GEN.846 | Store a rotation period for every star (stars have none, so MAP.328 leaves them standing still) | none | done, PR #3984 |
| GEN.847 | radiation.ground_dose_msv_yr passes a float rock_ppm to radiogenic.surface_content, which fails with "float object is not iterable" and turns 17 tests red on main (bug) | none | done, PR #3467 |
| GEN.848 | Store a day and night swing for every temperature the sheet shows (only the equator has one) (bug) | none | done, PR #3973 |
| GEN.849 | Re-fit the delivered-water rows b2/b3 against the bodies bulk water (the rows deliver 0.1 to 0.2 percent of mass against 0.01 to 0.03 percent of body water) | none | done, PR #3973 |
| GEN.850 | Wire radiative_convective.warming into generation behind a lookup table of the model answers, and delete the infrared_opacity.py plateau code | none | done, PR #4105 |
| GEN.851 | Wire haze into the air from aerosols.tau_haze (haze share and pressure are caller inputs today) | none | done, PR #4105 |
| GEN.852 | Infrared opacity leftovers: Venus clouds and wings (the IR3 far wing is a labelled stand-in), more gas tables, tropopause quantisation, and species without data in the cp over R path (OCS and others) | none | open |
| GEN.853 | No generated system reaches a technological civilization any more: tests/test_population.py cannot build one after 200 tries (bug) | none | done, PR #3663 |
| GEN.854 | Most class M planets do not look habitable: only a tiny share are, and most need a hardened pressure suit (bug) | none | open |
| GEN.855 | Biogenic oxygen: photosynthesis releases free oxygen into the air as the life timeline passes it, so habitable worlds can reach multicellular life (bug) | none | done, PR #3663 |
| GEN.856 | interior_heat._flux_w_m2 raises ZeroDivisionError and fails the albedo and tidal migration tests on clean main (bug) | none | done, PR #3665 |
| GEN.857 | Biogenic gases for non-oxygenic worlds: methane and haze from anoxygenic and chemoautotroph communities, donor drawdown, and registry entries for each stand-in | none | done, PR #3783 |
| GEN.858 | Stage-aware air: the atmosphere follows the world's life stage from the timeline (pre-life, microbial anoxic, oxygenation, ozone shield and its UV effect on the surface), so every planet carries the air its life type and stage would make | none | done, PR #3783 |
| GEN.859 | Microbial worlds without O2 keep a capped CO2 and CH4 air: the limit on how much CO2 and CH4 anoxic microbial life leaves in the air | none | done, PR #3783 |
| GEN.860 | Melanin radiotrophs weighted by radiation dose: the star weight times a dose factor from 0.1 at Earth dose to 1.0 at about 1,000 times Earth dose | none | dropped, absorbed by GEN.881 (Core 1, no melanin table exists in code) |
| GEN.861 | Chemoautotroph carbon base for retinal and melanin worlds: a primary producer that fixes the carbon those light-only chemistries need | none | dropped, absorbed by GEN.881 and the class store (Core 1) |
| GEN.862 | Design: the life class system, merged as docs/design/habitability/life-classes.md (nine classes, the stage ladder, the ledger, scoring and storage) | none | open |
| GEN.863 | Life class scores: every world stores its habitability per life class, shown on the object page, the API and the map filters | none | open |
| GEN.864 | Sentient life per life class: sentient species and civilizations draw from the life class scores of a world | none | open |
| GEN.865 | Headline of 9.5: sentient life draws from the life classes, so every sentient world is rated for the kind of life that built it | none | open |
| GEN.866 | A planet draw around an F5V star at an ecosphere orbit never returns when a large tidally forced moon is generated: the interior heat step size collapses (bug) | none | done, PR #3578 |
| GEN.867 | The star lists and the chemistry weights agree, and a carbon and water chemistry is never drawn on a non-water solvent | none | done, PR #3702 |
| GEN.868 | A world has a life chemistry only after abiogenesis has happened: the map icons, search facet and population code read life_stage, not life_chemical | none | done, PR #3702 |
| GEN.869 | Life descriptions and labels say what the chemistry does: Chlorophyll a builds the oxygen, Retinal is a supplement, Zinc-Bacteriochlorophyll is an aerobic acidophile, Melanin is not a pigment, and hosts without a life row are stated | none | done, PR #3702 |
| GEN.870 | Every life constant has a stand-in registry row with its tier and tweak level (LIFE_CHEMICALS gates and star weights, STAR_EVOLUTION, timelines, score thresholds, solvent factors, toxicity limits, O2 band, energy flux and shares) | none | done, PR #3581 |
| GEN.871 | Sub-ice microbial mats are scored as Flora and as edible plant food, and described with coastlines and light they do not have (bug) | none | done, PR #3702 |
| GEN.872 | The life pass draws from a per-body sub-stream, so a life-model change no longer re-rolls the names and flavour text of every later body in the system | none | done, PR #3702 |
| GEN.873 | Standing life chemistry tests: star weights match the star lists, every pace can reach complex cells inside the age of the universe, no life chemistry without a stage, no Flora under ice, no chemistry on a non-water solvent, every life constant is registered, an oxygenic world past photosynthesis has more O2 | none | done, PR #3918 |
| GEN.874 | Non-carbon and non-water life classes (silicon, hydrocarbon, sulfuric acid and ammonia solvents) in the life class system | none | folded into GEN.894 (strange-solvent life classes 5 to 9; Boss 9:11 and 9:20 PM EDT 2026-10-28) |
| GEN.875 | Design: adaptation traits drawn from a world's own stress (radiation, UV, sulfide, acid, brine, heat or cold, high CO2, low O2, toxins), each moving one named habitability threshold | none | done, PR #3583 |
| GEN.876 | Adaptation traits built: drawn in the life pass, stored in a text column, listed on the object page (adapted to: radiation, sulfide) and read by the habitability scores | none | open |
| GEN.877 | Animal oxygen floor by grade: route R1 replaces the human 8 to 50 kPa window, with a floor of about 0.1 to 0.8 kPa for sponge-grade animals rising with body size | none | done, PR #3702 |
| GEN.878 | Chemistry-specific energy: the pigment's absorption band times the star's spectrum is the multiplier on LIGHT_SHARE, instead of one bolometric curve for every metabolism | none | open |
| GEN.879 | A far-red oxygenic chemistry (chlorophyll d and f) that K and M stars can draw, so the most common stars can have oxygenic life and animals | none | open |
| GEN.880 | Chemolithoautotroph producers (H2 methanogens and acetogens, sulfur and iron oxidisers, radiolytic H2) replace Melanin as the producer under ice and in aquifers, with an oxidant budget that caps aphotic animals | none | open |
| GEN.881 | Retinal becomes a photoheterotroph supplement that needs a carbon fixer, and Melanin becomes the radiotolerance adaptation with dose in its draw, not a radiation food | none | open |
| GEN.882 | Route R2, chemosymbiotic animals: tubeworm, clam, mussel and cave-molly style communities on sulfide, methane or H2 with an oxidant, no light needed | none | done, PR #3702 |
| GEN.883 | Route R3, low-oxygen sessile animals: sponge-grade animals with Fauna capped low, on air above 0.5 percent of today's O2 or a CO2 or iron-rich ocean | none | done, PR #3702 |
| GEN.884 | Route R4, rare anoxic animals: a few percent of anoxic sulfidic worlds get animals that have lost aerobic respiration, Fauna capped low | none | done, PR #3702 |
| GEN.885 | Route R5, alternate-oxidant airs (nitrate-rich water, nitrous oxide, chlorate and perchlorate brines) for non-oxygen breathers | none | folded into GEN.887 and GEN.889 (classes 4b and 4c are Class 4 rows of the life class system; Boss 9:11 PM EDT 2026-10-28) |
| GEN.886 | The life class rule file and classifier: physics/data/life_classes.json, physics/life_classes.py, a shadow sweep that scores every body with no stored change, and a generated docs page | none | done, PR #3590 |
| GEN.887 | Life supply quantities: h2_supply, donor_supply, oxidant_supply, co2_phase and nh3_inventory, with stand-in rows and units registry entries | none | done, PR #3702 |
| GEN.888 | The acid medium model: a sulfuric acid cloud layer and acid pools from the temperature profile and the SO2 and water shares, filling the dead solvent() sulfuric branch | none | open |
| GEN.889 | The class switch: life_class, life_class_variant, life_class_best, the life_class_scores table and migration, life_chemical and the star weights and lists dropped, dominant class by net primary production | none | done, PR #3702 |
| GEN.890 | The complexity index and animal grades per class (sessile, small mobile, large active), and the Fauna rule for chemosymbiotic life (fauna never exceeds its class's producer base) | none | done, PR #3702 |
| GEN.891 | The class ledger and the class-keyed stage-aware air for Classes 2 to 9 (O2 rows come from GEN.855; the gas rows of GEN.857/GEN.858 carry class numbers from the start) | none | done, PR #3783 |
| GEN.892 | Class text: stage names and descriptions per class, and the eats-drinks-breathes-gives-off lines, replacing the chemistry text | none | done, PR #3595 |
| GEN.893 | Stand-in registry rows for every class window, supply and complexity threshold, with tier and tweak level | none | done, PR #3593 |
| GEN.894 | Strange-solvent life classes 5 to 9 (acid-world, dense CO2, cold hydrocarbon sea, ammonia-water, silicon-bearing in acid): windows, supplies, texts, stand-ins and tiers | none | done, PR #3973 |
| GEN.895 | A system draw crashes with ZeroDivisionError in interior_heat._flux_w_m2 (viscosity 0), reached through validation.reconcile_moved_planet for a moon of a moved planet (bug) | none | done, PR #3578 |
| GEN.896 | Since 4c (epoch 106) hot-zone planets and moons almost never score Yellow: the Yellow share fell to 0 for M2V, G2V and F5V (bug) | none | done, PR #3700 |
| GEN.897 | K dwarfs are 19.4% of stars against 11.6% to 12.9% observed (G 3.1% against 4.1% to 5.0%): the K class range reaches down to 3700 K (bug) | none | done, PR #3920 |
| GEN.898 | A and B stars are far too rare at the reference position: A is 0.11% of stars against 0.6% to 1.1% observed, B 0.01% against 0.04% (bug) | none | done, PR #3920 |
| GEN.899 | An M2V primary has a companion star in 15.8% of systems; the star-count table says 27% and Duchene and Kraus 2013 give 26% +/- 3% (bug) | none | done, PR #3920 |
| GEN.900 | Supernova remnants: the reference galaxy gets 260 while 294 are already catalogued (Green 2019) and the design target in tuning.py is 400 to 3000 (bug) | none | done, PR #3920 |
| GEN.901 | Magnetars: the reference galaxy gets 22 while 26 are known (McGill catalog) and tuning.py (GEN.278) aims at about 100 active (30 to 600) (bug) | none | done, PR #3920 |
| GEN.902 | Wolf-Rayet stars: the reference galaxy gets about 590 while Rosslowe and Crowther 2015 estimate 1200 (+300, -100) (bug) | none | done, PR #3920 |
| GEN.903 | Atmosphere possibility table e: rows e1 and e2 (bare rock, retained H/He envelope) are never reached for class I because group b has only row b4 for class I (low priority) (bug) | none | done, PR #3729 |
| GEN.904 | Daylight columns stay stale or empty after an orbit move: re-run daylight after generate_spin in the reconcile (bug) | none | done, PR #3920 |
| GEN.905 | Mantle redox and delta-IW are empty on class V ocean worlds (bug) | none | done, PR #3920 |
| GEN.906 | Icy-moon tide is stored as 0.0 where it should be NULL (bug) | none | done, PR #3783 |
| GEN.907 | Life census: only 0.08% of life worlds reach animals against the 15 to 35% band, because almost all life worlds are Melanin or Retinal and few ordinary worlds are oxygenic (bug) | none | done, diagnosed and decided (Boss 6:04 AM EDT 10-29: both in order, restate the bands first, then re-baseline after GEN.867); the work is TEST.435 |
| GEN.908 | Worlds labelled surface ocean are mostly frozen: at least half have a mean surface temperature below 273 K (bug) | none | done, PR #3920 |
| GEN.909 | Temperate liquid-water worlds (273 to 323 K, 10 kPa to 1 MPa, ocean) are rare: 10 in 120 systems at G2V, 6 at K5V, 0 at M2V (bug) | none | done, diagnosis only (cold by the physics, no code change; report gen909-cold-ecosphere-diagnosis-2026-10-29.md) |
| GEN.910 | K5V worlds are also cold in the inner third of the habitable zone (9 of 13 under 273 K): split the cold share by flux and by CO2 to find out why (bug) | none | done, diagnosis only (no K5V-specific fault; report gen910-k5v-cold-share-2026-10-29.md) |
| GEN.911 | Habitable-world requests and the habitable zone: the zone mostly holds frozen worlds, so decide where habitable requests aim (inner half of the zone or a temperate band from the solved climate) | none | done, PR #3973 |
| GEN.912 | The planet albedo draw does not follow the star's colour, so planets of cool stars come out a little too cold (stand-in) | none | done, PR #4385 |
| GEN.913 | A civilization is never forced: the intelligent_life request and the prevalence setting only tweak probabilities; remove forced placement, redraws and squeezed timelines (bug) | none | dropped (withdrawn: Boss never wants a civilization forced) |
| GEN.914 | G dwarfs are 3.1% of stars against 4.1% to 5.0% observed: the mass function from 0.5 to 1 solar mass (Kroupa slope 2.3) gives too few stars of 0.86 to 1.1 solar masses (bug) | none | superseded by GEN.922 (the whole star mass function chain, Boss 5:24 PM EDT 2026-10-29) |
| GEN.915 | Life depends only on the conditions, never on the planet class letter: every combination of conditions that can host life can evolve it on any class | none | open |
| GEN.916 | The habitable-world guarantee is judged on the conditions (a world that can hold life), not on the class letter M | none | done, PR #3973 |
| GEN.917 | The dominant surface material and the aerosol draw are chosen from the conditions, not from the class letter | none | done, PR #3920 |
| GEN.918 | The galaxy's star count follows its own star-formation history: the density model counts systems per parsec cubed while the history counts stars formed, so the galaxy holds 0.46 of the stars its history implies (root of GEN.900, GEN.901, GEN.902) | none | open |
| GEN.919 | White dwarfs are 8.1% of stars against 4.5% to 7% expected once the low-mass slope is 1.1; their share follows the old-population age spread, not the slope (bug) | none | absorbed by GEN.925 (white dwarfs fall out of the present-day populations) |
| GEN.920 | Giants around F, G and K stars are too rare after the class cut-over: 0.03 cold and 0.07 in all per star against 0.05 to 0.25 and 0.08 to 0.35 observed; retune the mass-law gas share (bug) | none | done, PR #4105 |
| GEN.921 | Moons per planet rose after the class cut-over (G2V 1.90 against 1.25): find out why and whether it is right (bug) | none | done (diagnosed by Bug 1, no defect: the planet mix and the 1e20 kg moon floor; report /mnt/project-files/reports/bug1/moons-per-planet-gen921-2026-10-30.md; the rocky-planet rate is GEN.962) |
| GEN.922 | The star mass function is built right: initial mass function and present-day populations match the observed curve segment by segment, for every mass range and every star-like object type (parent, replaces the one-number GEN.914 tweak) | none | open |
| GEN.923 | Observed star-population targets, with sources, for every mass segment and star-like object type (the table the star mass function must match) | none | dropped, replaced by TEST.454 and GEN.926 to GEN.932 (Bug 1 plan) |
| GEN.924 | A segmented initial mass function with sourced breakpoints and slopes replaces the single low-mass slope, for the whole mass range down to the brown dwarfs and up to the most massive stars | none | dropped, replaced by GEN.927 and GEN.928 (Bug 1 plan) |
| GEN.925 | Present-day star populations follow from the initial mass function, the star-formation history, main-sequence lifetimes and remnant masses: white dwarfs, neutron stars and black holes fall out, not tuned | none | dropped, replaced by GEN.929 and GEN.930 (Bug 1 plan) |
| GEN.926 | One initial mass function module (physics/imf.py) replaces the six copies of the power-law code, with no change of output | none | done, PR #3760 |
| GEN.927 | Single-star mass function rows: a bell below 1 solar mass, the field slope from 1 to 8, the sub-stellar tail and the massive end, with stand-in tiers and tweak levels and a design note | none | done, PR #3984 |
| GEN.928 | System primaries follow from the single-star curve plus multiplicity (a top-down solve), with brown dwarf companions and the brown dwarf desert; the star-count table agrees | none | done, PR #3984 |
| GEN.929 | Systems per stellar mass and the galaxy's star count come from the same construction (the mean system mass of item 4); supersedes GEN.918 | none | done, PR #3984 |
| GEN.930 | Remnants (white dwarfs, neutron stars, black holes) and giants follow from the curve and are checked against it; the pinned derived-density tests are re-derived, not loosened | none | done, PR #3984 |
| GEN.931 | The spiral-arm contrast at the Sun's place is set per stellar population (young arm populations strong, A stars weak), not by the single inter-arm minimum (the A-star lever of GEN.898) | none | done, PR #4105 |
| GEN.932 | The star mass function epoch: golden, design docs, census targets and tags, retire the old tags; items 3 to 7 land as one merge in one lock turn | none | open |
| GEN.933 | One mass-boundary table (brown dwarf/star, planet/brown dwarf, spectral-class mass edges, Chandrasekhar, neutron star/black hole) with sources and tiers, read by the mass function, companions, rogue bodies and the display ladders | none | done, PR #3777 |
| GEN.934 | Replace the hand-written mass conversions in about 47 modules with the unit registry helper (do not rename _sol/_solar) | none | done, PR #4281 |
| GEN.935 | Evolved stars from the real post-main-sequence tracks: read the MIST post-main-sequence rows at setup for subgiant, giant, core-helium-burning, asymptotic giant and supergiant luminosity, temperature, radius and duration; red supergiants stop at about 300,000 Lsun; a Cepheid and yellow-star flag | none | done, PR #3984 |
| GEN.936 | The 0.08 to 0.1 solar mass main sequence comes from the Baraffe 2015 table instead of an extrapolation, and the star/brown-dwarf boundary depends on [Fe/H] (no gap and no overlap) | none | done, PR #4105 |
| GEN.937 | Neutron star masses from the sourced Gaussians in tuning (two populations) and black hole masses from a progenitor-mass and metallicity law | none | absorbed by GEN.930 (remnant masses, Bug 1 triage 5:47 PM EDT 2026-10-29) |
| GEN.938 | One white dwarf initial-final mass relation, read by the star model, the multi-star roster, the disk history and mass_loss | none | absorbed by GEN.930 (one white dwarf initial-final relation, Bug 1 triage) |
| GEN.939 | Wolf-Rayet luminosity follows the present mass by the helium-star mass-luminosity relation (low priority) | none | done, PR #3984 |
| GEN.940 | Present mass for giants, horizontal-branch, asymptotic and massive stars: a Reimers-type loss for the red giant and asymptotic branches and a wind loss for massive stars | none | done, PR #3984 |
| GEN.941 | Forced star types draw their star from the model; the legacy letter tables and the dead random branch are removed; the class reference pages read the new ranges | none | done, PR #3984 |
| GEN.942 | Supported life paces come from the star's actual lifespan; the per-letter lifespan table (tuning.STAR_EVOLUTION) is removed | none | done, PR #4321 |
| GEN.943 | Spectral subclasses from a letter-and-subclass temperature lookup table (Pecaut and Mamajek style), the K/M boundary at about 3,900 K, and a decision on Yerkes class VI (subdwarf) | none | done, PR #3984 |
| GEN.944 | Later: an age-metallicity relation inside each population and a 10 to 13 Gy thick-disk window; variable-star flags (Cepheids, Mira and long-period variables, carbon and S stars) | none | open |
| GEN.945 | Audit every numeric constant the generator uses against the data we already download or hold: which are registered stand-ins, which come from a table, and which are unregistered hand-tuned numbers; ranked findings with a data source named for each | none | open |
| GEN.946 | Burn down the hand-tuned constants: replace each L3 or L4 registry row and each unregistered constant the audit ranks high with the data we hold, or record why no data exists | none | done, PRs #4319 and #4320; rest kept as no data |
| GEN.947 | Pre-main-sequence contraction times are calibrated to MIST (about 2 times too short at 1.5 solar masses and up and at 0.2 solar masses) (low priority) | none | done, PR #3984 |
| GEN.948 | Brown dwarf radius from an age-dependent table (Baraffe 2015) instead of nine hard-coded points | none | done, PR #4105 |
| GEN.949 | Wire the held methane absorption tables (AD1 Karkoschka k-tables, PDS, NMSU) into the visible colour of giant planets | none | open |
| GEN.950 | Water share of solids by distance from Lodders 2003 Table 11 and its condensation temperatures, replacing the fixed 170 K snow line (low priority) | none | done, PR #4105 |
| GEN.951 | Companion mass ratio depends on the primary mass and the orbital period, with the twin excess (Moe and Di Stefano 2017 Table 13) instead of a flat 0.1 to 1 (low priority) | none | open |
| GEN.952 | Store each scatter member's radius (scatter_members.radius_au) so a page can show the X-ray luminosity of an accreting pair without guessing | none | done, PR #3984 |
| GEN.953 | Birth spin retune: a scale on the rocky planet birth rotation so free days fall in the Solar System range | none | done, PR #3973 |
| GEN.954 | Tides stage 1a (lock-free): a rheology module giving k2 and Q from temperature, melt fraction and tide frequency | none | done, PR #3799 |
| GEN.955 | Tides stage 1b: wire the rheology into tidal_heating, interior_heat and circularisation | none | done, PR #4321 |
| GEN.956 | Tides stage 2a (lock-free): a torque-sum module for planet spin: star, moon and atmospheric thermal torques with floors and an attractor | none | done, PR #3801 |
| GEN.957 | Tides stage 2b: wire the torque sum into planet spin, with the moon-on-planet tide after the moons exist | none | done, PR #4321 |
| GEN.958 | Tides stage 3a (lock-free): a moon orbit evolution module: outward drift, spiral-in and the Roche fall | none | done, PR #3807 |
| GEN.959 | Tides stage 3b: wire moon orbit history, rings from fallen moons and the moons' tide into the planet's heat | none | done, PR #4321 |
| GEN.960 | Tides side thread: ice rheology for the hydrosphere ice shell and the icy branch of tidal_heating (c-ice) | none | done, PR #4321 |
| GEN.961 | Settle giants: after the star's age is final, recompute each giant's radius and redraw its interior, at the mean light of its real orbit (eccentricity-first folded in) | none | done, PR #3984 |
| GEN.962 | Rocky planets hold 1.5 to 2.0 moons of 1e20 kg and over each (the Solar System has 1 for 4 rocky planets): ground the moon rate in giant-impact and capture frequency | none | open |
| GEN.963 | The envelope-survival rule and the Hycean flag on the existing sub-Neptune classes (the part of GEN.350 that stays; no new class letter) | none | done, PR #4156 |
| GEN.964 | H-He equation of state tables (a): a setup catalogue entry and a reader with tests on Jupiter, Saturn, Uranus and Neptune | none | done, PR #4436 |
| GEN.965 | H-He equation of state tables (b): replace the n=1 polytrope in the giant interior solver with the tables | none | done, PR #4436 |
| GEN.966 | Dynamo field (a): a pure module giving field strength and lifetime from core heat flow and spin, with tests on Earth, Mercury, Ganymede, Jupiter and Mars (dead) | none | done, PR #3802 |
| GEN.967 | Dynamo field (b): the stored field strength, moment and lifetime come from the dynamo module and feed the radiation shielding | none | done, PR #4436 |
| GEN.968 | Neutron star structure (a): a TOV integration over the equation-of-state tables held on disk, tested against NICER radii | none | open |
| GEN.969 | Neutron star structure (b): wire the integrated radius and layers into the neutron star model, and add helium-core and oxygen-neon white dwarf layers | none | done, PR #3984 |
| GEN.970 | Stellar spectra (a): a spectral library entry (BT-Settl, MUSCLES for M dwarfs) and a lookup by temperature, gravity and activity, with tests on the Sun and an M dwarf | none | open |
| GEN.971 | Stellar spectra (b): the UV band shares and the map colour of a star come from the library spectrum, not from a blackbody with a seeded UV-C share | none | done, PR #4105 |
| GEN.972 | Climate feedback 1: ice-albedo on the surface temperature fixed point (ice cover grows as the planet cools and raises the albedo) | none | open |
| GEN.973 | Climate feedback 2: water vapour on the fixed point (a warmer ocean world holds more vapour and warms further, up to the runaway limit of GEN.639) | none | done, PR #4370 |
| GEN.974 | Climate feedback 3: cloud albedo on the fixed point (cloud cover changes the reflected share) | none | open |
| GEN.975 | Climate sensitivity sweep test: the fixed point converges on a cooling and a warming sweep and stays inside the published range of Earth's climate sensitivity | none | done, PR #3895 |
| GEN.976 | A giant's radius is read at the star's age before the age is moved: 28 percent of table-sized giants have a stored radius over 1 percent too large, up to 9.3 percent (bug) | none | done, PR #3984 |
| GEN.977 | Guard test: a seeded sweep fails if the generator ever makes a table-sized giant with light over 300 times Earth's and eccentricity over 0.3 | none | done, PR #3803 |
| GEN.978 | Birth spin derivation (a): a module for the last-impact mass ratio gamma by planet mass, from the N-body impact statistics | none | done, PR #3809 |
| GEN.979 | Birth spin derivation (b): the birth rotation of a rocky planet follows Kokubo and Genda eq. 4 with the gamma law, replacing the scale factor of GEN.953 | none | done, PR #4321 |
| GEN.980 | The stored class of 4 to 8 bodies differs from the class their final body gives (stored Q, conditions give D; stored P or K, derived D or P): the class is derived before the air is final (bug) | none | done, PR #3973 |
| GEN.981 | Induced magnetic moment part (b): conducting cores and planets in the star's field (the part of GEN.779 that remains) | none | done, PR #4436 |
| GEN.982 | Wire the photolysis lifetimes (physics/photolysis.py and photolysis_lifetimes.json) into the atmosphere chemistry | none | open |
| GEN.983 | Wire the radiation belts (physics/radiation_belts.py) into the surface dose and the habitability radiation score | none | done, PR #4436 |
| GEN.984 | Make the star's composition draw independent of its age so every planet reads the final composition (bug) | none | done, PR #3984 |
| GEN.985 | Circularise a close-in planet's eccentricity again at the final star age | none | done, PR #3984 |
| GEN.986 | Measure how far the rocky planets' dose, surface ages and magnetism move when read at the final star age | none | done, measured in reports/gen986-final-age-rocky-planets-2026-10-30.md |
| GEN.987 | Read a moon's tidal heat from the settled giant's radius | none | done, PR #3984 |
| GEN.988 | Give Ca, S and Ar a star-to-star scatter like the other generated [X/Fe] abundances (bug) | none | done, PR #3984 |
| GEN.989 | Lower generated [Mg/Fe] by 0.05 to 0.07 dex for [Fe/H] above -0.2 to match Brewer and Fischer 2016 (bug) | none | done, PR #3984 |
| GEN.990 | Raise the C, O and Si abundance scatter, now about half the Brewer and Fischer 2016 scatter (bug) | none | done, PR #3984 |
| GEN.991 | Raise generated [O/Fe] by 0.03 to 0.09 dex to the Brewer and Fischer 2016 table in every [Fe/H] bin | none | done, PR #3984 |
| GEN.992 | Re-finish every rocky planet and moon after the star's age is final, then re-judge the forced requests with the retry round | none | done, PR #4105 |
| GEN.993 | Widen the core light-element total in bulk_draw to 0.02 to 0.25 (0.25 dex) instead of the uniform 0.03 to 0.12 | none | open |
| GEN.994 | Wire the hot-giant regime edges (equilibrium_edges) into the atmosphere row for hot giants | none | open |
| GEN.995 | Draw the binary-stripped Wolf-Rayet channel so the Galaxy Wolf-Rayet count reaches the observed range (bug) | none | open |
| GEN.996 | Draw binary-made low-mass white dwarfs so the white dwarf mean mass reaches the accepted range (bug) | none | open |
| GEN.997 | A moon of a planet around a very young star (under 0.06 Gy) holds a life class: test_planets_respect_their_stars_age_and_history is red on main (bug) | none | done, PR #3887 |
| GEN.998 | The default galaxy disc gas depletion time is 0.85 Gy against the 4.2 Gy the Galaxy measures: check the value and its source (bug) | none | open |
| GEN.999 | Give the disc gas law a low-density threshold so a Milky Way mass disc depletes its gas in about 4.2 Gy, not 1.46 Gy | none | open |
| GEN.1000 | Draw the largest body of every belt from the small body size law (replaces the stand-in largest-body size in the GEN.477 belt mass) | none | open |
| GEN.1001 | A world whose star is younger than the first step of its life schedule is stored with a life class but no life stage, so the map markers, the object page line and the class filters show life that has not started (bug) | none | done, PR #3984 |
| GEN.1002 | Registered table rows that rest on a Solar System body are flagged stand-in with a replacement item, or re-sourced | none | done, PR #4038 |
| GEN.1003 | The gas equilibrium solver can stall with no convergence on a hot hydrogen-oxygen mixture (EquilibriumError at 725.97 K and 1.0 bar) (bug) | none | done, PR #4105 |
| GEN.1004 | A triple system whose third star cannot be found: it may sit inside another star or on an orbit that is not drawn (bug) | none | done, not reproducible on main; fixed by #4026 and MAP.255 (Aux 7) |
| GEN.1005 | Replace the Solar System-anchored redox, dry carbon and comet species rows with models that carry the physics (replaces the Solar System-resting rows of GEN.1002 that GEN.474 does not cover) | none | done, PR #4436 |
| GEN.1006 | Stars and bodies in one system never overlap by design: every member knows the others when any is rebuilt, then fit and adjust afterwards | none | open |
| GEN.1007 | Overlap-free stars, part a: pass the primary's [Fe/H] into the star tree (the tree evolves companions at [Fe/H] 0 while the build uses the real value) | none | dropped (folded into GEN.1006 P3 (MemberRecord), design 10:56 AM EDT 2026-10-30) |
| GEN.1008 | Overlap-free stars, part b: store a giant's radius draw on the member and read it at build, so the radius is drawn once | none | dropped (folded into GEN.1006 P3 (MemberRecord), design 10:56 AM EDT 2026-10-30) |
| GEN.1009 | Overlap-free stars P3: evolve each star once into a frozen MemberRecord that the tree, the build and the planets all read (replaces GEN.1007 and GEN.1008) | none | done, PR #4321 |
| GEN.1010 | Overlap-free stars P4: one keep-out rule at every depth of the star tree, replacing the Roche-or-extent split | none | done, PR #4321 |
| GEN.1011 | Overlap-free stars P5: truncated sampling of period and eccentricity, deleting the 60-try loop and the tree restart loop | none | done, PR #4321 |
| GEN.1012 | Overlap-free stars P6: a general history rule replacing the white-dwarf-only engulfment test in classify_pair | none | done, PR #4385 |
| GEN.1013 | Overlap-free stars P7: keep compact members and brown dwarfs in the tree for the orbit draw | none | done, PR #4385 |
| GEN.1014 | Overlap-free stars P8: a clearance ledger and planet, moon and belt sampling from it | none | open |
| GEN.1015 | Overlap-free stars P9: delete the GEN.1004 quick fix (built_tree_clears, the clearance branch of _companions_fit, _built_stars_clear), no shim | none | open |
| GEN.1016 | Overlap-free stars P10: per-level hierarchical stability margin and the 2+2 quadruple fit, if the verifier shows it is needed | none | open |
| GEN.1017 | Overlap-free stars P0: give the star tree the primary's [Fe/H] (and the build's age and mass) so the tree and the build evolve each companion the same way (replaces the dropped GEN.1007) | none | open |
| GEN.1018 | Overlap-free stars research part: port a 3+1 quadruple stability criterion and verify the sources from the papers (GEN.1006) | none | done, PR #4043 |
| GEN.1019 | Venus clouds and the CO2 far-wing line shape in the infrared model (the IR3 far wing is a labelled stand-in) | none | open |
| GEN.1020 | Infrared opacity of the acid cloud and the SO2, CO and H2O line wings, so a Venus-like air reaches its published surface temperature once the sourced CO2 far wing replaces IR3 | none | done, PR #4061 |
| GEN.1021 | Replace the 67P calibration of the comet species factors with a second comet or a laboratory trapping table | none | open |
| GEN.1022 | One latent ice formation temperature per comet drives its CO/H2O, N2 and Ar together (option B of the comet species model) | none | open |
| GEN.1023 | Two planets of a close pair in triple seed 208 sit 9.93 mutual Hill radii apart where the generator promises 10 (bug) | none | open |
| GEN.1024 | GEN.744 part B: the alpha-enhanced track shift through evolve_star and the scatter tables | none | done, PR #4436 |
| GEN.1025 | GEN.944 part 1b: the 10 to 13 Gyr thick-disc window in the galaxy density and age-law tables | none | done, PR #4436 |
| GEN.1026 | Carbon star and S star flags on the star sheet from the thermal-pulse and dredge-up count of the post-main-sequence track | none | done, PR #4399 |
| GEN.1027 | Ap and Am stars: a chemically peculiar branch of the A-star draw with rotation and magnetic field | none | done, PR #4436 |
| GEN.1028 | The Mardling and Aarseth 2001 stability criterion holds only for an outer mass ratio q_out of 5 or less and the code does not enforce it (bug) | none | done, PR #4276 |
| GEN.1029 | Research: port a 3+1 quadruple stability criterion without a machine-learning dependency (Vynatheya, Mardling and Hamers 2023 model files) | none | done, PR #4098 |
| GEN.1030 | Methane pseudo-adiabat in the giant planet column, so Uranus-like and Neptune-like one-bar temperatures match Voyager (GEN.852 tropopause and condensables family) | none | done, PR #4436 |
| GEN.1031 | Port the 2+2 quadruple stability network of Vynatheya, Mardling and Hamers 2023 the same way as the 3+1 one (about 8,300 numbers, plain numpy forward pass, no machine-learning dependency) | none | done, PR #4129 |
| GEN.1032 | Refit the greenhouse table on the far-wing and acid-cloud model that GEN.1019 and GEN.1020 landed, and move the tropopause-jump sweep (GEN.852 family) | none | done, PR #4436 |
| GEN.1033 | Ice-giant internal heat: the law goes through Jupiter and Neptune, so a Uranus-mass giant gets Neptune's 0.33 W/m2 instead of Uranus's 0.04 (the Uranus anomaly), and giants draw their Bond albedo from the rock range 0.12 to 0.35 | none | open |
| GEN.1034 | Ocean tide data: the open-ocean equilibrium tide height and range at a world's surface from each moon and the star, for worlds with a surface ocean | none | done, PR #4385 |
| GEN.1035 | Surface liquid inventory: which liquid a cold or hot world holds as a sea (methane-ethane-nitrogen, liquid nitrogen, liquid CO2, ammonia-water, sulfur) from its temperature, pressure, ice and air, so non-water oceans can have tides | none | done, PR #4321 |
| GEN.1036 | Oddity table and loader: physics/oddities.py, physics/data/oddities.json, the row schema (requires, rate with source, tier and tweak, named input effect, limits, text) and the GEN.725 registry rows, no wiring | none | done, PR #4126 |
| GEN.1037 | Confirm the oddity rates and thresholds marked [R] (W3, W4, W5, W10, W14, W15, W20 first) against papers, write the literature table with sources, register the stand-ins | none | done, PR #4145 |
| GEN.1038 | W2 module: the Venus-type retrograde slow rotator, the spin equilibrium of a thick-air planet from the thermal tide against core-mantle friction (Correia and Laskar 2001, 2003), pure function | none | done, PR #4133 |
| GEN.1039 | W1 module: a seasonal temperature swing from the stored obliquity in the regional temperature (tilt-driven climate), polar summer hotter than the equator, ice and wet-bulb limits by season | none | done, PR #4141 |
| GEN.1040 | W10 module: collisional mantle stripping to an iron-rich (Mercury type) core fraction, pure function (Marcus 2010, Benz 2007 scaling) | none | done, PR #4152 |
| GEN.1041 | W14 module: mass loss and lifetime of a disintegrating rocky planet (Perez-Becker and Chiang 2013), pure function | none | done, PR #4128 |
| GEN.1042 | W15 module: super-puff radius and its fade with age (Lee and Chiang 2016, Millholland 2019), pure function | none | done, PR #4143 |
| GEN.1043 | W5 module: detect resonant chains from the stored periods and compute the forced eccentricity and the resulting tidal heating (Goldreich and Soter, Murray and Dermott), pure function | none | done, PR #4144 |
| GEN.1044 | W19 module: the nightside collapse threshold of a thin-air locked planet (Wordsworth 2011, Leconte 2013), pure function | none | done, PR #4131 |
| GEN.1045 | W16 audit: does envelope_escape fractionate helium (hydrogen escapes first, leaving a helium-dominated air; Hu 2015, Malsky and Rogers 2020)? Add the fractionation if not | none | done, PR #4146 |
| GEN.1046 | W20 design: the two-branch climate output (snowball or warm, runaway or temperate) from the ice-albedo feedback and the CO2 exit rule, as a history draw between branches | none | done, PR #4290 |
| GEN.1047 | Wire the oddity table and the weirdness dial: epoch, stored oddities text column and numeric columns (db/column_units.py), the "Odd: ..." line with its cause sentence on the object page, the API field, an admin force-oddity like FORCE_CLASS, the dial in settings (default 1, 0 turns all off) | none | open |
| GEN.1048 | Wire W1 (tilt-driven seasonal climate), W2 (Venus-type spin), W10 (iron-rich cores) and W15 (super-puffs) into the bulk, spin and climate inputs | none | open |
| GEN.1049 | Wire W5 resonance heating and the W3 high-eccentricity and W4 misaligned-orbit tails in planet_architecture | none | open |
| GEN.1050 | Wire W14 (the dying-planet flag) and W9 (rings from decayed moons) | none | open |
| GEN.1051 | Wire the W19 nightside collapse and the W20 bistable branches (the history draw between branches) | none | open |
| GEN.1052 | Heavy-liquid layers and mercury (from GEN.1034): the deep global heavy-liquid path with the fluid Love number, the thin-layer refusal, and mercury and other liquid metals as rows of the surface liquid table | none | open |
| GEN.1053 | G1 module: the giant-impact tilt of an ice giant (Uranus 98 degrees), a tail row on the spin lever with its ring and moon plane following the equator (Slattery 1992, Morbidelli 2012, Kegerreis 2018), pure function | none | done, PR #4149 |
| GEN.1054 | Ring origin from history: a ring comes from a moon decayed inside the Roche limit (R1, Canup 2010, Hesselbrock and Minton 2017) or a comet passing inside it (R2, Dones 1991), pure function; the always-a-ring rule for giants is what it replaces | none | done, PR #4124 |
| GEN.1055 | Ring age and pollution darkening: a ring stores its age and its ice share and darkness follow from pollution over time (R6); ring shadow climate hook (R8) | none | open |
| GEN.1056 | Resonance rings of small bodies: rings beyond the Roche limit held by a spin-orbit resonance with a satellite or fast spin (Quaoar 7.4 radii, Haumea, Chariklo, Chiron; R4), pure function, plus R9 (no ring inside about 0.1 AU) | none | done, PR #4127 |
| GEN.1057 | Wire conditional rings: a ring follows from moon decay or comet passage and a giant with neither has none, so Saturn-like broad rings are the tail (epoch, replaces the BC20 draw), plus ring age, R4 and R9 | none | done, PR #4385 |
| GEN.1058 | Wire the G1 giant-impact tilt and the ring plane following the equator | none | open |
| GEN.1059 | The GEN.412 inversion of a stored giant row ignores the giant's own heat, so its stored albedo and greenhouse factor no longer reproduce its temperature from starlight alone (bug) | none | open |
| GEN.1060 | Belt research read: open the nine sources of design section 8 (Bottke 2006, Wisdom 1980 and 1982, Moons and Morbidelli 1995, Walsh 2011, Raymond and Izidoro 2017, Stern and Levison 2002, Soter 2006, Margot 2015, Kenyon and Bromley 2008, Debes and Sigurdsson 2002, Zuckerman 2010, Koester 2014, Ida and Guillot 2016, Gladman and Volk 2021), one register row per number used | none | done, PR #4180 |
| GEN.1061 | physics/belt_ledger.py: the size-law top-end order statistic, the individual-object draw by the cut rule with mass conserved, and the Dc iteration with belt_mass | none | done, PR #4180 |
| GEN.1062 | physics/resonance_gaps.py: the resonance and chaotic-zone map from the giants' orbits and the size-dependent Yarkovsky cut | none | done, PR #4217 |
| GEN.1063 | physics/body_status.py: roundness mass by composition, Soter mu, Stern-Levison Lambda, Margot Pi, and the status words (planet, dwarf planet, minor body, moon), both words and numbers stored | none | done, PR #4193 |
| GEN.1064 | physics/belt_history.py: early-depletion survival, instability sweep, and the star events (engulf partly, mass-loss destabilise share, companion truncation, white dwarf pollution link) | none | done, PR #4224 |
| GEN.1065 | Belt member designation: a pure function designation(belt_id, rank) printing the belt ID plus the rank in hex (no name, no salad name, no codec name), and its round trip through parse_id | none | done, PR #4161 |
| GEN.1066 | Design the system body table or view (system_bodies: id, parent, kind, mass, radius, density, a, e, i, status, tier) and the barycentre that includes the belts; decide view versus table | none | open |
| GEN.1067 | Per-member composition from the finished taxon kernel at the body's own a, differentiated bodies only if the 26Al ratio says so, density from the compaction taper | none | done, PR #4220 |
| GEN.1068 | Belt birth from the architecture: a belt is a gap or an unfinished tail in the planet chain, not a 10 percent slot roll, and sees the giants (GEN.728 slot loop, system.py) | none | open |
| GEN.1069 | Store the belt mass, Dc, age history text and formation snow line (DB.70 column, migration) and wire belt_mass into AsteroidBelt, replacing the three density words | none | open |
| GEN.1070 | Store the individual belt members (rows with parent link and rank), the designation display, the system map drawing them as dots at true orbits, and a body page per member; API version bump | none | open |
| GEN.1071 | Outer (Kuiper-type) disc built as a belt from the same chain, linked to the debris descriptor so the descriptor and the belt are one thing | none | open |
| GEN.1072 | White dwarf and giant-branch belt events wired into remnant_orbits and the engulfment step, replacing the uniform expansion factor of belts | none | open |
| GEN.1073 | Bring-to-epoch for belts: grinding recomputed from the stored age, members above the grinding scale kept | none | open |
| GEN.1074 | Captured asteroid moons: a giant near a belt captures bodies from it as irregular moons; research the capture mechanisms and rates, then plan the build | none | dropped (research done in the captured asteroid moons draft; the build items are GEN.1075 to GEN.1087) |
| GEN.1075 | Captured moons read list and register rows: Nesvorny et al. 2003 and 2007, Hamilton and Burns 1991, Agnor and Hamilton 2006, Bottke 2010, Sheppard and Jewitt 2003 and 2005, Granvik 2012, Rosenblatt 2016 | none | done, PR #4175 |
| GEN.1076 | physics/irregular_zone.py: Hill radius, stable limits (prograde and retrograde), the inner edge from the outermost regular moon, the Kozai exclusion | none | done, PR #4155 |
| GEN.1077 | physics/irregular_population.py: capture windows, the cumulative count N0(>R), grinding, and the order-statistic draw into individual, on-demand and population tiers | none | done, PR #4239 |
| GEN.1078 | physics/irregular_orbits.py: a/r_H, e, i, survival and families (parent breakup with a velocity spread) | none | done, PR #4164 |
| GEN.1079 | Collisional-saturation check for the irregular zone: reuse belt_mass_earth and show the invariance across the four giants, or replace the guess with a fitted scalar per planet registered at L3 | none | done, PR #4245 |
| GEN.1080 | Composition of captured bodies from the reservoir (the GEN.473 taxon kernel at the source belt or disc) and a heterogeneous parent | none | done, PR #4248 |
| GEN.1081 | Captured-origin weight depends on moon mass: the ORIGIN_TABLE row `captured` becomes a function of mass so a 1e20 kg or heavier moon is rarely captured, as a pure helper | none | done, PR #4184 |
| GEN.1082 | The architecture records its capture windows (the instability event, the number of giant-giant close encounters, the gas-disc lifetime already kept in DiskSystem) as a field on the system | none | open |
| GEN.1083 | Wire irregular moons into generate_moons as a second tier above the regular moons (MOON_MIN_MASS_KG left alone for the regular list) with the captured-origin weight | none | open |
| GEN.1084 | Store the individual irregular moons and the new moon columns (migration, API version), and call the stability test for irregulars in drop_unstable_moons and in bring-to-epoch when a planet moves | none | open |
| GEN.1085 | Belt ledger link: a captured member leaves its belt's mass total; until then captured bodies come from the outer-disc reservoir only | none | open |
| GEN.1086 | System map and pages for irregular moons (eccentric, inclined, retrograde orbits on the moon-orbit view, a body page per irregular) and the clearance ledger taking the irregular zone as occupied space | none | open |
| GEN.1087 | Census of Class D and dwarf-size bodies, moons of 50 to 500 km, and belts per system from the current generator, per star type | none | done, PR #4185 |
| GEN.1088 | physics/belt_origin.py: the cause table (gap, unfinished tail, snow-line bump, giant stirring, instability depletion, collisional debris, destroyed planet), the survival factors, the origin and formation record and the one-sentence text | none | done, PR #4168 |
| GEN.1089 | physics/scatter.py: the scatter step as pure functions of the giants and the disc (a rule over the architecture's events with stored mass shares, no N-body run) | none | done, PR #4160 |
| GEN.1090 | physics/dwarf_population.py: N(>D) from belt mass and size law, the order-statistic top body, the dwarf-planet line and Class D tagging | none | done, PR #4196 |
| GEN.1091 | Rewrite the Class D wording in tuning.PLANET_CLASSES["D"] and the docs to match the conditions rule plus the ledger scope (docs only) | none | done, PR #4188 |
| GEN.1092 | Wire the scatter phase and the belt origin: scatter tag on stored members, the outer disc fed from the scatter step, origin and formation shown on the belt page | none | open |
| GEN.1093 | physics/belt_density.py: the belt density factor (objects per square AU above a stated diameter, plus mass per square AU) from the size law, the belt annulus and the ledger | none | done, PR #4166 |
| GEN.1094 | Read Goldberg and Batygin 2023 for the librating share of near-resonant pairs (the share of pairs that are heating-capable, W5), and register the number | none | done, PR #4153 |
| GEN.1095 | Confirm the remaining [R] oddity rates (W6, W7, W9, W12, W16, W21 to W26, G1, G7) against papers and register the stand-ins | none | done, PR #4183 |
| GEN.1096 | planet_mass.gas_outcome returns max(runaway, core) and the core is capped at 0.01 of the star, not at the 13 Jupiter-mass edge: a 1.5 to 2 solar-mass host at twice solar metals draws 15 to 21 Jupiter masses (4 in 3,000 draws) (bug) | none | done, PR #4321 |
| GEN.1097 | Every object has a day-and-night state and counter: solar day length, day and night hours, and an explicit state for tidally locked, non-rotating, polar-day and no-solid-surface bodies, stored for planets, moons, belt bodies and (where meaningful) stars (bug) | none | done, PR #4385 |
| GEN.1098 | Interstellar fields research read: the sources of design section 10 (LS12 and the rest), one GEN.725 row per number used | none | done, PR #4194 |
| GEN.1099 | physics/rogue_collision.py: rogue-rogue collision rate (cross-section with focusing), Maxwellian speed window, Q_R against Q*_RD, shock pressure, unbound mass and ejecta dispersion | none | done, PR #4197 |
| GEN.1100 | physics/field_dispersal.py: expansion radius from age, age from radius, freeze-out radius per size, shear term, retirement age, and a body position from (seed, rank, age) | none | done, PR #4200 |
| GEN.1101 | physics/field_rates.py: the four source rates (collision, supernova-unbound belts, and the others), the hosted rule for unbound belts, the labelled floor dial keeping both figures | none | done, PR #4205 |
| GEN.1102 | Field composition from two differentiated parents (iron, silicate, ice shards) replacing draw_field("collisional") | none | done, PR #4209 |
| GEN.1103 | AsteroidField columns, radius derived from age, the field density leaving 0, registry text changed | none | open |
| GEN.1104 | Supernova host rule: a young supernova remnant and its remnant core draw the unbound belt as a field, like NEBULA_HOST_RULES | none | open |
| GEN.1105 | Field members (individual bodies by the belt cut rule) with `<field ID>.<rank hex>` designations, status flags, and the sector and field pages | none | open |
| GEN.1106 | Field render: a 3D view for the field page (today _NO_VIEW), the sector map drawing the cloud or annulus plus the individual points at v x age, the render seed from the uid not the database id | none | open |
| GEN.1107 | Bring-to-epoch for fields: age advances, radius follows, retirement past 1 ly marks or deletes the row (Boss to pick; default retire) | none | open |
| GEN.1108 | Climate by latitude band for every body in seven bands (the equator, three above, three below, the poles inside the outermost band each side), the temperature categories of regional_temperature.py extended to seven: day and night length and the temperature range by day and by night (bug) | none | done, PR #4436 |
| GEN.1109 | Review rogue planet generation so physics drives it: where the counts, masses and properties rest on stand-ins or flat rates, and what replaces each (ejection from the architecture, scattering, the mass function, the survival of moons and air) | none | done, PR #4176 |
| GEN.1110 | Average wind speeds by day and by night in each of the seven latitude bands: the day-to-night thermal wind from the pressure difference, zonal jets and trade winds from the equator-to-pole temperature difference and the Rossby number, faded out below about 1 kPa | none | open |
| GEN.1111 | Rogue rows created from the architecture's ejections: mass, composition, moons kept, ejection speed, age and the parent system, with the flat rate becoming a test target | none | open |
| GEN.1112 | Replace the rogue-planet defaults one by one with their physical derivations (mass function from ejection, gas-giant threshold from the envelope fraction, age from the ejecting system, envelope and water share from formation, tidal flux from kept moons, radiogenic inputs from the parent star), each a small pure module | none | done, PR #4210 |
| GEN.1113 | Rogue ejection kick research (D1): read the full texts, extract ejection speed against perturber mass and orbit, direction statistics and moon survival, replace the kick stand-in, test against Sayed et al. | none | done, PR #4218 |
| GEN.1114 | physics/rogue_drift.py (D2): the closed-form drift of an ejected planet (the wrap rule and the velocity), checked against a direct orbit integration to 1 kpc and against hill_drift.py | none | done, PR #4178 |
| GEN.1115 | Smoothed star-density function for rogues (D3): ring average plus radial smoothing, shared by the sector draw and the low-zoom map | none | done, PR #4190 |
| GEN.1116 | Regenerate a rogue's donor on demand from its seed (D5): the history sentence and the parent record shown on the object page | none | dropped (replaced by the real donor link of GEN.1111) |
| GEN.1117 | Register the rogue drift and kick stand-ins in GEN.725 (D7): the rows of design section 8, each with tier and tweak level and a burn-down line where L3 or L4 | none | done, PR #4218 |
| GEN.1118 | Design study: one galaxy-wide object store generated in region work units (bands within the sphere), with sectors as a navigation table that flags whether each is charted and where its uncharted part lies | none | dropped (study done) |
| GEN.1119 | RP1: classify every architecture loss as ejected, collided or bound by the Safronov number and return an ejection record (mass, composition, kept moons, orbit at ejection, time), counting the small bodies the chain deletes today | none | done, PR #4389 |
| GEN.1120 | Rogue planet abundance option: a None, Few, Normal, Many dropdown in the galaxy and the sector generation options, with Normal the target and the physics-based spread tweaked to the low side of the observed range | none | open |
| GEN.1121 | Designation history: an object may change its name and ID, keeps a stored history of its previous designations, and its current ID is always completely unique across the whole galaxy | none | done (split into GEN.1138 to GEN.1141 after the design note (#4234) and Boss answers) |
| GEN.1122 | Read Pollack 1979 and Cuk and Burns 2004 for the gas-drag capture window strength and register the number (CM14) | none | done, PR #4192 |
| GEN.1123 | Audit every table that carries sector_id and its ON DELETE rule, and every path that deletes or regenerates a sector, and report which objects would be lost when a sector is deleted though they drifted in | none | done, PR #report sector-delete-audit-gen1123-2026-10-30 |
| GEN.1124 | Resident-mass subtraction: a pure function giving a sector's remaining mass budget after what it already holds (arrived rogues, pre-placed bright stars, scattered phenomena, extended objects whose centre lies there) | none | done, PR #4207 |
| GEN.1125 | Rogue top-up planner: a pure module that, from a sector's mass shortfall and the Rogues level cap, says how many rogues to place, their masses by the existing mass bins (96 percent terrestrial, 4 percent gas giant), and gives each a placement inside the sector, a realistic velocity and an apparent origin marked "apparent" | none | done, PR #4204 |
| GEN.1126 | Top-up asteroid field kind: lighter interstellar asteroid fields made when the Rogues cap leaves sector mass unmatched, smaller bodies, lower density, a larger volume, no full body history, origin a registered stand-in "unbound remnant" | none | done, PR #4212 |
| GEN.1127 | chart_ledger table (tile, run, epoch, time, topped_up_at) with the charted state of a sector derived from it, written by the sector fill units and backfilled from sector_stats.filled_at; the ON DELETE CASCADE removed from position-owned objects (rogue_planets first) so deleting or regenerating a sector removes only objects whose birth tile it is | none | open |
| GEN.1128 | Region run command: a sphere of a given radius from a start point becomes whole sectors, cut into units of a block of sectors by the galaxy's own (shell, band) arcs, nearest bands first, with `--grow-shells K`, all writing the chart ledger | none | open |
| GEN.1129 | Fingerprint split: a sector's digest (GEN.58) becomes "made here" and "arrived here" so rogues landing from other sectors do not change the sector's own digest | none | open |
| GEN.1130 | The +intelligent_life and +habitable_world request options become probability tweaks only: they raise the odds of such a world in the generated system and never force a placement, a ladder of adaptations or a civilization | none | open |
| GEN.1131 | Class descriptions for N, K and Q agree with the bodies their rules select | none | done, PR #4213 |
| GEN.1132 | Assess retiring Class Q: what depends on it, which bodies it covers, and whether eccentricity belongs in a separate orbit or climate flag, with a recommendation | none | done (report core4-class-q-assessment-2026-10-30.md, Q retired by GEN.1135) |
| GEN.1133 | Rocky rogue planets of 10 to 16 Earth masses come out up to 17,600 km, past class S's 11,468 km top: apply Core 4's recommendation so the size range and the class agree | none | done, PR #4385 |
| GEN.1134 | Rogue radiogenic inputs from the ejecting star's composition: the long-lived heat (Th, U, K) and 26Al of an ejected planet from the parent star's stored composition, replacing the sigma 0.22 dex offsets and ROGUE_ACTIVE_HEAT_FLUX_W_M2 | none | done, PR #4215 |
| GEN.1135 | Retire Class Q: remove the row from planet_classes.json, tuning.PLANET_CLASSES and HABITABLE_PLANET_CLASSES, and show the orbit swing as the stored flux_variation plus an "At a glance" flag when it is very high | none | done, PR #4385 |
| GEN.1136 | Move the Class V gravity line from 1.25 g to 2 g so Earth-sized and super-Earth temperate worlds can land in Class M, and register the line in GEN.725 | none | done, PR #4385 |
| GEN.1137 | cooling_track (GEN.948) jumps for a mass just under the lowest mass the table holds: 0.0699 solar masses at 5 Gy gives 1220 K, 0.07 gives 1627 K, so the temperature cap never binds below the limit mass (bug) | none | open |
| GEN.1138 | Designation history part b: the designations table (retired IDs and names), galaxy/designations.py (rules, retire, renumber) and the three rename writers and _rank_chain writing history | none | open |
| GEN.1139 | Designation history part c: a read of a retired ID redirects to the current object, a write against it answers 409, a deleted object answers 410 (api/ids.py and the page routes) | none | open |
| GEN.1140 | Designation history part d: Formerly lines on the object page, the former_designations API field and search by a past name | none | open |
| GEN.1141 | Designation history part e: bodies, belt members and captured bodies follow a parent ID change, and tombstones on delete | none | open |
| GEN.1142 | Class lines patch for the train 6 epoch: T hydrogen-envelope test ahead of N, N Venus-like sentence and hot, not-hydrogen-dominated rule, K without the under 273 K condition and ordered after P | none | done, PR #4385 |
| GEN.1143 | Class count script: temperate-ocean bodies per class, hydrogen-rich N bodies under 15 Earth masses, and the E, H and L reach (prints only) | none | done, PR #4243 |
| GEN.1144 | Post-train-6 class re-sweep on the new main, then decide the E and F merge, the T versus N boundary and the H and L reach | none | open |
| GEN.1145 | Ocean-fraction spread of temperate air-bearing bodies, printed before the train 6 epoch (decides whether the H line can fill) | none | done, PR #4247 |
| GEN.1146 | Upstream water inventory and land share: temperate bodies with air nearly always keep their water, so no desert world (Class H) exists; decide whether the generator should produce dry temperate planets (only if GEN.1145 shows H stays empty) | none | done, PR #report-ocean-fraction-physics-gen1146 |
| GEN.1147 | Ocean-fraction law: replace the sqrt cover law with a two-level hypsometry plus the mantle water cap by mass (2 to 8 oceans), in one epoch | none | open |
| GEN.1148 | Asteroid belt share fell from 55 percent of systems to 3 percent and binary share rose from 24 percent to 31.5 percent between the census before trains 3 and 4 and current main (Bug 2 census); find out whether each is the intended change or a bug (bug) | none | done, report gen1148-belt-binary-shift-2026-10-31.md |
| GEN.1149 | Lagrange trojans T1: pure physics module physics/trojans.py: Gascheau-Routh stability test (host, secondary and trojan masses), libration period, tadpole and horseshoe limit, e and i screen, swarm mass from the scatter trojan share with leading and trailing split, plausibility harness against Jupiter, Neptune, Saturn moons and the exoplanet upper bound (WR10), GEN.725 registry rows, design note | none | done, PR #4273 |
| GEN.1150 | Lagrange trojans T2: the generator draws a trojan swarm at the L4 and L5 of each giant planet or brown dwarf from the scatter trojan share; new columns co_orbital_host_id, lagrange_point, libration amplitude and phase (db/column_units.py); orbit updater libration overlay; the object and sector pages and the designation rule "<host> L4 trojan region" | none | open |
| GEN.1151 | Lagrange trojans T3: whole planets and moons (Class M and Class J allowed by the mass rule) at L4 and L5 from the disc formation model only; no forced rate, the observed 5 percent upper bound (WR10) is a plausibility test only; folds in the weird-worlds W7 row | none | open |
| GEN.1152 | Mass literals that differ from constants.py (changing them changes output): moon_composition JUPITER_MASS_KG 1.89813e27, ring_origin star_mass_kg 1.989e30, giant_interior 317.8 and EARTH_MASS_JUPITER 317.83, stellar_abundances JUPITER_MASS_SOL, tuning 317.8; use the constants.py values with an epoch | none | open |
| GEN.1153 | Day and night state for stars and rogue planets (follow-up to GEN.1097, which covers planets and moons): a star stores the constant "no night"; a rogue planet gets the same state, night hours and solar days as an in-system planet from its own spin, with no star the answer is no solar day | none | open |
| GEN.1154 | Halo mass of galaxies above 3e10 solar masses of stars is drawn too heavy by GEN.368 (0.17 dex at 6e10, 0.27 at 1e11, 0.34 at 2e11 against a Tinker 2008 halo function): replace median_halo_mass plus the Gaussian with a precomputed Tinker-based conditional drawn by inverse CDF | none | done, PR #4385 |
| GEN.1155 | A REBOUND N-body sample of 2+2 quadruple systems to replace the 2+2 stability margin GM14 that has no N-body sample behind it (the 3+1 margin GM13 was tuned on the 331 systems of GEN.1029) | none | done, PR #4343 |
| GEN.1156 | Raise the belt share after the GEN.728 disc architecture (3 to 12 percent of systems) toward the documented 4 to 14 percent with a physics-driven change | none | open |
| GEN.1157 | Earth-type moons around gas giants: how often a giant has a moon the size of a planet, and how big, from the physics of moon formation and capture | none | open |
| GEN.1158 | Roundness engine: a pure module physics/roundness.py with yield strength from ice fraction, the round radius from density and strength by the Lineweaver and Norman law, the roundness ratio, the 10 percent sphericity rule and the binding energy of a sphere | none | done, PR #4338 |
| GEN.1159 | Store ice fraction, bulk density, roundness ratio, yield strength and binding energy on planets, moons, asteroids and belt or field members (columns, migration, API version) | none | open |
| GEN.1160 | Class D and C rules read the roundness ratio and the 0.1 Earth mass line, and the Asteroid object type (no class letter, moon-capable) for bodies under roundness ratio 1 | none | open |
| GEN.1161 | Asteroid field mass rule and a member ledger like the belts' (individual bodies above 1e19 kg stored, the rest drawn on demand) | none | open |
| GEN.1162 | Store bulk density and binding energy (and roundness ratio, yield strength and ice fraction where they apply) on stars, brown dwarfs, giant planets, white dwarfs, neutron stars and black holes; a black hole's binding energy is stored as -1, meaning infinite | none | open |
| GEN.1163 | Stars hotter than about 6300 K all take the hot-star spin branch (a median 180 km/s whatever the mass), so F stars rotate in under a day where the real Kraft break rises gradually (about 20 to 70 km/s) (bug) | none | done, PR #4385 |
| GEN.1164 | About 3 percent of moons (70 of 2,273) sit in a 3:2 or 2:1 spin-orbit state, all captured moons with eccentricity 0.2 to 0.5, where real moons have none apart from chaotic Hyperion: find out whether the tide damps the eccentricity before the spin settles | none | done, answered by Tides 1 (GEN.1172) |
| GEN.1165 | Giant planet obliquity never exceeds 60 degrees (312 giants) and locked planets have exactly 0 degrees tilt: draw a wider, physically derived obliquity and a small non-zero tilt for locked bodies | none | done, PR #4385 |
| GEN.1166 | A forced F5V star gets 1.24 to 1.69 solar masses at a fixed 6500 K: the forced-type temperature override overrules the mass, as GEN.365 found for hot types (bug) | none | dropped (withdrawn by Aux 7: a forced F5V is an exact model star; the 6500 K was rounding and the mass spread is main-sequence evolution) |
| GEN.1167 | The stellar-to-halo relation (Moster 2013) has 0.3 to 0.4 dex more stars than Behroozi 2019 in halos of 1e10.75 to 1e11.5 and 0.06 to 0.26 dex fewer in halos of 1e12.5 to 1e14.5: compare the published fits and pick the best-supported relation | none | open |
| GEN.1168 | No white dwarf is lighter than 0.56 solar masses in the generated population while field samples go down to 0.45: check whether the low-mass white dwarf channel is missing | none | done, report gen1168 (GEN.996 targets the shares) |
| GEN.1169 | Pre-main-sequence stars are 2.5 percent of living stars at the Sun's position (late M dwarfs, median age 0.45 Gy): check the share against the contraction times of the tracks and a published fraction | none | done, report gen1169 (GEN.1173, GEN.1174) |
| GEN.1170 | Spin of A and B stars follows mass: the fast-rotator mode rises from 160 km/s at 1.6 to 2 solar masses to 238 km/s at 3.15 to 3.85 (about 46 km/s per solar mass), where the generator keeps one 180 km/s median for every hot star from 1.8 solar masses up | none | done, PR #4385 |
| GEN.1171 | Ring lifetime by ring mass and planet instead of one fixed decay, because fallen-moon rings now persist on 13 percent of rocky planets | GEN.1057 | open |
| GEN.1172 | Spin-orbit resonance capture is a draw, not a certainty: 3:2 and 2:1 capture inside the window with a probability inside the published 0.07 to 0.55 bracket | none | done, PR #4385 |
| GEN.1173 | Time to the main sequence for 0.1 to 0.3 solar-mass stars has an unphysical step between 0.105 and 0.11 solar masses (718 to 1358 Myr) from the isochrone arrival grid | none | done, PR #4436 |
| GEN.1174 | Decide the pre-main-sequence time definition: the generator's time to the main sequence is 3 to 4 times BHAC15 for 0.1 to 0.3 solar masses and 1.5 to 2 times from 0.4 up | none | open |
| GEN.1175 | Porosity weakening of body strength from a sourced law instead of the placeholder exponential | none | open |
| GEN.1176 | Rubble-pile branch for the rounding engine: friction plus a little cohesion and a shape-set temperature instead of a size threshold | none | done, PR #4347 |
| GEN.1177 | Strength rows for the ice phases, methane clathrate, solid CO2 and ammonia ices | none | done, PR #4352 |
| GEN.1178 | Strength of a body sampled from its range so the round boundary is fuzzy, as in the real Solar System | GEN.1159 | open |
| GEN.1179 | Optional: Hoek-Brown geological strength index as a body parameter for rock | none | open |
| GEN.1180 | Follow-up strength research: peridotite and olivine, rock fracture toughness, solid N2 and CO ices, ices VII and X, methane hydrate constants at low temperature | none | done, report strength-followup-gen1180 (olivine flow constants read; no data exists for the rest and none is needed for rounding) |
| GEN.1181 | Set the molecular gas velocity dispersion to 7.0 km/s (range 5 to 9) from four read sources and record it in GEN.725 and the structure.py comment | none | dropped (Replaced by the item that keeps 6.5 km/s and records its sources: the young-star scale height sources (GEN.1182) put 66 pc inside the sourced Milky Way band and 73 pc (7.0 km/s) just outside it.) |
| GEN.1182 | Sourced scale height of young stars: the 'observed about 50 pc' figure has no source in anything read and no CO velocity dispersion reaches it | none | done, report gen1182 (sourced Milky Way band 35 to 70 pc; kept in GEN.1188) |
| GEN.1183 | Ice-type strength ratios (IT1 to IT3) become a ratio-by-temperature table from Durham 1997 Table 2 instead of constants | none | done, PR #4383 |
| GEN.1184 | Rock strength follows olivine creep in the roundness temperature rules, and the hot rubble strength RN9 is re-tiered to olivine creep at about 1300 K | none | done, PR #4359 |
| GEN.1185 | Interior activation energy range 260 to 390 kJ/mol instead of 200 to 300, with the reference viscosity retuned | none | open |
| GEN.1186 | Small evidence updates to registered stand-ins from the strength papers review | none | open |
| GEN.1187 | Strength constants still without a source: ammonia dihydrate 1993 constants, the melt factor, Andrade creep, ice grain-boundary-sliding laws and the ice phase boundary table | none | done, report gen1187 (GEN.1190) |
| GEN.1188 | Record the sources of the 6.5 km/s molecular gas dispersion and pin the Milky Way young-star scale height test to the sourced band of 35 to 70 pc | none | done, PR #4364 |
| GEN.1189 | Greenhouse refit method: the refit table is blind to the methane share of a hydrogen air and a dense carbon dioxide air is 2 to 8 percent warm, so find a better fit than one L-BFGS network | GEN.1032 | open |
| GEN.1190 | Replace unsourced rheology constants with the sourced ones: partial-melt factor, Andrade exponent and ratio, ice grain-boundary sliding and the ice melting curves | none | open |
| GEN.1191 | Satellite supernova remnants scale with the satellite's own star formation instead of the Milky Way's per-system rate | none | done, PR #4436 |
| GEN.1192 | Satellite gas: a gas disc longer than the stars and a longer depletion time at low density | none | done, PR #4436 |
| GEN.1193 | Stripping rule on a satellite's distance to its host so gas follows Grcevich and Putman | none | done, PR #4436 |
| GEN.1194 | Runaway greenhouse: the radiative-convective model has no plateau in outgoing flux, so water window opacity, continuum temperature dependence and wings must be reworked until it plateaus | none | open |
| GEN.1195 | Cluster maximum mass of the central ring from the gas, and cloud-shock disruption of its clusters | none | done, PR #4436 |
| GEN.1196 | Disc cluster maximum mass per ring, replacing the flat 1e6 cap | none | done, PR #4436 |
| GEN.1197 | Faint gas-bearing dwarfs far from the host: the star formation threshold stops anything under 1e7 solar masses, so Leo T-like dwarfs are never made | none | done, research only, no change (gen1197-faint-gas-dwarfs-2026-10-31.md) |
| GEN.1198 | The sector profile's remnant factor in a Magellanic Cloud cell uses the gas factor, so scatter and fill disagree with the satellite remnant count | none | done, PR #4436 |
| GEN.1199 | A small moon of a giant carries a hydrogen air and is called class T: can a 0.1 Earth-mass moon keep a hydrogen air at all? | none | done, PR #4436 |
| GEN.1200 | The galactic-centre ring inflow of the generator is 0.044 solar masses per year against the measured 0.07 | none | open |
| GEN.1201 | A per-body shape-setting temperature: the roundness engine uses the temperature a body took its shape at, but no body stores one, so strength is the cold strength today | none | open |
| GEN.1202 | The generator's disc gas surface density at 8 kpc is 4.3 solar masses per pc2 against 13 measured, which holds the disc cluster maximum mass too low | none | open |
| GEN.1203 | Polluted white dwarf census: 25 to 50 percent of white dwarf systems polluted once belt birth sees the giants | none | open |
| GEN.1204 | Dark members of an orbit tree (brown dwarf, neutron star, black hole in dark_members) store density, binding energy and the shape values the stars, rogues and comets store | none | open |
| GEN.1205 | Best-fit classes, step 1: every class gets its own full positive conditions (temperature, pressure, gravity, water, air chemistry) as new data, with shoulder widths registered in GEN.725 | none | done, PR #4469 |
| GEN.1206 | Best-fit classes, step 2: a read-only shadow scorer prints the confusion table (stored letter against best-fit letter), the share at 100 percent and the low-match list | none | done, PR #4473 |
| GEN.1207 | Best-fit classes, step 3: the classifier replaces the first-match tree, stores class_match and the runner-up letter and match, and the class pages show the percent match | none | open |
| GEN.1208 | Best-fit classes, wording: drop the age and stage words from the class tier labels and the "Checked after" and "leftover" wording from the class pages | none | open |
| GEN.1209 | Several life types coexist on one planet: a per-body list of living classes, each with its stage, animal grade, kinds (microbes, flora, fauna) and prevalence | none | open |
| GEN.1210 | Life readings: abundance, trend and security for each kind of life (microbes, flora, fauna, sentient), as one pure function over stored values with no new column | none | done, PR #4481 |
| GEN.1211 | A new class for hot thin-air rocks (about 1.6 percent of bodies: 1.4 to 2 g, 340 to 1100 K, 1 to 10 kPa, no ocean) is Class Z of the class code system (Boss confirmed) | none | open |
| GEN.1212 | Rewrite the stored system flavour lines (tuning.SYSTEM_FLAVOR, about 40 story lines) to the plain-words style guide, and stop storing text the pages can build | none | open |
| GEN.1213 | Atmospheric composition accounts for the whole mix of living classes on a planet: each class changes the air in its own way, summed in timeline order | none | open |
| GEN.1214 | Carbon dioxide and hydrocarbon (methane-ethane) seas are not being created: find out why on a census of the surface liquid and fix it (bug) | none | open |
| GEN.1215 | Store the acid and ammonia fractions of the surface liquid, and let habitability read them (an acid or ammonia sea can make a world uninhabitable) | none | open |
| GEN.1216 | Research: the composition of surface liquids from physics and from life: salts and brines, acids, ammonia-water, liquid CO2, methane-ethane-nitrogen mixes, sulfur and what life and chemistry do to each | none | done, research report written (reports/gen1216-surface-liquid-composition-2026-10-31.md) |
| GEN.1217 | Design: the surface liquid composition model: what is stored (fractions of water, salts, acid, ammonia, CO2, hydrocarbons), how it is derived, and how life changes it | none | open |
| GEN.1218 | Build the surface liquid composition model: the composition columns, the derivation, the effect of life on it and of it on life, pages and API | none | open |
| GEN.1219 | Body types: planet_type (terrestrial, transitional, gas) from the envelope share of the mass, with a shadow report first | none | open |
| GEN.1220 | Class code system: a letter plus a digit (class_variant column, 0 not written), the digit meanings and the letters of the oxygen-water, anoxic-water, CO2, acid, methane-ethane, ammonia-water and hot thin-air families | none | open |
| GEN.1221 | Class-by-domain matrix as data: six domains plus regional, day and night temperature and wet bulb, stored relative humidity and a groundwater flag, and a phase gate per solvent | none | open |
| GEN.1222 | Tiny airless moons hold a water surface ocean at 20 to 33 K: 60 of 2,071 bodies in one sample (39 class D, 20 class C, 1 class T), all moons (bug) | none | open |
| GEN.1223 | Research: the gas and Transitional class lists from sources (Sudarsky-style temperature classes, sub-Neptune, ice giant, sub-Saturn, water world) | none | done, research report written (reports/gen1223-giant-and-transitional-classes-2026-10-31.md) |
| GEN.1224 | Research: envelope share against mass and metallicity (Thorngren 2016 mass-metallicity relation) to replace the envelope-share guesses of the giant and transitional class rows | none | done, research report written (reports/gen1224-giant-heavy-element-and-envelope-share-2026-10-31.md) |
| GEN.1225 | Research: tune the ammonia life gate and the nitrogen ice amount from physics and chemistry (instead of guesses) | none | done, research report written (reports/gen1225-ammonia-gate-and-nitrogen-ice-2026-10-31.md) |
| GEN.1226 | Store the envelope share for EVERY body, so hydrogen-air dwarfs (about 400 class T bodies) are typed and classed correctly | none | open |
| GEN.1227 | Store the coma class and the dominant volatiles of an interstellar comet, so the map can colour its coma | none | open |
| GEN.1228 | Hematite (the oxidised ground type) reflects the same in all three bands, so a rusty world renders dark grey instead of red; give it a sourced spectrum or finer bands (bug) | none | open |
| GEN.1229 | Methane haze depth per unit of mass flux: the fall term in aerosols.haze_fall (eddy speed replacing fall speed) and the aggregate cross section give Titan its depth only through an oversized mass flux (bug) | none | open |
| GEN.1230 | The disc star formation history should decline to the present for a Milky Way mass galaxy: the model gives 0.49 to 0.55 of the ten Gy mean, measured is 0.20 to 0.3 (bug) | none | open |
| GEN.1231 | White dwarf pollution flag reads the belt's scattered rate (is_polluted, pollution_rate_g_s) and the census band 25 to 50 percent | none | open |
| GEN.1232 | Belt history shown on the belt page: history_text and the giant event / companion limit (giant_event) | none | open |
| GEN.1233 | Dust rings fed by small moons (Jupiter, Uranus, Neptune): a steady-state supply so giants with faint rings are not left ringless by the ring lifetime law | none | open |
| GEN.1234 | Four stored ozone columns are below the least the oxygen allows because the oxygen changed after the column was drawn (bug) | none | done, PR #4699 |
| GEN.1235 | The sweep rule does not know the 2000 K surface temperature clamp: Icilat Aobaris I has 2000.0 where albedo and greenhouse give 2165 K (bug) | none | done, PR #4699 |
| MAP.1 | Galaxy Map follow-ups (edge cases) | 11 (2026-09-30 16:44Z); 12 (2026-09-30 16:49Z to 18:09Z); 19 (2026-09-30 18:14Z to 2026-10-01 04:16Z) | done in 7.42.1, PR #168 |
| MAP.2 | Drill-down navigation (new parent) | none | done (all subitems shipped), PR #234 |
| MAP.3 | A bigger Galaxy Map with controls underneath | 63 (2026-10-01 01:44Z to 05:05Z) | done in 7.55.0, PR #178 |
| MAP.4 | System Map names never overlap | 62 (2026-09-30 20:01Z to 20:27Z); 49 (2026-09-30 20:08Z to 2026-10-01 04:37Z) | done in 7.21.1, PR #135 (see note 3) |
| MAP.5 | Rework the Galaxy Map | 9 (2026-09-24 01:32Z to 02:18Z); 6 (2026-09-24 01:57Z to 02:02Z); 4 (2026-09-24 02:25Z to 05:38Z); 3 (2026-09-24 02:53Z to 2026-09-30 17:58Z) | replaced on 2026-09-30 by MAP.31 to MAP.42 and MAP.1; bugs MAP.37 and MAP.43 fixed, PR #201 |
| MAP.6 | Unfilled-sector skeleton draws as a sphere | 5 (2026-09-24 01:32Z to 02:18Z); 2 (2026-09-24 01:57Z to 02:02Z) | done in 5.47.0, PR #72 (test added in 5.47.1, PR #73) |
| MAP.7 | Remove the large sphere marker for a star in a sector | 2 (2026-09-24 01:32Z to 02:18Z) | done in 5.47.1, PR #73 |
| MAP.8 | System Map orbital paths back, not over bodies | 3 (2026-09-24 01:32Z to 02:18Z) | done in 5.47.1, PR #73 |
| MAP.9 | Star glow renders as an opaque shell | 4 (2026-09-24 01:32Z to 02:18Z) | done in 5.47.1, PR #73 |
| MAP.10 | Draw the Measure distance path around obstacles | 7 (2026-09-24 01:32Z to 02:18Z); 4 (2026-09-24 01:57Z to 02:02Z); 2 (2026-09-24 02:25Z to 2026-09-30 18:09Z); 9 (2026-09-30 18:14Z to 23:48Z) | done in 7.22.0, PR #135 |
| MAP.11 | Every kind of phenomenon on the Sector Map, clickable | 8 (2026-09-24 01:32Z to 02:18Z); 5 (2026-09-24 01:57Z to 02:02Z); 3 (2026-09-24 02:25Z to 05:38Z) | done in 5.51.0, PR #81; bugs MAP.45 and MAP.46 fixed, PR #200 |
| MAP.12 | Arc-segment wireframe on the Sector Map | 60 (2026-09-30 20:01Z to 20:27Z); 47 (2026-09-30 20:08Z to 2026-10-01 00:34Z) | done in 7.29.1, PR #143 |
| MAP.13 | Sector Map star dots sized to giants and white dwarfs | 55 (2026-10-01 00:07Z to 00:34Z) | done in 7.29.2, PR #143 |
| MAP.14 | Bright stars on the Galaxy Map (new done item; never had a number) | none | done in 7.42.0, PR #160; bugs MAP.47, MAP.48 (PR #201) and MAP.51 (PR #214, #218) fixed |
| MAP.15 | Stars and glowing phenomena as points of light on the Sector Map | 108 (2026-10-01 05:44Z to 05:54Z) | done, PR #234 |
| MAP.16 | Drill-down stages | 66 (2026-10-01 02:24Z to 02:26Z); 72 (2026-10-01 02:27Z to 04:33Z) | done in 7.44.0, PR #171 (bug MAP.17 fixed, PR #208) |
| MAP.17 | No free camera: drill down from a top-down view by wedge, slice and block (bug) | 100 (2026-10-01 05:22Z to 05:29Z) | done, PR #208 |
| MAP.18 | The block under the pointer is too hard to see from above (bug) | 111 (2026-10-01 05:50Z to 05:56Z) | done, PR #208 |
| MAP.19 | Big wedge picks in the drill-down (bug) | none | done, PR #208 |
| MAP.20 | Generate from the sector level | 67 (2026-10-01 02:24Z to 02:26Z); 73 (2026-10-01 02:27Z to 05:29Z) | done: map buttons and radius dialog in 7.53.0 (PR #177), block and layer generate after 7.58.2 (PR #182, #183) |
| MAP.21 | Sector Map pick mode and Nav links | 68 (2026-10-01 02:24Z to 02:26Z); 74 (2026-10-01 02:27Z to 05:05Z) | done in 7.58.0, PR #178 |
| MAP.22 | NAV page picks on the map | 69 (2026-10-01 02:24Z to 02:26Z); 75 (2026-10-01 02:27Z to 05:29Z) | done, Bookmarks select on NAV, PR #234 |
| MAP.23 | Bookmarks | 70 (2026-10-01 02:24Z to 02:26Z); 76 (2026-10-01 02:27Z to 05:29Z) | done, per-browser bookmarks (static/bookmarks.js), PR #234 |
| MAP.24 | Address bar | 71 (2026-10-01 02:24Z to 02:26Z); 77 (2026-10-01 02:27Z to 04:41Z) | done in 7.50.0, PR #172 |
| MAP.25 | "Show on Galaxy Map" links | 72 (2026-10-01 02:24Z to 02:26Z); 78 (2026-10-01 02:27Z to 05:29Z) | done after 7.58.2, PR #188 (bug MAP.26 fixed, PR #208) |
| MAP.26 | "Show on Galaxy Map" opens at the sector; map Back and Forward (bug) | 101 (2026-10-01 05:24Z to 05:29Z) | done, PR #208 |
| MAP.27 | NAV course on the Galaxy Map | 73 (2026-10-01 02:24Z to 02:26Z); 79 (2026-10-01 02:27Z to 04:51Z) | done in 7.52.0, PR #176 |
| MAP.28 | Nested ladder geometry (drill-down section 3) | 64 (2026-10-01 02:24Z to 02:26Z); 70 (2026-10-01 02:27Z to 03:54Z) | done in 7.41.2, PR #160 |
| MAP.29 | Stage contents API (drill-down section 7) | 65 (2026-10-01 02:24Z to 02:26Z); 71 (2026-10-01 02:27Z to 03:54Z) | done in 7.41.3, PR #160 |
| MAP.30 | Slab list to the left of the map, and a 3:4 map | 109 (2026-10-01 05:50Z to 05:56Z) | done (slab slider to the right of a 4:3 map), PR #234 |
| MAP.31 | Spiral arms stand out in the density shading | 3 (2026-09-30 16:44Z to 18:09Z); 10 (2026-09-30 18:14Z to 22:03Z) | done in 7.9.0, PR #120 |
| MAP.32 | Scale readout in sectors, pc and ly | 4 (2026-09-30 16:44Z to 18:09Z); 11 (2026-09-30 18:14Z to 22:03Z) | done in 7.9.0, PR #120 |
| MAP.33 | Hybrid master-wedge slot rule (schema v35) | 5 (2026-09-30 16:44Z to 18:09Z); 12 (2026-09-30 18:14Z to 22:42Z) | done in 7.13.0, PR #129 |
| MAP.34 | Mega-blocks sized from the pixel scale | 6 (2026-09-30 16:44Z to 18:09Z); 13 (2026-09-30 18:14Z to 23:14Z) | done in 7.17.0, PR #134 |
| MAP.35 | Continuous blocks and a Slice control | 7 (2026-09-30 16:44Z to 18:09Z); 14 (2026-09-30 18:14Z to 22:03Z) | done in 7.9.0, PR #120 |
| MAP.36 | One solid of blocks for filled and unfilled sectors | 8 (2026-09-30 16:49Z to 18:09Z); 15 (2026-09-30 18:14Z to 2026-10-01 00:35Z) | done in 7.25.0, PR #142 |
| MAP.37 | Generated systems are hard to find on the map (bug) | 99 (2026-10-01 05:15Z to 05:29Z) | done, PR #201 |
| MAP.38 | Block info on click | 8 (2026-09-30 16:44Z); 9 (2026-09-30 16:49Z to 18:09Z); 16 (2026-09-30 18:14Z to 2026-10-01 00:35Z) | done in 7.25.0, PR #142 |
| MAP.39 | Smooth zooming | 9 (2026-09-30 16:44Z); 10 (2026-09-30 16:49Z to 18:09Z); 17 (2026-09-30 18:14Z to 2026-10-01 01:44Z) | done in 7.32.0, PR #147 |
| MAP.40 | Keep three.js; record why | 10 (2026-09-30 16:44Z); 11 (2026-09-30 16:49Z to 18:09Z); 18 (2026-09-30 18:14Z to 2026-10-01 02:31Z) | done, PR #153 (docs only, shipped with 7.36.0) |
| MAP.41 | Remove the server's leftover density sampling | 12 (2026-09-30 16:44Z); 13 (2026-09-30 16:49Z to 18:09Z); 20 (2026-09-30 18:14Z to 22:03Z) | done in 7.9.0, PR #120 |
| MAP.42 | Wedge lines from the center | 21 (2026-09-30 18:14Z to 22:03Z) | done in 7.9.0, PR #120 |
| MAP.43 | Wedge lines run past the galaxy's edge (bug) | 96 (2026-10-01 05:15Z to 05:29Z) | done, PR #201 |
| MAP.44 | Wedge lines and ring circles run far past a zoomed-in block (bug) | 110 (2026-10-01 05:50Z to 05:56Z) | done, PR #208 |
| MAP.45 | Rogue planets (and maybe other objects) drawn outside the sector's wireframe (bug) | 104 (2026-10-01 05:40Z to 05:54Z) | done, PR #200 |
| MAP.46 | Rogue planets are hard to find on the Sector Map (bug) | 107 (2026-10-01 05:44Z to 05:54Z) | done, PR #200 |
| MAP.47 | Bright stars vanish when zoomed out (bug) | 97 (2026-10-01 05:15Z to 05:29Z) | done, PR #201 |
| MAP.48 | Stars take a while to appear after a zoom (bug) | 98 (2026-10-01 05:15Z to 05:29Z) | done, PR #201 |
| MAP.49 | The System Map shows planet orbits inside an asteroid belt (bug) | none | done, PR #204 |
| MAP.50 | Names run off the edge of the map (bug) | none | done, PR #204 |
| MAP.51 | No stars drawn in filled sectors past certain zoom levels (bug, under MAP.14) | none | done, PR #214 and #218 |
| MAP.52 | Galaxy Map highlights the wrong area; pick a 40-degree wedge around the cursor (bug) | none | done, PR #369 |
| MAP.53 | Rotate a zoomed-in wedge, and zoom it to fit the window (bug) | none | done, PR #410 |
| MAP.54 | Slab leader lines instead of the slab slider (bug) | none | done, PR #410 |
| MAP.55 | Galaxy Map buttons: a menu, with only back, forward, up, reset and bookmark showing | none | done, PR #369 |
| MAP.56 | Drop the 3x3 block pick: select a slab, zoom in, select a segment (bug) | none | done, PR #408 |
| MAP.57 | The System Map writes NaN or infinite positions into its SVG (bug) | none | done, PR #405 |
| MAP.58 | Galaxy Map zoom limits: a short manual range on the galaxy wedge, locked below it | none | superseded by MAP.125 (Boss, 2026-10-09) |
| MAP.59 | Make it plain that a zoomed-in slab is a slab, not a wedge | none | dropped (Boss, 2026-10-16, 3:00 PM EDT; superseded: a picked region is a shell, band or wedge outline now) |
| MAP.60 | Galaxy Map scale readout: one scale line | none | done, PR #369 |
| MAP.61 | One map engine and control set for the Galaxy Map and the Sector Map | none | done, PR #629 |
| MAP.62 | A full 3D star system view with a free camera | none | done, PR #673 |
| MAP.63 | Shared map helpers in one module | none | done, PR #351 |
| MAP.64 | One camera and input controller | none | done, PR #351 |
| MAP.65 | One picking, hover and info-panel layer | none | done, PR #521 |
| MAP.66 | The sector as the drill-down's last stage, on the same page | none | done, PR #548 |
| MAP.67 | One URL and history scheme for every level | none | done, PR #673 |
| MAP.68 | Remove the old Sector Map code | none | done, PR #573 |
| MAP.69 | A system scene endpoint with 3D orbits | none | done, PR #668 |
| MAP.70 | Positions at any time | none | done, PR #670 |
| MAP.71 | Scale modes that keep everything visible | none | done, PR #673 |
| MAP.72 | Rendering at system scale | none | done, PR #673 |
| MAP.73 | Free camera on the shared engine | none | done, PR #673 |
| MAP.74 | The 3D view on the system page, the flat diagram kept | none | done, PR #673 |
| MAP.75 | The mini map as a second engine view | none | open |
| MAP.76 | Leader-line layout | none | done, PR #410 |
| MAP.77 | Galaxy Map draws block divisions inside a picked slab before zooming to it (bug) | none | done, PR #410 |
| MAP.78 | Zooming into a wedge must show the whole wedge at every drill-down level (bug) | none | done, PR #410 |
| MAP.79 | Rogue planets clog the Sector Map: dim them, and a show/hide button per kind of object (bug) | none | done, PR #577 |
| MAP.80 | Sector-level zoom on the Galaxy Map should show almost every star in the sector (bug) | none | done, PR #429 |
| MAP.81 | Ctrl+1 to Ctrl+9 bookmark keys clash with the browser's tab switching (bug) | none | done, PR #351 |
| MAP.82 | Unmarked rogue planets barely visible (bug) | none | done, PR #351 |
| MAP.83 | The "Mark rogue planets" button shows when it is on (bug) | none | done, PR #351 |
| MAP.84 | Marked rogue planets grow and become clickable; unmarked ones stay small (bug) | none | done, PR #351 |
| MAP.85 | The galaxy pick is an arc, on a 3D galaxy with no sector lines | none | done, PR #369 |
| MAP.86 | Sector and block colors from what is in them: filled sectors translucent (bug) | none | done, PR #432 |
| MAP.87 | Stars on the Sector Map and Galaxy Map need to be brighter, most of all the dim ones (bug) | none | done, PR #351 |
| MAP.88 | Parts of a star system run off the edge of the System Map (bug) | none | done, PR #405 |
| MAP.89 | System Map: space orbits with a fitted scale and a minimum ring gap instead of plain log | none | closed as done (Boss, 2026-10-09) |
| MAP.90 | The tile-level helper crashes on a subnormal view radius (bug) | none | done, PR #365 |
| MAP.91 | Hovering the map while picking a slab highlights the whole slab, not one cube | none | done, PR #395 |
| MAP.92 | The System Map's side panel leaves out a planet's or moon's radius and mass (bug) | none | done, PR #405 |
| MAP.93 | The Galaxy Map breadcrumb wraps onto several lines instead of collapsing its middle steps into a "…" menu (bug) | none | done, PR #399 |
| MAP.94 | On a phone the breadcrumb should give way to a round menu button between Back and Forward (bug) | none | done, PR #399 |
| MAP.95 | A "Forward to current" button next to the map's Back and Forward | none | done, PR #562 |
| MAP.96 | The Galaxy Map can't be turned freely: the tilt stops at straight down and at 80 degrees (bug) | none | done, PR #413 (open: keep the rotation in the URL and bookmarks? Not stored today) |
| MAP.97 | The Galaxy Map camera should go top-down for the galaxy and a slab, isometric for a block, at every zoom step (bug) | none | done, PR #413 (defaults: a manual turn does not carry to the next step; Reset view returns to the step's preset) |
| MAP.98 | Slab button lines should end at the nearest edge of their slab (bug) | none | done, PR #422 |
| MAP.99 | Slab buttons that don't fit the window split across both sides of the map, shrink, or give way to map picking (bug) | none | done, PR #422 |
| MAP.100 | Slab button labels on one line: "#N" and how much is charted (bug) | none | done, PR #422 |
| MAP.101 | Stars on the Galaxy Map take the click, so a dense sector can't be picked (bug) | none | done, PR #431 |
| MAP.102 | Galaxy Map streaming with a BVH and 3D tiles, and camera-relative rendering | none | done, PR #574 |
| MAP.103 | Nebulae don't show on the Galaxy Map or any other map (bug) | none | done, PR #592 |
| MAP.104 | Nebula shading is missing on unfilled sectors (bug) | none | done, PR #592 |
| MAP.105 | The nebula view should show its whole shape with the dimmed galaxy around it (bug) | none | done, PR #590 |
| MAP.106 | The breadcrumb trail falls out of sync with the map (bug) | none | done, PR #586 |
| MAP.107 | Selecting an empty slab near the core says "There is no layer x here" (bug) | none | done, PR #537 |
| MAP.108 | Empty slabs and wedges near the core can't be selected, and the side buttons block clicks (bug) | none | closed, not reproducible; covered by tests (#537, #541, #521) |
| MAP.109 | Zooming in and out loads slowly (bug) | none | done, PR #603 |
| MAP.110 | Slab button lines come out of numerical order (bug) | none | done, PR #520 |
| MAP.111 | "Generated only" should be "Charted only" and dim the stars too (bug) | none | done, PR #531 |
| MAP.112 | Nothing can be selected while "Generated only" is on (bug) | none | done, PR #541 |
| MAP.113 | A nebula covering the whole sector can't be unselected, and nebulae need a show/hide toggle (bug) | none | done, PR #584 |
| MAP.114 | The Sector Map's Reset view leaves the picture unchanged (bug) | none | done, PR #442 |
| MAP.115 | Comets, rogue planets and asteroid fields show above the sector level (bug) | none | merged into MAP.116; done, PR #601 |
| MAP.116 | Dense generated sectors crowd the Galaxy Map when zoomed out (bug) | none | done, PR #601 |
| MAP.127 | See into and step to the neighbouring blocks and slabs on the Galaxy Map | none | merged into MAP.121 (2026-10-07) |
| MAP.128 | Decide what the Sector Map tints by star age, density and luminosity, and test it | none | closed, decided: sector space is not tinted (Boss 2026-10-08 13:35Z); the Galaxy Map half was done in #556 |
| MAP.129 | Blocks colored from their sectors' statistics, weighted by star count | none | done, PR #556 |
| MAP.130 | Retire the average-star-color rule in the docs and tests | none | done, PR #556 |
| MAP.131 | A Color by switch on the Galaxy Map, with a legend | none | done, PR #692 |
| MAP.132 | Overlay markers for black holes, nebulae and habitable worlds | none | open |
| MAP.133 | The Galaxy Map shows no bright stars in the bulge, and only layers -121 to 121 (bug) | none | done, PR #722 |
| MAP.134 | Build the Galaxy Map's opening view ahead of time on every update | none | done, PR #764 |
| MAP.117 | Surface pressure missing from the planet and moon side panel (bug) | none | done, PR #457 |
| MAP.118 | The Galaxy Map shows an unfilled sector above the galaxy that can't be filled (bug) | none | done, PR #457 |
| MAP.119 | Expected star density editable by admins on the Galaxy Map | none | open |
| MAP.120 | Bright-star backfill from the Galaxy Map's block, slab and wedge menus | none | open |
| MAP.121 | Every map shows and steps to its neighbouring regions, on one map engine | none | done, PR #3430 |
| MAP.122 | A Select mode on every galaxy view: Galaxy (blocks and sectors) or Star | none | done, PR #3553 |
| MAP.123 | Show or hide star types and phenomena, and set the luminosity floor, on the Galaxy and Sector Maps | none | done, PR #753 |
| MAP.124 | The Galaxy Map opens zoomed to fit all charted space | none | done, PR #696 |
| MAP.125 | Infinite zoom: one 3D interface from the galaxy down to a moon | none | done, PR #683 |
| MAP.126 | Show the orbital trajectories of selected objects in their frame of reference | none | done, PR #686 |
| MAP.135 | Selecting the first slab or wedge shows its bounds (bug) | none | done, PR #788 |
| MAP.136 | Binary stars pick as one system and their 3D orbits are drawn clearly (bug) | none | done, PR #797 |
| MAP.137 | Rogue planets are easy to see, with the right default filters (bug) | none | done, PR #797 |
| MAP.138 | Recenter the camera in every 3D view (bug) | none | done, PR #797 |
| MAP.139 | The Galaxy View uses its spare space: an info box with a Details link, and menu items | none | done, PR #2468 |
| MAP.140 | Double-click on a selected object goes there and opens its information | none | open |
| MAP.141 | Context around the selection: faint neighbours, and the sectors above and below | none | open |
| MAP.142 | Nebulae have fuzzy, fading boundaries | none | open |
| MAP.143 | Color sectors by their number of habitable locations | none | done, PR #2742 |
| MAP.144 | Replace `THREE.Clock` with `THREE.Timer` in `phenomenonrender.js` (bug) | none | done, PR #865 |
| MAP.145 | A sky and Galaxy Map drawing rule for neighbour galaxies | none | open |
| MAP.146 | Fly through the galaxy: scroll-zoom, double-click flight, distance-based visibility and a see-through near field | none | done, PR #2103 |
| MAP.147 | The Galaxy Map wire format: the investigation (done) and the record of what was built from it | none | done, PR #1849 |
| MAP.148 | The star visibility law: apparent-magnitude opacity, flux-based brightness and an on-screen limit from a histogram | none | done, PR #1160 |
| MAP.149 | The near field: depth fade, a see-through focus tube, drawing from inside a container, and picking that matches what is drawn | none | done, PR #1141 |
| MAP.150 | The free camera: wheel zoom to the cursor, double-click flight, and the observer inside, with the container named from position | none | done, PR #1382 |
| MAP.151 | The region data layer: exact-centred frame, aligned cells, per-level aggregates and slot-wrap ranges | none | done, PR #1539 |
| MAP.152 | Scale hand-offs: galaxy, sector and system cross-fade with hysteresis, and per-tile camera-relative origins | none | done, PR #1832 and #1838 |
| MAP.153 | Stars fade in with the zoom: a birth radius from each star's rank in its tile list (first client stage) | none | done, PR #998 |
| MAP.154 | Nested bright-star lists on the server, so every parent list is a subset of its child's | none | done, PR #1303 |
| MAP.155 | Other objects fade in too: point objects from level 8, a size ramp for cloud sprites, and stars that grow from a faint dot | none | done, PR #1253 |
| MAP.156 | View one layer or a range of layers top-down from the galaxy view, as a secondary option | none | open |
| MAP.157 | Trim the Galaxy Map tile JSON and serve it from prebuilt, precompressed bytes | none | done, PR #1112 |
| MAP.158 | A gentler tile prefetch and an IndexedDB tile cache instead of localStorage | none | done, PR #1123 |
| MAP.159 | Packed binary Galaxy Map tiles (quantised planes) with the nested tile lists, on one cache stamp bump | none | done, PR #1305 |
| MAP.160 | Quantise the Galaxy Map GPU buffers (deferred) | none | open |
| MAP.161 | Load the Galaxy Map faster on a first visit: bundle or preload its scripts | none | done, PR #1129 |
| MAP.162 | A sector holding scattered objects but never generated can still be opened, marked uncharted | none | done, PR #1085 |
| MAP.163 | Galaxy Map brightness scale: floor 2,500 L_sun at full zoom, then min/max scaling per zoom | none | done, PR #1029 |
| MAP.164 | Galaxy Map shows phenomena: black holes purple, neutron stars dark blue, sized by mass, dark colors out-shine brighter stars | none | done, PR #1029 |
| MAP.165 | Scattered phenomena store a mass so the Galaxy Map sizes them exactly | none | done, PR #1110 |
| MAP.166 | Galaxy Map "Dimmest star shown" says "every star" only when the view is complete | none | done, PR #1070 |
| MAP.167 | Gravity map: a heat map of the gravitational field inside a sector (the seminal feature of 9.0) | none | open |
| MAP.168 | The Sector Map gravity layer: coloured zones for pull, well depth and tidal strength | none | open |
| MAP.169 | The System Map gravity layer: the orbital plane as a heat map with Lagrange points and Hill spheres | none | open |
| MAP.170 | The Galaxy Map gravity layer: the galaxy potential and region aggregates, coarse and optional | none | open |
| MAP.171 | The Galaxy Map keys a massive phenomenon's visibility to its map luminosity (replaces "always lit") | none | open |
| MAP.172 | Nebula cell aggregates in the region pyramid | none | done, PR #2431 |
| MAP.173 | The Galaxy Map tile carries a nebula layer: regions and singles, nested like the star lists | none | done, PR #3609 |
| MAP.174 | Nebula cover in the stage-view cell statistics, for filled and unfilled cells | none | done, PR #2694 |
| MAP.175 | The dark-family nebula fill is invisible on the dark theme (bug) | none | done, PR #1203 |
| MAP.176 | Nebula regions on the Galaxy Map: grouped soft sprites in the dominant family colour | none | done, PR #3623 |
| MAP.177 | Color by "Nebula cover" on the Galaxy Map | none | open |
| MAP.178 | Field nebula clouds drawn one by one from the nebula table | none | open |
| MAP.179 | Scattered nebulae as markers from level 8 | none | open |
| MAP.180 | Picking, hover and fly-to for nebula regions and cover | none | open |
| MAP.181 | Dust colour for the dark nebula family on both themes | none | open |
| MAP.182 | Optional, default do not build: hosted H II and reflection nebulae for pre-placed O and B stars, drawn before their sector exists | none | merged into GEN.241 (Boss, 2026-10-15) |
| MAP.183 | Optional, default not built for the first release: near-view P and Q nebulae from the seeded field below about 250 pc | none | open |
| MAP.184 | "Mark rogue planets" on the Sector Map is redundant (bug) | none | done, PR #1282 |
| MAP.185 | Black holes vanish in the Sector Map (bug) | none | done, PR #1279 |
| MAP.186 | A binary pair is drawn in different places on the Galaxy Map and the Sector Map (bug) | none | done, PR #1283 |
| MAP.187 | The see-through "tunnel" in front of the camera leaves a visible halo and hides stars (bug) | none | done, PR #1288 |
| MAP.188 | Sector and slab borders add up to a near-opaque cube grid when zooming in (bug) | none | done, PR #1289 |
| MAP.189 | The Galaxy Map opens flown in to the created sectors | none | open |
| MAP.190 | Several separate regions of created sectors: which one the Galaxy Map opens on | none | done, PR #3445 |
| MAP.191 | Galaxy Map default view: no density shading, only nebulae shade the stars behind them | none | open |
| MAP.192 | See the objects whose positions changed in the last orbital update | none | open |
| MAP.193 | Show every object's axis of rotation, direction of rotation and axial tilt | none | open |
| MAP.194 | Nebulae and globular clusters open their own zoom window from the Galaxy Map | none | open |
| MAP.195 | Galaxy Map finer tiles leave out the stars a coarser tile already sent | none | open |
| MAP.196 | galaxystages.js reads the system= and object= parts of a sector link as digits only, so object links may fail (bug) | none | done, PR #1331 |
| MAP.197 | Remove the obsolete drill container code left by the map flip | none | done, PR #1455 |
| MAP.198 | The Galaxy Map outlines and shades the whole block a click would zoom into (bug) | none | done, PR #1470 |
| MAP.199 | Galaxy Map zoom and rotation share one centre, and the user can choose a new one (bug) | none | done, PR #1513 |
| MAP.200 | Galaxy Map on-map controls: a zoom slider, a direction star and a centre indicator (bug) | none | done, PR #1556 |
| MAP.201 | A selected sector that is not filled stays highlighted on the map until something else is selected (bug) | none | done, PR #1494 |
| MAP.202 | Galaxy Map sector and slab borders are nearly invisible unless a block is being selected (bug) | none | done, PR #1490 |
| MAP.203 | Stars in the sector view look too big next to the Galaxy Map, and the sector map should match (bug) | none | done, PR #1519 |
| MAP.204 | Regions under the galactic plane are hard to highlight and select on the Galaxy Map (bug) | none | done, PR #1692 |
| MAP.205 | Galaxy Map block ladder for shells | none | done, PR #1612 |
| MAP.206 | Galaxy Map drawing on shells | none | done, PR #1612 |
| MAP.207 | The Galaxy Map draws satellites, clusters and wanderers as points beyond the disk outline | none | open |
| MAP.208 | Stars and all phenomena on the three maps draw as points of light until near sector level, phenomena as bright as a star of similar mass | none | done, PR #1774 |
| MAP.209 | Opening a sector from the Galaxy Map flashes one view, reloads into another, and leaves the details panel empty (bug) | none | done, PR #FOLDED |
| MAP.210 | Nothing on the map can be selected, even in sector view (bug) | none | done, PR #FOLDED |
| MAP.211 | Hypervelocity stars show no accessible information: show their type, and a button to the phenomenon or the star system (bug) | none | done, PR #1787 |
| MAP.212 | A region selection tool on the Galaxy Map and the Sector Map, for admin actions | none | open |
| MAP.213 | Scene, maps and the JavaScript twin for trees: per-node zoom, log-radial overview and positions_at in JavaScript | none | done, PR #4107 |
| MAP.214 | Galaxy Map default view: no density, age or luminosity shading, only stars and markers for known space | none | done, PR #3284 |
| MAP.215 | Nebulae tint the stars seen through them, with soft translucent boundaries | none | open |
| MAP.216 | Admin page of the objects that moved in the last update, with the movement vector and a map ghost of the old position | none | open |
| MAP.217 | Map layer interface: one registry the Galaxy stage view, the Sector Map and the System Map share for a toggled layer (menu entry, modes, URL key, legend, picker layer, fade rule) | none | done, PR #2006 |
| MAP.218 | Galaxy Map density control: a readout for everyone, an editor for admins over any region, with preview, confirm and undo | none | open |
| MAP.219 | Galaxy Map marker layer: a kind registry with per-kind toggles, the zoom budget and picking, shipping the black-hole kind | none | done, PR #2860 |
| MAP.220 | Markers for habitable worlds and chosen classes, optional per kind | none | done, PR #3547 |
| MAP.221 | Station and vessel icons on the Galaxy Map | none | open |
| MAP.222 | Nebula zoom window: choosing a nebula on the Galaxy Map flies into a view of that nebula alone | none | open |
| MAP.223 | Globular cluster zoom window: choosing a cluster flies into a view of its own stars | none | open |
| MAP.224 | Cluster pages and the Galaxy Map use the parent link | none | open |
| MAP.225 | The wheel zooms toward the point under the cursor, with clearance | none | done, PR #2031 |
| MAP.226 | Double-click flight on the Galaxy Map: fly to the thing clicked, or half way toward empty space | none | done, PR #1859 |
| MAP.227 | The camera names its container: the header and the breadcrumb come from its position | none | done, PR #2032 |
| MAP.228 | The camera in the address: the URL and bookmarks hold the camera's target, distance and orientation | none | done, PR #2032 |
| MAP.229 | Map label engine: names placed by priority, with no overlaps, shown and hidden by zoom | none | done, PR #2980 |
| MAP.230 | Names of stars, systems and sectors on the Galaxy Map | none | done, PR #2980 |
| MAP.231 | The Sector Map and the System Map use the label engine | none | done, PR #3001 |
| MAP.232 | Names for the big features: spiral arms, bar, bulge, halo, satellites and clusters | none | open |
| MAP.233 | Names for nebulae on the Galaxy Map | none | done, PR #2980 |
| MAP.234 | Galaxy Map key: a Key section in the help dialog and an on-map key that names everything the map draws | none | done, PR #3978 |
| MAP.235 | Superbubbles as a registered Galaxy Map layer: a toggle, a hover name and radius, and a key row | none | open |
| MAP.236 | Superbubbles need clearer bounds and slightly more transparency on the Galaxy Map | none | done, PR #3739 |
| MAP.237 | The Galaxy Map's "Charted only" view is invisible again, and charted regions are shown by colour and opacity, not outlined blocks (bug) | none | done, PR #FOLDED |
| MAP.238 | Charted regions shown by colour and opacity, not outlined blocks: drop the outlines everywhere charted space is marked | none | merged into MAP.237 (Boss, 2026-10-16, 11:58 PM EDT: consolidation sweep) |
| MAP.239 | Raise the Galaxy Map tile root to the 262 kpc reach: `TILE_ROOT` in `galaxy/viewport.py` and `galaxymap3d.js` | none | done, PR #2152 |
| MAP.240 | Browser buffers hold camera-relative coordinates: territory dots, superbubbles, made rings, block prisms, stage views, system view dots and belts | none | done, PR #3605 |
| MAP.241 | Stars on the star maps should glow brighter, not grow bigger: size stays small, glow and colour carry luminosity and temperature (bug) | none | done, PR #3412 |
| MAP.242 | The Galaxy Map and its key draw what the family has: arms and bar only where they exist, a spheroid shell, a lens or a clumpy body | none | open |
| MAP.243 | The Galaxy Map cannot select or enter anything in the halo: the halo's scatter objects are not pickable at galaxy scale and a zoom toward one stops without opening a view (bug) | none | done, PR #2199 |
| MAP.244 | The maps and pages use the atmosphere colour where they showed the planet colour | none | done, PR #2805 |
| MAP.245 | The wheel zoom's remaining parts: the column-density ray fallback, the clearance from the nearest pickable thing and the half-way double click | none | done, PR #2103 |
| MAP.246 | A loose group of stars far above the disc on a fresh galaxy cannot be reached by zoom or navigation and cannot be selected (bug) | none | folded into MAP.243 (Bug Hunt reproduced it on Boss's seed, 2026-10-17, 4:30 PM EDT as the lane reported) |
| MAP.247 | The wheel zoom's column-density ray fallback: where the cursor ray hits nothing pickable, the point is the first sample where the integrated column density passes a threshold | none | done, PR #2229 |
| MAP.248 | The maps and pages follow the plan configuration: with the halo or the satellite galaxies off, nothing is drawn, picked or offered beyond the configured boundaries | none | open |
| MAP.249 | The Galaxy Map is laggy and slow, and changing zoom level slogs the whole page (bug) | none | done, PR #3420 |
| MAP.250 | When a map chooses which stars to show, halo and satellite-galaxy stars are the lowest priority and are always cut first, before anything else | none | done, PR #2155 |
| MAP.251 | When a map has to drop stars or objects, drop them evenly across the viewing space and keep the most massive and brightest of the galaxy last (bug) | none | done, PR #2375 |
| MAP.252 | Viewing only the filled sectors does not dim the phenomena (bug) | none | done, PR #2252 |
| MAP.253 | Map and scene read the parent chain: systemscene, body_positions, the browser twin orbitpositions.js and the systemmap branches drop the three origin conventions | none | open |
| MAP.254 | One shared info box: selecting anything on any view opens the same box, with a View button that opens the item's own page | none | done, PR #2514 |
| MAP.255 | The browser twin orbitpositions.js draws a wide pair as a circle: move it onto the stored orbit tree (bug) | none | done, PR #4015 |
| MAP.256 | The map shows the information box twice when an object is selected: once in the window and once to its right (bug) | none | done, PR #2852 |
| MAP.257 | Infinite zoom is still clunky: the zoom zones and segments other than sectors get in the way (bug) | none | open |
| MAP.258 | N1: the sector-block ladder above 243 sectors and drillChildren in JavaScript | none | done, PR #2896 |
| MAP.259 | N2: node aggregates for the pre-placed bright stars and scattered phenomena | none | done, PR #2937 |
| MAP.260 | N3: node read and cache: GET /galaxy/nodes, the node cache, the change feed and the warm-up | none | done, PR #2937 |
| MAP.261 | N4: the node cut on the page: galaxynodecut.js replaces the tile cut | none | done, PR #3070 |
| MAP.262 | N5: retire the cube tiles | none | done, PR #3070 |
| MAP.263 | S1: the camera is the only state: stages follow the camera, with no gesture gate, preset flights, per-stage zoom range or pan clamp | none | done, PR #2888 |
| MAP.264 | S2: nothing is clipped when you zoom: drop setWedgeClip and markClippedStars | none | done, PR #2903 |
| MAP.265 | S3: cells are a tool, not a display stage: the hover and selected sector outline, an optional grid tool and a closed-form pick | none | done, PR #2941 |
| MAP.266 | S4: the unresolved glow: a soft light drawn from the node aggregates, switchable in the layers menu | none | done, PR #3088 |
| MAP.267 | S5: one star buffer in relative-to-eye pairs, filled node by node (MAP.240 folded in) | none | done, PR #3114 |
| MAP.268 | S6: labels and selection on the cut: MAP.229 to MAP.233 on the projected-size query, selection by object ID across scales | none | done, PR #2980 |
| MAP.269 | Stars on the maps do not glow: a star further than about 3 sectors is a 1.5 px point with no halo and a mostly white core (bug) | none | open |
| MAP.270 | Far stars on the maps: use the blackbody colour and a luminosity halo (bug) | none | done, PR #2838 |
| MAP.271 | The bodylook record: the server builds one data-look record per body from stored columns | none | done, PR #2789 |
| MAP.272 | Shared body painter module: seeded noise, texture cache, paint queue and size tiers | none | done, PR #2789 |
| MAP.273 | Ocean, land and ice painter on stored data: sea level from the ocean fraction, caps from temperature, land hues from temperature and class family, with an estimate line in the info panel | none | done, PR #2789 |
| MAP.274 | Land colour from the bulk composition: oxidation and rock type | none | done, PR #4543 |
| MAP.275 | Gas and ice giant painter: bands, storms and thermal glow from rotation, temperature and the atmosphere colour, rings with resonance gaps from the moons | none | open |
| MAP.276 | Star painter: blackbody colour shared with star_spectrum, class texture, spots and glow | none | done, PR #2793 |
| MAP.277 | Atmosphere glow from pressure and composition, and the day and night terminator | none | done, PR #4515 |
| MAP.278 | Phenomena painted from classes: nebula tints, remnant shells by age, asteroid field grains, black hole and pulsar geometry, rogue planet and comet with today's data | none | done, PR #2798 |
| MAP.279 | Small bodies: lumpy meshes for moons and asteroids under about 400 km, from the stored shape | none | done, PR #4572 |
| MAP.280 | Cloud layer, terminator and cover on solid worlds; cloud-deck worlds paint the deck as the surface | none | open |
| MAP.281 | Giant chemistry bands once the band strengths exist | none | open |
| MAP.282 | Life tint: the pigment colour of the life chemistry, faint, on worlds rated Yellow or better; life that is not green mixes in the colour of its own evolutionary chemistry | none | open |
| MAP.283 | Comet colour from the coma species | none | open |
| MAP.284 | Night-side thermal glow on hot worlds, with no city lights until 9.2 | none | open |
| MAP.285 | Close-up tier: fragment-shader detail over about 500 px, shared with the sky view | none | open |
| MAP.286 | Supernova remnant shells drawn by age on the maps | none | done, PR #3746 |
| MAP.287 | The 3D system view takes the body look record: systemscene.py still colours by class_color | none | done, PR #3130 |
| MAP.288 | Do away with the info panel at the side of the galaxy and sector views: the popup alone, in a map window that takes the full width | none | done, PR #2875 |
| MAP.289 | The warm-up count of warm sectors always reads one answer built, because /galaxy/profile is answered fresh on each warm run for addressed sectors (bug) | none | done, PR #3743 |
| MAP.290 | The System Map draws only the close pair of a triple or larger system, not a scene for each extra star | none | done, PR #3054 |
| MAP.291 | Galaxy view map window does not use the full page width or the page height (bug) | none | done, PR #3132 |
| MAP.292 | Galaxy Map scroll zoom gets stuck and never reaches a sector, a system or a star (bug) | none | done, PR #3104 |
| MAP.293 | Colour by Habitable worlds no longer fills any sector since MAP.267: the legend shows but the sectors keep their normal colour (bug) | none | done, PR #3167 |
| MAP.294 | The sector page map fails to compile its vertex shader at phone width: camHigh, starLow, camLow and viewRotation are undeclared (bug) | none | done, PR #3167 |
| MAP.295 | The system infobox on the map does not scale to the size of the map window (bug) | none | done, PR #3299 |
| MAP.296 | The system infobox on the map lists the nearest systems under the location, each name a link (bug) | none | done, PR #3299 |
| MAP.297 | Remove the sector boxes from the maps: clicking goes to a location, never to a sector block (bug) | none | done, PR #3284 |
| MAP.298 | Every visible item that has room for a name gets a name, and the name is clickable to fly to that location (bug) | none | done, PR #3316 |
| MAP.299 | In a dense area no name overlaps another and each keeps a buffer; when names must be cut, the lower mass or brightness goes (bug) | none | done, PR #3754 |
| MAP.300 | Zoom reaches far enough in to single out one item in even a dense cluster (bug) | none | done, PR #3751 |
| MAP.301 | Scroll zoom heads for the most populated layer in front of the zoom direction so it does not fly past it (bug) | none | done, PR #3756 |
| MAP.302 | At sector-level zoom every kind of item appears: comets, rogue planets and the rest, not only stars (bug) | none | done, PR #3945 |
| MAP.303 | Stars with a habitable world (Yellow and above) are labelled with an icon on the map (bug) | none | done, PR #3496 |
| MAP.304 | Any object or icon on the map can be clicked: a single click selects it and shows its infobox in the map window, a double-click flies to it (bug) | none | done, PR #3284 |
| MAP.305 | Systems of possible interest carry one generic point-of-interest mark; a click lists what the points of interest are, and in dense areas one mark stands for many (bug) | none | done, PR #3494 |
| MAP.306 | The map menu filters points of interest by kind, so habitable worlds, gas giants and other kinds can be pulled forward (bug) | none | done, PR #3501 |
| MAP.307 | Stations and facilities become points of interest: interest reasons for the map mark and its menu filter | none | open |
| MAP.308 | A planet with moons has its own system map: the lunar orbits, the planet drawn from its rendered look at the moons' scale, and the atmosphere glow in the viewing space (bug) | none | done, PR #3327 |
| MAP.309 | Star map filter for trade routes: draw the routes (all, or by polity, length or hop count) as paths on the Galaxy Map | none | open |
| MAP.310 | Routes away from a selected object: select a system, spaceport or world on the star map and show only the trade routes leaving it | none | open |
| MAP.311 | Map markers and a menu filter for worlds with life, with sentient life and with a civilization, beside the habitable-world markers | none | done, PR #3562 |
| MAP.312 | Resource locations become points of interest: interest reasons for valuable asteroid belts and fields, icy bodies and giants, each with its own row in the map menu filter | none | open |
| MAP.313 | Galaxy Map resource richness layer: a heat map of the resource score by sector, by resource kind, with the richest sectors listed | none | open |
| MAP.314 | Galaxy Map labels the charted regions by zoom: the family word with Reach at 81 and 243, the Region name at 27 and 9, the corner or infobox below that (regions chain, step 3) | none | done, PR #3432 |
| MAP.315 | The map infobox shows the region and its parent: "Region: Vang Zarh (in Vang Reach)", each a link that flies the camera to the block (regions chain, step 4) | none | done, PR #3403 |
| MAP.316 | Regions are the nodes of the galactic-level map: a soft sphere with its name, coloured by share charted, replacing the sector boxes there, and a sortable charted-regions list (regions chain, step 3a) | none | done, PR #3445 |
| MAP.317 | Domain-sliced label budget: a share for inhabited places, uninhabited bodies, stars, nebulae and regions, and only drawn markers are labelled (regions chain, step 3c) | none | done, PR #3395 |
| MAP.318 | Population is a metric slot: the inhabited slice and the charted-regions list read a body's population once it is stored, and rank habitable worlds by habitability score until then (regions chain, step 3d) | none | done, PR #3483 |
| MAP.319 | flyToFit: flying to a region, star system, planet system of moons or nebula places the camera so the target fills the window, with no sector boxes drawn (regions chain, step 3b) | none | done, PR #3406 |
| MAP.320 | Distinct vector point-of-interest icons: habitable world, uninhabitable world (brown or grey, not red), inhabited world (coloured by habitability score), gas giant, asteroid field, sector and region | none | done, PR #3415 |
| MAP.321 | Point-of-interest markers by level: a star carries its best world marker in the star views, the planet or moon carries its own marker at the system level, and regions show inhabited tokens | none | done, PR #3482 |
| MAP.322 | Inside the lowest-level region, the sectors shown and labelled are the most densely filled ones by total mass (regions chain, step 3e) | none | done, PR #3513 |
| MAP.323 | The legend on the Galaxy Map covers the scale text: stack them so the legend sits above the scale and they never overlap (bug) | none | done, PR #3331 |
| MAP.324 | The soft glow around stars on the Galaxy Map floods the view: at sector and sector-block zoom it makes big white blobs and haze that hide the stars, the black holes and the selection, so finding anything is nearly impossible (bug) | none | done, PR #3359 |
| MAP.325 | The Galaxy Map infobox is too dense for a readable font: a block card lists ten rows plus a help paragraph in tiny text, and it sits over the zoom controls (bug) | none | done, PR #3331 |
| MAP.326 | Zooming in from the Galaxy Map does not carry smoothly into the star system map: on no system could Boss zoom from the galaxy into the system view as one continuous move (bug) | none | done, PR #3345 |
| MAP.327 | Galactic view draws a locked wireframe round the current region with no dimming (the sector part, dimming and wall icons, is done in #3422) | none | done, PR #3445 |
| MAP.328 | Rendered objects slowly rotate on their axis in every view, and every star does the same wherever stars can be seen (not in the galactic views) | none | done, PR #3573 |
| MAP.329 | System Map stepping reaches systems across a sector's edge | none | done, PR #3789 |
| MAP.330 | Star glow halos grow huge when stars are blurred and cover most of the screen: fade the glow out strongly with halo size and cut it to zero long before that (bug) | none | done, PR #3454 |
| MAP.331 | The purple wireframe cube on the Galaxy Map has no legend entry or label, so nothing says it outlines the charted region (bug) | none | done, PR #3535 |
| MAP.332 | The Galaxy Map info card (Block 10 1 63 and similar) is cut off at the bottom, so its region line and buttons are hidden below the More row (bug) | none | done, PR #3533 |
| MAP.333 | Zoom redesign Z0: fly-tos are no longer clamped to the zone and do not snap back after a region flight; the slider shows absolute log distance over one global range and works in flights; plus, minus, slider and pinch cross levels like the wheel; the wheel never dead-ends (bug) | none | done, PR #3556 |
| MAP.334 | Zoom redesign Z1: one global zoom policy, with pace and clearance taken from the distance to the nearest thing; zoom floor 1e-10 pc (1e-12 pc after Z5) | none | done, PR #3569 |
| MAP.335 | Zoom redesign Z2: sector and system scenes are chosen from the camera position, with one new server read (generated sectors near a point); the window fetch stops being a display source | none | done, PR #3647 |
| MAP.336 | Zoom redesign Z3: one goTo(target) for every fly with the verb Go to, a find-a-place box, a vertical log ruler with named clickable stops, crumbs as Go to, and Up one level | none | done, PR #3778 |
| MAP.337 | Zoom redesign Z4: delete the stage code and stage URLs and rewrite the stage browser tests | none | done, PR #4423 |
| MAP.338 | Zoom redesign Z5: planet and moon frames, the camera as frame plus offset, automatic near and far planes, and the MAP.240 leftovers | none | done, PR #4444 |
| MAP.339 | Black holes: 237 of 240 drawn at a 10 pc zoom are far from the view; show them by a near-field rule (bug) | none | done, PR #3636 |
| MAP.340 | Black hole crowd badges: crowded black hole marks merge into one mark with a count badge and a click list (bug) | none | done, PR #3640 |
| MAP.341 | Markers outside an opened sector do not dim like the stars do (a MAP.327 gap) (bug) | none | done, PR #3652 |
| MAP.342 | Server stops listing the stellar black holes in nodes over 27 sectors (optional, low priority) | none | dropped (measured by Core 3, 2026-10-29: saves under 15 KB gzipped once per view) |
| MAP.343 | Sector Map and System Map signs for the reasons to look (life, sentient life, civilization) and the chosen-class markers: body-level markers know only inhabited, gas giant, habitable, uninhabitable and belt | none | done, PR #3576 |
| MAP.344 | The civilization marker gains a technology band: filter and mark worlds by the technology level of their species | none | done, PR #3659 |
| MAP.345 | Map and search filters: a world filter for 'supports Class N' and 'has Class N' with a minimum score, from the life_class_scores index | none | done, PR #3829 |
| MAP.346 | Map browser tests fail on main since MAP.240 (#3605): the wheel moves the target off its line of sight, the window does not fill the page, the info popup does not scroll (bug) | none | done, PR #3639 |
| MAP.347 | Galaxy Map info popup redesign: it is hard to read and scrolls off the screen even on a maximized 1080p display (bug) | none | done, PR #3940 |
| MAP.348 | Galaxy Map hover preview: pointing at a place shows a small preview before the click | none | done, PR #3990 |
| MAP.349 | Galaxy Map panels: a Go to button on every container line (the region, the parent, each listed place)  | none | done, PR #4006 |
| MAP.350 | The Galaxy Map window fills the screen with only about 20 pixels of buffer at the top and bottom, so the whole window shows with no page scrolling (bug) | none | done, PR #3905 |
| MAP.351 | Galaxy Map reason per life class: a point-of-interest reason for each of the 15 life classes with map menu entries, backfilled for stored worlds | none | done, PR #3897 |
| MAP.352 | Back to a system in a different sector fails when pressed before the first Back has finished (bug) | none | done, PR #3946 |
| MAP.353 | Galaxy Map is always 95 percent of the window height, recomputed on every resize (follow-up to MAP.350) | none | done, PR #4139 |
| MAP.354 | Lagrange trojans T4: System Map icon at the centre of each Lagrange point that has objects, labelled with the host and point (for example Jupiter-L4); shares the geometry of the MAP.169 gravity layer | none | open |
| MAP.355 | The map on a sector page zooms out to the whole galaxy and pans into neighbouring sectors; limit it to the sector (zoom-out capped at 1.1 times the sector fit, pan kept inside the sector block) with the Show on Galaxy Map link at the top (bug) | none | done, PR #4279 |
| MAP.356 | The map legend covers too much of the map; make it a small corner of the map (bug) | none | done, PR #4393 |
| MAP.357 | Boss cannot reliably find and navigate from the galaxy down to filled and charted sectors and regions (bug) | none | done, PR #4393 |
| MAP.358 | The View system and Expand buttons of the sector info display cannot be seen or reached (bug) | none | done, PR #4393 |
| MAP.359 | The sector map and system map should match the size of the galaxy map for the browser window (bug) | none | done, PR #4393 |
| MAP.360 | Sector view highlights only the sector being viewed and zooms in as close as the window allows (bug) | none | done, PR #4393 |
| MAP.361 | Sector view shows every star system that fits when there is room, so it shows more of the sector than the galaxy view would | MAP.360 | done, PR #4393 |
| MAP.362 | Superbubbles are hidden on the map at load, with a Menu button to show them (bug) | none | done, PR #4393 |
| MAP.363 | Object page: draw the lumpy silhouette of an unrounded small body, as the System Map and the 3D system view do (bodyview.js) | none | done, PR #4594 |
| NAV.1 | Courses in "bearing mark mark" on nested frames | 28 (2026-09-30 18:14Z); 33 (2026-09-30 18:39Z to 2026-10-01 02:57Z) | done in 7.14.0, PR #130 (see note 4) |
| NAV.2 | Warp and fold speeds | 29 (2026-09-30 18:14Z); 34 (2026-09-30 18:39Z to 21:54Z) | done in 7.8.0, PR #121 |
| NAV.3 | One shared picker for the Galaxy, Sector and System displays | none | done, PR #724 |
| NAV.4 | Save a course | none | open |
| NAV.5 | Show a course on the Galaxy Map | none | open |
| NAV.6 | Courses that steer clear of gravity wells | none | open |
| NAV.7 | One reference for every object, with its parents | none | done, PR #688 |
| NAV.8 | Pages and anchors for stars, planets, moons and belts | none | done, PR #1095 |
| NAV.9 | Search and locate return references for every kind | none | done, PR #1139 |
| NAV.10 | Routing that scales past a few thousand systems | none | done, PR #814 |
| NAV.11 | Travel times for the system-to-system route too | none | done, PR #1099 |
| NAV.12 | No maximum hop length: a route always reaches the nearest star it can, across any number of sectors | none | done, PR #841 |
| NAV.13 | A picker module: select, step out, step in, step sideways | none | done, PR #564 |
| NAV.14 | One breadcrumb for every level | none | done, PR #571 |
| NAV.15 | Pick mode everywhere | none | done, PR #566 |
| NAV.16 | NAV endpoints can be any object | none | done, PR #697 |
| NAV.17 | A saved course record with both forms | none | done, PR #2553 |
| NAV.18 | Saved courses list: open, rename, delete and the Direct/Route switch, over the saved-list store | none | open |
| NAV.19 | Saved courses in the account (after USR.7) | none | open |
| NAV.20 | Draw the direct line and the route apart | none | done, PR #711 |
| NAV.21 | Fit the view to the whole course | none | open |
| NAV.22 | Courses inside a sector and a system | none | open |
| NAV.23 | Open a saved course on the map | none | open |
| NAV.24 | A keep-out radius for every kind of object | none | done, PR #719 |
| NAV.25 | Find the obstacles along a path | none | open |
| NAV.26 | Bend the path around keep-out spheres | none | open |
| NAV.27 | Moving bodies inside a system | none | open |
| NAV.28 | Show and save the adjusted course | none | open |
| NAV.29 | Replace "Nav from here" and "Nav to here" with "Start Here" and "End Here" while picking (bug) | none | done, PR #615 |
| NAV.30 | Hide "View phenomenon" and "View system" links while picking a course (bug) | none | done, PR #351 |
| NAV.31 | Galaxy wedges don't highlight on the navigation screens (bug) | none | done, PR #399 |
| NAV.32 | Every Galaxy and Sector Map control works on the navigation screens (bug) | none | done, PR #631 |
| NAV.33 | After picking one end of a course, stay at that zoom level (bug) | none | done, PR #615 |
| NAV.34 | Courses between separately generated areas find no route: the route graph splits into islands (bug) | none | done, PR #427 |
| NAV.35 | Mark jumps through unknown space in the route | none | merged into NAV.12 (the per-hop flag), PR #346 |
| NAV.36 | Unknown-space jumps drawn red and glowing | none | open |
| NAV.37 | An optional ship range for routes (open question) | none | dropped: conflicts with NAV.12 (a route always reaches the nearest star), PR #346 |
| NAV.38 | Every sector a straight line passes through | none | done, PR #357 |
| NAV.39 | Saved courses remember their unknown-space jumps and check them again | none | open |
| NAV.40 | Bookmarks can't be used to find the start or destination once a course pick has begun (bug) | none | done, PR #399 |
| NAV.41 | The NAV page's course map is too small to read (bug) | none | done, PR #457 |
| NAV.42 | Each route stop shows the course and distance to the next stop | none | done, PR #1153 |
| NAV.43 | Find everything within a distance of a place: the query and the API | none | done, PR #804 |
| NAV.44 | A "What's nearby" page: pick a place, enter a distance in parsecs, list what is there | none | done, PR #804 |
| NAV.45 | "What's within N pc" from the Galaxy Map and Sector Map | none | done, PR #2001 |
| NAV.46 | The NAV picker can't click galaxy wedges to zoom in (bug) | none | closed, already fixed; covered by a browser test (#568) |
| NAV.47 | Unknown-space jumps stop at scattered stars, black holes, neutron stars and quasars | none | done, PR #1389 |
| NAV.48 | Offer to generate the uncharted sectors that block a course | none | done, PR #1195 |
| NAV.49 | Waypoints: pick objects in Star select mode and plot a course through them, kept on the map until cleared | none | open |
| NAV.50 | Pick any object down to a moon as a NAV endpoint | none | done, PR #724 |
| NAV.51 | Courses route around asteroid fields | none | open |
| NAV.52 | Port `join_islands` and the k-d tree to cKDTree | none | done, PR #1172 |
| NAV.53 | A `cells_touching_sphere` helper, and pad `store.sectors_reached_by` by one edge (bug) | none | done, PR #843 |
| NAV.54 | One keep-out table for every kind of object: fields, supermassive holes and quasars, moons, stars in their own system, nebulae and remnants | none | open |
| NAV.55 | A tuning block for the keep-out knobs | none | merged into NAV.54 (Boss, 2026-10-16, 11:58 PM EDT: consolidation sweep) |
| NAV.56 | Census of overlapping keep-out spheres in a generated galaxy | none | open |
| NAV.57 | The Intergalactic Frame in navigation-frames.md and `navigation.py` | none | open |
| NAV.58 | Edit a course: add, remove, reorder and drag its waypoints, then save the edit | none | open |
| NAV.59 | Headline of 8.5: courses that bend around gravity wells, with waypoints, editable and saved courses | none | open |
| NAV.60 | What is nearby reads by the cell key: one statement for every kind, including objects with no sector and unbuilt scatter rows | none | done, PR #1977 |
| NAV.61 | Nearest-first by a radius ladder, and a dense-core guard with an early-stop count | none | done, PR #3685 |
| NAV.62 | Gravity of the stars and objects a cell leaves out: work out how a course and the orbital update account for the unresolved mass, derived from the density model | none | open |
| NAV.63 | Find the nearest trade route: by spaceport, and by where someone could put themselves in the path of a route | none | open |
| NAV.64 | Find the nearest resource location from a position or along a course, by resource kind and minimum score | none | open |
| OPS.1 | Build the version number from the category counters (item 80's version-scheme questions) | none (split from 80 by the renumbering) | done in the version-from-todo-counters PR |
| OPS.2 | Apache OOM-killed on the production server | 1 (2026-09-24 01:32Z to 02:02Z) | done in 5.47.0, PR #72 |
| OPS.3 | PowerShell installers and macOS-safe bash scripts | 50 (2026-09-30 20:43Z to 2026-10-01 02:57Z) | done in 7.16.0, PR #125 (see note 4) |
| OPS.4 | Generate page jobs on native Windows | 54 (2026-09-30 20:48Z); 55 (2026-09-30 20:48Z to 22:12Z) | done in 7.9.2, PR #124 |
| OPS.5 | Install and update check the log locations and say how to fix them | none | done, PR #288 |
| OPS.6 | Admin scripts accept impossible `--mysql-port` values (bug) | none | done, PR #457 |
| OPS.7 | Update asks to fill a wiped database with population data (bug) | none | done, PR #457 |
| OPS.8 | Update reloads Apache itself when run as root | none | done, PR #808 |
| OPS.9 | Multi-line messages lose their prefix in the debug log (bug) | none | done, PR #448 |
| OPS.10 | The galaxy seed and version at the top of every generation log | none | done, PR #391 |
| OPS.11 | Define "the same galaxy" and which versions stay reproducible | none | done, PR #362 (docs/design/generation/reproducible-galaxies.md) |
| OPS.12 | `generate.py reproduce`: a version and a seed rebuild a galaxy and check it | none | dropped (Boss, 2026-10-09 20:42Z) |
| OPS.13 | Every update records the version key, keeping the last 10 | none | done, PR #810 |
| OPS.14 | A warning when the running version key differs from the galaxy's | none | done, PR #876 |
| OPS.15 | Each update says whether it changes generated output | none | done, PR #3015 |
| OPS.16 | A daily maintenance script for Linux and macOS, with its lock helper | none | open |
| OPS.17 | Install and update set up the daily maintenance schedule | none | open |
| OPS.18 | Settings JSON backups kept in 18 slots: 7 daily, 4 weekly, 6 monthly, 1 yearly | none | dropped (Boss, 2026-10-09 20:42Z) |
| OPS.19 | The Generate jobs folder is /var/lib/planetgen while the checkout is /var/lib/planetGen (bug) | none | done, PR #525 |
| OPS.20 | Move the code base from zero dependencies to third-party open-source libraries | none | done, PR #767 |
| OPS.21 | Pinned third-party dependencies and a Redis server in install, update and CI | none | done, PR #437 |
| OPS.22 | Reorganize the code into importable Python packages with shared utility libraries | none | done, PR #473 |
| OPS.23 | A package layout plan for the reorganization | none | done, PR #435 |
| OPS.24 | Move the code into the new package layout, one package per PR | none | done, PR #473 |
| OPS.25 | Update tries to drop bright_star_blocks, a table that no longer exists (bug) | none | done, PR #442 |
| OPS.26 | Installer mixes apt's NumPy-1 builds (astropy, erfa, scikit-image) with pip's NumPy 2, so the requirements probe fails (bug) | none | done, PR #467 |
| OPS.27 | The Windows installer and docs point at Redis in WSL, not Memurai | none | done, PR #475 |
| OPS.28 | Generator epoch and battery digest: say whether two checkouts generate the same galaxy | none | done, PR #2654 |
| OPS.29 | Update reload: cover the gunicorn units and non-Apache hosts, and say what a reload aborts | none | done, PR #3192 |
| OPS.30 | A lock helper for the maintenance run | none | merged into OPS.16 (Boss, 2026-10-16, 11:58 PM EDT: consolidation sweep) |
| OPS.31 | Lint every example plist, XML and service file in CI | none | done, PR #1715 |
| OPS.32 | `examples/macos/org.planetgen.update.plist` is not well-formed XML, so the update daemon silently fails to install (bug) | none | done, PR #865 |
| OPS.33 | `.gitattributes` has no LF pins for the lock files and word list, so hashes differ between a Windows and a Linux checkout (bug) | none | done, PR #843 |
| OPS.34 | Windows Redis in WSL: fix the keep-alive advice and add a Start-RedisInWsl remedy | none | done, PR #981 |
| OPS.35 | A vendored-version lock file for the static libraries | none | done, PR #2990 |
| OPS.36 | Space and size checks measure the boot drive, not the drive holding the database (bug) | none | done, PR #868 |
| OPS.37 | A Generator version number: one sequential integer, shown in admin and the API | none | done, PR #2950 |
| OPS.38 | Two Redis dump files (dump.rdb and src/dump.rdb) are committed to main and should be removed and git-ignored (bug) | none | done, PR #967 |
| OPS.39 | Remove Windows support; keep only a simple docs/deployment/WINDOWS.md | none | done, PR #981 |
| OPS.40 | update.sh step 8 fails: setup-debug-log.sh loads the deleted util/appconfig.py (bug) | none | done, PR #989 |
| OPS.41 | Put the site in an "updating" state while update.sh runs long database migrations | none | done, PR #1780 |
| OPS.42 | Install numba: Python 3.12 minimum, numba required, and the web app never imports its compiled code | none | done, PR #1291 |
| OPS.43 | A fresh install gets markdown 3.11.1, which removes a pattern the page renderer names, and 74 tests fail (bug) | none | done, PR #1307 |
| OPS.44 | The docs site build (GitHub Pages, Jekyll) fails on every release commit and nothing in our tests would have caught it (bug) | none | done, PR #1307 |
| OPS.45 | Install planetGen into its own virtual environment and clear the packages out of the system Python | none | done, PR #1333 |
| OPS.46 | Make the installer check the Apache mod_wsgi venv rules: venv root, process group, Apache user access, matching Python | none | done, PR #1359 |
| OPS.47 | update.sh stops at the database migrate step in the venv: cryptography is missing from the requirements (bug) | none | done, PR #1346 |
| OPS.48 | Redis server settings: a documented example config and an update.sh warning (append-only file, noeviction, memory limit, loopback) | none | done, PR #1395 |
| OPS.49 | update.sh stops at the database migrate step with 'NoneType' object has no attribute 'is_auth_switch_request' (bug) | none | done, PR #1387 |
| OPS.50 | Keep the pymysql login workaround and make it compatible with the upstream fix | none | done, PR #1425 |
| OPS.51 | The maintenance-mode flag file cannot be written by the web user on the server, so maintenance mode fails with a PermissionError (bug) | none | done, PR #1741 |
| OPS.52 | Everything under the planetGen install is owned by www-data:www-data, and the example Apache2 config restricts directories | none | done, PR #1711 |
| OPS.53 | Before it changes anything, update.sh lists running and queued work and stops unless the operator agrees to cancel it, then nothing else can start until the update ends | none | done, PR #1781 |
| OPS.54 | The orbital update refuses to run while an update is running, with tests | none | done, PR #1782 |
| OPS.55 | Nightly full backup: a backup step in the daily maintenance run that takes every part into one archive, with a schedule setting and a nightly retention setting | none | open |
| OPS.56 | The generator epoch counter and the PR rule that enforces it | none | done, PR #1877 |
| OPS.57 | generator_battery_digest(): a database-free digest over 6 to 8 canned cases | none | done, PR #2628 |
| OPS.58 | The version history row and the galaxy's internal record carry the epoch, the battery digest and the environment | none | done, PR #2654 |
| OPS.59 | `generate galaxy --shell 0` does not stop on SIGTERM: the busy `worker --burst` child keeps running and the parent is left a zombie (bug) | none | done, PR #2124 |
| OPS.60 | The admin commands that did two things by flag: `lockouts` lists or lifts, `version_history --list` repeats `planetgen versions` | none | done, PR #2227 |
| OPS.61 | Remove the one-off tools `cli.render_parity` and `cli.dedupe`: they check and fix a schema and a naming rule that migrated long ago | none | done, PR #2227 |
| OPS.62 | The cloud environment setup script fails at its verify step because it imports `markdown`, which nothing uses (bug) | none | done, PR #3094 |
| OPS.63 | The setup script works on Ubuntu 24.04 and 26.04 LTS | none | done, PR #2712 |
| OPS.64 | The version code carries the database schema, generator epoch and API versions as well | none | open |
| OPS.65 | Debug mode wrote almost nothing to stdout: planetgen system --debug printed 3 lines because the generation runs in queue workers whose stdout is /dev/null (bug) | none | done, PR #2723 |
| OPS.66 | Remove the unused older executor in src/planetgen/queue/work.py: RQExecutor, _rq_task and the --debug FILE worker plumbing | none | done, PR #3093 |
| OPS.67 | Detached helpers (the background planetgen warm a run starts, an auto-started worker's own lifecycle) write the debug log file only: decide whether they relay to stdout | none | done, PR #3093 |
| OPS.68 | update.sh looked hung on a big data set download: no progress and no stall check (bug) | none | done, PR #2786 |
| OPS.69 | Remove all wiki functionality: the upload and download code, the wikitext renderer, the wiki routes, UI, settings and columns; the export to Markdown stays | none | done, PR #2963 |
| OPS.70 | Gate resumes after a restart, and a lane-start health check runs before every gate | none | done, PR #3192 |
| OPS.71 | Epoch and migration numbers are allocated by a script during the merge rebase, never renumbered by hand | none | done, PR #3695 |
| OPS.72 | Delete the dead tuning knobs and uncalled helpers the unused-values audit found (G7) | none | done, PR #3938 |
| OPS.73 | Delete the dead tuning knobs (HABITABLE_ZONE_BUFFER_AU, ROUND_RADIUS_KM, RADIUS_KM_SCIENTIFIC_NOTATION_THRESHOLD, STELLAR_PROTO_HELIUM_OFFSET_DEX, MAX_FLAVOR_PLANET, NAV_COURSE_DECIMAL_PLACES, STELLAR_U_TH_UNCERTAINTY) (G7, step 1) | none | done, PR #3938 |
| OPS.74 | Delete the uncalled helpers (class_allowed_in_zone, build_sector_configs, shell_chunks, draw_mass, mantle_mass_kg, build_nodes, get_galaxy_phenomena, register_card_renderer, read_config, settings_hash, mass_limit_argv) (G7, step 2) | none | done, PR #3938 |
| OPS.75 | update.sh ends with "the update failed at step 9/9, checking that the web app imports" when the import check passed and the real failure is the health check after the Apache reload; the site stays closed with a wrong reason (bug) | none | done, PR #3363 |
| OPS.76 | update.sh warns "uncommitted local changes in /var/lib/planetGen -- these will be overwritten: ?? state/" on every update although state/ is expected (bug) | none | done, PR #3363 |
| OPS.77 | update.sh ends with the site closed because its final /api/health check gets a 403 (a reemergence of OPS.75) (bug) | none | done, PR #3525 |
| OPS.78 | update.sh leaves the migrated site closed when the /api/health probe fails after a clean reload (127.0.0.1 answers 403, the site address answers HTTP 000) (bug) | none | done, PR #4256 |
| OPS.79 | Web app stops booting when a downloaded dataset is missing: update.sh took 52 minutes in stage 3 with the new code already serving and every page returned 500 (DataMissing mamajek_dwarf_sequence) (bug) | none | done, PR #4258 |
| OPS.80 | Every page of the live site hung (even a 404): the page freshness check built the whole galaxy shape, which takes minutes per web process and passes the 60 second daemon limit, so the daemon was killed in a loop (same hal incident as OPS.78 and OPS.79) (bug) | none | done, PR #4284 |
| OPS.81 | The Galaxy Map first load after a restart is slow: the web process works out the mean system mass of every stored cluster age on its first request (follow-up to OPS.80); update.sh should work it out once and store it for the web app (bug) | none | done, PR #4287 |
| OPS.82 | update.sh reloads the site and probes /api/health before the caches are warm, so the first requests after an update are slow or fail; run the cache warm-up in the foreground before the reload and let the probe wait after it (bug) | none | done, PR #4297 |
| OPS.83 | update.sh waits on the whole cache warm-up before the reload; start it in the background once the code, database, permissions and Redis are ready, wait only for its opening phase (at most 20 minutes), then reload Apache, wait 10 seconds and check health (bug) | none | done, PR #4300 |
| OPS.84 | update.sh says what every wait is for and for how long, the cluster-mass step shows a progress bar, and there are pauses of 5 seconds before and 10 seconds after the reload | none | done, PR #4303 |
| OPS.85 | update.sh is silent while it records the version key after the migration, so it looks hung (bug) | none | done, PR #4658 |
| OPS.86 | The update uses the API address from config.json and never tries localhost unless that is the configured address (bug) | none | done, PR #4664 |
| OPS.87 | The update shows a progress bar while it records the version key | none | done, PR #4660 |
| OPS.88 | The update splits the cache warm-up around the Apache reload: the opening phase before it, the rest after it with the site up | none | done, PR #4669 |
| OPS.89 | The health API returns the reason for maintenance mode (bug) | none | done, PR #4670 |
| OPS.90 | The update warm-up fails on the Galaxy Map page with 302 while the site is held in maintenance mode (bug) | none | done, PR #4683 |
| OPS.91 | The site sets and keeps its own secret key when none is configured | none | done, PR #4685 |
| OPS.92 | The warm-up shows which chunk of how many it is on, with a timed progress bar on a terminal | none | done, PR #4687 |
| OPS.93 | Recording the version key is faster: the bright-star scatter case skips table bounds it never reads | none | done, PR #4688 |
| PERF.1 | Generation at scale | none | done, PR #425 |
| PERF.2 | Cache so pages don't hit the database every request | 6 (2026-09-24 01:32Z to 02:18Z); 3 (2026-09-24 01:57Z to 02:02Z); 1 (2026-09-24 02:25Z to 2026-09-30 18:09Z); 8 (2026-09-30 18:14Z to 2026-10-01 05:05Z) | done in 7.56.0, PR #178 |
| PERF.3 | Estimate size and time before bulk generation | 86 (2026-10-01 03:15Z to 05:29Z) | done, PR #238 (stats in control schema v6) |
| PERF.4 | Second progress bar for slow plan layers | 88 (2026-10-01 03:26Z to 05:29Z) | done, PR #258 |
| PERF.5 | Scatter bright stars in stages | 89 (2026-10-01 03:36Z to 05:29Z) | done, PR #229 |
| PERF.6 | Rate-limit SQL calls, do more per call | 90 (2026-10-01 03:46Z to 05:29Z) | done: investigation, then PR #222, #223 and #225 (PERF.8 caps the writers) |
| PERF.7 | Parallelize sector and system generation | 91 (2026-10-01 03:46Z to 05:29Z) | done, PR #225 and #227 |
| PERF.8 | Parallel background work queue in the API | 92 (2026-10-01 03:46Z to 05:29Z) | done, PR #225 and #227 |
| PERF.9 | Weight the bright-star ETA by the shape of the galaxy | 93 (2026-10-01 04:50Z to 05:29Z) | done, PR #258 |
| PERF.10 | Record generation speed across a log scale of densities | 94 (2026-10-01 04:58Z to 05:29Z) | done, PR #238 |
| PERF.11 | Store each sector's expected and actual density | 95 (2026-10-01 04:58Z to 05:29Z) | done, PR #425 |
| PERF.12 | Check the schema once per process during generation | none | done, PR #222 |
| PERF.13 | Write each sector in batches | none | done, PR #222 |
| PERF.14 | Reserve a sector's names in bulk, safe with several writers at once | none | done, PR #222 |
| PERF.15 | Fewer queries per web page | none | done, PR #223 |
| PERF.16 | Search names without scanning every row | none | done, PR #223 |
| PERF.17 | A time limit on web database statements | none | done, PR #223 |
| PERF.18 | Run the GEN.30 bright-star backfill in parallel on the work queue | none | merged into PERF.108 (Boss, 2026-10-16, 11:58 PM EDT: consolidation sweep) |
| PERF.19 | Everything the API or web site starts runs on the work queue (investigate) | none | done, PR #492 |
| PERF.20 | Short-term caching through the work queue and API | none | done, PR #2154 |
| PERF.21 | Generation works with any worker count: the parallel path is built, used and tested (bug) | none | done, PR #371 |
| PERF.22 | On Python 3.12 a run hangs forever when a worker process dies (bug) | none | done, PR #371 |
| PERF.23 | The bright-star progress bar can end at 101% (bug) | none | done, PR #371 |
| PERF.24 | The work queue and web jobs on Redis with RQ | none | done, PR #517 |
| PERF.25 | The page cache on cachetools; the tile cache stays | none | done, PR #478 |
| PERF.26 | Size estimates don't match what generation stores (bug) | none | done, PR #470 |
| PERF.27 | On Python 3.9, `--workers=--` comes back as a list and crashes generate.py's option checks (bug) | none | done, PR #442 |
| PERF.28 | The console's second progress bar (the bright-star backfill after a sector run) never updates its ETA (bug) | none | done, PR #454 |
| PERF.29 | Record which runs a partly filled sector still needs | none | dropped (Boss, 2026-10-16, 3:00 PM EDT; dropped with PERF.100: no automatic resume in the Redis plan) |
| PERF.30 | Finish an interrupted block or sector run on the next start | none | dropped (Boss, 2026-10-16, 3:00 PM EDT; dropped with PERF.100: no automatic resume in the Redis plan) |
| PERF.34 | The site stays responsive during heavy generation jobs (bug) | none | done, PR #811 |
| PERF.31 | Investigate: where generation spends its time, from the plan to a finished galaxy | none | done, PR #1179 |
| PERF.32 | Generation performance stats: rates recorded per run, deleted on every new version | none | done, PR #905 |
| PERF.33 | Progress bars and ETAs from measured performance | none | done, PR #1553 |
| PERF.35 | An interval or chunk ledger for untouched sectors once block-first backfill lands | none | open |
| PERF.36 | Memory and request guard: list at most about 50,000 candidate cells, stream the rest, and refuse enumerations above about 2 million cells in a web request | none | done, PR #2070 |
| PERF.37 | `DecayingRate` starts from the first single completion, so the ETA is up to twice too long early in a run (bug) | none | done, PR #865 |
| PERF.38 | Cache fixes for the Galaxy Map under a fill: single-flight tile builds, a busy rule for the page cache, a deletion epoch in place of COUNT(*) | none | done, PR #1144 |
| PERF.39 | Every API job costs 2.4 s and 195 MB: import lazily and cap the burst workers | none | done, PR #1151 |
| PERF.40 | Two shared queues, a reserved interactive worker and a real "cancel now" | none | done, PR #1506 |
| PERF.41 | Stamp the tile and page caches with a TILE_FORMAT constant instead of the version (optional) | none | done, PR #2316 |
| PERF.42 | Warm the RQ worker before the fork: pre-import generation modules and build the bright-star table once | none | done, PR #841 |
| PERF.43 | Lazy word-salad names for phenomena named by object ID | none | done, PR #861 |
| PERF.44 | Compute object uids in Python and write them with the row | none | done, PR #863 |
| PERF.45 | Nearest-system links and containment as one later pass | none | done, PR #863 |
| PERF.46 | Planets and moons: set the position once per body | none | done, PR #1186 |
| PERF.47 | The PERF.31 benchmark records the buffer pool, table sizes and worker start-up cost | none | done, PR #2050 |
| PERF.48 | Low priority: a numeric-only INSERT formatter or C driver for bright_stars and phenomenon_scatter | none | done, PR #2050 |
| PERF.49 | Batch system-name reservation: remove the quadratic scan and the long-held registry locks (re-measure first) | none | done, PR #870 |
| PERF.50 | A progress bar inside one sector's save: workers report their sub-steps to the main process | none | done, PR #922 |
| PERF.51 | One progress mechanism for every sub-step: a bar starts by itself when a step is predicted to take over 15 seconds | none | done, PR #919 |
| PERF.52 | Admin generation-stats table is wrong (bug) | none | done, PR #994 |
| PERF.53 | Timing stats are skewed by layers and sectors that generated nothing (bug) | none | done, PR #1015 |
| PERF.54 | Separate generation-stats rows for the mass pass and the luminosity pass of the star scatter | none | done, PR #1284 |
| PERF.55 | One global progress bar for generation jobs that run in phases, with an ETA across all phases | none | done, PR #1131 |
| PERF.56 | Record how long each stage of a staged job takes, with the settings it ran with | none | done, PR #1077 |
| PERF.57 | Skip empty stretches in a galactic scatter by combining layers into growing groups | none | done, PR #1101; removed by PERF.58 (PR #1108) |
| PERF.58 | Object-first star sampler for the mass and luminosity passes (top priority, replaces PERF.57) | none | done, PR #1108 |
| PERF.59 | Share the ring inputs between the phenomena pass and the backfill rings (top priority) | none | closed as superseded by PERF.58, PERF.61 and PERF.63 (no PR) |
| PERF.60 | Large-mean Poisson helper for per-layer counts (top priority) | none | done, PR #1106 |
| PERF.61 | Object-first sampler for the phenomena pass, own prototype first | none | done, PR #1113 (measured, no gain; not merged) |
| PERF.62 | Interim early stop in the current layer-walking scatters: 100 empty layers in a row (tuning.SCATTER_DRY_LAYERS) | none | done, PR #1087; removed by PERF.58 (PR #1108) |
| PERF.63 | Vectorise the candidate work of the phenomena scatter with numpy | none | done, PR #1174 |
| PERF.64 | The Galaxy Map and the Systems list stop with "Took too long" and the whole site is sluggish (bug) | none | done, PR #1161 |
| PERF.65 | The text output of a multi-step job counts its steps, not its tasks: "Step 1 of 4" when the job has 12 (bug) | none | done, PR #1197 |
| PERF.66 | The whole-job progress bar should track total elapsed time and estimate a stage with no performance data from the earlier stages (bug) | none | done, PR #1197 |
| PERF.67 | Record how long every operation takes, from each stage up to the whole job, and use the records for ETAs | none | done, PR #1249 |
| PERF.68 | Measure the Galaxy Map tile queries on a big galaxy and make them fit the time limit (bug) | none | done, PR #1193 |
| PERF.69 | Store the Planets, Moons and Phenomena table counts like the Systems and Sectors counts (bug) | none | done, PR #1206 |
| PERF.70 | Sorting the Systems list by Sector or Octant on millions of rows must not sort them all (bug) | none | done, PR #1193 (built as a group-by-group read: no stored sector-name column, exact and needs no upkeep) |
| PERF.71 | Research how to keep slow database calls on large data sets from timing out: queue, split, or answer in parts | none | done, PR #1175 |
| PERF.72 | Research the cost of a sector's gravity grid and where to cut between exact and aggregated sources | none | open |
| PERF.73 | Cut the cost of writing the phenomenon rows to the database (now the biggest part of a scatter run) | none | done, PR #1209 |
| PERF.74 | Store a per-sector system count so the Sectors list does not count every system on each request | none | done, PR #1193 |
| PERF.75 | The jump-by-value box for the big data tables (keyset paging is done for the Systems, Uncharted stars and Phenomena tables) | none | done, PR #2154 |
| PERF.76 | Give each Galaxy Map tile piece its own time budget and serve an "incomplete" tile | none | done, PR #1199 |
| PERF.77 | Capped counts: "10,000 or more" where no stored count exists for a filter | none | done, PR #1199 |
| PERF.78 | A reserved warm worker for long admin operations, and the poll pattern for them | none | done, PR #1506 |
| PERF.79 | Two workers make the plan scatter about 12 times slower than one (bug) | none | done, PR #1206 |
| PERF.80 | "Save the sector paths" takes 24% of a small run and does not speed up with workers | none | done, PR #1251 |
| PERF.81 | Nearest-neighbour search in sector generation: a spatial grid instead of 1.8 million distance calls | none | done, PR #1236 |
| PERF.82 | Cut the cost of inserting a sector (about 17% of the fill) | none | done, PR #1352 |
| PERF.83 | The phenomena scatter takes about 5 minutes to start (bug) | none | done, PR #1259 |
| PERF.84 | Neighbour linking stalls on a large run (bug) | none | done, PR #1274 |
| PERF.85 | Generation stalls after the phenomena scatter and the layer count on the progress bar is wrong (bug) | none | done, PR #1265 |
| PERF.86 | Redo the work queue: drop the task and subtask logic and use Redis simply | none | done, PR #3656 |
| PERF.87 | Placement loop: test a new position against its neighbours as one array, not 919,000 distance calls | none | done, PR #1300 |
| PERF.88 | Sector-path integrator on numba: saving the sector paths takes 13 s instead of 34 s on the benchmark fill | none | done, PR #1873 |
| PERF.89 | Bright and massive star scatters: population densities as NumPy arrays (about 4 s of 15 to 20 s) | none | done, PR #1308 |
| PERF.90 | Memoise the ring and slot count functions (3.2 and 6.5 million calls in the fill) | none | done, PR #1300 |
| PERF.91 | Bulk random draws for the scatter layers (4 to 12 times on batches above 1,000 draws) | none | done, PR #2070 |
| PERF.92 | Save generation timing per object and per object type, as well as per operation, stage and job | none | done, PR #1371 |
| PERF.93 | api_jobs.ensure_workers counts dead workers as alive, so a killed API worker is never replaced (bug) | none | done, PR #1364 |
| PERF.94 | One Redis client module: timeouts, health checks and a cached pool, with JSON instead of pickle | none | done, PR #1383 |
| PERF.95 | Work queue step 1: the queue spine, two queues, on-demand workers with an idle timeout and the wanted-worker count | none | done, PR #1506 |
| PERF.96 | Work queue step 2: the run record, stage counters, event stream and the ETA function, as a library | none | done, PR #1553 |
| PERF.97 | Work queue step 3: the mass scatter end to end on the new queue | none | done, PR #1999 |
| PERF.98 | Work queue step 4: the remaining stages and the sector and backfill units on the new queue | none | done, PR #3654 |
| PERF.99 | Work queue step 5: the Generate page, Queue page and job log read Redis; cancel now through the stop command | none | done, PR #3030 |
| PERF.100 | Work queue step 6: delete the old machinery and the work_* tables (Alembic migration) | none | done, PR #3656 |
| PERF.101 | No progress is shown while the phenomena indexes are built, and a phase with no updates shows no time remaining (bug) | none | done, PR #1620 |
| PERF.102 | The migration progress display shows no estimated time remaining on either bar (bug) | none | done, PR #1626 |
| PERF.103 | A per-night cost ceiling on the orbital update: warn with the backlog size and move the rest the next night | none | done, PR #1769 |
| PERF.104 | Cache warm-up for the whole site: build every list, map, sector, system and search index cache after an update, then rebuild only what a change invalidates | none | done, PR #1826 |
| PERF.105 | Slow page requests queue a job that builds the answer into Redis, with a status route (the mechanism) | none | done, PR #2066 |
| PERF.106 | The plan command on the run record, with stage timing, the killed-worker handler and job timeouts | none | done, PR #2278 |
| PERF.107 | The sector units and the stages around them: start, neighborhood, sectors, link, population and settle | none | done, PR #2305 |
| PERF.108 | The backfill as queue units: chunks of at least 10,000 cells, one enumeration (folds in PERF.18) | none | done, PR #2305 |
| PERF.109 | Cache the Generate page estimate | none | done, PR #2050 |
| PERF.110 | XFetch refresh for costly entries that expire by age | none | done, PR #2050 |
| PERF.111 | `ResponseCache.put` without a database never trims the page cache: the cache grows past `max_entries`, and the thread test fails when its last put is one of those (bug) | none | done, PR #2974 |
| PERF.112 | The uncharted-star count runs a full COUNT(*) on every request: keep the count so the page and the API do not scan bright_stars each time (performance) | none | done, PR #2154 |
| PERF.113 | The neighbour link of a galaxy run: measure it, and split it into chunks of sectors if the time warrants | none | done, PR #2846 |
| PERF.114 | The population pass of a galaxy run: measure it, and split it into chunks of sectors if the time warrants | none | done, PR #2846 |
| PERF.115 | The settle stage of a galaxy run: measure it, and split it into chunks of sectors if the time warrants | none | done, PR #2794 |
| PERF.116 | The habitability re-score (`plan --habitability`): measure it, and split it into chunks of systems if the time warrants | none | done, PR #2794 |
| PERF.117 | The map-luminosity fill (`plan --map-luminosity`): measure it, and split it into id-range chunks if the time warrants | none | done, PR #2794 |
| PERF.118 | Generation as commands: every stage is its own command, and one wrapper builds the command list and queues each command as one queue item | none | done, PR #2159 |
| PERF.119 | The Generate page, the API and the one-off sector and system commands start their runs through the wrapper | none | done, PR #3664 |
| PERF.120 | Every command does one thing: the umbrella for reshaping every command and entry point | none | done, PR #3703 |
| PERF.121 | `plan` does only the skeleton, and the scatters become one `scatter` command with exactly one kind per run (mass, luminosity, phenomena, nebula field, band) | none | done, PR #2441 |
| PERF.122 | The map-luminosity fill and the habitability re-score are their own commands, not `plan` flags | none | done, PR #2441 |
| PERF.123 | `galaxy` stops being a bundle: sector generation, the neighbour link, the backfill, the population pass and the settle are each their own command, and `galaxy` is a step list the wrapper builds | none | done, PR #2865 |
| PERF.124 | `sector` and `phenomenon` only generate: the population pass leaves `sector` and the settle leaves `phenomenon` | none | done, PR #3194 |
| PERF.125 | The population pass and the territories pass are separate commands | none | done, PR #3199 |
| PERF.126 | `check-db`: each check is its own command and the report is the list the wrapper runs | none | done, PR #2211 |
| PERF.127 | `warm_map`: each cache family is its own command and `update.sh` builds the list | none | done, PR #2166 |
| PERF.128 | `benchmark` is a step list handed to the wrapper, not a sequence of its own | none | done, PR #3703 |
| PERF.129 | Robust tile caching: the first zoom levels built ahead, tiles stored compressed, built proactively around every visit, kept until the data changes | none | open |
| PERF.130 | Build the first several zoom levels of the Galaxy Map ahead of the first visitor: the opening view and every level down to a chosen depth | none | done, PR #2276 |
| PERF.131 | Store the tiles compressed and serve the stored bytes with their encoding, so a request does no compression | none | skipped (Boss, "Skip it", 11:27 AM EDT 2026-10-21) |
| PERF.132 | Build proactively on a visit: the surrounding tiles and the next zoom levels are built and cached in the background | none | done, PR #2298 |
| PERF.133 | Keep every tile until the data changes: invalidate by what changed, rebuild what a change touched, and never evict a valid tile for size | none | done, PR #2316 |
| PERF.134 | The admin page shows the tile cache and can clear it: tiles per level, bytes, built and visited counts, the build queue, and clear and build-now actions | none | open |
| PERF.135 | The Uncharted stars and Phenomena tables get a jump-by-value box for their numeric and scatter-id sorts | none | done, PR #2163 |
| PERF.136 | Measure the fill with wanderers and their moons, and decide whether the rogue cost is on by default | none | open |
| PERF.137 | Audit the remaining `LIMIT ? OFFSET ?` queries: say which serve a big table and move those to keyset paging | none | done, report limit-offset-final-audit-perf137-2026-10-31.md |
| PERF.138 | The nuclear cluster term made the extent majorant about 1000 times looser: the shell scatter takes over 250 s, not 5.6 s (bug) | none | done, PR #2297 |
| PERF.139 | One cache size setting, `cache.max_mb`, default 1024: tiles, stored lists and indexes share one disk budget | none | done, PR #2316 |
| PERF.140 | The index builder and its automatic refresh: each kind is rebuilt from what changed whenever generation, an orbital update or an edit changes it | none | done, PR #2370 |
| PERF.141 | The warm-up runs in phases after every generation, update or change, as deep as the budget allows | none | done, PR #2413 |
| PERF.142 | A stale job left in the database or the queue blocks every new run and Boss cannot proceed (bug) | none | done, PR #2345 |
| PERF.143 | A command resets all jobs and the queue: clears every run record, the active-run key, the waiting and started units and the old job tree | none | done, PR #2345 |
| PERF.144 | Every job queued to RQ gets a time limit from the best available estimate of how long it takes, never no limit, so a killed worker's entry expires | none | done, PR #2366 |
| PERF.145 | The run output repeats the stage: drop the stage-list line and the 'Stage 4 of 13: ...' line the engine logs under the wrapper's '=== Stage 4 of 13: ... ===' banner | none | done, PR #2416 |
| PERF.146 | A performance profile of every generation option and scatter, with a repeatable benchmark command, a report and a ranked optimization list: one of the last three 8.1 stages | none | done, PR #2559 |
| PERF.147 | A per-task progress bar with an ETA for any task expected to take, or taking, longer than 15 seconds | none | open |
| PERF.148 | Boss's stalled-job fix: find why generation jobs stall and fix it, with units of equal expected work shown in the log line | none | done, PR #2397 |
| PERF.149 | A map visit queues the missing sector as a run instead of generating it in the web request, and the Galaxy View is closed while a galaxy generation runs | none | done, PR #3723 |
| PERF.150 | Maintenance lane: warm-up, tile prefetch and rebuild and the database check run on their own queue with its own worker cap, a guaranteed minimum and a lend rule (F1) | none | done, PR #2588 |
| PERF.151 | A reserved interactive worker outside `wanted` at nice 5, and a long-job lane (cap 1) for API jobs expected to run over 60 seconds (F2) | none | done, PR #2562 |
| PERF.152 | Pause and resume reach the worker queues: a Redis flag the workers honour, written by the Queue page (F3) | none | done, PR #2562 |
| PERF.153 | An API job that finds a live generation run waits (a bounded retry, shown as waiting for run X) instead of failing Busy, and a Generate stage waits for an API job holding the key (F4) | none | done, PR #2562 |
| PERF.154 | One waiting list for runs in Redis (up to 8 runs, a setting), replacing next.json, the file lock and Busy (F5) | none | done, PR #3727 |
| PERF.155 | A unit-length ceiling: no unit of a stage runs longer than a set time (default 5 minutes) unless the stage is one unit; longer stages are split (F6) | none | done, PR #2923 |
| PERF.156 | The warm-up counts the search check answering 202 (a queued job) as a failure, so planetgen warm stops after phase 4 on a large galaxy (bug) | none | done, PR #2731 |
| PERF.157 | The warm-up reads 95.9 million rows for 1,044 deep tiles, about 92,000 rows a tile | none | done, PR #2731 |
| PERF.158 | scatter --band writes 348,743 single-row UPDATE statements and 2.4 million rows in 5.2 minutes | none | done, PR #2731 |
| PERF.159 | scatter --mass and --luminosity expand 65,536 shells serially in the supervisor, so the machine runs at 1.3 to 1.7 of 4 cores | none | done, PR #2745 |
| PERF.160 | fill-map-luminosity issues 3.1 million single-row updates: 162 s server CPU against 45 s client | none | closed, no change (the fill is already one statement per id block; Performance lane, 7:34 PM EDT 2026-10-23) |
| PERF.161 | bright_stars._reach_spots still calls tier_at for every candidate: batch it like the component pass | none | done, PR #2763 |
| PERF.162 | Member-tree drawing is 53 percent of a band unit, and speeding it changes values, so it needs a generator epoch and Boss's word | none | open |
| PERF.163 | The Stage 4 plan takes 264 s against a 3.71 s skeleton: find where the other time goes (bug) | none | done, PR #2576 |
| PERF.164 | The generation performance report holds only the small size: measure medium and large after PERF.163 and re-rank | none | done, PR #2623 |
| PERF.165 | The plan stage still stores 22,267 components row by row and spends about 10 s outside its unit (plan_shape twice, math check, queue start-up) (bug) | none | done, PR #2636 |
| PERF.166 | The mean-system-mass lookup bisects per age (0.2 ms each); a faster form changes the last bits and needs an epoch and Boss's word (bug) | none | open |
| PERF.167 | The phenomenon scatter on a large galaxy runs over an hour and stores over 66 million rows without finishing | none | done, PR #2679 |
| PERF.168 | test_sql_portability fails on main: benchmark_suite.py runs SHOW DATABASES, a reserved word on MySQL 8.0 and 8.4 (bug) | none | done, PR #2632 |
| PERF.169 | The plan stage spends 3.5 s of a 12.4 s unit on its report lines | none | done, PR #2648 |
| PERF.170 | The math check (the TEST.68 gate) runs at every update and from the admin suite, and no longer before every generation | none | done, PR #2696 |
| PERF.171 | The bright-star scatter keeps its cell and due-time indexes on while its rows go in: drop and rebuild them as PERF.167 did for the phenomenon scatter | none | done, PR #2763 |
| PERF.172 | The galaxy as one wrapper run: the population pass and the phenomenon work are units of the step list (the D2 slice of PERF.123) | none | done, PR #3697 |
| PERF.173 | A new-galaxy Generate run dies in stage 5 with "Execution not found in Redis": the worker warm outlives the 90 s execution record (bug) | none | done, PR #3034 |
| PERF.174 | Every work unit's database connection runs the schema check again: skip it when the schema version is current | none | done, PR #3449 |
| PERF.175 | benchmark --suite accepts --keep but drops the scratch databases on success: make it honour --keep (bug) | none | done, PR #3449 |
| PERF.176 | Storing the galaxy shape takes 52 percent of the plan (20 s of 38.6 s): find why and speed it up | none | done, PR #2877 |
| PERF.177 | A galaxy run died at the phenomena stage (14 percent) on a Redis out-of-memory (maxmemory 512mb, noeviction), the later stages never ran, and the log is a 70-line double traceback (bug) | none | done, PR #2836 |
| PERF.178 | The plan estimate said 33 GB and 2 min 33 s but the neighbourhood took 4 min 28 s, and the line reads "Refused ... Going ahead" (bug) | none | done, PR #3766 |
| PERF.179 | Progress bars are cut off in the log: "elapsed 0:0... remaining 0:00:..." (bug) | none | done, PR #3766 |
| PERF.180 | Single-unit stages already past the unit ceiling on big galaxies: clearing the bright stars (211 s on the medium galaxy) and the cache-warm families (one unit each, minutes on a large galaxy) | none | done, PR #3773 |
| PERF.181 | CRITICAL (bug): MySQL is killed out of memory during New galaxy Stage 7 (scatter --mass), 65,536 shells and 19,795 components (bug) | none | done, PR #3031 |
| PERF.182 | Finish an incomplete run: resume a New galaxy run that died partway from the stage it stopped in | none | done, PR #3711 |
| PERF.183 | Nightly New galaxy smoke: stage times and peak memory per stage, fail above a budget | none | done, PR #3078 |
| PERF.184 | Time estimates break when the server has no benchmark run: estimate from what is known and update as the system records how long each star, layer or unit takes (bug) | none | done, PR #3176 |
| PERF.185 | A hard host-wide cap on the total number of planetgen worker processes, across every lane and command, as a second guard against memory exhaustion | none | done, PR #3040 |
| PERF.186 | CRITICAL (bug): Cancel does not cancel, python workers keep saturating the host; add a working Stop everything that cancels all runs, kills all workers and blocks any new start until the owner says so (bug) | none | done, PR #3078 |
| PERF.187 | The phenomena scatter holds 6.4 GB in its own process while every other generate step stays under 300 MB: find what it keeps and bound it (high priority) (bug) | none | done, PR #3125 |
| PERF.188 | The bright-star fraction table build (star_population._bright_windows) takes 60 to 130 seconds per test in pure Python | none | open |
| PERF.189 | The block light rebuild at the end of every scatter reports its progress (Step bar, a log line, phase lines) (bug) | none | done, PR #3125 |
| PERF.190 | Measure the speed of block_lights.rebuild on a New galaxy with the scatter on, and speed it up if it dominates | none | done, PR #4028 |
| PERF.191 | The warm-up after an epoch change fails: "a stage queues at most 4000 units, not 4051" (bug) | none | done, PR #3120 |
| PERF.192 | The nebula field slab stage queues one unit per 50 pc slab with no cap handling, so a large galaxy can exceed the 4000-unit limit (bug) | none | done, PR #3731 |
| PERF.193 | The nebula block totals pass keeps the level-3 cover as Python dicts, about 2 GB on a 12 pc galaxy: hold it as numpy arrays (bug) | none | done, PR #3449 |
| PERF.194 | The phenomena stage progress line still says "186M to place": the weight is printed as a number of phenomena instead of work (bug) | none | done, PR #3726 |
| PERF.195 | Generation timing records are shown or no longer recorded: object_timing_stats and object_timing_samples (G5) | none | open |
| PERF.196 | Decide: show the object timing records or stop recording them (G5, step 1) | none | done, decided in reports/perf196-object-timing-decision-2026-10-31.md |
| PERF.197 | Carry out the object timing decision: show the records with a test, or delete the recording, the tables, the migration and the fingerprint entry (G5, step 2) | none | open |
| PERF.198 | Speed up util/detmath.matmul and cache equal gas_equilibrium.solve calls (finding 6) | none | open |
| PERF.199 | The words "rough guess" are removed from every ETA, time left and time estimate display (bug) | none | done, PR #3271 |
| PERF.200 | Galaxy Map node cut bookkeeping: about half the busy time in the profile is string building, geometry recomputation and cache re-insertion | none | open |
| PERF.201 | The neighborhood job test (test_api_job_runs, a_neighborhood_job_runs_its_sectors_through_the_wrapper) takes over 150 s because generation physics is slow on main: find the slow step and speed it up with identical values (bug) | none | done, PR #3633 |
| PERF.202 | Move each generating verb's starter handlers and start-time checks into the step units, so a verb is expressible as units in one run | none | done, PR #3857 |
| PERF.203 | The Generate page and job_runner follow one wrapper run: one command list, one run, with the page's progress, state and log files kept | none | done, PR #4100 |
| PERF.204 | An admin's visit to an uncharted cell queues a one-sector generation run by itself, and the page says queued (the later option 2 of PERF.149 rule 2) | none | open |
| PERF.205 | Keyset paging for the species list and the built phenomena list (slice 2 of the PERF.137 OFFSET audit) | none | done, PR #3916 |
| PERF.206 | Keyset walk for the octant and binary system group queries (slice 3 of the PERF.137 OFFSET audit) | none | done, PR #3929 |
| PERF.207 | planetgen planets and moons take --limit and --after (keyset), replacing OFFSET listing (slice 4 of the PERF.137 OFFSET audit) | none | done, PR #3935 |
| PERF.208 | A stored, indexed sector system count, density and quadrant so the remaining Sectors list sorts page by key (slice 5 of the PERF.137 OFFSET audit) | none | done, PR #3949 |
| PERF.209 | The generator battery digest costs about 30 s in every test process (its docstring says 3 to 5 s); the bright-star table build alone is 18 s | none | done, PR #4029 |
| PERF.210 | Cut the cost of GEN.992: re-finish only the bodies whose result can change after the star age is final | none | done, PR #4150 |
| PERF.211 | Phenomena and species list sorts that still sort the whole set per page: phenomena by type, descriptor, radius, sector and placed; species by homeworld and by system | none | done, PR #4009 |
| PERF.212 | The bright-star, phenomena and shell scatters take their params when their stage opens, not when the run is built | none | done, PR #4033 |
| PERF.213 | Link and settle take the started_at of the sectors stage when their stage opens | none | done, PR #3985 |
| PERF.214 | The sectors planner (_plan_sectors) runs in a worker inside the run, not in the starting process | none | done, PR #4013 |
| PERF.215 | An unfiltered planets or moons list pages by key: an index on planets (star_system_id, orbital_index) so the name-then-orbit order is served by the index | none | done, PR #3949 |
| PERF.216 | Find where the 50 pc spherical region count spends its time: the key statement is 3x to 10x slower than the research figure | none | done, PR #3992 |
| PERF.217 | The near list (near._systems_in_reach) reads the full row of every system in reach: list ids by the covering index and read names for the page only | none | done, PR #4019 |
| PERF.218 | run_plan._draw_sector_bands recomputes the bright-star band table for every capped cell: memoise it per (shell, luminosity floor) | none | done, PR #4075 |
| PERF.219 | The Systems list with a star type filter and the sector or octant sort pages with OFFSET only (a G filter matches about a tenth of the galaxy); move it to keyset paging, with a MySQL measurement at a full galaxy against the 50 ms target | none | open |
| PERF.220 | Cap the workers of the bright-star and massive-star scatter stages by the memory available, so a run is not killed for running out of memory (bug) | none | done, PR #4560 (scatter workers capped at 60 percent of host memory at 600 MB each; Boss's kernel log confirmed the out of memory kill, 450 to 900 MB seen) |
| PERF.221 | A failed or killed unit names the wrong stage: the message says "unit 101 of stage 5" for a unit of stage 8 of 14 (bug) | none | done, PR #4576 |
| PERF.222 | Generation hosts: MySQL redo log is too small (InnoDB warns the checkpointer lags behind); setup recommends and checks innodb_redo_log_capacity (bug) | none | done, PR #4568 |
| PERF.223 | Protect MySQL (and Redis) from the kernel out of memory killer on a generation host: setup docs and update.sh recommend OOMScoreAdjust=-900 for mysql.service and a memory limit on the workers (bug) | none | done, PR #4570 |
| PERF.224 | Cheaper bright-table bound: fewer compositions per bin, so the star tables build about 1.7 times faster | none | dropped (Star 2 12:14 AM EDT as the lane reported: three compositions per bin is not a certified bound, the caches are not the lever, so as written it cannot be done) |
| PERF.225 | The galaxy shape's cluster field is built by a queued job in small saved chunks and kept between requests, so no web request rebuilds it (about 56 seconds) | none | done, PR #4641 |
| PERF.226 | Every generation stage has a measured worker memory, and the worker cap uses it, so no stage past the scatter can run the host out of memory | none | open |
| PERF.227 | Starting a single scatter part from the Generate console on a fresh download kills the web daemon: five 'Truncated or oversized response headers' errors, an error 500, and the job reads 'stopped without finishing' (bug) | none | open |
| PERF.228 | The rescore-habitability stage fails with a MySQL deadlock (error 1213) when four workers run its units at once (bug) | none | open |
| PERF.229 | Every supergiant and bright-giant star type times out at 900 s in the fuzz test (bug) | none | open |
| PERF.230 | Faster isochrone interpolation for the star tables: edge_point and end_state are about 60 percent of a table build, output identical | none | open |
| POP.1 | Government ownership of systems | 12 (2026-09-24 01:32Z to 02:18Z); 9 (2026-09-24 01:57Z to 02:02Z); 7 (2026-09-24 02:25Z to 05:38Z); 6 (2026-09-24 02:53Z to 2026-09-30 17:58Z); 15 (2026-09-30 16:44Z); 16 (2026-09-30 16:49Z); 17 (2026-09-30 16:51Z to 18:09Z); 34 (2026-09-30 18:14Z); 39 (2026-09-30 18:39Z to 18:41Z); 59 (2026-09-30 19:02Z to 19:17Z); 63 (2026-09-30 20:01Z to 20:27Z); 46 (2026-09-30 20:07Z); 50 (2026-09-30 20:08Z to 20:48Z); 51 (2026-09-30 20:43Z to 2026-10-01 04:37Z) | done in 7.49.0, PR #169 (optional since 7.58.2, PR #180) |
| POP.2 | Names for dominant species on living worlds | 13 (2026-09-24 01:32Z to 02:18Z); 10 (2026-09-24 01:57Z to 02:02Z); 8 (2026-09-24 02:25Z to 05:38Z); 7 (2026-09-24 02:53Z to 2026-09-30 17:58Z); 16 (2026-09-30 16:44Z); 17 (2026-09-30 16:49Z); 18 (2026-09-30 16:51Z to 18:09Z); 35 (2026-09-30 18:14Z); 40 (2026-09-30 18:39Z to 18:41Z); 60 (2026-09-30 19:02Z to 19:17Z); 64 (2026-09-30 20:01Z to 20:27Z); 47 (2026-09-30 20:07Z); 51 (2026-09-30 20:08Z to 20:48Z); 52 (2026-09-30 20:43Z to 2026-10-01 04:37Z) | done in 7.49.0, PR #169 |
| POP.3 | Database of spacefaring species | 14 (2026-09-24 01:32Z to 02:18Z); 11 (2026-09-24 01:57Z to 02:02Z); 9 (2026-09-24 02:25Z to 05:38Z); 8 (2026-09-24 02:53Z to 2026-09-30 17:58Z); 17 (2026-09-30 16:44Z); 18 (2026-09-30 16:49Z); 19 (2026-09-30 16:51Z to 18:09Z); 36 (2026-09-30 18:14Z); 41 (2026-09-30 18:39Z to 18:41Z); 61 (2026-09-30 19:02Z to 19:17Z); 65 (2026-09-30 20:01Z to 20:27Z); 48 (2026-09-30 20:07Z); 52 (2026-09-30 20:08Z to 20:48Z); 53 (2026-09-30 20:43Z to 2026-10-01 04:37Z) | done in 7.49.0, PR #169 |
| POP.4 | Younger and older civilizations | 15 (2026-09-24 01:32Z to 02:18Z); 12 (2026-09-24 01:57Z to 02:02Z); 10 (2026-09-24 02:25Z to 05:38Z); 9 (2026-09-24 02:53Z to 2026-09-30 17:58Z); 18 (2026-09-30 16:44Z); 19 (2026-09-30 16:49Z); 20 (2026-09-30 16:51Z to 18:09Z); 37 (2026-09-30 18:14Z); 42 (2026-09-30 18:39Z to 18:41Z); 62 (2026-09-30 19:02Z to 19:17Z); 66 (2026-09-30 20:01Z to 20:27Z); 49 (2026-09-30 20:07Z); 53 (2026-09-30 20:08Z to 20:48Z); 54 (2026-09-30 20:43Z to 2026-10-01 04:37Z) | done in 7.49.0, PR #169 |
| POP.5 | Population pages: species, polities, dominant species, territory (unnumbered in TODO.md) | none | done after 7.58.2, PR #184 |
| POP.6 | Territories overlay on the Galaxy Map (unnumbered in TODO.md) | none | done in 7.54.0 (PR #179) and 7.58.1 (PR #181) |
| POP.7 | Tech levels for technological species | none | done, PR #1461 |
| POP.8 | A tech-level design from the six domain indices | none | done, PR #1422 |
| POP.9 | A tech level generated for every technological species | none | done, PR #1461 |
| POP.10 | Facility types: programmable, picked from a dropdown, with affiliation and Green/Yellow/Red ratings | none | done, PR #1524 |
| POP.11 | Scatter space vessels and stations for the races in their territories | none | open |
| POP.12 | List every species that could live in a given environment, and give every species its environment limits | none | done, PR #1554 |
| POP.13 | A species editor: create and edit species and their environment limits | none | open |
| POP.14 | A race is only placed where it can survive: a placement check and a consistency check over the data | none | done, PR #2651 |
| POP.15 | A per-item technology level table for the species fit: which equipment item a species can field at which tech level | none | done, PR #2766 |
| POP.16 | Design: settlements and nation-states on a planet with a technological civilization (the tech gate, nation count and dominance laws, the city size law, the city to nation draw, space ports), with the sources and stand-in tiers | none | open |
| POP.17 | Which worlds get settlements: the tech gate and the planet inputs (habitable land area, carrying capacity) for a species whose civilization is dated | none | open |
| POP.18 | Nation-states on a planet: how many, their land shares, names, governments, power and which are dominant | none | open |
| POP.19 | Cities: the capital and the ten largest cities of a planet, their populations, names and positions on land | none | open |
| POP.20 | Cities belong to nation-states: a probability draw weighted by city size and nation dominance, with the capital in the dominant nation | none | open |
| POP.21 | Station and space port scatter: orbital stations, space ports and starbases placed for each nation and polity from the tech gate, as stored facilities | none | open |
| POP.22 | The settlements stage: its own command, its place in the galaxy run and the rescan rules | none | open |
| POP.23 | Pages and API for nations, cities and ports: the planet page sections, nation and city pages, species and polity links | none | open |
| POP.24 | Settlements become dynamic: nations and cities grow, borders and dominance shift and cities change hands over time | none | open |
| POP.25 | Option: how many space-faring governments the galaxy has (a number, or natural), on the Generate page, the CLI, the API and the stored plan | none | open |
| POP.26 | Generate the requested number of space-faring governments: force civilizations where chance gives too few, hold back extra ones where it gives too many, homeworlds spread apart | none | open |
| POP.27 | Territories of several worlds and sectors: every government holds a minimum, its homeworld is always inside, and the check is a test | none | open |
| POP.28 | Show the governments: the option and result on the Generate page and galaxy page, a government page with its worlds, systems and sectors, and the territories overlay legend | none | open |
| POP.29 | Design: the world population simulation from habitability scores and the trade routes between spaceports (partner selection, hop runs, the global route list, the map filters, the nearest-route query) | none | open |
| POP.30 | World population simulation from habitability scores: population, carrying capacity and growth state for every world with a civilization, as its own command and galaxy stage | none | open |
| POP.31 | Population on the pages and API: world population, capacity and share of capacity on the planet, species and polity pages | none | open |
| POP.32 | Trade partners: which spaceports trade with which, by the gravity model on population and distance, with a cap per port and hubs | none | open |
| POP.33 | Generate the trade routes: a hop-based run on the navigation graph between each pair of spaceports, saved in the one global route list | none | open |
| POP.34 | The global trade route list: a page and an API every account can read, with filters by polity, route length, hop count and spaceport | none | open |
| POP.35 | Dynamic world populations: populations step forward with the nightly update from the stored growth state, so cities and capacity follow the habitability scores over time | none | open |
| POP.36 | Dynamic trade routes: routes grow, shift and close as populations, ports and territories change, using the same hop pathfinding | none | open |
| POP.37 | Facility ratings read Suitability, the food rung and the life scores: Health and Housing use only PHI-4 today | none | done, PR #3276 |
| POP.38 | Species traits and facility needs read the homeworld and host: dose, pressure, air, thermal range, day length and the equipment kit (G1) | none | open |
| POP.39 | Facility constants are registered with a tier and tweak level, and a facility's richness reads the host's stored composition (G6) | none | open |
| POP.40 | Species traits take radiation tolerance, pressure and air tolerance, thermal range and day length from the homeworld (G1, step 1) | none | open |
| POP.41 | Facility ratings (health, housing, maintenance, resources) read the host's dose, pressure, air and temperature (G1, step 2) | none | open |
| POP.42 | Facility need and upkeep read the equipment kit the host requires (eq_body_*, eq_rad_suit*, eq_gas_filter, eq_co2_scrub, eq_dust, eq_supply*, eq_shelter*, outside_hours) (G1, step 3) | none | open |
| POP.43 | Settlement and facility pages say why: the species traits and the host's dose, air and kit behind each rating (G1, step 4) | none | open |
| POP.44 | Register the facility constants (_WEIGHTS, FACILITY_HOST_RICHNESS, FACILITY_KIND_SIZE, FACILITY_RATING_*) in the GEN.725 registry with tier and tweak level (G6, step 1) | none | open |
| POP.45 | A facility's richness reads the host's stored composition and resource scores instead of a per-type constant (G6, step 2) | none | open |
| SEC.1 | Lock out an IP after failed logins | 61 (2026-10-01 01:15Z to 05:29Z) | done, PR #220 |
| SEC.2 | Security audit findings of 2026-09-30 (new parent) | none | done in 7.5.0, PR #112, except SEC.17 to SEC.19 |
| SEC.3 | Seeded admin/password login claimable | 39 (2026-09-30 19:02Z to 20:27Z) | done in 7.5.0, PR #112 |
| SEC.4 | Web-user compromise can become root via the installer | 40 (2026-09-30 19:02Z to 20:27Z) | done in 7.5.0, PR #112 |
| SEC.5 | HTML pages have no rate limit | 41 (2026-09-30 19:02Z to 20:27Z) | done in 7.5.0, PR #112 |
| SEC.6 | /tmp fallback for jobs and tile cache can be hijacked | 42 (2026-09-30 19:02Z to 20:27Z) | done in 7.5.0, PR #112 |
| SEC.7 | Debug log is mode 0666 | 43 (2026-09-30 19:02Z to 20:27Z) | done in 7.5.0, PR #112 |
| SEC.8 | Login time reveals which usernames exist | 44 (2026-09-30 19:02Z to 20:27Z) | done in 7.5.0, PR #112 |
| SEC.9 | Changing credentials leaves other sessions logged in | 45 (2026-09-30 19:02Z to 20:27Z) | done in 7.5.0, PR #112 |
| SEC.10 | Sector wiki link accepts any scheme | 46 (2026-09-30 19:02Z to 20:27Z) | done in 7.5.0, PR #112 |
| SEC.11 | /api/databases lists the control schema | 47 (2026-09-30 19:02Z to 20:27Z) | done in 7.5.0, PR #112 |
| SEC.12 | Public endpoints return raw database errors | 48 (2026-09-30 19:02Z to 20:27Z) | done in 7.5.0, PR #112 |
| SEC.13 | CSRF token not tied to the login session | 49 (2026-09-30 19:02Z to 20:27Z) | done in 7.5.0, PR #112 |
| SEC.14 | New API key's value rides in the flash cookie | 50 (2026-09-30 19:02Z to 20:27Z) | done in 7.5.0, PR #112 |
| SEC.15 | Nothing sets config.json's permissions | 51 (2026-09-30 19:02Z to 20:27Z) | done in 7.5.0, PR #112 |
| SEC.16 | Hardening (HSTS, lock file, login backoff, input bounds) | 52 (2026-09-30 19:02Z to 20:27Z); 39 (2026-09-30 20:07Z to 22:53Z) | done: HSTS in 7.5.0 (PR #112); rest as SEC.17 to SEC.19 |
| SEC.17 | Per-username login backoff (cited as "39a") | 39a (2026-09-30 22:33Z to 22:51Z) | done in 7.14.1, PR #131 (see note 8) |
| SEC.18 | Upper bounds on admin generation inputs (cited as "39b") | 39b (2026-09-30 21:37Z to 22:27Z) | done in 7.10.2, PR #118 (see note 8) |
| SEC.19 | Hashed lock file for pip dependencies (cited as "39c") | 39c (2026-09-30 21:39Z to 21:54Z) | done in 7.7.0, PR #119 (see note 8) |
| SEC.20 | Log every failed and locked login with its address | none | done, PR #217 |
| SEC.21 | Keep the per-username backoff in the control database | none | done, PR #220 |
| SEC.22 | Trusted-device cookie so lockouts can't shut out the real admin | none | done, PR #221 |
| SEC.23 | Wrong current passwords on /account aren't counted (bug) | none | done, PR #221 |
| SEC.24 | Refuse common and breached passwords | none | done, PR #221 |
| SEC.25 | Check the password hashing cost and re-hash on login | none | done, PR #221 |
| SEC.26 | Two-factor sign-in (TOTP) for admins | none | done, PR #221 |
| SEC.27 | A fail2ban filter and jail for login brute force (was: a fail2ban recipe in the deployment docs) | none | done, PR #221 |
| SEC.28 | An always-on log in the standard log location | none | done, PR #217 |
| SEC.29 | Two-step sign-in on pyotp, QR codes on segno | none | done, PR #612 |
| SEC.30 | Login and request rate limits on Flask-Limiter with Redis storage | none | done, PR #643 |
| SEC.31 | Signing in as admin works but shows a "form expired" error (bug) | none | done, PR #457 |
| SEC.32 | Argon2id password hashing, a 1,024-character password limit and the `__Host-` cookie prefix | none | done, PR #2667 |
| SEC.33 | Release gate: a full security code audit of the whole code base, done before every major or minor version change (8.1 first) | none | open |
| SEC.34 | Upgrade Pillow from 12.3.0 to 13.0.0 in both lock files (37 advisories fixed in 13.0.0) | none | done, PR #3141 |
| SEC.35 | The picture upload opens the file with every format Pillow knows before refusing all but PNG, JPEG and ICO: restrict the formats on open | none | done, PR #3141 |
| SEC.36 | Avatars may be GIFs, animated GIFs included: a new avatar picture kind with its own format list and decoder limits | none | done, PR #3192 |
| USR.1 | User accounts | none | open |
| USR.2 | Accounts with roles: user, moderator, admin and Owner | 64 (2026-10-01 02:13Z to 05:29Z) | open |
| USR.3 | SMTP settings in the admin config | 65 (2026-10-01 02:13Z to 05:29Z) | open |
| USR.4 | Invite-only sign-up by unique link | 66 (2026-10-01 02:13Z to 05:29Z) | open |
| USR.5 | Email loop for setting and resetting passwords | 67 (2026-10-01 02:13Z to 05:29Z) | open |
| USR.6 | Owner transfer | 68 (2026-10-01 02:13Z to 05:29Z) | open |
| USR.7 | A user-level interface with bookmarks | 69 (2026-10-01 02:13Z to 05:29Z) | open |
| USR.8 | Every signed-in user can generate a one-off system | none | open |
| USR.9 | `seo.privacy_note` text, "download my data" and "delete my account" on the account page | none | open |
| USR.10 | Users add their own objects (ships and stations) and track their movement on the map | none | open |
| USR.11 | A "home" field on user accounts: where the Galaxy Map zooms to when it opens | none | open |
| USR.12 | Permission catalogue and one `can()` check: user, moderator, admin and Owner as data | none | open |
| USR.13 | Moderator role: assign and remove it on the admin accounts page, with the audit entry | none | open |
| USR.14 | Move every route, page and command from "is admin" checks to permission checks | none | open |
| USR.15 | Submit an edit request: a user proposes an edit and the plan engine computes what it would change | none | open |
| USR.16 | Request queue: who sees what | none | open |
| USR.17 | Approve, reject or send back a request, by role: moderators decide user requests, admins decide moderator requests | none | open |
| USR.18 | Approver screen: the requester's note, the plan preview and the decision buttons | none | open |
| USR.19 | Comments on a request between the requester and the approver | none | open |
| USR.20 | "My requests" in the user interface: state of each request and withdraw | none | open |
| USR.21 | Request lifecycle: stale, withdrawn and expired requests | none | open |
| USR.22 | Request limits: a per-user cap on open requests, a rate limit, and blocking a user from requesting | none | open |
| USR.23 | Email notices: a request was submitted, decided or sent back | none | open |
| USR.24 | Audit view for roles and requests: who changed a role, who asked, who decided | none | open |
| USR.25 | Headline of 10.0: roles and permissions, user edit requests and moderators | none | open |
| USR.26 | Roles in sessions, API keys and route guards; the Owner rules; the accounts CLI | none | open |
| USR.27 | Admin Users page: list accounts, promote, demote, disable and delete | none | open |
| USR.28 | Account tokens: one issue-and-redeem library for invite, first-password, reset, email-change and owner-transfer links | none | open |
| USR.29 | Vessel and station objects: add, move along a course and list, with ownership | none | open |
| USR.30 | Add a ship or station from a page and see it and its track on the maps | none | open |
| USR.31 | The Owner account stub: one Owner row and a role check, enough for the Owner-tier settings of ADM.86 in 8.1 | none | done, PR #3083 |
| UX.0 | Bugs and small fixes (standing item) | none | open while it holds bugs |
| UX.1 | Class reference pages | 56 (2026-10-01 01:19Z to 04:37Z) | done in 7.46.0, PR #167 |
| UX.2 | Menus sized to what they hold (bug) | 62 (2026-10-01 01:44Z to 05:29Z) | done, PR #528 |
| UX.3 | Warn every visitor while a background job changes the galaxy | 87 (2026-10-01 03:26Z to 05:29Z) | done, PR #2086 |
| UX.4 | Phenomenon pages (new parent) | none | done (UX.17 and UX.18) |
| UX.5 | Place facilities from the web interface | 31 (2026-09-30 18:14Z); 36 (2026-09-30 18:39Z to 2026-10-01 04:37Z) | done in 7.47.0, PR #167 |
| UX.6 | Every distance in its most meaningful unit | 1 (2026-09-30 18:14Z to 22:15Z) | done in 7.10.0, PR #122 |
| UX.7 | Planet list: type chip, habitable-moon chip, belt distances | 2 (2026-09-30 18:14Z to 22:36Z) | done in 7.12.0, PR #126 |
| UX.8 | A less dense top bar (gear menu) | 3 (2026-09-30 18:14Z to 22:36Z) | done in 7.11.0, PR #126 |
| UX.9 | Tag search: collapsible groups and phenomena | 4 (2026-09-30 18:14Z to 2026-10-01 00:34Z) | done in 7.22.1, PR #135, and 7.29.0, PR #143 |
| UX.10 | Timestamps in the viewer's own time zone | 14 (2026-09-30 16:51Z to 18:09Z); 22 (2026-09-30 18:14Z to 23:48Z) | done in 7.21.0, PR #135 |
| UX.11 | Paginated list of every system | 59 (2026-09-30 20:01Z to 20:27Z); 46 (2026-09-30 20:08Z to 23:48Z) | done in 7.20.0, PR #135 |
| UX.12 | System page: one ordered list of everything in orbit | 61 (2026-09-30 20:01Z to 20:27Z); 48 (2026-09-30 20:08Z to 22:36Z) | done in 7.12.0, PR #126 |
| UX.13 | One meaningful-unit ladder for speeds | 102 (2026-10-01 05:27Z to 05:29Z) | done, PR #234 |
| UX.14 | One meaningful-unit ladder for time periods | 103 (2026-10-01 05:27Z to 05:29Z) | done, PR #234 |
| UX.15 | Put an object's data beside its 3D render when there's room (bug) | 105 (2026-10-01 05:43Z to 05:54Z) | done, PR #195 |
| UX.16 | Always leave space between buttons (bug) | none | done, PR #195 |
| UX.17 | A view that suits each phenomenon | 25 (2026-09-30 18:14Z to 2026-10-01 05:05Z) | done in 7.57.0, PR #178 |
| UX.18 | Show phenomena's octant and nearest systems | 26 (2026-09-30 18:14Z to 2026-10-01 04:37Z) | done: storage in 7.33.0 (PR #149), pages in 7.48.0 (PR #167) |
| UX.19 | Asteroid belt rows: density, range, top minerals; no Zone column (bug) | none | done, PR #213 |
| UX.20 | Scientific notation past 4 digits before the decimal point (bug) | none | done, PR #216 |
| UX.21 | Clean up the web interface: overlapping buttons and dead controls (bug) | none | done, PR #707 |
| UX.22 | Meaningful units for every measurement | none | done, PR #1504 |
| UX.23 | A shared unit-ladder module | none | done (shipped as src/planetgen/util/ladders.py; its entry left TODO.md earlier) |
| UX.24 | Sector contents: rogue planets after systems and phenomena, and expanded rows the full table width (bug) | none | done, PR #487 |
| UX.25 | Rogue planets: octant and a small map symbol beside each name (bug) | none | done, PR #490 |
| UX.26 | Edit and admin actions as a button that opens a menu (bug) | none | done, PR #558 |
| UX.27 | System page: the system and navigation buttons on one row that doesn't overlap (bug) | none | done, PR #539 |
| UX.28 | Investigate icons instead of words on buttons | none | done, PR #490 |
| UX.29 | Every comet in a system shows its type as a link (bug) | none | done, PR #487 |
| UX.30 | Planet information without the Markdown render | none | absorbed into UX.105 |
| UX.31 | Editing a star system: an edit button with a quick menu, not a long panel (bug) | none | done, PR #533 |
| UX.32 | Planet rows show the class only, without the type and moon labels | none | merged into UX.100 (Boss, 2026-10-16, 11:58 PM EDT: consolidation sweep) |
| UX.33 | Filter phenomena by their classes and types (bug) | none | done, PR #606 |
| UX.34 | The sector summary calls white dwarfs "B-type" and "A-type" systems (bug) | none | done, PR #448 |
| UX.35 | The NAV page route shown horizontally, wrapping onto several lines on narrow screens | none | done, PR #879 |
| UX.36 | Scientific notation starts too early for whole numbers (bug) | none | done, PR #457 |
| UX.37 | A UX sweep: remove redundant and duplicate controls so the interface gets out of the way | none | done, PR #671 |
| UX.38 | The nebula and remnant diagrams' "-" button does nothing at the 1 ly limit (bug) | none | done, PR #590 |
| UX.39 | Markdown rendered by the markdown library | none | done, PR #759 |
| UX.40 | Buttons, menus and dialogs from Shoelace web components | none | done, PR #544 |
| UX.41 | Tables on TanStack Table and TanStack Virtual | none | done, PR #606 |
| UX.42 | In-universe wording across the interface | none | done, PR #1466 |
| UX.43 | A visual design built like a pilot's starmap and navigation console | none | open |
| UX.44 | Search: mutually exclusive tags should combine with OR, the rest with AND (bug) | none | done, PR #457 |
| UX.45 | Bookmark management | none | open |
| UX.46 | Wiping the galaxy wipes the bookmarks | none | open |
| UX.47 | A bookmark manager: list, rename, sort, group and delete | none | open |
| UX.48 | Charted regions as bookmarks that frame and outline the region | none | open |
| UX.49 | Form fields as Shoelace components (sl-input, sl-select, sl-checkbox) across the site | none | done, PR #3387 |
| UX.50 | One short hint per map, the rest behind a Help entry | none | done, PR #640 |
| UX.51 | Cards that repeat the page title | none | done, PR #663 |
| UX.52 | Sector name repeated on every Contents row | none | done, PR #654 |
| UX.53 | The same star shown three times on a system page | none | done, PR #654 |
| UX.54 | Search shows empty result groups and echoes the query | none | done, PR #671 |
| UX.55 | Home and Systems repeat other pages’ tables | none | done, PR #671 |
| UX.56 | Admin hub repeats the gear menu | none | done, PR #632 |
| UX.57 | Two controls both called Reset | none | done, PR #637 |
| UX.58 | Move Current into the Steps menu | none | done, PR #637 |
| UX.59 | One breadcrumb trail on map pages | none | done, PR #640 |
| UX.60 | The Slabs rail has nothing in it at the top level | none | done, PR #637 |
| UX.61 | System Map: Measure distance floats beside an empty gap | none | done, PR #640 |
| UX.62 | Map pages jump 64 px left | none | done, PR #637 |
| UX.63 | One action bar on every object page | none | done, PR #618 |
| UX.64 | Action buttons are all solid primary, with no hierarchy | none | done, PR #620 |
| UX.65 | Facts and links mixed in the Sector header chips | none | done, PR #635 |
| UX.66 | “Cube edge” on arc-shaped sectors | none | done, PR #617 |
| UX.67 | Wikitext and Markdown buttons | none | done, PR #648 |
| UX.68 | System Admin menu lists every planet and moon | none | done, PR #646 |
| UX.69 | Tables are cut off on phones with no cue | none | done, PR #659 |
| UX.70 | The result page repeats the route and puts the map last | none | done, PR #648 |
| UX.71 | NAV landing page | none | done, PR #648 |
| UX.72 | Gear menu mixes unrelated things; the admin hub is not the hub | none | done, PR #632 |
| UX.73 | Sector wiki link form on the Admin hub | none | done, PR #632 |
| UX.74 | Generate page shows an empty Current job card | none | done, PR #617 |
| UX.75 | Admin menus on the sector and phenomenon pages sit in the page body, not the action bar | none | done, PR #679 |
| UX.76 | Icons and highlights follow the light and dark theme (bug) | none | done, PR #800 |
| UX.77 | The class list is alphabetized (bug) | none | done, PR #797 |
| UX.78 | Unit preference: Automatic, Metric only or Customary | none | done, PR #1562 |
| UX.79 | Python and JS round half-way values differently (bug) | none | done, PR #843 |
| UX.80 | Negative values that round to zero print "-0" (bug) | none | done, PR #843 |
| UX.81 | Time symbols Gyr, Myr, kyr in place of Gy, My, ky; AU from 1,000,000 km; scientific text below mantissa 1e-3 | none | done, PR #1504 |
| UX.82 | Theme checks after PR #800: SVG currentColor, two Shoelace contrast failures, alpha in --bg-subtle | none | done, PR #4708 |
| UX.83 | Generation steps that run long show no progress bar of their own: linking new sectors to their neighbours, the phenomenon scatter and others (bug) | none | done, PR #895 |
| UX.84 | Every sub-step must show a progress bar that starts by itself when it is predicted to take over 15 seconds (bug) | none | done, PR #930 |
| UX.85 | Button menus open out of sight and make the user scroll to see them (bug) | none | done, PR #932 |
| UX.86 | The Galaxy Map controls take too much room: buttons too large and filters one character wide (bug) | none | done, PR #936 |
| UX.87 | The system list shows uncharted systems: every scattered star, with its location and a way to generate it | none | done, PR #1097 |
| UX.88 | Say "uncharted" instead of "unbuilt" or "not generated" in all user-facing text outside Generate and admin | none | done, PR #1049 |
| UX.89 | Staged jobs show wrong stage counts and numbers, and skipped stages are not listed (bug) | none | done, PR #1057 |
| UX.90 | Explain the habitability chips in the web interface: a visible legend with the equipment labels, colours, thresholds and scores | none | done, PR #1065 |
| UX.91 | Planet and moon description carries a full PHI-4 explanation, each colour factor and why | none | done, PR #1104 |
| UX.92 | A Bookmark button for planets, moons and belts | UX.45 | open |
| UX.93 | No TODO code (like PERF.67 or NAV.42) appears anywhere a user can see it, with a test that fails if one does | none | done, PR #1378 |
| UX.94 | Gravity map wording: the legend, the mode names and the Map help text | none | open |
| UX.95 | Tables show their rows first and fill the filter-menu counts a moment later | none | done, PR #1299 |
| UX.96 | A table that hits the statement limit says the database is busy and retries, instead of failing with a 502 | none | done, PR #1267 |
| UX.97 | Search results: each panel runs under its own time limit and is fetched on its own | none | done, PR #1315 |
| UX.98 | The Generate log has odd wording and mixed section headers (bug) | none | done, PR #1369 |
| UX.99 | Rogue planets do not show their surface gravity (bug) | none | done, PR #1355 |
| UX.100 | Planet and moon lists show habitability as symbols, not notes | none | done, PR #2784 |
| UX.101 | The PHI-6 symbols as vector shapes in the icon library, both themes, coloured by rating | none | done, PR #2706 |
| UX.102 | One logging standard for every output the user or an admin reads: the same timestamps, format, dividers and step markers everywhere | none | done, PR #4450 |
| UX.103 | Design note for the description layer | none | done, PR #1496 |
| UX.104 | Build the document model for stars, pairs, planets, moons, belts, comets and the overview | none | done, PR #1509 |
| UX.105 | The system page builds every body's detail from the document model as HTML tables | none | done, PR #2206 |
| UX.106 | CLI output as plain text by default, Markdown and wikitext on request | none | done, PR #1533 |
| UX.107 | Units in the description go through the ladders | none | done, PR #2240 |
| UX.108 | Describe phenomena, sectors, galaxy and facilities in the same Document | none | done, PR #2206 |
| UX.109 | Remove the old per-class description text | none | done, PR #2234 |
| UX.110 | The phenomena count on the overview page is wrong (bug) | none | done, PR #2124 |
| UX.111 | Reads, API and web wording on shells | none | done, PR #1607 |
| UX.112 | Maintenance mode: the site is closed to everyone but admin while the maps are rebuilt | none | done, PR #1607 |
| UX.113 | The sector page shows the profile of any address in the reach, and heavy cells and far neighbourhoods warn before they run | none | done, PR #1919 |
| UX.114 | Scatter-row pages list the stored members and the map shows a system mass | none | done, PR #1989 |
| UX.115 | Interface update: every page and control works with the spherical grid, 96-bit IDs, the new phenomena and the region queries (stage after the sphere and generation stages) | none | open |
| UX.116 | Sliders for every scatter and fill limit: one limit_slider with a number box, in three new segments under the scatter-limits row of the Generate page | none | done, PR #1911 |
| UX.117 | The Apply filters button on the phenomena list disappears while the list loads (bug) | none | done, PR #1778 |
| UX.118 | Slow pages show a progress bar while their data loads: the search page and the phenomena list pages | none | done, PR #2089 |
| UX.119 | An HTML renderer for the Document: facts, rating and composition tables, with unit spans, in both themes | none | done, PR #2206 |
| UX.120 | The system page and the one-off system page use it: every body's detail from the Document | none | done, PR #2206 |
| UX.121 | The logging layer and the written standard | none | done, PR #3049 |
| UX.122 | The command line, the job log files and the downloadable log use the layer | none | done, PR #3600 |
| UX.123 | The Generate page log, the admin log pages and update.sh output use the layer | none | done, PR #3540 |
| UX.124 | One form-field look and error slot, piloted on the search, nearby and list-filter fields | none | done, PR #2716 |
| UX.125 | The Generate page and the one-off system page on the field macro | none | done, PR #3184 |
| UX.126 | The admin pages and the edit dialogs on the field macro | none | done, PR #3387 |
| UX.127 | The Galaxy Map block menu offers Generate again | none | done, PR #3567 |
| UX.128 | Map side panels and other browser call sites take Quantity | none | done, PR #2801 |
| UX.129 | Saved lists in the browser and in the account: one store with a local and an account backend | none | done, PR #2544 |
| UX.130 | Console design brief: colours, type, panel shapes and iconography, approved by Boss | none | open |
| UX.131 | Console tokens, type and base template in both themes | none | open |
| UX.132 | Console look for the shared components: buttons, tables, panels, status chips, icons and the maps' chrome | none | open |
| UX.133 | `GET /polities` and `GET /species` answer 404 before the population pass has run: show an empty page instead (bug) | none | done, PR #2124 |
| UX.134 | Generate page: the galaxy type shows what the seed drew and can be pinned | none | open |
| UX.135 | One simple vector outline per galaxy family in the icon library: spiral, barred spiral, lenticular, elliptical, irregular and dwarf spheroid | none | open |
| UX.136 | The data tables' own row fetches (scroll and sort) wait on a page job when the answer is slow, like the pages do | none | done, PR #2999 |
| UX.137 | The Generate page, the settings and the command line set the halo and the satellite galaxies switches, each through its own command, and say what they are | none | open |
| UX.138 | Pages, symbols and descriptions for the new rows: moons of wanderers, free brown dwarfs as stars, unbound comets | none | open |
| UX.139 | Display reads the registry: ladder rungs take their symbol and factor from it, and the 21 format wrappers become one format_quantity entry | none | done, PR #2490 |
| UX.140 | Browser table: the registry generated to JS beside unitladders.js, and the repeated browser constants removed | none | done, PR #2571 |
| UX.141 | List pages for every kind: systems, stars, planets, moons, belts, comets and each phenomenon type, and the name search, read the object index with keyset paging | none | done, PR #2408 |
| UX.142 | The viewing windows scale to the browser window so the controls are always in view without scrolling the page | none | open |
| UX.143 | A dedicated page for every kind of item: star, planet, moon, belt, comet, each phenomenon, nebula, cluster, system, sector and region, reached from its Object ID | none | open |
| UX.144 | Every item page says where the item is: links up the chain (galaxy, sector, system, parent star or planet) and the nebula or other region that contains it | none | open |
| UX.145 | The admin menu opens past the edge of the visible screen on the Sector Map (and possibly on other pages) (bug) | none | done, PR #2416 |
| UX.146 | The top menu bar: Lists moves into the gear menu, the header search becomes an icon-only button, Nearby goes under Nav as Plot Course then Nearby Search (bug) | none | done, PR #2610 |
| UX.147 | The Generate page has no clear way to generate only a galaxy plan, with no scatter, to see the sector layout on the galaxy screen and add sectors by hand (bug) | none | done, PR #2729 |
| UX.148 | Log output from Generate is no longer in a fixed-width font (bug) | none | done, PR #2838 |
| UX.149 | Chemical names and formulas in the HTML output use subscripts and superscripts, as water (H2O) with a subscript 2 (bug) | none | done, PR #2767 |
| UX.150 | Pages, lists, search, filters, the API and the map markers show and filter by the derived class and the PHI-6 colour; the class wording and the descriptions follow the rule table | none | open |
| UX.151 | Icons for Microbial, Flora, Fauna and Sentient life, with a 0 to 9 digit, drawn as a small cluster next to the PHI-6 icons, both themes | none | done, PR #2733 |
| UX.152 | Planet and moon pages, lists, the system page and the API show the four life scores in the cluster; the old Microbial, Complex life and Human operability line goes | none | done, PR #3038 |
| UX.153 | Search offers a Needs tick per equipment item, the habitability levels help page and the describer sentences read the needs list, and no sentence says none | none | done, PR #3038 |
| UX.154 | Lists and the system page show the equipment lists as chips: solid for required, outlined for recommended | none | done, PR #2844 |
| UX.155 | The /lists planet and moon pages show the PHI-6 symbols in a Habitability column and the API fields carry them (UX.100 part 2) | none | done, PR #2784 |
| UX.156 | The planet and moon item pages show the PHI-6 symbols (UX.100 part 3) | none | done, PR #2784 |
| UX.157 | Every sprite icon draws as a solid black shape because the site's Content-Security-Policy refuses the sprite's own style element (bug) | none | done, PR #2890 |
| UX.158 | Info panel line: how this colour was worked out, which inputs and whether partial or full quality | none | done, PR #4511 |
| UX.159 | Remove the old /<name>.py CGI redirects (web/old_urls.py) and the old NAV from_id and to_id parameter translation, with their tests | none | done, PR #2873 |
| UX.160 | The search icon button opens a box over the page: it should go straight to the search page (bug) | none | done, PR #2884 |
| UX.161 | After a galaxy reset the counts on the main page do not change even though the database has been wiped (bug) | none | done, PR #2884 |
| UX.162 | The plan prints a "Where the plan's time went" table in every galaxy plan; it should show only in debug mode or with a --timing flag (bug) | none | done, PR #2877 |
| UX.163 | Every planet and moon page and the descriptions carry a Required line and a Recommended line for the equipment, on every run | none | done, PR #2883 |
| UX.164 | The wipe line reads "Wiped 75 tables (about 1 row)" (bug) | none | open |
| UX.165 | Class pages and names come from the rule rows | none | done, PR #2781 |
| UX.166 | The habitability icons display as intended in Chrome but draw as solid black shapes with tiny coloured dots in Firefox, even with a fully cleared cache (bug) | none | done, PR #2890 |
| UX.167 | Rewrite the equipment category row sentences from the needs list instead of the 10-15 wording | none | done, PR #2902 |
| UX.168 | Property sheet registry and GET /api/objects/<ref>/sheet, with a completeness test against column_units.py | none | done, PR #3003 |
| UX.169 | Children API GET /api/objects/<ref>/children and the sector-level tree feed | none | done, PR #3003 |
| UX.170 | Rebuild item.html: header, Location chain with the Inside remainder, siblings, property sheet, description and children | none | done, PR #3003 |
| UX.171 | bodyview.js: the page rendered example from data-look, with the colour note and a browser test | none | done, PR #3003 |
| UX.172 | Tree component: partial, script, keyboard and screen-reader behaviour, no-script fallback and narrow-screen layout | none | done, PR #3003 |
| UX.173 | System page: stars table and body list become the tree; remove the dropdown rows and move species, facilities and the admin menu to the object page | none | done, PR #3011 |
| UX.174 | Sector page: Contents rows expand into the tree, keeping facets, sort, jump and paging at the top level | none | done, PR #3011 |
| UX.175 | The same page anatomy for phenomena, facilities, clusters and regions (absorbs the UX.143 remainder) | none | done, PR #3065 |
| UX.176 | Link sweep: object_url, search, lists, describer references and the map popup View button all open /object/<ref>; delete body_anchor, bodyanchor.js and the row-highlight styles | none | done, PR #3011 |
| UX.177 | Tests for the object pages: a page-walk for every kind, tree browser tests, accessibility, no-script, a huge sector, a moon of a moon-less planet and an unknown ID | none | done, PR #3011 |
| UX.178 | Docs for the object pages: docs/design/interface/object-pages.md, the route list, body-rendering.md and the sheet registry rules | none | done, PR #3011 |
| UX.179 | Sitemap entries and an Open Graph image for object pages (optional, can wait) | none | open |
| UX.180 | A Help menu in the top bar on the search-icon side that holds Classes and Habitability; the separate Habitability levels links go away | none | done, PR #2867 |
| UX.181 | Every habitability symbol, life score and equipment or Habitable chip links straight to its section of the habitability page | none | done, PR #2867 |
| UX.182 | The standalone Planet class, Spectral class and Luminosity class links go away and the class text itself links to its section | none | done, PR #2867 |
| UX.183 | Facility pages: add the facility kind to the object index and resolver so /object/facility:<ID> works | none | done, PR #3060 |
| UX.184 | Star cluster and region pages: stored cluster and region objects and their object pages | none | done, PR #3065 |
| UX.185 | Phenomenon pages: paragraph text and neighbour links from the stored row | none | done, PR #3045 |
| UX.186 | Pages show the X-ray luminosity, the microquasar flag and the spider millisecond pulsar kind of an accreting compact binary | none | done, PR #3775 |
| UX.187 | Pages, lists, the API and the help text show every physics rectification output | none | open |
| UX.188 | Captain-first object pages, step C0: the GLANCE registry, the At a glance panel partial and its test (every key is in the sheet; the page renders the panel above the sheet) | none | done, PR #3116 |
| UX.189 | Captain-first object pages, step C1: planet and moon pages | none | done, PR #3116 |
| UX.190 | Captain-first object pages, step C2: star and star system pages | none | done, PR #3134 |
| UX.191 | Captain-first object pages, step C3: sector pages | none | done, PR #3134 |
| UX.192 | Captain-first object pages, step C4: phenomenon pages | none | done, PR #3134 |
| UX.193 | Captain-first object pages, step C5: belt, comet, facility, cluster and superbubble pages | none | done, PR #3134 |
| UX.194 | Captain-first object pages, step C6: galaxy page | none | done, PR #3134 |
| UX.195 | Every page has a proper heading outline: one h1 for the title, then h2 and below for each section and table, with no skipped level (bug) | none | done, PR #3196 |
| UX.196 | Every page gets an automatic table of contents built from its headings (bug) | none | done, PR #3196 |
| UX.197 | The habitability icons link to the habitability section of the planet information page, not to the help page (bug) | none | done, PR #3150 |
| UX.198 | The Atmosphere section of the planet property sheet shows 14 empty rows (bug) | none | done, PR #3196 |
| UX.199 | Pages show raw or placeholder values: "not stored" rows in At a glance, and a raw key such as Technological_civilization (bug) | none | done, PR #3196 |
| UX.200 | The Location block lists a system and its star separately when they share a name, shows "not inside any nebula" style lines, and links the Galaxy (bug) | none | done, PR #3188 |
| UX.201 | A Nav button on object pages opens one menu: navigate from here, navigate to here, nearby locations (bug) | none | done, PR #3188 |
| UX.202 | Remove the "How the composition was found" line from the object pages (bug) | none | done, PR #3188 |
| UX.203 | Pages state a missing thing in one sentence instead of a table of zeros: no life, no atmosphere (0 Pa), effectively no water | none | done, PR #3188 |
| UX.204 | Every collapsible section on an object page starts collapsed | none | done, PR #3188 |
| UX.205 | Object pages show the day and night time shares, moons with the extra night from their planet | none | done, PR #3755 |
| UX.206 | Link-preview (SEO and social) images for pages with a rendered image or map view are a cached snapshot of that view | none | open |
| UX.207 | Breadcrumb bar: the system link is labelled "System" and the star link carries the star(s) name (bug) | none | done, PR #3188 |
| UX.208 | Fill audit of the property sheet: sample bodies of every class, every sheet column filled or not applicable by a rule, no "not stored" for an applicable column | none | done, PR #3629 |
| UX.209 | Object pages order: At a glance, then About, then the Properties list, then the other details (bug) | none | done, PR #3196 |
| UX.210 | Map infoboxes and labels: the orbit infobox eccentricity and every other number on a map says what it is and its unit (scan) (bug) | none | done, PR #3347 |
| UX.211 | Pages, lists, search, Nearby and the admin pages: every number carries a unit, a noun or a label (scan) (bug) | none | done, PR #3210 |
| UX.212 | Admin pages show bare numbers: the statistics and queue tiles (stat-value 0, 18) and the admin settings values carry no unit or noun (bug) | none | done, PR #3423 |
| UX.213 | Generate page task option values show as bare numbers (1, 2, 4, 6): name what they count (bug) | none | done, PR #3725 |
| UX.214 | Search, object lists, Nearby and the API filter and sort by Suitability, the four life scores (Flora, Fauna, Sentient, food) and PHI-6: today they only draw the PHI-6 symbol column | none | done, PR #3518 |
| UX.215 | Object pages: the table of contents is a menu on the left side, one entry for each section heading, and the location links sit in a Location row of the at-a-glance table | none | done, PR #3573 |
| UX.216 | Object pages: the "Rendered example" panel gets an in-universe navigation title | none | done, PR #3573 |
| UX.217 | System and planet pages name the zones: "Ecosphere / Habitable Zone", "Hot Zone", "Cold Zone", each with its range in appropriate units, instead of "Zone E" | none | done, PR #3573 |
| UX.218 | The class letter on pages is a link to the page for that class (Class M and the rest) | none | done, PR #3573 |
| UX.219 | Planet pages: Atmosphere and Composition are each a heading followed by their table; the composition table is the planet's actual composition by mass | none | done, PR #3573 |
| UX.220 | Planet pages: a short summary paragraph generated from the data replaces "It is a temperate world"; it holds no details, class, surface conditions, safety equipment, PHI-6 summary, sustainability or life scores | none | done, PR #3573 |
| UX.221 | Planet pages: a Life paragraph (only when the world has life) and a Civilization / Political paragraph (only when it has a civilization), in paragraph form from the life scores and the polities | none | done, PR #3573 |
| UX.222 | Planet pages: the section order is About planet, Table, Summary, Life, Civ / Political, PHI-6 table and details with Radiation and Survival equipment, Atmosphere, Composition, each a heading in the TOC menu | none | done, PR #3573 |
| UX.223 | Resource pages: belts, fields, giants and icy bodies show their resources, and lists, search and the API filter and sort by resource and score | none | open |
| UX.224 | A planet's chips and description never show Habitable, Inhabited and Habitable class together: only the most specific one (bug) | none | done, PR #4452 |
| UX.225 | Search, object lists, Nearby and /api/search filter and sort on every stored physical number: gravity, temperature, pressure, orbit, day, star mass, dose, field, ocean depth, moon count and gases (G4) | none | open |
| UX.226 | Numeric range filters and sort keys for the stored physical values: one engine in db/query.py and the /api/search parameters (G4, step 1) | none | done, PR #4523 |
| UX.227 | The search page shows range controls for the physical values, with units, and the active filters as chips (G4, step 2) | none | done, PR #4552 |
| UX.228 | The object lists and Nearby filter and sort by the same physical values (G4, step 3) | none | open |
| UX.229 | Colour contrast of the habitability life scores and tree-leaf stats passes 4.5 to 1 in both themes (bug) | none | done, PR #3280 |
| UX.230 | The Generate system page options have a valid listbox parent so axe aria-required-parent passes (bug) | none | done, PR #3280 |
| UX.231 | Sector and object pages say "in Vang Zarh Region" in the at-a-glance block, hidden when the position is in no charted region (regions chain, step 5) | none | done, PR #3616 |
| UX.232 | Search results: the planet link in a result row goes to a bad link while the star system link in the same row works (bug) | none | done, PR #3518 |
| UX.233 | The Map symbols key on the object pages | none | done, PR #3611 |
| UX.234 | Searching for all M class planets with the search button says it is preparing the results, loops and never shows data (bug) | none | done, PR #3498 |
| UX.235 | Pages that already have a table of contents also show the automatic contents bar, so the contents appear twice (bug) | none | done, PR #3515 |
| UX.236 | User preferences: a screen and storage for per-user settings, starting with the Galaxy Map region landing (frame or densest) | none | done, PR #4448 |
| UX.237 | Planet page Life section: Supports Class N, also supports, Has Class N at the X stage with its animal grade, the eats-drinks-breathes-gives-off block, the Life and air chart and the adaptation traits; the class number on lists | none | open |
| UX.238 | The "On this page" table of contents on object pages goes under the tab bar, not off to the left (bug) | none | done, PR #3900 |
| UX.239 | Search for star systems by how many stars they hold: single, binary, trinary, quaternary and more | none | done, PR #4569 |
| UX.240 | Belt page data builder: the designated individual bodies (designation, mass, size, status, composition), the size table of the rest of the belt mass (counts, mass share and typical composition by diameter bin), and the density factor, as one pure structure | none | done, PR #4259 |
| UX.241 | Asteroid belt information page: route, template and API field showing the individually tracked bodies (by designation, each linked to its member page), the contents and sizes of the rest of the belt mass, and the density factor with its unit; link from the system page and the belt row | none | open |
| UX.242 | Systems tables: Stars, Objects and Habitable world columns replace the Binary column (All systems table and the search Systems panel) | none | done, PR #4216 |
| UX.243 | Day-and-night counter on every object page: a clock counting the object's own day at its real rate, with the state word, shown separately from the spin animation speed (spin.js compresses periods onto 20 s to 10 min per turn, stand-in SP1, L4) | none | open |
| UX.244 | Slow and tidally locked bodies look frozen on screen: spin.js clamps every period over 1000 h to one turn per 600 s (stand-in SP1, L4), so Eakkar IV (class V, rotation 46.8 days, locked to its orbit) never visibly turns (bug) | none | done, PR #4254 |
| UX.245 | Tidally locked moon faces its planet in the System Map, and the locked label is drawn on the sensor view picture itself (follow-up to UX.244) (bug) | none | done, PR #4266 |
| UX.246 | Search page Star: Spectral Class and Star: Luminosity Class tags list only the types already stored, so a small database offered 2 types; list every class the generator produces with a count of 0 for a type no star holds yet (bug) | none | done, PR #4261 |
| UX.247 | The "On this page" menu on object pages is a single sideways-scrolling row (UX.238); restore the grouped list of sections and subsections (left column on wide windows, list under the title on narrow ones) (bug) | none | done, PR #4263 |
| UX.248 | Abralinar (triple star: close A/B pair plus distant C): Abralinar I is drawn circling the system barycenter instead of the A/B pair it orbits, and nothing says the system is a close pair orbited by a third star (bug) | none | done, PR #4268 |
| UX.249 | The "On this page" menu should sit in one row across the top when every link fits without scrolling, and be the grouped list otherwise (follow-up to UX.247) (bug) | none | done, PR #4274 |
| UX.250 | Planet render: the terminator between night and day should be more of a crescent, showing more of the day side, and the "to star" arrow overlaps the render (bug) | none | done, PR #4462 |
| UX.251 | Planet and moon pages: a Life summary paragraph in At a glance, from the life readings | none | done, PR #4488 |
| UX.252 | Readability rewrite of ALL planet, moon and belt text: summary sentences, descriptions, notes, labels, class pages | none | open |
| UX.253 | Readability rewrite of ALL star, system and star-class text: summary sentences, descriptions, notes, labels, spectral class pages | none | done, PR #4496 |
| UX.254 | Show a star's chemistry on its page: thin or thick disc branch, alpha enhancement, metal fraction and the dredge-up shift, from stored values with stellar_composition.star_makeup | none | done, PR #4540 |
| UX.255 | Remove the light-theme branches in the map scripts that black map backgrounds left unused | none | open |
| UX.256 | Cancelling a job from its job page returns to the page the visitor came from, or the jobs page (bug) | none | done, PR #4656 |
| VIEW.1 | View from a planet | none | open |
| VIEW.2 | A starmap seen from a planet. RESEARCH WITH BOSS FIRST | 83 (2026-10-01 02:55Z to 05:29Z) | open |
| VIEW.3 | Render the view as a PNG, with constellations | 84 (2026-10-01 02:55Z to 05:29Z) | open |
| VIEW.4 | Constellation names in the name generator | 85 (2026-10-01 02:55Z to 05:29Z) | open |
| VIEW.5 | Light-travel positions: where an object appears to a distant observer | none | done, PR #719 |
| VIEW.6 | Settle the handedness of the generated galaxy before mapping real sky coordinates | none | open |
| VIEW.7 | Sky tiers and the backfill cost measurement | none | open |
| VIEW.8 | A dust model for the sky: the galactic extinction law and `nebulae.extinction_av` | none | open |
| VIEW.9 | Replace the BC_V table with a published relation (Flower 1996 or Torres 2010) | none | done, PR #4508 |
| VIEW.10 | A HYG calibration test: a generated sky has about 9,000 stars at V <= 6.5 | none | merged into VIEW.12 (Boss, 2026-10-16, 11:58 PM EDT: consolidation sweep) |
| VIEW.11 | Photometry library: apparent magnitude, bolometric correction, blackbody colour and the extinction law | none | open |
| VIEW.12 | Sky catalogue: the stars an observer sees, with direction, brightness and light-time | none | open |
| VIEW.13 | Galactic band: the glow of the unresolved stars and dust along each line of sight | none | open |
| VIEW.14 | Orientation and horizon: the planet's own sky, time of day and twilight | none | open |
| VIEW.15 | Whole-sky PNG: Hammer-Aitoff projection, written with numpy and zlib, cached and served by the API | none | open |
| VIEW.16 | Constellations: the seeded partition over the brightest stars, stored per system, named by the codec, with Bayer labels | none | open |
| VIEW.17 | Constellation overlay: SVG lines and names over the sky chart, and the planet-sky card | none | open |
| VIEW.18 | Headline of 9.1: the view from a planet, with its starmap, its picture and its constellations | none | open |

## Tree IDs to flat IDs

The dotted IDs used from PR #189 (about 06:05Z on 2026-10-01) to the
flat-ID PR, and what they are now. Top-level tree IDs (`MAP.2`, `GEN.4`)
kept their numbers; every subitem got the next number in its category,
in tree order. The tree had a standing `UX.0` "Bugs and small fixes"
parent; it is gone, and its bugs are top-level items (UX.15, UX.16).

| Tree ID | Flat ID |
|---|---|
| ADM.1.1 | ADM.5 |
| ADM.1.2 | ADM.6 |
| ADM.1.3 | ADM.7 |
| ADM.1.4 | ADM.8 |
| GEN.4.1 | GEN.10 |
| GEN.4.2 | GEN.11 |
| GEN.4.3 | GEN.12 |
| GEN.4.4 | GEN.13 |
| GEN.4.5 | GEN.14 |
| GEN.5.1 | GEN.15 |
| GEN.5.2 | GEN.16 |
| GEN.5.3 | GEN.17 |
| GEN.5.4 | GEN.18 |
| GEN.5.5 | GEN.19 |
| GEN.5.6 | GEN.20 |
| GEN.7.1 | GEN.21 |
| GEN.7.2 | GEN.22 |
| MAP.2.1 | MAP.16 |
| MAP.2.1.1 | MAP.17 |
| MAP.2.1.2 | MAP.18 |
| MAP.2.1.3 | MAP.19 |
| MAP.2.2 | MAP.20 |
| MAP.2.3 | MAP.21 |
| MAP.2.4 | MAP.22 |
| MAP.2.5 | MAP.23 |
| MAP.2.6 | MAP.24 |
| MAP.2.7 | MAP.25 |
| MAP.2.7.1 | MAP.26 |
| MAP.2.8 | MAP.27 |
| MAP.2.9 | MAP.28 |
| MAP.2.10 | MAP.29 |
| MAP.3.1 | MAP.30 |
| MAP.5.1 | MAP.31 |
| MAP.5.2 | MAP.32 |
| MAP.5.3 | MAP.33 |
| MAP.5.4 | MAP.34 |
| MAP.5.5 | MAP.35 |
| MAP.5.6 | MAP.36 |
| MAP.5.6.1 | MAP.37 |
| MAP.5.7 | MAP.38 |
| MAP.5.8 | MAP.39 |
| MAP.5.9 | MAP.40 |
| MAP.5.10 | MAP.41 |
| MAP.5.11 | MAP.42 |
| MAP.5.11.1 | MAP.43 |
| MAP.5.11.2 | MAP.44 |
| MAP.11.1 | MAP.45 |
| MAP.11.2 | MAP.46 |
| MAP.14.1 | MAP.47 |
| MAP.14.2 | MAP.48 |
| PERF.1.1 | PERF.3 |
| PERF.1.2 | PERF.4 |
| PERF.1.3 | PERF.5 |
| PERF.1.4 | PERF.6 |
| PERF.1.5 | PERF.7 |
| PERF.1.6 | PERF.8 |
| PERF.1.7 | PERF.9 |
| PERF.1.8 | PERF.10 |
| PERF.1.9 | PERF.11 |
| SEC.2.1 | SEC.3 |
| SEC.2.2 | SEC.4 |
| SEC.2.3 | SEC.5 |
| SEC.2.4 | SEC.6 |
| SEC.2.5 | SEC.7 |
| SEC.2.6 | SEC.8 |
| SEC.2.7 | SEC.9 |
| SEC.2.8 | SEC.10 |
| SEC.2.9 | SEC.11 |
| SEC.2.10 | SEC.12 |
| SEC.2.11 | SEC.13 |
| SEC.2.12 | SEC.14 |
| SEC.2.13 | SEC.15 |
| SEC.2.14 | SEC.16 |
| SEC.2.14.1 | SEC.17 |
| SEC.2.14.2 | SEC.18 |
| SEC.2.14.3 | SEC.19 |
| TEST.1 | Test category and suite markers | none | done, PR #251 (db/slow/browser markers; test_todo_tags reads bump_version.TODO_CATEGORIES) |
| TEST.2 | Parallel test runs | none | done, PRs #246 and #248 (pytest-xdist, per-worker control DB; the template schema was dropped, about 0.35 s per DB test) |
| TEST.3 | MariaDB in CI | none | done, PR #251 (CI legs MySQL 8.0, MySQL 8.4, MariaDB 11.4; 10.11 covered by local runs) |
| TEST.4 | Revive and widen the known-bug tests | none | done, PR #303 |
| TEST.5 | Real 4 pc in boundary tests | none | done, PR #303 |
| TEST.6 | SQL portability lint | none | done, PR #315 |
| TEST.7 | Strict sql_mode on both engines | none | done, PR #315 |
| TEST.8 | Migrate from real old schemas | none | done, PR #315 |
| TEST.9 | Migration crash and re-run | none | done, PR #315 |
| TEST.10 | Database newer than the code | none | done, PR #288 |
| TEST.11 | Every column round-trips | none | done, PR #315 |
| TEST.12 | Boundary values round-trip | none | done, PR #315 |
| TEST.13 | Collation collisions | none | done, PR #315 |
| TEST.14 | CHECK constraints enforced | none | done, PR #315 |
| TEST.15 | Sector save fails halfway | none | done, PR #315 |
| TEST.16 | Id blocks after reset and rollback | none | done, PR #315 |
| TEST.17 | Batched writes at the limits | none | done, PR #315 |
| TEST.18 | Full-text search edge cases | none | done, PR #315 |
| TEST.19 | Same galaxy at any worker count | none | done, PR #321 |
| TEST.20 | Work queue failure paths | none | done, PR #305 |
| TEST.21 | Cancelling a run | none | done, PR #305 |
| TEST.22 | Every bulk mode in parallel | none | done, PR #321 |
| TEST.23 | Resume after an interrupted fill | none | done, PR #303 |
| TEST.24 | Bright-star scatter edge cases | none | done, PR #303 |
| TEST.25 | Interrupted bright-star scatter | none | done, PR #303 |
| TEST.26 | `--force` scatter then fill | none | done, PR #303 |
| TEST.27 | Progress and ETA under bad clocks | none | done, PR #321 |
| TEST.28 | CLI errors by message | none | done, PR #303 |
| TEST.29 | Limits stay consistent | none | done, PR #303 |
| TEST.30 | Grid seams and the nucleus | none | done, PR #303 |
| TEST.31 | Sector placement exhaustion | none | done, PR #303 |
| TEST.32 | System builder internals | none | done, PR #303 |
| TEST.33 | Moon stability helpers | none | done, PR #303 |
| TEST.34 | Kepler solver extremes | none | done, PR #303 |
| TEST.35 | Star and evolution helpers | none | done, PR #303 |
| TEST.36 | Phenomenon class helpers | none | done, PR #303 |
| TEST.37 | Names under parallel saves | none | done, PR #321 |
| TEST.38 | Population incremental rescans | none | done, PR #321 |
| TEST.39 | Navigation graph | none | done, PR #321 |
| TEST.40 | Two admins start a job at once | none | done, PR #302 |
| TEST.41 | Job files damaged | none | done, PR #302 |
| TEST.42 | Pages fresh after a CLI write | none | done, PR #293 |
| TEST.43 | Auth sweep over every route | none | done, PR #293 |
| TEST.44 | What an API key may do | none | done, PR #293 |
| TEST.45 | More than one admin | none | done, PR #293 |
| TEST.46 | Trusted device and TOTP edge cases | none | done, PR #293 |
| TEST.47 | Oversized requests | none | done, PR #293 |
| TEST.48 | Security headers everywhere | none | done, PR #293 |
| TEST.49 | Thin API routes | none | done, PR #293 |
| TEST.50 | Galaxy URLs combined | none | done, PR #293 |
| TEST.51 | Page-number sweep gaps | none | done, PR #293 |
| TEST.52 | Old URLs and error codes | none | done, PR #293 |
| TEST.53 | Formatters with bad numbers | none | done, PR #293 |
| TEST.54 | Caches under threads | none | done, PR #293 |
| TEST.55 | Map buttons do something | none | done, PR #307 |
| TEST.56 | No overlapping controls | none | done, PR #307 |
| TEST.57 | Galaxy Map JavaScript logic | none | done, PR #307 |
| TEST.58 | Other map JavaScript | none | done, PR #307 |
| TEST.59 | Galaxy Map drill-down in a browser | none | done, PR #307 |
| TEST.60 | Admin script command lines | none | done, PR #288 |
| TEST.61 | SQLite import script | none | done, PR #288 |
| TEST.62 | update.sh against a real database | none | done, PR #288 |
| TEST.63 | Math check that runs first (new parent) | none | done, PRs #281 and #290 |
| TEST.64 | Reference values | none | done, PR #281 |
| TEST.65 | Identities and invariants | none | done, PR #281 |
| TEST.66 | Distributions match their targets | none | done, PR #281 |
| TEST.67 | Runs first in the suite and in CI | none | done, PR #281 |
| TEST.68 | Gate before bulk generation | none | done, PR #290 |
| TEST.69 | Intermittent failure in the colony test (bug) | none | done, PR #303 |
| TEST.70 | Tests for the map JavaScript | none | done, PR #351 |
| TEST.71 | Intermittent failure in the admin planet-regenerate test (bug) | none | done, PR #484 |
| TEST.72 | Intermittent failure in the two-step (2FA) sign-in test (bug) | none | done, PR #612 |
| TEST.73 | Intermittent failure in the parallel galaxy-run interrupt test (bug) | none | done, PR #371 |
| TEST.74 | Generation tests at more than one worker | none | done, PR #371 |
| TEST.75 | Tests for forcing and prevalence | none | done, PR #502 |
| TEST.76 | A bright-star test breaks on Python 3.9 and 3.10 (bug) | none | done, PR #371 |
| TEST.77 | A golden-seed regression test | none | done, PR #3152 |
| TEST.78 | A resume test's sector query fails under ONLY_FULL_GROUP_BY on MariaDB 10.11 (bug) | none | done, PR #484 |
| TEST.79 | Route edge cases, written before NAV.12 | none | done, PR #427 |
| TEST.80 | Intermittent failure in the admin change-star test (bug) | none | done, PR #442 |
| TEST.81 | Two processes reserving id blocks of one table can deadlock (bug) | none | done, PR #442 |
| TEST.82 | Intermittent failure in the orbit-ceiling trim test (bug) | none | done, PR #484 |
| TEST.83 | Rate-limit tests fail under parallel load (bug) | none | done, PR #643 |
| TEST.84 | The every-column round-trip test depends on whether a quasar got placed (bug) | none | done, PR #484 |
| TEST.85 | Name collisions can count -1 existing names and fail generation (bug) | none | done, PR #403 |
| TEST.86 | Intermittent failure in the concurrent-insert recovery test (bug) | none | done, PR #484 |
| TEST.87 | The two-process id-block test times out under full parallel load (bug) | none | done, PR #442 |
| TEST.88 | The facilities test fails when the drawn gas giant's sphere of influence is too small (bug) | none | done, PR #484 |
| TEST.89 | The Galaxy Map drill-down browser test fails intermittently (bug) | none | done, PR #484 |
| TEST.90 | test_old_jobs_are_pruned raises JobBusy under parallel tests (bug) | none | done, PR #484 |
| TEST.91 | The System Map drill-and-measure browser test can't click the first moon (bug) | none | done, PR #493 |
| TEST.92 | The web job runner's first-failure test reports the job as interrupted under full-suite load (bug) | none | done, PR #522 |
| TEST.93 | Timing tests fail and MariaDB drops connections under full-suite load (bug) | none | done, PR #522 |
| TEST.94 | test_old_jobs_are_pruned raises JobBusy again: the job lock outlives the finished job (bug) | none | done, PR #549 |
| TEST.95 | test_stabilize_lunar_system_respaces_crowded_moons keeps only one moon in the full suite (bug) | none | done, PR #549 |
| TEST.96 | test_galaxy_map_camera_presets_at_each_zoom_step fails in the full suite (bug) | none | done, PR #549 |
| TEST.97 | test_class_change_regenerates_surface_conditions raises StopIteration in the full suite (bug) | none | done, PR #549 |
| TEST.98 | test_scale_line_follows_the_zoom fails in the full suite (bug) | none | done, PR #549 |
| TEST.99 | test_slab_buttons_have_lines_that_follow_the_view fails in the full suite (bug) | none | done, PR #549 |
| TEST.100 | test_hover_while_picking_a_slab fails in the full suite (bug) | none | done, PR #549 |
| TEST.101 | test_interrupting_a_parallel_galaxy_run[2-True] fails in the full suite (bug) | none | done, PR #549 |
| TEST.102 | test_a_pass_removes_species_stored_without_a_civilization fails now and then in the full suite (bug) | none | done, PR #610 |
| TEST.103 | test_controls_do_not_overlap fails on clean main: controls overlap on the sector page at 600 and 820 px (bug) | none | done, PR #625 |
| TEST.104 | test_system_map_selection_drill_and_measure fails now and then under load (bug) | none | done, PR #634 |
| TEST.105 | test_web_a11y.py fails on the Admin hub: its links have 2.49:1 contrast (bug) | none | done, PR #648 |
| TEST.106 | test_controls_do_not_overlap fails on the Sector Map at 390, 600 and 1280 px: its buttons overlap the Contents filters (bug) | none | done, PR #650 |
| TEST.107 | test_controls_do_not_overlap fails on main for /galaxy at 390 px and /sector: buttons overlap each other (bug) | none | closed, not reproducible (passes on main) |
| TEST.108 | test_a_planet_holds_one_whose_system_frame_is_its_offset_from_the_star fails about one run in four (bug) | none | done, PR #690 |
| TEST.109 | test_sector_map_click_on_the_selected_nebula_clears_it fails on main since PR #605 (bug) | none | done, PR #690 |
| TEST.110 | Object ID tests: identical IDs on 1 and 4 workers, none reused, none missing | none | done, PR #1075 |
| TEST.111 | test_ensure_sector_generated_creates_then_reuses_the_same_sector fails in a busy parallel run (bug) | none | open |
| TEST.112 | test_regenerate_phenomenon_keeps_id_name_and_place fails now and then in the full suite (bug) | none | done, PR #955 |
| TEST.113 | tests/js/generatejobs.test.mjs fails on main: "asks for the job's status two seconds in" (bug) | none | done, PR #920 |
| TEST.114 | test_the_check_writes_nothing fails now and then in a parallel full run (bug) | none | done, PR #946 |
| TEST.115 | Three tests fail on the Windows CI leg in every recent run: Redis in WSL is unreachable (bug) | none | done, PR #981 |
| TEST.116 | test_bughunt_end_to_end stores M2V or M6V where it expects K2V on the MySQL 8.4 and MariaDB legs once (bug) | none | open |
| TEST.117 | generatejobs.test.mjs fails on main since PERF.33 (PR #910) (bug) | none | done, PR #920 |
| TEST.118 | test_star_scatter_passes.py fails twice on main since GEN.184 raised the luminosity floor to 2500 or more (bug) | none | done, PR #975 |
| TEST.119 | 15 browser map tests fail on main since the UX.86 menu regrouping (bug) | none | done, PR #993 |
| TEST.120 | test_sampled_stars_stay_inside_their_mass_range still fails on main after TEST.118 (bug) | none | done, PR #975 |
| TEST.121 | test_open_map_menus_hold_no_overlap fails under load: the Menu panel intercepts the close click (bug) | none | done, PR #1002 |
| TEST.122 | Browser map tests fail on plain main in Bugfixes lane 2's container (fixture maps, controls, system-page maps) (bug) | none | done, PR #1125 |
| TEST.123 | test_a_loaded_sector_knows_every_objects_cell_and_velocity fails once under full-suite load (bug) | none | done, PR #1121 |
| TEST.124 | Phenomena page tests fail with KeyError scattered in system_pages._phenomena_load (bug) | none | done, PR #1077 |
| TEST.125 | Three browser-map tests fail on main with "no generated system has a moon" (bug) | none | done, PR #1313 |
| TEST.126 | The heavy WebGL browser-map tests time out when four workers run them together (bug) | none | open |
| TEST.127 | test_nebula_shape_endpoint_serves_a_mesh failed once in a full parallel run (bug) | none | open |
| TEST.128 | test_spatial_position_db and test_web_db_fields fail under a parallel run on one MySQL (bug) | none | open |
| TEST.129 | Two orbit-update tests fail on the MySQL 8.4 leg: the stored last_updated_at rounds up (bug) | none | done, PR #1313 |
| TEST.130 | test_the_scene_positions_at_the_epoch_match_the_stored_ones fails when a random system has comets (bug) | none | done, PR #3017 |
| TEST.131 | Tests for the gravity field evaluator: two bodies, Lagrange points, softening and the rotation curve | none | merged into GEN.198 (Boss, 2026-10-16, 11:58 PM EDT: consolidation sweep) |
| TEST.132 | test_a_fill_after_a_layer_failed_mid_scatter_builds_no_leftover_star failed once in a full parallel run (bug) | none | done, PR #1578 |
| TEST.133 | A reusable big-galaxy query budget test: EXPLAIN every page and list query on 2,000,000 systems | none | done, PR #1405 |
| TEST.134 | test_many_processes_starting_at_once_start_one_job fails under a parallel run (bug) | none | done, PR #2124 |
| TEST.135 | Big-galaxy query-budget test for the nebula reads | none | open |
| TEST.136 | A nebula-rate audit test: the nebula counts and fillings stay inside the observed bands | none | done, PR #1239 |
| TEST.137 | test_a_neighborhood_from_the_generate_page[core] fails once under load (bug) | none | open |
| TEST.138 | test_open_map_menus_hold_no_overlap[web.galaxy-1280] fails once in a full parallel run (bug) | none | merged into TEST.126 (Boss, 2026-10-16, 11:58 PM EDT: consolidation sweep) |
| TEST.139 | Rework, optimise and rebuild the whole test suite so it runs faster and in modules | none | done, PR #1605 |
| TEST.140 | test_a_body_is_reached_by_its_id_and_a_row_number_names_nothing fails once in a full parallel run (bug) | none | open |
| TEST.141 | test_interrupting_a_parallel_galaxy_run_leaves_no_half_written_sector[15-True] fails once in a full parallel run (bug) | none | open |
| TEST.142 | test_a_table_keeps_its_sort_and_filter_across_a_reload fails on main: its fake get_systems lacks the new `after` argument (bug) | none | done, PR #1313 |
| TEST.143 | test_install_checkout_leaves_a_current_editable_install_alone fails when the suite runs from a git worktree (bug) | none | done, PR #1313 |
| TEST.144 | test_every_plain_button_does_something[web.sector] failed once in a full parallel run (bug) | none | merged into TEST.126 (Boss, 2026-10-16, 11:58 PM EDT: consolidation sweep) |
| TEST.145 | Subprocess tests must run the checkout under test, not the installed copy of another checkout (bug) | none | done, PR #1367 |
| TEST.146 | Rate-limit tests fail when the between-test reset wipes every limiter key on the shared Redis (bug) | none | done, PR #1367 |
| TEST.147 | Two nebula-field table tests fail under load and pass alone (bug) | none | done, PR #1405 |
| TEST.148 | Investigate resuming a test run where it left off when the build is the same | none | merged into TEST.208 (Boss, 2026-10-16, 11:58 PM EDT: consolidation sweep) |
| TEST.149 | Two physics tests fail under load and pass alone (bug) | none | done, PR #1405 |
| TEST.150 | Test run progress log and build id for resuming a run | none | done, PR #1431 |
| TEST.151 | scripts/test.py --resume: skip the tests that passed in a run of the same build | none | done, PR #1484 |
| TEST.152 | Clear leftover Redis keys at the start of a resumed run; merge the dependency map; document resuming | none | done, PR #1520 |
| TEST.153 | A single hung test can stall a whole run (bug) | none | done, PR #1451 |
| TEST.154 | A flaky test: the first facility's host is not the first planet of the random system in test_a_body_is_reached_by_its_id_and_a_row_number_names_nothing (bug) | none | done, PR #1486 |
| TEST.155 | Wording and index agree: property test | none | done, PR #2206 |
| TEST.156 | Rewrite the tests that name ring, layer and slot | none | done, PR #2050 |
| TEST.157 | A flaky test: test_system_detail_names_the_containing_nebula fails now and then (bug) | none | done, PR #1572 |
| TEST.158 | The math check tests the spherical shell, band and wedge grid | none | done, PR #1697 |
| TEST.159 | Whole-sphere tests: required stars per nebula kind, Wolf-Rayet counts, tiling at the reach, profile and empty-sector determinism, and the pre-placed share per population | none | done, PR #2504 |
| TEST.160 | Representation audit: seeded sectors in each region class show every kind at its expected frequency, with systems and required contents checked | none | done, PR #2504 |
| TEST.161 | Fill-stage and system-pass tests: stream independence, same seed same system, and the probability table | none | done, PR #3957 |
| TEST.162 | Multiplicity audit: companion states, row bookkeeping, stability, determinism, no planets in a scatter, fill reuse, motion and the census anchors | none | open |
| TEST.163 | test_nav_ends_picked_from_bookmarks_on_every_page fails: #course-heading count is 0 on the nav route page (bug) | none | done, PR #1801 |
| TEST.164 | Catalogue type audit: type counts, region gates, the N tail, interaction consistency and group decay | none | open |
| TEST.165 | Re-measure the region query on a real generated galaxy, warm and cold | none | done, PR #3980 |
| TEST.166 | The cell key follows an object across a cube boundary, and the point columns and the key round-trip bit for bit | none | done, PR #1881 |
| TEST.167 | test_random_start_mode_respects_min_start_density hit a 900 s timeout in a full run and passes alone (bug) | none | merged into TEST.111 (Boss, 2026-10-16, 11:58 PM EDT: consolidation sweep) |
| TEST.168 | Control sweep: a browser test drives every control on every page against a spherical galaxy with 96-bit IDs | none | open |
| TEST.169 | Admin page browser tests: every admin page and task works end to end after the streamline | none | open |
| TEST.170 | test_math_check::test_whole_run_is_fast and the Galaxy Map breadcrumb folding test fail under full-suite load and pass alone (bug) | none | done, PR #2124 |
| TEST.171 | test_gen_resume::test_sectors_a_forced_scatter_skipped_fill_correctly_afterwards failed once under gate load with no scattered system (0 > 0) and passes alone (bug) | none | done, PR #2124 |
| TEST.172 | test_galaxy_map_forward_to_current_jumps_to_the_newest_view failed once: clicking a Galaxy Map block with 0 sectors charted did nothing, and it passes alone (bug) | none | done, PR #2124 |
| TEST.173 | test_gen_bright_scatter_edges::test_a_fill_after_a_shell_failed_mid_scatter_builds_no_leftover_star failed once under gate load with an empty ring in _nearest (bug) | none | merged into TEST.171 (Boss, 2026-10-16, 11:58 PM EDT: consolidation sweep) |
| TEST.174 | test_web_a11y::test_page_is_accessible_and_fits[web.sector-desktop-light] failed once under gate load and passes alone (bug) | none | done, PR #2124 |
| TEST.175 | Edit everything end to end: the headline scenarios with approval and undo | none | open |
| TEST.176 | Plan engine properties: any sequence of edits leaves every system valid, and undo restores every row | none | open |
| TEST.177 | Edited sectors survive repair, regeneration and the golden seed | none | open |
| TEST.178 | Roles end to end: user, moderator and admin through the whole request chain | none | open |
| TEST.179 | Permission matrix test: every route and page checked for every role against the catalogue | none | open |
| TEST.180 | Backup and repair end to end: take a full archive, lose the control schema and the live parity file, restore and repair, and prove the galaxy and the logins come back | none | open |
| TEST.181 | The order-independence test: the same sectors at 1 and 4 workers and in two orders | none | done, PR #2628 |
| TEST.182 | Put one block-first backfill into the golden galaxy | none | open |
| TEST.183 | A test that fails when a new output path formats its own lines | none | done, PR #4450 |
| TEST.184 | Whole-sphere tests, nebula and Wolf-Rayet side: contents per class, the Wolf-Rayet count, tiling, determinism and the King solver | none | done, PR #2504 |
| TEST.185 | Whole-sphere tests, component side: pre-placed share per population, no double count of a cluster star, nothing hosted in a halo or old cluster cell, D and E members inside their radius | none | done, PR #2504 |
| TEST.186 | Census audit harness: fill seeded sectors per region class, count by kind and stage, compare to expectation within k Poisson sigma, print one report | none | done, PR #2040 |
| TEST.187 | Browser-test failure reports are buried in the werkzeug request log: set that logger to WARNING in the test server | none | done, PR #2050 |
| TEST.188 | The skimage marching-cubes call reads `.shape` in a way NumPy 2.5 deprecates: change the call and keep the warning out of the log | none | done, PR #2050 |
| TEST.189 | A test that runs `migrate` leaves maintenance mode on: restore `state/maintenance.json` in the fixture (bug) | none | done, PR #2124 |
| TEST.190 | A failing session fixture makes every dependent test wait its full 900 s timeout: fail the dependents at once (bug) | none | open |
| TEST.191 | The plan tests write real settings JSON into `/var/lib/planetGen/jobs/galaxy-settings` and log "No naming key drawn ... Unknown database": give them a temporary directory and a test database | none | done, PR #2070 |
| TEST.192 | The galaxy-grid fuzz test `test_skeleton_bound_is_an_upper_bound_at_every_angle` fails on a float cancellation of the halo subtraction far above the disk (bug) | none | done, PR #1883 |
| TEST.193 | A module fixture's `DROP DATABASE IF EXISTS` blocks for 853 s on a metadata lock that another connection holds, so a browser test hangs in teardown (bug) | none | done, PR #1997 |
| TEST.194 | `test_system_map_selection_drill_and_measure` failed once in the browser phase of a full run under load and passed alone: find the cause and make it robust (bug) | none | open |
| TEST.195 | `test_opening_a_sector_is_one_continuous_view_with_a_filled_panel` timed out at 30 s fetching `/sector/<id>/galaxy` inside a full gate run and passed alone (bug) | none | merged into TEST.193 (Boss, 2026-10-16, 11:58 PM EDT: consolidation sweep) |
| TEST.196 | Close-binary check: generated close-binary fractions by primary mass and [Fe/H] bin against the published ranges, with out-of-range bins flagged | none | done, PR #2642 |
| TEST.197 | `test_page_is_accessible_and_fits[web.sector-phone-light]` failed under load and passed alone: find the cause and make it robust (bug) | none | merged into TEST.174 (Boss, 2026-10-16, 11:58 PM EDT: consolidation sweep) |
| TEST.198 | `test_galaxy_map_does_not_prefetch_on_a_data_saving_connection` failed with a `/galaxy/tiles` net::ERR under load and passed alone: find the cause and make it robust (bug) | none | open |
| TEST.199 | `test_skeleton_bound_is_an_upper_bound_at_every_angle` fails on main for a far-tail case: the halo subtraction loses precision, 2.4e-9 against a 1e-9 tolerance (bug) | none | done, PR #2124 |
| TEST.200 | The nightly local offset limit is documented and pinned by `test_local_offset_after_a_rotation_is_within_one_ulp_of_the_absolute_point` | none | done, PR #3962 |
| TEST.201 | `block_stats` totals do not depend on the order of the rows: a test | none | done, PR #3959 |
| TEST.202 | `test_galaxy_map_drill_down_by_clicks` failed under load and passed alone: call 28 s to 52 s on a loaded worker (bug) | none | merged into TEST.172 (Boss, 2026-10-16, 11:58 PM EDT: consolidation sweep) |
| TEST.203 | `test_one_worker_leaves_the_runs_own_random_stream_alone` failed once under heavy load with WorkerDied "the task's job ended canceled" and passes alone (bug) | none | done, PR #2124 |
| TEST.204 | Family plausibility tests: masses, sizes, indices, bulge shares, bar fractions and arm pitch fall inside the published ranges | none | open |
| TEST.205 | `test_no_reserved_word_is_an_unquoted_identifier` fails on main: the DB.18 migration helpers use LOCK in `ALGORITHM=INPLACE, LOCK=NONE` (bug) | none | done, PR #1927 |
| TEST.206 | Superbubble density tests: rim excess, interior depletion, old stars unchanged, scatter against fill and determinism | none | open |
| TEST.207 | The page_targets fixture of the web accessibility tests errors 24 tests: web.admin_backup answers 302 and 404 with the sample database (bug) | none | withdrawn: fixed before it was filed (PR #1921, NO_SAMPLE_OBJECT, 2026-10-16, 9:26 PM EDT) |
| TEST.208 | `scripts/test.py --resume` reuses an earlier pass across builds when none of that test's own dependencies changed | none | done, PR #2124 |
| TEST.209 | One consistency sweep over every planet and moon class: composition, density, temperature, ozone, aerosols and text all agree | none | done, PR #3613 |
| TEST.210 | A bulk composition sweep over every class: fractions, density, mass, radius and gravity agree | none | done, PR #3613 |
| TEST.211 | A star composition sweep: every kind of star, remnant and brown dwarf, fractions, trends and determinism | none | done, PR #3625 |
| TEST.212 | A radiation sweep over every class: transmission, greenhouse, shielding and dose agree with the composition | none | open |
| TEST.213 | `test_a_systems_disc_fades_in_as_the_camera_comes_up_to_its_star` failed in a loaded browser run: the wheel loop opened the system before the disc reached full share (bug) | none | done, PR #2103 |
| TEST.214 | Every registered possibility table: for every class the gated weights add up to 100 percent, no class has an empty gate, and every row has a source and a licence | none | done, PR #3615 |
| TEST.215 | The Solar System calibration rows live in test files only: a flag fences them and a test fails when a generator table reads one | none | done, PR #3954 |
| TEST.216 | `test_systems_inside_a_nebula_point_at_it_and_the_innermost_wins` failed in the full parallel gate: the nested system's `inside_nebula_id` was None, and it passes alone (bug) | none | done, PR #2124 |
| TEST.217 | `test_every_column_is_written` found a NULL in every row of `moons.evolutionary_speed` in the full parallel gate, and passes alone (bug) | none | done, PR #1989 |
| TEST.218 | `test_first_control_connections_at_once` leaves 3 idle connections open on `planetgen_test_*` when its fixture ends: the test should close its own connections (bug) | none | done, PR #2899 |
| TEST.219 | Browser pages still have 1 to 7 requests running 8 s after the server stops: find what keeps them busy and make it faster | none | done, PR #4052 |
| TEST.220 | test_the_scatter_lists_each_shells_phenomena_by_kind_like_the_star_scatter finds `neutron-star` inside `magnetar neutron-star` (bug) | none | done, PR #2214 |
| TEST.221 | test_a_runaway_system_moves_beyond_the_rotation_curve_and_its_stars_and_bodies_share_it raises StopIteration when the first system's only planet entry is an asteroid belt (bug) | none | done, PR #2070 |
| TEST.222 | Remove the NumPy 2.5 deprecation filter around the skimage marching-cubes call once scikit-image stops raising it | none | open |
| TEST.223 | Two galaxy-map browser tests fail in a loaded full browser run and pass alone: the free camera from a window down, and forward to the current view (bug) | none | done, PR #2252 |
| TEST.224 | `test_systems_inside_a_nebula_point_at_it_and_the_innermost_wins` failed once in a full gate: the nested system's inside_nebula_id was None (bug) | none | done, PR #2124 |
| TEST.225 | `test_random_start_neighborhood_matches_the_real_skeleton_plan` generated 18,210 systems against a predicted 83,346.6 in one gate run (bug) | none | done, PR #2227 |
| TEST.226 | Tests for the plan configuration and its two switches: the interior is identical with each combination, off generates nothing beyond the boundaries, addresses outside them are refused, the 8.1 sphere is the all-on case | none | open |
| TEST.227 | test_a_failed_job_shows_its_error_and_try_again_starts_a_new_one[light-phone] failed in a loaded browser run and passes alone (bug) | none | done, PR #2261 |
| TEST.228 | test_stars_and_phenomena_fit_within_their_sectors_real_cells compares a dark cloud's centre with its sector's cell, but the cloud is filed by its centroid | none | open |
| TEST.229 | The browser phase of a gate stalls for minutes at one test with a headless_shell gpu process near 200% CPU (bug) | none | open |
| TEST.230 | test_o_stars_sit_in_h_ii_regions_and_cool_stars_light_nothing fails in about 10 of 12 runs since the stellar bundle merged (bug) | none | done, PR #2128 |
| TEST.231 | test_the_wheel_anchors_on_the_star_under_the_pointer_not_the_mesh_behind_it failed alone once in 8 runs (bug) | none | done, PR #2147 |
| TEST.232 | The a11y page_targets fixture fails 199 browser tests when /search answers 202 (slow-page job) under load (bug) | none | done, PR #2149 |
| TEST.233 | test_a_measured_feature_moves_by_the_percentage[asteroid_belt-30] failed once in a loaded gate and passed 13 of 13 alone: its seed changes every run (bug) | none | done, PR #2899 |
| TEST.234 | test_the_review_and_error_pages_are_accessible[admin] found 606 px of width on a 390 px page once, at the site gear dropdown (bug) | none | open |
| TEST.235 | test_stop_waiting_stops_the_polling_and_check_again_picks_it_up[light] failed once in a loaded gate and passes alone (bug) | none | done, PR #2218 |
| TEST.236 | test_galaxy_map_drill_down_by_clicks rejects a window address with a negative band, because its regex allows a minus sign only on the wedge (bug) | none | done, PR #2227 |
| TEST.237 | The wheel-anchor browser test failed alone at box load 9.7 (call 5.47 s) after the TEST.231 fix merged (bug) | none | open |
| TEST.238 | test_a_plan_run_stores_the_scatter_and_its_seed_and_a_rerun_replaces_it fails in a gate with no hypervelocity-star stored, and passes alone (bug) | none | done, PR #2252 |
| TEST.239 | Tests for the unified chain: determinism by host, representation counts per kind and level, rogue and comet rates, derived remnant shares | none | open |
| TEST.240 | Round-trip and registry tests for the mass object core: every kind's stored row round-trips, and no kind table is named in a list the registry should own | none | open |
| TEST.241 | Test: no absolute position is stored for anything inside a system, so every 8.1 item that stores a position stores it parent-relative | none | done, PR #2454 |
| TEST.242 | Tests: invariants for any system at any time: zero mass-weighted offset sum, add and remove returns exactly, a parent change keeps every absolute position, one member equals the binary formula bit for bit | none | open |
| TEST.243 | Two test_web_browser_units tests errored in setup, one after 1041 s in pymysql (bug) | none | closed, could not reproduce |
| TEST.244 | Tests: re-primary and reparent invariants, frame round trips, momentum conservation, a moon's velocity as the sum of its hops, a locked body's rate follows a primary change | none | open |
| TEST.245 | test_generate_sector_real_capacity_overflow_returns_a_partial_sector fails in about 1 run in 100 because it draws an unseeded sector, so the draw sometimes fits all 120 systems (bug) | none | done, PR #2227 |
| TEST.246 | The web edit menu tests time out in setup in a loaded gate: the accessibility site database backfill took 1024 s while other tests ran, and the tests pass alone in 46 to 49 s (bug) | none | done, PR #2286 |
| TEST.247 | The random-start neighbourhood test with --min-start-density 1000 (dense_core) takes 466 s alone and 554 s with its sibling, and times out at 900 s under gate load (bug) | none | done, PR #2272 |
| TEST.248 | The slow-pages test 'stop waiting stops the polling and check again picks it up' fails in a loaded browser phase: the status line still reads 'About 30% ... about 5 s left.' where 'Stopped waiting.' is expected, and it passes alone in 0.3 s (bug) | none | done, PR #2261 |
| TEST.249 | The sector and system page query-count test fails once in the non-browser phase: 12 statements for a lone system where 11 for one with planets, on unseeded systems (bug) | none | done, PR #2272 |
| TEST.250 | test_workers.py: two tests fail once in a full run because ensure_workers is the inline test stand-in although the module is marked real_queue (something left the inline queue patch active in that xdist worker) (bug) | none | done, PR #2267 |
| TEST.251 | test_orbital_facilities_advance_their_phase fails once in a gate: advance_facility_orbits returns 0 where 1 is expected for orbit_clock(period / 4) (bug) | none | done, PR #2272 |
| TEST.252 | test_a_double_click_on_a_star_flies_the_camera_to_it fails once in a loaded gate and passes alone (bug) | none | done, PR #2286 |
| TEST.253 | test_nebula_shape_endpoint_serves_a_mesh fails once under load: the largest vertex radius is 1.0288 against a bound of 1.0, and it passes alone (bug) | none | done, PR #2332 |
| TEST.254 | test_core_ball_batch_fill_makes_planets builds 4 systems instead of 3 for some random galaxy seeds (bug) | none | open |
| TEST.255 | test_workers: two real-Redis worker tests fail in loaded gates and pass alone, again after TEST.250's fix (bug) | none | done, PR #2332 |
| TEST.256 | test_black_hole_density: the black-hole stream test fails in loaded and full gates and passes alone (bug) | none | done, PR #2327 |
| TEST.257 | test_sector_page_and_tiles_show_a_phenomenon_the_cli_added asserts an empty cloud list on an unseeded galaxy and sometimes gets a nebula (bug) | none | done, PR #2416 |
| TEST.258 | Other tests that scatter a low shell of the toy shape may hold gigabytes after GEN.404: sweep them for resident memory | none | done, PR #2461 |
| TEST.259 | test_open_map_menus_hold_no_overlap audits the menu boxes while the open animation is still moving and reports an overlap that is gone a frame later (bug) | none | done, PR #2332 |
| TEST.260 | test_galaxy_absurd_radius_is_clean[200.0] hangs past its 10 s limit on main and leaves the run lock set, so every later test that starts a run fails (bug) | none | done, PR #2425 |
| TEST.261 | Release gate: a full bug hunt audit of the whole code base and its behaviour, done before every major or minor version change (8.1 first) | none | open |
| TEST.262 | test_a_double_click_flies_even_when_the_page_is_slow_to_answer_the_first_click fails once under load: the camera target stays at (3.56, -1.36, -0.29) instead of (0, 0, 0) (bug) | none | open |
| TEST.263 | test_every_check_rejects_a_bad_row fails on main: the seeded rich galaxy has no planet_reflection_spectrum row to break planet_reflection_spectrum_chk_1 on (bug) | none | done, PR #2416 |
| TEST.264 | test_web_a11y: every page test errors on main because the new /lists/<kind> route has no sample value in sample_params (bug) | none | done, PR #2425 |
| TEST.265 | A test database drop ended another worker's query: "Lost connection to MySQL server during query" and 503 "database unavailable" in parallel gates (bug) | none | done, PR #2442 |
| TEST.266 | test_interrupting_a_parallel_galaxy_run_leaves_no_half_written_sector[15-False] fails under load with "run ... is running" (Busy) and passes alone (bug) | none | done, PR #2466 |
| TEST.267 | test_object_list.py::test_a_walk_of_the_pages_visits_every_object_once_in_the_order_of_the_sort[False-mass] fails once with "the galaxy has too few planets for a page walk" on unseeded galaxy data (bug) | none | done, PR #2485 |
| TEST.268 | Four test_galaxy_gen tests for an address outside the outline fail on main: generate exits with status 1 (bug) | none | done, PR #2502 |
| TEST.269 | test_only_the_newest_archives_are_kept compares archive names as plain strings, so two saves in one second fail it (bug) | none | done, PR #2494 |
| TEST.270 | check-db reports a damaged revision during a full parallel gate and fails two db_check tests (flake) | none | done, PR #2573 |
| TEST.271 | test_drop_database_holders gentle_wait fails on a clean main with MySQL "Lock wait timeout exceeded" (bug) | none | done, PR #2441 |
| TEST.272 | test_a_star_in_an_opened_sector_is_selectable (browser) fails in a full gate on main (bug) | none | done, PR #2573 |
| TEST.273 | test_the_index_takes_its_size_off_the_tile_cache_budget fails once in a full run: max_cache_bytes 1073577984 against 1073741824 (flake) | none | done, PR #2899 |
| TEST.274 | test_orbital_limits.py::test_a_circular_body_far_from_its_epoch_is_where_the_remainder_puts_it fails on main: _relative_circular_km was removed by #2481 (bug) | none | done, PR #2493 |
| TEST.275 | test_fuzz_web_routes shape-report case runs past the 900 s pytest limit on main: /admin/generate/shape-report under hostile query strings (bug) | none | done, PR #2600 |
| TEST.275 | test_fuzz_web_routes shape-report case runs past the 900 s pytest limit on main: /admin/generate/shape-report under hostile query strings (bug) | none | done, PR #2600 |
| TEST.276 | test_clearing_a_scatter_with_nothing_stored_says_so_and_changes_nothing fails once under parallel load with empty captured output (flake) | none | done, PR #3732 |
| TEST.277 | test_a_dead_worker_never_hangs_the_run times out once each for [1] and [2] under heavy load (flake) | none | done, PR #2600 |
| TEST.278 | test_column_units fails on main: 36 columns of orbital_events and orbital_event_bodies are not declared in column_units.py (bug) | none | done, PR #2510 |
| TEST.279 | Eight phenomenon page tests answer 404 on main: test_web_system_phen.py and test_class_pages.py (bug) | none | done, PR #2540 |
| TEST.280 | Five tests are refused as Busy in a full gate because workers share the default Redis database (race) | none | done, PR #2573 |
| TEST.281 | test_a_star_in_an_opened_sector_is_selectable fails on main: no click on a star in an opened sector selects it (bug) | none | done, PR #2619 |
| TEST.282 | test_precision imports calculate_hill_sphere from generation.star, where GEN.565d moved it away: ImportError on main (bug) | none | done, PR #2547 |
| TEST.283 | Re-run the whole-sphere audits TEST.185, TEST.184, TEST.159 and TEST.160 on main after the PERF.146 speedups | none | done, audits run on main dd6d46d0 (9:57 PM EDT) |
| TEST.284 | test_an_uncharted_black_hole_page_names_its_class fails under parallel load in the changed gate: the page renders without the class name (bug) | none | open |
| TEST.285 | test_a_real_run_of_the_plan_writes_the_json_and_drops_its_databases fails on a machine without the static downloads: the plan step exits 1 and the test expects 0 (bug) | none | done, PR #2602 |
| TEST.286 | test_quasar_physics_follows_from_mass_and_eddington_ratio fails by chance: it builds 50 quasars with no seed and a mass draw below 1e5 breaks its assertion (bug) | none | done, PR #2612 |
| TEST.287 | test_web_browser_fixture_maps.py::test_galaxy_breadcrumb_always_matches_the_view_its_url_names times out in the browser phase under load (bug) | none | open |
| TEST.288 | test_work_queue_admin.py::test_a_paused_queue_node_stands_by_until_resumed[2] fails in the four-worker non-browser phase: the lease holder is another test worker's node (bug) | none | done, PR #2625 |
| TEST.289 | Two API tests fail once under the gate's four workers: test_a_neighborhood_job_runs_its_sectors_through_the_wrapper and test_a_quick_edit_answers_in_the_same_response (bug) | none | open |
| TEST.290 | test_work_units and test_galaxy_seed each fail once under the gate load: a galaxy run stores each stage once, and a new seed replaces the old only before any sector (bug) | none | open |
| TEST.291 | test_a_star_in_a_sector_far_from_the_core_is_selectable fails once in the gate and passes alone (bug) | none | done, PR #2669 |
| TEST.292 | test_web_browser_sector_reach.py::test_a_systems_disc_fades_in_as_the_camera_comes_up_to_its_star fails at 130 s under load 9.7 and takes 97 s alone (bug) | none | open |
| TEST.293 | test_comet_population.py::test_a_draw_is_deterministic_from_its_seed_at_one_and_four_workers fails once under the gate load and passes alone (bug) | none | done, PR #2663 |
| TEST.294 | full --check run on a main checkout has an empty change set and reports every failure as not selected (the dependency map itself selects test_validation.py::test_system_methods_still_delegate for the GEN.325 change) (bug) | none | done, PR #2663 |
| TEST.295 | test_parallel_generation.py::test_one_worker_draws_what_two_and_four_workers_draw hit the 900 s timeout waiting for its queue workers in a saturated gate (bug) | none | open |
| TEST.296 | test_species_limits.py::test_a_generated_galaxys_species_are_all_inside_their_limits fails once in a gate: population.run_pass makes 1 new species where 2 or more are expected (bug) | none | done, PR #2687 |
| TEST.297 | The 13 tests of test_work_queue_failures.py and test_web_browser_sector_open.py::test_opening_a_sector_is_one_continuous_view_with_a_filled_panel fail in a gate at load 6 to 8 and pass alone (bug) | none | open |
| TEST.298 | test_admin_script_cli.py: four fixture-finalizer errors under load in a gate (bug) | none | open |
| TEST.299 | test_api_species_tech_filters fails once in a full run: the first species of the tech-descending list has tech 6.55, not the 6.0 the test set (bug) | none | done, PR #2945 |
| TEST.300 | End-to-end browser tests for the habitability system: symbols on lists, sector colouring and habitable-world markers | none | open |
| TEST.301 | Tests for the condition-based classes: every planet and moon class from its conditions, class and PHI-6 agreement, reclassification after an edit, determinism | none | done, PR #3619 |
| TEST.302 | benchmark_suite.stop_workers() kills other tests' queue workers: test_benchmark_suite fails or breaks test_queue_pause, test_reset_jobs, test_workers and test_api_jobs with exit -9 when the data sets are fetched (bug) | none | done, PR #4042 |
| TEST.303 | test_update_database_step::test_a_failed_migration_stops_the_step_and_keeps_the_version fails on untouched main on the lane machines because the test database user lacks GRANT (bug) | none | dropped (Debug mode thread, 2:21 PM EDT 2026-10-24 as the lane reported: its own test MySQL user lacked GRANT OPTION) |
| TEST.304 | Nine Playwright test files are not marked browser, so they run in the all-core phase | none | done, PR #2945 |
| TEST.305 | The sector-map browser tests wait 120 s for the star list and fail on a fresh four-core box (bug) | none | closed, no change (not reproduced: the file passes 103 of 103 alone on a clean four-core machine, Bug Hunt 7:38 AM EDT 2026-10-24; load-bound, see TEST.302) |
| TEST.306 | Five physics modules are imported only by tests and no item names them | none | done, no code (2:59 AM EDT 2026-10-24: Bug Hunt checked, nothing deleted; GEN.197, GEN.110, GEN.479 and GEN.428 name the modules) |
| TEST.307 | test_systems.py::test_asteroid_belt_and_habitable_world_forced_together_both_succeed[A0V] fails about once in 5 to 8 runs (bug) | none | open |
| TEST.308 | Body rendering tests: look record stability, painter determinism, tier switching and the memory ceiling | none | dropped (Boss, 10:26 PM EDT 2026-10-23: tests ride with each body-rendering item) |
| TEST.309 | Delete the dead, duplicate and tombstone regression tests, the unused test imports and variables, and the unused galaxy_shape_visualizer_cli.py tool | none | done, PR #2771 |
| TEST.310 | Re-run the 16 tests that failed in the Blocker 5 gate on current main in an environment with MySQL and Redis, and report | none | done, PR #2837 |
| TEST.311 | test_hosted_nebula_rows::test_a_pre_placed_star_rolls_the_nebula_its_sector_builds fails on main (bug) | none | done, PR #2945 |
| TEST.312 | test_datasets_real::test_the_catalogue_reports_nothing_missing fails on main for the new ozone_uv dataset (bug) | none | done, PR #2810 |
| TEST.313 | test_systems.py test_asteroid_belt_and_habitable_world_forced_together_both_succeed[M5V] is a random flake: about 2 in 12 runs the belt is missing (bug) | none | done, PR #2918 |
| TEST.314 | Only classify writes the class column, and the class is recomputed when an input is edited | none | open |
| TEST.315 | TEST.210 reports the class shares and ranges as outputs; TEST.209 adds that the class agrees with the conditions | none | done, PR #3568 |
| TEST.316 | test_no_page_shows_a_formula_that_is_not_written_out fails on main: /classes/habitability shows a bare "CO2" (bug) | none | done, PR #2850 |
| TEST.317 | test_a_systems_disc_fades_in_as_the_camera_comes_up_to_its_star fails under heavy load (bug) | none | open |
| TEST.318 | test_web_generate_size::test_the_line_is_the_estimate_only_output_for_the_same_values fails on untouched main: refused True vs False (bug) | none | done, PR #3271 |
| TEST.319 | test_datasets_real::test_the_catalogue_reports_nothing_missing fails when the ozone_uv file is not downloaded (bug) | none | open |
| TEST.320 | test_db_batching::test_sector_rows_are_written_in_a_few_multi_row_inserts fails once in a full gate and passes alone (bug) | none | open |
| TEST.321 | test_web_browser_tables: the jump box test flakes under load (bug) | none | open |
| TEST.322 | test_no_todo_codes_shown fails on main: planetgen/db/point_columns.py:189 shows the item code DB.73 in a string the code shows (bug) | none | done, PR #2881 |
| TEST.323 | test_web_browser_fixture_maps.py: 41 tests time out waiting for .starmap-sr-list button or find nothing to click on the Galaxy Map (bug) | none | closed, no change (not reproduced: the file passes 103 of 103 alone on a clean four-core machine, Bug Hunt 7:38 AM EDT 2026-10-24; load-bound, see TEST.302) |
| TEST.324 | test_web_header_menus top-bar [1280x720-light] times out under load (bug) | none | open |
| TEST.325 | test_migration_habitable_worlds: the v135 galaxy count test fails on main because it downgrades to v138 while the DB.73 v140 marker still reads 140 (bug) | none | done, PR #2881 |
| TEST.326 | test_systems.py::test_asteroid_belt_and_habitable_world_forced_together_both_succeed[M5V] fails about 1 run in 12: a forced asteroid belt is not always made (bug) | none | dropped (duplicate of TEST.313, TODO 4:07 AM EDT 2026-10-24) |
| TEST.327 | test_web_browser_tables::test_search_planet_rows_draw_their_habitability_symbols_with_and_without_scripts fails on main: it counts 0 a.help-link elements and expects 7 since UX.180 to UX.182 (bug) | none | done, PR #2892 |
| TEST.328 | Tests that downgrade the schema by hand need an undo step for every new migration: derive the undo from one helper instead of three hand lists | none | done, PR #3693 |
| TEST.329 | test_web_a11y::test_page_is_accessible_and_fits[web.system-*] fails on main: the planet list class links have 1.06:1 contrast and sit as focusable links inside a summary (bug) | none | done, PR #3011 |
| TEST.330 | test_web_formulas::test_no_page_shows_a_formula_that_is_not_written_out fails on (/classes/habitability, CO2) (bug) | none | done, PR #2850 |
| TEST.331 | test_web_browser_tables search_planet_rows_draw_their_habitability_symbols fails: no a.help-link in the search planet rows (bug) | none | open |
| TEST.332 | test_object_rows::test_planets_moons_belts_and_comets_are_read_back_as_the_dicts_they_were_written_from fails on main: the comet is missing from the row-kind set (bug) | none | open |
| TEST.333 | test_systems forced belt and habitable world together still fails for M5V about 1 trial in 15 after the TEST.313 fix (bug) | none | open |
| TEST.334 | Queue and generation tests pass alone but hang or "DID NOT RAISE" in a full run on lane machines: test_gen_resume (7), test_gen_bright_scatter_edges (5), test_work_units (7 in sequence), test_worker_warm, test_work_queue parallel (bug) | none | dropped (Debug mode thread, 2:21 PM EDT 2026-10-24 as the lane reported: not red on main; its lane machine ran without PLANETGEN_TEST_REDIS_URL so real workers ran) |
| TEST.335 | test_order_independence::test_a_sector_reached_by_the_backfill_comes_out_the_same_in_either_order fails 3 runs in 8 on clean main (bug) | none | dropped (not reproduced on main by Bug 1, 8 of 8 pass, 7:14 PM EDT 2026-10-24 as the lane reported) |
| TEST.336 | test_migration_habitable_worlds::test_a_v135_galaxy_gets_its_counts_from_the_planets_and_moons_it_holds fails on main: 142 == 138 (bug) | none | open |
| TEST.337 | test_gen_cell_resolution::test_a_cell_above_the_ceiling_is_filled_down_to_the_cap_and_no_further fails on unmodified main: minimum luminosity 5.4e28 against the required 1.9e29 (bug) | none | dropped (not reproduced on main by Bug Hunt, 10:41 AM EDT 2026-10-24 as the lane reported) |
| TEST.338 | test_bright_star_scatter::test_plan_scatter_stores_stars_and_fill_builds_their_systems fails on main: the built primary differs from the stored row in 2 of 7 systems (bug) | none | done, PR #2986 |
| TEST.339 | test_systems::test_asteroid_belt_forcing_is_reliable fails about 1 run in 4: 2 of 60 trials without a belt, not seeded (bug) | none | dropped (not reproduced on main by Bug Hunt, 10:41 AM EDT 2026-10-24 as the lane reported) |
| TEST.340 | Three more tests fail on clean main, likely the same cause as TEST.338: test_gen_stored_trees fill_builds_each_row, test_hosted_nebula_rows pre_placed_star, test_emission_members pre_placed_hot_star (bug) | none | done, PR #2986 |
| TEST.341 | A spatial-position test failed once with "no system matching" in the lane gate and passed on rerun (bug) | none | open |
| TEST.342 | test_api::test_galaxy_tiles_carry_no_comets_rogue_planets_or_asteroid_fields fails on unmodified main (bug) | none | done, PR #3009 |
| TEST.343 | Three tests red on main since the habitable_systems table (MAP.220, #2945): point columns, parent-relative positions, and VALUES() in migration 0148 (bug) | none | done, PR #3071 |
| TEST.344 | Gate: browser and WebGL tests run as their own serial phase, with a worker cap by free memory | none | open |
| TEST.345 | test_gen_stored_trees::test_a_fill_builds_each_row_from_its_stored_tree fails on main on a lane machine (bug) | none | done, PR #4018 |
| TEST.346 | web.system a11y variants (desktop and phone, light and dark) fail on lane machines of Core 3 and Core 5 but pass for Bug Hunt (bug) | none | open |
| TEST.347 | test_admin_edits::test_system_page_regenerates_a_planet and ::test_system_page_changes_a_class fail alone on main: the system page lacks the 'Regenerated ...' and 'Change class' text (bug) | none | done, PR #3050 |
| TEST.348 | test_class_pages: test_phenomenon_page_links_its_class (5 cases) and test_phenomenon_unknown_class_stays_text fail on clean main (bug) | none | done, PR #3050 |
| TEST.349 | The browser phase of scripts/test.py changed fails in fixture setup with KeyError 'planets' before any test body runs, on every lane and on main (bug) | none | done, PR #3068 |
| TEST.350 | Tests red on clean main: test_web_facilities (18 tests), test_web_population (test_system_page_names_species_and_territory, test_real_population) and test_web_browser_controls test_every_plain_button_does_something[web.system] (bug) | none | done, PR #3077 |
| TEST.351 | test_property_sheet::test_every_column_is_on_the_sheet_or_says_why_it_is_not fails on main for planet and moon: GEN.699's columns internal_heat_flux_w_m2, mantle_activation_kj_mol and mantle_eta0_pa_s are neither on the sheet nor in NOT_SHOWN (bug) | none | done, PR #3059 |
| TEST.352 | test_db_column_round_trip::test_every_column_is_written fails on main: the star columns disk_kind, disk_mass_solar, disk_inner_au and disk_outer_au are now written, so they must leave NULL_IN_THIS_GALAXY (bug) | none | done, PR #3103 |
| TEST.353 | test_interior_heat::test_every_generated_rocky_body_stores_its_draws_and_an_urey_ratio_that_is_an_output is red on main: a generated body has urey_ratio 2.09 and the test caps it at 1.5 (bug) | none | done, PR #3103 |
| TEST.354 | test_api.py::test_the_description_splits_into_every_body is red on clean main: the facts section appears twice (facts, prose, composition, ratings, facts) (bug) | none | done, PR #3098 |
| TEST.355 | test_class_pages: the phenomenon link (nebula D) is red on clean main (bug) | none | done, PR #3098 |
| TEST.356 | test_web_search: the star result link is red on clean main (bug) | none | done, PR #3098 |
| TEST.357 | changes/map-unresolved-glow.minor.md has a '### Tests' heading, so bump_version.py --check-pr fails on main (bug) | none | done, PR #3103 |
| TEST.358 | test_property_sheet planet and moon "every column shown or excused" fails on main for the interior columns of GEN.642 and DB.89 (bug) | none | done, PR #3123 |
| TEST.359 | test_object_page_walk fixtures error on main: the system detail has no "stars" (bug) | none | done, PR #3118 |
| TEST.360 | test_nebula_blocks_db.py::test_a_plan_makes_the_totals fails on main (bug) | none | dropped (not reproduced on main by Bug 1, 6:35 PM EDT 2026-10-24 as the lane reported) |
| TEST.361 | test_no_todo_codes_shown.py fails on main: MAP.220 is named in a shown description in db/point_columns.py line 177 (bug) | none | dropped (not reproduced on main by Bug 1, 6:35 PM EDT 2026-10-24 as the lane reported) |
| TEST.362 | test_object_rows read-back fails on main with KeyError partial_pa (bug) | none | done, PR #3155 |
| TEST.363 | test_protection describer disagrees with habitability_explain on main (bug) | none | done, PR #3155 |
| TEST.364 | Clean up the unused imports, unused variables and redefinitions that ruff and vulture report (82 ruff F findings, 12 vulture hits) | none | done, PR #3154 |
| TEST.365 | Six phenomenon-page tests in test_class_pages.py fail on clean main (bug) | none | done, PR #3142 |
| TEST.366 | test_memwatch::test_a_step_that_holds_far_more_memory_than_its_budget_fails_the_smoke_run fails on clean main: it counts the host MySQL and returns three rows instead of one (bug) | none | dropped (passes 9 of 9 on main since #3125 changed the memory budgets, Bug 1, 7:20 PM EDT 2026-10-24 as the lane reported) |
| TEST.367 | Four admin-page tests fail on origin/main: two system-page edit tests and two admin-page redirect tests (bug) | none | done, PR #3291 |
| TEST.368 | A fixture galaxy with a dense filled region for the map browser tests: zoom, labels, picking and the system infobox (bug) | none | done, PR #3171 |
| TEST.369 | test_age_needs (two tests) and test_eta_estimator fail on clean main (bug) | none | done, PR #3155 |
| TEST.370 | test_web_pages in_process_transport fails: DatabaseLoginError, there is no planetgen database account on the lane machine (bug) | none | done, PR #3293 |
| TEST.371 | test_workers::a_worker_that_is_killed fails under load and passes alone (bug) | none | done, PR #3291 |
| TEST.372 | scripts/test.py and docs/reference/testing.md call `full` "the merge gate", which lanes read as a whole-suite gate for every merge (bug) | none | done, PR #3192 |
| TEST.373 | test_protection::test_there_is_one_describer_of_the_six_categories fails on main: habitability_explain.py exists again after GEN.445 (bug) | none | done, PR #3155 |
| TEST.374 | scripts/bump_version.py --check-pr fails on main-side change notes with a "### Tests" heading: gen398, gen691, gen719 (bug) | none | dropped (not reproduced on main: bump_version --check passes and the notes with "### Tests" were stamped out, Bug 1, 10:31 PM EDT 2026-10-24 as the lane reported) |
| TEST.375 | Tests Lane L reports red on clean main: test_api_auth_sweep (2), test_class_pages phenomenon links (6), test_db_column_round_trip, test_save_progress sector save, test_web_generate generate_anyway (2), test_memwatch smoke run (bug) | none | done, PR #3219 |
| TEST.376 | scripts/test.py checks MySQL, Redis and the test environment before and during every run, and restarts and resumes after a container kill | none | done, PR #3192 |
| TEST.377 | test_cache_warmup::test_a_big_unit_stops_part_way_when_the_cache_fills_and_counts_the_rest_as_left_out is red on clean main (bug) | none | done, PR #3291 |
| TEST.378 | test_no_todo_codes_shown fails on main: the composition_origin reason in describe/sheet.py names UX.202 in a shown string (bug) | none | done, PR #3196 |
| TEST.379 | Twelve unrelated tests fail on main: description_splits_into_every_body, eta_estimator queue page, protection describers, wheel scroll, web search star link (bug) | none | done, PR #3291 |
| TEST.380 | First live test pass of species, civilizations, polities and territories on a generated galaxy, before the 8.1 stamp | none | open |
| TEST.381 | First live test pass of nations, cities and ports on a generated galaxy, before the 8.4 stamp | none | open |
| TEST.382 | First live test pass of the world population and the trade routes on a generated galaxy, before the 8.5 stamp | none | open |
| TEST.383 | test_web_browser_sector_reach::test_the_wheel_anchor_on_empty_space_is_where_the_column_of_stars_passes_the_threshold fails alone on main (bug) | none | done, PR #3385 |
| TEST.384 | Failures red on a clean copy of main: habitable_systems point columns and sheet tests, mantle columns, facilities pages, the VALUES() in migration 0148, test_prevalence binaries (bug) | none | closed, does not reproduce on main (Bug 1, 1:37 AM EDT 2026-10-25) |
| TEST.385 | Browser tests error at setup on the 50 per hour rate limit of the API fixture: 259 setup errors in one gate (bug) | none | done, PR #3273 |
| TEST.386 | First live test pass of the resource locations on a generated galaxy, before the 8.4 stamp | none | open |
| TEST.387 | test_db_column_round_trip::test_every_column_is_written fails when a random batch has no disk star: stars.disk_inner_au, disk_kind, disk_mass_solar and disk_outer_au are NULL in every row (bug) | none | open |
| TEST.388 | test_reproducible_draws fails on math.hypot in physics/impact_history.py:140 and physics/iron_branch_tuning.py:79 (bug) | none | done, PR #3665 |
| TEST.389 | test_golden_seed one-versus-four-workers reports a sector-structure difference after epoch 101 (GEN.512); unconfirmed, the clean-main run timed out at 900 s (bug) | none | open |
| TEST.390 | More tests red on a clean copy of main: runner self-tests (on_child), admin queue log route 403, cache warm-up, test_stages, test_web_fresh_after_cli, test_workers killed worker, phenomenon diagram zoom tests, test_web_browser_tables life scores (bug) | none | done, PR #3293 |
| TEST.391 | Two tests red on a clean main: the surface water a planet keeps (runaway_history) and the 95 percent abiotic oxygen bound (atmosphere_mixture_escape) (bug) | none | open |
| TEST.392 | test_prevalence.py takes 66 minutes: it draws 2000 full systems per check and each check runs 4 to 12 minutes, the longest single block of any changed gate (bug) | none | done, PR #3528 |
| TEST.393 | Each gate leaves about 50 orphaned planetgen.cli.worker processes that can fill 16 GB of RAM and get mysqld OOM-killed: tests that start detached workers do not reap them (bug) | none | done, PR #3308 |
| TEST.394 | Six fixture-map browser tests are red on clean main (sector-click flow), and test_web_browser_fixture_maps.py is order-sensitive in one process because dense_site leaks extra worlds into the shared fixture (bug) | none | done, PR #3385 |
| TEST.395 | A planet distribution test with observed-rate bands, run after the class cut-over on the final classes and on galaxy-placed stars: giants and hot Jupiters per host mass, planets per mass bin, 1 to 20 Earth-mass planets per star, habitable-class bodies per star, and a rule for classes O, P and Q | none | open |
| TEST.396 | test_fuzz_web_routes is red on main: the admin queue log route answers 302 where the test expects 403 (bug) | none | done, PR #3329 |
| TEST.397 | test_api_auth_sweep is red on main (bug) | none | done, PR #3418 |
| TEST.398 | The test preflight kills stale planetgen.cli.worker processes, which pick up the units of later runs (bug) | none | done, PR #3690 |
| TEST.399 | test_golden_seed::test_the_pinned_galaxy_is_made_at_one_worker_and_at_four is red on main: the pinned galaxy differs between one worker and four (bug) | none | done, PR #3504 |
| TEST.400 | test_habitability_fill::test_the_command_re_scores_every_system is red on main: Namespace has no attribute 'command' in generate.validate_fill_args (bug) | none | done, PR #3358 |
| TEST.401 | test_atmosphere_mixture_escape (abiotic oxygen) and test_atmosphere_mixture_rocky (delivered-water class) are red on main (bug) | none | done, PR #3381 |
| TEST.402 | test_db_column_round_trip fails once under load: the stars' disk columns read back NULL (bug) | none | done, PR #3689 |
| TEST.403 | Three scatter tests are far too slow: test_a_galaxy_with_no_scatter_still_rolls_its_own runs over 4.5 minutes, and two others hang for an hour under -n 4 (bug) | none | done, PR #4041 |
| TEST.404 | golden_update's first run wrote galaxy.json with the sector section made before the epoch bump, so test_golden_seed failed until a second run on the same tree (bug) | none | done, PR #3358 |
| TEST.405 | The system page tests and the migration 0148 habitable_systems test are red on main (bug) | none | done, PR #3689 |
| TEST.406 | Four tests are red on main: test_surface_air stored-partials sum (1e-9 relative) and test_resume_across_builds, test_progress_log (2), test_redis_cleanup fail with an unexpected keyword on_child (bug) | none | done, PR #3689 |
| TEST.407 | golden_update takes about 55 minutes because its pinned-galaxy test builds the galaxy twice (one worker, then four), so it blocks every generation merge (bug) | none | done, PR #3371 |
| TEST.408 | test_gen_fill_stages::test_a_saved_fill_stores_the_pulsar_systems_and_marks_every_scatter_row_built hangs more than 20 minutes at 98 percent CPU on clean main and the 900 second timeout does not stop it (bug) | none | open |
| TEST.409 | test_web_population::test_real_population errors in setup: "could not generate a system with a civilization", 7.75 minutes, every run on clean main (bug) | none | open |
| TEST.410 | test_planet_physics_fixes::test_class_p_is_colder_than_class_m times out at 900 s: it generates 600 planets at about 5 s each since 4c, 14 of 16 seconds in the air solve (bug) | none | done, PR #3405 |
| TEST.411 | Three map browser tests flake on a loaded 4-core box after a long run and pass alone: the indexeddb revisit test and two sector_reach double-click flight tests (bug) | none | open |
| TEST.412 | src/tests/test_web_ui_fields.py hangs on a lane machine at test_a_bad_distance_is_shown_on_the_distance_field_not_above_the_form, on the nearby page, after 6 passes (bug) | none | open |
| TEST.413 | test_property_sheet fails on main: the GEN.441 ground_*_fraction columns are not on the property sheet (bug) | none | done, PR #3499 |
| TEST.414 | test_habitability_end_to_end_browser errors: the database fixture times out at 900 s before any page loads (8 errors) (bug) | none | open |
| TEST.415 | test_cache_warmup::test_a_big_unit_stops_part_way fails on unmodified main (bug) | none | done, PR #3418 |
| TEST.416 | test_habitability_explain fails ten tests on main: the tests still call the surface_content API as it was before GEN.441 (bug) | none | done, PR #3467 |
| TEST.417 | Whole-sector and whole-galaxy tests take 20 to 45 minutes: make them smaller or cached so a gate run finishes (bug) | none | open |
| TEST.418 | Seven map browser tests are red on clean main (bug) | none | done, PR #3690 |
| TEST.419 | test_opening_a_sector_is_one_continuous_view_with_a_filled_panel failed once in a full file run and passed on the next two (bug) | none | done, PR #3428 |
| TEST.420 | Runner self-tests and the limit-slider follow-again browser test are red on main (bug) | none | done, PR #3690 |
| TEST.421 | test_migration_habitable_worlds: the filled-sector test finds no system with a Blue or Green world in 120 seeds, its _fill_one_sector helper no longer draws one (bug) | none | done, PR #3467 |
| TEST.422 | test_generation_size_estimate (3 tests error in setup) and test_phenomenon_scatter (scatter_lists_each_shells_phenomena_by_kind, a_fill_builds_its_rows_at_their_stored_points) are red on clean main (bug) | none | done, PR #3665 |
| TEST.423 | test_admin_edits::test_recommended_class_change_saves_and_validates fails on clean main: the stand-in body has no eccentricity since migration 0163 (bug) | none | done, PR #3486 |
| TEST.424 | test_property_sheet::test_every_column_is_on_the_sheet_or_says_why_it_is_not[moon|planet] fails on main: the new ground_* columns have no sheet row or NOT_SHOWN excuse (bug) | none | done, PR #3499 |
| TEST.425 | test_wolf_rayet and test_field_units[RoguePlanet] are red on main (unit g/cm2 unknown) (bug) | none | done, PR #3525 |
| TEST.426 | test_galaxy_map_window_follows_the_browser_window_when_it_is_resized fails in long full-file runs and passes alone (bug) | none | done, PR #3687 |
| TEST.427 | Three column registry tests are red on main: aerosol_draw has no unit, tidal_flux_w_m2 is written but listed as NULL, host_eccentricity is not serializable (bug) | none | done, PR #3920 |
| TEST.428 | test_benchmark_suite (real run) and the test_work_units batch fills are red on main (bug) | none | done, PR #3665 |
| TEST.429 | Four test groups are red on clean main: three test_api_job_runs tests, and test_web_seo for web.galaxy_system_nearest (bug) | none | done, PR #3681 |
| TEST.430 | Class census: a standing seeded report of every planet and moon class share and overall colour per star type and zone, with stored expected ranges and anomaly flags | none | open |
| TEST.431 | System census: full systems per star type with planets, moons, belts and comets per system, the zero-planet share, binary share, habitable-world share, and life stage and civilization shares (prevalence suite part B) | none | done, PR #3585 |
| TEST.432 | Star census: star type mix, white dwarfs, subgiants, giants, pre-main-sequence stars, brown dwarfs and multiplicity by type against RECONS and Kroupa (prevalence suite part C) | none | done, PR #3588 |
| TEST.433 | Phenomena census: black holes, neutron stars, magnetars, pulsars, Wolf-Rayet stars, planetary nebulae, supernova remnants, nebula classes, asteroid fields, clusters, rogue planets and interstellar comets against researched densities (prevalence suite part D) | none | done, PR #3621 |
| TEST.434 | Tests red on clean main after the habitable-world markers went on by default: the POI icons browser test, test_no_todo_codes_shown and seven node-JS test files (bug) | none | done, PR #3681 |
| TEST.435 | Prevalence tests for the life classes: the share of worlds in each life class and stage inside researched ranges | none | open |
| TEST.436 | Prevalence tests for the chemistries, the adaptation traits and the animal routes: shares per chemistry per star letter, per route and per trait inside researched ranges | none | open |
| TEST.437 | test_climate_draws test_every_class_stores_the_draws_its_temperature_is_made_of fails for 20 classes: the stored surface temperature differs from the draws' formula by 0.03 K (365.639 vs 365.607) (bug) | none | done, PR #3558 |
| TEST.438 | test_planet_classes test_only_the_known_writers_assign_a_class fails: db/interests.py is a new writer of planet_class and is not in CLASS_WRITERS (bug) | none | done, PR #3562 |
| TEST.439 | test_web_browser_map_perf.py fails on main: the climb moves stars over the buffer at up 1 and up 4, and teardown raises page errors on indexOf and scale (bug) | none | done, PR #3674 |
| TEST.440 | test_sector_map_arrow_keys_turn_the_view_and_keep_the_scale fails on main: the scale text reads 0.25 sectors where 2.5 is expected (bug) | none | done, PR #3571 |
| TEST.441 | Tests red on main from the albedo-greenhouse, moon-tidal-flux and mass-backfill migrations; the white dwarf models dataset is also missing on a fresh test machine (bug) | none | done, PR #3665 |
| TEST.442 | Browser map tests red on main: menu_filters_sit_side_by_side (2 tests), test_system_map_selection_drill_and_measure and test_a_star_in_a_sector_far_from_the_core_is_selectable (bug) | none | done, PR #3639 |
| TEST.443 | Reclassification after an edit: an edited planet or moon changes class when its conditions cross a rule, and keeps it when they do not | none | open |
| TEST.444 | Browser map tests still red on main after TEST.442: the sector_reach double-click tests and the POI icon key test (bug) | none | done, PR #3643 |
| TEST.445 | Four map browser tests in test_web_browser_fixture_maps.py red on main: point phenomena links, black hole markers toggle and legend, nebula click clears the selection, wide map card keeps the page edge (bug) | none | done, PR #3650 |
| TEST.446 | Two tests red on main: the generated config.schema.json out of date, and the batch-fill units listing rescore-habitability first (bug) | none | done, PR #3665 |
| TEST.448 | test_a_sector_is_coloured_by_its_habitable_worlds_and_the_legend_holds_the_stored_count is red on main: blocks show no tint in the new sector view (bug) | none | done, PR #3672 |
| TEST.449 | test_sheet_fill "no empty column that applies" is red on main: life_age_complex_cells_gyr is empty (bug) | none | done, PR #3681 |
| TEST.450 | test_web_browser_sector_open.py::test_opening_a_sector_is_one_continuous_view_with_a_filled_panel is red on main: 80 wheel ticks at a sector centre no longer reach open=1 (bug) | none | done, PR #3879 |
| TEST.451 | Ten reds Core 5 saw on main at its PERF.172 merge (9:59 AM EDT 2026-10-29): re-run each on current main, close the ones already fixed, file the real ones (bug) | none | done, PR #3707 |
| TEST.452 | A sweep finds the life-capable conditions on every planet class and proves the generator lets life evolve there (a class-blind life test) | none | done, PR #3728 |
| TEST.453 | Five tests are red on a clean main: the query near test (NAV.61, expects exit 2), two item_pages planet-page tests, the moon tidal flux migration test, and the fuzz O-star timeout in gas_equilibrium (bug) | none | done, PR #3920 |
| TEST.454 | The observed star mass-function curve as data, with a per-bin harness that measures the generator against it (segments from 0.01 to 150 solar masses by object type) | none | done, PR #3764 |
| TEST.455 | A completeness gate: every numeric constant in the generation and physics modules has a stand-in registry row or is read from a table; the list of unregistered constants may only shrink | none | open |
| TEST.456 | A test that the generated [X/Fe] abundances match the Brewer and Fischer 2016 table, not only its trend and scatter (low priority) | none | done, PR #3817 |
| TEST.457 | test_no_todo_codes_shown fails on main: util/settings.py shows the item code PERF.154 to users (bug) | none | done, PR #3761 |
| TEST.458 | test_sector_map_scale_line_follows_the_zoom fails on main (bug) | none | done, PR #3771 |
| TEST.459 | The point_phenomena links browser test is red on clean main (bug) | none | done, already green on main (repaired by MAP.349 and MAP.337 step 2) |
| TEST.460 | test_sweep_life has 5 setup errors on main: KeyError 'h' because tuning.PLANET_CLASSES has no class h since the class cut-over (bug) | none | dropped (Bug 2: the KeyError was fixed in PR 3747 at 6:00 PM EDT as the lane reported; a gate on an older base shows it, a rebase clears it) |
| TEST.461 | test_db_column_round_trip fails on planets.porosity: the column is never written (bug) | none | open |
| TEST.462 | test_planets_respect_their_stars_age_and_history fails in a batch run: a class E planet in the habitable list around a star under 0.1 Gy (bug) | none | done, PR #3794 |
| TEST.463 | Reds on main reported by Core 5 at 10:33 PM EDT 2026-10-29: triage each one and fix or assign it (bug) | none | done, PR #4307 |
| TEST.464 | Cache band_fractions between calls so test_gen_component_pass and the whole-sphere census fixtures finish inside the gate (bug) | none | done, PR #4075 |
| TEST.465 | test_admin_script_cli: query_db exits 1 instead of 2 for an incomplete command line when the database is not reachable (bug) | none | done, PR #3932 |
| TEST.466 | test_generate_py_system_output_writes_a_file_and_no_database: the plain system output starts with a newline (bug) | none | done, PR #3932 |
| TEST.467 | Two tests red on the base branch seen by Star 1: the wide-binary planet-class check and the planet serialization allowlist (formation_distance_au, host_eccentricity) (bug) | none | dropped (Bug 2 rechecked on main c4b6bd167: the serialization fuzz and generation fuzz files pass and every wide-binary test passes; Star 1 named no failing test. Refile with the test name if it recurs after rebasing) |
| TEST.468 | test_work_units::test_the_batch_fills_are_units_of_id_ranges_recorded_as_their_commands and test_generation_stats::test_an_unmeasured_server_estimates_from_the_defaults are red on a clean main (bug) | none | done, PR #4425 |
| TEST.469 | test_interior_heat::test_a_generated_body_has_a_urey_ratio_above_one_only_while_its_mantle_is_still_warming fails on main: a generated body has urey_ratio 1.85 with a mantle that is no longer warming (bug) | none | done, already green on main |
| TEST.470 | test_orbit_elements fails on main: the moon-light test (relative error 6e-6) and the ellipse tests (a lunar-range problem) (bug) | none | done, already fixed on main (58714c8de, 80543090d) |
| TEST.471 | tests/test_near.py::test_the_command_line_asks_for_the_nearest_few fails: argparse exits with code 2 so the pytest.raises(SystemExit, match="--distance") never sees the text (bug) | none | done, PR #4019 |
| TEST.472 | Four migration tests are red on clean main: aerosols same-seed draws, albedo_greenhouse saved-system temperature, and two gas_tables tests (AsteroidBelt has no attribute moons) (bug) | none | done, PR #4421 |
| TEST.473 | test_star_census_raises_no_new_flag is red on main: G stars 2.90 percent against the 3.13 to 6.99 percent range, and the A0V companion rate 0.7456 against up to 0.7367 (bug) | none | done, PR #4031 |
| TEST.474 | tests/test_web_a11y.py::_add_bright_stars fails with pymysql 1136 "Column count doesn't match value count" and errors every site_app browser test in setup (bug) | none | done, PR #4023 |
| TEST.475 | tests/test_web_nearby.py: four of eight tests error in fixture setup at 900 s because the client fixture generates two sectors per test (bug) | none | done, PR #4049 |
| TEST.476 | A test that measures how often items are redrawn or regenerated to make them fit: stars in a tree, planets, moons and belts | none | done, PR #4199 |
| TEST.477 | The overlap verifier: audit() of a generated system and a seed-sweep test over forced triples, quadruples and ordinary systems (GEN.1006 P2) | none | done, PR #4071 |
| TEST.478 | tests/test_api.py _bright_row (and its row near line 1491) is one value short since GEN.846 added bright_stars.rotation_period_hours (bug) | none | done, PR #4056 |
| TEST.479 | test_radiative_convective_real::test_earth_sends_out_what_it_absorbs fails on the train 4 branch: the Earth tropopause temperature is 168.33 K against the test bound of 170 to 230 K (bug) | none | open |
| TEST.480 | test_gen_builder_internals::test_slot_planet_honours_class_and_moon_count gets class J where it forced class I at a cold slot (bug) | none | done, already green on main (fixed by the class epochs) |
| TEST.481 | Three tests red on clean main: the Hill sphere moon drop, the halo evolutionary flux counts and the large component phenomena unit (bug) | none | done, PR #4096 |
| TEST.482 | test_stages::test_plan_is_the_skeleton_and_each_scatter_and_fill_is_one_stage_of_its_own is red on main: fill-map-luminosity lists the block-lights stage (bug) | none | done, PR #4093 |
| TEST.483 | test_disregard_missing_data (6 tests) and test_rq_time_limits::test_the_generate_pages_job_gets_a_limit_from_the_time_its_steps_are_recorded_to_take (3 params) are red on clean main (bug) | none | done, PR #4415 |
| TEST.484 | Four tests red on clean main: climate draws [c], habitability world scored consistently, life sweep class D "supports 3b, drew no life", and the moon-helper Hill sphere test (bug) | none | done, PR #4120 |
| TEST.485 | test_planets h-0: atmospheric_pressure 0.99999999993 sits under the 1.0 floor (bug) | none | done, PR #4136 |
| TEST.486 | test_planets c-2: a life class is kept but life_class_scores is None when no class scores above the line after the oxygen build (bug) | none | done, PR #4321 |
| TEST.487 | test_no_todo_codes_shown fails twice: a GEN.474 string and the PERF.203 help text show TODO codes to users (bug) | none | done, PR #4140 |
| TEST.488 | test_migration_albedo_greenhouse is red on clean main (bug) | none | dropped (duplicate of TEST.472 (albedo_greenhouse migration is named there)) |
| TEST.489 | Reds on clean main found by Core 3: test_cache_warmup (nodes_to_warm, full_cache_leaves, nodes_family: insert_bright_stars column count mismatch), test_satellites (2), test_activity (stored field and exposure), test_sector_gen research rates per star (black-hole rate 0.00123 over the 0.0012 cap) (bug) | none | done, PR #4230 |
| TEST.490 | Oddity plausibility harness (src/tests/oddity_plausibility_cli.py): force each oddity on 200 bodies with a fixed seed and assert the general invariants (energy, structure, stability, air, class derived, no NaN) and each row's own | none | done, PR #4147 |
| TEST.491 | W21 audit: check end to end that M dwarf water loss can leave an oxygen air and that the life scores never read oxygen alone (Luger and Barnes 2015, Wordsworth and Pierrehumbert 2014) | none | done, PR #4135 |
| TEST.492 | Census re-baseline for every oddity share and the prevalence suite rows (one pass after each of trains 12, 13 and 14) | none | open |
| TEST.493 | Giant oddity audit: the super-Jupiter mass ceiling (G9), the hot Neptune desert survivors (G7, LTT 9779 b), and which G rows are derived outputs only (G2 to G6, G8, G10, G11); file bugs for any that are not | none | done, PR #4138 |
| TEST.494 | Ring invariants in the plausibility harness: ring mass under 1e-4 of the planet, ring inside the Roche limit unless R4 names its resonance, optical depth an output, no ring inside 0.1 AU, ring plane in the equator | none | done, PR #4147 |
| TEST.495 | Designation for captured bodies: the belt-member designation `<belt ID>.<rank hex>` is kept with a parent link, outside-reservoir bodies print `<planet ID>.<rank hex>`; both pure and round-trip through parse_id | none | done, PR #4163 |
| TEST.496 | test_the_wheel_carries_the_map_from_the_whole_galaxy_into_a_system fails on main in its wheel-out leg (the URL keeps sector= after wheeling out) (bug) | none | done, PR #4444 |
| TEST.497 | test_planet_physical_properties e-3 fails on an unseeded star (air pressure 0.99999 Pa against a floor of 1.0), the same cause as TEST.486 (bug) | none | done, PR #4165 |
| TEST.498 | A rogue planet is functionally and programmatically identical to a bound planet: the same columns, the same moon and composition generation, the same classes, differing only in the rogue note and having no host star | none | open |
| TEST.499 | RP9: rogue yield calibration against microlensing and the terrestrial scale: rogues per star against Sumi 2023 (central 21, lower 1 sigma 8 over the probed range), the Mroz 2017 Jupiter limit and Barclay 2017's 2.5 terrestrial bodies per star | none | open |
| TEST.500 | test_moon_composition::test_a_disk_moon_gains_ice_past_its_disks_ice_line fails on a clean main (reported by Aux 4 on PR #4184) (bug) | none | done, PR #4230 |
| TEST.501 | Order-independence test: chart region A then B and B then A and compare the union of real-donor objects and each sector's made-here digest; chart a region twice and nothing changes; the top-up checked as "same stored state in, same objects out" | none | open |
| TEST.502 | test_orbital_motion::test_moon_having_planet_reflex_offset_matches_manual_formula fails on a clean main: settle_orbit moves captured moons without refreshing the planet offset (bug) | none | done, PR #4385 |
| TEST.503 | Four tests are red on a clean main: test_db_column_round_trip::test_every_column_is_written (moons.magnetopause_rp is NULL in every row), the test_object_rows star read-back (columns state and decay_age_yr), test_gen_hosting_wiring directed moon count, and test_reproducible_draws hypot/BLAS (bug) | none | done, PR #4385 |
| TEST.504 | test_runaway_history::test_the_surface_water_a_planet_keeps_is_what_its_class_gave_less_what_it_lost fails on main (the M3V case sees no loss) (bug) | none | done, PR #4283 |
| TEST.505 | test_greenhouse_table agreement with the model fails for the Venus and 5 bar CO2 airs (2 cases) on main, with and without the detmath swap (bug) | none | dropped (duplicate of GEN.1032) |
| TEST.506 | test_no_todo_codes_shown fails on main: GEN.1128 shows in the --help text of planetgen/cli/generate.py (line 696) and of the sectors and galaxy commands (bug) | none | done, PR #4333 |
| TEST.507 | Add Canes Venatici I (and decide on Canis Major) to the reference satellite table so it holds all 12 Milky Way satellites above 1e5 solar masses | none | done, PR #4387 |
| TEST.508 | Leo I fixture row in test_satellites has an age of 5 Gy so it forms gas; change the age or label it a stand-in | none | done, PR #4371 |
| TEST.509 | The climate sweep twin sits at 313 K, far above the 288 K of the Earth it stands for: make it a real Earth twin so the sweep tests the feedbacks where they act | none | done, PR #4376 |
| TEST.510 | test_no_new_hand_written_mass_conversion is red on main: new hand-written multiplies or divides by a mass constant in sheet.py, member_record, rogue, belt_mass, binding_energy, planet_architecture and tuning | none | done, PR #4403 |
| TEST.511 | test_every_column_is_on_the_sheet_or_says_why_it_is_not is red on main for belt and asteroid_field: chondrite_group (GEN.475) is neither on the sheet nor excused | none | done, PR #4410 |
| TEST.512 | test_rogue_planet_parity::test_the_gap_between_a_rogue_and_a_bound_planet_never_grows is red on main | none | done, PR #4405 |
| TEST.513 | test_strings_the_code_shows_have_no_todo_code is red on main: TODO codes in shown strings (sheet.py GEN.1142, small_body_composition.py GEN.475/GEN.473, surface_composition.py GEN.474 twice) | none | done, PR #4412 |
| TEST.514 | test_habitability_world::test_the_score_round_trips_and_a_cloud_rescores_the_dose is red on main: phi4_temperature 0.75 against 0.8075 after the cloud refresh (bug) | none | done, PR #4437 |
| TEST.515 | test_planet_zones::test_the_share_of_pair_nodes_with_a_p_type_window_is_the_measured_one is red on main for 2 and 7 stars after the dark members (GEN.1013) changed the tree draw (bug) | none | done, already fixed by PR #4442 (TEST.515 re-measure); not reproduced on main |
| TEST.516 | tests/test_readme_external_files.py has 3 failures on main: the README external-files section differs from the dataset catalogue (bug) | none | done, PR #4457 |
| TEST.517 | test_class_pages_show_the_constants is red on main: the brown dwarf mass range text on the class page differs from the test (bug) | none | done, PR #4458 |
| TEST.518 | test_a_body_with_a_stored_day_turns_at_the_shared_rate_and_a_tilted_one_turns_backwards fails in a whole-file run of test_web_browser_bodyview.py but passes alone (timing sensitive) (bug) | none | done, PR #4465 |
| TEST.519 | tests/test_web_a11y.py::test_page_is_accessible_and_fits[web.galaxy-*] fails the axe rule scrollable-region-focusable on .map-key-lines (all four size and theme runs) on clean main (bug) | none | done, PR #4531 |
| TEST.520 | Extend the wiring guard to stored columns and shadow functions, and wire, show or drop every column only the Details sheet reads | none | open |
| TEST.521 | Wire, register or delete the 27 pure functions nothing calls in 13 modules (Core 1's share of the TEST.520 shadow functions) | none | done, PR #4545 |
| TEST.522 | Wire, register or delete the 29 pure functions nothing calls in 21 modules (Core 3's share of the TEST.520 shadow functions) | none | done, PR #4528 |
| TEST.523 | Wire, register or delete the 17 pure functions nothing calls in 12 modules (Star 2's share of the TEST.520 shadow functions) | none | done, PR #4519 |
| TEST.524 | Wire, register or delete the 38 pure functions nothing calls in 21 modules (Aux 4's share of the TEST.520 shadow functions) | none | done, PR #4537 |
| TEST.525 | Wire, register or delete the 26 pure functions nothing calls in 14 modules (Aux 3's share of the TEST.520 shadow functions) | none | done, PR #4533 |
| TEST.526 | Wire, register or delete the 8 pure functions nothing calls in 4 modules (Tides 1's share of the TEST.520 shadow functions) | none | done, PR #4613 |
| TEST.527 | Wire, register or delete the 12 pure functions nothing calls in 5 modules (Bug 1's share of the TEST.520 shadow functions) | none | done, PR #4546 |
| TEST.528 | tests/test_systemmap.py::test_a_planet_and_its_moon_glow_in_the_colour_of_their_air_under_their_own_star is red on main (bug) | none | done, PR #4520 |
| TEST.529 | test_the_job_strip_shows_the_running_job_on_every_admin_page fails on clean main (bug) | none | done, PR #4556 |
| TEST.530 | Three tests are red on clean main: two in test_roche_interaction (post-common-envelope classifier) and the tuning pair-order check (REFERENCE_BLACK_HOLE_DRAWS 0.55 against 0.5) (bug) | none | done, PR #4535 |
| TEST.531 | Delete the sector-settle queue chain that no route reaches any more (api_jobs.settle_sectors and the test helper queue_settle) | none | done, PR #4563 |
| TEST.532 | Give the bright-star backfill unit path (backfill_bright_stars_unit) its own direct tests, since the serial form the tests cover is no longer the path a run takes | none | done, PR #4566 |
| TEST.533 | test_disk_history::test_the_present_rate_over_the_mean_is_within_a_factor_of_two_of_the_milky_ways is red on main (0.505 against a bound of 0.393, 0.554 after GEN.999): the model's recent star formation rate over its mean is about 2.5 times the Milky Way's (bug) | none | done, PR #4575 |
| TEST.534 | test_near::test_the_search_counts_what_is_not_generated_yet_and_never_reads_the_sectors_to_search is red on clean main: its _scatter_rows insert fails with column count does not match value count (bug) | none | done, PR #4580 |
| TEST.535 | The two test_planet_zones P-type share tests are red on clean main (dark members) (bug) | none | done, already fixed by PR #4442 (TEST.515 re-measure); not reproduced on main |
| TEST.536 | test_deploy_scripts::test_scripts_only_name_repo_paths_that_exist is red on main (bug) | none | done, PR #4578 |
| TEST.537 | The stored-trees fill test counted dark members as star rows (red on main since GEN.1013) (bug) | none | done, PR #4592 |
| TEST.538 | test_atmosphere_standin_registry fails on main: stand-in rows AM28 and AM29 name a test that no longer exists (bug) | none | done, PR #4609 |
| TEST.539 | test_accreting_compact_binaries: the black-widow and redback spider cases fail on main (post-common-envelope against detached) (bug) | none | done, PR #4616 |
| TEST.540 | test_orbit_elements: an eccentric planet is as warm as a circular one, and a moon gets the light of its parent's orbit, both fail on main (bug) | none | done, PR #4618 |
| TEST.541 | Four star tests fail on clean main: the population model, the hot white dwarf classes O5VII and B0VII, and the forced subdwarf (bug) | none | done, PR #4704 |
| TEST.542 | The life sweep expects 20 class letters and the rules give 19, since GEN.1135 retired Class Q (bug) | none | done, PR #4677 |
| TEST.543 | The atmosphere sweep audit reads planetgen.physics.ozone.STRUCTURE_RANGE, which is now STRUCTURE_KNOTS, so three sweep tests fail (bug) | none | done, PR #4699 |
| TEST.544 | The class sweep finds five rocky bodies classed T, a gas class, so the totality test fails (bug) | none | done, PR #4677 |
| TEST.545 | The table sweep raises temperature must be positive for small-body-taxa/asteroid-field at 0 K, and two rows are never reached (bug) | none | done, PR #4714 |
| TEST.546 | The object timing test finds a unit saved with the kind star where it expects system (bug) | none | open |
| TEST.547 | Dry-class carbon bracket test is red on main after GEN.1005 (bug) | none | done, PR #4701 |
| TEST.548 | Two phenomenon-scatter count tests fail on main: remnant counts moved past 1.5 times the model (bug) | none | open |
| TEST.549 | test_star_population passes the removed bounds_only argument to _bright_table and is red on main (bug) | none | done, PR #4704 |
| TEST.550 | Five census tests stay red on main until the census re-tag for epochs 129 and 130 (bug) | none | open |
| TEST.551 | test_burn_down_register::BR20 is red on main: row 6.95:BR20 has a line but is no longer L3 or L4 (bug) | none | done, PR #4710 |
| TEST.552 | test_plan_outline expects the outline crest at 3.32 times and the density model now gives 3.18 (bug) | none | open |
| USR.1.1 | USR.2 |
| USR.1.2 | USR.3 |
| USR.1.3 | USR.4 |
| USR.1.4 | USR.5 |
| USR.1.5 | USR.6 |
| USR.1.6 | USR.7 |
| UX.0.1 | UX.15 |
| UX.0.2 | UX.16 |
| UX.4.1 | UX.17 |
| UX.4.2 | UX.18 |
| VIEW.1.1 | VIEW.2 |
| VIEW.1.2 | VIEW.3 |
| VIEW.1.3 | VIEW.4 |

## Where old numbers are cited

Each citation below resolves to one item by its date.

### CHANGELOG.md

| Release | Text | New ID |
|---|---|---|
| 7.6.0 | "`docs/TODO.md` item 54 tracks the admin Generate page's POSIX-only job" | OPS.4 (written on the PR #115 branch at 20:48Z, where the item was 54; it became 55 on merge) |
| 7.6.1 | "TODO item 45" | GEN.20 |
| 7.9.1 | "TODO item 44" | GEN.19 |
| 7.13.1 | "TODO item 42" | GEN.17 |
| 7.16.1 | "TODO item 43" | GEN.18 |
| 7.18.1 | "TODO items 40 and 41" | GEN.15, GEN.16 |
| 5.46.34, 5.46.35 | describe the 2026-09-24 compaction and first numbering | (no item) |
| 7.49.0 | "Design: docs/design/habitability/population-and-politics.md (TODO 51-54)" | POP.1 to POP.4 |

### Commit and PR titles

| Citation | Date (UTC) | New ID |
|---|---|---|
| "item 1's cube-tile fetching" (commit 5a89ff4) | 2026-09-24 02:02Z | OPS.2 |
| "new item 8", "item 8" (commits 822f0ed, c26685a, 7605b86) | 2026-09-30 16:49Z to 17:59Z | MAP.36 |
| "item 14" (commit 5e3f9d4) | 2026-09-30 16:51Z | UX.10 |
| "TODO item 61" (commit d292b45) | 2026-09-30 20:03Z | UX.12 |
| "Windows jobs item becomes 55" (commit 6851fec) | 2026-09-30 20:48Z | OPS.4 |
| TODO 45 (PR #116) | 2026-09-30 21:29Z | GEN.20 |
| TODO 39b (PR #118), 39c (commit 51b3543, PR #119), 39a (PR #131) | 2026-09-30 21:37Z to 22:54Z | SEC.18, SEC.19, SEC.17 |
| TODO item 1 (PR #122) | 2026-09-30 21:38Z | UX.6 |
| TODO 20, TODO 14 (commits in PR #120) | 2026-09-30 21:40Z, 21:45Z | MAP.41, MAP.35 |
| TODO 34 (PR #121) | 2026-09-30 21:43Z | NAV.2 |
| TODO 44 (PR #123) | 2026-09-30 21:33Z | GEN.19 |
| TODO 55 (commit 2bf3ef3, PR #124) | 2026-09-30 21:43Z | OPS.4 |
| TODO 12 (PR #129), "local item 12" (commit 3d83a1b) | 2026-09-30 21:53Z, 22:18Z | MAP.33 |
| TODO 33 (PR #130) | 2026-09-30 21:47Z | NAV.1 |
| TODO 42 (PR #128) | 2026-09-30 21:39Z | GEN.17 |
| TODO 50 (commit 5d776ce, PR #125) | 2026-09-30 22:15Z | OPS.3 |
| TODO 2, 3, 48 (PR #126) | 2026-09-30 22:42Z | UX.7, UX.8, UX.12 |
| TODO 43 (PR #133) | 2026-09-30 21:43Z | GEN.18 |
| TODO 7 (PR #132) | 2026-09-30 23:07Z | GEN.3 |
| TODO 40, 41 (PR #137) | 2026-09-30 21:48Z | GEN.15, GEN.16 |
| TODO 13 (PR #134) | 2026-09-30 23:21Z | MAP.34 |
| TODO 5 and 6 (PR #136) | 2026-09-30 23:25Z | GEN.1, GEN.2 |
| TODO 28 and 31 (PR #138) | 2026-09-30 23:45Z | GEN.11, GEN.14 |
| TODO 29 (PR #140) | 2026-09-30 23:48Z | GEN.12 |
| TODO 15, 16 (PR #142) | 2026-10-01 00:34Z | MAP.36, MAP.38 |
| TODO 27 (PR #144, commits b2e4da3 and fa7655f, PR #153) | 2026-10-01 01:04Z to 02:35Z | GEN.10 |
| "add TODO 56" (commit e2f6d57) | 2026-10-01 01:19Z | UX.1 |
| items 57-61 (PR #146) | 2026-10-01 01:21Z | ADM.5 to ADM.8, SEC.1 |
| TODO 30 (PR #145) | 2026-10-01 01:29Z | GEN.13 |
| TODO 17 (PR #147) | 2026-10-01 01:48Z | MAP.39 |
| TODO 26 (PR #149, commit 92de010) | 2026-10-01 01:57Z | UX.18 (storage part) |
| items 62-63 (PR #150) | 2026-10-01 01:58Z | UX.2, MAP.3 |
| TODO 24 (PR #151) | 2026-10-01 02:10Z | ADM.3 |
| TODO 35 (PR #152) | 2026-10-01 02:22Z | DB.1 |
| items 64-69 (PR #154) | 2026-10-01 02:25Z | USR.2 to USR.7 |
| TODO 64-73 (PR #155, commit a18bc46) | 2026-10-01 02:24Z to 02:26Z | drill-down: 64 MAP.28, 65 MAP.29, 66 MAP.16, 67 MAP.20, 68 MAP.21, 69 MAP.22, 70 MAP.23, 71 MAP.24, 72 MAP.25, 73 MAP.27 (note 6) |
| "renumber the drill-down items to 70-79" (PR #156) | 2026-10-01 02:31Z | MAP.28, MAP.29, MAP.16 to MAP.27 |
| TODO 18 (PR #153, commit e9191ad) | 2026-10-01 02:04Z | MAP.40 |
| TODO 32 (PR #157) | 2026-10-01 02:35Z | GEN.6 |
| TODO 56 (PR #159) | 2026-10-01 02:41Z | GEN.22 (uncertain, note 2) |
| items 80-85 (PR #158), 81-82, 83-85, 80 (commits) | 2026-10-01 02:41Z to 02:56Z | DOC.1, DOC.2, DOC.3, VIEW.2 to VIEW.4 |
| "shipped items 33 and 50" (PR #162) | 2026-10-01 03:00Z | NAV.1, OPS.3 (note 4) |
| item 86 (PR #163) | 2026-10-01 03:16Z | PERF.3 |
| items 87-92 (PR #164), 89 (commit fac0da9) | 2026-10-01 03:26Z to 03:49Z | UX.3, PERF.4 to PERF.8 |
| TODO 70, 71 (PR #160) | 2026-10-01 03:58Z | MAP.28, MAP.29 |
| TODO 19 (PR #168) | 2026-10-01 04:16Z | MAP.1 (note 10) |
| TODO 51-54 (commit f352b1e, design doc; commit 084d875, PR #169) | 2026-10-01 04:09Z to 04:37Z | POP.1 to POP.4 |
| TODO 56 (commit 6199719, PR #167) | 2026-10-01 01:52Z | UX.1 (note 2) |
| TODO 26 (commit 76cff6f, PR #167) | 2026-10-01 02:05Z | UX.18 (pages) |
| TODO 36 (commit 5158d7e, PR #167) | 2026-10-01 02:42Z | UX.5 |
| "Delete shipped TODO 49" (commit 9c47d5d, PR #167) | 2026-10-01 04:06Z | MAP.4 (note 3) |
| TODO 25 (commit abd6965, PR #178) | 2026-10-01 04:15Z | UX.17 |
| TODO 37, 38 (PR #170) | 2026-10-01 04:24Z | API.1, API.2 |
| TODO 63 (commit 0889d91, PR #178) | 2026-10-01 04:21Z | MAP.3 |
| TODO 74 (commit 292540d, PR #178) | 2026-10-01 04:24Z | MAP.21 |
| TODO 8 (commit 75448bf, PR #178) | 2026-10-01 04:30Z | PERF.2 |
| TODO 72 (PR #171) | 2026-10-01 04:33Z | MAP.16 |
| TODO 77 (PR #172) | 2026-10-01 04:41Z | MAP.24 |
| item 93 (commit dcbf0fb, PR #175) | 2026-10-01 04:50Z | PERF.9 |
| TODO 79 (PR #176) | 2026-10-01 04:51Z | MAP.27 |
| TODO 73, map side (PR #177) | 2026-10-01 04:57Z | MAP.20 (part) |
| items 94, 95 (commit 907e4e4, PR #175) | 2026-10-01 04:58Z | PERF.10, PERF.11 |
| items 96-99 (merge commit 8767450), TODO 99 (commit 8bac5cb), PR #175 | 2026-10-01 05:15Z | MAP.43, MAP.47, MAP.48, MAP.37 (note 13) |
| item 100 (commit 83f1c98, PR #175) | 2026-10-01 05:22Z | MAP.17 |
| item 101 (commit bb62014, PR #175) | 2026-10-01 05:24Z | MAP.26 |
| items 102, 103 (commit b33b5b1, PR #175) | 2026-10-01 05:27Z | UX.13, UX.14 |

### Docs and code comments

| Where | Text | Date (UTC) | New ID |
|---|---|---|---|
| [navigation-frames.md](../design/orbits/navigation-frames.md) | "TODO items 33 and 34" (first written as "items 28 and 29" at 18:14Z) | 2026-09-30 18:39Z | NAV.1, NAV.2 |
| [interstellar-object-rates.md](../design/systems/interstellar-object-rates.md) | "TODO items 5 (rates) and 6", "TODO item 27" | 2026-09-30 20:27Z | GEN.1, GEN.2, GEN.10 |
| [nebula-and-asteroid-field-classes.md](../design/systems/nebula-and-asteroid-field-classes.md) | "items 27-31", "item 28" ... "item 31" | 2026-09-30 18:39Z | GEN.10 to GEN.14 |
| [galaxy-drilldown-navigation.md](../design/interface/galaxy-drilldown-navigation.md) | "items 64-73", then "items 70-79" and the build-order table; "TODO 19"; "TODO 63" | 2026-10-01 02:24Z, 02:27Z | MAP.16 to MAP.29; MAP.1; MAP.3 |
| Windows and macOS hosting guides | "item 50"; "item 54" | 2026-09-30 20:48Z to 20:50Z | OPS.3; OPS.4 |
| Installer scripts | "See docs/TODO.md item 50" | 2026-09-30 20:43Z | OPS.3 |
| `program_constants.py`, `generate.py` | "TODO item 28", "TODO item 31", "TODO item 27", "docs/TODO.md item 5" | 2026-09-30 18:14Z to 2026-10-01 01:04Z | GEN.11, GEN.14, GEN.10, GEN.1 |
| `utils.py`, `test_distance_format.py` | "docs/TODO.md item 1" | 2026-09-30 18:14Z to 21:38Z | UX.6 |
| `test_web_generate.py` | "docs/TODO.md item 39" | 2026-09-30 21:37Z | SEC.16 (SEC.18) |
| System page code and `test_web_system_phen.py` | "TODO items 2 and 48", "TODO 48" | 2026-09-30 22:42Z | UX.7, UX.12 |
| Galactic motion code and tests | "TODO item 32" | 2026-10-01 02:35Z | GEN.6 |
| `generate.py`, `cometData.py` | "TODO item 27", "TODO item 30" | 2026-10-01 01:04Z, 01:29Z | GEN.10, GEN.13 |
| `test_local_time.py` | "docs/TODO.md item 22" | 2026-10-01 00:07Z | UX.10 |
| `test_star_type_rendering.py` | "TODO 55" | 2026-10-01 00:08Z | MAP.13 (note 5) |
| `test_bright_star_sampling.py` | "TODO items 55 and 56" | 2026-09-30 23:46Z | GEN.21, GEN.22 (note 2) |
| [population-and-politics.md](../design/habitability/population-and-politics.md), `population.py`, `api/population.py`, `_db.py`, `schema.sql`, `program_constants.py`, `generate.py`, `test_population.py`, `docs/reference/database-schema.md` | "TODO 51-54", "TODO items 51-54" | 2026-10-01 04:09Z to 05:28Z | POP.1 to POP.4 |
| `lib/pagecache.py`, `web/__init__.py`, `test_page_cache.py` | "TODO 8" | 2026-10-01 04:30Z | PERF.2 |
| `lib/phenomenonrender.py`, `test_phenomenon_views.py` | "TODO 25" | 2026-10-01 04:15Z | UX.17 |
| `static/style.css` | "TODO 63" | 2026-10-01 04:21Z | MAP.3 |
| `test_web_sector_nav.py` | "TODO 74" | 2026-10-01 04:24Z | MAP.21 |
| `test_web_facilities.py` | "TODO 36" | 2026-10-01 02:42Z | UX.5 |

### Code tags `TODO(<area> #N)`

Tags were renumbered along with the file, so read them by the date of the
commit that wrote them. Only the last row is still in the code.

| Written (UTC) | Tags | New IDs |
|---|---|---|
| 2026-09-30 16:44Z | galaxy-map #3 to #12 | #3 MAP.31, #4 MAP.32, #5 MAP.33, #6 MAP.34, #7 MAP.35, #8 MAP.38, #9 MAP.39, #10 MAP.40, #11 MAP.1, #12 MAP.41 |
| 2026-09-30 16:49Z to 16:58Z | galaxy-map #8 to #13 (#3 to #7 unchanged) | #8 MAP.36, #9 MAP.38, #10 MAP.39, #11 MAP.40, #12 MAP.1, #13 MAP.41 |
| 2026-09-30 18:14Z | galaxy-map #10 to #21 | #10 MAP.31, #11 MAP.32, #12 MAP.33, #13 MAP.34, #14 MAP.35, #15 MAP.36, #16 MAP.38, #17 MAP.39, #18 MAP.40, #19 MAP.1, #20 MAP.41, #21 MAP.42 |
| 2026-09-30 18:14Z | distances #1, system-list #2, site-header #3, search #4, phenomena #5 to #7, sector-map #23, #24, phenomena #25, #26 | UX.6, UX.7, UX.8, UX.9, GEN.1 to GEN.3, ADM.2, ADM.3, UX.17, UX.18 |
| 2026-09-30 18:14Z to 18:39Z | orbits #27, nav #28, nav #29, facilities #30, facilities #31 | GEN.6, NAV.1, NAV.2, DB.1, UX.5 |
| 2026-09-30 18:39Z | phenomena #27 to #31, orbits #32, nav #33, nav #34, facilities #35, facilities #36 | GEN.10 to GEN.14, GEN.6, NAV.1, NAV.2, DB.1, UX.5 |
| 2026-09-30 19:02Z to 20:07Z | security #39 to #52; physics #53 to #58 | SEC.3 to SEC.16; GEN.15 to GEN.20 |
| 2026-09-30 20:01Z to 20:08Z | web-pages #59 to #62 | UX.11, MAP.12, UX.12, MAP.4 |
| 2026-09-30 20:07Z | physics #40 to #45 | GEN.15 to GEN.20 |
| 2026-09-30 20:08Z | web-pages #46 to #49 | UX.11, MAP.12, UX.12, MAP.4 |
| 2026-09-30 20:43Z, 20:48Z | installers #50; windows #54, then #55 | OPS.3; OPS.4 |
| 2026-10-01 00:07Z | sector-map #55 | MAP.13 |
| Still in the code at 7.40.0 | facilities #36, galaxy-map #19, installers #50 (a test asserting it is gone), phenomena #25, phenomena #26, web-pages #56 | UX.5, MAP.1, OPS.3, UX.17, UX.18, UX.1 |
| Still in the code at 7.58.2 | installers #50 only (the test asserting it is gone); the others left with their items | OPS.3 |

## Notes and ambiguities

1. **2026-09-24, two numberings at once.** From 01:57Z to 02:27Z the
   branch for PR #73 numbered its copy 1 OOM, 2 skeleton sphere, 3 cache,
   4 Measure path, 5 phenomena on the Sector Map, 6 Galaxy Map rework,
   7 and 8 the API items, 9 to 12 Population, while `main` still had the
   first numbering (1 OOM, 2 sphere marker, 3 orbit lines, 4 star glow,
   5 skeleton sphere, 6 cache ...). Both rows are in the table. The only
   citation from that window ("item 1") means OPS.2 in both.
2. **56 for pre-placed bright stars (uncertain).** PR #159 ("TODO 56"),
   its PR text ("TODO 56 and the population-ages item are removed") and
   `test_bright_star_sampling.py` ("TODO items 55 and 56") use 55 and 56
   for the star population work. `docs/TODO.md` never had a 56 for it:
   only 55 ("Use population ages at sector fill") was written there. The
   56 most likely came from the working thread's own plan. It is mapped
   to GEN.22. At the same time, another branch wrote 56 as "Class
   reference pages" (2026-10-01 01:19Z, on `main` from 03:54Z), which is
   UX.1; that branch's commit e2f6d57 ("add TODO 56") means UX.1.
3. **49, System Map names.** Shipped in 7.21.1 (PR #135, deleted at
   2026-10-01 00:07Z), then the item came back in PR #159's merge
   (02:41Z), almost certainly by accident. It was still open in TODO.md
   at 7.40.0. PR #167 deleted it as "shipped in #135" (commit 9c47d5d,
   on `main` at 2026-10-01 04:37Z), so MAP.4 is recorded as done in
   7.21.1.
4. **33 and 50 came back.** NAV.1 (33) shipped in 7.14.0 and OPS.3 (50)
   in 7.16.0, but stale merges in PR #157 (33) and PR #159 (50) put the
   items back, and PR #162 removed them again at 2026-10-01 03:00Z.
   Citations of 33 and 50 between 22:51Z on 2026-09-30 and 03:00Z on
   2026-10-01 still mean NAV.1 and OPS.3.
5. **55, three items.** 55 meant OPS.4 (Windows jobs, 2026-09-30 20:48Z
   to 22:14Z), then on 2026-10-01 both MAP.13 (star dots, from PR #135
   at 00:07Z until PR #143 at 01:00Z) and GEN.21 (population ages, from
   PR #141 at 00:13Z until PR #159 at 02:41Z). Between 00:13Z and 01:00Z
   `main` had two items numbered 55. Tell them apart by content: star
   dots or rendering means MAP.13; star population, ages or S7 means
   GEN.21. Earlier, during the security numbering (19:02Z to 20:07Z),
   55 was GEN.17.
6. **64 to 73 collision.** PR #154 put user accounts at 64 to 69 at
   02:25Z; PR #155 put the drill-down at 64 to 73 at 02:26Z, so `main`
   briefly had two 64s to 69s. PR #156 moved the drill-down to 70 to 79
   at 02:31Z. "64-73" in PR #155's title and in the design doc's first
   version means the drill-down; "64-69" in PR #154 means user accounts.
7. **The security numbering (2026-09-30 19:02Z to 20:23Z on `main`).**
   For about an hour, 39 to 51 were individual security findings, 52 was
   Hardening, 53 to 58 the generation bugs and 59 to 62 (then 63 to 66)
   Population. PR #112 fixed the findings and renumbered. Each finding
   has its own ID under SEC.2 so a citation from that hour still
   resolves; none is known outside the code tags.
8. **39a, 39b, 39c.** These letters were never in TODO.md. They are the
   three parts of item 39 (Hardening) as named in commit and PR titles:
   a = per-username login backoff, b = upper bounds on admin generation
   inputs, c = hashed lock file. The fourth part, HSTS, shipped with the
   security fixes in 7.5.0.
9. **2026-09-30 18:14Z to 18:39Z.** A numbering seen only on the branch
   for PR #109 (commit 885a21c), replaced 25 minutes later before it
   merged. In it, 27 was the correlative update, 28 courses, 29 warp
   speeds, 30 starbases, 31 facilities on the web, 32 and 33 the API
   items, 34 to 37 Population. Code tags and one design doc were written
   with these numbers.
10. **After the first snapshot (7.40.0).** PR #168 (04:16Z) finished
   TODO 19 (MAP.1), and PR #167 (04:37Z) finished 56 (UX.1), 26
   (UX.18) and 36 (UX.5) and deleted 49 (MAP.4). Both are now recorded
   as done.
11. **Item 4 and item 24 shipped in two parts.** UX.9 (tag search) was
   done by 7.22.1 and 7.29.0; ADM.3 (Generate buttons) by 7.27.0 for the
   Sector Map and 7.34.0 for the Galaxy Map, after its title changed to
   "Generate buttons on the Galaxy Map's unfilled sectors".
12. **The Galaxy Map rework (MAP.5).** The 2026-09-24 item "Rework the
   Galaxy Map" was not finished as such: on 2026-09-30 it was replaced by
   the Galaxy Map plan (MAP.31 to MAP.42, and MAP.1 for the
   follow-ups), so those are its subitems.

13. **93 to 103 (2026-10-01 04:50Z to 05:27Z).** All written on the
   PR #175 branch and merged to `main` at 05:27Z. 96 to 99 came in a
   merge commit (8767450, "Merge main; add Galaxy Map bugs"), not a
   plain commit, so `git log -S` on the item text misses them. They
   were never renumbered.
14. **73 is half done.** PR #177 ("TODO 73, map side") shipped the
   map's Generate buttons and the light-year radius dialog (7.53.0);
   the item stays open in TODO.md for Generate this layer or slab. A
   citation of 73 after 04:57Z means the remaining part.
15. **Population after 7.40.0.** 51 to 54 left TODO.md with PR #169
   (04:37Z), replaced by an unnumbered "Population and Politics" note
   that lists the pages (POP.5) and a Galaxy Map territory overlay as
   still open. The overlay shipped in 7.54.0 (PR #179) and is POP.6,
   but TODO.md at 7.58.2 still lists it as open. The pass became
   optional and off by default in 7.58.2 (PR #180). The Plan section
   of TODO.md also still says "Extend the cache (8)" though 8 shipped
   in 7.56.0.

## Notes for the renumbering run

New IDs created here for finished items, recorded here as taken (they
never appear in TODO.md, since finished items are deleted):

- New parents: MAP.5 (Galaxy Map rework, with MAP.31 to MAP.42),
  SEC.2 (security audit findings, with SEC.3 to SEC.16 and
  SEC.17 to SEC.19), GEN.4 (nebulae, remnants and asteroid
  fields, GEN.10 to GEN.14), GEN.5 (known generation bugs, GEN.15 to
  GEN.20), GEN.7 (star population and bright stars, GEN.21, GEN.22).
- Subitems under an open parent: MAP.28 (nested ladder geometry) and
  MAP.29 (stage contents API), the first two drill-down pieces. In the
  tree IDs they were MAP.2.9 and MAP.2.10, numbered after the open
  drill-down items that already had MAP.2.1 to MAP.2.8.
- Single items: UX.6 to UX.12, MAP.6 to MAP.13, GEN.1 to GEN.3, GEN.6,
  NAV.1, NAV.2, DB.1, ADM.2, ADM.3, OPS.2 to OPS.4.
- Added in the 7.58.2 update: MAP.14 (bright stars on the Galaxy Map,
  done, never numbered) as the parent of bugs MAP.47 and MAP.48;
  POP.5 (population pages, open) and POP.6 (territories overlay, done),
  both unnumbered in TODO.md; and bug subitems under finished items
  (MAP.17, MAP.37, MAP.43) and under an open one (MAP.26).
- Added in the renumbering run (after `main` at 05:54Z, PR #188): 104
  to 108 became MAP.45 and MAP.46 (bugs under the finished MAP.11),
  UX.15 (a bug with no open item it breaks), GEN.8 and MAP.15; OPS.1
  took item 80's version-scheme questions; DOC.1 to DOC.3 (80 to 82)
  finished in the docs-refresh PR; MAP.20 (73), MAP.25 (78) and POP.5
  finished after 7.58.2. Then 109 to 111 (after `main` at 05:56Z,
  PR #187) became MAP.30 (under the finished MAP.3), MAP.44 and
  MAP.18 (bugs).
- Added after the renumbering: UX.16 (bug, button spacing); MAP.19
  (bug, big wedge picks), GEN.9 (multiple galaxies plan) and ADM.4
  (collapsible Generate page sections).
- Added in the bugs-and-security plan (after PR #192): SEC.20 to SEC.27
  (login logging, the per-username backoff in the control database, a
  trusted-device cookie, the `/account` password-guessing bug, a password
  blocklist, hashing cost, two-factor sign-in and a fail2ban recipe),
  then MAP.49 (bug, planet orbits drawn inside an asteroid belt) and
  MAP.50 (bug, names running off the edge of the map), then UX.19,
  UX.20 and ADM.9 (bugs: belt rows, scientific notation, the facility
  form), then SEC.28 (the always-on log; SEC.20 moved under it).
- Flat IDs (after PR #190): every dotted ID above was replaced by the
  next number in its category (see "Tree IDs to flat IDs"). A bug with no
  open item it breaks is now a top-level "(bug)" item (UX.15, UX.16)
  instead of a child of a `.0` item.
- The next free IDs are in "Next free IDs" at the top of this document.
- Finished bugs got their own IDs (for example GEN.15 to GEN.20, MAP.6
  to MAP.9) so they can be cited.
- Notes 3 and 10 are settled: MAP.4 and MAP.1 are done. TODO.md's
  stale Population and Plan text (note 15) was rewritten with the IDs.
