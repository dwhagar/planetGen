# TODO number map: old numbers and tree IDs to category IDs

## Next free IDs

A new item, subitem or bug takes its category's next free ID, and the
same PR moves that row up by one. The release scripts read this table:
each category's counter is its next free number minus one (the highest
ID ever issued in it, which is also its item count), and the version's
third number is the sum of the counters (see `changes/README.md`).
`scripts/bump_version.py --check` fails if `docs/TODO.md` uses an ID at
or past a category's next free one, so a stale row is caught before a
release is stamped.

| Category | Next free ID |
|---|---|
| UX | UX.92 |
| MAP | MAP.167 |
| NAV | NAV.58 |
| GEN | GEN.197 |
| PERF | PERF.58 |
| DB | DB.22 |
| API | API.24 |
| ADM | ADM.50 |
| SEC | SEC.33 |
| TEST | TEST.124 |
| USR | USR.10 |
| OPS | OPS.41 |
| DOC | DOC.17 |
| VIEW | VIEW.11 |
| POP | POP.11 |

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
| 87 | 2026-10-01 03:26Z to 05:29Z | UX.3 | Warn visitors while a background job changes the galaxy | open |
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
| ADM.15 | Change the worker count from the Queue page, with a "Ludicrous Speed" mode | none | open |
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
| ADM.43 | A full configuration page under Admin | none | open |
| ADM.44 | Web, Open Graph and SEO settings | none | open |
| ADM.45 | Prevalence fields take the override share directly and must total 100% | none | done, PR #991 |
| ADM.46 | Generate page: a progress line and per-layer counts instead of one line per sector | none | dropped |
| ADM.47 | Generating a neighbourhood from the Generate page shows no per-sector stats and looks slow or silent (bug) | none | done, PR #846 |
| ADM.48 | Two test_api_auth_sweep tests fail: /admin/stats/galaxy-settings/<name> answers 302, not 403, to an unauthorised caller (bug) | none | done, PR #863 |
| ADM.49 | Galaxy shape density settings: the user changes the density range of the spiral arms, the inter-arm space, the core and the bulge | none | open |
| API.1 | Create a system inside an existing sector | 10 (2026-09-24 01:32Z to 02:18Z); 7 (2026-09-24 01:57Z to 02:02Z); 5 (2026-09-24 02:25Z to 05:38Z); 4 (2026-09-24 02:53Z to 2026-09-30 17:58Z); 13 (2026-09-30 16:44Z); 14 (2026-09-30 16:49Z); 15 (2026-09-30 16:51Z to 18:09Z); 32 (2026-09-30 18:14Z); 37 (2026-09-30 18:39Z to 2026-10-01 04:24Z) | done in 7.43.0, PR #170 |
| API.2 | Edit a system's generated content | 11 (2026-09-24 01:32Z to 02:18Z); 8 (2026-09-24 01:57Z to 02:02Z); 6 (2026-09-24 02:25Z to 05:38Z); 5 (2026-09-24 02:53Z to 2026-09-30 17:58Z); 14 (2026-09-30 16:44Z); 15 (2026-09-30 16:49Z); 16 (2026-09-30 16:51Z to 18:09Z); 33 (2026-09-30 18:14Z); 38 (2026-09-30 18:39Z to 2026-10-01 04:24Z) | done in 7.43.0, PR #170 |
| API.3 | Remote generate: generate on a local machine, upload through the API | none | open |
| API.4 | API compatibility data in the docs | none | open |
| API.5 | API version and compatibility checking | none | open |
| API.6 | User-level API keys, owned by the account that created them, that can read but not upload | none | open |
| API.7 | Investigate and plan upload limits | none | open |
| API.8 | Verify uploaded data before it is finalized | none | open |
| API.9 | Key scopes | none | done, PR #1033 |
| API.10 | Reservations: claimed sectors and id blocks per run | none | open |
| API.11 | Staging tables | none | open |
| API.12 | The download: seed, skeleton and name state | none | open |
| API.13 | Generation without a database | none | open |
| API.14 | Upload routes, compressed, in batches | none | open |
| API.15 | Log every API call with its user, how it came in, and its HTTP response code | none | open |
| API.16 | The API reports the galaxy's seed, version and run history | none | dropped (Boss, 2026-10-09 20:42Z) |
| API.17 | Remote generation reproduces what the server would make | none | open |
| API.18 | Generate by recipe: JSON for sectors, systems, planets, moons and phenomena | none | open |
| API.19 | Galaxy-scale recipes: build a whole galaxy, piece by piece, from JSON | none | open |
| API.20 | `require_json_body` returns 500 for a deeply nested JSON body (bug) | none | done, PR #865 |
| API.21 | Flask-Limiter puts `Retry-After` on successful responses (bug) | none | done, PR #865 |
| API.22 | An API version number: one sequential integer, shown in admin and in the status response | none | open |
| API.23 | The object ID as the public reference: pages, URLs, the API, wiki links and objectref use it in place of row ids | none | open |
| DB.1 | Starbases, colonies and outposts in the database | 30 (2026-09-30 18:14Z); 35 (2026-09-30 18:39Z to 2026-10-01 02:24Z) | done in 7.35.0, PR #152 |
| DB.2 | Asteroid field and comet composition rows are written but never read (bug) | none | done, PR #347 |
| DB.3 | resetDb while another process holds id blocks can duplicate primary keys (bug) | none | done, PR #347 |
| DB.4 | A database with an emptied schema_migrations table is treated as current (bug) | none | done, PR #342 |
| DB.5 | Several first connections to an empty database race to create the schema (bug) | none | done, PR #342 |
| DB.6 | Store the galaxy's 128-bit seed, the version that made it, and every generation run | none | done, PR #387 |
| DB.7 | The version that generated each sector, and a warning for mixed-version galaxies | none | done, PR #813 |
| DB.8 | Check a galaxy database and say whether it is damaged | none | done, PR #893 |
| DB.9 | Repair a damaged galaxy database from a parity file | none | open |
| DB.10 | Repair reads the newest settings JSON and the pending deltas | none | dropped (Boss, 2026-10-09 20:42Z) |
| DB.11 | The database layer and migrations on SQLAlchemy and Alembic | none | done, PR #751 |
| DB.12 | MariaDB reads a stored -0.0 back as +0.0 (bug) | none | done, PR #442 |
| DB.13 | Every stored value in its own column, not in JSON blocks, and indexed for search | none | done, PR #766 |
| DB.14 | Sector stats keep the raw star statistics (mean age, summed luminosity, star count) and no baked color | none | done, PR #556 |
| DB.15 | A migration progress bar with the time remaining | none | done, PR #924 |
| DB.16 | Store the generator epoch and run id on each sector instead of four version text columns | none | open |
| DB.17 | Repair by regenerating a damaged sector from its seed when parity cannot rebuild it | none | open |
| DB.18 | Migration helpers for slow DDL: online indexes, instant columns and batched updates | none | open |
| DB.19 | Compact or derive phenomenon rows (needed only if the mass cut is lowered to 10 solar masses or less) | none | closed, not needed at the 20 solar mass cut (PR #866); reopen if the cut is lowered to 10 solar masses or less |
| DB.20 | Object IDs in the schema: uid becomes BINARY(10), unique on its own, plus an id_counters table | none | open |
| DB.21 | A deep pass for the database check: validate every star system, with the estimated time shown first | none | open |
| DOC.1 | Number TODO items by category (this renumbering) | 80 (2026-10-01 02:41Z to 05:29Z) | done in the docs-refresh PR (version-scheme questions moved to OPS.1) |
| DOC.2 | Architecture document | 81 (2026-10-01 02:51Z to 05:29Z) | done in the docs-refresh PR |
| DOC.3 | Design documents current, with reasons | 82 (2026-10-01 02:51Z to 05:29Z) | done in the docs-refresh PR |
| DOC.4 | Correct the stale statements the research found in docs, docstrings and comments | none | open |
| DOC.5 | Rewrite the object ID docs: object-ids.md, database-schema.md and api.md | none | open |
| DOC.6 | A static help section in the web interface: page template, index, per-page help links and a coverage test | none | open |
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
| GEN.1 | Real-world rates for interstellar objects | 5 (2026-09-30 18:14Z to 23:21Z) | done in 7.18.0, PR #136 |
| GEN.2 | Rogue planet mass bins | 6 (2026-09-30 18:14Z to 23:21Z) | done in 7.18.0, PR #136 |
| GEN.3 | A supermassive black hole in every galaxy | 7 (2026-09-30 18:14Z to 22:54Z) | done in 7.15.0, PR #132 |
| GEN.4 | Nebulae, remnants and asteroid fields (new parent) | none | done |
| GEN.5 | Known generation bugs of 2026-09-30 (new parent) | none | done |
| GEN.6 | Finish the correlative update (galactic motion) | 27 (2026-09-30 18:14Z); 32 (2026-09-30 18:39Z to 2026-10-01 02:35Z) | done in 7.37.0, PR #157 |
| GEN.7 | Star population and bright stars (new parent) | none | done |
| GEN.8 | Give rogue planets a planet class, with a rogue flag in the class constants | 106 (2026-10-01 05:43Z to 05:54Z) | done, PR #245 |
| GEN.9 | Plan for more than one galaxy in the database | none | open |
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
| GEN.24 | Generate the galactic core on layer 0 | none | open |
| GEN.25 | A moon reclassified after its planet moves can be too large for its planet (bug) | none | done, PR #350 |
| GEN.26 | Rogue planet surface conditions | none | done, PR #263 (schema v48; design docs/design/rogue-planet-surface.md) |
| GEN.27 | Class P (glaciated world) only in the habitable zone, and fitting there | none | open |
| GEN.28 | Seven new planet classes in the letter gaps (R, S, U, W, X, Y, Z) | none | open |
| GEN.29 | Sweep every planet class for sense once the new ones are in (bug) | none | open |
| GEN.30 | Bright-star thresholds: 1000 L_sun galaxy-wide, tiered backfill around generated sectors | none | done, PRs #295, #296 |
| GEN.31 | A point just under layer 0's top face lands in layer 1 (bug) | none | done, PR #353 |
| GEN.32 | Re-running an interrupted bright-star band draws it twice (bug) | none | done, PR #371 |
| GEN.33 | One class per PR, each with its tests | none | open |
| GEN.34 | Gas and ice giants come out too light, so there are no super-Jupiters (bug) | none | done, PR #350 |
| GEN.35 | Rocky planets only ever get Class D moons (bug) | none | done, PR #350 |
| GEN.36 | Moon regeneration can produce gas-giant or blacklisted moon classes (bug) | none | done, PR #350 |
| GEN.37 | 97% of planets land in the cold zone (bug) | none | done, PR #350 |
| GEN.38 | Rocky rogue planets over 10,000 km are still classed C (bug) | none | done, PR #415 |
| GEN.39 | The same seed can't reproduce the same galaxy (bug) | none | done, PR #381 |
| GEN.40 | Weed out sectors by star density before the bright-star backfill | none | open |
| GEN.41 | Investigate: how much backfill work a density pre-pass would save | none | done, PR #812 |
| GEN.42 | A pass that drops sectors from a region by probability | none | open |
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
| GEN.55 | Same seed, same data: a sector's contents depend only on the seed, the version and its address (internal) | none | open |
| GEN.56 | Every random draw in generation comes from the derived seeds | none | done, PR #791 |
| GEN.57 | A sector's contents depend only on the seed, the version and its address | none | open |
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
| GEN.90 | Refactor the planet classes around the habitability index | none | open |
| GEN.91 | Classes like S and V in the hot and cold zones | none | open |
| GEN.92 | Life and its highest stage follow the habitability score | none | open |
| GEN.93 | Nebula conditions in planet generation | none | open |
| GEN.94 | Feasibility study: can planets form in each nebula class, and what changes | none | done, PR #823 |
| GEN.95 | Nebula conditions applied when planets and surfaces are generated | none | open |
| GEN.96 | Generation directives for a sector (an override button) | none | done, PR #915 |
| GEN.97 | Generate N random neighborhoods | none | done, PR #950 |
| GEN.98 | Bright-star backfill from the farthest generated boundary outward | none | done, PR #795 |
| GEN.99 | Nebula volume backfill with the star types the nebula needs | none | open |
| GEN.100 | Phenomena placed galaxy-wide first and kept when sectors fill | none | done, PR #793 |
| GEN.101 | Fill order: nearest sectors first along a pruned Hilbert octree curve | none | open |
| GEN.102 | Investigate filling all near-zero-density void space at once | none | open |
| GEN.103 | Research where each star type and phenomenon belongs in the galaxy's structure | none | open |
| GEN.104 | A spin vector and a realistic axial tilt for every rotating object | none | done, PR #826 |
| GEN.105 | Orbital updates | none | open |
| GEN.106 | Movement thresholds and a next-update-due column | none | done, PR #802 |
| GEN.107 | The update reports how many objects moved, changed sector, or entered or left a nebula | none | done, PR #817 |
| GEN.108 | Orbital math limits: where each method breaks down and what happens there | none | done, PR #815 |
| GEN.109 | N-body influence from every object inside the largest nearby Hill sphere plus the galactic gradient, with a Hill-radius warning | none | open |
| GEN.110 | Rogue planet collisions: asteroid fields, merged giants and new stars | none | open |
| GEN.111 | Email the admin when two objects are inside each other's Hill radius | none | open |
| GEN.112 | Plan asteroid fields and belts as object systems for rendering | none | open |
| GEN.113 | Analyze the anomaly docs: which anomalies to add and how | none | open |
| GEN.114 | Add the chosen anomalies to the starmap | none | open |
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
| GEN.128 | Design: multi-star hierarchies and compact-object primaries | none | open |
| GEN.129 | Multi-star systems of up to seven stars | none | open |
| GEN.130 | Exotic star systems: a black hole, neutron star or similar at the center | none | open |
| GEN.131 | Bright-star scatter logs how many stars it added to each layer, by type | none | done, PR #805 |
| GEN.132 | Per-sector regional rates for neutron stars and black holes, extended system-wide where the research supports it | none | done, PR #807 |
| GEN.133 | Analysis of every star type's rate against its distance from the galactic core | none | done, PR #807 |
| GEN.134 | Tune the star populations to the observed star-formation profile by galactic radius | none | open |
| GEN.135 | Deterministic math helpers for every stored float, with a lint and frozen constants | none | open |
| GEN.136 | Fingerprint encoding: floats to 9 significant digits, a stored leaf digest per sector and a ring-and-layer tree | none | open |
| GEN.137 | Placed objects move along a galactic orbit but never by velocity times time, and phenomenon_scatter has no plan time (bug) | none | done, PR #843 |
| GEN.138 | Moon `hill_radius_km` uses the star's mass, so moon spacing and the orbit slider are wrong (bug) | none | done, PR #865 |
| GEN.139 | Orbit-update thresholds: per-object epoch, path-length rule and what the 0.01 mpc applies to (GEN.106 built) | none | open |
| GEN.140 | Orbital math guards the edge-case table adds (GEN.108 built) | none | open |
| GEN.141 | Faster Kepler solver (Mikkola or Markley) with brentq as fallback | none | open |
| GEN.142 | Peculiar velocity for rogue planets and asteroid fields | none | open |
| GEN.143 | A collision_events table, an admin report and a test that runs the whole collision path | none | open |
| GEN.144 | One distance for the Sun from the galactic centre across the constants, the density model and the design docs | none | open |
| GEN.145 | Class S atmosphere rule: S keeps air unless the shoreline ratio is over 30 | none | open |
| GEN.146 | Teff-dependent habitable zone from the Kopparapu table | none | open |
| GEN.147 | Classes N and Q carry a life chemical and an uncapped life timeline though they are lifeless (bug) | none | done, PR #865 |
| GEN.148 | Habitability index follow-ups from the research (GEN.84 built) | none | open |
| GEN.149 | Planetary-nebula central stars: 0.5 to 0.7 Msun, 1e2 to 1e4 Lsun, up to 2e5 K | none | open |
| GEN.150 | H II region radius and density from the ionizing photon rate, and IMF-based nebula hosts | none | open |
| GEN.151 | Supernova remnant sizes from the density-dependent Sedov-Taylor law (GEN.10 follow-up) | none | open |
| GEN.152 | Nebula cloud field is 10 to 40 times too full; lower it to the observed filling (GEN.47 rate check) | none | open |
| GEN.153 | Magnetar subtype of neutron star, and an age-dependent pulsar fraction | none | open |
| GEN.154 | Show the Einstein radius on compact-object pages | none | open |
| GEN.155 | A nuclear-cluster object for the Sgr A* sector (optional) | none | open |
| GEN.156 | Pin astropy to CODATA 2018 and IAU 2015 constants before its first import (GEN.66 follow-up) | none | open |
| GEN.157 | The `neighbor_galaxies` table and a verified data file of about 25 real galaxies | none | open |
| GEN.158 | Add a metallicity value to stars | none | open |
| GEN.159 | Globular clusters: cluster table, King tables and the Milky Way catalogue | none | open |
| GEN.160 | Cluster density in the sector gate, with a "cluster" population | none | open |
| GEN.161 | Bright-first fill for cluster sectors | none | open |
| GEN.162 | Planet cull and blue stragglers in clusters | none | open |
| GEN.163 | Type-B pulsar planets in globular clusters (GEN.130 follow-on) | none | open |
| GEN.164 | Synthetic globular-cluster systems for generated galaxies | none | open |
| GEN.165 | A comet's time-0 scene position disagrees with its stored position in a rare random system (bug) | none | done, PR #948 |
| GEN.166 | A mass_range argument on NeutronStar and BlackHole, and intermediate-mass black holes as their own kind | none | done, PR #866 |
| GEN.167 | A lowest-mass option for the phenomenon scatter: --phenomenon-min-mass, default 20 solar masses | none | done, PR #866 |
| GEN.168 | The sector fill draws the phenomena below the scatter cut | none | done, PR #866 |
| GEN.169 | Decide the phenomenon scatter rates: regional factors and the 0.1% intermediate-mass black holes | none | open |
| GEN.170 | Object ID layout: an 80-bit ID of birth sector, serial and body number, with pack, unpack, format and parse functions | none | done, PR #1017 |
| GEN.171 | The sector fill gives object IDs by generation rank | none | open |
| GEN.172 | Run-time births get object IDs from the counters | none | open |
| GEN.173 | Deleting a body and then adding one fails with IntegrityError 1062 on uq_planets_uid (bug) | none | done, PR #987 |
| GEN.174 | Bodies an admin adds are saved with a NULL uid (bug) | none | done, PR #987 |
| GEN.175 | Regenerating a phenomenon sets its uid to NULL (bug) | none | done, PR #960 |
| GEN.176 | A nebula or remnant is born in the sector holding the centre of the space it occupies | none | open |
| GEN.177 | Planetary magnetic fields: a stagnant-lid factor | none | open |
| GEN.178 | Magnetic fields: the induced field of an ocean moon | none | open |
| GEN.179 | Store each sector's generation directive and attempt record with the sector | none | open |
| GEN.180 | Directives: a forced fill (met_forced) after K failed draws | none | open |
| GEN.181 | Directives: refuse impossible requests up front from the compound-Poisson tables | none | open |
| GEN.182 | Some comets' orbits do not bring them back: they are ejected into space (bug) | none | done, PR #938 |
| GEN.183 | A mass cut the user sets: preset values on a slider from 8 to 20 solar masses | none | done, PR #969 |
| GEN.184 | A luminosity floor the user sets: presets from 2500 to 4 million solar luminosities, default 3000 | none | done, PR #965 |
| GEN.185 | The scatter in five passes: mass-limit objects, mass-limit stars, brightest-star sector marks, the luminosity pass, then other phenomena | none | done, PR #953 |
| GEN.186 | Random neighborhoods: an option to keep away from filled space | none | open |
| GEN.187 | Bright-star back scatter by mass: rings of 1, 2, 5 and 8 solar masses around filled space | none | open |
| GEN.188 | Mass limit default is 8, and the mass slider and luminosity dropdown sit side by side on Generate and New galaxy | none | done, PR #1007 |
| GEN.189 | Gamma-ray burst and AGN ozone loss for the ozone_loss_flag (low priority) | none | open |
| GEN.190 | The Phenomena table is empty after a web-generated galaxy: the phenomena scatter pass never runs (bug) | none | done, PR #1008 |
| GEN.191 | New galaxy ignores the mass limit slider: the plan step does not store the limit, so the scatter uses the default whatever the form says (bug) | none | done, PR #1008 |
| GEN.192 | Phenomena scatter log is thin and the Phenomena table page needs checking after a run (bug) | none | done, PR #1013 |
| GEN.193 | Phenomena table stays empty after the scatter: scattered unbuilt phenomena are not listed (bug) | none | done, PR #1019 |
| GEN.194 | Default mass limit 14 solar masses and default luminosity floor 9,000 solar luminosities | none | done, PR #1051 |
| GEN.195 | A separate mass limit for neutron stars and black holes, and the central black hole or quasar always created | none | open |
| GEN.196 | One "Redo scatters" box on Generate: choose which scatters to redo, with new settings for each | none | open |
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
| MAP.59 | Make it plain that a zoomed-in slab is a slab, not a wedge | none | open |
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
| MAP.121 | Every map shows and steps to its neighbouring regions, on one map engine | none | open |
| MAP.122 | A Select mode on every galaxy view: Galaxy (blocks and sectors) or Star | none | open |
| MAP.123 | Show or hide star types and phenomena, and set the luminosity floor, on the Galaxy and Sector Maps | none | done, PR #753 |
| MAP.124 | The Galaxy Map opens zoomed to fit all charted space | none | done, PR #696 |
| MAP.125 | Infinite zoom: one 3D interface from the galaxy down to a moon | none | done, PR #683 |
| MAP.126 | Show the orbital trajectories of selected objects in their frame of reference | none | done, PR #686 |
| MAP.135 | Selecting the first slab or wedge shows its bounds (bug) | none | done, PR #788 |
| MAP.136 | Binary stars pick as one system and their 3D orbits are drawn clearly (bug) | none | done, PR #797 |
| MAP.137 | Rogue planets are easy to see, with the right default filters (bug) | none | done, PR #797 |
| MAP.138 | Recenter the camera in every 3D view (bug) | none | done, PR #797 |
| MAP.139 | The Galaxy View uses its spare space: an info box with a Details link, and menu items | none | open |
| MAP.140 | Double-click on a selected object goes there and opens its information | none | open |
| MAP.141 | Context around the selection: faint neighbours, and the sectors above and below | none | open |
| MAP.142 | Nebulae have fuzzy, fading boundaries | none | open |
| MAP.143 | Color sectors by their number of habitable locations | none | open |
| MAP.144 | Replace `THREE.Clock` with `THREE.Timer` in `phenomenonrender.js` (bug) | none | done, PR #865 |
| MAP.145 | A sky and Galaxy Map drawing rule for neighbour galaxies | none | open |
| MAP.146 | Fly through the galaxy: scroll-zoom, double-click flight, distance-based visibility and a see-through near field | none | open |
| MAP.147 | The Galaxy Map wire format: the investigation (done) and the record of what was built from it | none | open |
| MAP.148 | The star visibility law: apparent-magnitude opacity, flux-based brightness and an on-screen limit from a histogram | none | open |
| MAP.149 | The near field: depth fade, a see-through focus tube, drawing from inside a container, and picking that matches what is drawn | none | open |
| MAP.150 | The free camera: wheel zoom to the cursor, double-click flight, and the observer inside, with the container named from position | none | open |
| MAP.151 | The region data layer: exact-centred frame, aligned cells, per-level aggregates and slot-wrap ranges | none | open |
| MAP.152 | Scale hand-offs: galaxy, sector and system cross-fade with hysteresis, and per-tile camera-relative origins | none | open |
| MAP.153 | Stars fade in with the zoom: a birth radius from each star's rank in its tile list (first client stage) | none | done, PR #998 |
| MAP.154 | Nested bright-star lists on the server, so every parent list is a subset of its child's | none | open |
| MAP.155 | Other objects fade in too: point objects from level 8, a size ramp for cloud sprites, and stars that grow from a faint dot | none | open |
| MAP.156 | View one layer or a range of layers top-down from the galaxy view, as a secondary option | none | open |
| MAP.157 | Trim the Galaxy Map tile JSON and serve it from prebuilt, precompressed bytes | none | open |
| MAP.158 | A gentler tile prefetch and an IndexedDB tile cache instead of localStorage | none | open |
| MAP.159 | Packed binary Galaxy Map tiles (quantised planes) with the nested tile lists, on one cache stamp bump | none | open |
| MAP.160 | Quantise the Galaxy Map GPU buffers (deferred) | none | open |
| MAP.161 | Load the Galaxy Map faster on a first visit: bundle or preload its scripts | none | open |
| MAP.162 | A sector holding scattered objects but never generated can still be opened, marked uncharted | none | open |
| MAP.163 | Galaxy Map brightness scale: floor 2,500 L_sun at full zoom, then min/max scaling per zoom | none | done, PR #1029 |
| MAP.164 | Galaxy Map shows phenomena: black holes purple, neutron stars dark blue, sized by mass, dark colors out-shine brighter stars | none | done, PR #1029 |
| MAP.165 | Scattered phenomena store a mass so the Galaxy Map sizes them exactly | none | open |
| MAP.166 | Galaxy Map "Dimmest star shown" says "every star" only when the view is complete | none | open |
| NAV.1 | Courses in "bearing mark mark" on nested frames | 28 (2026-09-30 18:14Z); 33 (2026-09-30 18:39Z to 2026-10-01 02:57Z) | done in 7.14.0, PR #130 (see note 4) |
| NAV.2 | Warp and fold speeds | 29 (2026-09-30 18:14Z); 34 (2026-09-30 18:39Z to 21:54Z) | done in 7.8.0, PR #121 |
| NAV.3 | One shared picker for the Galaxy, Sector and System displays | none | done, PR #724 |
| NAV.4 | Save a course | none | open |
| NAV.5 | Show a course on the Galaxy Map | none | open |
| NAV.6 | Courses that steer clear of gravity wells | none | open |
| NAV.7 | One reference for every object, with its parents | none | done, PR #688 |
| NAV.8 | Pages and anchors for stars, planets, moons and belts | none | open |
| NAV.9 | Search and locate return references for every kind | none | open |
| NAV.10 | Routing that scales past a few thousand systems | none | done, PR #814 |
| NAV.11 | Travel times for the system-to-system route too | none | open |
| NAV.12 | No maximum hop length: a route always reaches the nearest star it can, across any number of sectors | none | done, PR #841 |
| NAV.13 | A picker module: select, step out, step in, step sideways | none | done, PR #564 |
| NAV.14 | One breadcrumb for every level | none | done, PR #571 |
| NAV.15 | Pick mode everywhere | none | done, PR #566 |
| NAV.16 | NAV endpoints can be any object | none | done, PR #697 |
| NAV.17 | A saved course record with both forms | none | open |
| NAV.18 | Save, list, open, rename and delete, per browser | none | open |
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
| NAV.42 | Each route stop shows the course and distance to the next stop | none | open |
| NAV.43 | Find everything within a distance of a place: the query and the API | none | done, PR #804 |
| NAV.44 | A "What's nearby" page: pick a place, enter a distance in parsecs, list what is there | none | done, PR #804 |
| NAV.45 | "What's within N pc" from the Galaxy Map and Sector Map | none | open |
| NAV.46 | The NAV picker can't click galaxy wedges to zoom in (bug) | none | closed, already fixed; covered by a browser test (#568) |
| NAV.47 | Unknown-space jumps stop at scattered stars, black holes, neutron stars and quasars | none | open |
| NAV.48 | Offer to generate the uncharted sectors that block a course | none | open |
| NAV.49 | Waypoints: pick objects in Star select mode and plot a course through them, kept on the map until cleared | none | open |
| NAV.50 | Pick any object down to a moon as a NAV endpoint | none | done, PR #724 |
| NAV.51 | Courses route around asteroid fields | none | open |
| NAV.52 | Port `join_islands` and the k-d tree to cKDTree | none | open |
| NAV.53 | A `cells_touching_sphere` helper, and pad `store.sectors_reached_by` by one edge (bug) | none | done, PR #843 |
| NAV.54 | Keep-out radii for asteroid fields, supermassive holes, moons and nebulae (NAV.24 built) | none | open |
| NAV.55 | A tuning block for the keep-out knobs | none | open |
| NAV.56 | Census of overlapping keep-out spheres in a generated galaxy | none | open |
| NAV.57 | The Intergalactic Frame in navigation-frames.md and `navigation.py` | none | open |
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
| OPS.11 | Define "the same galaxy" and which versions stay reproducible | none | done, PR #362 (docs/design/reproducible-galaxies.md) |
| OPS.12 | `generate.py reproduce`: a version and a seed rebuild a galaxy and check it | none | dropped (Boss, 2026-10-09 20:42Z) |
| OPS.13 | Every update records the version key, keeping the last 10 | none | done, PR #810 |
| OPS.14 | A warning when the running version key differs from the galaxy's | none | done, PR #876 |
| OPS.15 | Each update says whether it changes generated output | none | open |
| OPS.16 | A daily maintenance script for Linux, macOS and Windows | none | open |
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
| OPS.28 | Generator epoch and battery digest: say whether two checkouts generate the same galaxy | none | open |
| OPS.29 | Update reload: cover the gunicorn units and non-Apache hosts, and say what a reload aborts | none | open |
| OPS.30 | A lock helper for the maintenance run | none | open |
| OPS.31 | Lint every example plist, XML and service file in CI | none | open |
| OPS.32 | `examples/macos/org.planetgen.update.plist` is not well-formed XML, so the update daemon silently fails to install (bug) | none | done, PR #865 |
| OPS.33 | `.gitattributes` has no LF pins for the lock files and word list, so hashes differ between a Windows and a Linux checkout (bug) | none | done, PR #843 |
| OPS.34 | Windows Redis in WSL: fix the keep-alive advice and add a Start-RedisInWsl remedy | none | done, PR #981 |
| OPS.35 | A vendored-version lock file for the static libraries | none | open |
| OPS.36 | Space and size checks measure the boot drive, not the drive holding the database (bug) | none | done, PR #868 |
| OPS.37 | A Generator version number: one sequential integer, shown in admin and the API | none | open |
| OPS.38 | Two Redis dump files (dump.rdb and src/dump.rdb) are committed to main and should be removed and git-ignored (bug) | none | done, PR #967 |
| OPS.39 | Remove Windows support; keep only a simple docs/WINDOWS.md | none | done, PR #981 |
| OPS.40 | update.sh step 8 fails: setup-debug-log.sh loads the deleted util/appconfig.py (bug) | none | done, PR #989 |
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
| PERF.18 | Run the GEN.30 bright-star backfill in parallel on the work queue | none | open |
| PERF.19 | Everything the API or web site starts runs on the work queue (investigate) | none | done, PR #492 |
| PERF.20 | Short-term caching through the work queue and API (needs planning) | none | open |
| PERF.21 | Generation works with any worker count: the parallel path is built, used and tested (bug) | none | done, PR #371 |
| PERF.22 | On Python 3.12 a run hangs forever when a worker process dies (bug) | none | done, PR #371 |
| PERF.23 | The bright-star progress bar can end at 101% (bug) | none | done, PR #371 |
| PERF.24 | The work queue and web jobs on Redis with RQ | none | done, PR #517 |
| PERF.25 | The page cache on cachetools; the tile cache stays | none | done, PR #478 |
| PERF.26 | Size estimates don't match what generation stores (bug) | none | done, PR #470 |
| PERF.27 | On Python 3.9, `--workers=--` comes back as a list and crashes generate.py's option checks (bug) | none | done, PR #442 |
| PERF.28 | The console's second progress bar (the bright-star backfill after a sector run) never updates its ETA (bug) | none | done, PR #454 |
| PERF.29 | Record which runs a partly filled sector still needs | none | open |
| PERF.30 | Finish an interrupted block or sector run on the next start | none | open |
| PERF.34 | The site stays responsive during heavy generation jobs (bug) | none | done, PR #811 |
| PERF.31 | Investigate: where generation spends its time, from the plan to a finished galaxy | none | open |
| PERF.32 | Generation performance stats: rates recorded per run, deleted on every new version | none | done, PR #905 |
| PERF.33 | Progress bars and ETAs from measured performance | none | open |
| PERF.35 | An interval or chunk ledger for untouched sectors once block-first backfill lands | none | open |
| PERF.36 | Memory and request guard: never list more than about 50,000 candidate cells, and refuse huge enumerations in a web request | none | open |
| PERF.37 | `DecayingRate` starts from the first single completion, so the ETA is up to twice too long early in a run (bug) | none | done, PR #865 |
| PERF.38 | Cache fixes for the Galaxy Map under a fill: single-flight tile builds, a busy rule for the page cache, a deletion epoch in place of COUNT(*) | none | open |
| PERF.39 | Every API job costs 2.4 s and 195 MB: import lazily and cap the burst workers | none | open |
| PERF.40 | Two shared queues, a reserved interactive worker and a real "cancel now" | none | open |
| PERF.41 | Stamp the tile and page caches with a TILE_FORMAT constant instead of the version (optional) | none | open |
| PERF.42 | Warm the RQ worker before the fork: pre-import generation modules and build the bright-star table once | none | done, PR #841 |
| PERF.43 | Lazy word-salad names for phenomena named by object ID | none | done, PR #861 |
| PERF.44 | Compute object uids in Python and write them with the row | none | done, PR #863 |
| PERF.45 | Nearest-system links and containment as one later pass | none | done, PR #863 |
| PERF.46 | Planets and moons: set the position once per body | none | open |
| PERF.47 | The PERF.31 benchmark records the buffer pool, table sizes and worker start-up cost | none | open |
| PERF.48 | Low priority: a numeric-only INSERT formatter or C driver for bright_stars and phenomenon_scatter | none | open |
| PERF.49 | Batch system-name reservation: remove the quadratic scan and the long-held registry locks (re-measure first) | none | done, PR #870 |
| PERF.50 | A progress bar inside one sector's save: workers report their sub-steps to the main process | none | done, PR #922 |
| PERF.51 | One progress mechanism for every sub-step: a bar starts by itself when a step is predicted to take over 15 seconds | none | done, PR #919 |
| PERF.52 | Admin generation-stats table is wrong (bug) | none | done, PR #994 |
| PERF.53 | Timing stats are skewed by layers and sectors that generated nothing (bug) | none | done, PR #1015 |
| PERF.54 | Separate generation-stats rows for the mass pass and the luminosity pass of the star scatter | none | open |
| PERF.55 | One global progress bar for generation jobs that run in phases, with an ETA across all phases | none | open |
| PERF.56 | Record how long each stage of a staged job takes, with the settings it ran with | none | open |
| PERF.57 | Stop a layer-walking scatter early once the last 100 layers produced no stars | none | open |
| POP.1 | Government ownership of systems | 12 (2026-09-24 01:32Z to 02:18Z); 9 (2026-09-24 01:57Z to 02:02Z); 7 (2026-09-24 02:25Z to 05:38Z); 6 (2026-09-24 02:53Z to 2026-09-30 17:58Z); 15 (2026-09-30 16:44Z); 16 (2026-09-30 16:49Z); 17 (2026-09-30 16:51Z to 18:09Z); 34 (2026-09-30 18:14Z); 39 (2026-09-30 18:39Z to 18:41Z); 59 (2026-09-30 19:02Z to 19:17Z); 63 (2026-09-30 20:01Z to 20:27Z); 46 (2026-09-30 20:07Z); 50 (2026-09-30 20:08Z to 20:48Z); 51 (2026-09-30 20:43Z to 2026-10-01 04:37Z) | done in 7.49.0, PR #169 (optional since 7.58.2, PR #180) |
| POP.2 | Names for dominant species on living worlds | 13 (2026-09-24 01:32Z to 02:18Z); 10 (2026-09-24 01:57Z to 02:02Z); 8 (2026-09-24 02:25Z to 05:38Z); 7 (2026-09-24 02:53Z to 2026-09-30 17:58Z); 16 (2026-09-30 16:44Z); 17 (2026-09-30 16:49Z); 18 (2026-09-30 16:51Z to 18:09Z); 35 (2026-09-30 18:14Z); 40 (2026-09-30 18:39Z to 18:41Z); 60 (2026-09-30 19:02Z to 19:17Z); 64 (2026-09-30 20:01Z to 20:27Z); 47 (2026-09-30 20:07Z); 51 (2026-09-30 20:08Z to 20:48Z); 52 (2026-09-30 20:43Z to 2026-10-01 04:37Z) | done in 7.49.0, PR #169 |
| POP.3 | Database of spacefaring species | 14 (2026-09-24 01:32Z to 02:18Z); 11 (2026-09-24 01:57Z to 02:02Z); 9 (2026-09-24 02:25Z to 05:38Z); 8 (2026-09-24 02:53Z to 2026-09-30 17:58Z); 17 (2026-09-30 16:44Z); 18 (2026-09-30 16:49Z); 19 (2026-09-30 16:51Z to 18:09Z); 36 (2026-09-30 18:14Z); 41 (2026-09-30 18:39Z to 18:41Z); 61 (2026-09-30 19:02Z to 19:17Z); 65 (2026-09-30 20:01Z to 20:27Z); 48 (2026-09-30 20:07Z); 52 (2026-09-30 20:08Z to 20:48Z); 53 (2026-09-30 20:43Z to 2026-10-01 04:37Z) | done in 7.49.0, PR #169 |
| POP.4 | Younger and older civilizations | 15 (2026-09-24 01:32Z to 02:18Z); 12 (2026-09-24 01:57Z to 02:02Z); 10 (2026-09-24 02:25Z to 05:38Z); 9 (2026-09-24 02:53Z to 2026-09-30 17:58Z); 18 (2026-09-30 16:44Z); 19 (2026-09-30 16:49Z); 20 (2026-09-30 16:51Z to 18:09Z); 37 (2026-09-30 18:14Z); 42 (2026-09-30 18:39Z to 18:41Z); 62 (2026-09-30 19:02Z to 19:17Z); 66 (2026-09-30 20:01Z to 20:27Z); 49 (2026-09-30 20:07Z); 53 (2026-09-30 20:08Z to 20:48Z); 54 (2026-09-30 20:43Z to 2026-10-01 04:37Z) | done in 7.49.0, PR #169 |
| POP.5 | Population pages: species, polities, dominant species, territory (unnumbered in TODO.md) | none | done after 7.58.2, PR #184 |
| POP.6 | Territories overlay on the Galaxy Map (unnumbered in TODO.md) | none | done in 7.54.0 (PR #179) and 7.58.1 (PR #181) |
| POP.7 | Tech levels for technological species | none | open |
| POP.8 | A tech-level design from the six domain indices | none | open |
| POP.9 | A tech level generated for every technological species | none | open |
| POP.10 | Facility types: programmable, picked from a dropdown, with affiliation and Green/Yellow/Red ratings | none | open |
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
| SEC.32 | Argon2id password hashing, a 1,024-character password limit and the `__Host-` cookie prefix | none | open |
| USR.1 | User accounts | none | open |
| USR.2 | Accounts with roles: user, admin and Owner | 64 (2026-10-01 02:13Z to 05:29Z) | open |
| USR.3 | SMTP settings in the admin config | 65 (2026-10-01 02:13Z to 05:29Z) | open |
| USR.4 | Invite-only sign-up by unique link | 66 (2026-10-01 02:13Z to 05:29Z) | open |
| USR.5 | Email loop for setting and resetting passwords | 67 (2026-10-01 02:13Z to 05:29Z) | open |
| USR.6 | Owner transfer | 68 (2026-10-01 02:13Z to 05:29Z) | open |
| USR.7 | A user-level interface with bookmarks | 69 (2026-10-01 02:13Z to 05:29Z) | open |
| USR.8 | Every signed-in user can generate a one-off system | none | open |
| USR.9 | `seo.privacy_note` text, "download my data" and "delete my account" on the account page | none | open |
| UX.0 | Bugs and small fixes (standing item) | none | open while it holds bugs |
| UX.1 | Class reference pages | 56 (2026-10-01 01:19Z to 04:37Z) | done in 7.46.0, PR #167 |
| UX.2 | Menus sized to what they hold (bug) | 62 (2026-10-01 01:44Z to 05:29Z) | done, PR #528 |
| UX.3 | Warn every visitor while a background job changes the galaxy | 87 (2026-10-01 03:26Z to 05:29Z) | open |
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
| UX.22 | Meaningful units for every measurement | none | open |
| UX.23 | A shared unit-ladder module | none | open |
| UX.24 | Sector contents: rogue planets after systems and phenomena, and expanded rows the full table width (bug) | none | done, PR #487 |
| UX.25 | Rogue planets: octant and a small map symbol beside each name (bug) | none | done, PR #490 |
| UX.26 | Edit and admin actions as a button that opens a menu (bug) | none | done, PR #558 |
| UX.27 | System page: the system and navigation buttons on one row that doesn't overlap (bug) | none | done, PR #539 |
| UX.28 | Investigate icons instead of words on buttons | none | done, PR #490 |
| UX.29 | Every comet in a system shows its type as a link (bug) | none | done, PR #487 |
| UX.30 | Planet information without the Markdown render | none | open |
| UX.31 | Editing a star system: an edit button with a quick menu, not a long panel (bug) | none | done, PR #533 |
| UX.32 | Planet rows show the class only, without the type and moon labels | none | open |
| UX.33 | Filter phenomena by their classes and types (bug) | none | done, PR #606 |
| UX.34 | The sector summary calls white dwarfs "B-type" and "A-type" systems (bug) | none | done, PR #448 |
| UX.35 | The NAV page route shown horizontally, wrapping onto several lines on narrow screens | none | done, PR #879 |
| UX.36 | Scientific notation starts too early for whole numbers (bug) | none | done, PR #457 |
| UX.37 | A UX sweep: remove redundant and duplicate controls so the interface gets out of the way | none | done, PR #671 |
| UX.38 | The nebula and remnant diagrams' "-" button does nothing at the 1 ly limit (bug) | none | done, PR #590 |
| UX.39 | Markdown rendered by the markdown library | none | done, PR #759 |
| UX.40 | Buttons, menus and dialogs from Shoelace web components | none | done, PR #544 |
| UX.41 | Tables on TanStack Table and TanStack Virtual | none | done, PR #606 |
| UX.42 | In-universe wording across the interface | none | open |
| UX.43 | A visual design built like a pilot's starmap and navigation console | none | open |
| UX.44 | Search: mutually exclusive tags should combine with OR, the rest with AND (bug) | none | done, PR #457 |
| UX.45 | Bookmark management | none | open |
| UX.46 | Wiping the galaxy wipes the bookmarks | none | open |
| UX.47 | A bookmark manager: list, rename, sort, group and delete | none | open |
| UX.48 | Charted regions as bookmarks that frame and outline the region | none | open |
| UX.49 | Form fields as Shoelace components (sl-input, sl-select, sl-checkbox) across the site | none | open |
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
| UX.78 | Unit preference: Automatic, Metric only or Customary | none | open |
| UX.79 | Python and JS round half-way values differently (bug) | none | done, PR #843 |
| UX.80 | Negative values that round to zero print "-0" (bug) | none | done, PR #843 |
| UX.81 | Time symbols Gyr, Myr, kyr in place of Gy, My, ky; AU from 1,000,000 km; scientific text below mantissa 1e-3 | none | open |
| UX.82 | Theme checks after PR #800: SVG currentColor, two Shoelace contrast failures, alpha in --bg-subtle | none | open |
| UX.83 | Generation steps that run long show no progress bar of their own: linking new sectors to their neighbours, the phenomenon scatter and others (bug) | none | done, PR #895 |
| UX.84 | Every sub-step must show a progress bar that starts by itself when it is predicted to take over 15 seconds (bug) | none | done, PR #930 |
| UX.85 | Button menus open out of sight and make the user scroll to see them (bug) | none | done, PR #932 |
| UX.86 | The Galaxy Map controls take too much room: buttons too large and filters one character wide (bug) | none | done, PR #936 |
| UX.87 | The system list shows uncharted systems: every scattered star, with its location and a way to generate it | none | open |
| UX.88 | Say "uncharted" instead of "unbuilt" or "not generated" in all user-facing text outside Generate and admin | none | done, PR #1049 |
| UX.89 | Staged jobs show wrong stage counts and numbers, and skipped stages are not listed (bug) | none | open |
| UX.90 | Explain the habitability chips in the web interface: a visible legend with the equipment labels, colours, thresholds and scores | none | open |
| UX.91 | Planet and moon description carries a full PHI-4 explanation, each colour factor and why | none | open |
| VIEW.1 | View from a planet | none | open |
| VIEW.2 | A starmap seen from a planet. RESEARCH WITH BOSS FIRST | 83 (2026-10-01 02:55Z to 05:29Z) | open |
| VIEW.3 | Render the view as a PNG, with constellations | 84 (2026-10-01 02:55Z to 05:29Z) | open |
| VIEW.4 | Constellation names in the name generator | 85 (2026-10-01 02:55Z to 05:29Z) | open |
| VIEW.5 | Light-travel positions: where an object appears to a distant observer | none | done, PR #719 |
| VIEW.6 | Settle the handedness of the generated galaxy before mapping real sky coordinates | none | open |
| VIEW.7 | Sky tiers and the backfill cost measurement | none | open |
| VIEW.8 | A dust model for the sky: the galactic extinction law and `nebulae.extinction_av` | none | open |
| VIEW.9 | Replace the BC_V table with a published relation (Flower 1996 or Torres 2010) | none | open |
| VIEW.10 | Declare `pillow` in `setup.py` if it is used, and add a HYG calibration test | none | open |

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
| TEST.77 | A golden-seed regression test | none | open |
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
| TEST.110 | Object ID tests: identical IDs on 1 and 4 workers, none reused, none missing | none | open |
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
| TEST.122 | Browser map tests fail on plain main in Bugfixes lane 2's container (fixture maps, controls, system-page maps) (bug) | none | open |
| TEST.123 | test_a_loaded_sector_knows_every_objects_cell_and_velocity fails once under full-suite load (bug) | none | open |
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
| 7.49.0 | "Design: docs/design/population-and-politics.md (TODO 51-54)" | POP.1 to POP.4 |

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
| [navigation-frames.md](navigation-frames.md) | "TODO items 33 and 34" (first written as "items 28 and 29" at 18:14Z) | 2026-09-30 18:39Z | NAV.1, NAV.2 |
| [interstellar-object-rates.md](interstellar-object-rates.md) | "TODO items 5 (rates) and 6", "TODO item 27" | 2026-09-30 20:27Z | GEN.1, GEN.2, GEN.10 |
| [nebula-and-asteroid-field-classes.md](nebula-and-asteroid-field-classes.md) | "items 27-31", "item 28" ... "item 31" | 2026-09-30 18:39Z | GEN.10 to GEN.14 |
| [galaxy-drilldown-navigation.md](galaxy-drilldown-navigation.md) | "items 64-73", then "items 70-79" and the build-order table; "TODO 19"; "TODO 63" | 2026-10-01 02:24Z, 02:27Z | MAP.16 to MAP.29; MAP.1; MAP.3 |
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
| [population-and-politics.md](population-and-politics.md), `population.py`, `api/population.py`, `_db.py`, `schema.sql`, `program_constants.py`, `generate.py`, `test_population.py`, `docs/database-schema.md` | "TODO 51-54", "TODO items 51-54" | 2026-10-01 04:09Z to 05:28Z | POP.1 to POP.4 |
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
