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
| UX | UX.23 |
| MAP | MAP.57 |
| NAV | NAV.3 |
| GEN | GEN.30 |
| PERF | PERF.18 |
| DB | DB.2 |
| API | API.9 |
| ADM | ADM.14 |
| SEC | SEC.29 |
| TEST | TEST.70 |
| USR | USR.8 |
| OPS | OPS.6 |
| DOC | DOC.4 |
| VIEW | VIEW.5 |
| POP | POP.7 |

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
| 62 | 2026-10-01 01:44Z to 05:29Z | UX.2 | Menus sized to what they hold | open |
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
| ADM.4 | Collapsible Generate page sections; pick the center sector | none | open |
| ADM.5 | Central validate module | 57 (2026-10-01 01:15Z to 05:29Z) | done, PR #235 |
| ADM.6 | Override a planet's or moon's class | 58 (2026-10-01 01:15Z to 05:29Z) | done, PR #260 |
| ADM.7 | Override a star | 59 (2026-10-01 01:15Z to 05:29Z) | done, PR #260 |
| ADM.8 | Delete and regenerate buttons, sector down | 60 (2026-10-01 01:15Z to 05:29Z) | done, PR #244 |
| ADM.9 | "Place a facility": host by placement, log orbit slider, moving belt facilities (bug) | none | done, PR #211 |
| ADM.10 | Admin page to view and manage the work queue | none | open |
| ADM.11 | Jobs keep running after the browser closes | none | open |
| ADM.12 | Jobs as a tree, with timing for every node | none | open |
| ADM.13 | Incomplete uploads page | none | open |
| API.1 | Create a system inside an existing sector | 10 (2026-09-24 01:32Z to 02:18Z); 7 (2026-09-24 01:57Z to 02:02Z); 5 (2026-09-24 02:25Z to 05:38Z); 4 (2026-09-24 02:53Z to 2026-09-30 17:58Z); 13 (2026-09-30 16:44Z); 14 (2026-09-30 16:49Z); 15 (2026-09-30 16:51Z to 18:09Z); 32 (2026-09-30 18:14Z); 37 (2026-09-30 18:39Z to 2026-10-01 04:24Z) | done in 7.43.0, PR #170 |
| API.2 | Edit a system's generated content | 11 (2026-09-24 01:32Z to 02:18Z); 8 (2026-09-24 01:57Z to 02:02Z); 6 (2026-09-24 02:25Z to 05:38Z); 5 (2026-09-24 02:53Z to 2026-09-30 17:58Z); 14 (2026-09-30 16:44Z); 15 (2026-09-30 16:49Z); 16 (2026-09-30 16:51Z to 18:09Z); 33 (2026-09-30 18:14Z); 38 (2026-09-30 18:39Z to 2026-10-01 04:24Z) | done in 7.43.0, PR #170 |
| API.3 | Remote generate: generate locally, upload through the API | none | open |
| API.4 | API compatibility data in the docs | none | open |
| API.5 | API version and compatibility checking | none | open |
| API.6 | Admin-created user-level API keys that can read but not upload | none | open |
| API.7 | Investigate and plan upload limits | none | open |
| API.8 | Verify uploaded data before it is finalized | none | open |
| DB.1 | Starbases, colonies and outposts in the database | 30 (2026-09-30 18:14Z); 35 (2026-09-30 18:39Z to 2026-10-01 02:24Z) | done in 7.35.0, PR #152 |
| DOC.1 | Number TODO items by category (this renumbering) | 80 (2026-10-01 02:41Z to 05:29Z) | done in the docs-refresh PR (version-scheme questions moved to OPS.1) |
| DOC.2 | Architecture document | 81 (2026-10-01 02:51Z to 05:29Z) | done in the docs-refresh PR |
| DOC.3 | Design documents current, with reasons | 82 (2026-10-01 02:51Z to 05:29Z) | done in the docs-refresh PR |
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
| GEN.25 | A moon reclassified after its planet moves can be too large for its planet (bug) | none | open |
| GEN.26 | Rogue planet surface conditions | none | done, PR #263 (schema v48; design docs/design/rogue-planet-surface.md) |
| GEN.27 | Class P (glaciated world) only in the habitable zone, and fitting there | none | open |
| GEN.28 | Seven new planet classes in the letter gaps (R, S, U, W, X, Y, Z) | none | open |
| GEN.29 | Sweep every planet class for sense once the new ones are in | none | open |
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
| MAP.52 | Galaxy Map highlights the wrong area; pick a 40-degree wedge around the cursor (bug) | none | open |
| MAP.53 | Rotate a zoomed-in wedge, and zoom it to fit the window (bug) | none | open |
| MAP.54 | Slab leader lines instead of the slab slider (bug) | none | open |
| MAP.55 | Galaxy Map buttons: a menu, with only back, forward, up, reset and bookmark showing | none | open |
| MAP.56 | Drop the 3x3 block pick: select a slab, zoom in, select a segment (bug) | none | open |
| NAV.1 | Courses in "bearing mark mark" on nested frames | 28 (2026-09-30 18:14Z); 33 (2026-09-30 18:39Z to 2026-10-01 02:57Z) | done in 7.14.0, PR #130 (see note 4) |
| NAV.2 | Warp and fold speeds | 29 (2026-09-30 18:14Z); 34 (2026-09-30 18:39Z to 21:54Z) | done in 7.8.0, PR #121 |
| OPS.1 | Build the version number from the category counters (item 80's version-scheme questions) | none (split from 80 by the renumbering) | done in the version-from-todo-counters PR |
| OPS.2 | Apache OOM-killed on the production server | 1 (2026-09-24 01:32Z to 02:02Z) | done in 5.47.0, PR #72 |
| OPS.3 | PowerShell installers and macOS-safe bash scripts | 50 (2026-09-30 20:43Z to 2026-10-01 02:57Z) | done in 7.16.0, PR #125 (see note 4) |
| OPS.4 | Generate page jobs on native Windows | 54 (2026-09-30 20:48Z); 55 (2026-09-30 20:48Z to 22:12Z) | done in 7.9.2, PR #124 |
| OPS.5 | Install and update check the log locations and say how to fix them | none | open |
| PERF.1 | Generation at scale (new parent) | none | open |
| PERF.2 | Cache so pages don't hit the database every request | 6 (2026-09-24 01:32Z to 02:18Z); 3 (2026-09-24 01:57Z to 02:02Z); 1 (2026-09-24 02:25Z to 2026-09-30 18:09Z); 8 (2026-09-30 18:14Z to 2026-10-01 05:05Z) | done in 7.56.0, PR #178 |
| PERF.3 | Estimate size and time before bulk generation | 86 (2026-10-01 03:15Z to 05:29Z) | done, PR #238 (stats in control schema v6) |
| PERF.4 | Second progress bar for slow plan layers | 88 (2026-10-01 03:26Z to 05:29Z) | done, PR #258 |
| PERF.5 | Scatter bright stars in stages | 89 (2026-10-01 03:36Z to 05:29Z) | done, PR #229 |
| PERF.6 | Rate-limit SQL calls, do more per call | 90 (2026-10-01 03:46Z to 05:29Z) | done: investigation, then PR #222, #223 and #225 (PERF.8 caps the writers) |
| PERF.7 | Parallelize sector and system generation | 91 (2026-10-01 03:46Z to 05:29Z) | done, PR #225 and #227 |
| PERF.8 | Parallel background work queue in the API | 92 (2026-10-01 03:46Z to 05:29Z) | done, PR #225 and #227 |
| PERF.9 | Weight the bright-star ETA by the shape of the galaxy | 93 (2026-10-01 04:50Z to 05:29Z) | done, PR #258 |
| PERF.10 | Record generation speed across a log scale of densities | 94 (2026-10-01 04:58Z to 05:29Z) | done, PR #238 |
| PERF.11 | Store each sector's expected and actual density | 95 (2026-10-01 04:58Z to 05:29Z) | open |
| PERF.12 | Check the schema once per process during generation | none | done, PR #222 |
| PERF.13 | Write each sector in batches | none | done, PR #222 |
| PERF.14 | Reserve a sector's names in bulk, safe with several writers at once | none | done, PR #222 |
| PERF.15 | Fewer queries per web page | none | done, PR #223 |
| PERF.16 | Search names without scanning every row | none | done, PR #223 |
| PERF.17 | A time limit on web database statements | none | done, PR #223 |
| POP.1 | Government ownership of systems | 12 (2026-09-24 01:32Z to 02:18Z); 9 (2026-09-24 01:57Z to 02:02Z); 7 (2026-09-24 02:25Z to 05:38Z); 6 (2026-09-24 02:53Z to 2026-09-30 17:58Z); 15 (2026-09-30 16:44Z); 16 (2026-09-30 16:49Z); 17 (2026-09-30 16:51Z to 18:09Z); 34 (2026-09-30 18:14Z); 39 (2026-09-30 18:39Z to 18:41Z); 59 (2026-09-30 19:02Z to 19:17Z); 63 (2026-09-30 20:01Z to 20:27Z); 46 (2026-09-30 20:07Z); 50 (2026-09-30 20:08Z to 20:48Z); 51 (2026-09-30 20:43Z to 2026-10-01 04:37Z) | done in 7.49.0, PR #169 (optional since 7.58.2, PR #180) |
| POP.2 | Names for dominant species on living worlds | 13 (2026-09-24 01:32Z to 02:18Z); 10 (2026-09-24 01:57Z to 02:02Z); 8 (2026-09-24 02:25Z to 05:38Z); 7 (2026-09-24 02:53Z to 2026-09-30 17:58Z); 16 (2026-09-30 16:44Z); 17 (2026-09-30 16:49Z); 18 (2026-09-30 16:51Z to 18:09Z); 35 (2026-09-30 18:14Z); 40 (2026-09-30 18:39Z to 18:41Z); 60 (2026-09-30 19:02Z to 19:17Z); 64 (2026-09-30 20:01Z to 20:27Z); 47 (2026-09-30 20:07Z); 51 (2026-09-30 20:08Z to 20:48Z); 52 (2026-09-30 20:43Z to 2026-10-01 04:37Z) | done in 7.49.0, PR #169 |
| POP.3 | Database of spacefaring species | 14 (2026-09-24 01:32Z to 02:18Z); 11 (2026-09-24 01:57Z to 02:02Z); 9 (2026-09-24 02:25Z to 05:38Z); 8 (2026-09-24 02:53Z to 2026-09-30 17:58Z); 17 (2026-09-30 16:44Z); 18 (2026-09-30 16:49Z); 19 (2026-09-30 16:51Z to 18:09Z); 36 (2026-09-30 18:14Z); 41 (2026-09-30 18:39Z to 18:41Z); 61 (2026-09-30 19:02Z to 19:17Z); 65 (2026-09-30 20:01Z to 20:27Z); 48 (2026-09-30 20:07Z); 52 (2026-09-30 20:08Z to 20:48Z); 53 (2026-09-30 20:43Z to 2026-10-01 04:37Z) | done in 7.49.0, PR #169 |
| POP.4 | Younger and older civilizations | 15 (2026-09-24 01:32Z to 02:18Z); 12 (2026-09-24 01:57Z to 02:02Z); 10 (2026-09-24 02:25Z to 05:38Z); 9 (2026-09-24 02:53Z to 2026-09-30 17:58Z); 18 (2026-09-30 16:44Z); 19 (2026-09-30 16:49Z); 20 (2026-09-30 16:51Z to 18:09Z); 37 (2026-09-30 18:14Z); 42 (2026-09-30 18:39Z to 18:41Z); 62 (2026-09-30 19:02Z to 19:17Z); 66 (2026-09-30 20:01Z to 20:27Z); 49 (2026-09-30 20:07Z); 53 (2026-09-30 20:08Z to 20:48Z); 54 (2026-09-30 20:43Z to 2026-10-01 04:37Z) | done in 7.49.0, PR #169 |
| POP.5 | Population pages: species, polities, dominant species, territory (unnumbered in TODO.md) | none | done after 7.58.2, PR #184 |
| POP.6 | Territories overlay on the Galaxy Map (unnumbered in TODO.md) | none | done in 7.54.0 (PR #179) and 7.58.1 (PR #181) |
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
| USR.1 | User accounts (new parent) | none | open |
| USR.2 | Roles: user, admin, Owner | 64 (2026-10-01 02:13Z to 05:29Z) | open |
| USR.3 | SMTP settings | 65 (2026-10-01 02:13Z to 05:29Z) | open |
| USR.4 | Invite-only sign-up | 66 (2026-10-01 02:13Z to 05:29Z) | open |
| USR.5 | Email loop for passwords | 67 (2026-10-01 02:13Z to 05:29Z) | open |
| USR.6 | Owner transfer | 68 (2026-10-01 02:13Z to 05:29Z) | open |
| USR.7 | User-level interface with bookmarks | 69 (2026-10-01 02:13Z to 05:29Z) | open |
| UX.0 | Bugs and small fixes (standing item) | none | open while it holds bugs |
| UX.1 | Class reference pages | 56 (2026-10-01 01:19Z to 04:37Z) | done in 7.46.0, PR #167 |
| UX.2 | Menus sized to what they hold | 62 (2026-10-01 01:44Z to 05:29Z) | open |
| UX.3 | Warn visitors while a background job changes the galaxy | 87 (2026-10-01 03:26Z to 05:29Z) | open |
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
| UX.21 | Clean up the web interface: overlapping buttons and dead controls | none | open |
| UX.22 | Meaningful units for every measurement | none | open |
| VIEW.1 | View from a planet (new parent) | none | open |
| VIEW.2 | Starmap seen from a planet | 83 (2026-10-01 02:55Z to 05:29Z) | open |
| VIEW.3 | Render the view as a PNG with constellations | 84 (2026-10-01 02:55Z to 05:29Z) | open |
| VIEW.4 | Constellation names in the name generator | 85 (2026-10-01 02:55Z to 05:29Z) | open |

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
| TEST.4 | Revive and widen the known-bug tests | none | open |
| TEST.5 | Real 4 pc in boundary tests | none | open |
| TEST.6 | SQL portability lint | none | open |
| TEST.7 | Strict sql_mode on both engines | none | open |
| TEST.8 | Migrate from real old schemas | none | open |
| TEST.9 | Migration crash and re-run | none | open |
| TEST.10 | Database newer than the code | none | open |
| TEST.11 | Every column round-trips | none | open |
| TEST.12 | Boundary values round-trip | none | open |
| TEST.13 | Collation collisions | none | open |
| TEST.14 | CHECK constraints enforced | none | open |
| TEST.15 | Sector save fails halfway | none | open |
| TEST.16 | Id blocks after reset and rollback | none | open |
| TEST.17 | Batched writes at the limits | none | open |
| TEST.18 | Full-text search edge cases | none | open |
| TEST.19 | Same galaxy at any worker count | none | open |
| TEST.20 | Work queue failure paths | none | open |
| TEST.21 | Cancelling a run | none | open |
| TEST.22 | Every bulk mode in parallel | none | open |
| TEST.23 | Resume after an interrupted fill | none | open |
| TEST.24 | Bright-star scatter edge cases | none | open |
| TEST.25 | Interrupted bright-star scatter | none | open |
| TEST.26 | `--force` scatter then fill | none | open |
| TEST.27 | Progress and ETA under bad clocks | none | open |
| TEST.28 | CLI errors by message | none | open |
| TEST.29 | Limits stay consistent | none | open |
| TEST.30 | Grid seams and the nucleus | none | open |
| TEST.31 | Sector placement exhaustion | none | open |
| TEST.32 | System builder internals | none | open |
| TEST.33 | Moon stability helpers | none | open |
| TEST.34 | Kepler solver extremes | none | open |
| TEST.35 | Star and evolution helpers | none | open |
| TEST.36 | Phenomenon class helpers | none | open |
| TEST.37 | Names under parallel saves | none | open |
| TEST.38 | Population incremental rescans | none | open |
| TEST.39 | Navigation graph | none | open |
| TEST.40 | Two admins start a job at once | none | open |
| TEST.41 | Job files damaged | none | open |
| TEST.42 | Pages fresh after a CLI write | none | open |
| TEST.43 | Auth sweep over every route | none | open |
| TEST.44 | What an API key may do | none | open |
| TEST.45 | More than one admin | none | open |
| TEST.46 | Trusted device and TOTP edge cases | none | open |
| TEST.47 | Oversized requests | none | open |
| TEST.48 | Security headers everywhere | none | open |
| TEST.49 | Thin API routes | none | open |
| TEST.50 | Galaxy URLs combined | none | open |
| TEST.51 | Page-number sweep gaps | none | open |
| TEST.52 | Old URLs and error codes | none | open |
| TEST.53 | Formatters with bad numbers | none | open |
| TEST.54 | Caches under threads | none | open |
| TEST.55 | Map buttons do something | none | open |
| TEST.56 | No overlapping controls | none | open |
| TEST.57 | Galaxy Map JavaScript logic | none | open |
| TEST.58 | Other map JavaScript | none | open |
| TEST.59 | Galaxy Map drill-down in a browser | none | open |
| TEST.60 | Admin script command lines | none | open |
| TEST.61 | SQLite import script | none | open |
| TEST.62 | update.sh against a real database | none | open |
| TEST.63 | Math check that runs first (new parent) | none | open |
| TEST.64 | Reference values | none | done, PR #281 |
| TEST.65 | Identities and invariants | none | done, PR #281 |
| TEST.66 | Distributions match their targets | none | done, PR #281 |
| TEST.67 | Runs first in the suite and in CI | none | done, PR #281 |
| TEST.68 | Gate before bulk generation | none | open |
| TEST.69 | Intermittent failure in the colony test | none | open |
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
