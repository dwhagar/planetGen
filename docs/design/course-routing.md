# Course routing

How the NAV page finds a route from one system to another, today and as
planned. Recorded 2026-10-02 from Boss's decision on hop length (NAV.12,
01:53Z) and the hop-length study (report in the project's shared files
under `nav-hop-length/`). The frames a course is described in (bearing
and mark, warp and fold speeds) are in
[navigation-frames.md](navigation-frames.md).

**Status (2026-10-02, checked against main after PR #359):** section 1
is built, and so is NAV.38's line-to-sectors helper. Everything else
from section 2 on is planned and not built; each piece names its TODO
item and phase.

## 1. Today

- `queryDb.nav_between` picks a scope. Two systems in the same sector
  route over that sector's own systems only (the Sector Local Frame);
  anything else routes over every placed system in every generated
  sector (`_galaxy_frame_positions`, the Galactic Standard Frame).
- `navGraph.build_knn_adjacency(positions, k=6)` builds a k-d tree and
  links each system to its 6 nearest, made two-way. It is rebuilt on
  every request.
- `navGraph.shortest_path` runs Dijkstra on straight-line hop lengths.
- Only star systems are stops between the ends. Phenomena can be the two
  ends but never a stop; bright stars from the galaxy scatter are not
  stops; ungenerated sectors add nothing.
- Nothing limits hop length and nothing reports it: the route returns
  its stops and total length. Warp and fold times are given for the
  direct distance only.
- The NAV page lists the route as a vertical list (`<ol class="nav-route">`
  in `nav.html`) and draws it on one square map at one scale
  (`lib/navmap.py`).

What goes wrong (from the study):

- **Islands.** A separately generated area with 7 or more systems had no
  link out, because each of its systems' 6 nearest are inside it. 2,000
  generated sectors spread over the disk split into 714 islands, and a
  course across the galaxy found no route ("No route via adjacent
  systems could be found"). Fixed by NAV.34 (PR #427):
  `navGraph.join_islands` links each island to its `NAV_ISLAND_LINKS`
  (6) nearest islands, in rounds until one island is left, so a route
  always exists.
- **Hidden long hops.** A lone system or a small cluster always links
  out, however far: a system 2 kpc above the disk was reached by a route
  whose last hop was 6,504 ly.
- **Speed.** The rebuild dominates: 0.5 s at 10,000 systems, 16 s at
  200,000. Dijkstra itself takes 0.7 s at 200,000.

## 2. No maximum hop length (NAV.12, phase 1)

Boss (2026-10-02 01:53Z):

> routs will always find the nearest star they can even if it crosses
> sector boundaries even across multiple sectors. This is for a game
> mechanic I need in place and I also want a jump through unknown space
> is marked in red and glows to draw attention to it. Also UX change
> here to display the path horizontally and find a way to split it
> between multiple lines for mobile or limited displays.

This drops the study's optional "ship range" field; there is no hop
limit of any kind.

- **A route always exists** between two placed endpoints. The islands
  of the route graph are joined, each linked to its nearest few islands
  (NAV.34, built in PR #427). In the study, linking each island to its 6 nearest
  joined all 714 in 0.2 s and gave a cross-galaxy route 1.29 times the
  direct distance.
- **Same-sector routes may leave the sector**, so the nearest star is
  used whichever sector it is in.
- **The longest hop is shown** with the route.
- **Unknown-space jumps.** A hop is flagged when its straight line
  crosses one or more unfilled (ungenerated) sectors; that is the
  default reading of "unknown space". The sectors a line crosses come
  from `galaxyGeometry.sectors_along_segment` (NAV.38, built in PR #357),
  exact at sector faces and edges, with its browser twin
  `galaxyprisms.sectorsAlongSegment`. `/api/nav` returns the flag per hop.
- Route edge cases from the study (islands, lone systems, the halo, the
  dense core, endpoints that are phenomena) are written as tests first
  (TEST.79, built in PR #427).

## 3. Routing that scales (NAV.10, phase 1)

Built with NAV.12, which changes the same code.

- The search loads only the systems in a corridor around the direct
  line (a box query on indexed sector centers, widened when no route is
  found), or reads the stored `nearest_systems` table instead of
  rebuilding the graph.
- A* with the straight-line distance as its heuristic.
- A query for every system, star and phenomenon within a distance of a
  line segment (used later by NAV.6 to steer around gravity wells).
- A galaxy schema migration adds the position indexes.
- Timed on a 100,000-sector database.

Travel times for the route itself, per hop (unknown-space jumps
included) and in total, with the same warp and fold tables as the
direct distance, follow (NAV.11, phase 1), the total assuming a stop at
every system on the route (Boss, 2026-10-02 04:19Z).

## 4. Showing the route

- **A readable course map (NAV.41, phase 0, bug).** The NAV page's map
  (`lib/navmap.py`) becomes a wide panel across the usable width of the
  device, with labels at least body-text size, and the route list stays
  below it (Boss, 2026-10-02 04:19Z).
- **Course and distance per stop (NAV.42, phase 1).** Each stop shows the
  course to the next in the existing notation, for example
  `045 mark 012, 3.2 ly`, worked out in that hop's frame.
- **Horizontal, wrapping (UX.35, phase 1, alongside NAV.12).** The stops
  run left to right, each a link, with the hop distance between them.
  On phones and narrow panels the strip wraps onto several lines (a
  container query, not a device check) and never splits a stop across
  lines. Screen readers still get an ordered list.
- **Red and glowing (NAV.36, phase 2).** A flagged hop is drawn red with
  a glow in the route strip, on the NAV page's map and on the Galaxy Map
  course, with a legend entry. Under `prefers-reduced-motion` it stays
  red without the pulse; it reads in both themes; the route list also
  labels it in text, so it isn't shown by colour alone.
- **Saved courses (NAV.39, phase 2).** A saved course (NAV.17) keeps
  which hops were unknown-space jumps, and opening it checks again,
  since sectors may have been filled since; a hop that is now known
  shows as ordinary and the course says what changed.

## 5. Order

| ID | Piece | Phase | Needs |
|---|---|---|---|
| TEST.79 | Route edge cases, written first | 0 | **Built, PR #427** |
| NAV.34 | Join the route graph's islands (bug) | 0 | **Built, PR #427** |
| NAV.38 | Every sector a straight line passes through | 0 | **Built, PR #357** |
| NAV.10 | Corridor search, A*, position indexes | 1 | queues behind PERF.11's schema migration |
| NAV.12 | No hop limit, longest hop shown, unknown-space flag | 1 | NAV.34, NAV.38, TEST.79, NAV.10 |
| UX.35 | Horizontal, wrapping route strip | 1 | alongside NAV.12 |
| NAV.11 | Travel times per hop and for the route | 1 | NAV.10, NAV.12 |
| NAV.36 | Unknown-space jumps red and glowing | 2 | NAV.12, UX.35, NAV.20 |
| NAV.39 | Saved courses keep and re-check their unknown-space jumps | 2 | NAV.12, NAV.17 |

## 6. Why it works this way

- **No cap.** With most of the galaxy ungenerated, every capped test in
  the study that crossed ungenerated space found no route, at every cap
  from 5 ly to 5,000 ly, and the answer would change as more sectors are
  generated. At top speed (warp 9.995, fold 8.5) even a hop across the
  whole galaxy takes years, so no hop is impossible under the current
  rules.
- **Why the route has stops at all.** With straight-line costs the
  direct edge is always shortest; the route has stops only because the
  graph leaves long edges out. The 6-nearest rule is today's implicit
  limit (about 8 ly at the Sun's radius, 34 ly 1 kpc above the disk).
- **Flag instead of forbid.** A long hop across ungenerated space hides
  a gap behind a list that looks like a chain of nearby stars, and it
  makes NAV.6's gravity-well steering costlier. Flagging it tells the
  user without refusing the course.
