# Lane designations

Boss (5:31 PM EST on 2026-11-07): every active thread carries a number and its job. Docs, items and labels say "Lane NN" (for example "Lane 11"); the thread title carries the same number and job. Numbers are stable and are never reused. This page maps the numbers to the old names, so history stays traceable: older text in the design notes, reports, memory and pull requests still uses the old names.

Numbering: 0x runs the system, 1x is core generation and product families, 2x is composition and physics lanes.

## Active lanes

| Lane | Job | Old name | Notes |
|---|---|---|---|
| 01 | TODO: items, issues, docs, merges | TODO | files and retires items, merges non-failing pull requests |
| 02 | Train Manager: merge trains | Train Manager (was Bug 2) | builds and merges the trains |
| 03 | Gate: run gates, fix reds | Gate 1 | started 2026-11-07, runs gate shards and fixes reds for lane 02 |
| 04 | Bugs and test health | Bug 1 | the current Bug 1 thread (an older Bug 1 thread ran in October 2026 and is retired) |
| 05 | Bugs and test health, second lane | none (new 2026-11-08) | Boss, 8:10 PM EST 2026-11-07: a second bug lane to assist lane 04 with the 8.1 bug-hunt queue (GEN.1390, TEST.695 and their splits) |
| 11 | Oddities, census, ore design | Core 1 | |
| 12 | Storage, remnants, volatility law | Core 2 | |
| 13 | Galaxy and admin | Core 3 | |
| 14 | Atmosphere classes and sweeps | Core 4 | |
| 15 | Air, climate and performance | Core 5 | |
| 21 | Planet composition chain | Planets 1 | |
| 22 | Belts, small bodies, mineral rows | Aux 5 | the current Aux 5 thread (an older Aux 5 ran in October 2026 and is retired) |
| 23 | Tides, spin, moons | Tides 1 | |
| 24 | Materials | none yet | reserved for the real-material items that wait on their prerequisites (GEN.1377, GEN.1378, GEN.1381 to GEN.1383, TEST.686; GEN.1380 is with lane 22) |

## Retired lanes (no number; the old names stay)

Aux 1, Aux 2, Aux 3, Aux 4, Aux 6, Aux 7, Aux 8 to Aux 16, Sources 1, Sources 2, Star 1, Trojans, Math 1, the units lane, the search-filter lane, the admin redesign thread, and the older Bug Hunt, Bugfix, Blocker, Habitability, Composition, Map nodes and Foundations lanes listed in `lanes.md`. Where an old name equals a number's old name (Core 1 to Core 5, Bug 1, Aux 5), the text in older files means whichever thread held the name then; the number above is the current one.

## How the rest of the repository reads

- `docs/TODO.md`, `docs/plan/` and the generated pages use the new numbers in owner lines, hand-offs and notes (rewritten 2026-11-07, 6:35 PM EST).
- Design notes under `docs/design/`, reports, handoffs and memory files keep the old names; read them through this table.
- GitHub carries no lane labels, so none were renamed; goal issues keep their numbers.
