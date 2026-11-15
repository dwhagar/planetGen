# Lane designations

Boss (5:31 PM EST on 2026-11-07): every active thread carries a number and its job. Docs, items and labels say "Lane NN" (for example "Lane 11"); the thread title carries the same number and job. Numbers are stable and are never reused. This page maps the numbers to the old names, so history stays traceable: older text in the design notes, reports, memory and pull requests still uses the old names.

Numbering (Boss, 8:46 and 8:48 PM EST on 2026-11-14): 00 is TODO; 01 to 08 are the generator lanes (lock items that ride generator trains), 01 being Train Manager; 12 and up are non-lock lanes (job engine, tests, gate, storage, admin). Before that date the scheme was 0x system, 1x core, 2x composition; text written before 9:00 PM EST on 2026-11-14 uses the old numbers.

## Active lanes

| Lane | Job | Number before 2026-11-14 | Old name |
|---|---|---|---|
| 00 | TODO: items, issues, docs, merges | 01 | TODO |
| 01 | Train Manager: merge trains | 02 | Train Manager (was Bug 2) |
| 02 | Air and climate | 15 | Core 5 |
| 03 | Atmosphere classes and sweeps | 14 | Core 4 |
| 04 | Planet composition chain | 21 | Planets 1 |
| 05 | Oddities, census, ore design | 11 | Core 1 |
| 06 | Belts, small bodies, mineral rows | 22 | Aux 5 |
| 07 | Tides, spin, moons | 23 | Tides 1 |
| 08 | Materials (reserved) | 24 | none |
| 12 | Gate: run gates, fix reds | 03 | Gate 1 |
| 13 | Job engine (was Gate 2) | 06 | none (new 2026-11-08) |
| 14 | Bugs 2: bugs and test health | 05 | none (new 2026-11-08) |
| 15 | Storage and volatility | 12 | Core 2 |
| 16 | Galaxy and admin | 13 | Core 3 |

Closed 2026-11-14: old 04 (Bugs and test health, Bug 1) and old 07 (Dev, builds ahead of trains); their open items went to Lane 14. Numbers 09 to 11 are unused.

## Retired lanes (no number; the old names stay)

Aux 1, Aux 2, Aux 3, Aux 4, Aux 6, Aux 7, Aux 8 to Aux 16, Sources 1, Sources 2, Star 1, Trojans, Math 1, the units lane, the search-filter lane, the admin redesign thread, and the older Bug Hunt, Bugfix, Blocker, Habitability, Composition, Map nodes and Foundations lanes listed in `lanes.md`. Where an old name equals a number's old name (Core 1 to Core 5, Bug 1, Aux 5), the text in older files means whichever thread held the name then; the number above is the current one.

## How the rest of the repository reads

- `docs/TODO.md`, `docs/plan/` and the generated pages use the new numbers in owner lines, hand-offs and notes (rewritten 2026-11-07, 6:35 PM EST).
- Design notes under `docs/design/`, reports, handoffs and memory files keep the old names; read them through this table.
- GitHub carries no lane labels, so none were renamed; goal issues keep their numbers.
