**Headline:** The All systems list opens any page as fast as the first when sorted by sector, octant or binary, as it already did by name.

### Changed
- **The All systems table pages by key in every sort.** Sorted by sector, octant or binary, the page links carry the system the page before ended on and the list carries on from it inside its group, so a page deep in the list costs what the first does. Before, a page counted and skipped every system before it in its group: with 3,000,000 systems, a page 300,000 deep in an octant took 1.5 seconds and now takes 3 milliseconds. `GET /api/systems?after=<system ID>` now works in those sorts as well as by name.
