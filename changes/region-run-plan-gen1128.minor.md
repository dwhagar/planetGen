**Headline:** The plan of a region run is written: a sphere around a start point becomes whole sectors, cut into units along the galaxy's own shell and band arcs, nearest first, with resume and the stage barrier.

### Added
- **The region run plan (GEN.1128, first slice).** `planetgen/galaxy/region_run.py`: `plan` (a sphere of a radius from a start point, plus `grow_shells` whole shells, as the whole sectors it touches; units of at most 27 consecutive wedges of one (shell, band) arc; nearest first; sectors already charted taken out so a stopped run or an overlapping region builds only what is missing), `next_units` (every generate unit before any top-up unit, so a stopped run resumes in the generate stage until the whole set is charted), `unit_items` (a unit as the triples a sector stage takes) and `cost` (sectors, units, shells, the sphere's volume in sectors, seconds only when a per-sector time is given). The plan reads no stream and no database, and no draw reads it, so the same seed gives the same sectors whatever the start, the radius, the cut or the order. The command that runs the units and the chart ledger write follow (GEN.1127).

No generator output changes.

generation-output: unchanged

TODO items resolved: none (GEN.1128 stays open for the command).
