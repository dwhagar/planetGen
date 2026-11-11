**Headline:** The solids a giant planet accretes now come from the whole disk outside it, so a giant that formed warm still gets argon, krypton and xenon from the planetesimals that froze them.

### Changed

- `disk_condensation.provenance_condensed` (GEN.1519): the planetesimals in a giant's envelope are the mean of the species frozen over the annuli between its place and the disk's edge (equal per AU), not only what froze at its own place. `envelope_c_over_o` and `envelope_enrichments` take it (`solids_edge_k`); `mixture_giant` sets it to the disk's edge temperature.
- Giants that formed warmer than a noble gas's line now hold it: on a grid of stars, argon 1 to 2, krypton 1.7 to 2.5 and xenon 2.9 to 3.6 times the star's (they were 1 before this stack and 0.03 to 0.05 with the tapered pebbles alone; Jupiter has 2.5, 2.7 and 2.6). Nitrogen and the median carbon over oxygen are unchanged; sulfur is 3.1 to 3.3 (probe 2.75). Stand-ins PRV1 to PRV3 are registered with the falsifiers.

generation-output: changed
generation-affects: fill
