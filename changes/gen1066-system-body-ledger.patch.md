**Headline:** A system can now be checked and summed as one list of every body with a mass, including a barycentre that has the belts in it and a measure of how much the belts move the star.

### Changed
- **The system body ledger (GEN.1066).** `physics/system_bodies.py` holds the pure half of the plan to see a whole system at once: one row per star, planet, moon, belt, belt member, remainder of a belt, comet and ring (mass, radius, density, orbit, status, tier, and a place relative to its parent). It checks that the rows are one tree (one root, parents of an allowed kind, no cycle, no mass on a belt itself), places every body in the system frame, sums the mass of a body and everything that orbits it, compares a belt's rows with the belt's own mass, and takes the barycentre with the belts in it. The belt wobble is how far the belts move the star from where the planets alone put it. The design note decides view against table: a view over the existing tables, no new storage. Nothing reads it yet; the view and its pages are the wiring step.

generation-output: unchanged

TODO items resolved: none (GEN.1066 stays open until the view is wired).
