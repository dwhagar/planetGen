**Headline:** A scattered multi-star system moves as one point and its members are placed from their stored orbits.

### Changed

- A scattered system's barycentre is its scatter point and moves with the nightly update like every other point; the update writes no member or node row.
- A member's position at any time is the barycentre plus the offsets of the stored orbit tree, each orbit's phase read from its epoch and period (`planetgen/physics/system_tree.py`).
- `planetgen/db/scatter_systems.py` stores, reads and places a drawn tree under a `bright_stars` or `phenomenon_scatter` row, and gives the mass a point carries for the influence rule (the system's mass where it has members).

TODO items resolved: GEN.265. Added: none.
