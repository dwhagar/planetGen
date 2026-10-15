**Headline:** Every object now has a 96-bit ID, every scatter point moves with the galaxy, and the grid reaches 262 kpc.

### Changed
- Object IDs are 96 bits (GEN.230): 19 bits of shell, 18 of band, 19 of wedge, 28 of serial and 12 of body. Every `uid` column is `BINARY(12)`, the sector designation is 56 bits, and the layout is fixed for the whole galaxy instead of fitted to the outline. Migration 0093 renumbers every object in place, in chunks, and can be resumed.
- The sector grid now reaches shell 65,535 (262 kpc, GEN.233). The disk outline (`GalaxyBounds.qualifies`) still limits what `generate` builds on its own; `GalaxyBounds.contains` answers for the whole reach. `generate --reach-shells` sets the reach.
- The nightly orbit update moves every scatter point (GEN.223): bright stars and every kind of phenomenon scatter (the nucleus stays at the center). Points on the rotation curve turn about the galactic axis, hypervelocity stars fly a straight line, and a point that crosses a cell boundary is refiled there. Each point keeps its position, velocity, `epoch_unix` and `next_update_due` in its own columns.
- New tables `scatter_members` and `scatter_nodes` hold a multi-star point's members and orbit tree (GEN.260), with the columns the system-type catalogue needs (member state and flags, node interaction, a group's offsets and velocities, up to 12 members), and the scatter rows carry `member_count` and `system_mass_solar`; the point is the barycentre.
- Every table that holds a point in the galaxy frame gets a `cell_key` (DB.34, DB.35): a generated column with the Morton key of the object's 0.25 pc cube, and a covering index `(cell_key, x, y, z)`, the base of the region queries. A star system also keeps its point in the galaxy frame (`galaxy_x/y/z_pc`) now, since its own position is relative to its sector. `galaxy/cellkey.py` is the same function in Python.

TODO items resolved: GEN.223, GEN.230, GEN.233, DB.34, DB.35, GEN.260 (tables and columns). Added: none.
