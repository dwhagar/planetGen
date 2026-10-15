**Headline:** Every object now has a 96-bit ID, every scatter point moves with the galaxy, and the grid reaches 262 kpc.

### Changed
- Object IDs are 96 bits (GEN.230): 19 bits of shell, 18 of band, 19 of wedge, 28 of serial and 12 of body. Every `uid` column is `BINARY(12)`, the sector designation is 56 bits, and the layout is fixed for the whole galaxy instead of fitted to the outline. Migration 0093 renumbers every object in place, in chunks, and can be resumed.
- The sector grid now reaches shell 65,535 (262 kpc, GEN.233). The disk outline (`GalaxyBounds.qualifies`) still limits what `generate` builds on its own; `GalaxyBounds.contains` answers for the whole reach. `generate --reach-shells` sets the reach.
- The nightly orbit update moves every scatter point (GEN.223): bright stars and every kind of phenomenon scatter (the nucleus stays at the center). Points on the rotation curve turn about the galactic axis, hypervelocity stars fly a straight line, and a point that crosses a cell boundary is refiled there. Each point keeps its position, velocity, `epoch_unix` and `next_update_due` in its own columns.
- New tables `scatter_members` and `scatter_nodes` hold a multi-star point's members and orbit tree (GEN.260), and the scatter rows carry `member_count` and `system_mass_solar`; the point is the barycentre.

TODO items resolved: GEN.223, GEN.230, GEN.233, GEN.260 (tables). Added: none.
