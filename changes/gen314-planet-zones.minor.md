**Headline:** One function now works out where a planet can stay in a hierarchy of stars and how much light it gets there.

### Added
- **Planet zones in a hierarchy, and the flux a planet receives from every member (GEN.314).** `physics/planet_zones.py` gives an S-type window round every member (the Holman and Wiegert limit or the Hill radius against the sibling subtree as a point mass, whichever is tighter, from the member's engulfment floor) and a P-type window round every pair node (the circumbinary fit or the highest floor of a member below it, up to the pair's own S-type limit against its sibling), all capped at 100 AU (`PLANET_ZONE_DISK_CAP_AU`). For two stars the windows are the ones `Star.a_crit_au` and `StarSystem._orbit_floor_au` hold. `flux_at` sums the host's own light and every other member's at the time-averaged inverse square of the node they meet at. On a seeded draw of the generator's trees the share of stars with an S-type window above 0.3 AU falls from 73% at two stars to 15% at seven, the design table's row within 7 points. GEN.264 stage C, GEN.130 and GEN.275 call it and work out no zones of their own.

TODO items resolved: GEN.314. Added: none.
