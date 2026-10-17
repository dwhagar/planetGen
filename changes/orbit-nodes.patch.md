**Headline:** A built star system can store its tree of orbits, its stars' places and velocities, and its disks.

### Changed
- **Orbit tree storage (DB.53).** A built multi-star system gets `orbit_nodes` (the same columns as the scatter's nodes, plus a disk descriptor), each star gets its node, offset, velocity, state, flags, decay age and epoch, and a circumbinary planet gets its host node. Existing close and wide pairs convert to one node each from their stored orbit. Migration 0104.

TODO items resolved: DB.53. Added: none.
