**Headline:** The to-do list places the capture and ejection items in the Phase 2 and Phase 3 plan.

### Changed

- GEN.289 to GEN.293 (planets and moons that can be captured or ejected) now carry Boss's design of 11:06 PM and 11:10 PM EDT on 2026-10-15: every planet and moon is the universal point-in-space object with its point mass and movement vector stored, and each step sums the point-mass pull of the host, the system's own planets and moons and passing bodies. GEN.289, GEN.109, GEN.290 and GEN.291 run in that order in Phase 2 (GEN.105 now waits on GEN.289 to GEN.291); GEN.292 and GEN.293 move to a new Phase 3 thread. No lane has them yet.
- Filed GEN.294 (Boss, 11:12 PM EDT): every object in the galaxy is held as one SpatialPosition3D with its position, point mass and movement vector stored as columns. It is a prerequisite of GEN.109 and shares one migration with GEN.289.

TODO items resolved: none / Added: GEN.294
