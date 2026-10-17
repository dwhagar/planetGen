**Headline:** The galaxy-type and planet-class items no longer redo work the galaxy laws and the composition chain already build.

### Changed

- GEN.380, GEN.382 and GEN.387 (galaxy types) now extend GEN.400, GEN.401, GEN.402 and GEN.403 instead of drawing the masses, the shape, the black hole and the cluster count a second time; GEN.402 and GEN.403 take the spheroid mass and a family factor as inputs from the start; GEN.164's older cluster count is superseded by GEN.403's halo-mass law; GEN.386 is the family scaling of GEN.405 and GEN.406.
- The planet classes (GEN.350 to GEN.356, GEN.27, GEN.145, GEN.29) wait for GEN.428, GEN.411 and the gas table and draw their albedo, mixture and composition priors through the composition chain; GEN.29 becomes the judgement pass on top of TEST.209 and TEST.210; the four composition proofs share one fixture on TEST.186's harness.

TODO items resolved: none / Added: none
