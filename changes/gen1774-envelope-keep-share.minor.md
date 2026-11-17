**Headline:** A rocky planet of 1 to 3 Earth masses now keeps its hydrogen envelope only some of the time, a share that rises smoothly with the mass, so an Earth-mass planet in the habitable zone is bare again instead of keeping a thin envelope with no ground.

### Changed

- GEN.1774: `envelope_capture.KEEPS_FROM_MASS_EARTH` (row EC10) was a step: every body of 1.0 Earth masses or more kept the envelope it drew and none under 1.0 did, so the Earth twin (1.00003 Earth masses) kept 0.3 percent, was typed transitional with a radius of 8,678 km and classed T, in every seed. The share of bodies that keep their envelope is now `envelope_capture.keep_share`, the smooth step from 0 at 1 Earth mass (`KEEPS_FROM_MASS_EARTH`) to 1 at 3 (`KEEPS_ALL_FROM_MASS_EARTH`), decided by one draw from a stream of the body's own (`keeps_envelope`, no draw of the run's); a body that does not keep it is bare, as a body under 1 Earth mass is. A tweak of EC10 with no source for 1 to 3 Earth masses (tweak L4), not a new modifier. Bodies of 3 Earth masses or more and of 1 or less come out as they did. Reuse: `envelope_capture.start_stream` key style, `planets.apply_envelope`, the Class M tier's `bare_share` dial (already bare under a tier).

generation-output: changed
