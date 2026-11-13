**Headline:** when the life of a world changes its temperature or pressure, the water is now made again at the new air, so a world the life heated past the boiling point no longer keeps a liquid ocean, and the class is read after the thaw path too.

### Changed

- `life._rebuild(planet, remake_water)`: when the life pass moved the surface temperature or the pressure, the hydrosphere, ice, ocean cover and ground are made again with the existing `planets.refresh_at_final_age` (which also makes the habitability score) on the body's own stream, before the class is applied. A world the life did not move keeps what it had, draw for draw. A liquid ocean stored at 432 K under 63 kPa becomes vapor and a dry surface.
- The thaw path (`_remake_ground_if_thawed`) now applies the class after its refresh, so the class no longer reads the cover the thaw replaced.
- New `scripts/census_stale_water.py` counts bodies of a class census that store a liquid above its boiling point (reports only). Not done here: reading the liquid cover instead of the ice-or-liquid maximum where the class asks for ocean cover (GEN.1567 owns the liquid and ice split).
- Changes the stored water and class of worlds whose life moved their air, so it is a lock item: the golden and epoch are the train's.

generation-output: changed
generation-affects: fill, population
