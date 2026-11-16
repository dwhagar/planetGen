**Headline:** a sector being generated now shows a count while it works, so a run's bar no longer stands still for the whole time a sector takes.

### Changed

- OPS.266: a sector is one piece of its unit, and its stars, planets and moons were drawn inside that piece with nothing moving the unit's count (18 to 42 s for a sector of one system on a small galaxy, minutes for a crowded one). The sector's second-level count is now the system under way ("System 3 of 40: bodies made", "Pre-placed star 2 of 5: bodies made") counting each planet and moon as it is made, the save of the sector is a count of its own ("Saving the sector"), and every timed phase of generation that ends starts the quiet clock again. The same holds for `planetgen sector` and `planetgen system`. Measured on the silence test's small galaxy: the longest quiet stretch of a sector unit fell from 18 to 42 s to 3 to 7 s. The galaxy is the same.
