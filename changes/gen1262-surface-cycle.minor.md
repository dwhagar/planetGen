**Headline:** One calculation now gives the sunlight on every place of a body through its year and its turns, and the day and night search of the latitude bands runs on it.

### Added

- `physics/surface_cycle.py` (GEN.1262): the cosine of the angle between the local vertical and the star for each cell of a 12 by 12 grid of equal-area cells and each step of time (a 48th of a solar day, a 3000th of the year or a quarter of the surface's thermal time constant, whichever is least, up to 20,000 steps), from the body's stored tilt, rotation period and orbital period alone. A locked world, a world that turns backwards, a 28-day spinner near resonance and Earth all come out of the same formula; an optional transmission (the ring's shade) multiplies the light where the star is up.

### Changed

- The search for the longest unbroken day and night and the mean day and night of the latitude bands (GEN.1260) now takes its geometry from the new module. The stored values are the same as before to the last digit; only the code path moved.
