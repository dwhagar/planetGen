**Headline:** the orbit ring drawn for a moon in the 3D system view now lies in the plane the moon moves in.

### Changed

- MAP.425: a moon's stored inclination and ascending node are measured in its planet's equator (GEN.1352), but its stored position and velocity are in the star's frame. The system scene took the ring's plane from the columns and the moon's place from the state, so the highlighted ring of a moon (Agnowee VIb, VIc and VIg) did not pass through it. The scene now takes the plane from the state too, as it already took the eccentricity and the periapsis. Display only; nothing stored changes.
