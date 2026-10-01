### Fixed

- System page, "Place a facility" (ADM.9): the form asks for the name,
  then the placement, then the host, and the host list holds only the
  stars, planets, moons or belts that take that placement.
- An orbital facility's distance is now a logarithmic slider from just
  above the host's surface to the edge of its sphere of influence (a
  planet's or moon's Hill sphere, a star's heliosphere), with the
  distance, period and speed read out as it moves. It shows only for "in
  orbit"; a surface facility has no distance control. The API refuses an
  orbit outside the host's sphere of influence.
- A facility in an asteroid belt gets a random spot in the belt and a
  circular orbit around its star from there, and `updateOrbits.py` moves
  it along with orbital facilities. No schema change; belt facilities
  saved before this have no orbit and stay where they are.
