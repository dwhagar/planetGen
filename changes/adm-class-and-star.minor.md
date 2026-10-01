### Added
- **Change a planet's or moon's class, or a system's star (ADM.6,
  ADM.7).** On the system page's Edit panel an admin can pick a new class
  for any planet or moon: the menu lists the classes that fit where it
  is without moving any other planet first, and every other class under
  "Force", which is kept even where it couldn't form. A single-star
  system can have its star replaced by any spectral type: planets,
  moons and belts keep their classes and their orbits scale with the new
  star's light, and the page lists anything that had to move or that no
  longer had a stable orbit and was removed. New API routes: `GET
  /api/systems/<id>/class-options`, `POST /api/planets/<id>/class`,
  `POST /api/moons/<id>/class` and `POST /api/systems/<id>/star`.
