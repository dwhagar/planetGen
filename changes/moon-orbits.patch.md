### Fixed
- **Moons could orbit outside their planet's Hill sphere, or inside the
  planet itself.** Moons now orbit between the planet's surface (plus
  room for the largest moon it can hold and its atmosphere) and the
  prograde stability limit of about half the Hill radius (Domingos,
  Winter & Yokoyama 2006). A planet whose class is regenerated after a
  move drops the moons that no longer fit (TODO items 40 and 41).
