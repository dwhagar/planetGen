### Fixed
- **An asteroid belt could overlap the next planet out.** When spacing
  pushed a planet into the habitable zone and it was reclassified there,
  some classes redrew its distance anywhere in the zone, which could put
  it back inside the belt it had just been moved past. It happened most
  often around giant stars, whose habitable zones are wide, and made
  `test_systems.py` fail at random. A reclassified planet now keeps the
  distance it was moved to.
