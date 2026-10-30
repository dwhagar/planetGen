**Headline:** The irregular moons of the four giants are checked for saturation: they end within a factor 2.5 of each other, and the repo's belt law is shown not to fit a moon zone.

### Changed
- **Saturation check of the irregular zones (GEN.1079).** New tests pin how far the grinding of `physics/irregular_population.py` saturates (counts above 5 km and above the storage floor and the mass left within 2.5 to 3 times across Jupiter, Saturn, Uranus and Neptune), that the stored counts match the observed 8, 5, 6 and 8 to within 1.5 times, that a single scalar on the collision probability cannot fix the 5 km counts without breaking those, and that `belt_mass_earth` leaves 400 to 4400 times too little for a zone round a planet. No scalar is added. Test-only.

TODO items resolved: GEN.1079.

No generator output changes.

generation-output: unchanged
