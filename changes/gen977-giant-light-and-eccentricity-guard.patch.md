**Headline:** A new test keeps the generator from making a giant whose radius is read at the wrong light.

### Changed
- **Guard for giants that are bright and eccentric (GEN.977).** A giant's radius is read from the tables at the light of its orbit's semi-major axis, because the eccentricity is drawn after the radius. That is within 1% of the radius at the orbit's mean light unless the giant gets more than 300 times Earth's light and an eccentricity over 0.3 together. A seeded sweep now fails if the generator ever makes one, and the bound is written into the stand-in registry (OE12). Nothing you see on a planet changes.
