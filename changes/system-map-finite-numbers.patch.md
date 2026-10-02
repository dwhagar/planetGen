### Fixed
- **The System Map never draws a broken number (MAP.57).** A NaN or infinite value stored for a star, planet, moon, belt or facility no longer ends up in the map's SVG. A planet or moon whose position was lost is drawn at its orbit distance, due east of what it orbits, with a note in its info panel; one with no distance either is left out.
