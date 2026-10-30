**Headline:** The rule for where a captured asteroid can stay in orbit round a planet now exists as a tested calculation, with the Hill radii and the inner edges of the four giants as checks.

### Changed
- **The irregular zone (GEN.1076, first part).** A new module (`physics/irregular_zone`) gives, for a planet and its star, the Hill radius (0.35, 0.43, 0.47 and 0.77 AU for Jupiter, Saturn, Uranus and Neptune), the distance out to which a prograde orbit (0.4 of it) and a retrograde one (0.55 of it) stay, the inclination band where the star drives a body into the planet, and the inner edge set by the outermost regular moon (0.013 AU at Jupiter, 0.023 at Saturn, 0.004 at Uranus). A hot giant's zone comes out empty. This is the first of the pieces for moons that were once asteroids; nothing reads it yet.

generation-output: unchanged
