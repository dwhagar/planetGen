**Headline:** Planets and moons now travel on real ellipses, with seasons from the changing distance to their star.

### Changed
- **Planets and moons have an eccentricity and a mutual inclination (GEN.727).** A planet's eccentricity is drawn from the
  Kepler laws (Van Eylen et al. 2019: a half-Gaussian of width 0.083 for a system of several planets, a mixture peaking
  near 0.26 for a lone one), a regular moon's forced eccentricity from 0.001 to 0.01 (kept for good, as Io's is by its
  neighbours) and a captured moon's from the thermal law, all cut where the orbit would leave the stable region (periapsis
  outside the Roche limit, apoapsis inside the Hill sphere). Tidal circularisation over the star's age pulls close planets
  and captured moons back toward a circle. Inclination is Rayleigh about a
  plane shared with the others (mutual inclination 1.6 degrees over the square root of two).
- **The state vector is the ellipse.** Position and velocity come from the orbital elements, the stored distance is the
  semi-major axis, and the orbit update steps the state on its conic. The new `eccentricity` column (planets and moons,
  migration 0156) is stored; older rows read as a circle.
- **Light follows the orbit.** A body's mean flux is the circular value over the square root of one minus e squared, its
  temperature range gains the swing between periapsis and apoapsis, and a lone star's `flux_variation` is that swing.
- **Spacing, belts, moon stability and the system map use periapsis and apoapsis.** Orbit reach, the gap past a belt and
  between planets, and the moon stability check use the ellipse; the system map draws an ellipse for each orbit and the
  object page has an Eccentricity row.
- **Stand-ins recorded.** The drawn laws and their source tiers and tweak levels are rows OE1 to OE11 in the stand-in
  registry (GEN.725).
generation-output: changed
