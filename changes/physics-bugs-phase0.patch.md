### Fixed

- Planets no longer pile up in the cold zone (GEN.37). The first planet now
  sits at 5-40% of the habitable zone's inner edge and each next one 1.3-1.8
  times farther out, stopping at the protoplanetary disk's outer edge, so a
  red dwarf's planets stay close in. About 30% of planets are now hot, 6% in
  the ecosphere and 64% cold (it was about 1%, 2% and 96%).
- Gas and ice giants have real densities, and super-Jupiters occur (GEN.34).
  A giant's mass is drawn first (dN/dlogM ~ M^-0.31 within its class) and its
  radius follows the Chen & Kipping (2017) mass-radius relation, so the median
  giant is about 1.2 g/cm^3 (it was 0.25) and about 8% of giants pass a
  Jupiter mass. A giant is at most 1% of its star's mass, so red dwarfs don't
  get super-Jupiters.
- Rocky planets get the moon classes their size and zone allow, not only
  Class D (GEN.35): a class only has to fit at its smallest size, and a moon's
  radius is capped so it stays under the planet's moon size and a tenth of its
  mass.
- A moon regenerated after its planet moves only gets a class and size a moon
  of that planet may have, never a gas giant or a blacklisted class (GEN.36),
  and is never too large for its planet (GEN.25). Moons are re-spaced after one
  is regenerated, and a moon heavier than a tenth of its planet is dropped
  when the planet is reclassified.
- A planet reclassified after a spacing push keeps the Hill sphere of its real
  orbit, not of the spot its new class first drew.
- Rogue planets follow one mass function, dN/dlogM ~ M^-0.65 from 0.1 Earth
  masses to 13 Jupiter masses (GEN.45): about 96% terrestrial and 4% gas
  giants (it was 91% and 9%), still 6.5 per star.
