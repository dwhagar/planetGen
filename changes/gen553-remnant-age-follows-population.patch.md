**Headline:** A neutron star's or stellar black hole's age now follows its population: the time since the star that made it collapsed.

### Changed
- `stellar_evolution.collapse_age_gy` draws a progenitor from the IMF over the kind's progenitor range, a birth age from the population's window and age law (as any star's), and returns the birth age less the time the progenitor took to collapse. Neutron stars and stellar black holes use it in place of an age uniform over 10 Gy; a young population's remnants are under 0.1 Gy old, a bulge's or halo's 8 to 13 Gy. The sector fill gives scattered rows, below-cut remnants, ejected remnants and rolled remnants their population as it gives a brown dwarf.
- `youngest_progenitor_sol` moved from `multistar` into `stellar_evolution` beside it, and the three progenitor ranges moved into `tuning`; no stored field changes.
- Intermediate-mass and supermassive black holes, and a remnant given an age (a supernova remnant's core, a magnetar row), keep their ages.
- `docs/design/remnant-age.md` records the rule.

generation-output: changed

TODO items resolved: none (GEN.553 stays open for the remnant row, [Fe/H], the composition hook and the free brown dwarf's table) / Added: none
