**Headline:** The galaxy's mass chain is merged and retired: every galaxy draws its stellar mass, halo mass, virial radius and concentration from its seed.

### Changed

- GEN.400 (#2011, generator epoch 10) is retired: the stellar mass comes from the Schechter function, the halo mass, virial radius and concentration follow, and the Hill spheres read the virial mass; the fixed Milky Way mass is gone.
- GEN.401, GEN.403 and GEN.399 no longer wait on it.
- GEN.406 now carries the stellar halo constants (`HALO_MASS_SOLAR`, `HALO_SYSTEM_FRACTION`): they are the stellar halo, not the dark matter halo, so GEN.400 left them, and the laws table has no stellar halo law until GEN.406 researches one.

TODO items resolved: GEN.400 / Added: none
