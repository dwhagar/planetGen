**Headline:** Every star now has a stored metallicity, and the sector map shows it.

### Changed
- Migration 0095 adds `metallicity_fe_h` (DOUBLE, nullable) to `stars` and `bright_stars`: the star's [Fe/H] in dex against the Sun. A scattered star draws its own at scatter time from its population and its distance from the galactic centre, on a stream of its own seed, so no other scatter draw moves; the fill that builds its system keeps it, and the two stars of a binary share one. Rows made before 0095 hold NULL, and so does a cluster star until the cluster table (GEN.159) gives it its cluster's value; a star built from a NULL draws its own as it did before.
- The info panel of a star on the Sector Map and in an opened sector now has a "Metallicity [Fe/H]" row (for example "−0.12 dex") when one is stored.
- `store.BRIGHT_STAR_COLUMNS` ends with `metallicity_fe_h`, so a scatter row has one more value.

TODO items resolved: GEN.158 (the stored column and the display) / Added: none
