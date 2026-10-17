**Headline:** Every galaxy now draws its own stellar mass, halo mass, virial radius and concentration from its seed, and the Hill spheres read that mass in place of the fixed Milky Way mass.

### Changed
- **The galaxy's mass chain (GEN.400).** `galaxy/mass_chain.py` draws the stellar mass from the double Schechter function (Baldry et al. 2012) inside 3e9 to 2e11 solar masses, the halo mass from the stellar-to-halo relation (Moster et al. 2013) by the approximation Boss chose, and the virial radius and NFW concentration (Dutton and Maccio 2014, 0.11 dex scatter), all from the galaxy seed: the same seed gives the same galaxy.
- **No fixed Milky Way mass.** `MILKY_WAY_MASS` is removed. A star's, a pair's and a compact object's Hill sphere (the system perimeter and the keep-out radius) read the drawn galaxy's virial mass, so the perimeters differ with the galaxy. `keepout.galactic_hill_radius_km` and `compact_keep_out` take the central mass as an argument.
- **Stays.** `HALO_MASS_SOLAR` and `HALO_SYSTEM_FRACTION` are the stellar halo (a fraction of the visible mass), not the dark matter halo; their law is GEN.406's research.

generation-output: changed

TODO items resolved: GEN.400. Added: none.
