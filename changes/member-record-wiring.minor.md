**Headline:** The companion stars of a system are now the same stars its orbits were drawn for, so a star can no longer end up wider than its own orbit after the build.

### Changed
- **One evolution per star (GEN.1009).** The tree that draws a system's orbits and the build that makes its stars used to evolve each companion separately, with different random draws and a different metallicity, so the built star could differ by several per cent in radius from the one the orbit was drawn around. Each companion is now evolved once from its own seed, and the build makes the star from that record. The primary's real radius is what the tree clears.
- **Stripped helium companions are left out of a built system (GEN.1009).** A stripped helium star was built as an ordinary star of its birth mass, a different star from the one drawn; it now stays in the stored tree like a neutron star or brown dwarf and is not built.
- **Largest radius and mass of a stored member (GEN.1009).** `scatter_members` gains `radius_max_au` and `mass_max_solar` (schema v187), the red-giant tip and the birth mass, for the history rule to come.
- generation-output: changed
