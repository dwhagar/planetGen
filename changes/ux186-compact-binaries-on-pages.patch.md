**Headline:** A scattered star or phenomenon with a neutron star or black hole in a close pair now says what kind of pair it is.

### Changed
- **Compact binaries on the pages and in the API (UX.186).** The page of an uncharted star or phenomenon that holds an X-ray binary or a spider pulsar gains a Compact binaries panel: the kind (high-mass X-ray binary, Be/X-ray binary, low-mass X-ray binary, black widow or redback), the compact object, its companion, and whether the accretor is a microquasar. The same list is in the API as `compact_pairs`, and the stored `flags` of members and the `interaction` and `flags` of orbits are now returned too. The X-ray luminosity is not shown yet: it needs each member's radius, which a stored tree does not keep (GEN.952).
