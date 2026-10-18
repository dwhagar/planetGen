**Headline:** The mass of a white dwarf, neutron star or stellar black hole comes from one place in the code instead of six; no seed's output changes.

### Changed
- New `physics/remnant_mass.py` holds the initial-final relation for white dwarfs, the uniform neutron star and black hole mass draw, and the moments the placement terms use. The star model, the multi-star member draw, the scatter's planetary-nebula central star, the compact objects, the phenomenon scatter's mass law and the placement terms call it; the six copies are gone with no shim.
- Every draw uses the same stream the same way, so no seed gives a different object (`tests/test_remnant_mass.py` compares each rule with the line it replaced); the numbers stay in `tuning`.

TODO items resolved: GEN.562 / Added: none
