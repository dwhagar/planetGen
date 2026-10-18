**Headline:** Any star's composition can now be worked out from what the generator already stores for it.

### Added
- **A star's composition, derived and not stored (GEN.440).** `physics/stellar_composition.py` gives the composition of any star, brown dwarf or remnant from its kind, its stored [Fe/H], its age and mass, and a key. The star's own random numbers come from a stream of their own keyed on the star, always drawn in the same order, so the same star has the same composition however many workers ran the galaxy, editing a star's metallicity or age works it out again with nothing to repair, and the members of a multi-star system or a cluster share their sequence and element pattern while each keeps its own phase. A red giant has its first dredge-up, a brown dwarf its deuterium, lithium and hydrogen limits at its mass and metal content, a white dwarf its core, atmosphere and radius, a neutron star its outer crust, a black hole its horizon quantities, and every kind lists what the design leaves open. No column is added and no stored value changes.

TODO items resolved: GEN.440. Added: none.
