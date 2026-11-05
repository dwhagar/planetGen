**Headline:** A supernova remnant's ejecta and swept-up gas, and its progenitor's mass, can now be worked out from the physics (not yet stored or shown).

### Added

- `remnant_composition.py` gives a remnant two parts: the interstellar gas swept up inside its radius, and the ejecta of the star that exploded, element by element, with the progenitor's initial mass. A Type Ia remnant throws out the 1.37 solar masses of the W7 model (mostly iron, no hydrogen or helium); a core-collapse remnant throws out oxygen-rich ejecta that grow with the progenitor's mass, from a progenitor drawn from the windows of the core it left (a neutron star or a black hole).
- A progenitor above 30 solar masses is stripped of its hydrogen, and a third of the lighter ones are; the remnant's core and the ejecta add up to the star's mass at death.
- Checked against Cassiopeia A: the swept-up mass at 350 years and 3 particles per cm^3 is 8.0 solar masses (7.9 measured) and the ejecta of a stripped 16 solar mass star is 2.4 (2.2 measured).
- The tables are the published ones (a Type Ia and a Type II table, read as `remnant_yields.json`); what scales them with the progenitor's mass, strips the star and bounds the sums are stand-ins with rows SNC1 to SNC7 in the registry.
- Nothing reads the model yet: GEN.1276 writes the values at generation and GEN.1272 shows them. The swept-up gas is the star's own birth gas until the shared nebula model (GEN.1269) gives the gas of the place.

generation-output: unchanged
