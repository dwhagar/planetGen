**Headline:** A gas giant's or brown dwarf's magnetic field now comes from its own metallic layer, its heat and its spin, so giants of the same mass differ in field with their interiors and a slowly turning giant has a multipolar one.

### Changed
- **A giant's field is its dynamo (GEN.809, GEN.801).** The moment of a gas giant or brown dwarf was read off three anchors (Neptune, Saturn and Jupiter) with a random factor. It is now the dynamo of the metallic layer the body's interior was drawn with (`giant_dynamo.giant_field`): the layer's radii and mass, the share of the heat the giant radiates that drives it (DY2), and the spin, whose Rossby number sets the dipole's share. A slow spin makes the field multipolar. An ice giant keeps the anchor law until its layers carry radii (the structure of GEN.733).
- generation-output: changed
