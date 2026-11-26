**Headline:** The brown dwarfs, neutron stars, black holes and stripped helium stars that hang in a multi-star system's tree now store their bulk density and binding energy like the single ones do.

### Added

- `dark_members` gains `density_g_cm3` and `binding_energy_j` (-1 for the infinite binding energy of a black hole), worked out by `body_strength.dark_member_values` from the member's mass and the radius the tree already drew it with (`bulk_values`, the code a single brown dwarf, neutron star or black hole goes through, GEN.1162). A system made before has them empty until its sector is regenerated; a migration does not draw.

generation-output: changed
generation-affects: fill
