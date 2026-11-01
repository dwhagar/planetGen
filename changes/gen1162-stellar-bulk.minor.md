**Headline:** Stars, white dwarfs, neutron stars, black holes, rogue planets and interstellar comets now keep their mean density and how tightly they are bound.

### Added
- `stars`, `black_holes` and `neutron_stars` keep `density_g_cm3` and `binding_energy_j` (GEN.1162, schema v200, migration 0200). The energy is each kind's own formula (`physics/binding_energy.py`): the n = 3 polytrope for a star over a third of a solar mass, n = 3/2 for a convective star, brown dwarf or white dwarf, [R2573] for a neutron star. A black hole stores -1, meaning infinite, and the property sheet says so in words; its density is the mean inside the event horizon.
- `rogue_planets` keep the density and binding energy; a terrestrial one also keeps a roundness ratio and a strength from its rock and the water of its shell and ocean (a giant or a brown dwarf is a fluid and has none). `interstellar_comets` keep what `comets` do: density, ice share, roundness ratio, strength and binding energy.
- Not stored yet: the dark members of an orbit tree (a brown dwarf, neutron star or black hole in `dark_members`) have no radius of their own to work from.

generation-output: changed
