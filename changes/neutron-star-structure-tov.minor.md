**Headline:** A neutron star's radius and layers can now be computed from the 21 equation-of-state tables by integrating the general-relativity structure equations.

### Added
- **Neutron star structure from an equation of state (GEN.968).** `physics/neutron_star_structure.py` integrates the Tolman-Oppenheimer-Volkoff equations over the 21 Malik, Pais and Providencia 2024 tables: a star's mass, radius, compactness, crust and neutron-drip radii, each table's mass-radius curve and maximum mass (2.0 to 2.75 solar masses), and the radius band the tables span at a mass (R of 1.4 solar masses: 11.8 to 13.8 km against NICER's 12.45 +- 0.65). Nothing in generation calls it yet; the wiring is GEN.969. Three stand-ins registered (NS1 to NS3).
