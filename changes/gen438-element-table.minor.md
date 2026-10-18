**Headline:** The generator now has one table of the elements and their isotopes, with their sources.

### Added
- **One table of the elements and their isotopes (GEN.438).** `physics/elements.py` holds the atomic weights of the 118 elements (CIAAW's Atomic Weights 2021), the masses and natural abundances of the isotopes (NIST), the half-lives, decay modes and branchings of the radioactive ones (NUBASE2020), the energy each releases per decay (Lugaro and others 2018) and the decay chains of uranium-238, uranium-235 and thorium-232 down to lead, with the energy of every step. `molar_mass("CaSO4.2H2O")` gives any compound's molar mass from its atoms, and the gas, materials and stellar tables will take theirs from it. The radiogenic heat, activity and mass fraction of an isotope are there for the heat and dose items to use; the heat of uranium and thorium reproduces Dye's published table to 0.3 percent. Every data file states its sources and licence and the loader refuses one that does not.

TODO items resolved: GEN.438. Added: none.
