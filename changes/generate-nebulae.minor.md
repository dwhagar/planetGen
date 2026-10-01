### Added

- **Nebulae and supernova remnants are generated with the stars they
  need.** Sectors now make molecular clouds (dark classes M-Q) at the
  research rate. Every planetary nebula comes with its own new hot white
  dwarf system at its center. Every O star, and half the B0-B2 stars,
  sits in an H II region (classes C-E or G), and a few later B and A
  stars light a reflection nebula (`NEBULA_HOST_RULES`). A core-collapse
  remnant's neutron star or black hole has drifted off-center by its
  birth kick (`SUPERNOVA_KICK_SPEED_RANGE_KMS`) times the remnant's age.
  Diffuse gas (classes A-B) is the background and isn't generated.
- The sector's nearby phenomena now carry each nebula, supernova remnant and asteroid field's letter `class`.
