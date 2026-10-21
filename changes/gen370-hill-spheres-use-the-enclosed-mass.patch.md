**Headline:** Every galaxy's Hill spheres use the mass the galaxy holds inside the object's own distance from the centre, not the whole halo's mass.

### Changed
- **A star's, binary's, compact object's and rogue's Hill sphere takes the galaxy's mass inside its own radius** (`galaxy/enclosed_mass.py`, GEN.370). The central mass was the mass chain's whole virial mass at every radius (`mass_chain.central_mass_kg`), but most of the halo lies outside a star's orbit, so the Hill spheres were several times too small (a drawn galaxy of 9.5e12 solar masses of halo and 9.0e10 of stars encloses 1.6e10 inside 1 kpc and 1.2e11 inside 8 kpc). The enclosed mass is the density model's stars (the smooth components integrated on a cylindrical grid, the arms averaged round each circle), the central black hole and the chain's own NFW dark halo (its halo mass less its stars, its virial radius and concentration). The gas is left out until GEN.338's potential gives the enclosed mass.
  - `Star.calculate_system_perimeter`, `BinaryStarProxy._calculate_system_perimeter_static` and the keep-out query (`db/query.py` `keep_out_radius`) read it at the object's own distance; `mass_chain.central_mass_kg`, `current_central_mass_kg` and `MassChain.halo_mass_kg` are gone.
  - A unit reads the galaxy's stored shape (`seeded(..., shape=...)`, bound by `run_galaxy`, `nebula_field` and `superbubble_field`), or the shape its seed draws with every option at its default (`shape_draw.default_shape`; `plan` and `generate` take the same default arm and nuclear options from it).
  - A galaxy twice the size encloses less at one radius, and so has wider Hill spheres (`tests/test_enclosed_mass.py`).
  - **The placement fill is rechecked** (`docs/design/course-avoidance.md` section 4.1): the Hill spheres of a mean-mass system fill 0.3 to 0.9 of the volume over a disk's inner and middle radii (the old value never passed 0.03), and in a dense galaxy the placer now fits 60 to 70% of the requested systems of a 4 pc cell at 1 kpc and 35 to 45% at 3 kpc, where it fitted all of them; sparser galaxies are unchanged except 1 system in 16 at 3.7 kpc.
- The generator epoch is 42.
generation-output: changed

### Added
- `tests/test_enclosed_mass.py`: the mass grows with the radius and ends at the galaxy's whole mass, the black hole is most of the centre, the stars agree with a Monte Carlo of the density model within 6%, the dark halo follows NFW, a galaxy twice the size changes the Hill radius, and a star's and a pair's perimeter read the mass inside their own radius from the shape the unit is given.
