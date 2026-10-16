**Headline:** Stars, black holes, neutron stars, rogue planets, comets, asteroid fields, nebulae and remnants keep a movement vector and a mass.

### Changed
- **One movement vector per free object (GEN.336).** A black hole, neutron star, rogue planet, interstellar comet, asteroid field, nebula or supernova remnant is stored with a velocity vector (the tangent of the galaxy's rotation at its place, at its speed) instead of a scalar speed, and the nightly orbit update turns the vector with the position. A star and a close pair take their system's vector. Migration 0103 drops the scalar speed columns.
- **Masses are derived.** A comet nucleus is its sphere at the density of comet 67P, and a nebula or remnant is its bounding sphere at its particle density with the helium counted; an asteroid field, whose mass is not known, is a stated zero.
- **Pages.** The galactic speed shown for an object is the length of its vector.

TODO items resolved: GEN.336. Added: none.
