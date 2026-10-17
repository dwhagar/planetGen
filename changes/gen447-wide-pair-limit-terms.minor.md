**Headline:** The widest orbit a wide pair keeps now uses the relative velocity of the passing stars and the mass-weighted mass of everything that passes, remnants included.

### Changed
- The encounter limit of a wide pair (Jiang and Tremaine 2010, equation 1) takes the dispersion of the relative velocity, root 2 times a population's own, where it used the population's own, and the second moment of the perturbers' mass over its first, `<M^2> / <M>`, where it used the mean system mass (0.4 solar masses). The mass function is the generator's own initial mass function aged to each population's age range, with the white dwarfs, neutron stars and black holes of what has left the main sequence. For the stars alone it gives 0.59 solar masses against Jiang and Tremaine's 0.64; with the remnants 1.46 for the whole disk, and the black holes alone hold about twice the second moment of all the stars. The widest axis a wide pair keeps is about 2.5 times narrower than before in every place (11,400 AU at the Sun's distance for two solar masses at 10 Gyr, against 29,000).
- `place_terms.mass_weighted_perturber_mass_solar(populations, remnants=True)` (the stars alone with `remnants=False`) and `place_terms.RELATIVE_VELOCITY_FACTOR`; `PlaceTerms.perturber_mass_solar` defaults to the whole disk's value; `from_profile` mixes the populations a cell holds.
- `docs/design/placement-by-population.md` 3.3 records the derivation, the paper's definitions and the new table.

generation-output: changed

TODO items resolved: GEN.447 / Added: whether the dark remnants should count as perturbers (the black holes do most of the second moment), for Boss and TEST.196
