**Headline:** A helium-core white dwarf lighter than the lightest cooling sequence (0.2 solar masses) now has a radius, where it had none.

### Changed
- `stellar_composition._grid_radius_solar` reads the lightest sequence for a white dwarf under 0.2 solar masses and scales it by the zero-temperature mass-radius relation's ratio of the two masses, the way a star above 1.30 solar masses already reads the heaviest. The composition sweep (`test_star_composition_sweep`) found such a star with no radius in the train 11 draws.
