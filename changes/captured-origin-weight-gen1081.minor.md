**Headline:** A moon's chance of having been captured now falls with its mass as a tested calculation, so a big moon is rarely a captured one; the moon generator does not use it yet.

### Added
- **Captured-origin weight by moon mass (GEN.1081).** Two pure functions in `physics/moon_composition` (`captured_weight`, `origin_probabilities`) weigh the `captured` origin of a moon against the others by its mass: full weight at 1e16 kg or less, then falling as the mass to the power of -1/3 from the size law of the irregular moons that survive (0.046 at 1e20 kg, 0.008 for a Triton). The registry has a new row (CM17). Nothing in the generator calls them yet: the irregular moons item (GEN.1083) wires them on a later generation train.

generation-output: unchanged
