**Headline:** How much is left in an asteroid belt is now worked out from the disc and the collisions that grind it down (not yet used by the generator).

### Changed

- New library `physics/belt_mass.py`: the mass of a planetesimal belt after collisional evolution (Wyatt et al. 2007), from the disc's solids over the belt, its width, the star's mass, the age and the size of its largest body, and the belt's density word from its mass per area. Nothing in generation reads it yet; the wiring step (the belt's stored mass and density word) follows with the epoch train.
