**Headline:** A giant's layered interior is now saved with the world and read back with it.

### Changed
- **Giant interior storage** (DB.89): each giant planet and moon now keeps the kind of interior it was drawn with, its compact and dilute core masses, the pressure where its hydrogen turns metallic, its helium depletion and whether it has a helium shell, in new columns on the planet and moon rows. Its layers (kind, state, pressure and radius at the top and bottom, and the mass of each material) are saved one row per layer in two new tables, `planet_interior_layers` and `moon_interior_layers`. Worlds saved before this version have no interior until the galaxy is regenerated at the generator epoch. Database version 154.
