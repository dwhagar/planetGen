**Headline:** More than one kind of life can now live on the same world, each on its own schedule.

### Added
- **Several life classes on one planet or moon (GEN.1209).** Any class the conditions support (one for each class number, any solvent) can take hold beside the dominant one: it draws its own schedule on a stream of its own, lives only if its first step is inside the star's age and the air that all the classes make does not poison it, and has its own stage, animal grade, complexity and four kind scores. The body's microbial, flora, fauna and sentient scores are the highest of each over its classes. `planet_life` and `moon_life` (schema v197, migration 0197) keep one row for each class, the dominant one too, and `life_living` carries the list on the body and in the system file.

### Changed
- **The wait for complex cells follows the oxygen of all the classes.** A class that breathes oxygen starts it when the oxygen of the whole mix has risen (GEN.1213), not at its own photosynthesis plus a delay.

generation-output: changed
