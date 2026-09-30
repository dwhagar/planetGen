### Added

- **Schema v38: letter classes for nebulae, supernova remnants and
  asteroid fields.** Nebulae are classed A-Q and supernova remnants R-W
  (`program_constants.NEBULA_CLASSES`), each with what fills it: dominant
  species, particle density, temperature and optical extinction. A
  remnant's class follows its progenitor and core (a Type Ia remnant is
  always W; a pulsar wind nebula needs a pulsar). Asteroid fields get a
  class made of a composition-and-density letter and a size digit, like
  `C3`, with composition drawn from real asteroid families
  (`ASTEROID_FIELD_COMPOSITIONS`). Nebulae gain a `diffuse` family.
  Existing rows get the class they most likely are.
