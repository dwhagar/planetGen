**Headline:** A number with a unit is now one `Quantity` object that keeps what it was given and answers in every other unit of its kind.

### Added
- **`planetgen.physics.quantity.Quantity`** (GEN.569). `Quantity("power", 3.828e26, "W")` keeps the number and unit as
  given and never converts them; `.kW`, `.Lsun` and `.as_("erg_s")` give every other unit of its kind, worked out from
  the truth on first read and kept, and `set` drops them so every form is in step. A relation of the registry answers
  a quantity of another kind (`related`: a mass's gravitational parameter, a radius's diameter, a period's frequency).
  A value is checked finite and inside its kind's bounds when made or set, and a unit of another kind is refused by
  name. `to_dict` gives `{kind, value, unit}`. It is for the edges; the generator's kernels stay on plain floats.

### Changed
- `unit_registry.convert` multiplies or divides by the power of ten itself when two units' factors are a power of ten
  apart (`unit_registry.scale`): 1.5 km is exactly 1500.0 m, where the two-step product rounded twice for about a
  third of all values.
