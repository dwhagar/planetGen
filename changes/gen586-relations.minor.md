**Headline:** The physical laws the generator writes by hand about fifty times are now declared once in `physics/relations.py`, each with a fast SI kernel and a Quantity form, and a dimension test in the math check.

### Added
- **`physics/relations.py`** (GEN.586): 21 named laws (Schwarzschild radius, surface gravity, escape, circular and
  vis-viva speeds, Kepler's third law and mean motion, break-up rate, rotation period and rate, Hill radius, Roche
  limit, sphere volume, mean density and radius, black-body luminosity and effective temperature, inverse-square flux,
  equilibrium temperature, mass flow). Each is one function: `relations.surface_gravity(mass_kg, radius_m)` takes
  plain SI floats and is the kernel; `relations.quantities.surface_gravity(mass_quantity, radius_quantity)` runs the
  same formula on `Quantity` values and answers a `Quantity` of the declared kind, bit for bit the kernel's answer.
  Constants come from the registry's table (`c.G`), never typed. The vector laws (obliquity, galactic orbit) are left
  to GEN.584.
- **A dimension test**: each formula is run once on probe objects that carry a dimension, and the answer must have its
  declared kind's dimension; the new `relations_have_their_dimension` invariant in the math check runs it before
  anything else, so a missing square, a missing constant, a wrong root or a length added to a mass fails the suite.
- Kinds area and volume (m², km², cm², m³, km³, cm³, liters) in the unit registry.
- Nothing is replaced yet: GEN.587 moves the hand-written sites onto the laws, file group by file group, wherever the
  golden set stays bit-identical.
