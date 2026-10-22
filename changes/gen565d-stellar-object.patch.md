**Headline:** Stars, binary proxies and compact remnants now share one `StellarObject` family class, and the `_skip_property_init` workaround is gone.

### Changed
- **`generation/stellar_object.py`** (GEN.565d): `StellarObject` holds what every star-like object shares: the persisted field
  list and class-level defaults, the save and load of `habitable_zone` and `lifespan`, the age rules for planets
  (`planet_age_requirement_gy`, `age_ceiling_gy`, `adjust_age_for_planets`), `calculate_system_perimeter` and
  `calculate_galactic_orbit`, and `get_star_spectral_class` and `get_star_evolutionary_profile` (moved from `star.py`,
  their importers updated, no copy left). `Star`, `CompactRemnant` and `BinaryStarProxy` subclass it, so `Star.__init__`
  no longer takes `_skip_property_init` and `Star.from_params` builds on the base class directly. A black hole or neutron
  star is no longer an instance of `Star`. Every stored value, field name and draw is unchanged.
