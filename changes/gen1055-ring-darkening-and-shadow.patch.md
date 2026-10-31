**Headline:** A pure module now gives how bright a ring is for its age and dust, and how much sunlight a ring shades from a latitude; nothing uses it yet.

### Added
- **Ring brightness and ring shadow (GEN.1055, rows R6 and R8).** `physics/ring_darkening.py` has `pollution_flux`, `surface_density_kg_m2`, `pollution_gain`, `albedo`, `ring_brightness` and `shadow_transmission`, with four registry stand-ins (RGD1 to RGD4, section 6.148 of the gap-filling plan) and a burn-down line for the one hand-tuned row. It is not wired into `rings.apply` or the climate yet.

No generator output changes.

generation-output: unchanged

TODO items resolved: GEN.1055 pure part (the wiring into rings and the climate stays open).
