**Headline:** A belt body that became a planet now has the radius the planet structure model gives its makeup, so no body of the sweep breaks the radius rule.

### Changed
- **Promoted belt bodies under 1e21 kg take the structure model's radius (GEN.1692).** `bulk_draw._take_belt_member` kept the member's own radius whatever its porosity, so a member denser than its grains' zero-pressure mix (porosity none) had a radius 0.1 to 2 percent under what the structure model gives its composition and mass; four such bodies failed the radius rule of the sweep (`tests/audit/sweep_bulk.py`, a millionth for a planet). The porosity is still what the member's radius leaves void, and the radius is now the structure model's for it, as every planet's. On 20 sample systems, 1 body of the 35 promoted under 1e21 kg moves (radius up to 2.0 percent, density down to 5.8 percent); the rest and every heavier body are as they were.

generation-output: changed
generation-affects: fill
