**Headline:** The overlap checker now also holds every irregular moon to the limits a planet's Hill sphere allows, and the clearance ledger gives the irregular zone as a span.

### Changed

- GEN.1086 (slice 3): `physics/clearance.audit` checks the stored irregular moons (`irregular-hill`: past the stable share of the planet's Hill radius, 0.4 prograde and 0.55 retrograde; `irregular-kozai`: inclined between 55 and 130 degrees; `irregular-inner`: closest approach inside the outermost regular moon's orbit), all hard, from what the system stores (`stored_rows` now reads `irregular_moons`). `clearance_ledger.irregular_bounds_km` gives the span the irregular moons hold, from the outermost regular moon out to the stable share of the Hill radius.
- A neighbor already keeps out of that space: another planet stays ten mutual Hill radii away and a belt five Hill radii, so no placement changes; the generated systems of the clearance sweep hold none of the new violations.

generation-output: unchanged
