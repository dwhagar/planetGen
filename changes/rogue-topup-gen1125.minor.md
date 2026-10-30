**Headline:** The plan for topping a sector up with rogue planets now exists as a tested calculation: how many to add, their masses, places, speeds and an "apparent" origin; no world uses it yet.

### Added
- **Rogue top-up planner (GEN.1125).** A new module (`physics/rogue_topup`) takes a sector's mass shortfall, the Rogues level (None, Few, Normal or Many, with ceilings of 0, 2.5, 4.3 and 10 per star, kept in `tuning.py`) and the rogues already there, and says how many rogue planets to add (never past the ceiling), their masses by the existing bins (about 96 percent terrestrial, 4 percent gas giant), a uniform place inside the sector, a velocity through the Galaxy and an origin marked "apparent" that never links to a real system. The placed mass plus the mass left over equals the shortfall exactly, and the leftover is returned for the light asteroid fields (GEN.1126). Five stand-in rows (TU1 to TU5) are registered. Nothing reads it yet: the region run wires it (GEN.1120).

generation-output: unchanged
