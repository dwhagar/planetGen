**Headline:** The generator can now say how many small captured bodies a planet holds, from the capture windows its history gave it, instead of every giant getting a few by chance.

### Added
- **`physics/irregular_population.py` (GEN.1077).** The captured small bodies of a planet from its windows: the instability (its own close encounters with another giant, scaled by its Hill volume), the gas era (gas giants only, bodies of 1 to 500 km), a binary exchange for a massive capture at an ice giant, and a temporary-capture flag. The count is ground down by the collisions of the planet's age, and the largest bodies are drawn as objects (every body above 1e16 kg, then the next heaviest until 90 percent of the mass, at most 100). A Jupiter-like planet gets about one body above 40 km and a few tens above 5 km, the four giants come out the same order, and an Earth gets none. Nothing reads it yet; the captured-moon wiring (GEN.1083) does.
- **Stand-in rows IP1 to IP14** in the gap-filling plan (section 6.118), the stand-in list, the burn-down register and the wiring register.

generation-output: unchanged
