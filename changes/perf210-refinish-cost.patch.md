**Headline:** Building a system is about a fifth to a third faster, with every planet, seed and number exactly as before.

### Changed
- **Cheaper atmosphere escape history and haze fall (PERF.210).** The thermal-escape history of an air no longer works out a gas's diffusion limit again at every step when the air has not changed, nor for a gas that is losing nothing, and the haze fall works out the air's viscosity once for each temperature of the profile instead of once for each level. Twelve Sun-like systems went from 36 s to 29 s and twelve of any star from 17 s to 12 s, and a checksum of every stored value of every body is the same before and after, so no seed gives a different world.
