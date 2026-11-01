**Headline:** The equilibrium temperature can now take a body's own heat, which is what sets the temperature of a planet with no star: the next piece of building a rogue through the planet code.

### Changed
- **The body's own heat in the equilibrium temperature (GEN.1111, TEST.498).** `climate_draws.equilibrium_temperature_k` takes an optional `internal_flux_w_m2` and adds it to the absorbed starlight as a flux (`T^4 = (1 - A) F / 4 sigma + F_int / sigma`), the way the giants' effective temperature already adds theirs. With none it is the old value exactly, so no body changes; under a star a rocky planet's heat moves it by less than a part in a thousand. A body with no light is then its heat's temperature, as the rogue chain has it, instead of the 2.725 K floor the bound chain returned for a body far from a star. Nothing passes it yet; the starless planet's temperature reads it in the lock wiring.

generation-output: unchanged

TODO items resolved: none (GEN.1111 stays open).
