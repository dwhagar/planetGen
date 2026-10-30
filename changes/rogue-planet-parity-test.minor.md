**Headline:** A test now holds a rogue planet to being a bound planet without a star, and says how far from that it is today.

### Added
- **The rogue planet parity contract (TEST.498).** `tests/test_rogue_planet_parity.py`: the columns a bound planet stores that a rogue planet does not (148 of the 228 under the same name, apart from the host-only ones: the orbit, position and velocity round a host, the zone, the Hill radius and the reflex offset), and the 13 names a rogue uses for what a planet already stores (`mass_kg` for `mass`, `radius_km`, `surface_temperature_k`, `planet_type` for `body_type` and the like), as a ceiling that can only fall. Five strict expected failures state what GEN.1111 and the rogue train build: every planet column on a rogue, no column of its own beyond the rogue note, galactic motion and ejecting system, a rogue that is a `Planet` (no rogue-only generator), a moon system of `Planet` moons instead of a flag, and the class set by the one classifier. Each mark comes off when the rogue is built through the planet code.

No generator output changes.

generation-output: unchanged

TODO items resolved: none (TEST.498 stays open until GEN.1111 turns the five expected failures into passes).
