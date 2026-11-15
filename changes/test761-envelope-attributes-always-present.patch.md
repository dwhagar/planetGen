**Headline:** a planet that kept no envelope now carries the same private attributes as one that did, so a rogue and a bound planet of one mass have the same attributes again.

### Fixed

- TEST.761: `apply_envelope` set `_envelope_helium_y` and `_tidal_radius_km` only on a planet that kept an envelope, so a rogue planet (no star's X-rays to strip it) that kept a small one had two attributes a bound planet of the same mass, which lost it, never got (`test_a_rogue_and_a_bound_planet_carry_the_same_attributes_and_their_moons_do_too`, 4 cases, red on main). Both are now set to `None` with the other start-of-draw resets in `generate_planet_properties`, so every planet has them and a planet drawn again in place carries nothing of its last draw's helium share or tide radius. Every reader already treated a missing attribute as `None`, and the pinned golden galaxy is unchanged.
generation-output: unchanged
