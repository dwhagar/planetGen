**Headline:** The bare stage of a body now has its own scale height and surface temperature, so the forward build can place moons and measure tide decay before any air exists.

### Changed

- GEN.1548 (built ahead for the forward build): new `physics/bare_body.py` with `scale_height_km` (a bare body's atmospheric scale height from its light, gravity and the mean molar mass of the gas it could hold; stand-in BB1, registry 6.273), `surface_k` (the body's own surface temperature, else the temperature its light gives it with no air, in place of the fixed 273 K of `interior_heat.surface_of`), `equilibrium_k`, `reference_albedo` and `thermal_tide_is_known`. Nothing calls them yet; the reordered build switches `moon_orbit_bounds_km` and the tide decay to them.

generation-output: unchanged
generation-affects: none
