**Headline:** The Compact binaries panel of an unbuilt system now gives the X-ray luminosity of each accreting pair.

### Changed
- **Stored member radius and X-ray luminosity (GEN.952).** Every member of a stored companion tree keeps its radius (`scatter_members.radius_au`, schema v184). An accreting neutron star or black hole with a companion now shows the power it radiates in X-rays in the Compact binaries panel of the system page and as `x_ray_luminosity_w` in the API, worked out from the masses, radii and orbit of the pair with the same law the generator uses. Systems made before this change have no radius and show a dash until the galaxy is regenerated.

generation-output: changed
