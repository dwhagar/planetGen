**Headline:** A hot body whose water would not saturate its air now holds it all as vapor, so a steam-start Venus twin no longer ends with a 586 K ocean.

### Fixed
- `hydrosphere`: a body at a temperature under water's critical point whose whole surface water makes less pressure than water's saturation pressure there has no liquid at the surface (`vapor`), as the air's own fixed point already said (`Air.ocean_present` False). Before, only the boiling point against the total pressure was read, so 1.4e20 kg of water at 586 K under 118 bar of air was called a `surface ocean` (GEN.1773, found by TEST.786).

generation-output: changed
