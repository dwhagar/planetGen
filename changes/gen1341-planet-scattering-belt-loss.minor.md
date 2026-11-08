**Headline:** A belt past the outermost planet loses mass to planets that scattered each other, by the Safronov number of that planet (GEN.1341).

### Changed

- `belt_origin` gains the shaping cause `planet-scattering`: when a system has an instability level above 0 and two giants or more, the part of the disk that starts at or past the outermost giant and lies inside its reach (apoapsis plus 3 chaotic half-widths) loses the share `instability x (1 - exp(-Theta / 10))`, Theta being the Safronov number of that giant (`scatter.safronov_number`, `THETA_EJECT`). Nothing is tuned to the observed cold share. Registry rows BG20 (reach rule) and BG21 (radius of a giant with none stored); sources R5138 to R5141.
- `belt_origin.Body` carries a radius, `belt_birth.chain_of` passes the planet radius, and `belt_birth.build` passes the instability level to `belt_origin.build`. A belt with one giant or none, or outside the reach, is as before.
- Measured: 600 seeds, outer belts (10 AU and past) keep their median mass (1.75 Earth masses); 2 of 48 change (32.3 to 23.6 and 14.7 to 12.0 Earth masses). The share of systems with a cold belt moves from 62 to 61 percent over 240 systems (F5V, G2V, K2V, 80 each). See the report `gen1341-planet-scattering-belt-loss-2026-11-07.md`. Every system with a scattered outer belt changes (epoch).
