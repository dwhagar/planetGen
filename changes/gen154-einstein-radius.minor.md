**Headline:** Black hole, neutron star and rogue planet pages now show their Einstein radius against the galactic centre.

### Added
- **Einstein radius on compact-object pages.** A black hole, neutron star or rogue planet page has a new row, "Einstein Radius (seen from Sol, against the galactic centre)": `sqrt(4GM/c^2 * d * (1 - d/D))` for a lens on the line between Sol and the centre, worked out from the mass and the galactic radius (nothing is stored). A 10 solar-mass black hole 1 kiloparsec from Sol shows 8.4 AU. A lens farther from the centre than Sol is has no row.
- `planetgen/physics/lensing.py` holds the formula.

TODO items resolved: GEN.154 / Added: none
