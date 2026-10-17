**Headline:** A neutron star's radius now follows from its mass and the measured radius band, no longer a flat draw from 10 to 13 km.

### Changed
- **The neutron star radius is the band of Miller et al. 2021 (GEN.463).** It was drawn uniformly from 10 to 13 km, a floor the NICER and gravitational-wave radii do not support (they give 11.6 to 13.1 km). It is now the centre of the band at the star's mass, 12.45 km at 1.4 solar masses and 12.35 km at 2.08 (the line between the two, carried on to 1.1 and 2.2), plus a normal draw of the spread (0.65 and 0.75 km) clipped at two deviations, so the radius lies in 11.1 to 13.8 km at 1.4 solar masses and no star is under 10.8 km. The measured radii of PSR J0030+0451 (12.71 km) and PSR J0740+6620 (13.7 km) lie inside the band (`physics/neutron_star.py`).
- A neutron star that is a companion in a multiple system takes the centre radius of its mass in place of a fixed 12 km, and the plausibility check reads the band at the star's mass. `NEUTRON_STAR_RADIUS_RANGE_KM` is gone.

generation-output: changed

TODO items resolved: GEN.463. Added: none.
