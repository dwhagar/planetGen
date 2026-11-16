**Headline:** A star's own disk dispersal time is now a law with sources, an exponential of 2.5 Myr at 0.5 solar masses shortened by a close companion, ready to replace the fixed 3 Myr step.

### Added

- GEN.1551 slice 0 (unwired, no stored value moves): `physics/disk_dispersal.py` holds the law the star's disk time will be drawn from. The share of young stars with a primordial disk falls about exponentially with age, with a characteristic time of 2.5 Myr (half-life 1.7 Myr) [R5298], 2 to 3 Myr in 22 clusters of 1 to 100 Myr [R5296] and a half-life of 1.3 to 2 Myr in 69 clusters [R5297]; the characteristic time falls from 3 Myr for brown dwarfs to 1 Myr over 1.3 solar masses. A star's time is `-tau ln u` from a draw of its own stream (`draw_dispersal_time_yr`, keyed to the star and not to its age), and a companion inside 47 AU shortens it with the same draw (an exponential of 0.91 Myr, two thirds gone within 1 Myr [R5225]), so GEN.1475 and this item are one draw. Stand-ins DZ1 to DZ3 say so. The shares a lane remembered (0.8 at 1 Myr, half at 3 Myr, 0.1 at 5 to 6) are not an exponential; the law gives 0.67, 0.30 and 0.11 at 0.5 solar masses, and the choice is left to Boss. Slice 1 wires it.

generation-output: unchanged
