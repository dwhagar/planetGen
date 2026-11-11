**Headline:** A planet whose orbit the secular exchange tilts now keeps its spin axis with the orbit, so its stored axial tilt matches its axis again.

### Changed

- The secular exchange of GEN.1450 changes a planet's orbital inclination and node, but left its spin axis and axial tilt as they were drawn, so `axial_tilt_deg` no longer matched the angle between the axis and the orbit normal (`test_every_generated_system_body_spins` red on main). The axis is now turned with the orbit normal (`planets.carry_spin_with_orbit`, `spin.carry_axis`): the obliquity and the angle round the normal are kept, because a planet's spin precesses about its normal far faster than the plane's nodes regress. A locked or resonant planet keeps its Cassini-state share of the new inclination (its tilt moves in proportion to the folded inclination). The sky of the other suns follows the new normal, the day and night shares follow a changed tilt, and the moons are set in the new equator again. No draws.
- Registry row SC6 (section 6.244).

generation-output: changed
generation-affects: fill
