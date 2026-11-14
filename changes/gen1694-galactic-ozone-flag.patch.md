**Headline:** A world's ozone is no longer held down by the galaxy's rate of nearby supernovae and bursts: the flag is now the chance that an event fell inside the ozone's recovery time, drawn from the galaxy's own scales, and a system with no place in a galaxy has none.

### Changed

- **The galactic ozone-loss flag was always on (GEN.1694, bug).** `radiation.ozone_lost` flagged every world with ozone whose star had more than ten lethal galactic events per Gyr, and the rate was worked out with the Milky Way's constants (a Sun radius of 8,000 pc and a supernova scale length of 1,450 pc) on any galaxy and, for a system made with no galaxy, at the default distance, where it is 165 per Gyr. All 199 oxygen-rich Earth copies of a probe carried the flag (ozone about 29 DU against 300, ground ultraviolet up to 50 times Earth's) and the world's Radiation score could never be Blue.
- **The rate is a chance, not a state.** A galactic event strips the ozone for years to decades, so the world is flagged only if an event fell inside `radiation.OZONE_RECOVERY_YR` (100 years, a stand-in) before today: `galactic_event_in_recovery_time` draws against `1 - exp(-rate * time)` on the star's own stream (the same event for every world of the star, no run draw moves). Even at the nucleus's rates this flags about one world in 25,000. The flare branch of the flag is unchanged. `LETHAL_EVENT_LIMIT_PER_GYR` is gone.
- **The rates use the galaxy's own scales.** The radius the supernova and burst rates are quoted at is `scales.reference_radius_pc()` and the supernova scale length is 0.56 of the galaxy's thin disk scale length (`activity.SN_SCALE_LENGTH_DISK_FRACTION`, the Milky Way's 1.45 kpc over 2.6 kpc); the nucleus's 2.9 kpc stays, it is set by the fluence of an episode. `SUN_GALACTIC_RADIUS_PC` and `SN_SCALE_LENGTH_PC` are gone.
- **A star placed in no galaxy has no rate** (`lethal_event_rate_per_gyr` is empty, shown as none on its page) and so no galactic ozone loss; it was rated at the distance a system with no place is taken at.
- Tests: the rate profile in the galaxy's scales, the chance and its determinism, the star with no place, the flag from a chance of 1.

generation-output: changed
generation-affects: fill
