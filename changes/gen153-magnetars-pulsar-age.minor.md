**Headline:** Neutron stars can now be magnetars, and only young ones pulse, as in the real Galaxy.

### Added
- **Magnetars.** A neutron star younger than 10,000 years has a 10% chance of being born a magnetar: a field of 1e14 to 1e15 gauss, a spin of 2 to 12 seconds, listed as "Neutron Star (Magnetar)". Older ones have decayed to ordinary neutron stars.
- `PULSAR_RATE_SCALE` multiplies the pulsar share (1.0 is the physical rate, 0 turns pulsars off).

### Changed
- **Migration 0096** widens the CHECK on `neutron_stars.pulsar_type` to allow `magnetar` (no row is touched).
- **The pulsar fraction followed no physics (a flat chance for every age).** A neutron star now pulses with probability `exp(-age / 10 million years)`, plus a small share of recycled millisecond pulsars among the old ones, so almost all neutron stars are dark, as observed. The core of a supernova remnant is as old as the remnant, so the young remnants still hold pulsars and magnetars.
- Compact-remnant ages may start at 100 years (was 1 million) so young neutron stars exist.
- `NEUTRON_STAR_PULSAR_CHANCE` and `PULSAR_MILLISECOND_CHANCE` are gone, replaced by `PULSAR_ACTIVE_LIFETIME_YR`, `PULSAR_RECYCLED_CHANCE` and the `MAGNETAR_*` settings in `tuning.py`.

TODO items resolved: GEN.153 / Added: none
