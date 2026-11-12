**Headline:** The two regions of a trojan swarm, where each lies and how it wanders about its point, can now be drawn for any pair; nothing reads it yet.

### Changed

- GEN.1150 (slice 0): new `physics/trojan_swarm.py`, a pure draw that follows `trojans.swarm_split`. `draw_regions(host_mass, secondary_mass, trojan_mass, rng)` gives the leading (L4) and trailing (L5) `Region`: its share of the mass, the libration amplitude of its bodies (a truncated Rayleigh below the 23.9 degree tadpole limit), the libration phase, and its mean eccentricity and inclination (truncated Rayleigh below 0.2 and 35 degrees). A pair that cannot hold the swarm gives none, and the number of draws is fixed so no other draw moves. `longitude_offset_deg` and `longitude_deg` give the libration overlay the orbit updater adds to the secondary's mean longitude, `+-60 + A sin(2 pi t / P_lib + phase)` with the slow libration period of `trojans.libration_period`. The three scales are stand-ins TJ9 to TJ11 (tier T4, tweak L4) in section 6.121 of the gap-filling plan.
- Nothing in the generator reads it yet (wiring register): the stored swarm and its migration wait for GEN.1006, and GEN.1151 draws the whole bodies from these regions.

generation-output: unchanged
