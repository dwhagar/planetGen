**Headline:** A system's pair is read from its orbit node, not kept in a second copy, so a pair has one place to be wrong.

### Changed
- **`StarSystem.pair_orbit` replaces `WideBinaryPair` and `system.wide_binary` (GEN.640, slice b).** The pair's axis, eccentricity, period, periapsis, apoapsis, angles, speed and the places of its two stars are now worked out from the pair's node in the orbit tree (`generation/pair_orbit.PairOrbit`). On 300 seeded wide pairs of two to five stars every number equals what the old class held, to 1e-9. The wide-pair draws stay in `generation/wide_binary.py`.
- **Serialization schema 10.** A saved system no longer carries a `wide_binary` key; it is read from `orbit_nodes`.
