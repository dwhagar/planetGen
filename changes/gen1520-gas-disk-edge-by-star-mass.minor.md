**Headline:** A star's protoplanetary disk now ends at a radius that follows the star's mass (100 au for the Sun, mass to the half power, a factor of 2 of scatter per star) instead of 18 snow lines, so Sun-like stars grow cold giants and belts farther out and small stars get larger disks than before.

### Changed

- The outer edge of the gas disk a star's planets, belts and comets formed in is `physics/disk_size.edge_au`: median 100 au x (M/Msun)^0.5 (45 au at 0.2 Msun, 71 au at 0.5, 118 au at 1.4), drawn once per star with a 0.3 dex log-normal scatter cut at 2 sigma from a stream keyed by the star (no stored column, no other draw moves). It replaces `tuning.DISK_OUTER_RADIUS_SNOWLINE_MULTIPLIER` (18 snow lines, 49 au for the Sun, 5.8 au for an M2V dwarf), which followed the star's light rather than its mass. Evidence and the decision are in `reports/gen1520-disk-size-by-mass-2026-11-11.md` (Boss approved 3:04 PM EST 11-11); stand-ins DS1 to DS4 in `gap-filling-plan.md` section 6.274, design note `docs/design/generation/gas-disk-edge.md`.
- It sets where random planet slots end and the oligarch chain's limit (`system.py`), the outer-disk belt and the scatter shares (`belt_birth.py`), the zone a comet's ice formed in (`comet.py`) and where the pebbles end and which volatile ices evaporate on the way in (`disk_drift.py` through `mixture_subject` and `mixture_giant`, now taking the edge as an argument; a larger edge is a colder one, so more of the CO2, CO and N2 lines lie inside the disk).
- New redraw kind `disk_edge.scatter`. Planet, giant, belt and comet counts, the outer planet and belt distances, comet reservoirs and giant C/O change for every star; the before and after sweep (150 systems, same seeds) is in the report.

generation-output: changed
generation-affects: plan, scatter, fill, population
