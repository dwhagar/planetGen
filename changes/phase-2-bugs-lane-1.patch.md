### Fixed
- A moon's Hill sphere (and so its minimum orbit spacing and the facility orbit slider) is measured about its planet, not its star (GEN.138); the Hill sphere uses the pair's total mass.
- Class N (a Venus analog) no longer carries a life chemical or a life timeline; Class Q is a habitable class capped at microbial life (GEN.147).
- The progress ETA is a ratio of decayed sums, so the early estimate of a run on several workers is no longer up to twice too long (PERF.37).
- The phenomenon render uses `THREE.Timer` in place of the deprecated `THREE.Clock` (MAP.144).
- The macOS update daemon's plist is well-formed XML again (OPS.32).
- The settings-file download answers a plain 403 to non-admins like the other script endpoints.
