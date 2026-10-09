### Fixed
- Hypervelocity stars now move in a straight line (position plus velocity times time) when the orbit update runs, keeping their velocity, instead of being turned about the galactic axis like bound objects (GEN.137). `phenomenon_scatter` gains `epoch_unix` (schema v69), the time a scattered star's position holds at, and the built system inherits it.
- The sector search reached by a radius now covers every cell that touches the sphere, not only cells whose centers lie inside it (`cells_touching_sphere`, NAV.53); the center-based listing's misses are documented.
