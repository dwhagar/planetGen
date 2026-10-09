### Added
- A route now reports its longest hop and flags each hop whose line crosses sectors that have not been generated as unknown space (NAV.12). `/api/nav` returns `route.hops` and `route.longest_hop_ly`, and the NAV page states the longest hop and the unknown-space jumps.

### Changed
- Two systems in one sector are routed by their nearest stars even when those lie in the sector next door.
