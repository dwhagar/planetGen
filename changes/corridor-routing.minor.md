### Changed
- A route across sectors now searches only the systems near the straight line between its ends (a corridor that widens when the route has long hops) with A*, instead of rebuilding the nearest-neighbour graph over every placed system, so routing scales to large galaxies (NAV.10). `scripts/bench_nav.py` times it on a synthetic galaxy of any size.

### Added
- `corridor.objects_near_segment`: every generated system, star and phenomenon within a distance of a line segment, in order along it (for NAV.6's course steering).
