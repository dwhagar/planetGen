### Added

- `galaxyGeometry.sectors_along_segment` lists every sector a straight
  segment passes through, in the order the segment enters them, exact at
  sector faces and edges (NAV.38). Its browser twin is
  `galaxyprisms.sectorsAlongSegment`, and tests check the two agree.
  Course planning (NAV.12) uses it to flag hops through unknown space.
