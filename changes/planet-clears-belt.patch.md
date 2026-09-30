### Fixed
- **A planet's Hill sphere could overlap the asteroid belt inside it.**
  A planet after a belt now keeps 5 Hill radii clear of the belt's outer
  edge, the same rule a belt after a planet already followed, and a
  planet moved to fix spacing now gets its Hill radius recomputed for
  its new distance (TODO item 43).
