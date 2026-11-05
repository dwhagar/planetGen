**Headline:** The generator can now say how many of the stars a crowded cell leaves unmade a straight course would pass through.

### Added
- **The expected unresolved keep-out spheres of a segment.** `galaxy/unresolved_crossings.py` (`UnresolvedField.crossings`) gives, for a straight segment, the mean number of keep-out spheres of the unresolved stars it crosses, the length it spends in capped cells and their unresolved systems. A star's keep-out is its galactic Hill radius from the mean unresolved system mass and the galaxy's own mass inside the point. A segment that stays beyond the capped region's reach (`capped_reach_pc`) reads no cell. The course and the Navigation page do not use it yet (NAV.75).
