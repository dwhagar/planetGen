**Headline:** A star system with hundreds of bodies no longer runs out of recursion when its planets are given life.

### Fixed

- GEN.1734: the trial copy of a planet that the life pass makes (`life._trial_copy`) deep-copied the moons the planet was born with (`_born_moons`), and each body holds the siblings of its host, so one copy walked every body of the system in a single chain: a planetary nebula's system of a few hundred bodies raised `RecursionError` (`test_census_component_wholesphere` on main). A list, set or dict that holds bodies is now shared with the copy, like `moons` and the star already were; one with no body in it is still copied. No generated value changes.
