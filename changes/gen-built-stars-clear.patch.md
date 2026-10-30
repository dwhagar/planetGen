**Headline:** A system of three or more stars no longer puts one star inside another.

### Fixed
- **A third (or later) star could sit inside another star.** A drawn tree of stars is built again from each member's birth mass and the system's age, so a member the tree had evolved (a white dwarf in a post-common-envelope pair) could come out a giant wider than its own orbit, with its partner or a third star inside it. `StarSystem` now checks the stars as built against the orbits (`multistar.built_tree_clears`) and draws the tree again when one does not clear (about 1 forced triple in 150, 1 forced quadruple in 50).
  generation-output: changed
