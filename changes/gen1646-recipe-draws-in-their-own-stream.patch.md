**Headline:** A world's second climate branch now holds the same volatile inventory as its first, as GEN.1646 meant.

### Fixed

- GEN.1646: the two branches of a bistable world (oddity W20) repeated the design constants and the late accretion but not the volatile inventory, in 6 of 11 branch worlds on main: the recipe rows (the c-after row is "after cooling" for a warm ocean world and a fallback row for a frozen one) draw a different number of values on each branch, and the inventory constants (water, carbon dioxide, nitrogen, sulfur, chlorine) are drawn after them from the same stream, so they shifted. The recipe now draws in a stream of its own (one value taken from the mixture stream for its seed), so what follows it is the same draw on both branches. Every mixture draw after the recipe moves.

generation-output: changed
generation-affects: fill
