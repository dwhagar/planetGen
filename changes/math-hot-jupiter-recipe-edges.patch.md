**Headline:** The temperature edges between a hot gas giant's atmosphere regimes now agree with the equilibrium solver the rest of the atmosphere chain uses.

### Changed
- **The edges of the hot-Jupiter recipe rows h1 to h3 are refreshed from the one equilibrium solver (GEN.723).** They were computed with the first equilibrium module, deleted in GEN.672: the water edge sat 9 to 19 K above the solver (0.4 to 0.6 percent) and the methane edge 1 to 2 K below it. A hot gas giant near an edge now takes one regime or the other with a slightly different probability, so a few drawn atmospheres change; the 100 K softness of the edges stays as it was.

generation-output: changed
