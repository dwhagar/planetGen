**Headline:** The four-star systems the generator draws are now the ones the quadruple stability network calls stable: before, 44 percent of the 3+1 trees and 22 percent of the 2+2 trees it drew were refused by it.

### Changed

- The tree draw (`generation/multistar.py`, `_Builder.hierarchy_stable`) asks `hierarchy_stability.tree_is_stable` of every drawn tree of four stars or more, after the overlap screen, and draws the tree again (`tree.restart_unstable_hierarchy`) when a level refuses: each triple level at its own criterion, a 3+1 (an inner pair, a star, a star) and a 2+2 (two pairs) by the network of the 2023 quadruple paper [R4575] with its axis margin. The pair draw tested each orbit only as a triple against the orbit it goes round, which cannot see the outer star pumping the middle orbit. The bound on the restarts is the existing one: a shape no draw of which passes loses its lightest loose member.
- Measured on 160 four-star systems drawn on main: of 97 3+1 trees 43 were refused and of 63 2+2 trees 14; after the change 0 of 78 and 1 of 82 (the one left has its close pair's eccentricity set to zero after the tree is drawn, `system.py` line 1156, which the draw does not see; filed as a follow-up).
- Systems of three stars or fewer are not asked, and draw the same orbits as before except where the shared random stream moved.
- The wiring register loses the GEN.1016 lines for `quadruple_stability`, `hierarchy_stability` and their data files, which are read now, and the `clearance` line, which is the independent verifier by design.

generation-output: changed
