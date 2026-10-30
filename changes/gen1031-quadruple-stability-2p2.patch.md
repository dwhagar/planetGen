**Headline:** A 2+2 quadruple star system can now be tested for stability with a published classifier, as the 3+1 one already can.

### Changed
- `physics/quadruple_stability.py` (GEN.1031) also holds the 2+2 quadruple stability network of Vynatheya, Mardling and Hamers 2023 as plain numbers and a numpy forward pass (functions with the suffix `_2p2`, no machine-learning library). On my N-body check of 185 systems drawn like the paper's it agreed with the integration 91 percent of the time against 82 percent for the nested-triple test. It is lenient near the stability boundary (it calls some weakly chaotic systems stable), which the design note says. Nothing in generation calls it yet (GEN.1016 wires it), so no generated world changes.
- `scripts/build_quadruple_stability.py` now takes `3p1` or `2p2` and the model file; registry rows QT1 to QT4 (section 6.66) and the design note `docs/design/systems/quadruple-stability-2p2.md`.
